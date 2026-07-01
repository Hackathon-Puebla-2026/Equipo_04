"""Solver por etapas / horizonte recedente (chunking) sobre el QUBO de Falcon.

Reimplementacion limpia del enfoque "por bloques" (parte T en tramos que entran en el
simulador) usando NUESTRO QAOA que funciona (`falcon_qaoa.run_qaoa`, statevector) o el DP
exacto por bloque, con las mejoras aprendidas (ver docs/SPEC_IMPLEMENTACION_QUBO.md E2):

- **Parametros oficiales GLOBALES** (Δu, niveles, pesos) sobre la ventana T completa, usados
  en TODOS los bloques (no se recomputan por bloque).
- **Carry de storage** entre bloques (`S0` del bloque = storage final del anterior).
- **C_smooth linkeado** entre bloques (`u_prev`/`k_prev_init`) si `link_smooth`.
- **Balance por `eta_local`** (`B_blk = η·ΣR_blk`), que garantiza el balance global por
  desigualdad triangular; opcional `global_greedy`.
- **Sin fallback silencioso**: cada bloque reporta su metodo y factibilidad; nada se enmascara.
- El SRS y la factibilidad finales se evaluan GLOBALMENTE con los modulos canonicos.

QAOA usa `balance="soft"` por bloque para que `T_blk·L` entre en statevector (p.ej. L5 y
bloques de 5 = 25 qubits); el balance global se verifica al final. El DP por bloque usa el
balance exacto de slack. Correr QAOA con `.venv-quantum/bin/python`.
"""
from __future__ import annotations

import numpy as np

import falcon_baselines as bl
import falcon_config as fcfg
import falcon_constants as fc
import falcon_qubo as fq
import falcon_srs as srs
import falcon_storage as st


def _block_budget(balance_split, eta, R_blk, B_global, running_abs):
    """Presupuesto de balance del bloque segun la estrategia de reparto."""
    if balance_split == "eta_local":
        return eta * float(R_blk.sum())
    if balance_split == "global_greedy":
        return max(0.0, B_global - running_abs)   # lo que quede del presupuesto global
    raise ValueError(f"balance_split invalido: {balance_split!r}")


def staged_solve(weekly_release_m3, deltaS_all, S0: float, T: int, *,
                 S_min: float, S_max: float, L: int, block_size: int,
                 eta: float = 0.10, solver: str = "qaoa",
                 balance_split: str = "eta_local", link_smooth: bool = True,
                 qaoa_kwargs: dict | None = None) -> dict:
    """Resuelve por etapas y devuelve u global + SRS/factibilidad global + tabla por bloque.

    `weekly_release_m3` es el record completo de release (para Δu oficial de la ventana T);
    `deltaS_all` el ΔS_obs completo. Se usan las primeras T semanas.
    `solver`: "qaoa" (nuestro statevector) o "dp" (exacto por bloque, oraculo del chunking).
    """
    qaoa_kwargs = qaoa_kwargs or {}
    Rall = np.asarray(weekly_release_m3, dtype=float)
    R = Rall[:T]
    dS = np.asarray(deltaS_all, dtype=float)[:T]
    half = (L - 1) // 2

    # --- parametros oficiales GLOBALES (sobre la ventana T completa) ---
    pr = fc.instance_params(Rall, T, L)
    delta_u, u_max, levels = pr["delta_u"], pr["u_max"], np.asarray(pr["levels"], dtype=float)
    weights = fc.compute_weights(T, S_min, u_max)
    B_global = eta * float(R.sum())

    block_balance = "soft" if solver == "qaoa" else "slack"   # soft para caber en statevector

    S_cur = float(S0)
    k_prev = None                    # nivel entero de la ultima semana del bloque anterior
    u_all: list[float] = []
    running = 0.0                    # Σu acumulado (para global_greedy)
    blocks = []

    for start in range(0, T, block_size):
        end = min(start + block_size, T)
        Rb, dSb, Tb = R[start:end], dS[start:end], end - start
        B_blk = _block_budget(balance_split, eta, Rb, B_global, abs(running))
        u_prev = (k_prev * delta_u) if (link_smooth and k_prev is not None) else None

        cfg = fcfg.FalconConfig(T=Tb, L=L, balance=block_balance)
        Q, const, vi, meta = fq.build_qubo(
            cfg, Rb, dSb, S_cur, S_min=S_min, delta_u=delta_u, levels=levels,
            weights=weights, B=B_blk, u_prev=u_prev)

        if solver == "qaoa":
            import falcon_qaoa as qa     # perezoso (requiere qiskit / .venv-quantum)
            out = qa.run_qaoa(Q, const, meta, vi, levels, Rb, dSb, S_cur, S_min=S_min,
                              S_max=S_max, u_max=u_max, B=B_blk, weights=weights, half=half,
                              **qaoa_kwargs)
            d = out["decoded"]
            u_blk = np.asarray(d["u"], dtype=float)
            row = {"method": "qaoa", "feasible": bool(d["feasible"]),
                   "n_qubits": meta["n_qubits"], "prob": d["prob"],
                   "betas": out["betas"], "gammas": out["gammas"],
                   "selected_levels": [int(x) for x in d["lv"]]}
        elif solver == "dp":
            dp = bl.dp_optimal(Rb, dSb, S_cur, S_min=S_min, S_max=S_max, delta_u=delta_u,
                               L=L, weights=weights, B=B_blk,
                               k_prev_init=(k_prev if link_smooth else None))
            u_blk = np.asarray(dp["u_star"], dtype=float)
            row = {"method": "dp", "feasible": bool(dp["feasible"]), "n_qubits": None,
                   "selected_levels": [int(round(x / delta_u)) + half for x in u_blk]}
        else:
            raise ValueError(f"solver invalido: {solver!r}")

        # eval del bloque (SRS local, informativo) + carry
        S_blk = st.simulate_storage(S_cur, dSb, u_blk)
        row.update({"start": start, "end": end, "Tb": Tb, "B_blk": B_blk,
                    "block_SRS": srs.compute_srs(srs.compute_costs(S_blk, u_blk, S_min), weights),
                    "u": u_blk.tolist()})
        blocks.append(row)
        u_all.extend(u_blk.tolist())
        running += float(u_blk.sum())
        S_cur = float(S_blk[-1])
        k_prev = int(round(u_blk[-1] / delta_u)) if len(u_blk) else k_prev

    # --- evaluacion GLOBAL con los modulos canonicos (sin doble conteo de frontera) ---
    u = np.array(u_all, dtype=float)
    S = st.simulate_storage(S0, dS, u)
    costs = srs.compute_costs(S, u, S_min)
    SRS = srs.compute_srs(costs, weights)
    chk = st.check_constraints(R, u, S, S_max, u_max, B_global)
    return {"u": u, "S": S, "SRS": SRS, "costs": costs,
            "feasible": bool(chk["feasible"]), "violations": chk["violations"],
            "blocks": blocks, "params": pr, "weights": weights, "B": B_global,
            "solver": solver, "block_size": block_size, "balance_split": balance_split,
            "link_smooth": link_smooth,
            "n_qaoa_blocks": sum(1 for b in blocks if b["method"] == "qaoa" and b["feasible"]),
            "n_blocks": len(blocks)}


def gap_vs_full(staged_result, R, dS, S0, *, S_min, S_max, L):
    """Gap del resultado por etapas contra el DP exacto GLOBAL (oraculo)."""
    w = staged_result["weights"]
    pr = staged_result["params"]
    dp_full = bl.dp_optimal(np.asarray(R)[:len(staged_result["u"])],
                            np.asarray(dS)[:len(staged_result["u"])], S0,
                            S_min=S_min, S_max=S_max, delta_u=pr["delta_u"], L=L,
                            weights=w, B=staged_result["B"])
    return {"SRS_full_dp": dp_full["SRS_star"], "SRS_staged": staged_result["SRS"],
            "gap_vs_full": dp_full["SRS_star"] - staged_result["SRS"], "dp_full": dp_full}
