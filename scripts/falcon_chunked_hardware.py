"""Planificación (en simulador) de la QAOA chunked para correr en hardware IBM.

`plan_blocks` replica el loop por etapas (como `falcon_chunked.staged_solve`, pero SIN editar
ese módulo de la otra instancia) usando XY-mixer, y CAPTURA por bloque todo lo necesario para
armar/enviar el circuito de hardware: `Q, const, meta, vi, groups, betas, gammas, S0_blk, B_blk`.

Idea: el simulador es barato → entrena β/γ y fija las fronteras de storage (carry) por bloque.
Luego cada bloque se envía a hardware con esos params (jobs independientes → paralelos). El SRS
global final se evalúa exacto sobre el `u` decodificado del HARDWARE (ver runner).
"""
from __future__ import annotations

import numpy as np

import falcon_config as fcfg
import falcon_constants as fc
import falcon_qaoa as qa
import falcon_qubo as fq
import falcon_storage as st


def plan_blocks(weekly_release_m3, deltaS_all, S0, T, *, S_min, S_max, L, block_size,
                eta=0.10, p=1, restarts=3, maxiter=100, seed=42, link_smooth=True):
    """Entrena la chunked QAOA (XY-mixer) en sim y devuelve los builds por bloque + params globales."""
    Rall = np.asarray(weekly_release_m3, dtype=float)
    R = Rall[:T]
    dS = np.asarray(deltaS_all, dtype=float)[:T]
    half = (L - 1) // 2

    pr = fc.instance_params(Rall, T, L)
    delta_u, u_max, levels = pr["delta_u"], pr["u_max"], np.asarray(pr["levels"], dtype=float)
    weights = fc.compute_weights(T, S_min, u_max)
    B_global = eta * float(R.sum())

    S_cur = float(S0)
    k_prev = None
    blocks = []
    for start in range(0, T, block_size):
        end = min(start + block_size, T)
        Rb, dSb, Tb = R[start:end], dS[start:end], end - start
        B_blk = eta * float(Rb.sum())                      # eta_local (balance global garantizado)
        u_prev = (k_prev * delta_u) if (link_smooth and k_prev is not None) else None

        # XY-mixer: sin penalización one-hot ni slacks -> n = Tb*L, todo one-hot válido
        cfg = fcfg.FalconConfig(T=Tb, L=L, balance="soft", onehot="xy_mixer")
        Q, const, vi, meta = fq.build_qubo(cfg, Rb, dSb, S_cur, S_min=S_min, delta_u=delta_u,
                                           levels=levels, weights=weights, B=B_blk, u_prev=u_prev)
        groups = [[vi.idx(t, l) for l in range(L)] for t in range(Tb)]

        out = qa.run_qaoa(Q, const, meta, vi, levels, Rb, dSb, S_cur, S_min=S_min, S_max=S_max,
                          u_max=u_max, B=B_blk, weights=weights, half=half, p=p, restarts=restarts,
                          maxiter=maxiter, seed=seed, mixer="xy")
        u_blk = np.asarray(out["decoded"]["u"], dtype=float)   # decode del SIM -> fija boundary
        S_blk = st.simulate_storage(S_cur, dSb, u_blk)

        blocks.append({"start": start, "end": end, "Tb": Tb, "B_blk": B_blk,
                       "Q": Q, "const": const, "meta": meta, "vi": vi, "groups": groups,
                       "betas": out["betas"], "gammas": out["gammas"],
                       "S0_blk": S_cur, "Rb": Rb, "dSb": dSb,
                       "sim_feasible": bool(out["decoded"]["feasible"]),
                       "sim_lv": [int(x) for x in out["decoded"]["lv"]]})
        S_cur = float(S_blk[-1])
        k_prev = int(round(u_blk[-1] / delta_u)) if len(u_blk) else k_prev

    return {"blocks": blocks, "params": pr, "weights": weights, "levels": levels,
            "delta_u": delta_u, "u_max": u_max, "half": half, "B_global": B_global,
            "R": R, "dS": dS, "S0": float(S0), "T": T, "L": L, "block_size": block_size}
