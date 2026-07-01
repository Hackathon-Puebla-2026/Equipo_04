"""Solvers sobre el QUBO: exhaustivo (valida vs DP), simulated annealing, y
decode+verify (post-seleccion de factibilidad).

Trabajan sobre el espacio estructurado one-hot (indices de nivel por semana) para
garantizar one-hot valido; el slack de balance se fija al valor optimo (s=M_cap-M)
antes de evaluar la energia del QUBO. Asi la energia normalizada * scale == J == -SRS
para factibles (validado por el gate de energia).
"""
from __future__ import annotations

import numpy as np

import falcon_qubo as fq
import falcon_srs as srs
import falcon_storage as st


def _x_from_levels(lv, vi, meta, half):
    """Construye el bitvector (bits de decision + slack de balance optimo) desde indices de nivel.

    Los bits de decision se arman con el encoding (one-hot o binary) via `vi.encode_levels`.
    """
    x = np.zeros(meta["n_qubits"], dtype=float)
    x[:vi.n] = vi.encode_levels([int(l) for l in lv])
    if meta["n_balance_slack"] > 0:
        M = int(sum(int(l) - half for l in lv))
        s = meta["M_cap"] - M
        s = min(max(s, 0), 2 * meta["M_cap"])          # clamp por si infactible
        for r in range(meta["n_balance_slack"]):
            x[vi.n + r] = float((s >> r) & 1)
    return x


def _srs_of_levels(lv, levels, dS, S0, S_min, weights):
    u = levels[np.asarray(lv, dtype=int)]
    S = st.simulate_storage(S0, dS, u)
    costs = srs.compute_costs(S, u, S_min)
    return srs.compute_srs(costs, weights), u, S, costs


def decode_and_verify(lv, vi, levels, R_obs, dS, S0, *, S_min, S_max, u_max, B, weights):
    """Decodifica indices de nivel -> u, verifica factibilidad y calcula SRS."""
    srs_val, u, S, costs = _srs_of_levels(lv, levels, dS, S0, S_min, weights)
    chk = st.check_constraints(R_obs, u, S, S_max, u_max, B)
    return {"u": u, "SRS": srs_val, "costs": costs, "feasible": chk["feasible"],
            "violations": chk["violations"]}


def exhaustive_qubo(Q, const, vi, meta, levels, R_obs, B, half,
                    max_combos=2_000_000, chunk=50_000):
    """Minimiza la energia del QUBO sobre one-hot factibles (release+balance), vectorizado.

    Enumera `L^T` combinaciones de nivel; filtra factibilidad; evalua energia por lotes
    (`X @ Q`). Devuelve el mejor (min energia == max SRS). None si `L^T > max_combos`.
    """
    T, L = vi.T, vi.L
    N = L ** T
    if N > max_combos:
        return None
    levels = np.asarray(levels, dtype=float)
    R_obs = np.asarray(R_obs, dtype=float)

    best = {"energy": np.inf, "lv": None}
    idx_all = np.arange(N)
    for start in range(0, N, chunk):
        idx = idx_all[start:start + chunk]
        b = len(idx)
        # combos: (b, T) indices de nivel via contador base-L
        combos = np.empty((b, T), dtype=np.int64)
        tmp = idx.copy()
        for t in range(T):
            combos[:, t] = tmp % L
            tmp //= L
        u = levels[combos]                        # (b,T)
        # factibilidad release y balance (storage nunca ata)
        feas = ((R_obs[None, :] + u) >= -1e-6).all(axis=1) & (np.abs(u.sum(axis=1)) <= B + 1e-6)
        if not feas.any():
            continue
        combos_f = combos[feas]
        bf = len(combos_f)
        # construir X (one-hot + slack) y energia por lote
        X = np.zeros((bf, meta["n_qubits"]))
        rows = np.arange(bf)[:, None]
        cols = combos_f + np.arange(T)[None, :] * L   # idx(t,l)=t*L+l
        X[rows, cols] = 1.0
        if meta["n_balance_slack"] > 0:
            M = (combos_f - half).sum(axis=1)
            s = np.clip(meta["M_cap"] - M, 0, 2 * meta["M_cap"]).astype(np.int64)
            for r in range(meta["n_balance_slack"]):
                X[:, vi.n + r] = ((s >> r) & 1).astype(float)
        energies = np.einsum("bi,ij,bj->b", X, Q, X) + const
        j = int(np.argmin(energies))
        if energies[j] < best["energy"]:
            best = {"energy": float(energies[j]), "lv": combos_f[j].copy()}

    srs_val = -best["energy"] * meta["scale"]     # energia_norm*scale = J = -SRS
    return {"lv": best["lv"], "SRS_star": srs_val, "energy": best["energy"]}


def simulated_annealing_qubo(Q, const, vi, meta, half, *, n_iter=40_000,
                             T0=1.0, Tf=1e-3, seed=42):
    """SA sobre indices de nivel (cambia el nivel de una semana por paso), min energia QUBO."""
    rng = np.random.default_rng(seed)
    T, L = vi.T, vi.L
    lv = rng.integers(0, L, size=T)
    x = _x_from_levels(lv, vi, meta, half)
    e = fq.qubo_energy(Q, x, const)
    best_lv, best_e = lv.copy(), e
    alpha = (Tf / T0) ** (1.0 / max(n_iter, 1))
    temp = T0
    for _ in range(n_iter):
        t = int(rng.integers(0, T))
        new_l = int(rng.integers(0, L))
        if new_l == lv[t]:
            temp *= alpha
            continue
        old_l = lv[t]
        lv[t] = new_l
        x = _x_from_levels(lv, vi, meta, half)
        e_new = fq.qubo_energy(Q, x, const)
        if e_new <= e or rng.random() < np.exp(-(e_new - e) / max(temp, 1e-12)):
            e = e_new
            if e < best_e:
                best_e, best_lv = e, lv.copy()
        else:
            lv[t] = old_l
        temp *= alpha
    return {"lv": best_lv, "energy": float(best_e), "SRS_star": -best_e * meta["scale"]}
