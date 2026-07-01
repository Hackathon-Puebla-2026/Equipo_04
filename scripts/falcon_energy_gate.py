"""Gate de energia: valida que el QUBO codifica el objetivo correcto.

Para bitstrings FACTIBLES (one-hot valido, R>=0, |Σu|<=B con el slack bien puesto),
las penalizaciones valen 0, asi que debe cumplirse:
    qubo_energy(Q,x,const)*scale  ==  J  ==  -SRS
donde J = w1 Ccrit + w2 Cdev + w3 Csmooth (via falcon_srs sobre el u decodificado).
Se prueba en el optimo del DP y en varios calendarios factibles aleatorios.

Uso: python scripts/falcon_energy_gate.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

import falcon_baselines as bl
import falcon_config as fcfg
import falcon_constants as fc
import falcon_data as fd
import falcon_qubo as fq
import falcon_srs as srs
import falcon_storage as st


def _set_balance_slack(x, vi, meta, half):
    """Fija los bits de slack de balance para que M+s=M_cap (penalizacion 0)."""
    if meta["n_balance_slack"] == 0:
        return x
    lvls = vi.decode_levels(x[:vi.n])
    M = int(sum(l - half for l in lvls))
    s = meta["M_cap"] - M                       # en [0, 2*M_cap] si factible
    for r in range(meta["n_balance_slack"]):
        x[vi.n + r] = float((s >> r) & 1)
    return x


def gate(T=12, L=3, n_random=8, tol=1e-4):
    df = fd.build_weekly_benchmark(write=False)
    c = fc.load_official_constants()
    R = df["R_obs_m3_week"].to_numpy()[:T]
    dS = df["DeltaS_obs_m3"].to_numpy()[:T]
    S0 = df.attrs["S0_m3"]
    p = fc.instance_params(df["R_obs_m3_week"].to_numpy(), T, L)
    w = fc.compute_weights(T, c["S_min_m3"], p["u_max"])
    levels = p["levels"]; du = p["delta_u"]; half = (L - 1) // 2
    B = c["eta"] * float(R.sum())

    cfg = fcfg.FalconConfig(T=T, L=L)   # defaults: onehot, storage drop, c_crit soft, balance slack
    Q, const, vi, meta = fq.build_qubo(cfg, R, dS, S0, S_min=c["S_min_m3"],
                                       delta_u=du, levels=levels, weights=w, B=B)
    n = meta["n_qubits"]
    print(f"QUBO T={T} L={L}: n_qubits={n} (dec={meta['n_decision']}, "
          f"bal_slack={meta['n_balance_slack']}), scale={meta['scale']:.3e}, "
          f"forbidden={meta['n_forbidden_levels']}")

    # candidatos factibles: optimo DP + aleatorios
    dp = bl.dp_optimal(R, dS, S0, S_min=c["S_min_m3"], S_max=c["S_max_m3"],
                       delta_u=du, L=L, weights=w, B=B)
    dp_lvls = np.round(dp["u_star"] / du).astype(int) + half
    candidates = [dp_lvls]
    rng = np.random.default_rng(0)
    while len(candidates) < 1 + n_random:
        lv = rng.integers(0, L, size=T)
        u = levels[lv]
        if (R + u < -1e-6).any():                 # R>=0
            continue
        if abs(u.sum()) > B + 1e-6:               # balance
            continue
        candidates.append(lv)

    max_err = 0.0
    for k, lv in enumerate(candidates):
        x = vi.encode_levels(lv)
        x = np.concatenate([x, np.zeros(meta["n_balance_slack"])])
        x = _set_balance_slack(x, vi, meta, half)
        # energia del QUBO (des-normalizada)
        e = fq.qubo_energy(Q, x, const) * meta["scale"]
        # objetivo directo (oraculo)
        u = vi.decode_u(x[:vi.n], levels)
        S = st.simulate_storage(S0, dS, u)
        costs = srs.compute_costs(S, u, c["S_min_m3"])
        J = w["w1"] * costs["Ccrit"] + w["w2"] * costs["Cdev"] + w["w3"] * costs["Csmooth"]
        srs_val = srs.compute_srs(costs, w)
        err = abs(e - J)
        max_err = max(max_err, err)
        rel = err / max(abs(J), 1e-30)
        tag = "DP*" if k == 0 else f"rand{k}"
        assert rel < tol, f"[{tag}] energia {e:.6e} != J {J:.6e} (rel {rel:.2e})"
        # ademas energia == -SRS
        assert abs(e + srs_val) < tol * max(abs(srs_val), 1e-30) + 1e-12
    print(f"GATE OK: energia == J == -SRS en {len(candidates)} bitstrings factibles "
          f"(max |Δ| = {max_err:.2e})")


if __name__ == "__main__":
    gate(12, 3)
    gate(26, 5)
    print("\nENERGY GATE PASSED")
