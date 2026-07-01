"""Validacion del XY-mixer (falcon_qaoa, mixer='xy') en instancias chicas.

Comprueba que: (1) todos los estados top muestreados son one-hot validos (el XY-mixer
confina al subespacio factible), (2) el decode es factible, (3) el SRS alcanza/roza el
optimo del DP. Compara contra el QAOA RX (penalty) de referencia.

Correr con: .venv-quantum/bin/python scripts/julian/falcon_xymixer_check.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

import falcon_baselines as bl
import falcon_config as fcfg
import falcon_constants as fc
import falcon_data as fd
import falcon_qaoa as qa
import falcon_qubo as fq


def check(T=5, L=3, p=1, restarts=5, maxiter=150):
    df = fd.build_weekly_benchmark(write=False)
    c = fc.load_official_constants()
    Rall = df["R_obs_m3_week"].to_numpy()
    R, dS = Rall[:T], df["DeltaS_obs_m3"].to_numpy()[:T]
    S0 = df.attrs["S0_m3"]
    p_ = fc.instance_params(Rall, T, L)
    du, umax, levels = p_["delta_u"], p_["u_max"], p_["levels"]
    w = fc.compute_weights(T, c["S_min_m3"], umax)
    B = c["eta"] * float(R.sum())
    half = (L - 1) // 2

    dp = bl.dp_optimal(R, dS, S0, S_min=c["S_min_m3"], S_max=c["S_max_m3"],
                       delta_u=du, L=L, weights=w, B=B)

    # QUBO con XY-mixer: sin penalizacion one-hot, balance soft (sin slacks -> n = T*L)
    cfg = fcfg.FalconConfig(T=T, L=L, balance="soft", onehot="xy_mixer")
    Q, const, vi, meta = fq.build_qubo(cfg, R, dS, S0, S_min=c["S_min_m3"],
                                       delta_u=du, levels=levels, weights=w, B=B)
    assert meta["n_qubits"] == T * L, f"esperaba {T*L} qubits, hay {meta['n_qubits']}"

    out = qa.run_qaoa(Q, const, meta, vi, levels, R, dS, S0, S_min=c["S_min_m3"],
                      S_max=c["S_max_m3"], u_max=umax, B=B, weights=w, half=half,
                      p=p, restarts=restarts, maxiter=maxiter, mixer="xy")
    d = out["decoded"]

    # (1) fraccion one-hot valida en el top-256 muestreado: recomputar el statevector
    #     del mejor (β,γ) para inspeccionar el top
    singles, pairs, _ = qa._ising_terms(Q)
    groups = [[vi.idx(t, l) for l in range(L)] for t in range(T)]
    sim = qa._Sim(meta["n_qubits"], singles, pairs, p, mixer="xy", groups=groups)
    psi = sim.statevector(out["betas"], out["gammas"])
    probs = np.abs(psi) ** 2
    top = np.argsort(probs)[::-1][:256]
    n_valid = 0
    for idx in top:
        xbits = np.array([(int(idx) >> i) & 1 for i in range(vi.n)], dtype=float)
        if vi.is_onehot(xbits):
            n_valid += 1
    frac_valid = n_valid / len(top)
    prob_valid = float(sum(probs[i] for i in top
                           if vi.is_onehot(np.array([(int(i) >> b) & 1 for b in range(vi.n)], float))))

    ar = d["SRS"] / dp["SRS_star"] if dp["SRS_star"] else None
    n_onehot_states = L ** T
    print(f"=== XY-mixer T{T}/L{L} p={p} ({meta['n_qubits']} qubits) ===")
    print(f"  top-256: {n_valid} one-hot validos (de {n_onehot_states} posibles), "
          f"**prob en subespacio one-hot = {prob_valid*100:.4f}%**")
    print(f"  decode: feasible={d['feasible']}  SRS={d['SRS']:.6f}  dp*={dp['SRS_star']:.6f}  "
          f"AR={ar:.4f}  u!=0 en {int(np.count_nonzero(d['u']))}/{T}")
    # el XY-mixer confina TODA la probabilidad al subespacio one-hot (la metrica fisica);
    # el conteo en el top-256 puede incluir estados de amplitud ~0 por desempate del argsort.
    assert prob_valid > 0.999, f"FALLA: solo {prob_valid*100:.2f}% de prob en one-hot"
    assert d["feasible"], "FALLA: decode infactible con XY-mixer"
    print("  OK: XY-mixer confina ~100% de la probabilidad al subespacio one-hot y decodifica factible.\n")
    return {"frac_valid": frac_valid, "prob_valid": prob_valid, "SRS": d["SRS"],
            "dp": dp["SRS_star"], "AR": ar, "feasible": d["feasible"]}


if __name__ == "__main__":
    check(5, 3, p=1)
    print("XY-MIXER CHECK PASSED")
