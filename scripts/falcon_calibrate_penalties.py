"""Barrido de sensibilidad de penalties para QAOA (spec §8, docs/SPEC_IMPLEMENTACION_QUBO.md).

Motivo: con penalties = 10·J_scale el landscape queda mal condicionado (el objetivo ~0.1
frente a penalties ~1 tras normalizar) y QAOA de poca profundidad no resuelve el subespacio
factible (ver docs/ANALISIS_DP_Y_RESULTADOS.md, pitfall Georgia). Este script barre un
multiplicador de penalty sobre la instancia chica (debug T5/L3), mide para QAOA:
- SRS_qaoa y AR vs el optimo exacto (DP),
- factibilidad del decodificado y del estado MAS probable,
- probabilidad total en el subespacio factible (cuanto concentra QAOA en factibles).

Se corre una vez en la instancia chica; el mejor multiplicador se fija como default y se
reusa en instancias grandes (spec §8). Usar el venv cuantico: .venv-quantum/bin/python.
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
import falcon_qaoa as qa
import falcon_qubo as fq


def _base_penalties(S0, dS, S_min, weights, T, B):
    """Penalties base (mult=1) = las auto-escaladas del builder (10·J_scale)."""
    Jscale = fq._historical_Jscale(S0, dS, S_min, weights, T)
    return {
        "P_onehot": 10.0 * Jscale,
        "P_R": 10.0 * Jscale,
        "P_bal_soft": 10.0 * Jscale / max(B * B, 1e-30),
        "P_bal_slack": 10.0 * Jscale,
    }


def _feasible_mass(probs, vi, levels, R_obs, dS, S0, *,
                   S_min, S_max, u_max, B, weights):
    """Prob total en estados factibles + factibilidad del estado mas probable."""
    import falcon_solvers as sv

    # Evaluar factibilidad es caro por estado; muestreamos el top-256 por probabilidad.
    order = np.argsort(probs)[::-1][:256]
    mass = 0.0
    top_feasible = None
    for rank, idx in enumerate(order):
        lv = np.empty(vi.T, dtype=int)
        ok_onehot = True
        for t in range(vi.T):
            wk = [(int(idx) >> vi.idx(t, l)) & 1 for l in range(vi.L)]
            if sum(wk) != 1:
                ok_onehot = False
            lv[t] = int(np.argmax(wk))
        dv = sv.decode_and_verify(lv, vi, levels, R_obs, dS, S0, S_min=S_min, S_max=S_max,
                                  u_max=u_max, B=B, weights=weights)
        feas = ok_onehot and dv["feasible"]
        if feas:
            mass += float(probs[idx])
        if rank == 0:
            top_feasible = feas
    return mass, top_feasible


def sweep(T=5, L=3, mults=(0.1, 0.3, 1.0, 3.0, 10.0), p=1, restarts=5, seed=42, maxiter=150):
    df = fd.build_weekly_benchmark(write=False)
    c = fc.load_official_constants()
    Rall = df["R_obs_m3_week"].to_numpy()
    R, dS = Rall[:T], df["DeltaS_obs_m3"].to_numpy()[:T]
    S0 = df.attrs["S0_m3"]
    pr = fc.instance_params(Rall, T, L)
    w = fc.compute_weights(T, c["S_min_m3"], pr["u_max"])
    B = c["eta"] * float(R.sum())
    half = (L - 1) // 2
    base = _base_penalties(S0, dS, c["S_min_m3"], w, T, B)

    dp = bl.dp_optimal(R, dS, S0, S_min=c["S_min_m3"], S_max=c["S_max_m3"],
                       delta_u=pr["delta_u"], L=L, weights=w, B=B)

    print(f"instancia debug T{T}/L{L}  DP*={dp['SRS_star']:.6f}  (mult=1 == default 10·Jscale)\n")
    print(f"{'mult':>5} {'SRS_qaoa':>11} {'AR':>7} {'feas':>5} {'topFeas':>7} "
          f"{'P(feas)':>8} {'prob*':>9} {'iters':>5}")
    rows = []
    for m in mults:
        pen = {k: m * v for k, v in base.items()}
        cfg = fcfg.FalconConfig(T=T, L=L, penalties=pen)
        Q, const, vi, meta = fq.build_qubo(cfg, R, dS, S0, S_min=c["S_min_m3"],
                                           delta_u=pr["delta_u"], levels=pr["levels"],
                                           weights=w, B=B)
        out = qa.run_qaoa(Q, const, meta, vi, pr["levels"], R, dS, S0, S_min=c["S_min_m3"],
                          S_max=c["S_max_m3"], u_max=pr["u_max"], B=B, weights=w, half=half,
                          p=p, restarts=restarts, seed=seed, maxiter=maxiter)
        d = out["decoded"]
        # reconstruir probs del mejor set de params para medir masa factible
        singles, pairs, _ = qa._ising_terms(Q)
        sim = qa._Sim(meta["n_qubits"], singles, pairs, p)
        psi = sim.statevector(out["betas"], out["gammas"])
        probs = np.abs(psi) ** 2
        mass, topf = _feasible_mass(probs, vi, pr["levels"], R, dS, S0,
                                    S_min=c["S_min_m3"], S_max=c["S_max_m3"], u_max=pr["u_max"],
                                    B=B, weights=w)
        ar = d["SRS"] / dp["SRS_star"] if dp["SRS_star"] else float("nan")
        print(f"{m:>5.1f} {d['SRS']:>11.6f} {ar:>7.4f} {str(d['feasible']):>5} "
              f"{str(topf):>7} {mass:>8.4f} {d['prob']:>9.2e} {out['n_iter']:>5}")
        rows.append({"mult": m, "SRS": d["SRS"], "AR": ar, "feasible": d["feasible"],
                     "top_feasible": topf, "feasible_mass": mass, "prob": d["prob"]})
    return rows, dp["SRS_star"]


if __name__ == "__main__":
    sweep()
