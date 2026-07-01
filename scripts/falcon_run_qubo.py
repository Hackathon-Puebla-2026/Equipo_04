"""Corre los solvers del QUBO y los registra con el writer estandarizado.

- small (T12/L3): exhaustivo sobre el QUBO (debe == dp).
- medium (T26/L5): simulated annealing sobre el QUBO (factible, heuristico).
Registra en results/ con los campos cuanticos del esquema (encoding, n_qubits, penalties).

Uso: python scripts/falcon_run_qubo.py
"""
from __future__ import annotations

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import falcon_baselines as bl
import falcon_config as fcfg
import falcon_constants as fc
import falcon_data as fd
import falcon_qubo as fq
import falcon_results as res
import falcon_solvers as sv
import falcon_storage as st


def run(T, L, method):
    df = fd.build_weekly_benchmark(write=False)
    c = fc.load_official_constants()
    Rall = df["R_obs_m3_week"].to_numpy()
    R, dS = Rall[:T], df["DeltaS_obs_m3"].to_numpy()[:T]
    S0 = df.attrs["S0_m3"]
    p = fc.instance_params(Rall, T, L)
    w = fc.compute_weights(T, c["S_min_m3"], p["u_max"])
    B = c["eta"] * float(R.sum())
    half = (L - 1) // 2
    cfg = fcfg.FalconConfig(T=T, L=L)
    Q, const, vi, meta = fq.build_qubo(cfg, R, dS, S0, S_min=c["S_min_m3"],
                                       delta_u=p["delta_u"], levels=p["levels"], weights=w, B=B)

    # referencias (historical, threshold pura, dp)
    dp = bl.dp_optimal(R, dS, S0, S_min=c["S_min_m3"], S_max=c["S_max_m3"],
                       delta_u=p["delta_u"], L=L, weights=w, B=B)
    uh = bl.historical(T); Sh = st.simulate_storage(S0, dS, uh)
    import falcon_srs as srs
    srs_h = srs.compute_srs(srs.compute_costs(Sh, uh, c["S_min_m3"]), w)
    ut = bl.threshold_rule(R, dS, S0, c["S_min_m3"], p["delta_u"])
    srs_t = srs.compute_srs(srs.compute_costs(st.simulate_storage(S0, dS, ut), ut, c["S_min_m3"]), w)
    refs = {"historical": srs_h, "threshold": srs_t, "dp": dp["SRS_star"]}

    t0 = time.perf_counter()
    if method == "exhaustive":
        out = sv.exhaustive_qubo(Q, const, vi, meta, p["levels"], R, B, half)
        energy = out["energy"]
    else:
        out = sv.simulated_annealing_qubo(Q, const, vi, meta, half, n_iter=60_000)
        energy = out["energy"]
    rt = time.perf_counter() - t0

    dv = sv.decode_and_verify(out["lv"], vi, p["levels"], R, dS, S0, S_min=c["S_min_m3"],
                              S_max=c["S_max_m3"], u_max=p["u_max"], B=B, weights=w)
    S = st.simulate_storage(S0, dS, dv["u"])
    ar = dv["SRS"] / dp["SRS_star"] if dp["SRS_star"] != 0 else None   # ratio vs optimo (dp)
    solver = {"encoding": meta["encoding"], "n_qubits": meta["n_qubits"], "p_depth": None,
              "energy": float(energy), "approximation_ratio": ar, "seed": 42,
              "simulator": f"numpy_{method}", "penalties": meta["penalties"]}

    res.record_run(method=f"qubo_{method}", instance=res.instance_label(T, L), T=T, L=L,
                   params=p, weights=w, constants=c, B=B, u=dv["u"], S=S, costs=dv["costs"],
                   srs=dv["SRS"], feasibility={"feasible": dv["feasible"], "violations": dv["violations"]},
                   runtime_seconds=rt, references=refs, solver=solver)
    print(f"[{res.instance_label(T,L):7s} T={T} L={L}] qubo_{method}: SRS={dv['SRS']:.6f} "
          f"dp={dp['SRS_star']:.6f} feasible={dv['feasible']} n_qubits={meta['n_qubits']} ({rt:.1f}s)")


def main():
    run(12, 3, "exhaustive")     # debe == dp
    run(26, 5, "sa")             # factible, heuristico
    print(f"\nRegistrado en {res.RESULTS_ROOT}/runs_summary.csv")


if __name__ == "__main__":
    main()
