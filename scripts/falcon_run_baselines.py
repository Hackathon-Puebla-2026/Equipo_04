"""Corre y registra los baselines (historical, threshold, DP exacto) para las
instancias estandar, dejando la SRS de referencia (SRS_hist) guardada en results/.

Uso: python scripts/falcon_run_baselines.py
Escribe results/runs_summary.csv + results/runs/*.json (raiz compartida).
"""
from __future__ import annotations

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import falcon_baselines as bl
import falcon_constants as fc
import falcon_data as fd
import falcon_results as res
import falcon_srs as srs
import falcon_storage as st

STANDARD_INSTANCES = [(12, 3), (26, 5), (52, 5)]


def _evaluate(u, S0, dS, S_min, weights):
    S = st.simulate_storage(S0, dS, u)
    costs = srs.compute_costs(S, u, S_min)
    return S, costs, srs.compute_srs(costs, weights)


def run_instance(df, constants, T: int, L: int) -> None:
    R = df["R_obs_m3_week"].to_numpy()
    dS = df["DeltaS_obs_m3"].to_numpy()
    S0 = df.attrs["S0_m3"]
    R_T, dS_T = R[:T], dS[:T]

    p = fc.instance_params(R, T, L)
    w = fc.compute_weights(T, constants["S_min_m3"], p["u_max"])
    B = constants["eta"] * float(R_T.sum())
    S_min, S_max = constants["S_min_m3"], constants["S_max_m3"]
    instance = res.instance_label(T, L)

    # 1) Historical (define SRS_hist = referencia)
    t0 = time.perf_counter()
    uh = bl.historical(T)
    dt_h = time.perf_counter() - t0
    Sh, ch, srs_h = _evaluate(uh, S0, dS_T, S_min, w)

    # 2) Threshold (regla pura del spec)
    t0 = time.perf_counter()
    ut = bl.threshold_rule(R_T, dS_T, S0, S_min, p["delta_u"])
    dt_t = time.perf_counter() - t0
    Stt, ct, srs_t = _evaluate(ut, S0, dS_T, S_min, w)

    # 3) DP exacto (baseline clasico fuerte / ground truth)
    t0 = time.perf_counter()
    dp = bl.dp_optimal(R_T, dS_T, S0, S_min=S_min, S_max=S_max,
                       delta_u=p["delta_u"], L=L, weights=w, B=B)
    dt_dp = time.perf_counter() - t0
    u_dp = dp["u_star"]
    Sdp = st.simulate_storage(S0, dS_T, u_dp)

    refs = {"historical": srs_h, "threshold": srs_t, "dp": dp["SRS_star"]}

    common = dict(instance=instance, T=T, L=L, params=p, weights=w,
                  constants=constants, B=B, references=refs)

    res.record_run(method="historical", u=uh, S=Sh, costs=ch, srs=srs_h,
                   feasibility=st.check_constraints(R_T, uh, Sh, S_max, p["u_max"], B),
                   runtime_seconds=dt_h, **common)
    res.record_run(method="threshold", variant="pure", u=ut, S=Stt, costs=ct, srs=srs_t,
                   feasibility=st.check_constraints(R_T, ut, Stt, S_max, p["u_max"], B),
                   runtime_seconds=dt_t, **common)
    res.record_run(method="dp", u=u_dp, S=Sdp, costs=dp["costs"], srs=dp["SRS_star"],
                   feasibility=st.check_constraints(R_T, u_dp, Sdp, S_max, p["u_max"], B),
                   runtime_seconds=dt_dp, **common)

    print(f"[{instance:7s} T={T:2d} L={L}] "
          f"hist={srs_h:.4e}  thr={srs_t:.4e} (dSRS {srs_t-srs_h:+.2e})  "
          f"dp={dp['SRS_star']:.4e} (dSRS {dp['SRS_star']-srs_h:+.2e})")


def main() -> None:
    df = fd.build_weekly_benchmark(write=False)
    constants = fc.load_official_constants()
    print(f"Dataset: {df.attrs['n_weeks']} semanas | S_min oficial = {constants['S_min_m3']:,.0f} m^3\n")
    for T, L in STANDARD_INSTANCES:
        run_instance(df, constants, T, L)
    print(f"\nResultados en {res.RESULTS_ROOT}/runs_summary.csv  y  {res.RESULTS_ROOT}/runs/")


if __name__ == "__main__":
    main()
