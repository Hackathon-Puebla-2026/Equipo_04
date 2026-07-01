"""Entrena β/γ en simulador y corre las instancias en HARDWARE real IBM Quantum.

Instancias: T5/L3 (one-hot, 17q) y T12/L3 (binary, 24q). Fixed-params: se entrena en el
statevector (`falcon_qaoa.run_qaoa`) y se envía UN job por instancia (auto least-busy).
Se envían AMBOS jobs antes de esperar (para encolar en paralelo). Registra
`method="qaoa_hardware"`. Correr con `.venv-quantum/bin/python`.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import falcon_baselines as bl
import falcon_config as fcfg
import falcon_constants as fc
import falcon_data as fd
import falcon_qaoa as qa
import falcon_qaoa_hardware as hw
import falcon_qubo as fq
import falcon_results as res
import falcon_srs as srs
import falcon_storage as st


def _setup(T, L, encoding, balance):
    df = fd.build_weekly_benchmark(write=False)
    c = fc.load_official_constants()
    Rall = df["R_obs_m3_week"].to_numpy()
    R, dS = Rall[:T], df["DeltaS_obs_m3"].to_numpy()[:T]
    S0 = df.attrs["S0_m3"]
    pr = fc.instance_params(Rall, T, L)
    w = fc.compute_weights(T, c["S_min_m3"], pr["u_max"])
    B = c["eta"] * float(R.sum())
    half = (L - 1) // 2
    cfg = fcfg.FalconConfig(T=T, L=L, encoding=encoding, balance=balance)
    Q, const, vi, meta = fq.build_qubo(cfg, R, dS, S0, S_min=c["S_min_m3"],
                                       delta_u=pr["delta_u"], levels=pr["levels"], weights=w, B=B)
    dp = bl.dp_optimal(R, dS, S0, S_min=c["S_min_m3"], S_max=c["S_max_m3"],
                       delta_u=pr["delta_u"], L=L, weights=w, B=B)
    uh = bl.historical(T)
    srs_h = srs.compute_srs(srs.compute_costs(st.simulate_storage(S0, dS, uh), uh, c["S_min_m3"]), w)
    ut = bl.threshold_rule(R, dS, S0, c["S_min_m3"], pr["delta_u"])
    srs_t = srs.compute_srs(srs.compute_costs(st.simulate_storage(S0, dS, ut), ut, c["S_min_m3"]), w)
    return dict(T=T, L=L, enc=encoding, R=R, dS=dS, S0=S0, c=c, pr=pr, w=w, B=B, half=half,
                Q=Q, const=const, vi=vi, meta=meta,
                refs={"historical": srs_h, "threshold": srs_t, "dp": dp["SRS_star"]},
                dp=dp["SRS_star"])


def main(shots=4096, p=1):
    instances = [
        dict(T=5, L=3, encoding="onehot", balance="slack"),
        dict(T=12, L=3, encoding="binary", balance="soft"),
    ]
    service = hw.get_service()
    jobs = []
    for ins in instances:
        s = _setup(**ins)
        print(f"\n[{res.instance_label(s['T'],s['L'])} T{s['T']}/L{s['L']} {s['enc']}] "
              f"entrenando β/γ en simulador (n={s['meta']['n_qubits']} qubits)...")
        out = qa.run_qaoa(s["Q"], s["const"], s["meta"], s["vi"], s["pr"]["levels"], s["R"], s["dS"],
                          s["S0"], S_min=s["c"]["S_min_m3"], S_max=s["c"]["S_max_m3"],
                          u_max=s["pr"]["u_max"], B=s["B"], weights=s["w"], half=s["half"],
                          p=p, restarts=3, seed=42, maxiter=100)
        s["betas"], s["gammas"] = out["betas"], out["gammas"]
        print(f"  β={[round(b,3) for b in s['betas']]} γ={[round(g,3) for g in s['gammas']]} "
              f"(SRS sim={out['decoded']['SRS']:.6f})")
        s["handle"] = hw.submit_job(s["Q"], s["meta"], s["betas"], s["gammas"],
                                    service=service, shots=shots, optimization_level=1)
        jobs.append(s)

    print("\n=== esperando resultados de hardware (timeout 1200s c/u) ===")
    for s in jobs:
        dec = hw.await_and_decode(s["handle"], s["Q"], s["const"], s["vi"], s["meta"],
                                  s["pr"]["levels"], s["R"], s["dS"], s["S0"],
                                  S_min=s["c"]["S_min_m3"], S_max=s["c"]["S_max_m3"],
                                  u_max=s["pr"]["u_max"], B=s["B"], weights=s["w"], half=s["half"])
        lbl = f"{res.instance_label(s['T'],s['L'])} T{s['T']}/L{s['L']} {s['enc']}"
        if dec["status"] != "OK":
            print(f"[{lbl}] {dec['status']} job_id={dec['job_id']} backend={dec['backend']} "
                  f"(retirar luego con QiskitRuntimeService().job('{dec['job_id']}').result())")
            continue
        d = dec["decoded"]
        if d is None:
            print(f"[{lbl}] hardware: 0 muestras factibles de {dec['counts_distinct']} bitstrings "
                  f"distintos (ruido); job_id={dec['job_id']} backend={dec['backend']}")
            continue
        S = st.simulate_storage(s["S0"], s["dS"], d["u"])
        ar = d["SRS"] / s["dp"] if s["dp"] else None
        solver = {"encoding": s["meta"]["encoding"], "n_qubits": s["meta"]["n_qubits"], "p_depth": p,
                  "energy": float(d["energy"]), "approximation_ratio": ar, "seed": 42,
                  "simulator": f"ibmq:{dec['backend']}", "penalties": s["meta"]["penalties"],
                  "beta": s["betas"], "gamma": s["gammas"],
                  "job_id": dec["job_id"], "shots": dec["shots"],
                  "transpiled_depth": dec["transpiled_depth"],
                  "n_feasible_samples": dec["n_feasible_samples"]}
        res.record_run(method="qaoa_hardware", instance=res.instance_label(s["T"], s["L"]),
                       T=s["T"], L=s["L"], params=s["pr"], weights=s["w"], constants=s["c"],
                       B=s["B"], u=d["u"], S=S, costs=d["costs"], srs=d["SRS"],
                       variant=f"{s['enc']}_p{p}",
                       feasibility={"feasible": d["feasible"], "violations": d["violations"]},
                       runtime_seconds=dec["runtime_seconds"], references=s["refs"], solver=solver)
        print(f"[{lbl}] HW SRS={d['SRS']:.6f} dp={s['dp']:.6f} AR={ar:.4f} feasible={d['feasible']} "
              f"prob={d['prob']:.3f} | {dec['backend']} job={dec['job_id']} "
              f"({dec['n_feasible_samples']} muestras factibles, {dec['runtime_seconds']:.0f}s)")
    print(f"\nRegistrado en {res.RESULTS_ROOT}/runs_summary.csv (method=qaoa_hardware)")


if __name__ == "__main__":
    main()
