"""Corre QAOA (Qiskit Aer statevector) sobre el QUBO y lo registra + compara con baselines.

- debug (T5/L3): one-hot = 15 qubits (+ slack de balance) -> entra en statevector.
- small (T12/L3): requiere encoding compacto (binary, Fase 2) para entrar en statevector
  (one-hot serian 36 qubits). Se corre solo si el encoding compacto esta disponible.

Registra en results/ con method="qaoa" y los campos cuanticos del schema
(encoding, n_qubits, p_depth, beta, gamma, optimizer_iterations, approximation_ratio).

Uso: python scripts/falcon_run_qaoa.py    (con el venv cuantico: .venv-quantum/bin/python)
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
import falcon_qaoa as qa
import falcon_qubo as fq
import falcon_results as res
import falcon_srs as srs
import falcon_storage as st


def run(T, L, *, p=1, restarts=5, seed=42, maxiter=200, encoding="onehot", balance="slack"):
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

    # referencias exactas (dp = optimo global) + baselines
    dp = bl.dp_optimal(R, dS, S0, S_min=c["S_min_m3"], S_max=c["S_max_m3"],
                       delta_u=pr["delta_u"], L=L, weights=w, B=B)
    uh = bl.historical(T)
    srs_h = srs.compute_srs(srs.compute_costs(st.simulate_storage(S0, dS, uh), uh, c["S_min_m3"]), w)
    ut = bl.threshold_rule(R, dS, S0, c["S_min_m3"], pr["delta_u"])
    srs_t = srs.compute_srs(srs.compute_costs(st.simulate_storage(S0, dS, ut), ut, c["S_min_m3"]), w)
    refs = {"historical": srs_h, "threshold": srs_t, "dp": dp["SRS_star"]}

    t0 = time.perf_counter()
    out = qa.run_qaoa(Q, const, meta, vi, pr["levels"], R, dS, S0, S_min=c["S_min_m3"],
                      S_max=c["S_max_m3"], u_max=pr["u_max"], B=B, weights=w, half=half,
                      p=p, restarts=restarts, seed=seed, maxiter=maxiter)
    rt = time.perf_counter() - t0

    d = out["decoded"]
    S = st.simulate_storage(S0, dS, d["u"])
    ar = d["SRS"] / dp["SRS_star"] if dp["SRS_star"] != 0 else None
    solver = {"encoding": meta["encoding"], "n_qubits": meta["n_qubits"], "p_depth": p,
              "energy": float(d["energy"]), "approximation_ratio": ar, "seed": seed,
              "simulator": "qiskit_aer_statevector", "penalties": meta["penalties"],
              "beta": out["betas"], "gamma": out["gammas"],
              "optimizer_iterations": out["n_iter"]}

    variant = f"{encoding}_p{p}"
    res.record_run(method="qaoa", instance=res.instance_label(T, L), T=T, L=L,
                   params=pr, weights=w, constants=c, B=B, u=d["u"], S=S, costs=d["costs"],
                   srs=d["SRS"], variant=variant,
                   feasibility={"feasible": d["feasible"], "violations": d["violations"]},
                   runtime_seconds=rt, references=refs, solver=solver)
    print(f"[{res.instance_label(T,L):5s} T={T} L={L}] qaoa {variant}: "
          f"SRS={d['SRS']:.6f}  dp={dp['SRS_star']:.6f}  hist={srs_h:.6f}  thr={srs_t:.6f}  "
          f"AR={ar:.4f}  feasible={d['feasible']}  n_qubits={meta['n_qubits']}  "
          f"prob={d['prob']:.2e}  iters={out['n_iter']}  ({rt:.1f}s)")
    return {"instance": res.instance_label(T, L), "T": T, "L": L, "encoding": encoding,
            "SRS_qaoa": d["SRS"], "dp": dp["SRS_star"], "hist": srs_h, "thr": srs_t,
            "dSRS_vs_hist": d["SRS"] - srs_h, "AR": ar, "feasible": d["feasible"],
            "n_qubits": meta["n_qubits"], "p": p, "runtime": rt}


def _print_table(rows):
    print("\n=== Benchmark QAOA vs baselines (ΔSRS = SRS_qaoa - SRS_baseline) ===")
    hdr = f"{'instance':7} {'enc':6} {'p':>2} {'nq':>3} {'SRS_qaoa':>11} {'DP*':>11} " \
          f"{'hist':>11} {'thr':>11} {'ΔvsHist':>9} {'AR':>7} {'feas':>5} {'s':>6}"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(f"{r['instance']:7} {r['encoding'][:6]:6} {r['p']:>2} {r['n_qubits']:>3} "
              f"{r['SRS_qaoa']:>11.6f} {r['dp']:>11.6f} {r['hist']:>11.6f} {r['thr']:>11.6f} "
              f"{r['dSRS_vs_hist']:>9.5f} {r['AR']:>7.4f} {str(r['feasible']):>5} {r['runtime']:>6.1f}")


def main():
    rows = []
    # debug T5/L3 (one-hot, 17 qubits): p=1 y p=2
    rows.append(run(5, 3, p=1, encoding="onehot"))
    rows.append(run(5, 3, p=2, encoding="onehot"))
    # small T12/L3 (binary + balance soft -> 24 qubits, entra en statevector): p=1
    rows.append(run(12, 3, p=1, encoding="binary", balance="soft", restarts=3, maxiter=100))
    _print_table(rows)
    print(f"\nRegistrado en {res.RESULTS_ROOT}/runs_summary.csv")


if __name__ == "__main__":
    main()
