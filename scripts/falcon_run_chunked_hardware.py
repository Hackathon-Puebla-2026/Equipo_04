"""QAOA chunked en HARDWARE real IBM, con todos los trucos, sobre varias instancias.

Por instancia: entrena la chunked QAOA (XY-mixer) en simulador (`plan_blocks`), envía UN job
por bloque a hardware least-busy (XY measured + DD + twirling + opt_level=3), decodifica con
post-selección del mejor sample FACTIBLE, concatena y evalúa SRS/factibilidad GLOBAL + gap vs DP.

Mejora sobre la corrida de ivan (infactible, sin mitigación, params heurísticos, depth ~2000):
params optimizados en sim, XY-mixer (subespacio factible), mitigación, bloques chicos (depth<<).

Correr con `.venv-quantum/bin/python` (largo → conviene background).
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

import falcon_baselines as bl
import falcon_chunked_hardware as chw
import falcon_constants as fc
import falcon_data as fd
import falcon_qaoa_hardware as hw
import falcon_results as res
import falcon_srs as srs
import falcon_storage as st

INSTANCES = [
    dict(T=12, L=3, block_size=4),
    dict(T=26, L=5, block_size=5),
    dict(T=52, L=5, block_size=5),
]


def main(shots=4096, p=1):
    df = fd.build_weekly_benchmark(write=False)
    c = fc.load_official_constants()
    Rall = df["R_obs_m3_week"].to_numpy()
    dSall = df["DeltaS_obs_m3"].to_numpy()
    S0 = df.attrs["S0_m3"]
    Smin, Smax = c["S_min_m3"], c["S_max_m3"]

    service = hw.get_service()
    backend = hw.pick_backend(service, 25)              # least-busy que cubra hasta 25q
    print(f"[HW] backend elegido: {backend.name} ({backend.num_qubits}q)")

    # --- 1) plan en sim + submit de TODOS los bloques (paralelo) ---
    plans = []
    for ins in INSTANCES:
        lbl = res.instance_label(ins["T"], ins["L"])
        print(f"\n[{lbl} T{ins['T']}/L{ins['L']} blk{ins['block_size']}] entrenando en sim (XY)...")
        plan = chw.plan_blocks(Rall, dSall, S0, ins["T"], S_min=Smin, S_max=Smax, L=ins["L"],
                               block_size=ins["block_size"], p=p, restarts=2, maxiter=60,
                               eta=c["eta"])
        for b in plan["blocks"]:
            b["handle"] = hw.submit_job(b["Q"], b["meta"], b["betas"], b["gammas"],
                                        service=service, backend=backend, shots=shots,
                                        optimization_level=3, mixer="xy", groups=b["groups"],
                                        mitigation=True, num_randomizations=5)
        plan["instance"] = ins
        plans.append(plan)

    # --- 2) await + decode (post-selección) + eval global por instancia ---
    print("\n=== esperando resultados de hardware (timeout 1800s/bloque) ===")
    rows = []
    for plan in plans:
        ins = plan["instance"]
        T, L = ins["T"], ins["L"]
        lbl = res.instance_label(T, L)
        w = plan["weights"]
        u_all, job_ids, depths, per_block = [], [], [], []
        for b in plan["blocks"]:
            dec = hw.await_and_decode(b["handle"], b["Q"], b["const"], b["vi"], b["meta"],
                                      plan["levels"], b["Rb"], b["dSb"], b["S0_blk"],
                                      S_min=Smin, S_max=Smax, u_max=plan["u_max"], B=b["B_blk"],
                                      weights=w, half=plan["half"], timeout=1800)
            job_ids.append(dec.get("job_id"))
            if dec["status"] != "OK" or dec.get("decoded") is None:
                print(f"[{lbl}] bloque {b['start']+1}-{b['end']}: {dec['status']} "
                      f"job={dec.get('job_id')} (usando u del SIM para continuar - marcado)")
                u_all.extend([plan["levels"][l] for l in b["sim_lv"]])   # continuidad (no silencioso)
                per_block.append({"block": f"{b['start']+1}-{b['end']}", "source": "sim_fallback",
                                  "job_id": dec.get("job_id")})
                continue
            d = dec["decoded"]
            u_all.extend([float(x) for x in d["u"]])
            depths.append(dec["transpiled_depth"])
            per_block.append({"block": f"{b['start']+1}-{b['end']}", "source": "hw",
                              "feasible": d["feasible"], "prob": d["prob"],
                              "n_feasible_samples": dec["n_feasible_samples"],
                              "job_id": dec["job_id"], "depth": dec["transpiled_depth"]})

        # --- eval GLOBAL sobre el u del hardware ---
        u = np.array(u_all, dtype=float)
        R = plan["R"]
        S = st.simulate_storage(S0, plan["dS"], u)
        costs = srs.compute_costs(S, u, Smin)
        SRS = srs.compute_srs(costs, w)
        chk = st.check_constraints(R, u, S, Smax, plan["u_max"], plan["B_global"])

        # referencias
        dp = bl.dp_optimal(R, plan["dS"], S0, S_min=Smin, S_max=Smax, delta_u=plan["delta_u"],
                           L=L, weights=w, B=plan["B_global"])
        uh = bl.historical(T)
        srs_h = srs.compute_srs(srs.compute_costs(st.simulate_storage(S0, plan["dS"], uh), uh, Smin), w)
        ar = SRS / dp["SRS_star"] if dp["SRS_star"] else None
        nq = max(b["meta"]["n_qubits"] for b in plan["blocks"])
        solver = {"encoding": "onehot", "n_qubits": nq, "p_depth": p, "energy": None,
                  "approximation_ratio": ar, "seed": 42, "simulator": f"ibmq:{backend.name}",
                  "penalties": None, "beta": plan["blocks"][0]["betas"],
                  "gamma": plan["blocks"][0]["gammas"], "job_ids": job_ids,
                  "shots": shots, "mixer": "xy", "mitigation": "DD+twirling",
                  "transpiled_depth_max": max(depths) if depths else None, "per_block": per_block}
        res.record_run(method="qaoa_chunked_hardware", instance=lbl, T=T, L=L, params=plan["params"],
                       weights=w, constants=c, B=plan["B_global"], u=u, S=S, costs=costs, srs=SRS,
                       variant=f"blk{ins['block_size']}_xy_p{p}",
                       feasibility={"feasible": chk["feasible"], "violations": chk["violations"]},
                       runtime_seconds=0.0, references={"historical": srs_h, "dp": dp["SRS_star"]},
                       solver=solver)
        rows.append(dict(lbl=lbl, T=T, L=L, SRS=SRS, dp=dp["SRS_star"], hist=srs_h, ar=ar,
                         feas=chk["feasible"], nq=nq, depth=(max(depths) if depths else None)))
        print(f"[{lbl}] HW-chunked SRS={SRS:.6f} dp={dp['SRS_star']:.6f} hist={srs_h:.6f} "
              f"AR={ar:.4f} feasible={chk['feasible']} nq={nq} depth_max={max(depths) if depths else '-'}")

    print("\n=== Resumen QAOA chunked en hardware (todos los trucos) ===")
    print(f"{'inst':7} {'nq':>3} {'SRS_hw':>11} {'DP*':>11} {'hist':>11} {'AR':>7} {'feas':>5} {'depth':>6}")
    for r in rows:
        print(f"{r['lbl']:7} {r['nq']:>3} {r['SRS']:>11.6f} {r['dp']:>11.6f} {r['hist']:>11.6f} "
              f"{r['ar']:>7.4f} {str(r['feas']):>5} {str(r['depth']):>6}")
    print(f"\nRegistrado en {res.RESULTS_ROOT}/runs_summary.csv (method=qaoa_chunked_hardware)")
    print("Comparar vs ivan (T26/L5 IBM): SRS=-0.3322, INFACTIBLE, depth~2000, sin mitigacion.")


if __name__ == "__main__":
    main()
