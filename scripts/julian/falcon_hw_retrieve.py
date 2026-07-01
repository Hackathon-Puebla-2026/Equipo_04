"""Recupera de IBM las corridas chunked-QAOA por job_id, decodifica y registra.

`plan_blocks` (seed 42) es determinista → reconstruye los MISMOS bloques (Q/meta/vi/boundaries)
que en el submit, así que cada job completo se decodifica por su id sin re-enviar. Post-selecciona
el mejor sample FACTIBLE por bloque, concatena el `u` global, evalúa SRS/factibilidad global y
`gap_vs_full` vs el DP exacto, y registra `method="qaoa_chunked_hardware"`.

Además vuelca `results/julian/hw_analysis_<instance>.json` con, por bloque: fracción one-hot válida en
HW real, nº de muestras factibles, SRS del mejor, y la distribución (SRS, count) de muestras factibles
-> insumo de las figuras (`falcon_hw_analysis_figures.py`).

Jobs aún no DONE se saltan (re-correr luego). Correr: .venv-quantum/bin/python scripts/julian/falcon_hw_retrieve.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
JULIAN = os.path.join(ROOT, "results", "julian")
MANIFEST = os.path.join(JULIAN, "hw_jobs_manifest.json")

import falcon_baselines as bl
import falcon_chunked_hardware as chw
import falcon_constants as fc
import falcon_data as fd
import falcon_qaoa_hardware as hw
import falcon_results as res
import falcon_srs as srs
import falcon_storage as st

# Instancias y sus job_ids sembrados (small ya corrió; medium/large vienen del manifest).
SEED_JOBS = {
    ("small", 12, 3, 4): ["d92jq4d958jc73brt6lg", "d92jq4nd07jc73dvdj8g", "d92jq4vd07jc73dvdj90"],
}


def _load_manifest():
    jobs = {}
    for (label, T, L, bs), ids in SEED_JOBS.items():
        jobs[(label, T, L, bs)] = list(ids)
    if os.path.exists(MANIFEST):
        for e in json.load(open(MANIFEST)):
            key = (res.instance_label(e["T"], e["L"]), e["T"], e["L"], e["block_size"])
            jobs.setdefault(key, [])
            # respetar orden de block_idx
            while len(jobs[key]) <= e["block_idx"]:
                jobs[key].append(None)
            jobs[key][e["block_idx"]] = e["job_id"]
    return jobs


def _feasible_distribution(counts, vi, levels, R_obs, dS, S0, *,
                           S_min, S_max, u_max, B, weights):
    """Fracción one-hot válida (en shots) y lista (SRS, count) de muestras factibles."""
    import falcon_solvers as sv
    total = sum(counts.values())
    n = vi.n
    onehot_shots = 0
    dist = []
    for bs, cnt in counts.items():
        x = np.array([int(ch) for ch in str(bs).replace(" ", "")[::-1]], dtype=float)
        xdec = x[:n]
        if not bool(vi.is_onehot(xdec)):
            continue
        onehot_shots += cnt
        lv = vi.decode_levels(xdec)
        dv = sv.decode_and_verify(lv, vi, levels, R_obs, dS, S0, S_min=S_min, S_max=S_max,
                                  u_max=u_max, B=B, weights=weights)
        if dv["feasible"]:
            dist.append((float(dv["SRS"]), int(cnt)))
    return onehot_shots / total, dist


def retrieve_instance(svc, label, T, L, block_size, job_ids):
    df = fd.build_weekly_benchmark(write=False)
    c = fc.load_official_constants()
    Rall = df["R_obs_m3_week"].to_numpy()
    dSall = df["DeltaS_obs_m3"].to_numpy()
    S0 = df.attrs["S0_m3"]
    Smin, Smax = c["S_min_m3"], c["S_max_m3"]

    plan = chw.plan_blocks(Rall, dSall, S0, T, S_min=Smin, S_max=Smax, L=L,
                           block_size=block_size, p=1, restarts=2, maxiter=60, eta=c["eta"])
    w = plan["weights"]
    if len(job_ids) != len(plan["blocks"]):
        print(f"[{label}] ⚠ {len(job_ids)} job_ids vs {len(plan['blocks'])} bloques; se procesan los disponibles")

    u_all, per_block, ok_backend = [], [], None
    for i, b in enumerate(plan["blocks"]):
        jid = job_ids[i] if i < len(job_ids) else None
        if jid is None:
            print(f"[{label}] bloque {i}: sin job_id -> saltado")
            return None
        job = svc.job(jid)
        stt = str(job.status())
        if stt not in ("DONE", "JobStatus.DONE"):
            print(f"[{label}] bloque {i} job {jid}: {stt} (aún no listo) -> instancia incompleta")
            return None
        ok_backend = job.backend().name if callable(getattr(job, "backend", None)) else "?"
        r = job.result()
        data = r[0].data
        creg = next(iter(vars(data)))
        counts = getattr(data, creg).get_counts()
        dec = hw.decode_counts(counts, b["Q"], b["const"], b["vi"], b["meta"], plan["levels"],
                               b["Rb"], b["dSb"], b["S0_blk"], S_min=Smin, S_max=Smax,
                               u_max=plan["u_max"], B=b["B_blk"], weights=w, half=plan["half"])
        oh_frac, dist = _feasible_distribution(counts, b["vi"], plan["levels"], b["Rb"], b["dSb"],
                                               b["S0_blk"], S_min=Smin, S_max=Smax,
                                               u_max=plan["u_max"], B=b["B_blk"], weights=w)
        d = dec["decoded"]
        if d is None:
            print(f"[{label}] bloque {i}: 0 factibles en HW -> uso u del SIM (marcado)")
            u_blk = [plan["levels"][l] for l in b["sim_lv"]]
            src = "sim_fallback"
        else:
            u_blk = [float(x) for x in d["u"]]
            src = "hw"
        u_all.extend(u_blk)
        per_block.append({"block_idx": i, "job_id": jid, "n_qubits": b["meta"]["n_qubits"],
                          "one_hot_frac_hw": oh_frac, "n_feasible_samples": len(dist),
                          "best_SRS": (d["SRS"] if d else None), "best_feasible": (d["feasible"] if d else False),
                          "source": src, "feasible_srs_dist": dist})
        print(f"[{label}] bloque {i} ({b['meta']['n_qubits']}q): one-hot={oh_frac:.3f} "
              f"feasibles={len(dist)} best_SRS={d['SRS']:.4f}" if d else
              f"[{label}] bloque {i}: sin factible")

    # --- eval GLOBAL ---
    u = np.array(u_all, dtype=float)
    R = plan["R"]
    S = st.simulate_storage(S0, plan["dS"], u)
    costs = srs.compute_costs(S, u, Smin)
    SRS = srs.compute_srs(costs, w)
    chk = st.check_constraints(R, u, S, Smax, plan["u_max"], plan["B_global"])
    dp = bl.dp_optimal(R, plan["dS"], S0, S_min=Smin, S_max=Smax, delta_u=plan["delta_u"],
                       L=L, weights=w, B=plan["B_global"])
    uh = bl.historical(T)
    srs_h = srs.compute_srs(srs.compute_costs(st.simulate_storage(S0, plan["dS"], uh), uh, Smin), w)
    ar = SRS / dp["SRS_star"] if dp["SRS_star"] else None
    nq = max(b["meta"]["n_qubits"] for b in plan["blocks"])

    solver = {"encoding": "onehot", "n_qubits": nq, "p_depth": 1, "energy": None,
              "approximation_ratio": ar, "seed": 42, "simulator": f"ibmq:{ok_backend}",
              "penalties": None, "mixer": "xy", "mitigation": "DD+twirling",
              "job_ids": job_ids, "one_hot_frac": [round(b["one_hot_frac_hw"], 4) for b in per_block]}
    res.record_run(method="qaoa_chunked_hardware", instance=label, T=T, L=L, params=plan["params"],
                   weights=w, constants=c, B=plan["B_global"], u=u, S=S, costs=costs, srs=SRS,
                   variant=f"blk{block_size}_xy_p1",
                   feasibility={"feasible": chk["feasible"], "violations": chk["violations"]},
                   runtime_seconds=0.0, references={"historical": srs_h, "dp": dp["SRS_star"]},
                   solver=solver)

    out = {"instance": label, "T": T, "L": L, "block_size": block_size, "backend": ok_backend,
           "SRS": SRS, "dp": dp["SRS_star"], "historical": srs_h, "AR": ar,
           "feasible": bool(chk["feasible"]), "gap_vs_full": dp["SRS_star"] - SRS,
           "n_qubits": nq, "u": u.tolist(), "per_block": per_block}
    json.dump(out, open(os.path.join(JULIAN, f"hw_analysis_{label}.json"), "w"), indent=2)
    print(f"[{label}] GLOBAL HW-chunked SRS={SRS:.6f} dp={dp['SRS_star']:.6f} hist={srs_h:.6f} "
          f"AR={ar:.4f} feasible={chk['feasible']} gap={dp['SRS_star']-SRS:.4f}\n")
    return out


def main():
    from qiskit_ibm_runtime import QiskitRuntimeService
    svc = QiskitRuntimeService()
    jobs = _load_manifest()
    results = []
    for (label, T, L, bs), ids in sorted(jobs.items(), key=lambda kv: kv[0][1]):
        if not ids or any(j is None for j in ids):
            print(f"[{label}] manifest incompleto ({ids}) -> saltado")
            continue
        print(f"=== recuperando {label} T{T}/L{L} ({len(ids)} bloques) ===")
        r = retrieve_instance(svc, label, T, L, bs, ids)
        if r:
            results.append(r)
    print("=== resumen ===")
    for r in results:
        print(f"{r['instance']:7} SRS={r['SRS']:.6f} dp={r['dp']:.6f} feasible={r['feasible']} "
              f"AR={r['AR']:.4f} gap={r['gap_vs_full']:.4f}")


if __name__ == "__main__":
    main()
