"""Corre el solver por etapas (chunking) y registra + descompone el gap.

Objetivo: T26/L5 con bloques de 5 semanas = 25 qubits/bloque (entra en statevector),
resuelto con QAOA real (`falcon_chunked.staged_solve`, solver="qaoa"). Descompone el gap:

    full DP (optimo global)  >=  DP-chunked (perdida por trocear)  >=  QAOA-chunked (perdida QAOA)

Registra el QAOA-chunked con method="qaoa_chunked". Correr con `.venv-quantum/bin/python`.
"""
from __future__ import annotations

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import falcon_baselines as bl
import falcon_chunked as ch
import falcon_constants as fc
import falcon_data as fd
import falcon_results as res
import falcon_srs as srs
import falcon_storage as st


def run(T=26, L=5, block_size=5, *, p=1, restarts=2, maxiter=60):
    df = fd.build_weekly_benchmark(write=False)
    c = fc.load_official_constants()
    Rall = df["R_obs_m3_week"].to_numpy()
    dSall = df["DeltaS_obs_m3"].to_numpy()
    S0 = df.attrs["S0_m3"]
    Smin, Smax = c["S_min_m3"], c["S_max_m3"]
    R = Rall[:T]

    # QAOA-chunked (real)
    t0 = time.perf_counter()
    qc = ch.staged_solve(Rall, dSall, S0, T, S_min=Smin, S_max=Smax, L=L, block_size=block_size,
                         eta=c["eta"], solver="qaoa", link_smooth=True,
                         qaoa_kwargs={"p": p, "restarts": restarts, "maxiter": maxiter})
    rt = time.perf_counter() - t0
    w = qc["weights"]

    # DP-chunked (oraculo del troceo) y full DP (oraculo global)
    dpc = ch.staged_solve(Rall, dSall, S0, T, S_min=Smin, S_max=Smax, L=L, block_size=block_size,
                          eta=c["eta"], solver="dp", link_smooth=True)
    g = ch.gap_vs_full(qc, Rall, dSall, S0, S_min=Smin, S_max=Smax, L=L)
    dp_full = g["dp_full"]["SRS_star"]

    # baselines de referencia
    uh = bl.historical(T)
    srs_h = srs.compute_srs(srs.compute_costs(st.simulate_storage(S0, dSall[:T], uh), uh, Smin), w)
    ut = bl.threshold_rule(R, dSall[:T], S0, Smin, qc["params"]["delta_u"])
    srs_t = srs.compute_srs(srs.compute_costs(st.simulate_storage(S0, dSall[:T], ut), ut, Smin), w)
    refs = {"historical": srs_h, "threshold": srs_t, "dp": dp_full}

    nq = max(b["n_qubits"] for b in qc["blocks"] if b["n_qubits"])
    ar = qc["SRS"] / dp_full if dp_full != 0 else None
    solver = {"encoding": "onehot", "n_qubits": nq, "p_depth": p,
              "energy": None, "approximation_ratio": ar, "seed": 42,
              "simulator": "qiskit_aer_statevector", "penalties": None,
              "optimizer_iterations": None}
    res.record_run(method="qaoa_chunked", instance=res.instance_label(T, L), T=T, L=L,
                   params=qc["params"], weights=w, constants=c, B=qc["B"], u=qc["u"], S=qc["S"],
                   costs=qc["costs"], srs=qc["SRS"], variant=f"blk{block_size}_p{p}",
                   feasibility={"feasible": qc["feasible"], "violations": qc["violations"]},
                   runtime_seconds=rt, references=refs, solver=solver)

    # --- reporte ---
    print(f"\n=== QAOA por etapas (chunking) T{T}/L{L}, bloques de {block_size} ({nq} qubits/bloque) ===")
    print(f"{'metodo':22} {'SRS':>11} {'gap_vs_full':>12} {'feasible':>8}")
    print("-" * 56)
    print(f"{'full DP (optimo)':22} {dp_full:>11.6f} {0.0:>12.2e} {'True':>8}")
    print(f"{'DP-chunked':22} {dpc['SRS']:>11.6f} {dp_full-dpc['SRS']:>12.2e} {str(dpc['feasible']):>8}")
    print(f"{'QAOA-chunked':22} {qc['SRS']:>11.6f} {g['gap_vs_full']:>12.2e} {str(qc['feasible']):>8}")
    print(f"{'historico u=0':22} {srs_h:>11.6f}")
    print(f"{'umbral (pura)':22} {srs_t:>11.6f}")
    print(f"\nQAOA real: {qc['n_qaoa_blocks']}/{qc['n_blocks']} bloques factibles por QAOA (sin fallback silencioso)")
    print(f"runtime QAOA-chunked: {rt:.1f}s")
    print("\nTabla por bloque (QAOA):")
    for b in qc["blocks"]:
        print(f"  bloque {b['start']+1:>2}-{b['end']:<2}: {b['method']:5} feasible={str(b['feasible']):5} "
              f"nq={b['n_qubits']} block_SRS={b['block_SRS']:.4f} niveles={b['selected_levels']}")
    print(f"\nRegistrado en {res.RESULTS_ROOT}/runs_summary.csv (method=qaoa_chunked)")
    return qc, dpc, dp_full


if __name__ == "__main__":
    run()
