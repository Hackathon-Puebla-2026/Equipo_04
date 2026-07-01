"""Envía a hardware IBM los bloques de medium/large (chunked-QAOA XY) y escribe el manifest.

DESACOPLADO del retrieve: solo entrena en sim (determinista, seed 42), envía cada bloque y anota
`{instance,T,L,block_size,block_idx,job_id,backend,shots,depth,betas,gammas}` en
`results/julian/hw_jobs_manifest.json`. NO espera resultados → `falcon_hw_retrieve.py` los procesa
cuando estén DONE (resiliente a la cola; el proceso puede morir sin perder nada).

Training barato (restarts=1, maxiter=40) porque el statevector XY a 25 qubits es caro (fue lo que colgó
el runner monolítico). Correr en background: .venv-quantum/bin/python scripts/julian/falcon_hw_submit.py
"""
from __future__ import annotations

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
JULIAN = os.path.join(ROOT, "results", "julian")
MANIFEST = os.path.join(JULIAN, "hw_jobs_manifest.json")

import falcon_chunked_hardware as chw
import falcon_constants as fc
import falcon_data as fd
import falcon_qaoa_hardware as hw

INSTANCES = [
    dict(T=26, L=5, block_size=5),
    dict(T=52, L=5, block_size=5),
]


def _append_manifest(entry):
    data = json.load(open(MANIFEST)) if os.path.exists(MANIFEST) else []
    data.append(entry)
    os.makedirs(JULIAN, exist_ok=True)
    json.dump(data, open(MANIFEST, "w"), indent=2)


def submit_instance(svc, backend, T, L, block_size, *, shots=4096, restarts=1, maxiter=40):
    df = fd.build_weekly_benchmark(write=False)
    c = fc.load_official_constants()
    Rall = df["R_obs_m3_week"].to_numpy()
    dSall = df["DeltaS_obs_m3"].to_numpy()
    S0 = df.attrs["S0_m3"]
    print(f"[{T}/{L}] entrenando en sim (XY, restarts={restarts}, maxiter={maxiter})...")
    plan = chw.plan_blocks(Rall, dSall, S0, T, S_min=c["S_min_m3"], S_max=c["S_max_m3"], L=L,
                           block_size=block_size, p=1, restarts=restarts, maxiter=maxiter,
                           eta=c["eta"])
    for i, b in enumerate(plan["blocks"]):
        h = hw.submit_job(b["Q"], b["meta"], b["betas"], b["gammas"], service=svc, backend=backend,
                          shots=shots, optimization_level=3, mixer="xy", groups=b["groups"],
                          mitigation=True, num_randomizations=5)
        _append_manifest({"T": T, "L": L, "block_size": block_size, "block_idx": i,
                          "job_id": h["job_id"], "backend": h["backend"], "shots": shots,
                          "depth": h["transpiled_depth"],
                          "betas": [float(x) for x in b["betas"]],
                          "gammas": [float(x) for x in b["gammas"]]})
    print(f"[{T}/{L}] {len(plan['blocks'])} bloques enviados y anotados en el manifest.\n")


def main():
    svc = hw.get_service()
    backend = hw.pick_backend(svc, 25)
    print(f"[HW] backend: {backend.name} ({backend.num_qubits}q)")
    for ins in INSTANCES:
        submit_instance(svc, backend, **ins)
    print(f"Manifest: {MANIFEST}")
    print("Recuperar cuando estén DONE: .venv-quantum/bin/python scripts/julian/falcon_hw_retrieve.py")


if __name__ == "__main__":
    main()
