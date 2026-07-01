"""Correr QAOA de Falcon en hardware REAL de IBM Quantum (fixed-params, un job).

Patrón (reusado de Georgia `run_qaoa_hardware.py`, MODE A): los parámetros β/γ se
entrenan en el simulador (`falcon_qaoa.run_qaoa`) y aquí solo se **envía UN job** con
`SamplerV2(mode=backend)` a un backend real (sin Session; sirve en el plan Open). No hay
loop variacional en hardware (serían cientos de jobs).

Ordenamiento: variable i del QUBO == qubit i. Al medir con `measure_all`, el bitstring de
qiskit viene MSB-first → se **revierte** a LSB para que el bit i sea el qubit i (consistente
con el resto del pipeline). Solo se decodifican los primeros `vi.n` bits (decisión).

Requiere `qiskit-ibm-runtime` y una cuenta guardada válida (`QiskitRuntimeService()`).
Correr con `.venv-quantum/bin/python`.
"""
from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FuturesTimeout

import numpy as np

import falcon_qaoa as qa
import falcon_qubo as fq
import falcon_solvers as sv


def _measured_circuit(n, singles, pairs, betas, gammas):
    """QAOA con (β,γ) numéricos y measure_all() (para muestreo en hardware).

    Misma estructura que `falcon_qaoa.build_qaoa_circuit` pero SIN save_statevector.
    """
    from qiskit import QuantumCircuit

    qc = QuantumCircuit(n)
    qc.h(range(n))
    for beta, gamma in zip(betas, gammas):
        for q, h in singles:
            qc.rz(2.0 * float(gamma) * h, q)
        for q, r, J in pairs:
            qc.rzz(2.0 * float(gamma) * J, q, r)
        for q in range(n):
            qc.rx(2.0 * float(beta), q)
    qc.measure_all()
    return qc


def _measured_circuit_xy(n, singles, pairs, betas, gammas, groups):
    """QAOA con XY-mixer (subespacio one-hot) + (β,γ) numéricos + measure_all().

    Espeja `falcon_qaoa.build_qaoa_circuit_xy` (init W-state por semana; costo RZ/RZZ;
    mixer XY en anillo RXX+RYY) para que los params entrenados en sim sean consistentes.
    """
    from qiskit import QuantumCircuit

    qc = QuantumCircuit(n)
    for g in groups:                                   # init one-hot (estado W) por semana
        qc.initialize(qa._w_state_vector(len(g)), list(g))
    for beta, gamma in zip(betas, gammas):
        for q, h in singles:
            qc.rz(2.0 * float(gamma) * h, q)
        for q, r, J in pairs:
            qc.rzz(2.0 * float(gamma) * J, q, r)
        for g in groups:
            L = len(g)
            edges = [(g[i], g[(i + 1) % L]) for i in range(L)] if L > 2 else [(g[0], g[1])]
            for a, b in edges:
                qc.rxx(2.0 * float(beta), a, b)
                qc.ryy(2.0 * float(beta), a, b)
    qc.measure_all()
    return qc


def default_mitigation_options(num_randomizations=5):
    """Opciones de mitigación de error para SamplerV2: DD (XY4) + gate/measure twirling."""
    return {
        "dynamical_decoupling": {"enable": True, "sequence_type": "XY4"},
        "twirling": {"enable_gates": True, "enable_measure": True,
                     "num_randomizations": num_randomizations},
    }


def get_service():
    """Carga el servicio desde la cuenta guardada; error accionable si falla."""
    from qiskit_ibm_runtime import QiskitRuntimeService

    try:
        return QiskitRuntimeService()
    except Exception as e:  # noqa: BLE001
        raise RuntimeError(
            "No se pudo cargar QiskitRuntimeService desde la cuenta guardada. "
            "Guardá credenciales una vez con:\n"
            "  QiskitRuntimeService.save_account(channel='ibm_cloud', token='<API>', "
            "instance='<CRN>', set_as_default=True, overwrite=True)\n"
            f"Detalle: {type(e).__name__}: {e}") from e


def pick_backend(service, n_qubits: int):
    """Backend real menos ocupado con >= n_qubits qubits."""
    return service.least_busy(operational=True, simulator=False, min_num_qubits=n_qubits)


def decode_counts(counts, Q, const, vi, meta, levels, R_obs, dS, S0, *,
                  S_min, S_max, u_max, B, weights, half, topk=10):
    """Decodifica el histograma de hardware: mejor factible (menor energía QUBO)."""
    total = sum(counts.values())
    is_binary = hasattr(vi, "bits")
    best = None
    scored = []
    for bs, cnt in counts.items():
        s = str(bs).replace(" ", "")
        xfull = np.array([int(b) for b in s[::-1]], dtype=float)   # revertir MSB->LSB (bit i = qubit i)
        xdec = xfull[:vi.n]
        lv = vi.decode_levels(xdec)
        valid = all(0 <= int(l) < vi.L for l in lv) if is_binary else bool(vi.is_onehot(xdec))
        if not valid:
            continue
        dv = sv.decode_and_verify(lv, vi, levels, R_obs, dS, S0, S_min=S_min, S_max=S_max,
                                  u_max=u_max, B=B, weights=weights)
        e = fq.qubo_energy(Q, sv._x_from_levels(lv, vi, meta, half), const)
        cand = {"lv": lv, "prob": cnt / total, "count": int(cnt), "energy": e, **dv}
        scored.append(cand)
        if best is None or (cand["feasible"], -cand["energy"]) > (best["feasible"], -best["energy"]):
            best = cand
    scored.sort(key=lambda c: c["prob"], reverse=True)
    return {"decoded": best, "n_feasible_samples": len(scored),
            "top_samples": [{"lv": [int(x) for x in c["lv"]], "prob": c["prob"],
                             "SRS": c["SRS"], "feasible": c["feasible"]} for c in scored[:topk]]}


def submit_job(Q, meta, betas, gammas, *, service=None, backend=None, shots=4096,
               optimization_level=3, mixer="x", groups=None, mitigation=True,
               num_randomizations=5):
    """Transpila y ENVÍA el job fixed-params (no espera). Devuelve handle + metadata.

    `mixer="xy"` usa el circuito XY (requiere `groups`); `mitigation=True` activa DD (XY4)
    + gate/measure twirling en SamplerV2.options.
    """
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
    from qiskit_ibm_runtime import SamplerV2

    n = meta["n_qubits"]
    service = service or get_service()
    if backend is None:
        backend = pick_backend(service, n)
    elif isinstance(backend, str):
        backend = service.backend(backend)

    singles, pairs, _ = qa._ising_terms(Q)
    if mixer == "xy":
        if not groups:
            raise ValueError("mixer='xy' requiere groups (qubits one-hot por semana)")
        qc = _measured_circuit_xy(n, singles, pairs, betas, gammas, groups)
    else:
        qc = _measured_circuit(n, singles, pairs, betas, gammas)
    pm = generate_preset_pass_manager(backend=backend, optimization_level=optimization_level)
    tqc = pm.run(qc)

    sampler = SamplerV2(mode=backend)
    if mitigation:
        opts = default_mitigation_options(num_randomizations)
        sampler.options.dynamical_decoupling.enable = opts["dynamical_decoupling"]["enable"]
        sampler.options.dynamical_decoupling.sequence_type = opts["dynamical_decoupling"]["sequence_type"]
        sampler.options.twirling.enable_gates = opts["twirling"]["enable_gates"]
        sampler.options.twirling.enable_measure = opts["twirling"]["enable_measure"]
        sampler.options.twirling.num_randomizations = opts["twirling"]["num_randomizations"]

    job = sampler.run([(tqc,)], shots=shots)
    jid = job.job_id()
    print(f"[HW] job {jid} -> {backend.name} | n={n}q, mixer={mixer}, "
          f"depth transpilado={tqc.depth()}, 2q-gates={tqc.num_nonlocal_gates()}, shots={shots}, "
          f"mitig={'DD+twirl' if mitigation else 'none'}")
    print(f"[HW] retirar luego: QiskitRuntimeService().job('{jid}').result()")
    return {"job": job, "job_id": jid, "backend": backend.name, "n_qubits": n,
            "shots": shots, "transpiled_depth": int(tqc.depth()), "mixer": mixer,
            "t_submit": time.time()}


def await_and_decode(handle, Q, const, vi, meta, levels, R_obs, dS, S0, *,
                     S_min, S_max, u_max, B, weights, half, timeout=1200):
    """Espera el resultado del job (con timeout) y decodifica el histograma."""
    job = handle["job"]
    with ThreadPoolExecutor(max_workers=1) as pool:
        fut = pool.submit(job.result)
        try:
            result = fut.result(timeout=timeout if timeout and timeout > 0 else None)
        except FuturesTimeout:
            return {"status": "TIMEOUT", "job_id": handle["job_id"],
                    "backend": handle["backend"], "n_qubits": handle["n_qubits"],
                    "shots": handle["shots"]}

    data = result[0].data
    creg = next(iter(vars(data)))                       # nombre del registro clásico (robusto)
    counts = getattr(data, creg).get_counts()
    dec = decode_counts(counts, Q, const, vi, meta, levels, R_obs, dS, S0, S_min=S_min,
                        S_max=S_max, u_max=u_max, B=B, weights=weights, half=half)
    dec.update({"status": "OK", "job_id": handle["job_id"], "backend": handle["backend"],
                "shots": handle["shots"], "runtime_seconds": time.time() - handle["t_submit"],
                "transpiled_depth": handle["transpiled_depth"], "counts_distinct": len(counts)})
    return dec


def run_on_hardware(Q, const, meta, vi, levels, R_obs, dS, S0, betas, gammas, *,
                    S_min, S_max, u_max, B, weights, half,
                    backend=None, shots=4096, timeout=1200, optimization_level=1):
    """Conveniencia: enviar UN job y esperar/decodificar (submit_job + await_and_decode)."""
    handle = submit_job(Q, meta, betas, gammas, backend=backend, shots=shots,
                        optimization_level=optimization_level)
    return await_and_decode(handle, Q, const, vi, meta, levels, R_obs, dS, S0, S_min=S_min,
                            S_max=S_max, u_max=u_max, B=B, weights=weights, half=half, timeout=timeout)
