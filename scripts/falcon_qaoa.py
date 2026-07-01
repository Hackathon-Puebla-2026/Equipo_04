"""QAOA sobre el QUBO de Falcon (simulador statevector, Qiskit + Aer).

Flujo: `Q,const` (de `falcon_qubo.build_qubo`) -> `to_quadratic_program` -> Ising
(`qp.to_ising()`, delega a qiskit el algebra x->z) -> circuito QAOA manual
(H inicial; por capa: capa de costo RZ/RZZ + mixer RX) -> statevector en Aer.

Convenciones (reuso de Georgia, `docs/georgia_qubo_snippets.md`): seed=42, restart
r usa seed+r, >=5 restarts si Q mal condicionada, se trackea `const` aparte.

Robustez de ordenamiento: la variable i del QUBO == qubit i == bit i del indice del
statevector (little-endian en qiskit). Por eso la energia esperada y el muestreo se
calculan desde el ARRAY del statevector y una diagonal QUBO precomputada, indexados
igual por construccion; nunca se leen labels de bitstrings de qiskit (evita reversos).

Solo apto para instancias chicas (statevector: n<=~24 qubits). Instancias grandes
requieren MPS/sampling (fuera de este modulo).
"""
from __future__ import annotations

import numpy as np

import falcon_qubo as fq
import falcon_solvers as sv

MAX_QUBITS_STATEVECTOR = 26   # tope de seguridad (2^26 * 16B ~ 1 GB)


def precompute_diagonal(Q: np.ndarray, const: float, chunk: int = 1 << 20) -> np.ndarray:
    """Energia QUBO H(x)=xᵀQx+const para los 2ⁿ estados base (indexados como el statevector).

    x[i] = bit i del indice (little-endian) == qubit i == variable i del QUBO.
    Se procesa por lotes (`chunk`) para acotar memoria (la matriz de bits (N,n) seria
    enorme en n=24: ~3 GB). Devuelve un array de 2ⁿ (134 MB en n=24).
    """
    n = Q.shape[0]
    N = 1 << n
    ar = np.arange(n)
    out = np.empty(N, dtype=np.float64)
    for start in range(0, N, chunk):
        idx = np.arange(start, min(start + chunk, N), dtype=np.int64)
        bits = ((idx[:, None] >> ar) & 1).astype(np.float64)        # (b, n)
        out[start:start + len(idx)] = np.einsum("bi,ij,bj->b", bits, Q, bits)
    out += const
    return out


def _ising_terms(Q: np.ndarray):
    """Q -> lista de terminos Ising [(qubits, coef)], via qiskit (algebra x->z correcta).

    Devuelve (singles, pairs, offset): singles=[(q,h)], pairs=[(q,r,J)], offset escalar.
    """
    qp = fq.to_quadratic_program(Q)
    op, offset = qp.to_ising()
    singles, pairs = [], []
    for pauli, coeff in op.to_list():
        c = float(np.real(coeff))
        # pauli: label de n chars, el de mas a la derecha es el qubit 0
        zs = [q for q, ch in enumerate(reversed(pauli)) if ch == "Z"]
        if len(zs) == 1:
            singles.append((zs[0], c))
        elif len(zs) == 2:
            pairs.append((zs[0], zs[1], c))
        # len 0 (identidad) -> ya contemplado en offset
    return singles, pairs, float(np.real(offset))


def build_qaoa_circuit(n, singles, pairs, p):
    """Circuito QAOA PARAMETRICO: |+>ⁿ; por capa: e^{-iγ H_cost} (RZ/RZZ) y e^{-iβ ΣX} (RX).

    RZ(θ)=e^{-iθZ/2} => e^{-iγ·c·Z}=RZ(2γc); RZZ(θ)=e^{-iθZZ/2} => e^{-iγ·c·ZZ}=RZZ(2γc).
    Construccion manual (evita PauliEvolutionGate, que puede colgar; ver Georgia).
    Devuelve (qc, betas[list Parameter], gammas[list Parameter]) para bindear luego.
    """
    from qiskit import QuantumCircuit
    from qiskit.circuit import Parameter

    betas = [Parameter(f"b{i}") for i in range(p)]
    gammas = [Parameter(f"g{i}") for i in range(p)]
    qc = QuantumCircuit(n)
    qc.h(range(n))
    for beta, gamma in zip(betas, gammas):
        for q, h in singles:
            qc.rz(2.0 * gamma * h, q)
        for q, r, J in pairs:
            qc.rzz(2.0 * gamma * J, q, r)
        for q in range(n):
            qc.rx(2.0 * beta, q)
    qc.save_statevector()
    return qc, betas, gammas


class _Sim:
    """AerSimulator statevector con el circuito parametrico transpilado UNA vez.

    Cada evaluacion bindea (β,γ) y corre; evita re-transpilar por iteracion (clave
    para que T12/L3 = 24 qubits sea viable).
    """

    def __init__(self, n, singles, pairs, p):
        from qiskit import transpile
        from qiskit_aer import AerSimulator

        self.backend = AerSimulator(method="statevector")
        qc, self.betas, self.gammas = build_qaoa_circuit(n, singles, pairs, p)
        self.tqc = transpile(qc, self.backend)
        self.p = p

    def statevector(self, betas, gammas) -> np.ndarray:
        binding = {self.betas[i]: float(betas[i]) for i in range(self.p)}
        binding.update({self.gammas[i]: float(gammas[i]) for i in range(self.p)})
        bound = self.tqc.assign_parameters(binding)
        res = self.backend.run(bound).result()
        return np.asarray(res.get_statevector().data)


def run_qaoa(Q, const, meta, vi, levels, R_obs, dS, S0, *, S_min, S_max, u_max, B,
             weights, half, p=1, restarts=5, seed=42, maxiter=200):
    """Optimiza QAOA (COBYLA) sobre <H> y devuelve el mejor bitstring decodificado.

    <H> = Σ_x |ψ(x)|² · diag[x] (diagonal QUBO, mismo indexado que el statevector).
    Minimiza sobre (β,γ). Reusa `decode_and_verify` para factibilidad/SRS.
    """
    from scipy.optimize import minimize

    n = meta["n_qubits"]
    if n > MAX_QUBITS_STATEVECTOR:
        raise ValueError(f"n={n} qubits excede el tope statevector ({MAX_QUBITS_STATEVECTOR})")

    diag = precompute_diagonal(Q, const)          # (2ⁿ,) energias QUBO por estado base
    e_min, e_max = float(diag.min()), float(diag.max())
    singles, pairs, _offset = _ising_terms(Q)
    sim = _Sim(n, singles, pairs, p)

    def energies_of_params(params):
        psi = sim.statevector(params[:p], params[p:])
        probs = np.abs(psi) ** 2
        return probs, float(probs @ diag)

    best = {"expval": np.inf, "params": None, "probs": None, "curve": None}
    for r in range(restarts):
        # semilla reproducible por restart (Georgia: seed+r)
        rng_r = np.random.default_rng(seed + r)
        x0 = rng_r.uniform(0, np.pi, size=2 * p)
        curve = []

        def obj(params):
            _, ev = energies_of_params(params)
            curve.append(ev)
            return ev

        res = minimize(obj, x0, method="COBYLA", options={"maxiter": maxiter})
        probs, ev = energies_of_params(res.x)
        if ev < best["expval"]:
            best = {"expval": ev, "params": res.x, "probs": probs, "curve": curve,
                    "n_iter": len(curve)}

    # --- muestreo: top estados por probabilidad, decodificar+verificar ---
    probs = best["probs"]
    order = np.argsort(probs)[::-1][:256]          # top candidatos por probabilidad
    is_binary = hasattr(vi, "bits")
    best_dec = None
    for idx in order:
        xbits = np.array([(int(idx) >> i) & 1 for i in range(vi.n)], dtype=float)  # solo decision
        lv = vi.decode_levels(xbits)
        # validez del codeword: one-hot exacto por semana / nivel binario en rango
        if is_binary:
            valid = all(0 <= int(l) < vi.L for l in lv)
        else:
            valid = bool(vi.is_onehot(xbits))
        if not valid:
            continue                               # infactible por construccion; siguiente
        dv = sv.decode_and_verify(lv, vi, levels, R_obs, dS, S0, S_min=S_min,
                                  S_max=S_max, u_max=u_max, B=B, weights=weights)
        # energia QUBO con slack optimo (mismo criterio que exhaustive/SA)
        x = sv._x_from_levels(lv, vi, meta, half)
        e_dec = fq.qubo_energy(Q, x, const)
        cand = {"lv": lv, "prob": float(probs[idx]), "energy": e_dec, **dv}
        # preferir factible; entre factibles, menor energia
        if best_dec is None or (cand["feasible"], -cand["energy"]) > \
                (best_dec["feasible"], -best_dec["energy"]):
            best_dec = cand

    if best_dec is None:
        # ningun codeword valido en el top: reportar el mas probable, clamp e infactible
        idx = int(order[0])
        xbits = np.array([(idx >> i) & 1 for i in range(vi.n)], dtype=float)
        lv = np.clip(vi.decode_levels(xbits), 0, vi.L - 1)
        dv = sv.decode_and_verify(lv, vi, levels, R_obs, dS, S0, S_min=S_min,
                                  S_max=S_max, u_max=u_max, B=B, weights=weights)
        dv["feasible"] = False
        x = sv._x_from_levels(lv, vi, meta, half)
        best_dec = {"lv": lv, "prob": float(probs[idx]),
                    "energy": fq.qubo_energy(Q, x, const), **dv}

    betas = best["params"][:p].tolist()
    gammas = best["params"][p:].tolist()
    return {
        "expval": best["expval"], "e_min": e_min, "e_max": e_max,
        "betas": betas, "gammas": gammas, "n_iter": best["n_iter"],
        "convergence_curve": best["curve"], "p": p, "restarts": restarts, "seed": seed,
        "decoded": best_dec,
    }
