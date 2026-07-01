# Georgia VQ-MAR — buffered QUBO/QAOA code snippets

**Origin:** adapted reference from the sibling **Georgia VQ-MAR** repo
(`/Users/jmelmer/Documents/Personal/Jeff-QTeam/Repository-aquifer-recharge-georgia_data`, read-only).
These are **examples to adapt, not source of truth**. They are buffered here so we (almost) never
need to open the external repo. Original locations are noted per block for rare deep dives.

The problem-agnostic pieces (`qubo_energy`, `build_quadratic_program`, `verify_energy`,
`precompute_diagonal`, `build_qaoa_circuit`, `brute_force`) port to Falcon **unchanged** — they only
care about the matrix `Q`. Only `assemble_qubo` is domain-specific; for Falcon, `Q` is already built
by `scripts/falcon_qubo.py::build_qubo` (one-hot/binary; see `docs/SPEC_IMPLEMENTACION_QUBO.md`), then
feed it to everything below.

Convention used everywhere: `H(x) = xᵀQx + const`, minimize over `x ∈ {0,1}ⁿ`. For Falcon, bits are
one-hot `x_{t,ℓ}=1 ⟺ u(t)=aₗ`, so `n = T·L` (36 for T12/L3, 130 for T26/L5).

---

## 1. QUBO assembly — the diagonal/off-diagonal split

`xᵢ² = xᵢ` for binaries, so every squared linear penalty `(c + Σaᵢxᵢ)²` folds its linear-from-square
part onto the **diagonal** and its cross terms onto the **off-diagonal**; the scalar `c²` is a
**constant offset** tracked separately (it does not change the argmin but is needed for energy /
approximation-ratio reporting).

Georgia's site-selection version (`unified_pipeline.py:328`) — read as a template for the algebra:

```python
def assemble_qubo(wi, Ci, ei, M, budget, lam, mu):
    n = len(wi)
    Q = np.zeros((n, n))
    # Diagonal: linear coeff + linear-from-square of the budget penalty + exclusion penalty
    for i in range(n):
        Q[i, i] = -wi[i] + lam * (Ci[i]**2 - 2*budget*Ci[i]) + mu * ei[i]
    # Off-diagonal: pairwise coupling + cross terms from the squared budget penalty
    for i in range(n):
        for j in range(i + 1, n):
            Q[i, j] = 2 * lam * Ci[i] * Ci[j] + M[i, j]
            Q[j, i] = Q[i, j]
    const = lam * budget**2          # <- track separately, never fold into Q
    return Q, const
```

**Falcon mapping** (build the analogous fills from §11-§15):
- one-hot penalty `P(Σₗ x_{t,ℓ} − 1)²` → `−P` on each in-week diagonal, `+2P` on each in-week pair.
- `C_dev = Σ wₗ aₗ²` → diagonal `+w₂ aₗ²` per bit.
- storage `C_crit` and `C_smooth` are squares of expressions linear in bits (`S_t = H_t − Σ_{k<t}Σ_ℓ aₗ x_{k,ℓ}`, `uₜ − u_{t−1}`) → expand exactly like the budget square above; `C_smooth` couples **consecutive weeks**.
- release-balance `P_bal(Σuₜ)²`, non-negativity `P_R x_{t,ℓ}` where `R_obs+aₗ<0` (forbidden-state penalty, like Georgia's `μeᵢ`).

```python
def qubo_energy(x, Q):            # unified_pipeline.py:348 — ports unchanged
    return float(x @ Q @ x)
```

---

## 2. Exhaustive ground truth (small instance)

For Falcon T=12/L=3 the feasible space is `3¹² ≈ 531k` (one-hot reduces the `2³⁶` raw space); enumerate
to validate the QUBO. Georgia enumerates `2ⁿ` directly (`unified_pipeline.py:352`), auto-skipped above
N≈24 — for Falcon enumerate over the `Lᵀ` one-hot assignments instead, scoring each with `qubo_energy`.

```python
def brute_force(Q, feasible_iter):
    best_e, best_x = np.inf, None
    for x in feasible_iter:                 # iterate one-hot-valid bit vectors
        e = qubo_energy(x, Q)
        if e < best_e:
            best_e, best_x = e, x.copy()
    return {"energy": float(best_e), "x": best_x}
```

---

## 3. Q → Qiskit `QuadraticProgram` (the exact conversion)

`unified_pipeline`/`qiskit_qubo.py:150`. Diagonal → linear, off-diagonal (upper triangle, doubled
because Q is symmetric) → quadratic.

```python
from qiskit_optimization import QuadraticProgram

def build_quadratic_program(Q):
    n = Q.shape[0]
    qp = QuadraticProgram("falcon_qubo")
    for i in range(n):
        qp.binary_var(f"x{i}")
    linear = {f"x{i}": float(Q[i, i]) for i in range(n)}        # diagonal -> linear
    quadratic = {}
    for i in range(n):
        for j in range(i + 1, n):
            val = float(Q[i, j] + Q[j, i])                       # = 2*Q[i,j]
            if val != 0.0:
                quadratic[(f"x{i}", f"x{j}")] = val
    qp.minimize(linear=linear, quadratic=quadratic)
    return qp
```

---

## 4. Energy verification — DO THIS BEFORE QAOA

`qiskit_qubo.py:185`. Re-evaluate `xᵀQx` on a known-good bitstring (brute-force/SA optimum) and
confirm it matches the stored classical cost. Catches Q-assembly sign/scale bugs early
(Falcon spec §19 item 10, §22 item 7).

```python
def verify_energy(Q, bitstring, expected_energy, tol=1e-6):
    x = np.array([int(b) for b in bitstring], dtype=float)
    computed = float(x @ Q @ x)
    ok = abs(computed - expected_energy) < tol
    if not ok:
        print("WARNING: energy mismatch — check Q sign/scale convention.")
    return ok
```

---

## 5. Precomputed diagonal — the ~100× QAOA speedup

`run_qaoa_native_diagonal.py:95`. Compute all `2ⁿ` QUBO energies once (vectorized), then each COBYLA
step is just `dot(|ψ|², diagonal)` instead of a per-call DiagonalEstimator. Feasible memory only for
small n (n=20 → ~80 MB). For Falcon use on the small instance; for large n, sample instead.

```python
def precompute_diagonal(Q):
    n = Q.shape[0]; N = 2 ** n
    idx  = np.arange(N, dtype=np.uint32)
    bits = ((idx[:, None] >> np.arange(n, dtype=np.uint32)) & 1).astype(np.float32)
    qb   = bits @ Q.astype(np.float32)               # (N, n)
    return np.sum(bits * qb, axis=1).astype(np.float64)   # (N,)  energy per basis state
```

---

## 6. Manual QAOA circuit (avoids PauliEvolutionGate hangs)

`run_qaoa_native_diagonal.py:131`. Build the cost layer gate-by-gate from the Ising op; `QAOAAnsatz`/
`PauliEvolutionGate` cause multi-minute hangs in statevector sim. Param order matters: name them so
`beta` sorts before `gamma` alphabetically, matching `circuit.parameters` order.

```python
from qiskit import QuantumCircuit
from qiskit.circuit import ParameterVector

def build_qaoa_circuit(ising_op, reps):
    n = ising_op.num_qubits
    gamma = ParameterVector('gamma', reps)
    beta  = ParameterVector('beta', reps)            # 'beta' < 'gamma' -> stable sorted order
    z_terms, zz_terms = [], []
    for label, coeff in zip(ising_op.paulis.to_labels(), ising_op.coeffs):
        pos = [i for i, p in enumerate(reversed(label)) if p == "Z"]   # rightmost char = qubit 0
        if len(pos) == 1:   z_terms.append((pos[0], float(coeff.real)))
        elif len(pos) == 2: zz_terms.append((pos[0], pos[1], float(coeff.real)))
    qc = QuantumCircuit(n)
    qc.h(range(n))
    for layer in range(reps):
        g, b = gamma[layer], beta[layer]
        for qi, h in z_terms:           qc.rz(2.0 * g * h, qi)         # exp(-i g Z)
        for qi, qj, J in zz_terms:                                      # exp(-i g ZZ)
            qc.cx(qi, qj); qc.rz(2.0 * g * J, qj); qc.cx(qi, qj)
        for qi in range(n):             qc.rx(2.0 * b, qi)             # mixer exp(-i b X)
    qc.save_statevector()               # STRIP this for MPS / noisy sims (see pitfalls)
    return qc, gamma, beta
```

After solving, always **decode** samples back to `u(t)` and **verify feasibility** (one-hot per week,
release ≥ 0, balance) before reporting (spec §22 item 7).

---

## 7. Conventions worth copying

- `seed = 42` base; restart *r* uses `seed + r`. Run **≥5 restarts** whenever `Q` is ill-conditioned
  (Georgia: κ(Q_real)=214.8 vs 9.1 flat; single cold start reliably fails — expect this once Falcon
  adds slack-heavy constraints).
- Unified JSON result schema (use `null`, never omit): `run_id, variant, cost_model, n_sites,
  simulator, p_depth, energy, approximation_ratio, beta, gamma, top10_bitstrings, convergence_curve,
  seed, runtime_seconds`.
- One shared loader for Q/meta; never re-implement per script.
- Parallelism block Georgia uses (tune thread count to the machine):

```python
import multiprocessing, os
n_cores = multiprocessing.cpu_count()
for v in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(v, str(n_cores))
# AerSimulator(method="statevector", max_parallel_threads=n_cores, statevector_parallel_threshold=12)
```

---

## 8. Code-level pitfalls (qiskit)

- **Track `const` separately** from `Q` (Georgia: `const=λB²`) — required for correct energy/AR.
- **MPS / noisy Aer**: `build_qaoa_circuit` appends `save_statevector()`, which is incompatible.
  Strip it from `circuit.data` (skip `save_statevector` instructions, rebuild), then re-add as needed.
- **`SamplerV2`** (qiskit-ibm-runtime ≥0.45): `SamplerV2(mode=backend)`, not `SamplerV2(backend=backend)`.
- **DataBin register name**: don't hardcode `data_bin.meas`; use
  `creg = next(iter(vars(data_bin))); counts = getattr(data_bin, creg).get_counts()`.
- Start on a **simulator with the small instance** (T=12, L=3); never debug on hardware.

---

## Deep-reference index (open the external repo only if a snippet above is insufficient)

`CLAUDE.md` (QUBO formulation, params, pitfalls) · `scripts/CLAUDE.md` (script rules, shared-function
index) · `scripts/unified_pipeline.py` (`assemble_qubo:328`, `qubo_energy:348`, `brute_force:352`,
`greedy_solve:388`, `sa_solve:418`, `ga_solve:459`, `sanity_check:516`) · `scripts/qaoa/qiskit_qubo.py`
(`build_quadratic_program:150`, `verify_energy:185`, `load_qubo:59`) · `scripts/qaoa/run_qaoa_native_diagonal.py`
(`precompute_diagonal:95`, `build_qaoa_circuit:131`, `run_qaoa:191`, `decode_result:327`, `compute_metrics:367`).
