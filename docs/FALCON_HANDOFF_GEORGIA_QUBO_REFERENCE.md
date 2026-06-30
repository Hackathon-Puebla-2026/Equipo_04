# Georgia VQ-MAR — QUBO/QAOA Reference for the Falcon Challenge

**Purpose.** This document points the FalconChallenge implementation at a sibling project that
already solves the same class of problem (constrained binary optimization → QUBO → Ising → QAOA).
Use it as a *reference implementation library*. Do not reinvent QUBO assembly, classical
baselines, or the QAOA harness — adapt what already works.

## Where the reference code lives

```
/Users/jmelmer/Documents/Personal/Jeff-QTeam/Repository-aquifer-recharge-georgia_data
```

This is the **Georgia VQ-MAR** repo. It is on the same machine. You may **read** any file in it
for reference. **Never modify, run, or commit anything in it** — it is a separate project. All
Falcon data-safety rules still apply: write only inside the FalconChallenge `data/processed/`,
`results/`, and `src/` directories.

Start by reading these in the Georgia repo:
- `CLAUDE.md` (root) — QUBO formulation, parameters, schema, pitfalls
- `scripts/CLAUDE.md` — script-level rules and shared-function index
- `scripts/unified_pipeline.py` — QUBO assembly + all classical solvers
- `scripts/qaoa/qiskit_qubo.py` — QUBO → Qiskit `QuadraticProgram`
- `scripts/qaoa/run_qaoa_native_diagonal.py` — production QAOA solver

## The one transferable idea

Both projects build a single matrix `Q` such that `H(x) = xᵀQx + const`, then minimize over
`x ∈ {0,1}ⁿ`. Every cost term and every constraint is written as a **linear expression in the
binary variables, then squared**. Because `xᵢ² = xᵢ` for binaries, the square's linear-from-square
part folds onto the **diagonal** of Q and the cross terms become **off-diagonal** entries. This is
exactly the `add_square_of_linear_expression` helper the Falcon spec (§12) asks you to write — the
Georgia repo already implements the same algebra inline.

Reference: `scripts/unified_pipeline.py:328` (`assemble_qubo`).
```
Diagonal:     Q[i,i] = (linear coeff for x_i)  +  (constraint penalties on x_i)
Off-diagonal: Q[i,j] = (pairwise coupling)     +  (cross terms from squared penalties)
const:        scalar energy offset, tracked separately (matters for energy/AR comparisons)
```

## Reusable building blocks (with exact locations)

| Falcon need (from Falcon CLAUDE.md) | Georgia reference | Location |
|---|---|---|
| `add_square_of_linear_expression` helper (§12) | the `(c+Σaᵢxᵢ)²` expansion onto diagonal/off-diagonal | `unified_pipeline.py:328` `assemble_qubo`, `:339` cross terms |
| QUBO energy eval `qubo_energy(Q,x)` (§17.5) | identical helper | `unified_pipeline.py:348` |
| Exhaustive search for small instance (§9.3, T=12 L=3) | brute-force enumeration, auto-skipped above N=24 | `unified_pipeline.py:352` `brute_force` |
| Simulated annealing baseline (§9.3) | geometric-cooling SA, multi-restart | `unified_pipeline.py:418` `sa_solve` |
| Evolutionary algorithm baseline (§9.3) | genetic algorithm, tournament selection | `unified_pipeline.py:459` `ga_solve` |
| Greedy/threshold-style baseline (§9.2) | ratio-ranked greedy | `unified_pipeline.py:388` `greedy_solve` |
| Decode + verify feasibility after solving (§17.6, §22.7) | selection/cost/exclusion/energy-match checks | `unified_pipeline.py:516` `sanity_check` |
| QUBO → Qiskit `QuadraticProgram` (§15) | diagonal→linear, off-diagonal→quadratic | `qiskit_qubo.py:150` `build_quadratic_program` |
| QUBO → Ising + QAOA on simulator (§15) | full QAOA build/run/decode | `run_qaoa_native_diagonal.py` (`build_qaoa_circuit`, `precompute_diagonal`, `decode_result`, `compute_metrics`) |
| Validate QUBO energy == classical cost (§19.10) | energy verification with tolerance | `qiskit_qubo.py:185` `verify_energy` |

## Mapping Georgia patterns → Falcon's one-hot formulation

The structures differ; the machinery is the same.

| | Georgia (site selection) | Falcon (release schedule) |
|---|---|---|
| Bit meaning | `xᵢ = 1` ⟺ select site *i* | `x_{t,ℓ} = 1` ⟺ `u(t) = aₗ` (one-hot) |
| Bit count | N sites (20, 50) | T×L (36 for T12/L3, 130 for T26/L5) |
| Diagonal source | `−wᵢ + λ(Cᵢ²−2BCᵢ) + μeᵢ` | one-hot `−P` per bit, `+wₗ aₗ²` from `C_dev`, storage/`C_crit` linear-from-square |
| Off-diagonal source | spatial `Mᵢⱼ = η e^(−dᵢⱼ/ℓ)` + budget coupling `2λCᵢCⱼ` | one-hot pair `+2P`, `C_smooth` consecutive-week coupling, storage cross terms |
| Equality constraint | budget `λ(ΣCᵢxᵢ−B)²` | one-hot `P(Σₗ x_{t,ℓ}−1)²`, release balance `P_bal(Σuₜ)²` |
| Forbidden states | exclusion penalty `μΣeᵢxᵢ` | invalid release `P_R x_{t,ℓ}` where `R_obs+aₗ < 0` |
| Ground truth | brute force `2^N` (N≤24) | exhaustive `L^T` (3¹²=531k feasible; 5²⁶ not) |
| MVP vs faithful | flat vs real cost model | soft-storage (Opt. A) vs deficit/surplus slacks (Opt. B) |

**Concrete porting notes:**
- Falcon's `S_t = H_t − Σ_{k<t} Σ_ℓ aₗ x_{k,ℓ}` is linear in bits → any squared storage cost is
  QUBO-able exactly the way Georgia squares the budget term `(ΣCᵢxᵢ − B)²`.
- Falcon's one-hot penalty `P(Σₗ x_{t,ℓ}−1)²` expands to `−P` on each diagonal entry and `+2P` on
  each in-week pair — the same diagonal/off-diagonal split Georgia uses for its budget penalty.
- `C_smooth = Σ(uₜ−u_{t−1})²` produces couplings between *consecutive weeks* — Georgia's nearest
  analogue is the spatial coupling `Mᵢⱼ`; reuse the off-diagonal-fill loop pattern.

## QAOA harness conventions worth copying

The Georgia QAOA scripts converged on conventions that will save Falcon time (see
`scripts/CLAUDE.md`):
- Base `seed = 42`; restart *r* uses `seed + r`.
- One **shared loader** (`qiskit_qubo.load_qubo`) — never re-implement Q/meta loading per script.
- A **unified JSON output schema** (mandatory fields: `run_id, variant, cost_model, n_sites,
  simulator, p_depth, energy, approximation_ratio, beta, gamma, top10_bitstrings,
  convergence_curve, seed, runtime_seconds`); use `null`, never omit.
- Always **decode and verify feasibility** of QAOA samples before reporting (mirrors Falcon §22.7).
- Available simulator variants to model your own after: native-diagonal (production, ~100×
  speedup via precomputed diagonal), CVaR-QAOA, MPS (`χ=64`), noisy Aer (FakeBrisbane), and real
  IBM hardware. Start on a simulator with the small instance (T=12, L=3), never hardware.

## Pitfalls that transfer

- **Penalty weights** (Falcon §14): too small → infeasible optima; too large → wrecks
  conditioning. Georgia learned this the hard way — its real cost model has κ(Q)=214.8 (vs 9.1
  flat) and *single cold-start QAOA reliably fails*; it always uses ≥5 restarts. Expect the same
  once Falcon adds slack-heavy constraints (Option B), and scale restarts accordingly.
- **Track the constant offset** separately from Q (Georgia: `const = λB²`). It does not affect the
  argmin but is required for correct energy / approximation-ratio reporting.
- **Validate energy == decoded classical cost** on several bitstrings before trusting QAOA
  (Falcon §19.10 and Georgia `verify_energy`). Do this first; it catches Q-assembly bugs early.
