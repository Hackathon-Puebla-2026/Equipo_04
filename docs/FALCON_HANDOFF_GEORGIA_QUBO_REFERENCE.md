# Georgia VQ-MAR — QUBO/QAOA reference for Falcon (buffered locally)

A sibling project (**Georgia VQ-MAR**) already solves the same class of problem: constrained binary
optimization → QUBO → Ising → QAOA. Reuse its machinery instead of reinventing QUBO assembly,
classical baselines, or the QAOA harness. The reusable parts are **buffered locally** so you rarely
need the external repo:

- **Code** → `docs/georgia_qubo_snippets.md` (on-demand): QUBO assembly pattern, `qubo_energy`,
  `build_quadratic_program`, `verify_energy`, `precompute_diagonal`, manual `build_qaoa_circuit`,
  conventions, pitfalls.
- **External repo** (deep reference only, read-only, never bulk-load):
  `/Users/jmelmer/Documents/Personal/Jeff-QTeam/Repository-aquifer-recharge-georgia_data`. Index at
  the bottom of the snippets file.

These are **examples to adapt, not source of truth** — the Falcon spec (`docs/FalconChallenge_V6.md`),
`docs/SPEC_IMPLEMENTACION_QUBO.md`, and the implemented `scripts/falcon_qubo.py` define the actual formulation.

## The one transferable idea

Build a single matrix `Q` with `H(x) = xᵀQx + const`, minimize over `x ∈ {0,1}ⁿ`. Each cost term and
constraint is written as a **linear expression in the bits, then squared**. Since `xᵢ²=xᵢ`, the
square's linear part lands on the **diagonal** of `Q`, cross terms on the **off-diagonal**, and the
scalar `c²` is a **constant offset** tracked separately (matters for energy / approximation ratio,
not for the argmin). This is exactly the `add_square_of_linear_expression` helper the Falcon spec
(§12) asks for.

## Georgia → Falcon mapping (same machinery, different structure)

| | Georgia (site selection) | Falcon (release schedule) |
|---|---|---|
| Bit meaning | `xᵢ=1` ⟺ select site *i* | `x_{t,ℓ}=1` ⟺ `u(t)=aₗ` (one-hot) |
| Bit count | N sites (20, 50) | T×L (36 for T12/L3, 130 for T26/L5) |
| Diagonal | `−wᵢ + λ(Cᵢ²−2BCᵢ) + μeᵢ` | one-hot `−P` per bit, `+w₂ aₗ²` from `C_dev`, storage/`C_crit` linear-from-square |
| Off-diagonal | spatial `Mᵢⱼ` + budget coupling `2λCᵢCⱼ` | one-hot pair `+2P`, `C_smooth` consecutive-week coupling, storage cross terms |
| Equality constraint | budget `λ(ΣCᵢxᵢ−B)²` | one-hot `P(Σₗ x_{t,ℓ}−1)²`, release balance `P_bal(Σuₜ)²` |
| Forbidden states | exclusion `μΣeᵢxᵢ` | invalid release `P_R x_{t,ℓ}` where `R_obs+aₗ<0` |
| Ground truth | brute force `2ⁿ` (N≤24) | exhaustive `Lᵀ` (3¹²≈531k feasible; 5²⁶ not) |

**Porting notes:** Falcon storage `S_t = H_t − Σ_{k<t}Σ_ℓ aₗ x_{k,ℓ}` is linear in bits, so any squared
storage cost is QUBO-able exactly as Georgia squares `(ΣCᵢxᵢ−B)²`. The one-hot penalty expands to `−P`
on each in-week diagonal and `+2P` on each in-week pair (same split as Georgia's budget penalty).
`C_smooth = Σ(uₜ−u_{t−1})²` couples consecutive weeks — reuse Georgia's off-diagonal-fill loop.

## Conventions and pitfalls that transfer

- `seed=42` base; restart *r* uses `seed+r`. Use **≥5 restarts** when `Q` is ill-conditioned — Georgia's
  real cost model has κ(Q)=214.8 (vs 9.1) and single cold-start QAOA reliably fails; expect the same
  once Falcon adds slack-heavy constraints (Option B).
- **Validate energy == decoded classical cost first** (Falcon §19 item 10) — catches Q-assembly bugs
  before any QAOA. Then always **decode + verify feasibility** of samples before reporting (§22 item 7).
- **Track the constant offset** separately from `Q`.
- Unified JSON result schema (use `null`, never omit): `run_id, variant, cost_model, n_sites, simulator,
  p_depth, energy, approximation_ratio, beta, gamma, top10_bitstrings, convergence_curve, seed,
  runtime_seconds`.
- Start on a **simulator with the small instance** (T=12, L=3), never hardware. Penalty weights too
  small → infeasible optima; too large → wrecks conditioning (Falcon §14).

## Data-safety

Write only inside the Falcon repo (per `docs/GUIDELINES.md`: shared code in root `scripts/`, per-person
folders, `results/`, `data/`). Never modify/run/commit anything in the Georgia repo. GUIDELINES.md is
authoritative for file placement (the Georgia repo's own paths like `src/`, `data/processed/` differ).
