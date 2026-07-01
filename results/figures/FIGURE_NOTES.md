# Figure notes — Falcon presentation

> Generado por `scripts/julian/falcon_figures.py`. Cada figura: qué prueba, rúbrica, talking point y pro/con. PNGs en `results/figures/`.

## A2 — Storage vs physical bounds  (`A2...png`)
- **Proves:** Storage stays ~10-20% of S_max all year; min≈317Mm³ vs 0, max≈649Mm³ vs 3289Mm³.
- **Rubric:** Quantum Impl. (25%) — justifies storage_bounds='drop' (0 qubits)
- **Talking point:** Dropping 0≤S≤S_max costs zero optimality here.
- **Pro/con:** Pro: saves O(T) slacks. Con: must re-enable if a wetter regime approaches the bounds.

## A1/A3 — Which constraints bind + how relaxed  (`A1_A3...png`)
- **Proves:** Balance binds at 97–99% where the optimizer acts (medium/large); debug/small are u*=0 (drought) so balance is slack there. Storage bounds never bind; R≥0 forbids only 0–20 levels. (*table '~100%' = the binding regime.)
- **Rubric:** Quantum Impl. (25%) + Benchmarking (20%) — the relaxation rationale
- **Talking point:** Each relaxation is chosen by whether the constraint actually binds in our data (measured, not assumed).
- **Pro/con:** Pro: minimal qubits, provably lossless for the dropped ones. Con: data-regime-specific.

## A4 — Qubit count by encoding  (`A4...png`)
- **Proves:** one-hot medium=130q, large=260q; binary halves it but still >26 → chunking needed.
- **Rubric:** Quantum Impl. (25%) — justifies encoding + chunking choices
- **Talking point:** Compact encodings and per-block chunking are what make statevector QAOA feasible.
- **Pro/con:** Pro: binary fits small in statevector. Con: binary exact only for L≤4; medium needs chunking.

## A5 — Penalty magnitude is not the lever  (`A5...png`)
- **Proves:** Exhaustive-QUBO SRS is constant (=-0.29626) across 0.1–100× penalties (feasible).
- **Rubric:** Quantum Impl. (25%) — tuning insight
- **Talking point:** Penalty multiplier doesn't change the optimum; depth and the XY-mixer are the real levers.
- **Pro/con:** Pro: no brittle penalty tuning. Con: penalties still must be 'large enough' (very small breaks feasibility).

## B1 — ΔS cross-validated vs official IBWC  (`B1...png`)
- **Proves:** Weekly derived == −official: corr=+1.000, err≈0; official used raw drives storage negative.
- **Rubric:** Baseline (20%) + Problem Formulation (25%) — data validity
- **Talking point:** Our ΔS matches the official series once the sign convention is fixed; results aren't preliminary.
- **Pro/con:** Pro: data is validated + sign trap documented. Con: official daily series has gaps/opposite sign.

## B2 — Drought regime frames the problem  (`B2...png`)
- **Proves:** Storage 10–20% of S_max; 0 of 2182 days (2020–2026) above S_min → C_crit dominates the SRS.
- **Rubric:** Problem Formulation (25%) — societal/hydrological context
- **Talking point:** The reservoir is critically low all year, so avoiding deeper shortfall (C_crit) is the priority.
- **Pro/con:** Pro: explains why u=0 optima appear at coarse L. Con: window is a drought; wetter years may differ.

## C1 — SRS comparison per instance  (`C1...png`)
- **Proves:** DP ≥ all feasible baselines; threshold-pure/clamped are infeasible (balance violated).
- **Rubric:** Baseline (20%) + Benchmarking (20%)
- **Talking point:** The strong classical baseline (exact DP) is the honest bar the quantum method is compared to.
- **Pro/con:** Pro: valid metric + feasibility shown. Con: naive threshold rule looks good only if you ignore feasibility.

## D1 — Quantum vs baseline on the official benchmark  (`D1...png`)
- **Proves:** full DP=-0.2904 ≥ DP-chunked=-0.3115 ≥ QAOA-chunked=-0.3563 (infeasible); chunking gap vs QAOA gap separated.
- **Rubric:** Benchmarking (20%) + Quantum Impl. (25%) — the headline comparison
- **Talking point:** QAOA solves each block feasibly; the gap to full-DP is the price of chunking (a scaling limit), per spec §7.
- **Pro/con:** Pro: honest, feasible quantum result at 25q/block. Con: eta_local chunking caps ΔSRS; no quantum advantage (expected).

## D2 — QAOA quality + XY-mixer + hardware  (`D2...png`)
- **Proves:** debug AR=1.0 on sim AND IBM hardware; XY-mixer keeps 100% probability in the one-hot subspace.
- **Rubric:** Quantum Impl. (25%)
- **Talking point:** QAOA reaches the optimum on debug (incl. real hardware); XY-mixer removes the one-hot penalty cleanly.
- **Pro/con:** Pro: validated on hardware; XY-mixer principled. Con: AR degrades at p=1 as size grows.

## D4 — Window robustness (E1)  (`D4...png`)
- **Proves:** DP SRS varies modestly across windows; ranking of methods is stable.
- **Rubric:** Benchmarking (20%)
- **Talking point:** Results aren't an artifact of the start-of-year window; the benchmark is representative.
- **Pro/con:** Pro: robustness shown across 3 windows. Con: all windows are within the same drought year.

## E1 — Scaling: search space vs DP tractability  (`E1...png`)
- **Proves:** L^T from 243 (debug) to ~10⁴³ (large L7); DP runtime stays <0.3 s (polynomial O(T²L²)).
- **Rubric:** Benchmarking (20%) — scaling analysis
- **Talking point:** Exact DP makes the combinatorial problem tractable and gives ground truth at every size.
- **Pro/con:** Pro: DP is the strong scalable baseline. Con: QAOA can't match DP's scaling (statevector limit).

## E2 — ΔSRS headroom by instance  (`E2...png`)
- **Proves:** ΔSRS_vs_hist: debug/small≈0 (drought+coarse L), medium=+0.021, large=+0.051/+0.075.
- **Rubric:** Benchmarking (20%) + Problem Formulation (25%)
- **Talking point:** The debug instances are trivially u=0; the real optimization payoff is at L=5, longer horizons.
- **Pro/con:** Pro: shows where the method matters. Con: small ΔSRS overall (drought caps achievable gains).

## F1 — Optimal policy (medium)  (`F1...png`)
- **Proves:** DP conserves in ~11/26 weeks (u=−2Δu), lifting storage while respecting the balance budget.
- **Rubric:** Problem Formulation (25%) + Presentation (10%)
- **Talking point:** The optimizer redistributes releases in time (not just withholding) to reduce critical shortfall.
- **Pro/con:** Pro: interpretable, operationally meaningful. Con: gains are modest given the drought.
