# Figure notes — Falcon presentation

> Generado por `scripts/julian/falcon_figures.py`. Cada figura: qué prueba, rúbrica, talking point y pro/con. PNGs en `results/figures/`.

## Feasibility explainer (para C1 / D1 / D4)

Una política `u(t)` es **factible** solo si cumple las **4 restricciones oficiales**: `R(t)=R_obs+u≥0`, `|u(t)|≤u_max`, `0≤S(t)≤S_max`, y **`|Σu(t)|≤B=η·ΣR_obs`** (presupuesto de balance). En nuestros datos las primeras tres casi nunca atan: **la que decide factibilidad es el balance**. *Factibilidad y SRS son ejes independientes*: una política puede ser factible y peor que el histórico (p.ej. threshold-balanced en debug/small — legítimo, en sequía `u=0` ya es óptimo), y una **infactible puede tener un SRS más alto solo por gastar de más el presupuesto** (p.ej. threshold-pure en medium, −0.289, pero `|Σu|−B>0`). **Regla de lectura: comparar SRS solo entre barras factibles.**

## Definición de ventanas (para D4)

*first* = semanas `[0,T)`; *middle* = start `(52−T)//2`; *stress* = la ventana de T semanas con **menor storage medio** (déficit más profundo, la más exigente).

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

## A4 — Qubit count by encoding (one-hot / domain-wall / binary)  (`A4...png`)
- **Proves:** medium T26/L5: one-hot=130q, domain-wall=104q, binary=78q — all >30 → chunking. Ivan's IBM hardware run used domain-wall (4 bits/week) in blocks 7+7+7+5 (≤30q).
- **Rubric:** Quantum Impl. (25%) — justifies encoding + chunking choices
- **Talking point:** Compact encodings + per-block chunking are what make QAOA fit a real device / statevector.
- **Pro/con:** Pro: domain-wall enabled the actual IBM-hardware run; binary fits small in statevector. Con: binary exact only for L≤4; every full instance still needs chunking.

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

## C1 — SRS comparison per instance (feasibility-aware)  (`C1...png`)
- **Proves:** Among FEASIBLE methods, DP ≥ historical ≥ threshold-balanced; threshold-pure/clamped look higher (e.g. medium −0.289) but are INFEASIBLE (|Σu|−B>0). Feasible-but-worse-than-hist (debug/small threshold-balanced) is legitimate — u=0 is optimal in the drought.
- **Rubric:** Baseline (20%) + Benchmarking (20%)
- **Talking point:** Compare only feasible bars; exact DP is the honest strong baseline for the quantum comparison.
- **Pro/con:** Pro: valid metric with feasibility made explicit. Con: naive threshold 'wins' only if you ignore the balance constraint it violates.

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

## D4 — Window robustness (E1) — DP vs historical per window  (`D4...png`)
- **Proves:** ΔSRS(DP−hist): debug/small ≈ 0 in ALL windows (u=0 optimal, drought+coarse L); medium >0 in every window → the optimum is genuinely non-trivial there, not a start-of-year artifact.
- **Rubric:** Benchmarking (20%)
- **Talking point:** Windows: first=weeks[0,T); middle=start (52−T)//2; stress=T-week window of lowest mean storage (deepest deficit). Comparing DP to historical per window shows WHERE optimization actually helps.
- **Pro/con:** Pro: separates 'found optimum' from 'optimum beats baseline' across windows. Con: all windows are within one drought year.

## E1 — Scaling: L^T search space vs DP states explored (the collapse)  (`E1...png`)
- **Proves:** medium L^T=1.49e18 → only 6,505 DP states (0.02s); large L7 8.81e43 → 49,739 states (0.20s). DP state=(t, C_t=Σk_j, k_prev): storage depends only on the integer cumulative sum C_t, so exponentially many schedules fuse into O(T²·L²) states — exact, lossless (brute==dp).
- **Rubric:** Benchmarking (20%) — scaling analysis
- **Talking point:** Exploiting optimal substructure turns a 10⁴³ search into ~10⁴ states: exact ground truth at every size.
- **Pro/con:** Pro: DP is exact + sub-second at all scales. Con: it's a classical baseline; QAOA can't match its scaling.

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
