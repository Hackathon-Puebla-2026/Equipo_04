# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Hackathon Puebla (Equipo 4) entry for the **Falcon Challenge** (`FalconChallenge/FalconChallenge_V6.pdf`): a guided quantum-computing challenge on transboundary water systems. The task is to find a reservoir **release-adjustment schedule** for the International Falcon Reservoir (Mexico–US Rio Grande / Río Bravo basin) that improves storage resilience, formulated as a discrete optimization problem suitable for QAOA / quantum annealing / hybrid quantum-classical methods.

The repo is in an early/setup stage: it contains the challenge spec (PDF), data-acquisition tooling, and quantum-execution smoke tests. There is no solution implementation yet — most work happens in Jupyter notebooks.

## Files

- `Ejemplos de ejecuciones.ipynb` — minimal smoke tests for the quantum stack: the same 4-qubit parametric circuit run on PennyLane (`lightning.gpu`), Qiskit Aer (`device="GPU"`), and CUDA Quantum (`nvidia` target). Use these to validate the environment and to fall back to CPU devices (`default.qubit` / `device="CPU"`) when no GPU is present.
- `FalconChallenge/Descargas.ipynb` — data downloader for challenge points 8 & 9. Pulls IBWC/USIBWC sources into `data/raw/`, scrapes link indexes, and writes a reproducible `data/falcon_download_manifest.json`. Intentionally avoids `beautifulsoup4`/`lxml`/`openpyxl` (uses stdlib `html.parser` + `requests`).
- `FalconChallenge/Smax_search.ipynb` — resolves `S_max` (Falcon conservation capacity) two independent ways: scraping the live IBWC reports and deriving it from `storage / (Conservation% / 100)`. Writes `data/falcon_reservoir_constants.json`.
- `FalconChallenge/data/falcon_reservoir_constants.json` — resolved model constants: `S_max ≈ 3,288,726 TCM` (3,288.726 MCM), `S_min = 0.25·S_max ≈ 822,181.5 TCM`, flood capacity 3923 MCM (context only). Both derivation methods agree exactly.
- `FalconChallenge/src/inspect_falcon_data.py` — helper to inspect the exported CSV datasets.
- `FalconChallenge/FalconChallenge_V6.pdf` — authoritative challenge spec. Read this before changing any math; the formulas and benchmark constants below come from it.
- `requirements.txt` — base Python deps; GPU variants (`lightning.gpu`, Aer-GPU, `cudaq`) are optional and machine-dependent.
- `notebooks/ scripts/ results/` — per-person work folders (`julian/ emilio/ ivan/`) plus shared code at the root of `scripts/`. See `docs/GUIDELINES.md` for ownership rules.

### Docs (`docs/`)

- `docs/FalconChallenge_V6.md` — markdown transcription of the spec PDF (cheaper to load and greppable). **Auto-imported.** The PDF remains the source of truth on any disagreement.
- `docs/GUIDELINES.md` — team collaboration guidelines: frozen `FalconChallenge/`, per-person folder ownership, where shared code lives, results/data conventions. **Auto-imported.** Follow it when creating or placing files.
- `docs/FALCON_HANDOFF_GEORGIA_QUBO_REFERENCE.md` — lean, self-contained buffer of the transferable QUBO/QAOA patterns from the sibling Georgia VQ-MAR project (the one transferable idea, the Georgia→Falcon one-hot mapping, conventions, pitfalls). **Auto-imported.**
- `docs/georgia_qubo_snippets.md` — buffered, ready-to-adapt code from Georgia (QUBO assembly, `precompute_diagonal`, manual QAOA circuit, energy verification). **On-demand** — read it when implementing the QUBO/QAOA pieces.
- `docs/CLAUDE_Falcon_QUBO_Input.md` — large implementation guide (EDA, classical baselines, full QUBO/Ising formulation §1-23, helper code, expected outputs). **On-demand** (not auto-imported, to save tokens) — read it when doing QUBO/EDA/baseline work. Guidance/example, not binding.
- `docs/HALLAZGOS_CLAVE.md` — decision-critical digest of what ivan & emilio have done (data ready, preliminary-vs-official `S_max`/`S_min`, `Δu` convention, one-hot encodings, no QAOA yet). **Auto-imported** (always loaded).
- `docs/RESUMEN_HALLAZGOS_EQUIPO.md` — full detailed summary of the team's EDA + preliminary QUBO work, with numbers, discrepancy table, and reusable-asset inventory. **On-demand** — read for the full picture.
- `docs/SPEC_IMPLEMENTACION_QUBO.md` — **the active implementation spec**: locked conventions, module architecture (`scripts/falcon_*.py`), phased todo checklist with done-criteria, exact-lattice DP, rubric→tasks mapping. **On-demand** — read before/while coding any pipeline part.
- `docs/ANALISIS_DP_Y_RESULTADOS.md` — DP complexity/timing analysis (states-explored vs `Lᵀ`, why sub-second is expected, lossless pruning) + current results digest (baselines/DP across instances, ΔSRS, feasibility, cross-validation vs ivan) + §8 qubit-scaling & droppable constraints. **On-demand.** Raw data in `results/runs_summary.csv`.
- `docs/RENDIMIENTO_M4.md` — how to run performantly on this laptop (**CPU-first M4; GPU server unreliable**): backend choice (`lightning.qubit`), thread tuning, ~28-30 qubit memory ceiling, compact encoding + chunking, vectorized classical, env setup. **On-demand** — read before Fase 3 / any simulation.

@docs/FalconChallenge_V6.md
@docs/GUIDELINES.md
@docs/FALCON_HANDOFF_GEORGIA_QUBO_REFERENCE.md
@docs/HALLAZGOS_CLAVE.md

## Reference material (examples, not source of truth)

Keep this hierarchy in mind — most docs are *guidance and worked examples*, not binding spec:

- **Source of truth:** `FalconChallenge/FalconChallenge_V6.pdf` + its transcription `docs/FalconChallenge_V6.md`, plus the official benchmark constants in this file. On any disagreement, the spec wins.
- **On-demand local guidance:** `docs/CLAUDE_Falcon_QUBO_Input.md` — a *suggested* QUBO formulation and file layout. Adopt what helps; deviate when the spec or our data says otherwise. Read it when implementing QUBO/EDA/baselines (it is not auto-imported, to keep context light).
- **Buffered QUBO/QAOA patterns:** `docs/FALCON_HANDOFF_GEORGIA_QUBO_REFERENCE.md` (auto-imported, lean) + `docs/georgia_qubo_snippets.md` (on-demand code). These buffer the reusable machinery from the sibling Georgia VQ-MAR project locally, so you should rarely need the external repo.
- **External deep reference (rarely needed):** the Georgia repo at `/Users/jmelmer/Documents/Personal/Jeff-QTeam/Repository-aquifer-recharge-georgia_data`. **Read-only**: never modify, run, or commit anything there. Read **specific** files on demand only — never bulk-load the repo into context. Treat its code as a pattern to adapt, not a spec to follow.

For file placement, `docs/GUIDELINES.md` is authoritative — note the Georgia handoff assumes paths like `FalconChallenge/src/` and `data/processed/` that differ from our actual layout; follow GUIDELINES.

## The optimization problem (from the spec)

Maximize the **Storage Resilience Score**: `SRS = -(w1·C_crit + w2·C_dev + w3·C_smooth)`.

- Decision variable is the release adjustment `u(t)`; the optimized release is `R(t) = R_obs(t) + u(t)`, so `u(t)=0` reproduces historical operation.
- Storage dynamics: `S_opt(t+1) = S_opt(t) + ΔS_obs(t) - u(t)`, using observed storage change from `Discharge.Total.Change-in-Storage@08461200`.
- Cost terms: `C_crit` penalizes storage below `S_min`; `C_dev = Σ u(t)²`; `C_smooth = Σ (u(t)-u(t-1))²`.
- Constraints: `R(t) ≥ 0`; `|u(t)| ≤ u_max`; `0 ≤ S_opt(t) ≤ S_max`; cumulative balance `|Σ u(t)| ≤ η·Σ R_obs(t)`.

**Official benchmark constants** (use these for the scored result; report any alternatives separately):
- `L = 5` adjustment levels: `u(t) ∈ {-2Δu, -Δu, 0, Δu, 2Δu}`, with `Δu = 0.25·median weekly observed release`, `u_max = 2Δu`.
- Weights: `w1 = 1/((T+1)·S_min²)`, `w2 = 0.1/(T·u_max²)`, `w3 = 0.1/((T-1)·(2·u_max)²)`.
- `η = 0.10`; `S_min = 0.25·S_max` (S_max = Falcon total conservation storage capacity).

**Baselines to compare against**: historical replay (`u=0`), the threshold conservation rule (`u_rule = -Δu` when `S < S_min`, else `0`), and optionally a stronger classical optimizer (DP / MILP / SA / evolutionary). The reported quantity is `ΔSRS = SRS_opt - SRS_baseline`, plus runtime and scaling behavior.

**Scaling instances**: small `T=12, L=3` (debug); medium `T=26, L=5` (official benchmark); large `T=52, L=5 or 7` (scaling analysis). Candidate schedules scale as `L^T`.

A quantum/hybrid solution must discretize `u(t)`, encode it as binary variables (QUBO/Ising), and justify the structure. The benchmark is specified but the solution is not provided.

## Data

Datasets come from the IBWC water-data portal. Two stations matter:
- `08461200` (International Falcon Reservoir) — storage, elevation, lake area, evaporation, storage change.
- `08461300` (Rio Grande Below Falcon Dam) — `Discharge.Best Available@08461300` is the observed historical release `R_obs`.

The official benchmark dataset is provided in a shared SharePoint folder (see `Descargas.ipynb` `source_urls`); downloading from IBWC directly is optional/for scaling extensions. If discharge is in m³/s, convert to volume over the time step before storage-balance modeling. `Descargas.ipynb` writes everything under `data/` (gitignored alongside `venv/`, `__pycache__/`, checkpoints).

## Compute environments

**CPU-first: the primary target is the M4 laptop.** The GPU server has availability problems, so do not depend on it. Full performance guidance: `docs/RENDIMIENTO_M4.md`.

- **Local — MacBook Pro M4 (Apple Silicon), primary.** No NVIDIA GPU. Use **CPU** simulators: prefer PennyLane `lightning.qubit` (native ARM C++); Qiskit Aer `device="CPU"` only if wheels install (may need a Python 3.11/3.12 venv). `lightning.gpu`/CUDA Quantum do NOT run locally. Many-core CPU is good for statevector-CPU sim and the precomputed-diagonal trick. Tune thread env vars; avoid oversubscription (see `docs/RENDIMIENTO_M4.md`).
- **Remote GPU — WCentroid cluster (unreliable, optional).** NVIDIA GPUs incl. a **T4 (16 GB)** where `lightning.gpu`/Aer-GPU/CUDA-Q *could* run, but availability is spotty — treat as a bonus, not a dependency. Keep device behind a flag so the same code runs CPU here and GPU there.

**Memory reality (drives encoding + instance choices):** dense statevector uses `2ⁿ × 16 bytes`. A 16 GB GPU (T4) tops out around **~30 qubits**; the precomputed-diagonal array (`2ⁿ × 8 bytes`) is similar. Implications:
- Small instance T12/L3: **one-hot = 36 qubits → too big for statevector**; domain-wall or binary = 24 qubits → fits (~256 MB). So compact encodings aren't just nicer, they're what makes exact statevector QAOA feasible.
- Medium (130 one-hot) / large: far beyond any statevector. Use **MPS / tensor-network** simulators (the 1-D temporal chain suits MPS well), sampling, or classical annealing — not brute statevector.

Portability: keep device selection behind a flag/config so the same code runs CPU locally and GPU on WCentroid. Always include a warmup run before timing, as the smoke-test notebook does.

## Conventions

- Notebook prose and comments are in Spanish; keep new content consistent with the surrounding language.
- Avoid heavy parsing/Excel dependencies in the data tooling — the downloader deliberately sticks to the standard library plus `pandas`/`requests`.
