# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Hackathon Puebla (Equipo 4) entry for the **Falcon Challenge** (`FalconChallenge/FalconChallenge_V6.pdf`): a guided quantum-computing challenge on transboundary water systems. The task is to find a reservoir **release-adjustment schedule** for the International Falcon Reservoir (Mexico–US Rio Grande / Río Bravo basin) that improves storage resilience, formulated as a discrete optimization problem suitable for QAOA / quantum annealing / hybrid quantum-classical methods.

The repo is in an early/setup stage: it contains the challenge spec (PDF), data-acquisition tooling, and quantum-execution smoke tests. There is no solution implementation yet — most work happens in Jupyter notebooks.

## Files

- `Ejemplos de ejecuciones.ipynb` — minimal smoke tests for the quantum stack: the same 4-qubit parametric circuit run on PennyLane (`lightning.gpu`), Qiskit Aer (`device="GPU"`), and CUDA Quantum (`nvidia` target). Use these to validate the environment and to fall back to CPU devices (`default.qubit` / `device="CPU"`) when no GPU is present.
- `FalconChallenge/Descargas.ipynb` — data downloader for challenge points 8 & 9. Pulls IBWC/USIBWC sources into `data/raw/`, scrapes link indexes, and writes a reproducible `data/falcon_download_manifest.json`. Intentionally avoids `beautifulsoup4`/`lxml`/`openpyxl` (uses stdlib `html.parser` + `requests`).
- `FalconChallenge/FalconChallenge_V6.pdf` — authoritative challenge spec. Read this before changing any math; the formulas and benchmark constants below come from it.

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

## Quantum environment

GPU-targeted by default. If running without an NVIDIA GPU or the GPU plugins:
- PennyLane: swap `qml.device("lightning.gpu", ...)` → `default.qubit` or `lightning.qubit`.
- Qiskit Aer: swap `device="GPU"` → `device="CPU"` (check `backend.available_devices()`).
- CUDA Quantum: requires a compatible NVIDIA install; no CPU fallback in the example.

Always include a warmup run before timing, as the smoke-test notebook does.

## Conventions

- Notebook prose and comments are in Spanish; keep new content consistent with the surrounding language.
- Avoid heavy parsing/Excel dependencies in the data tooling — the downloader deliberately sticks to the standard library plus `pandas`/`requests`.
