"""Escritura estandarizada de resultados de experimentos Falcon.

Toda corrida (baselines, DP, QUBO, QAOA, ...) se registra con `record_run`, que
escribe en la raiz COMPARTIDA `results/` (no en carpetas personales):
  - JSON por-corrida con el detalle completo (incluye el schedule `u`):
    results/runs/{run_id}.json
  - una fila en el CSV maestro de columnas FIJAS (facil de comparar):
    results/runs_summary.csv

Asi cualquier metodo es comparable: mismas columnas, mismas metricas, mismos
nombres. Campos cuanticos van null para metodos clasicos. Convencion de `null`
y schema inspirados en el repo Georgia (docs/georgia_qubo_snippets.md).
"""
from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_ROOT = REPO_ROOT / "results"

# Orden FIJO de columnas del CSV maestro (no reordenar; solo agregar al final).
CSV_COLUMNS = [
    "run_id", "timestamp", "owner", "method", "variant", "instance", "T", "L",
    "official_status",
    "SRS", "Ccrit", "Cdev", "Csmooth",
    "dSRS_vs_historical", "dSRS_vs_threshold", "dSRS_vs_dp",
    "min_storage", "weeks_below_Smin", "release_balance_error", "feasible",
    "runtime_seconds",
    "delta_u", "u_max", "eta", "S_max", "S_min", "delta_u_full_ref",
    "w1", "w2", "w3",
    "encoding", "n_qubits", "p_depth", "energy", "approximation_ratio", "seed",
    "simulator", "penalties",
]


def instance_label(T: int, L: int) -> str:
    """Etiqueta estandar de instancia."""
    if (T, L) == (12, 3):
        return "small"
    if (T, L) == (26, 5):
        return "medium"
    if (T, L) == (52, 5) or (T, L) == (52, 7):
        return "large"
    return f"customT{T}L{L}"


def make_run_id(instance: str, T: int, L: int, method: str,
                variant: str | None, timestamp: str) -> str:
    """run_id estandar: {instance}_T{T}_L{L}_{method}[_{variant}]_{timestamp}."""
    parts = [instance, f"T{T}", f"L{L}", method]
    if variant:
        parts.append(variant)
    parts.append(timestamp)
    return "_".join(parts)


def record_run(*, method: str, instance: str, T: int, L: int,
               params: dict, weights: dict, constants: dict, B: float,
               u, S, costs: dict, srs: float, feasibility: dict,
               runtime_seconds: float,
               variant: str | None = None, official_status: str = "preliminary",
               references: dict | None = None, solver: dict | None = None,
               owner: str = "", results_root: Path = RESULTS_ROOT,
               timestamp: str | None = None) -> dict:
    """Registra una corrida: escribe JSON por-corrida y agrega fila al CSV maestro.

    `references` = {"historical": srs_h, "threshold": srs_t, "dp": srs_dp} (opcional)
    llena las columnas dSRS_vs_*. `solver` = campos cuanticos (opcional).
    Devuelve el record completo.
    """
    if timestamp is None:
        timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    u = np.asarray(u, dtype=float)
    S = np.asarray(S, dtype=float)
    references = references or {}
    solver = solver or {}

    run_id = make_run_id(instance, T, L, method, variant, timestamp)

    def dsrs(ref_key):
        ref = references.get(ref_key)
        return None if ref is None else float(srs - ref)

    S_min = constants["S_min_m3"]
    record = {
        "run_id": run_id,
        "timestamp": timestamp,
        "owner": owner,
        "method": method,
        "variant": variant,
        "instance": instance,
        "T": int(T),
        "L": int(L),
        "official_status": official_status,
        "SRS": float(srs),
        "Ccrit": float(costs["Ccrit"]),
        "Cdev": float(costs["Cdev"]),
        "Csmooth": float(costs["Csmooth"]),
        "dSRS_vs_historical": dsrs("historical"),
        "dSRS_vs_threshold": dsrs("threshold"),
        "dSRS_vs_dp": dsrs("dp"),
        "min_storage": float(S.min()),
        "weeks_below_Smin": int(np.sum(S < S_min)),
        "release_balance_error": float(max(0.0, abs(u.sum()) - B)),
        "feasible": bool(feasibility["feasible"]),
        "runtime_seconds": float(runtime_seconds),
        "delta_u": float(params["delta_u"]),
        "u_max": float(params["u_max"]),
        "eta": float(constants["eta"]),
        "S_max": float(constants["S_max_m3"]),
        "S_min": float(S_min),
        "delta_u_full_ref": float(params.get("delta_u_full_ref", params["delta_u"])),
        "w1": float(weights["w1"]),
        "w2": float(weights["w2"]),
        "w3": float(weights["w3"]),
        # Campos cuanticos (null si clasico)
        "encoding": solver.get("encoding"),
        "n_qubits": solver.get("n_qubits"),
        "p_depth": solver.get("p_depth"),
        "energy": solver.get("energy"),
        "approximation_ratio": solver.get("approximation_ratio"),
        "seed": solver.get("seed"),
        "simulator": solver.get("simulator"),
        "penalties": solver.get("penalties"),
        # Detalle solo en JSON (no en CSV plano):
        "u": u.tolist(),
        "B": float(B),
        "violations": feasibility.get("violations"),
        "references": {k: (None if v is None else float(v)) for k, v in references.items()},
    }

    # 1) JSON por-corrida
    runs_dir = Path(results_root) / "runs"
    runs_dir.mkdir(parents=True, exist_ok=True)
    (runs_dir / f"{run_id}.json").write_text(json.dumps(record, indent=2), encoding="utf-8")

    # 2) Fila en el CSV maestro (orden de columnas fijo)
    _append_summary_row(Path(results_root) / "runs_summary.csv", record)
    return record


def _append_summary_row(csv_path: Path, record: dict) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not csv_path.exists()
    row = {}
    for col in CSV_COLUMNS:
        val = record.get(col)
        if col == "penalties" and isinstance(val, dict):
            val = json.dumps(val, sort_keys=True)
        row[col] = "" if val is None else val
    with csv_path.open("a", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        if write_header:
            writer.writeheader()
        writer.writerow(row)
