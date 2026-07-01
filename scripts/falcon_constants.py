"""Constantes oficiales y parametros por instancia del benchmark Falcon.

Todo en Sistema Internacional (m^3). Las constantes oficiales se leen de
``FalconChallenge/data/falcon_reservoir_constants.json`` (resueltas en Smax_search.ipynb);
NO se usa el maximo observado como S_max (eso era preliminar).

Referencias de formulas: docs/SPEC_IMPLEMENTACION_QUBO.md secciones 2 y 5.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

# Raiz del repo = dos niveles arriba de este archivo (scripts/ -> repo).
REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONSTANTS_PATH = REPO_ROOT / "FalconChallenge" / "data" / "falcon_reservoir_constants.json"

# Factor de conservacion critico: S_min = CRIT_FRACTION * S_max  (spec).
CRIT_FRACTION = 0.25
# Tolerancia de balance acumulado: |sum u| <= ETA * sum R_obs  (spec).
ETA = 0.10


def load_official_constants(path: str | Path = DEFAULT_CONSTANTS_PATH) -> dict:
    """Carga las constantes oficiales del embalse y las devuelve en m^3.

    Returns dict con: ``S_max_m3``, ``S_min_m3`` (= 0.25 * S_max),
    ``flood_capacity_m3`` (solo contexto, no usado por el SRS) y ``source``.
    """
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    s_max_m3 = float(data["s_max_m3"])
    s_min_m3 = CRIT_FRACTION * s_max_m3
    return {
        "S_max_m3": s_max_m3,
        "S_min_m3": s_min_m3,
        "flood_capacity_m3": float(data["flood_capacity_mcm"]) * 1e6,
        "crit_fraction": CRIT_FRACTION,
        "eta": ETA,
        "source": data.get("source", ""),
    }


def adjustment_levels(delta_u: float, L: int) -> np.ndarray:
    """Niveles discretos de ajuste a_l (multiplos enteros de delta_u, simetricos).

    L debe ser impar: L=3 -> {-1,0,1}*du, L=5 -> {-2..2}*du, L=7 -> {-3..3}*du.
    """
    if L % 2 == 0:
        raise ValueError(f"L debe ser impar (1-of-L simetrico); recibido L={L}")
    half = (L - 1) // 2
    ks = np.arange(-half, half + 1)  # enteros -half..+half
    return ks * float(delta_u)


def compute_delta_u(weekly_release_m3, T: int | None = None) -> float:
    """delta_u = 0.25 * mediana del release semanal sobre la ventana.

    Convencion DECIDIDA: por instancia. Si T se pasa, usa las primeras T semanas;
    si T es None, usa la ventana completa (valor de referencia).
    """
    r = np.asarray(weekly_release_m3, dtype=float)
    if T is not None:
        r = r[:T]
    return 0.25 * float(np.median(r))


def instance_params(weekly_release_m3, T: int, L: int) -> dict:
    """Parametros de una instancia: delta_u (por instancia), u_max, niveles.

    Tambien incluye ``delta_u_full_ref`` = delta_u sobre la ventana completa,
    solo como referencia de sanidad cross-instancia (no se optimiza con el).
    """
    delta_u = compute_delta_u(weekly_release_m3, T=T)
    half = (L - 1) // 2
    # u_max = mayor magnitud de ajuste = half*delta_u (=max|niveles|). Da 2*delta_u
    # para L=5 (oficial, spec ec. 11) y generaliza a L=3 (1*du) y L=7 (3*du),
    # consistente con |u(t)| <= u_max y con la definicion de niveles.
    return {
        "T": T,
        "L": L,
        "delta_u": delta_u,
        "u_max": half * delta_u,
        "levels": adjustment_levels(delta_u, L),
        "delta_u_full_ref": compute_delta_u(weekly_release_m3, T=None),
    }


def compute_weights(T: int, S_min_m3: float, u_max: float) -> dict:
    """Pesos oficiales del SRS (spec seccion 5). S_scale = S_min."""
    s_scale = S_min_m3
    w1 = 1.0 / ((T + 1) * s_scale ** 2)
    w2 = 0.1 / (T * u_max ** 2)
    w3 = 0.1 / ((T - 1) * (2.0 * u_max) ** 2)
    return {"w1": w1, "w2": w2, "w3": w3}
