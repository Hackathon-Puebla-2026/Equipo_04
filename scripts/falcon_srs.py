"""Storage Resilience Score (SRS) y sus terminos de costo (spec seccion 3).

SRS = -(w1*Ccrit + w2*Cdev + w3*Csmooth). Maximizar SRS = minimizar el costo.
Indices del spec:
  Ccrit   = sum_{t=0}^{T}   max(0, S_min - S(t))^2     (T+1 terminos; S tiene T+1 puntos)
  Cdev    = sum_{t=0}^{T-1}  u(t)^2
  Csmooth = sum_{t=1}^{T-1} (u(t) - u(t-1))^2
"""
from __future__ import annotations

import numpy as np


def compute_costs(S, u, S_min: float) -> dict:
    """Devuelve {Ccrit, Cdev, Csmooth} a partir de S (len T+1) y u (len T)."""
    S = np.asarray(S, dtype=float)
    u = np.asarray(u, dtype=float)
    deficit = np.maximum(0.0, S_min - S)          # solo cuenta por debajo de S_min
    Ccrit = float(np.sum(deficit ** 2))
    Cdev = float(np.sum(u ** 2))
    Csmooth = float(np.sum(np.diff(u) ** 2)) if len(u) >= 2 else 0.0
    return {"Ccrit": Ccrit, "Cdev": Cdev, "Csmooth": Csmooth}


def compute_srs(costs: dict, weights: dict) -> float:
    """SRS = -(w1*Ccrit + w2*Cdev + w3*Csmooth)."""
    return -(weights["w1"] * costs["Ccrit"]
             + weights["w2"] * costs["Cdev"]
             + weights["w3"] * costs["Csmooth"])
