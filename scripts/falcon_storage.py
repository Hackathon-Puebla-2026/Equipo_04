"""Modelo de almacenamiento y chequeo de restricciones (todo en m^3).

Dinamica simplificada del spec: S(t+1) = S(t) + DeltaS_obs(t) - u(t).
Restricciones oficiales: R(t)>=0, |u(t)|<=u_max, 0<=S(t)<=S_max, |sum u|<=B.
"""
from __future__ import annotations

import numpy as np


def simulate_storage(S0: float, deltaS, u) -> np.ndarray:
    """Trayectoria de almacenamiento. Devuelve S de longitud T+1 (S[0]=S0)."""
    deltaS = np.asarray(deltaS, dtype=float)
    u = np.asarray(u, dtype=float)
    T = len(u)
    if len(deltaS) < T:
        raise ValueError(f"deltaS (len {len(deltaS)}) mas corto que u (len {T})")
    S = np.empty(T + 1, dtype=float)
    S[0] = float(S0)
    for t in range(T):
        S[t + 1] = S[t] + deltaS[t] - u[t]
    return S


def compute_release(R_obs, u) -> np.ndarray:
    """Release optimizado R(t) = R_obs(t) + u(t)."""
    return np.asarray(R_obs, dtype=float) + np.asarray(u, dtype=float)


def check_constraints(R_obs, u, S, S_max: float, u_max: float, B: float,
                      tol: float = 1e-6) -> dict:
    """Verifica las 4 restricciones oficiales y reporta violaciones.

    - R(t) = R_obs+u >= 0
    - |u(t)| <= u_max
    - 0 <= S(t) <= S_max   (S incluye T+1 puntos)
    - |sum u| <= B         (B = eta * sum R_obs)
    """
    R_obs = np.asarray(R_obs, dtype=float)
    u = np.asarray(u, dtype=float)
    S = np.asarray(S, dtype=float)
    R = R_obs + u

    n_release_neg = int(np.sum(R < -tol))
    n_umax = int(np.sum(np.abs(u) > u_max + tol))
    n_storage_oob = int(np.sum((S < -tol) | (S > S_max + tol)))
    balance = float(abs(u.sum()))
    balance_ok = balance <= B + tol

    violations = {
        "release_negative": n_release_neg,
        "u_exceeds_umax": n_umax,
        "storage_out_of_bounds": n_storage_oob,
        "balance_excess": max(0.0, balance - B),
        "min_release": float(R.min()),
        "min_storage": float(S.min()),
        "max_storage": float(S.max()),
        "abs_balance": balance,
    }
    feasible = (n_release_neg == 0 and n_umax == 0 and n_storage_oob == 0 and balance_ok)
    return {"feasible": bool(feasible), "violations": violations}
