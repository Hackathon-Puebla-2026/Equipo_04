"""Abstraccion de encoding para el builder QUBO (empezando por one-hot).

Los terminos de costo del QUBO consumen solo "expresiones lineales en los bits":
un dict `{"constant": c, "linear": {var_index: coef}}`. Estas 3 funciones
(build_var_index, linear_expr_u, linear_expr_storage) son lo UNICO que cambia al
cambiar de encoding (one-hot -> domain-wall / binary), no los terminos de costo.

One-hot: x_{t,l}=1 <=> u(t)=a_l. n = T*L. u_t = sum_l a_l x_{t,l}.
Storage: S_t = H_t - sum_{k<t} u_k, con H_t = S0 + sum_{j<t} DeltaS_obs[j] (lineal en bits).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class VarIndex:
    """Mapea (semana t, nivel l) -> indice global de bit, para one-hot."""
    T: int
    L: int

    @property
    def n(self) -> int:
        return self.T * self.L

    def idx(self, t: int, l: int) -> int:
        return t * self.L + l

    def encode_levels(self, level_indices) -> np.ndarray:
        """Secuencia de indices de nivel (len T) -> vector de bits one-hot (len n)."""
        x = np.zeros(self.n, dtype=float)
        for t, l in enumerate(level_indices):
            x[self.idx(t, int(l))] = 1.0
        return x

    def decode_levels(self, x) -> np.ndarray:
        """Vector de bits -> indice de nivel por semana (argmax tolerante a one-hot imperfecto)."""
        x = np.asarray(x, dtype=float).reshape(self.T, self.L)
        return np.argmax(x, axis=1)

    def decode_u(self, x, levels) -> np.ndarray:
        """Vector de bits -> u(t) en m^3 (usa los valores de nivel a_l)."""
        levels = np.asarray(levels, dtype=float)
        return levels[self.decode_levels(x)]

    def is_onehot(self, x, tol: float = 1e-9) -> bool:
        """True si cada semana tiene exactamente un bit encendido."""
        xr = np.asarray(x, dtype=float).reshape(self.T, self.L)
        return bool(np.all(np.abs(xr.sum(axis=1) - 1.0) < tol))


def build_var_index(T: int, L: int) -> VarIndex:
    return VarIndex(T=T, L=L)


def linear_expr_u(t: int, levels, vi: VarIndex) -> dict:
    """u(t) = sum_l a_l x_{t,l} como expresion lineal en bits."""
    levels = np.asarray(levels, dtype=float)
    return {"constant": 0.0,
            "linear": {vi.idx(t, l): float(levels[l]) for l in range(vi.L)}}


def linear_expr_storage(t: int, S0: float, deltaS, levels, vi: VarIndex) -> dict:
    """S_t = H_t - sum_{k<t} sum_l a_l x_{k,l}  (lineal en bits). t en 0..T."""
    deltaS = np.asarray(deltaS, dtype=float)
    levels = np.asarray(levels, dtype=float)
    H_t = float(S0) + float(np.sum(deltaS[:t]))          # t=0 -> S0
    linear: dict[int, float] = {}
    for k in range(t):                                    # semanas previas
        for l in range(vi.L):
            linear[vi.idx(k, l)] = linear.get(vi.idx(k, l), 0.0) - float(levels[l])
    return {"constant": H_t, "linear": linear}


def add_exprs(*exprs) -> dict:
    """Suma expresiones lineales (p.ej. u(t) - u(t-1) para C_smooth)."""
    out: dict = {"constant": 0.0, "linear": {}}
    for e in exprs:
        out["constant"] += e.get("constant", 0.0)
        for var, coef in e.get("linear", {}).items():
            out["linear"][var] = out["linear"].get(var, 0.0) + coef
    return out


def scale_expr(expr: dict, factor: float) -> dict:
    """Multiplica una expresion lineal por un escalar (p.ej. -1 para restar)."""
    return {"constant": expr.get("constant", 0.0) * factor,
            "linear": {v: c * factor for v, c in expr.get("linear", {}).items()}}
