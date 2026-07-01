"""Abstraccion de encoding para el builder QUBO (one-hot y binary).

Los terminos de costo del QUBO consumen solo "expresiones lineales en los bits":
un dict `{"constant": c, "linear": {var_index: coef}}`. Cambiar de encoding solo
toca la capa de encoding (index, expresion de u, penalties de validez / niveles
prohibidos, expresion de balance), no los terminos de costo (`build_qubo`).

- One-hot: x_{t,l}=1 <=> u(t)=a_l. n = T*L. u_t = sum_l a_l x_{t,l}.
- Binary: l en {0..L-1} codificado en b=ceil(log2 L) bits/semana (n=T*b). Requiere
  niveles equiespaciados (a_l = a_0 + l*step): u_t = a_0 + step*sum_b 2^b y_{t,b}
  (lineal en bits). Codewords l>=L son invalidos -> se penalizan (producto de bits).
  Exacto en QUBO solo si b<=2 (L<=4); L=5 binary requiere manejo aparte (no cabe en
  statevector de todos modos). Storage se deriva de u (agnostico al encoding).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

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


@dataclass
class BinaryVarIndex:
    """Codifica el nivel l en {0..L-1} con b=ceil(log2 L) bits/semana. n=T*b.

    idx(t,r) = bit r (2^r) de la semana t. Requiere b<=2 (L<=4) para penalties
    exactas en QUBO (producto de bits de grado <=2).
    """
    T: int
    L: int
    bits: int = field(init=False)

    def __post_init__(self):
        self.bits = max(1, int(math.ceil(math.log2(self.L))))

    @property
    def n(self) -> int:
        return self.T * self.bits

    def idx(self, t: int, r: int) -> int:
        return t * self.bits + r

    def num_codewords(self) -> int:
        return 1 << self.bits

    def invalid_levels(self):
        return [l for l in range(self.num_codewords()) if l >= self.L]

    def bit_pattern(self, l: int):
        """Bits (r=0..b-1) del nivel l (little-endian)."""
        return [(int(l) >> r) & 1 for r in range(self.bits)]

    def encode_levels(self, level_indices) -> np.ndarray:
        x = np.zeros(self.n, dtype=float)
        for t, l in enumerate(level_indices):
            for r, bit in enumerate(self.bit_pattern(int(l))):
                x[self.idx(t, r)] = float(bit)
        return x

    def decode_levels(self, x) -> np.ndarray:
        """Vector de bits -> nivel por semana (puede ser >=L si el codeword es invalido)."""
        x = np.asarray(x, dtype=float)
        lv = np.empty(self.T, dtype=int)
        for t in range(self.T):
            lv[t] = sum((1 if x[self.idx(t, r)] >= 0.5 else 0) << r for r in range(self.bits))
        return lv

    def is_valid_week(self, l: int) -> bool:
        return 0 <= int(l) < self.L


def build_var_index(T: int, L: int) -> VarIndex:
    return VarIndex(T=T, L=L)


def build_binary_var_index(T: int, L: int) -> BinaryVarIndex:
    vi = BinaryVarIndex(T=T, L=L)
    if vi.bits > 2:
        raise NotImplementedError(
            f"binary encoding exacto solo con b<=2 (L<=4); L={L} necesita {vi.bits} bits")
    return vi


def linear_expr_u(t: int, levels, vi) -> dict:
    """u(t) como expresion lineal en bits (one-hot o binary)."""
    levels = np.asarray(levels, dtype=float)
    if isinstance(vi, BinaryVarIndex):
        a0 = float(levels[0])
        step = float(levels[1] - levels[0]) if vi.L > 1 else 0.0
        return {"constant": a0,
                "linear": {vi.idx(t, r): step * float(1 << r) for r in range(vi.bits)}}
    return {"constant": 0.0,
            "linear": {vi.idx(t, l): float(levels[l]) for l in range(vi.L)}}


def linear_expr_storage(t: int, S0: float, deltaS, levels, vi) -> dict:
    """S_t = H_t - sum_{k<t} u_k  (lineal en bits). Agnostico al encoding (usa linear_expr_u)."""
    deltaS = np.asarray(deltaS, dtype=float)
    H_t = float(S0) + float(np.sum(deltaS[:t]))          # t=0 -> S0
    expr = {"constant": H_t, "linear": {}}
    for k in range(t):                                    # semanas previas
        expr = add_exprs(expr, scale_expr(linear_expr_u(k, levels, vi), -1.0))
    return expr


def balance_M_expr(vi, half: int) -> dict:
    """Expresion lineal de M = sum_t (l_t - half) (entero), para el slack de balance.

    one-hot: M = sum_{t,l} (l-half) x_{t,l}.  binary: M = sum_{t,r} 2^r y_{t,r} - T*half.
    """
    if isinstance(vi, BinaryVarIndex):
        linear = {vi.idx(t, r): float(1 << r) for t in range(vi.T) for r in range(vi.bits)}
        return {"constant": -float(vi.T * half), "linear": linear}
    linear = {vi.idx(t, l): float(l - half) for t in range(vi.T) for l in range(vi.L)}
    return {"constant": 0.0, "linear": linear}


def add_codeword_penalty(Q, const: float, vi: BinaryVarIndex, t: int, l: int,
                         P: float) -> float:
    """Suma P * [semana t toma el codeword l] a Q (binary). Devuelve el nuevo const.

    Indicador = producto de factores (y_r si el bit r de l es 1, si no 1-y_r).
    Exacto en QUBO para bits<=2 (grado <=2). y_r²=y_r (diagonal = coef lineal).
    """
    pat = vi.bit_pattern(l)
    # factor r: (c_r, var_r, a_r) representa c_r + a_r * y_r
    factors = []
    for r in range(vi.bits):
        i = vi.idx(t, r)
        if pat[r] == 1:
            factors.append((0.0, i, 1.0))     # y_r
        else:
            factors.append((1.0, i, -1.0))    # 1 - y_r
    if len(factors) == 1:
        c0, i0, a0 = factors[0]
        Q[i0, i0] += P * a0
        return const + P * c0
    # len == 2
    (c0, i0, a0), (c1, i1, a1) = factors
    const += P * c0 * c1
    Q[i0, i0] += P * a0 * c1
    Q[i1, i1] += P * a1 * c0
    v = P * a0 * a1 / 2.0
    Q[i0, i1] += v
    Q[i1, i0] += v
    return const


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
