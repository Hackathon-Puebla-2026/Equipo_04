"""Builder QUBO configurable para Falcon.

H(x) = xᵀ Q x + const, x en {0,1}ⁿ. Cada termino de costo / restriccion se escribe
como expresion lineal en bits y se eleva al cuadrado (`add_square_of_linear_expression`).
Convencion: Q simetrica; `xᵀQx = Σ_i Q_ii x_i + 2 Σ_{i<j} Q_ij x_i x_j` (x binario).

Defaults (FalconConfig): one-hot; storage_bounds="drop" (nunca ata); c_crit="soft"
(EXACTO aqui, S<S_min siempre); balance="slack" exacto de 1 slack; release="prohibit".
Con esos defaults, para bitstrings factibles `energia == J == -SRS`.
Todo en m^3; `normalize` reescala Q a O(1) para condicionamiento (meta["scale"]).
"""
from __future__ import annotations

import math

import numpy as np

import falcon_encodings as fe


def add_square_of_linear_expression(Q: np.ndarray, const: float, expr: dict,
                                    weight: float) -> float:
    """Suma weight*(c + Σ aᵢ xᵢ)² a Q (in place) y devuelve el nuevo const.

    (c+Σaᵢxᵢ)² = c² + Σ_i (2c aᵢ + aᵢ²) xᵢ + 2 Σ_{i<j} aᵢaⱼ xᵢxⱼ  (xᵢ²=xᵢ).
    Diagonal Q_ii += w(aᵢ²+2c aᵢ); off-diag Q_ij=Q_ji += w aᵢaⱼ; const += w c².
    """
    c = float(expr.get("constant", 0.0))
    items = [(int(i), float(a)) for i, a in expr.get("linear", {}).items() if a != 0.0]
    for i, ai in items:
        Q[i, i] += weight * (ai * ai + 2.0 * c * ai)
    for p in range(len(items)):
        i, ai = items[p]
        for q in range(p + 1, len(items)):
            j, aj = items[q]
            val = weight * ai * aj
            Q[i, j] += val
            Q[j, i] += val
    return const + weight * c * c


def qubo_energy(Q: np.ndarray, x, const: float = 0.0) -> float:
    """H(x) = xᵀQx + const (en las unidades en que este Q; ver meta['scale'])."""
    x = np.asarray(x, dtype=float)
    return float(x @ Q @ x) + float(const)


def _historical_Jscale(S0, deltaS, S_min, weights, T) -> float:
    """Magnitud tipica del objetivo (costo del historico u=0) para auto-escalar penalties."""
    deltaS = np.asarray(deltaS, dtype=float)
    H = np.empty(T + 1)
    H[0] = S0
    H[1:] = S0 + np.cumsum(deltaS[:T])
    Ccrit = float(np.sum(np.maximum(0.0, S_min - H) ** 2))
    return max(weights["w1"] * Ccrit, 1e-30)


def build_qubo(cfg, R_obs, deltaS, S0: float, *, S_min: float, delta_u: float,
               levels, weights: dict, B: float) -> tuple[np.ndarray, float, fe.VarIndex, dict]:
    """Arma (Q, const, var_index, meta) segun cfg. Solo one-hot por ahora."""
    if cfg.encoding != "onehot":
        raise NotImplementedError(f"encoding {cfg.encoding!r} aun no implementado (Fase 2)")
    T, L = cfg.T, cfg.L
    half = (L - 1) // 2
    levels = np.asarray(levels, dtype=float)
    R_obs = np.asarray(R_obs, dtype=float)
    w1, w2, w3 = weights["w1"], weights["w2"], weights["w3"]
    vi = fe.build_var_index(T, L)
    n_dec = vi.n

    # --- slack de balance exacto (1 var log-encoded) ---
    balance_slack_bits: list[int] = []
    M_cap = int(math.floor(B / delta_u))
    if cfg.balance == "slack":
        smax = 2 * M_cap
        nbits = max(1, int(math.ceil(math.log2(smax + 1)))) if smax > 0 else 0
        balance_slack_bits = [n_dec + r for r in range(nbits)]
    n = n_dec + len(balance_slack_bits)

    Q = np.zeros((n, n), dtype=float)
    const = 0.0

    # --- penalties (auto-escala desde J_scale si no se pasan) ---
    Jscale = _historical_Jscale(S0, deltaS, S_min, weights, T)
    P = {
        "onehot": cfg.penalties.get("P_onehot", 10.0 * Jscale),
        "R": cfg.penalties.get("P_R", 10.0 * Jscale),
        "bal_soft": cfg.penalties.get("P_bal_soft", 10.0 * Jscale / max(B * B, 1e-30)),
        "bal_slack": cfg.penalties.get("P_bal_slack", 10.0 * Jscale),
    }

    S_target = S_min  # soft C_crit con S_target=S_min == C_crit oficial (S<S_min siempre)

    # --- C_dev = w2 Σ_t u_t² ---
    if cfg.c_dev:
        for t in range(T):
            const = add_square_of_linear_expression(Q, const, fe.linear_expr_u(t, levels, vi), w2)

    # --- C_smooth = w3 Σ_{t=1}^{T-1} (u_t - u_{t-1})² ---
    if cfg.c_smooth:
        for t in range(1, T):
            expr = fe.add_exprs(fe.linear_expr_u(t, levels, vi),
                                fe.scale_expr(fe.linear_expr_u(t - 1, levels, vi), -1.0))
            const = add_square_of_linear_expression(Q, const, expr, w3)

    # --- C_crit (soft, exacto aqui) = w1 Σ_{t=0}^{T} (S_t - S_target)² ---
    if cfg.c_crit == "soft":
        for t in range(T + 1):
            se = fe.linear_expr_storage(t, S0, deltaS, levels, vi)
            se = {"constant": se["constant"] - S_target, "linear": se["linear"]}
            const = add_square_of_linear_expression(Q, const, se, w1)
    else:
        raise NotImplementedError("c_crit='deficit_slack' (Opcion B) es Fase 2")

    # --- one-hot: P Σ_t (Σ_l x_{t,l} - 1)² ---
    if cfg.onehot == "penalty":
        for t in range(T):
            expr = {"constant": -1.0, "linear": {vi.idx(t, l): 1.0 for l in range(L)}}
            const = add_square_of_linear_expression(Q, const, expr, P["onehot"])
    # xy_mixer: sin penalizacion (se maneja en el mixer del QAOA, Fase 3)

    # --- release no negativo: prohibir niveles con R_obs+a_l<0 (penalizacion lineal) ---
    n_forbidden = 0
    if cfg.release_nonneg == "prohibit":
        for t in range(T):
            for l in range(L):
                if R_obs[t] + levels[l] < -1e-9:
                    Q[vi.idx(t, l), vi.idx(t, l)] += P["R"]
                    n_forbidden += 1

    # --- balance |Σu| ≤ B ---
    if cfg.balance == "soft":
        expr = {"constant": 0.0, "linear": {}}
        for t in range(T):
            expr = fe.add_exprs(expr, fe.linear_expr_u(t, levels, vi))
        const = add_square_of_linear_expression(Q, const, expr, P["bal_soft"])
    elif cfg.balance == "slack":
        # M = Σ_t Σ_l k_l x_{t,l} (k_l = l-half, entero). M + s - M_cap = 0.
        lin = {}
        for t in range(T):
            for l in range(L):
                lin[vi.idx(t, l)] = float(l - half)
        for r, var in enumerate(balance_slack_bits):
            lin[var] = float(2 ** r)
        expr = {"constant": -float(M_cap), "linear": lin}
        const = add_square_of_linear_expression(Q, const, expr, P["bal_slack"])

    # --- normalizacion interna (condicionamiento) ---
    scale = 1.0
    if cfg.normalize == "maxabs":
        m = float(np.max(np.abs(Q)))
        scale = m if m > 0 else 1.0
    elif cfg.normalize == "delta_u2":
        scale = float(delta_u * delta_u)
    if scale != 1.0:
        Q = Q / scale
        const = const / scale

    meta = {
        "n_qubits": n, "n_decision": n_dec, "n_balance_slack": len(balance_slack_bits),
        "encoding": cfg.encoding, "scale": scale, "S_target": S_target,
        "M_cap": M_cap, "n_forbidden_levels": n_forbidden, "penalties": P,
        "levels": levels.tolist(),
    }
    return Q, const, vi, meta


def to_quadratic_program(Q: np.ndarray, name: str = "falcon_qubo"):
    """Q -> Qiskit QuadraticProgram (diag->lineal, off-diag i<j -> 2*Q[i,j]).

    Import perezoso de qiskit_optimization (solo Fase 3; no requerido en Fase 1).
    """
    from qiskit_optimization import QuadraticProgram
    n = Q.shape[0]
    qp = QuadraticProgram(name)
    for i in range(n):
        qp.binary_var(f"x{i}")
    linear = {f"x{i}": float(Q[i, i]) for i in range(n)}
    quadratic = {}
    for i in range(n):
        for j in range(i + 1, n):
            val = float(Q[i, j] + Q[j, i])
            if val != 0.0:
                quadratic[(f"x{i}", f"x{j}")] = val
    qp.minimize(linear=linear, quadratic=quadratic)
    return qp
