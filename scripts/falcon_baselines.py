"""Baselines clasicos: historico, regla de umbral y DP exacto (ground truth).

El DP es exacto porque, al discretizar u en niveles enteros de Delta_u, el
almacenamiento vive en un lattice entero: S(t) = H_t - Delta_u * C_t, con
C_t = suma entera acumulada de los niveles. Estado (t, C_t, k_prev); k_prev solo
hace falta para C_smooth. Complejidad O(T^2 * L^2). Ver spec seccion 7.

Todo en m^3. Reusa falcon_storage / falcon_srs para reconstruir y reportar.
"""
from __future__ import annotations

import itertools

import numpy as np

import falcon_srs as srs
import falcon_storage as st


def historical(T: int) -> np.ndarray:
    """u(t) = 0 para todo t (operacion historica)."""
    return np.zeros(int(T), dtype=float)


def threshold_rule(R_obs, deltaS, S0: float, S_min: float, delta_u: float,
                   clamp_release: bool = False, B: float | None = None,
                   tol: float = 1e-6) -> np.ndarray:
    """Regla de conservacion por umbral del spec (ecs. 18-19):

        u_rule(t) = -delta_u  si S_rule(t) <  S_min
                  =  0        si S_rule(t) >= S_min

    S_rule(t) es la trayectoria bajo la propia regla (auto-consistente). Por
    defecto es FIEL al PDF (sin clamp), por lo que puede producir R(t)<0 en
    semanas de muy bajo caudal. Con ``clamp_release=True`` se aplica el limite
    fisico ``u >= -R_obs(t)`` (variante factible, opt-in).

    Si se pasa ``B`` (presupuesto de balance = eta*sum(R_obs)), se activa la
    variante **balanceada** (de ivan): detiene la reduccion cuando aplicar otra
    ``-delta_u`` haria ``|sum u| > B``, quedando factible en balance.
    """
    R_obs = np.asarray(R_obs, dtype=float)
    deltaS = np.asarray(deltaS, dtype=float)
    T = len(R_obs)
    u = np.zeros(T, dtype=float)
    S = float(S0)
    running = 0.0
    for t in range(T):
        ut = -delta_u if S < S_min else 0.0
        if clamp_release:
            ut = max(ut, -R_obs[t])      # mantener R(t) >= 0 (opt-in)
        if B is not None and abs(running + ut) > B + tol:
            ut = 0.0                     # balanceada: no reducir si rompe el balance
        u[t] = ut
        running += ut
        S = S + deltaS[t] - ut
    return u


def dp_optimal(R_obs, deltaS, S0: float, *, S_min: float, S_max: float,
               delta_u: float, L: int, weights: dict, B: float,
               tol: float = 1e-6) -> dict:
    """Optimo exacto del problema discretizado por DP sobre el lattice entero.

    Devuelve {u_star, SRS_star, costs, feasible}. Maximiza SRS = minimiza
    J = w1*Ccrit + w2*Cdev + w3*Csmooth, con factibilidad dura (R>=0, cotas de
    storage, balance |sum u| <= B).
    """
    R_obs = np.asarray(R_obs, dtype=float)
    deltaS = np.asarray(deltaS, dtype=float)
    T = len(R_obs)
    half = (L - 1) // 2
    ks = list(range(-half, half + 1))
    w1, w2, w3 = weights["w1"], weights["w2"], weights["w3"]

    # H_t = S0 + sum_{j<t} deltaS[j],  t = 0..T
    H = np.empty(T + 1, dtype=float)
    H[0] = S0
    H[1:] = S0 + np.cumsum(deltaS)
    C_cap = B / delta_u  # |C_T| <= C_cap (balance)

    def crit(t: int, C: int) -> float:
        S_t = H[t] - delta_u * C
        if S_t < -tol or S_t > S_max + tol:
            return float("inf")          # fuera de cotas fisicas -> estado podado
        deficit = max(0.0, S_min - S_t)
        return w1 * deficit * deficit

    # stages[t]: dict (C, k_prev) -> (cost_acumulado, parent_key, k_chosen)
    c0 = crit(0, 0)
    stages = [{(0, None): (c0, None, None)}]

    for t in range(T):
        cur = stages[t]
        nxt: dict = {}
        for (C, k_prev), (g, _, _) in cur.items():
            if g == float("inf"):
                continue
            for k in ks:
                if R_obs[t] + k * delta_u < -tol:    # R(t) >= 0
                    continue
                C2 = C + k
                cadd = crit(t + 1, C2)
                if cadd == float("inf"):
                    continue
                add = w2 * (delta_u * k) ** 2 + cadd
                if t >= 1:
                    add += w3 * (delta_u * (k - k_prev)) ** 2
                g2 = g + add
                key = (C2, k)
                prev = nxt.get(key)
                if prev is None or g2 < prev[0]:
                    nxt[key] = (g2, (C, k_prev), k)
        stages.append(nxt)

    # Terminal: aplicar balance |Delta_u * C_T| <= B
    terminal = stages[T]
    feasible = True
    candidates = {key: v for key, v in terminal.items()
                  if abs(key[0]) <= C_cap + tol}
    if not candidates:                    # ningun terminal respeta el balance
        feasible = False
        candidates = terminal
    if not candidates:
        raise RuntimeError("DP sin estados terminales (instancia infactible).")

    best_key = min(candidates, key=lambda kk: candidates[kk][0])

    # Reconstruir secuencia de niveles k via backpointers
    ks_path = []
    key = best_key
    for t in range(T, 0, -1):
        _, parent, k = stages[t][key]
        ks_path.append(k)
        key = parent
    ks_path.reverse()
    u_star = np.array(ks_path, dtype=float) * delta_u

    # Reconstruir costos/SRS con los modulos compartidos (cross-check)
    S = st.simulate_storage(S0, deltaS, u_star)
    costs = srs.compute_costs(S, u_star, S_min)
    srs_val = srs.compute_srs(costs, weights)
    return {"u_star": u_star, "SRS_star": srs_val, "costs": costs, "feasible": bool(feasible)}


def brute_force_optimal(R_obs, deltaS, S0: float, *, S_min: float, S_max: float,
                        delta_u: float, L: int, weights: dict, B: float,
                        max_combos: int = 600_000, tol: float = 1e-6) -> dict | None:
    """Optimo por enumeracion exhaustiva de las L^T secuencias (arbitro neutral).

    Filtra factibilidad dura (R>=0, 0<=S<=S_max, |sum u|<=B) y maximiza SRS.
    Devuelve la misma forma que `dp_optimal`, o `None` si `L**T > max_combos`
    (no enumerable). Sirve para validar el DP en instancias chicas.
    """
    R_obs = np.asarray(R_obs, dtype=float)
    deltaS = np.asarray(deltaS, dtype=float)
    T = len(R_obs)
    if L ** T > max_combos:
        return None
    half = (L - 1) // 2
    ks = range(-half, half + 1)

    best_srs, best_u = -np.inf, None
    best_infeas_srs, best_infeas_u = -np.inf, None
    for combo in itertools.product(ks, repeat=T):
        u = np.array(combo, dtype=float) * delta_u
        S = st.simulate_storage(S0, deltaS, u)
        R = R_obs + u
        feasible = (R >= -tol).all() and (S >= -tol).all() and (S <= S_max + tol).all() \
            and abs(u.sum()) <= B + tol
        val = srs.compute_srs(srs.compute_costs(S, u, S_min), weights)
        if feasible:
            if val > best_srs:
                best_srs, best_u = val, u
        elif val > best_infeas_srs:
            best_infeas_srs, best_infeas_u = val, u

    if best_u is not None:
        u_star, srs_val, feasible = best_u, best_srs, True
    else:                                   # no hay factible: reporta la mejor infactible
        u_star, srs_val, feasible = best_infeas_u, best_infeas_srs, False

    S = st.simulate_storage(S0, deltaS, u_star)
    costs = srs.compute_costs(S, u_star, S_min)
    return {"u_star": u_star, "SRS_star": srs_val, "costs": costs, "feasible": bool(feasible)}
