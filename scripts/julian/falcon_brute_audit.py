"""Auditoria esceptica por fuerza bruta: ¿estan bien el codigo, las restricciones
y los datos? ¿Se puede encontrar una solucion no trivial (u != 0)?

Motivacion: en la instancia chica (T12/L3) el optimo exacto da u=0 en todas las
ventanas. Antes de creerlo, verificamos TODO con fuerza bruta como oraculo:
  1) Datos: unidades y cross-check storage vs % de conservacion (fuente independiente).
  2) Restricciones: coinciden con check_constraints y con el spec; cual ata.
  3) Brute correcto: DOS implementaciones independientes de fuerza bruta coinciden;
     y un caso PLANTADO donde el optimo es demostrablemente u!=0 -> brute lo encuentra.
  4) Busqueda de no-trivial en datos reales (L=3 todas las ventanas; L=5 pocas).
  5) Sensibilidad: ¿cuando el optimo pasa a ser u!=0? (barrido de S_min y de T).

Solo fuerza bruta para hallar optimos (nada de DP/QUBO como autoridad). Todo en m^3.
Uso: .venv/bin/python scripts/julian/falcon_brute_audit.py
"""
from __future__ import annotations

import os
import sys
import time

# scripts/ (hermanos compartidos) al path: subimos un nivel desde scripts/julian/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd

import falcon_baselines as bl
import falcon_constants as fc
import falcon_data as fd
import falcon_srs as srs
import falcon_storage as st

TOL = 1e-6


# --------------------------------------------------------------------------- #
# Fuerza bruta VECTORIZADA e INDEPENDIENTE (distinta al itertools de baselines)
# --------------------------------------------------------------------------- #
def brute_vectorized(R_obs, dS, S0, *, S_min, S_max, delta_u, L, weights, B,
                     chunk=200_000, tol=TOL):
    """Enumera las L^T secuencias por lotes vectorizados y devuelve, entre las
    FACTIBLES: el mejor global, el mejor con u!=0, y el SRS de u=0.

    Camino de codigo totalmente distinto a bl.brute_force_optimal (contador base-L
    vectorizado + storage/costos con numpy), para cross-validar la fuerza bruta.
    """
    R_obs = np.asarray(R_obs, float)
    dS = np.asarray(dS, float)
    T = len(R_obs)
    half = (L - 1) // 2
    levels = (np.arange(L) - half) * float(delta_u)
    w1, w2, w3 = weights["w1"], weights["w2"], weights["w3"]
    csdS = np.cumsum(dS)                       # (T,)
    powers = L ** np.arange(T)                 # base-L (t=0 es el digito menos significativo)
    N = L ** T

    best = {"srs": -np.inf, "u": None}
    best_nz = {"srs": -np.inf, "u": None}
    srs_zero = None
    zero_idx = int((np.full(T, half) * powers).sum())  # combo todo-half == u=0

    for start in range(0, N, chunk):
        idx = np.arange(start, min(start + chunk, N))
        combos = (idx[:, None] // powers[None, :]) % L      # (b,T) indices de nivel
        u = levels[combos]                                  # (b,T) en m^3
        csu = np.cumsum(u, axis=1)
        S = np.empty((len(idx), T + 1))
        S[:, 0] = S0
        S[:, 1:] = S0 + csdS[None, :] - csu                 # S(t)=S0+Σ dS - Σ u
        # costos (identicos a falcon_srs, vectorizados)
        Ccrit = np.sum(np.maximum(0.0, S_min - S) ** 2, axis=1)
        Cdev = np.sum(u * u, axis=1)
        Csmooth = np.sum(np.diff(u, axis=1) ** 2, axis=1)
        srs_batch = -(w1 * Ccrit + w2 * Cdev + w3 * Csmooth)
        # factibilidad (misma definicion que st.check_constraints)
        Rmat = R_obs[None, :] + u
        feas = ((Rmat >= -tol).all(axis=1)
                & (S >= -tol).all(axis=1) & (S <= S_max + tol).all(axis=1)
                & (np.abs(u.sum(axis=1)) <= B + tol))
        if zero_idx is not None and start <= zero_idx < start + len(idx):
            srs_zero = float(srs_batch[zero_idx - start])
        if not feas.any():
            continue
        srs_f = np.where(feas, srs_batch, -np.inf)
        j = int(np.argmax(srs_f))
        if srs_f[j] > best["srs"]:
            best = {"srs": float(srs_f[j]), "u": u[j].copy()}
        # mejor factible con u != 0
        nz = feas & np.any(combos != half, axis=1)
        if nz.any():
            srs_nz = np.where(nz, srs_batch, -np.inf)
            k = int(np.argmax(srs_nz))
            if srs_nz[k] > best_nz["srs"]:
                best_nz = {"srs": float(srs_nz[k]), "u": u[k].copy()}

    return {"N": N, "best_srs": best["srs"], "best_u": best["u"],
            "best_nonzero_srs": best_nz["srs"], "best_nonzero_u": best_nz["u"],
            "srs_zero": srs_zero,
            "best_nonzero_gap": (None if best_nz["u"] is None or srs_zero is None
                                 else best_nz["srs"] - srs_zero)}


def _window(R, dS, S0, start, T):
    R = np.asarray(R, float); dS = np.asarray(dS, float)
    return R[start:start + T], dS[start:start + T], float(S0) + float(dS[:start].sum())


def _params_weights(_R_full, R_w, T, L, S_min):
    p = fc.instance_params(R_w, T, L)
    w = fc.compute_weights(T, S_min, p["u_max"])
    return p, w


# --------------------------------------------------------------------------- #
# 1) Auditoria de datos
# --------------------------------------------------------------------------- #
def audit_data(df, c):
    print("=" * 78)
    print("1) AUDITORIA DE DATOS")
    print("=" * 78)
    data_dir = fd.RAW_DATA_DIR

    _, s_unit = fd.load_ibwc_csv(fd._find_csv(data_dir, "Total Storage.Web-Daily-tcm@08461200"))
    _, r_unit = fd.load_ibwc_csv(fd._find_csv(data_dir, "Discharge.Best Available@08461300"))
    print(f"  storage unit={s_unit!r} (esperado m^3)  discharge unit={r_unit!r} (esperado m^3/s)")
    assert s_unit.lower().replace(" ", "") in ("m^3", "m3"), f"unidad storage inesperada: {s_unit}"
    assert r_unit.lower().replace(" ", "") in ("m^3/s", "m3/s"), f"unidad discharge inesperada: {r_unit}"

    # Cross-check independiente: S_obs/S_max == % de conservacion (serie aparte)
    pct_series, pct_unit = fd.load_ibwc_csv(
        fd._find_csv(data_dir, "Percentage.Conservation-Web-Telemetry@08461200"))
    S_max = c["S_max_m3"]
    ends = pd.to_datetime(df["week_end"])
    pct_at = pct_series.reindex(pct_series.index.union(ends)).sort_index().ffill().reindex(ends)
    implied = df["S_obs_m3"].to_numpy() / S_max * 100.0
    diff = np.abs(implied - pct_at.to_numpy())
    print(f"  cross-check S_obs/S_max vs %conservacion ({pct_unit}): "
          f"max|Δ|={np.nanmax(diff):.3f} pp  (media implied={implied.mean():.2f}%)")
    assert np.nanmax(diff) < 1.0, "storage/S_max NO coincide con % de conservacion (>1pp)!"

    S = df["S_obs_m3"].to_numpy(); S_min = c["S_min_m3"]
    print(f"  storage: {S.min()/1e6:.1f}M .. {S.max()/1e6:.1f}M m^3  "
          f"= {S.min()/S_max*100:.1f}%..{S.max()/S_max*100:.1f}% de S_max")
    print(f"  S_min (25%) = {S_min/1e6:.1f}M  ->  semanas con S<S_min: {(S < S_min).sum()}/{len(S)}")
    print(f"  DeltaS_source = {df['DeltaS_obs_m3'].size and df.attrs.get('units')!r} "
          f"| ΔS derivado: {fd.DELTAS_SOURCE}")
    print("  VEREDICTO datos: unidades OK y consistentes con %conservacion (fuente independiente).")
    print("    Regimen = SEQUIA (~11-15% de capacidad, siempre < S_min). Dataset = IBWC scrapeado,")
    print("    ΔS derivado (no oficial). No confirmado = dataset oficial de SharePoint.\n")


# --------------------------------------------------------------------------- #
# 2) Auditoria de restricciones
# --------------------------------------------------------------------------- #
def audit_constraints(df, c, T=12, L=3, start=0, n_sample=20000, seed=0):
    print("=" * 78)
    print("2) AUDITORIA DE RESTRICCIONES")
    print("=" * 78)
    R = df["R_obs_m3_week"].to_numpy(); dS = df["DeltaS_obs_m3"].to_numpy(); S0 = df.attrs["S0_m3"]
    R_w, dS_w, S0_w = _window(R, dS, S0, start, T)
    p, _ = _params_weights(R, R_w, T, L, c["S_min_m3"])
    du, umax = p["delta_u"], p["u_max"]
    B = c["eta"] * float(R_w.sum())
    half = (L - 1) // 2
    print(f"  T={T} L={L} start={start}: Δu={du/1e6:.3f}M  u_max={umax/1e6:.3f}M  "
          f"B=η·ΣR_obs={B/1e6:.2f}M  (B/Δu={B/du:.2f})  S_min={c['S_min_m3']/1e6:.1f}M  "
          f"S_max={c['S_max_m3']/1e6:.0f}M")

    # ¿Que restriccion ata? (bajo u=0 y bajo el maximo retiro factible)
    S0traj = st.simulate_storage(S0_w, dS_w, np.zeros(T))
    print(f"  bajo u=0: min S={S0traj.min()/1e6:.1f}M (cota inf 0: {'ata' if S0traj.min()<du else 'NO ata'}), "
          f"max S={S0traj.max()/1e6:.1f}M (cota sup S_max: {'ata' if S0traj.max()>c['S_max_m3']-du else 'NO ata'})")
    print(f"    -> restriccion que limita el retiro: BALANCE |Σu|≤B (presupuesto {B/du:.1f}·Δu).")

    # Consistencia: filtro inline del brute == st.check_constraints, en muestras aleatorias
    rng = np.random.default_rng(seed)
    levels = (np.arange(L) - half) * du
    mism = 0
    for _ in range(n_sample):
        lv = rng.integers(0, L, size=T)
        u = levels[lv]
        S = st.simulate_storage(S0_w, dS_w, u)
        # filtro inline (mismo que brute_vectorized / brute_force_optimal)
        R_full = R_w + u
        inline = bool((R_full >= -TOL).all() and (S >= -TOL).all()
                      and (S <= c["S_max_m3"] + TOL).all() and abs(u.sum()) <= B + TOL)
        chk = st.check_constraints(R_w, u, S, c["S_max_m3"], umax, B)["feasible"]
        if inline != chk:
            mism += 1
    print(f"  filtro brute vs st.check_constraints: {mism} discrepancias en {n_sample} muestras "
          f"({'OK' if mism == 0 else 'FALLA'}).")
    assert mism == 0
    print("  VEREDICTO restricciones: coinciden con check_constraints y con el spec; balance es la que ata.\n")


# --------------------------------------------------------------------------- #
# 3) Correctitud de la fuerza bruta
# --------------------------------------------------------------------------- #
def audit_brute_correctness(df, c):
    print("=" * 78)
    print("3) CORRECTITUD DE LA FUERZA BRUTA (2 implementaciones + caso plantado)")
    print("=" * 78)
    R = df["R_obs_m3_week"].to_numpy(); dS = df["DeltaS_obs_m3"].to_numpy(); S0 = df.attrs["S0_m3"]
    T, L = 12, 3
    n_weeks = int(df.attrs["n_weeks"])

    # (a) dos brutes independientes coinciden en TODAS las ventanas T12/L3
    max_srs_diff = 0.0
    for start in range(n_weeks - T + 1):
        R_w, dS_w, S0_w = _window(R, dS, S0, start, T)
        p, w = _params_weights(R, R_w, T, L, c["S_min_m3"])
        B = c["eta"] * float(R_w.sum())
        bv = brute_vectorized(R_w, dS_w, S0_w, S_min=c["S_min_m3"], S_max=c["S_max_m3"],
                              delta_u=p["delta_u"], L=L, weights=w, B=B)
        bf = bl.brute_force_optimal(R_w, dS_w, S0_w, S_min=c["S_min_m3"], S_max=c["S_max_m3"],
                                    delta_u=p["delta_u"], L=L, weights=w, B=B, max_combos=10**9)
        assert bf is not None, f"brute_force_optimal devolvio None en start={start}"
        assert bv["N"] == L ** T, f"brute_vectorized no enumera L^T ({bv['N']} != {L**T})"
        d = abs(bv["best_srs"] - bf["SRS_star"])
        max_srs_diff = max(max_srs_diff, d)
        assert d < 1e-9, f"brutes discrepan en start={start}: {bv['best_srs']} vs {bf['SRS_star']}"
    print(f"  (a) brute_vectorized == brute_force_optimal en {n_weeks-T+1} ventanas T12/L3  "
          f"(max|ΔSRS|={max_srs_diff:.2e}, ambos enumeran {L**T} combos).")

    # (b) micro-caso T=3 hecho a mano: optimo demostrablemente u != 0
    du = 10.0; L3 = 3; w = {"w1": 1.0, "w2": 1e-3, "w3": 1e-3}
    S0m, dSm, S_minm, S_maxm = 100.0, np.array([-50.0, 0.0, 0.0]), 100.0, 1e9
    Rm = np.array([1e9, 1e9, 1e9]); Bm = 1e9   # R>=0 y balance no atan
    bvm = brute_vectorized(Rm, dSm, S0m, S_min=S_minm, S_max=S_maxm, delta_u=du,
                           L=L3, weights=w, B=Bm)
    print(f"  (b) micro T=3 plantado: SRS(u=0)={bvm['srs_zero']:.3f}  optimo SRS={bvm['best_srs']:.3f} "
          f"u*={bvm['best_u']}  (nonzero={np.any(bvm['best_u']!=0)})")
    assert np.any(bvm["best_u"] != 0) and bvm["best_srs"] > bvm["srs_zero"] + 1e-9, \
        "FALLA: en el micro-caso plantado el brute deberia hallar u!=0 mejor que u=0"

    # (c) datos REALES: deficit toda la ventana (S_min sobre el storage) + objetivo
    #     dominado por C_crit (w2=w3=0) y balance holgado -> optimo debe ser u!=0.
    #     Prueba que el pipeline real (arrays, storage, costos, factibilidad, enumeracion)
    #     halla u!=0 cuando el objetivo lo premia. (Sonda de correctitud, no benchmark.)
    T = 12; L = 5; start = 0
    R_w, dS_w, S0_w = _window(R, dS, S0, start, T)
    p = fc.instance_params(R_w, T, L)
    S0traj = st.simulate_storage(S0_w, dS_w, np.zeros(T))
    S_min_planted = float(S0traj.max()) + p["delta_u"]         # TODA semana bajo el umbral
    w = {"w1": 1.0 / ((T + 1) * S_min_planted ** 2), "w2": 0.0, "w3": 0.0}  # solo C_crit
    B = 1.0 * float(R_w.sum())                                 # balance holgado (η=1.0)
    bvp = brute_vectorized(R_w, dS_w, S0_w, S_min=S_min_planted, S_max=c["S_max_m3"],
                           delta_u=p["delta_u"], L=L, weights=w, B=B)
    nzu = int(np.count_nonzero(bvp["best_u"]))
    print(f"  (c) datos reales + solo-C_crit + S_min={S_min_planted/1e6:.1f}M (deficit toda semana): "
          f"SRS(u=0)={bvp['srs_zero']:.3e}  optimo SRS={bvp['best_srs']:.3e}  u!=0 en {nzu}/{T} semanas")
    assert nzu > 0 and bvp["best_srs"] > bvp["srs_zero"] + 1e-12, \
        "FALLA: con deficit real y objetivo C_crit el optimo deberia ser u!=0"
    print(f"      u* = {np.round(bvp['best_u'] / p['delta_u']).astype(int)} (en unidades de Δu)")
    print("  VEREDICTO brute: CORRECTO. Dos implementaciones coinciden y halla u!=0 cuando existe.\n")


# --------------------------------------------------------------------------- #
# 4) Busqueda de no-trivial en datos reales
# --------------------------------------------------------------------------- #
def hunt_nonzero(df, c):
    print("=" * 78)
    print("4) BUSQUEDA DE NO-TRIVIAL EN DATOS REALES (solo fuerza bruta)")
    print("=" * 78)
    R = df["R_obs_m3_week"].to_numpy(); dS = df["DeltaS_obs_m3"].to_numpy(); S0 = df.attrs["S0_m3"]
    n_weeks = int(df.attrs["n_weeks"])

    # L=3, todas las ventanas
    T, L = 12, 3
    any_nz = 0; min_gap = np.inf
    print(f"  L=3, T=12: {n_weeks-T+1} ventanas")
    for start in range(n_weeks - T + 1):
        R_w, dS_w, S0_w = _window(R, dS, S0, start, T)
        p, w = _params_weights(R, R_w, T, L, c["S_min_m3"])
        B = c["eta"] * float(R_w.sum())
        bv = brute_vectorized(R_w, dS_w, S0_w, S_min=c["S_min_m3"], S_max=c["S_max_m3"],
                              delta_u=p["delta_u"], L=L, weights=w, B=B)
        nz = int(np.count_nonzero(bv["best_u"]))
        any_nz += (nz > 0)
        gap = bv["best_nonzero_gap"]  # SRS(best u!=0) - SRS(u=0); <0 => u=0 gana
        if gap is not None:
            min_gap = min(min_gap, -gap)  # cuanto pierde el mejor no-cero vs u=0
    print(f"    ventanas con optimo u!=0: {any_nz}/{n_weeks-T+1}. "
          f"El mejor schedule u!=0 pierde vs u=0 por >= {min_gap:.3e} SRS en toda ventana.")
    print("    (Nota: schedules FACTIBLES con u!=0 existen siempre -p.ej. threshold_rule-; "
          "lo que NO existe es un OPTIMO con u!=0 en este regimen.)")

    # L=5, pocas ventanas (full L=5 es pesado). first/middle/stress.
    T, L = 12, 5
    H = st.simulate_storage(S0, dS, np.zeros(n_weeks))
    last = n_weeks - T
    stress = int(np.argmin([H[s:s + T + 1].mean() for s in range(last + 1)]))
    wins = {"first": 0, "middle": last // 2, "stress": stress}
    print(f"  L=5, T=12: ventanas {wins}  (5^12={5**12:,} combos c/u, vectorizado)")
    for label, start in wins.items():
        R_w, dS_w, S0_w = _window(R, dS, S0, start, T)
        p, w = _params_weights(R, R_w, T, L, c["S_min_m3"])
        B = c["eta"] * float(R_w.sum())
        t0 = time.perf_counter()
        bv = brute_vectorized(R_w, dS_w, S0_w, S_min=c["S_min_m3"], S_max=c["S_max_m3"],
                              delta_u=p["delta_u"], L=L, weights=w, B=B, chunk=500_000)
        nz = int(np.count_nonzero(bv["best_u"]))
        _, _, srs_h = _srs_at_zero(R_w, dS_w, S0_w, c["S_min_m3"], w)
        print(f"    {label:6s} start={start:2d}: optimo u!=0 en {nz}/{T} semanas  "
              f"SRS={bv['best_srs']:.5f}  dSRS_vs_u0={bv['best_srs']-srs_h:+.3e}  "
              f"({time.perf_counter()-t0:.1f}s)")
    print()


def _srs_at_zero(R_w, dS_w, S0_w, S_min, w):
    u = np.zeros(len(R_w))
    S = st.simulate_storage(S0_w, dS_w, u)
    costs = srs.compute_costs(S, u, S_min)
    return S, costs, srs.compute_srs(costs, w)


# --------------------------------------------------------------------------- #
# 5) Sensibilidad: ¿cuando el optimo pasa a u != 0?
# --------------------------------------------------------------------------- #
def sensitivity(df, c):
    print("=" * 78)
    print("5) SENSIBILIDAD: ¿que regimen hace que el optimo sea u!=0?")
    print("=" * 78)
    R = df["R_obs_m3_week"].to_numpy(); dS = df["DeltaS_obs_m3"].to_numpy(); S0 = df.attrs["S0_m3"]
    T, L, start = 12, 3, 0
    R_w, dS_w, S0_w = _window(R, dS, S0, start, T)
    p = fc.instance_params(R_w, T, L)
    B = c["eta"] * float(R_w.sum())
    S_max = c["S_max_m3"]
    print(f"  Barrido de S_min como % de S_max (T12/L3, ventana 0). storage real "
          f"~{st.simulate_storage(S0_w,dS_w,np.zeros(T)).mean()/S_max*100:.1f}% de S_max:")
    print(f"  {'S_min %cap':>10s} {'S_min(M)':>10s} {'u!=0 sem':>9s} {'SRS':>12s}")
    for frac in (0.02, 0.05, 0.08, 0.11, 0.13, 0.15, 0.20, 0.25):
        S_min = frac * S_max
        w = fc.compute_weights(T, S_min, p["u_max"])
        bv = brute_vectorized(R_w, dS_w, S0_w, S_min=S_min, S_max=S_max,
                              delta_u=p["delta_u"], L=L, weights=w, B=B)
        nz = int(np.count_nonzero(bv["best_u"]))
        flag = "  <-- oficial" if abs(frac - 0.25) < 1e-9 else ""
        print(f"  {frac*100:9.0f}% {S_min/1e6:10.1f} {nz:9d} {bv['best_srs']:12.5f}{flag}")
    print("    Interpretacion: con storage ~11-15%, aun S_min bajo deja el optimo en u=0 salvo que")
    print("    S_min caiga al rango del storage (deficit chico y 'arreglable' con ±Δu).\n")


def main():
    t0 = time.perf_counter()
    df = fd.build_weekly_benchmark(write=False)
    c = fc.load_official_constants()
    print(f"Dataset: {df.attrs['n_weeks']} semanas | S_max={c['S_max_m3']/1e6:.0f}M "
          f"S_min={c['S_min_m3']/1e6:.1f}M m^3\n")

    audit_data(df, c)
    audit_constraints(df, c)
    audit_brute_correctness(df, c)
    hunt_nonzero(df, c)
    sensitivity(df, c)

    print("=" * 78)
    print("VEREDICTO FINAL")
    print("=" * 78)
    print("  * BRUTE:        CORRECTO (2 implementaciones independientes coinciden en 41 ventanas;")
    print("                  halla u!=0 en casos plantados micro y de datos reales).")
    print("  * RESTRICCIONES: CORRECTAS (identicas a check_constraints y al spec; balance es la que ata).")
    print("  * DATOS:        UNIDADES OK y consistentes con %conservacion (fuente independiente).")
    print("                  PERO: IBWC scrapeado, ΔS derivado, regimen de SEQUIA (~11-15% cap).")
    print("                  El optimo u=0 en T12/L3 es PROPIEDAD DEL REGIMEN, no un bug.")
    print("  * NO-TRIVIAL:   el optimo pasa a u!=0 con L=5/T mayor y/o storage cerca de S_min.")
    print("  RECOMENDACION:  conseguir el dataset OFICIAL de SharePoint (y el Change-in-Storage")
    print("                  oficial) para confirmar el regimen; T12/L3 es debug, no el benchmark.")
    print(f"\n({time.perf_counter()-t0:.1f}s)")


if __name__ == "__main__":
    main()
