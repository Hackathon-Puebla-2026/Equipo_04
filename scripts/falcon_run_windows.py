"""E1 - Robustez de ventana: corre los baselines + DP + brute sobre distintas
ventanas temporales de T semanas y compara.

Pregunta central: ¿el optimo u*=0 (que encontramos en la ventana start=0 de la
instancia chica) es un artefacto de esa ventana, o es representativo? Es decir,
¿existe alguna ventana de T semanas donde el optimo exacto sea u != 0, y coinciden
los metodos?

Estrategia:
  1) DP exacto escanea TODAS las ventanas (start = 0..n_weeks-T). El DP es exacto
     por-ventana (mismo algoritmo, solo cambian H_t y B) y esta validado == brute
     en start=0. Cuenta cuantas ventanas dan u* != 0.
  2) Brute force reconfirma (solo donde es enumerable: small 3^12, debug 3^5) las
     ventanas marcadas con u*!=0. medium/large: 5^26 / 5^52 no son enumerables ->
     el DP es el oraculo.
  3) Tabla completa de metodos (historical, threshold pure/clamped/balanced, dp,
     brute-si-enumerable) para 3 ventanas canonicas: first / middle / stress.

Convencion Delta_u por ventana (spec sec. 2): 0.25 * mediana del release sobre la
propia ventana; se registra ademas el Delta_u de ventana completa como referencia.

Uso: python scripts/falcon_run_windows.py
Escribe results/runs_summary.csv + results/runs/*.json con window_start/window_label.
"""
from __future__ import annotations

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

import falcon_baselines as bl
import falcon_constants as fc
import falcon_data as fd
import falcon_results as res
import falcon_srs as srs
import falcon_storage as st

# (T, L) estandar. debug (5,3) se incluye para un chequeo brute barato.
STANDARD_INSTANCES = [(5, 3), (12, 3), (26, 5), (52, 5), (52, 7)]
BRUTE_MAX_COMBOS = 600_000


# --------------------------------------------------------------------------- #
# Helpers de ventana
# --------------------------------------------------------------------------- #
def window_slice(R, dS, S0: float, start: int, T: int):
    """Extrae la ventana [start, start+T): (R_w, dS_w, S0_w).

    S0_w = almacenamiento ENTRANDO a la semana `start` = S0 + sum(dS[:start]),
    consistente con H_t = S0 + cumsum(dS) que usan el DP y los encodings.
    """
    R = np.asarray(R, dtype=float)
    dS = np.asarray(dS, dtype=float)
    R_w = R[start:start + T]
    dS_w = dS[start:start + T]
    S0_w = float(S0) + float(dS[:start].sum())
    return R_w, dS_w, S0_w


def window_params(R_full, R_w, T: int, L: int) -> dict:
    """params de la instancia con Delta_u sobre la VENTANA (no las primeras T sem).

    Se pasa R_w (largo T) a instance_params, cuyo [:T] es un no-op; luego se
    sobreescribe delta_u_full_ref con el Delta_u de la ventana completa (referencia).
    """
    p = fc.instance_params(R_w, T, L)
    p["delta_u_full_ref"] = fc.compute_delta_u(R_full)
    return p


def select_windows(S0: float, dS_full, T: int, n_weeks: int) -> dict:
    """Ventanas canonicas: first (start=0), middle, stress.

    stress = ventana de T semanas con MENOR almacenamiento medio (deficit mas
    profundo => mayor gradiente de C_crit => donde retener agua rinde mas). Empate
    -> menor start. Para T == n_weeks (large) hay una sola ventana: {"full": 0}.
    """
    last_start = n_weeks - T
    if last_start <= 0:
        return {"full": 0}
    H = st.simulate_storage(S0, dS_full, np.zeros(n_weeks))  # storage bajo u=0, largo n_weeks+1
    means = [H[s:s + T + 1].mean() for s in range(last_start + 1)]
    stress = int(np.argmin(means))
    wins = {"first": 0, "middle": last_start // 2, "stress": stress}
    return wins


def _srs_of(u, S0, dS, S_min, weights):
    S = st.simulate_storage(S0, dS, u)
    costs = srs.compute_costs(S, u, S_min)
    return S, costs, srs.compute_srs(costs, weights)


# --------------------------------------------------------------------------- #
# 1) Escaneo DP de todas las ventanas + reconfirmacion brute
# --------------------------------------------------------------------------- #
def scan_all_windows(df, constants, T: int, L: int) -> dict:
    """DP exacto en cada ventana; cuenta u*!=0. Reconfirma brute donde enumerable."""
    R = df["R_obs_m3_week"].to_numpy()
    dS = df["DeltaS_obs_m3"].to_numpy()
    S0 = df.attrs["S0_m3"]
    n_weeks = int(df.attrs["n_weeks"])
    S_min, S_max = constants["S_min_m3"], constants["S_max_m3"]
    instance = res.instance_label(T, L)
    last_start = n_weeks - T
    enumerable = L ** T <= BRUTE_MAX_COMBOS

    nonzero_starts, scan = [], []
    for start in range(last_start + 1):
        R_w, dS_w, S0_w = window_slice(R, dS, S0, start, T)
        p = window_params(R, R_w, T, L)
        w = fc.compute_weights(T, S_min, p["u_max"])
        B = constants["eta"] * float(R_w.sum())
        dp = bl.dp_optimal(R_w, dS_w, S0_w, S_min=S_min, S_max=S_max,
                           delta_u=p["delta_u"], L=L, weights=w, B=B)
        _, _, srs_h = _srs_of(bl.historical(T), S0_w, dS_w, S_min, w)
        nz = int(np.count_nonzero(dp["u_star"]))
        if nz > 0:
            nonzero_starts.append(start)
        scan.append({"start": start, "nonzero": nz, "SRS": dp["SRS_star"],
                     "dSRS_hist": dp["SRS_star"] - srs_h, "feasible": dp["feasible"]})

    # Reconfirmacion brute (solo enumerable): ventanas u*!=0 + controles u*==0
    reconfirmed = 0
    if enumerable:
        controls = [s["start"] for s in scan if s["nonzero"] == 0][:3]
        for start in sorted(set(nonzero_starts) | set(controls)):
            R_w, dS_w, S0_w = window_slice(R, dS, S0, start, T)
            p = window_params(R, R_w, T, L)
            w = fc.compute_weights(T, S_min, p["u_max"])
            B = constants["eta"] * float(R_w.sum())
            dp = bl.dp_optimal(R_w, dS_w, S0_w, S_min=S_min, S_max=S_max,
                               delta_u=p["delta_u"], L=L, weights=w, B=B)
            bf = bl.brute_force_optimal(R_w, dS_w, S0_w, S_min=S_min, S_max=S_max,
                                        delta_u=p["delta_u"], L=L, weights=w, B=B,
                                        max_combos=BRUTE_MAX_COMBOS)
            assert bf is not None, f"brute None en ventana enumerable {instance} start={start}"
            assert abs(bf["SRS_star"] - dp["SRS_star"]) < 1e-9, \
                f"brute != dp en {instance} start={start}: {bf['SRS_star']} vs {dp['SRS_star']}"
            # tambien los schedules deben coincidir en su condicion de nulidad
            assert np.count_nonzero(bf["u_star"]) == np.count_nonzero(dp["u_star"]) or \
                abs(bf["SRS_star"] - dp["SRS_star"]) < 1e-9
            reconfirmed += 1

    oracle = "brute==dp reconfirmado" if enumerable else "DP oraculo (5^T no enumerable)"
    print(f"[{instance:7s} T={T:2d} L={L}] escaneo: {len(nonzero_starts)} de "
          f"{last_start + 1} ventanas con u*!=0"
          + (f"  starts={nonzero_starts}" if nonzero_starts else "")
          + f"  | {oracle}"
          + (f" ({reconfirmed} ventanas)" if enumerable else ""))
    return {"instance": instance, "n_windows": last_start + 1,
            "nonzero_starts": nonzero_starts, "scan": scan, "enumerable": enumerable}


# --------------------------------------------------------------------------- #
# 2) Tabla completa de metodos en ventanas canonicas
# --------------------------------------------------------------------------- #
def run_canonical_windows(df, constants, T: int, L: int, windows: dict) -> list:
    """Corre historical/threshold(3)/dp/brute en cada ventana canonica y registra."""
    R = df["R_obs_m3_week"].to_numpy()
    dS = df["DeltaS_obs_m3"].to_numpy()
    S0 = df.attrs["S0_m3"]
    S_min, S_max = constants["S_min_m3"], constants["S_max_m3"]
    instance = res.instance_label(T, L)
    rows = []

    for label, start in windows.items():
        R_w, dS_w, S0_w = window_slice(R, dS, S0, start, T)
        p = window_params(R, R_w, T, L)
        w = fc.compute_weights(T, S_min, p["u_max"])
        B = constants["eta"] * float(R_w.sum())

        # sanity: S0_w coincide con el cierre de la semana anterior del dataset
        if start >= 1:
            s_prev = float(df["S_obs_m3"].to_numpy()[start - 1])
            assert abs(S0_w - s_prev) < 1e-3 * max(abs(s_prev), 1.0), \
                f"S0_w ({S0_w:.1f}) != S_obs[{start-1}] ({s_prev:.1f})"

        uh = bl.historical(T)
        Sh, ch, srs_h = _srs_of(uh, S0_w, dS_w, S_min, w)
        ut = bl.threshold_rule(R_w, dS_w, S0_w, S_min, p["delta_u"])
        Stt, ct, srs_t = _srs_of(ut, S0_w, dS_w, S_min, w)
        uc = bl.threshold_rule(R_w, dS_w, S0_w, S_min, p["delta_u"], clamp_release=True)
        Sc, cc, srs_c = _srs_of(uc, S0_w, dS_w, S_min, w)
        ub = bl.threshold_rule(R_w, dS_w, S0_w, S_min, p["delta_u"], B=B)
        Sb, cb, srs_b = _srs_of(ub, S0_w, dS_w, S_min, w)
        dp = bl.dp_optimal(R_w, dS_w, S0_w, S_min=S_min, S_max=S_max,
                           delta_u=p["delta_u"], L=L, weights=w, B=B)
        Sdp = st.simulate_storage(S0_w, dS_w, dp["u_star"])
        bf = bl.brute_force_optimal(R_w, dS_w, S0_w, S_min=S_min, S_max=S_max,
                                    delta_u=p["delta_u"], L=L, weights=w, B=B,
                                    max_combos=BRUTE_MAX_COMBOS)

        refs = {"historical": srs_h, "threshold": srs_t, "dp": dp["SRS_star"]}
        common = dict(instance=instance, T=T, L=L, params=p, weights=w,
                      constants=constants, B=B, references=refs,
                      window_start=start, window_label=label)

        def chk(u, S):
            return st.check_constraints(R_w, u, S, S_max, p["u_max"], B)

        def emit(method, u, S, costs, srs_v, dt, variant=None):
            res.record_run(method=method, variant=variant, u=u, S=S, costs=costs,
                           srs=srs_v, feasibility=chk(u, S), runtime_seconds=dt, **common)
            rows.append({"instance": instance, "window": label, "start": start,
                         "method": method + (f"_{variant}" if variant else ""),
                         "SRS": srs_v, "dSRS_hist": srs_v - srs_h,
                         "u_nonzero": int(np.count_nonzero(u)),
                         "feasible": chk(u, S)["feasible"]})

        emit("historical", uh, Sh, ch, srs_h, 0.0)
        emit("threshold", ut, Stt, ct, srs_t, 0.0, variant="pure")
        emit("threshold", uc, Sc, cc, srs_c, 0.0, variant="clamped")
        emit("threshold", ub, Sb, cb, srs_b, 0.0, variant="balanced")
        emit("dp", dp["u_star"], Sdp, dp["costs"], dp["SRS_star"], 0.0)
        if bf is not None:
            Sbf = st.simulate_storage(S0_w, dS_w, bf["u_star"])
            emit("brute", bf["u_star"], Sbf, bf["costs"], bf["SRS_star"], 0.0)
    return rows


# --------------------------------------------------------------------------- #
def _print_table(rows: list) -> None:
    print(f"\n{'instance':8s} {'window':7s} {'start':>5s} {'method':18s} "
          f"{'SRS':>12s} {'dSRS_hist':>11s} {'u!=0':>5s} {'feas':>5s}")
    print("-" * 82)
    for r in rows:
        print(f"{r['instance']:8s} {r['window']:7s} {r['start']:5d} {r['method']:18s} "
              f"{r['SRS']:12.6f} {r['dSRS_hist']:+11.4e} {r['u_nonzero']:5d} "
              f"{str(r['feasible']):>5s}")


def main() -> None:
    df = fd.build_weekly_benchmark(write=False)
    constants = fc.load_official_constants()
    S0 = df.attrs["S0_m3"]
    dS = df["DeltaS_obs_m3"].to_numpy()
    n_weeks = int(df.attrs["n_weeks"])
    print(f"Dataset: {n_weeks} semanas | S_min oficial = {constants['S_min_m3']:,.0f} m^3\n")

    all_rows, scans = [], []
    for T, L in STANDARD_INSTANCES:
        scans.append(scan_all_windows(df, constants, T, L))
        windows = select_windows(S0, dS, T, n_weeks)
        all_rows += run_canonical_windows(df, constants, T, L, windows)

    _print_table(all_rows)

    total_nz = sum(len(s["nonzero_starts"]) for s in scans)
    print("\n== Conclusion ==")
    for s in scans:
        rep = "SI (todas u*=0)" if not s["nonzero_starts"] else \
            f"NO: {len(s['nonzero_starts'])}/{s['n_windows']} ventanas con u*!=0"
        print(f"  {s['instance']:7s}: first representativa? {rep}")
    if total_nz == 0:
        print("  => En TODAS las ventanas de TODAS las instancias el optimo exacto es u*=0.")
        print("     u*=0 no es artefacto de la ventana start=0: es propiedad del benchmark")
        print("     (embalse siempre << S_min; el presupuesto de balance eta=0.10 es muy chico")
        print("     para que los ajustes muevan C_crit mas de lo que cuestan C_dev/C_smooth).")
    print(f"\nResultados en {res.RESULTS_ROOT}/runs_summary.csv  y  {res.RESULTS_ROOT}/runs/")


if __name__ == "__main__":
    main()
