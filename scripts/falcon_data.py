"""Carga de los CSV del IBWC y construccion canonica del dataset semanal Falcon.

Convenciones DECIDIDAS (docs/SPEC_IMPLEMENTACION_QUBO.md seccion 2):
- Unidades en m^3 (SI). La unidad se autodetecta del header ``Value (...)``.
- R_obs semanal = integracion nativa del caudal 15-min: ``sum(caudal * 900 s)``.
- S_obs = valor de cierre de semana (almacenamiento en la frontera de fin de semana).
- Semanas = bins de 7 dias anclados al primer timestamp del storage; se descarta el
  remanente parcial -> 52 semanas completas.
- DeltaS_obs = S(t+1) - S(t), driver del modelo. Ahora CROSS-VALIDADO contra el
  dataset oficial `Discharge.Total.Change-in-Storage@08461200`: a nivel semanal,
  DeltaS_obs == -(suma semanal del oficial) EXACTO (corr +1.000, err 0). OJO: el
  oficial usa la CONVENCION DE SIGNO OPUESTA a la ecuacion de balance del spec; si se
  usa crudo, la trayectoria de storage se va al reves (hacia negativo). Por eso el
  driver sigue siendo el derivado (== oficial con el signo corregido). Ver
  `validate_deltaS_vs_official` y `docs/AUDITORIA_DATOS_Y_BRUTE.md`.

NO escribe en FalconChallenge/ (congelado): el dataset procesado va a data/processed/.
"""
from __future__ import annotations

import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW_DATA_DIR = REPO_ROOT / "FalconChallenge" / "data"
PROCESSED_DIR = REPO_ROOT / "data" / "processed"

WEEK_DAYS = 7
DELTAS_SOURCE = "derived_from_total_storage_not_official"
DELTAS_SOURCE_VALIDATED = "derived_from_total_storage_validated_vs_official_change_in_storage"


def _find_csv(data_dir: Path, substring: str) -> Path:
    hits = [p for p in glob.glob(str(Path(data_dir) / "*.csv")) if substring in p]
    if not hits:
        raise FileNotFoundError(f"No se encontro CSV con '{substring}' en {data_dir}")
    return Path(sorted(hits)[0])


def load_ibwc_csv(path: str | Path) -> tuple[pd.Series, str]:
    """Lee un export del IBWC y devuelve (serie indexada por timestamp, unidad).

    Maneja el header especial (linea 1 ``#Data Set Export...``; linea 2
    ``Timestamp (...),Value (<unidad>)``) y descarta la fila final de disclaimer.
    """
    header = pd.read_csv(path, skiprows=1, nrows=0)
    value_col = str(header.columns[1]) if len(header.columns) > 1 else ""
    m = re.search(r"\(([^)]+)\)", value_col)
    unit = m.group(1).strip() if m else ""

    df = pd.read_csv(path, skiprows=1)
    df.columns = ["ts", "value"]
    df = df[df["ts"].astype(str).str.match(r"\d{4}-\d{2}-\d{2}")].copy()
    df["ts"] = pd.to_datetime(df["ts"])
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    s = df.dropna().set_index("ts")["value"].sort_index()
    return s, unit


def _to_m3(series: pd.Series, unit: str) -> pd.Series:
    """Normaliza una serie de volumen a m^3 (acepta m^3 / TCM / MCM)."""
    u = unit.lower().replace(" ", "")
    if u in ("m^3", "m3"):
        return series
    if u == "tcm":
        return series * 1e3
    if u == "mcm":
        return series * 1e6
    raise ValueError(f"Unidad de volumen no reconocida: {unit!r}")


def _weekly_official_deltaS(boundaries, data_dir: Path):
    """DeltaS semanal (m^3) desde el dataset oficial `Change-in-Storage@08461200`,
    con el SIGNO CORREGIDO para la ecuacion de balance del spec.

    El oficial reporta la variacion con la convencion de signo OPUESTA (su suma anual
    es -272M mientras el storage subio +273M). Aca se suma el diario dentro de cada
    semana [b0,b1) y se NIEGA, de modo que coincide con el derivado S(t+1)-S(t).
    Devuelve un np.ndarray de largo n_weeks, o None si el CSV oficial no esta.
    """
    try:
        path = _find_csv(data_dir, "Change-in-Storage@08461200")
    except FileNotFoundError:
        return None
    dS_daily, unit = load_ibwc_csv(path)
    dS_daily = _to_m3(dS_daily, unit)
    out = []
    for i in range(len(boundaries) - 1):
        b0, b1 = boundaries[i], boundaries[i + 1]
        seg = dS_daily[(dS_daily.index >= b0) & (dS_daily.index < b1)]
        out.append(-float(seg.sum()))          # signo corregido
    return np.asarray(out, dtype=float)


def build_weekly_benchmark(data_dir: str | Path = RAW_DATA_DIR,
                           write: bool = True) -> pd.DataFrame:
    """Construye el dataset semanal canonico en m^3 (52 semanas).

    Columnas: week, week_start, week_end, S_obs_m3, R_obs_m3_week, DeltaS_obs_m3,
    DeltaS_source. ``S_obs_m3`` es el cierre de semana (frontera de fin); el
    almacenamiento inicial S0 es la frontera de inicio de la semana 1
    (atributo ``df.attrs['S0_m3']``).
    """
    data_dir = Path(data_dir)

    storage_raw, s_unit = load_ibwc_csv(_find_csv(data_dir, "Total Storage.Web-Daily-tcm@08461200"))
    storage = _to_m3(storage_raw, s_unit)  # diario, m^3

    rel_raw, r_unit = load_ibwc_csv(_find_csv(data_dir, "Discharge.Best Available@08461300"))
    if r_unit.lower().replace(" ", "") not in ("m^3/s", "m3/s"):
        raise ValueError(f"Se esperaba caudal en m^3/s; header dice {r_unit!r}")
    release = rel_raw  # m^3/s, 15-min

    # Fronteras semanales ancladas al primer timestamp diario del storage.
    first = pd.Timestamp(storage.index[0])
    boundaries = []
    k = 0
    while True:
        b = first + pd.Timedelta(days=WEEK_DAYS * k)
        if b in storage.index:
            boundaries.append(b)
            k += 1
        else:
            break
    n_weeks = len(boundaries) - 1  # semanas completas
    if n_weeks < 1:
        raise RuntimeError("No se pudieron formar semanas completas de 7 dias.")

    rows = []
    for i in range(n_weeks):
        b0, b1 = boundaries[i], boundaries[i + 1]
        s_close = float(storage.loc[b1])            # cierre de semana
        s_open = float(storage.loc[b0])
        # Integracion nativa del caudal sobre la ventana: volumen = caudal_medio *
        # duracion de la semana. Para semanas uniformes (672 muestras a 15-min)
        # equivale a sum(caudal * 900 s); es robusto a huecos de muestreo (los
        # trata como el caudal promedio). Unit-safe (no depende de la resolucion
        # de datetime de pandas).
        win = release[(release.index >= b0) & (release.index < b1)]
        week_seconds = (b1 - b0).total_seconds()
        r_vol = float(win.mean() * week_seconds) if len(win) else 0.0  # m^3
        rows.append({
            "week": i + 1,
            "week_start": b0.date().isoformat(),
            "week_end": b1.date().isoformat(),
            "S_obs_m3": s_close,
            "R_obs_m3_week": r_vol,
            "DeltaS_obs_m3": s_close - s_open,      # derivado
            "DeltaS_source": DELTAS_SOURCE,
        })

    df = pd.DataFrame(rows)
    df.attrs["S0_m3"] = float(storage.loc[boundaries[0]])
    df.attrs["units"] = "m^3"
    df.attrs["n_weeks"] = n_weeks

    # Cross-validacion contra el oficial (si esta presente): agrega columna con el
    # oficial signo-corregido y marca la fuente como validada. El DRIVER
    # (DeltaS_obs_m3 = derivado) NO cambia -> numeros identicos.
    dS_official = _weekly_official_deltaS(boundaries, Path(data_dir))
    if dS_official is not None:
        df["DeltaS_official_m3"] = dS_official
        derived = df["DeltaS_obs_m3"].to_numpy()
        err = np.abs(derived - dS_official)
        df["DeltaS_source"] = DELTAS_SOURCE_VALIDATED
        df.attrs["deltaS_validation"] = {
            "max_abs_err_m3": float(err.max()),
            "median_abs_err_m3": float(np.median(err)),
            "sum_derived_m3": float(derived.sum()),
            "sum_official_signcorrected_m3": float(dS_official.sum()),
            "sign_note": "oficial usa signo opuesto; DeltaS_official_m3 ya viene negado",
        }

    if write:
        PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
        out = PROCESSED_DIR / "falcon_weekly_benchmark.csv"
        df.to_csv(out, index=False)
        df.attrs["path"] = str(out)
    return df


def validate_deltaS_vs_official(data_dir: str | Path = RAW_DATA_DIR,
                                tol: float = 1.0) -> dict:
    """Invariante verificable: el DeltaS derivado (driver) == oficial signo-corregido.

    A nivel semanal deben coincidir EXACTO (err ~0). Afirma
    `max|DeltaS_obs_m3 - DeltaS_official_m3| < tol` y devuelve estadisticas. Sirve
    para blindar la convencion de signo del dataset oficial (GUIDELINES seccion 9.3).
    Lanza si el CSV oficial no esta presente.
    """
    df = build_weekly_benchmark(data_dir=data_dir, write=False)
    if "DeltaS_official_m3" not in df.columns:
        raise FileNotFoundError(
            "Falta el dataset oficial Change-in-Storage@08461200 para validar DeltaS.")
    derived = df["DeltaS_obs_m3"].to_numpy()
    official = df["DeltaS_official_m3"].to_numpy()
    err = np.abs(derived - official)
    corr = float(np.corrcoef(derived, official)[0, 1])
    stats = {
        "n_weeks": int(df.attrs["n_weeks"]),
        "max_abs_err_m3": float(err.max()),
        "median_abs_err_m3": float(np.median(err)),
        "corr": corr,
        "sum_derived_m3": float(derived.sum()),
        "sum_official_signcorrected_m3": float(official.sum()),
    }
    assert err.max() < tol, (
        f"DeltaS derivado NO coincide con el oficial signo-corregido "
        f"(max|Δ|={err.max():.3e} m^3 >= tol={tol}).")
    return stats


if __name__ == "__main__":
    df = build_weekly_benchmark()
    print(f"Semanas: {df.attrs['n_weeks']} | S0 = {df.attrs['S0_m3']:,.0f} m^3")
    print(df.head(3).to_string(index=False))
    print("...")
    print(df.tail(2).to_string(index=False))
    if "DeltaS_official_m3" in df.columns:
        v = validate_deltaS_vs_official()
        print(f"\nValidacion DeltaS vs oficial (signo-corregido): "
              f"max|Δ|={v['max_abs_err_m3']:.3e} m^3  median|Δ|={v['median_abs_err_m3']:.3e}  "
              f"corr={v['corr']:+.4f}")
        print(f"  sum derived={v['sum_derived_m3']/1e6:+.1f}M  "
              f"sum oficial(signo-corr)={v['sum_official_signcorrected_m3']/1e6:+.1f}M  -> INVARIANTE OK")
    else:
        print("\n(dataset oficial Change-in-Storage ausente: sin cross-validacion)")
