"""Carga de los CSV del IBWC y construccion canonica del dataset semanal Falcon.

Convenciones DECIDIDAS (docs/SPEC_IMPLEMENTACION_QUBO.md seccion 2):
- Unidades en m^3 (SI). La unidad se autodetecta del header ``Value (...)``.
- R_obs semanal = integracion nativa del caudal 15-min: ``sum(caudal * 900 s)``.
- S_obs = valor de cierre de semana (almacenamiento en la frontera de fin de semana).
- Semanas = bins de 7 dias anclados al primer timestamp del storage; se descarta el
  remanente parcial -> 52 semanas completas.
- DeltaS_obs = S(t+1) - S(t), DERIVADO (falta dataset oficial) -> marcado preliminar.

NO escribe en FalconChallenge/ (congelado): el dataset procesado va a data/processed/.
"""
from __future__ import annotations

import glob
import re
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW_DATA_DIR = REPO_ROOT / "FalconChallenge" / "data"
PROCESSED_DIR = REPO_ROOT / "data" / "processed"

WEEK_DAYS = 7
DELTAS_SOURCE = "derived_from_total_storage_not_official"


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

    if write:
        PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
        out = PROCESSED_DIR / "falcon_weekly_benchmark.csv"
        df.to_csv(out, index=False)
        df.attrs["path"] = str(out)
    return df


if __name__ == "__main__":
    df = build_weekly_benchmark()
    print(f"Semanas: {df.attrs['n_weeks']} | S0 = {df.attrs['S0_m3']:,.0f} m^3")
    print(df.head(3).to_string(index=False))
    print("...")
    print(df.tail(2).to_string(index=False))
