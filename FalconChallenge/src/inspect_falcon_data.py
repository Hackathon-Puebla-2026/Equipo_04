from __future__ import annotations

import csv
import json
import re
from datetime import datetime
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
MANIFEST_PATH = DATA_DIR / "falcon_download_manifest.json"
OWNERSHIP_REPORTS_PATH = DATA_DIR / "ibwc_weekly_ownership_reports.txt"

DATE_NAME_HINTS = (
    "date",
    "fecha",
    "time",
    "timestamp",
    "datetime",
    "updated",
    "last_updated",
    "day",
)
DATE_FORMATS = (
    "%Y-%m-%d",
    "%Y/%m/%d",
    "%m/%d/%Y",
    "%d/%m/%Y",
    "%m-%d-%Y",
    "%d-%m-%Y",
    "%Y-%m-%d %H:%M:%S",
    "%Y/%m/%d %H:%M:%S",
    "%m/%d/%Y %H:%M:%S",
    "%d/%m/%Y %H:%M:%S",
    "%m/%d/%y",
    "%d/%m/%y",
    "%b %d, %Y",
    "%B %d, %Y",
)
DATE_PATTERN = re.compile(
    r"\b("
    r"\d{4}[-/]\d{1,2}[-/]\d{1,2}"
    r"|"
    r"\d{1,2}[-/]\d{1,2}[-/]\d{2,4}"
    r")\b"
)


def read_text(path: Path) -> str:
    for encoding in ("utf-8-sig", "utf-8", "latin-1"):
        try:
            return path.read_text(encoding=encoding)
        except UnicodeDecodeError:
            continue
    return path.read_text(errors="replace")


def truncate(value: Any, limit: int = 120) -> str:
    text = "" if value is None else str(value)
    text = text.replace("\n", " ").replace("\r", " ")
    if len(text) <= limit:
        return text
    return text[: limit - 3] + "..."


def parse_date(value: Any) -> datetime | None:
    if value is None:
        return None

    text = str(value).strip()
    if not text:
        return None

    normalized = text.replace("T", " ").replace("Z", "")
    try:
        return datetime.fromisoformat(normalized)
    except ValueError:
        pass

    for fmt in DATE_FORMATS:
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            continue

    match = DATE_PATTERN.search(text)
    if match:
        return parse_date(match.group(1))

    return None


def detect_date_columns(rows: list[dict[str, Any]], columns: list[str]) -> dict[str, tuple[datetime, datetime]]:
    detected: dict[str, tuple[datetime, datetime]] = {}

    for column in columns:
        values = [row.get(column) for row in rows if row.get(column) not in (None, "")]
        if not values:
            continue

        parsed_dates = [parsed for value in values if (parsed := parse_date(value))]
        has_name_hint = any(hint in column.lower() for hint in DATE_NAME_HINTS)
        enough_dates = len(parsed_dates) >= max(1, min(3, len(values)))

        if parsed_dates and (has_name_hint or enough_dates):
            detected[column] = (min(parsed_dates), max(parsed_dates))

    return detected


def print_section(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def print_rows(rows: list[dict[str, Any]]) -> None:
    if not rows:
        print("  (sin filas)")
        return

    for index, row in enumerate(rows[:5], start=1):
        compact_row = {key: truncate(value) for key, value in row.items()}
        print(f"  {index}. {compact_row}")


def print_report(
    path: Path,
    row_count: int,
    columns: list[str],
    first_rows: list[dict[str, Any]],
    date_ranges: dict[str, tuple[datetime, datetime]],
) -> None:
    print_section(str(path.relative_to(PROJECT_ROOT)))
    print(f"Nombre: {path.name}")
    print(f"Numero de filas: {row_count}")
    print(f"Numero de columnas: {len(columns)}")
    print(f"Nombres de columnas: {columns if columns else '(no detectadas)'}")
    print("Primeras 5 filas:")
    print_rows(first_rows)

    print("Posibles columnas de fecha:")
    if date_ranges:
        for column, (start, end) in date_ranges.items():
            print(f"  - {column}")
            print(f"    Rango de fechas: {start.date().isoformat()} a {end.date().isoformat()}")
    else:
        print("  (no detectadas)")
        print("Rango de fechas: (no detectado)")


def read_csv_file(path: Path) -> tuple[list[str], list[dict[str, Any]]]:
    text = read_text(path)
    sample = text[:4096]

    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",;\t|")
    except csv.Error:
        dialect = csv.excel

    rows = list(csv.DictReader(text.splitlines(), dialect=dialect))
    columns = list(rows[0].keys()) if rows else []
    return columns, rows


def inspect_csv(path: Path) -> None:
    columns, rows = read_csv_file(path)
    print_report(
        path=path,
        row_count=len(rows),
        columns=columns,
        first_rows=rows[:5],
        date_ranges=detect_date_columns(rows, columns),
    )


def json_to_rows(data: Any) -> tuple[list[str], list[dict[str, Any]]]:
    if isinstance(data, list):
        if all(isinstance(item, dict) for item in data):
            columns = sorted({key for item in data for key in item.keys()})
            rows = [{column: item.get(column) for column in columns} for item in data]
            return columns, rows
        return ["value"], [{"value": item} for item in data]

    if isinstance(data, dict):
        rows = [{"key": key, "value": value} for key, value in data.items()]
        return ["key", "value"], rows

    return ["value"], [{"value": data}]


def inspect_json(path: Path) -> None:
    data = json.loads(read_text(path))
    columns, rows = json_to_rows(data)
    print_report(
        path=path,
        row_count=len(rows),
        columns=columns,
        first_rows=rows[:5],
        date_ranges=detect_date_columns(rows, columns),
    )


def inspect_text(path: Path) -> None:
    lines = read_text(path).splitlines()
    rows = [{"line": line} for line in lines]
    print_report(
        path=path,
        row_count=len(rows),
        columns=["line"],
        first_rows=rows[:5],
        date_ranges=detect_date_columns(rows, ["line"]),
    )


def csv_paths() -> list[Path]:
    paths: list[Path] = []
    if DATA_DIR.exists():
        paths.extend(sorted(DATA_DIR.glob("*.csv")))
    if RAW_DIR.exists():
        paths.extend(sorted(RAW_DIR.glob("*.csv")))
    return paths


def main() -> None:
    for path in csv_paths():
        inspect_csv(path)

    if MANIFEST_PATH.exists():
        inspect_json(MANIFEST_PATH)
    else:
        print(f"\nNo existe {MANIFEST_PATH.relative_to(PROJECT_ROOT)}")

    if OWNERSHIP_REPORTS_PATH.exists():
        inspect_text(OWNERSHIP_REPORTS_PATH)
    else:
        print(f"\nNo existe {OWNERSHIP_REPORTS_PATH.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
