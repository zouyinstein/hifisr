#!/usr/bin/env python3
"""Initialize shell variables for the manual HiFiSR workflow from TSV files."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import shlex
import sys


ENV_REQUIRED_COLUMNS = ["env", "hifisr_root", "python", "soft_paths", "threads"]
SAMPLE_REQUIRED_COLUMNS = ["sample", "raw_reads", "mito_ref", "plastid_ref"]

EXPORT_MAP = {
    "env": "HIFISR_ENV",
    "hifisr_root": "HIFISR_ROOT",
    "results_dir": "RESULTS_DIR",
    "sample": "SAMPLE",
    "sample_dir": "SAMPLE_DIR",
    "python": "PY",
    "soft_paths": "SOFT_PATHS",
    "raw_reads": "RAW_READS",
    "mito_ref": "MITO_REF",
    "plastid_ref": "PLASTID_REF",
    "threads": "THREADS",
    "mito_read_limit": "MITO_READ_LIMIT",
    "plastid_read_limit": "PLASTID_READ_LIMIT",
    "manual_log_dir": "MANUAL_LOG_DIR",
    "hifisr_tmpdir": "HIFISR_TMPDIR",
    "mplconfigdir": "MPLCONFIGDIR",
}

PATH_KEYS = [
    "hifisr_root",
    "results_dir",
    "sample_dir",
    "python",
    "soft_paths",
    "raw_reads",
    "mito_ref",
    "plastid_ref",
    "manual_log_dir",
    "hifisr_tmpdir",
    "mplconfigdir",
]

EXISTING_INPUT_KEYS = [
    "hifisr_root",
    "python",
    "soft_paths",
    "raw_reads",
    "mito_ref",
    "plastid_ref",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read platform and sample TSV files, validate paths, create runtime "
            "directories, and emit shell exports for the manual HiFiSR workflow."
        )
    )
    parser.add_argument("env_tsv", help="Headered platform/environment TSV.")
    parser.add_argument("samples_tsv", help="Headered sample TSV.")
    parser.add_argument("--env", default=None, help="Environment row to select.")
    parser.add_argument("--sample", default=None, help="Sample row to select.")
    parser.add_argument(
        "--shell-out",
        default=None,
        help="Write shell exports to this file. Exports are always printed to stdout.",
    )
    parser.add_argument(
        "--create-dirs",
        action="store_true",
        help="Create sample_dir, manual_log_dir, hifisr_tmpdir, and mplconfigdir.",
    )
    parser.add_argument(
        "--no-validate-existing",
        action="store_true",
        help="Skip existence checks for input files and tools.",
    )
    return parser.parse_args()


def read_rows(tsv_path: Path, required_columns: list[str]) -> list[dict[str, str]]:
    with tsv_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(str(tsv_path) + " is missing a header row.")
        missing = [column for column in required_columns if column not in reader.fieldnames]
        if missing:
            raise ValueError(
                str(tsv_path) + " is missing required column(s): " + ", ".join(missing)
            )
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    if not rows:
        raise ValueError(str(tsv_path) + " has a header but no data rows.")
    return rows


def select_row(rows: list[dict[str, str]], key: str, value: str | None, label: str) -> dict[str, str]:
    if value is None:
        if len(rows) != 1:
            raise ValueError(f"{label} TSV has multiple rows; pass --{key}.")
        return rows[0]
    matches = [row for row in rows if row.get(key) == value]
    if not matches:
        raise ValueError(f"{label} not found in TSV: {value}")
    if len(matches) > 1:
        raise ValueError(f"{label} appears more than once in TSV: {value}")
    return matches[0]


def apply_defaults(values: dict[str, str]) -> dict[str, str]:
    values = dict(values)
    values.setdefault("results_dir", "")
    values.setdefault("sample_dir", "")
    values.setdefault("manual_log_dir", "")
    values.setdefault("hifisr_tmpdir", "")
    values.setdefault("mplconfigdir", "")
    values.setdefault("mito_read_limit", "50000")
    values.setdefault("plastid_read_limit", "50000")

    if not values["results_dir"]:
        values["results_dir"] = "{hifisr_root}/results"
    if not values["sample_dir"]:
        values["sample_dir"] = "{results_dir}/{sample}"
    if not values["manual_log_dir"]:
        values["manual_log_dir"] = "{sample_dir}/logs/manual"
    if not values["hifisr_tmpdir"]:
        values["hifisr_tmpdir"] = "{results_dir}/.tmp"
    if not values["mplconfigdir"]:
        values["mplconfigdir"] = "{results_dir}/.matplotlib"
    return values


def expand_templates(values: dict[str, str]) -> dict[str, str]:
    values = apply_defaults(values)
    for _ in range(10):
        changed = False
        for key, value in list(values.items()):
            try:
                expanded = value.format(**values)
            except KeyError as exc:
                raise ValueError(f"Unknown template key in column {key}: {exc}") from exc
            if expanded != value:
                values[key] = expanded
                changed = True
        if not changed:
            break

    for key in PATH_KEYS:
        if values.get(key):
            values[key] = str(Path(values[key]).expanduser())
    return values


def validate(values: dict[str, str], validate_existing: bool) -> None:
    for column in ENV_REQUIRED_COLUMNS + SAMPLE_REQUIRED_COLUMNS:
        if not values.get(column):
            raise ValueError("Required column has an empty value: " + column)
    for key in PATH_KEYS:
        value = values.get(key)
        if value and not os.path.isabs(value):
            raise ValueError(f"{key} must be an absolute path: {value}")
    if validate_existing:
        missing = [key for key in EXISTING_INPUT_KEYS if not Path(values[key]).exists()]
        if missing:
            details = "\n".join(f"  {key}: {values[key]}" for key in missing)
            raise FileNotFoundError("Missing required input path(s):\n" + details)
    int(values["threads"])
    int(values["mito_read_limit"])
    int(values["plastid_read_limit"])


def create_dirs(values: dict[str, str]) -> None:
    for key in ["sample_dir", "manual_log_dir", "hifisr_tmpdir", "mplconfigdir"]:
        Path(values[key]).mkdir(parents=True, exist_ok=True)


def shell_exports(values: dict[str, str]) -> str:
    lines = [
        "# Generated by scripts/init_manual_workflow.py",
        "# Source this file before running the manual HiFiSR workflow.",
    ]
    for key, env_name in EXPORT_MAP.items():
        value = values.get(key, "")
        if value:
            lines.append(f"export {env_name}={shlex.quote(value)}")
    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()
    try:
        env_rows = read_rows(Path(args.env_tsv), ENV_REQUIRED_COLUMNS)
        sample_rows = read_rows(Path(args.samples_tsv), SAMPLE_REQUIRED_COLUMNS)
        env_values = select_row(env_rows, "env", args.env, "Environment")
        sample_values = select_row(sample_rows, "sample", args.sample, "Sample")
        values = expand_templates({**env_values, **sample_values})
        validate(values, validate_existing=not args.no_validate_existing)
        if args.create_dirs:
            create_dirs(values)
        exports = shell_exports(values)
        sys.stdout.write(exports)
        if args.shell_out:
            shell_out = Path(args.shell_out)
            shell_out.parent.mkdir(parents=True, exist_ok=True)
            shell_out.write_text(exports)
    except Exception as exc:
        print("init_manual_workflow.py: " + str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
