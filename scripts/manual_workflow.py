#!/usr/bin/env python3
"""Shared helpers for the manual single-sample HiFiSR runners."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys
import time


SCRIPT_PATH = Path(__file__).resolve()
DEFAULT_HIFISR_ROOT = SCRIPT_PATH.parents[1]


def add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--hifisr-root", default=str(DEFAULT_HIFISR_ROOT))
    parser.add_argument("--env-tsv", default=None)
    parser.add_argument("--samples-tsv", default=None)
    parser.add_argument("--env", default="macOS")
    parser.add_argument("--sample", default="W3-5-2")
    parser.add_argument(
        "--output-sample",
        default=None,
        help="Result folder/sample index to write. Defaults to --sample.",
    )
    parser.add_argument("--results-dir", default=None)
    parser.add_argument("--sample-dir", default=None)
    parser.add_argument("--python", default=None)
    parser.add_argument("--soft-paths", default=None)
    parser.add_argument("--raw-reads", default=None)
    parser.add_argument("--mito-ref", default=None)
    parser.add_argument("--plastid-ref", default=None)
    parser.add_argument("--threads", default=None)
    parser.add_argument("--mito-read-limit", default=None)
    parser.add_argument("--plastid-read-limit", default=None)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without running analysis scripts.",
    )
    parser.add_argument(
        "--no-create-dirs",
        action="store_true",
        help="Do not create sample/log/cache directories before running.",
    )


def import_init_helpers(hifisr_root: Path):
    scripts_dir = hifisr_root / "scripts"
    if str(scripts_dir) not in sys.path:
        sys.path.insert(0, str(scripts_dir))
    import init_manual_workflow as init_helpers  # type: ignore

    return init_helpers


def _set_if_present(values: dict[str, str], key: str, value: str | None) -> None:
    if value is not None and value != "":
        values[key] = value


def load_context(args: argparse.Namespace) -> dict[str, str]:
    hifisr_root = Path(args.hifisr_root).expanduser().resolve()
    init_helpers = import_init_helpers(hifisr_root)

    env_tsv = Path(args.env_tsv or hifisr_root / "workflow/config/manual_env_macOS.tsv")
    samples_tsv = Path(args.samples_tsv or hifisr_root / "workflow/config/w3_5_2_samples.tsv")

    env_rows = init_helpers.read_rows(env_tsv, init_helpers.ENV_REQUIRED_COLUMNS)
    sample_rows = init_helpers.read_rows(samples_tsv, init_helpers.SAMPLE_REQUIRED_COLUMNS)
    env_values = dict(init_helpers.select_row(env_rows, "env", args.env, "Environment"))
    sample_values = dict(init_helpers.select_row(sample_rows, "sample", args.sample, "Sample"))

    env_values["hifisr_root"] = str(hifisr_root)
    _set_if_present(env_values, "results_dir", args.results_dir)
    _set_if_present(env_values, "python", args.python)
    _set_if_present(env_values, "soft_paths", args.soft_paths)
    _set_if_present(env_values, "threads", args.threads)
    _set_if_present(sample_values, "raw_reads", args.raw_reads)
    _set_if_present(sample_values, "mito_ref", args.mito_ref)
    _set_if_present(sample_values, "plastid_ref", args.plastid_ref)
    _set_if_present(sample_values, "mito_read_limit", args.mito_read_limit)
    _set_if_present(sample_values, "plastid_read_limit", args.plastid_read_limit)

    values = init_helpers.expand_templates({**env_values, **sample_values})
    values["input_sample"] = values["sample"]

    output_sample = args.output_sample or args.sample
    values["sample"] = output_sample
    if args.sample_dir:
        values["sample_dir"] = str(Path(args.sample_dir).expanduser().resolve())
    else:
        values["sample_dir"] = str(Path(values["results_dir"]) / output_sample)
    values["manual_log_dir"] = str(Path(values["sample_dir"]) / "logs/manual")

    init_helpers.validate(values, validate_existing=True)
    if not args.no_create_dirs:
        init_helpers.create_dirs(values)

    init_env = Path(values["manual_log_dir"]) / "init_env.sh"
    init_env.write_text(init_helpers.shell_exports(values))
    write_context_tsv(Path(values["manual_log_dir"]) / "run_context.tsv", values)
    return values


def write_context_tsv(path: Path, values: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt") as handle:
        print("key\tvalue", file=handle)
        for key in sorted(values):
            print(key + "\t" + str(values[key]), file=handle)


def runner_env(values: dict[str, str]) -> dict[str, str]:
    env = os.environ.copy()
    env["HIFISR_ROOT"] = values["hifisr_root"]
    env["RESULTS_DIR"] = values["results_dir"]
    env["SAMPLE"] = values["sample"]
    env["SAMPLE_DIR"] = values["sample_dir"]
    env["PY"] = values["python"]
    env["SOFT_PATHS"] = values["soft_paths"]
    env["RAW_READS"] = values["raw_reads"]
    env["MITO_REF"] = values["mito_ref"]
    env["PLASTID_REF"] = values["plastid_ref"]
    env["THREADS"] = values["threads"]
    env["MITO_READ_LIMIT"] = values["mito_read_limit"]
    env["PLASTID_READ_LIMIT"] = values["plastid_read_limit"]
    env["MANUAL_LOG_DIR"] = values["manual_log_dir"]
    env["HIFISR_TMPDIR"] = values["hifisr_tmpdir"]
    env["MPLCONFIGDIR"] = values["mplconfigdir"]
    return env


def sample_path(values: dict[str, str], *parts: str) -> str:
    return str(Path(values["sample_dir"], *parts))


def script_path(values: dict[str, str], script_name: str) -> str:
    return str(Path(values["hifisr_root"]) / "analysis_scripts" / script_name)


def run_logged(
    values: dict[str, str],
    log_name: str,
    command: list[str],
    *,
    dry_run: bool = False,
) -> None:
    log_path = Path(values["manual_log_dir"]) / log_name
    printable = " ".join(command)
    print("+ " + printable)
    if dry_run:
        print("  log: " + str(log_path))
        return

    log_path.parent.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    with log_path.open("wt") as log:
        print("+ " + printable, file=log)
        print("cwd: " + values["results_dir"], file=log)
        print("", file=log)
        completed = subprocess.run(
            command,
            cwd=values["results_dir"],
            env=runner_env(values),
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )
        elapsed = time.monotonic() - start
        print("", file=log)
        print("exit_status\t" + str(completed.returncode), file=log)
        print("elapsed_seconds\t" + f"{elapsed:.3f}", file=log)
    if completed.returncode != 0:
        raise RuntimeError("Command failed; see log: " + str(log_path))


def require_files(paths: list[str], *, allow_empty: bool = False) -> None:
    missing = [path for path in paths if not Path(path).exists()]
    empty = [
        path
        for path in paths
        if Path(path).exists()
        and not allow_empty
        and Path(path).is_file()
        and Path(path).stat().st_size == 0
    ]
    if missing or empty:
        lines = []
        if missing:
            lines.append("Missing required file(s):")
            lines.extend("  " + path for path in missing)
        if empty:
            lines.append("Empty required file(s):")
            lines.extend("  " + path for path in empty)
        raise FileNotFoundError("\n".join(lines))


def print_next_files(title: str, paths: list[str]) -> None:
    print("\n" + title)
    for path in paths:
        status = "OK" if Path(path).exists() and Path(path).stat().st_size > 0 else "NEEDED"
        print("  [" + status + "] " + path)
