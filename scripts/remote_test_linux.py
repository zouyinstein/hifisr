#!/usr/bin/env python3
# HiFiSR module guide:
# - base: command, file, and soft_paths helpers; import hifisr_functions.base as hfbase
# - reads: read extraction, filtering, sampling, and correction; import hifisr_functions.reads as hfreads
# - references: reference rotation, assembly, polishing, and alignment; import hifisr_functions.references as hfref
# - variants: read-variant calling, grouping, and frequency analysis; import hifisr_functions.variants as hfvar
# - transfer: organelle/nuclear transfer-fragment analysis; import hifisr_functions.transfer as hftrans
# - annotations: annotation tables and feature-level summaries; import hifisr_functions.annotations as hfanno
# - reports: read statistics, plots, Excel tables, and report outputs; import hifisr_functions.reports as hfrps

"""Sync this repository to a Linux server and run tests remotely."""

from __future__ import annotations

import argparse
import shlex
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run(command: list[str]) -> None:
    print("+", " ".join(command))
    subprocess.run(command, cwd=ROOT, check=True)


def remote_join(commands: list[str]) -> str:
    return " && ".join(commands)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("remote_host", help="Remote SSH target, for example user@host.")
    parser.add_argument(
        "remote_dir",
        nargs="?",
        default="~/hifisr-dev",
        help="Remote directory used for the synced repository.",
    )
    parser.add_argument(
        "--python",
        default="python3",
        help="Remote Python executable used to create the virtual environment.",
    )
    parser.add_argument(
        "--pytest-args",
        default="",
        help="Extra pytest arguments run on the remote server.",
    )
    args = parser.parse_args()

    remote_dir = args.remote_dir
    quoted_dir = shlex.quote(remote_dir)

    run(["ssh", args.remote_host, f"mkdir -p {quoted_dir}"])
    run(
        [
            "rsync",
            "-az",
            "--exclude",
            ".git/",
            "--exclude",
            ".venv/",
            "--exclude",
            "__pycache__/",
            "--exclude",
            ".pytest_cache/",
            "./",
            f"{args.remote_host}:{remote_dir}/",
        ]
    )

    remote_commands = [
        f"cd {quoted_dir}",
        f"{shlex.quote(args.python)} -m venv .venv",
        ". .venv/bin/activate",
        "python -m pip install --upgrade pip",
        "python -m pip install -r requirements-dev.txt",
        f"python -m pytest {args.pytest_args}".strip(),
    ]
    run(["ssh", args.remote_host, remote_join(remote_commands)])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

