#!/usr/bin/env python3
# HiFiSR module guide:
# - base: command, file, and soft_paths helpers; import hifisr_functions.base as hfbase
# - reads: read extraction, filtering, sampling, and correction; import hifisr_functions.reads as hfreads
# - references: reference rotation, assembly, polishing, and alignment; import hifisr_functions.references as hfref
# - variants: read-variant calling, grouping, and frequency analysis; import hifisr_functions.variants as hfvar
# - transfer: organelle/nuclear transfer-fragment analysis; import hifisr_functions.transfer as hftrans
# - annotations: annotation tables and feature-level summaries; import hifisr_functions.annotations as hfanno
# - reports: read statistics, plots, Excel tables, and report outputs; import hifisr_functions.reports as hfrps

"""Run local hifisr tests with the repository virtual environment if present."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def default_python() -> Path | str:
    if os.name == "nt":
        candidate = ROOT / ".venv" / "Scripts" / "python.exe"
    else:
        candidate = ROOT / ".venv" / "bin" / "python"
    return candidate if candidate.exists() else sys.executable


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pytest_args", nargs="*", help="Arguments passed to pytest.")
    parser.add_argument(
        "--python",
        default=str(default_python()),
        help="Python executable used to invoke pytest.",
    )
    args = parser.parse_args()

    command = [args.python, "-m", "pytest", *args.pytest_args]
    print("+", " ".join(command))
    return subprocess.run(command, cwd=ROOT, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())

