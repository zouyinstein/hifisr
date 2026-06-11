#!/usr/bin/env python3
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

