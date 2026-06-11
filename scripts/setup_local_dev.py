#!/usr/bin/env python3
# HiFiSR module guide:
# - base: command, file, and soft_paths helpers; import hifisr_functions.base as hfbase
# - reads: read extraction, filtering, sampling, and correction; import hifisr_functions.reads as hfreads
# - references: reference rotation, assembly, polishing, and alignment; import hifisr_functions.references as hfref
# - variants: read-variant calling, grouping, and frequency analysis; import hifisr_functions.variants as hfvar
# - transfer: organelle/nuclear transfer-fragment analysis; import hifisr_functions.transfer as hftrans
# - annotations: annotation tables and feature-level summaries; import hifisr_functions.annotations as hfanno
# - reports: read statistics, plots, Excel tables, and report outputs; import hifisr_functions.reports as hfrps

"""Create a local hifisr development environment."""

from __future__ import annotations

import argparse
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOFT_PATH_TOOLS = [
    "minimap2",
    "samtools",
    "seqkit",
    "mecat",
    "blastn",
    "bcftools",
    "bamtools",
    "meryl",
    "winnowmap",
    "pigz",
    "bandage",
    "hifiasm",
    "flye",
    "canu",
]

LOCAL_SOFT_PATHS = {
    "mecat": ROOT / "deps" / "MECAT2-20190314" / "Darwin-arm64" / "bin" / "mecat.pl",
    "meryl": ROOT / "deps" / "canu-2.3" / "bin" / "meryl",
    "winnowmap": ROOT / "deps" / "Winnowmap-2.03" / "bin" / "winnowmap",
    "bandage": ROOT
    / "deps"
    / "Bandage-0.9.0"
    / "Bandage.app"
    / "Contents"
    / "MacOS"
    / "Bandage",
    "hifiasm": ROOT / "deps" / "hifiasm-0.25.0" / "hifiasm",
    "canu": ROOT / "deps" / "canu-2.3" / "bin" / "canu",
}


def venv_python(venv_dir: Path) -> Path:
    if os.name == "nt":
        return venv_dir / "Scripts" / "python.exe"
    return venv_dir / "bin" / "python"


def run(command: list[str]) -> None:
    print("+", " ".join(command))
    subprocess.run(command, cwd=ROOT, check=True)


def soft_paths_header() -> list[str]:
    system = platform.system() or "unknown"
    machine = platform.machine() or "unknown"
    if system == "Darwin":
        note = "macOS local development only; replace this file on Ubuntu/Linux."
    elif system == "Linux":
        note = "Linux runtime paths."
    else:
        note = "local runtime paths."
    return [
        f"# hifisr soft_paths generated for {system} {machine}",
        f"# {note}",
        "# Format: software_name<TAB>absolute_path_to_executable",
    ]


def resolve_tool(tool: str) -> str | None:
    local_path = LOCAL_SOFT_PATHS.get(tool) if platform.system() == "Darwin" else None
    if local_path is not None and local_path.exists():
        return str(local_path)
    return shutil.which(tool)


def write_soft_paths(venv_dir: Path, soft_paths_file: Path) -> list[str]:
    soft_paths_file.parent.mkdir(parents=True, exist_ok=True)
    missing = []
    with soft_paths_file.open("w") as handle:
        for line in soft_paths_header():
            handle.write(f"{line}\n")
        handle.write(f"python\t{venv_python(venv_dir)}\n")
        for tool in SOFT_PATH_TOOLS:
            path = resolve_tool(tool)
            if path is None:
                missing.append(tool)
                handle.write(f"# missing\t{tool}\n")
            else:
                handle.write(f"{tool}\t{path}\n")
    return missing


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python executable used to create the virtual environment.",
    )
    parser.add_argument(
        "--venv",
        default=".venv",
        help="Virtual environment directory relative to the repository root.",
    )
    parser.add_argument(
        "--skip-install",
        action="store_true",
        help="Only create the virtual environment; do not install requirements.",
    )
    parser.add_argument(
        "--soft-paths",
        default="deps/soft_paths_macOS.txt",
        help="Path written with local third-party executable locations.",
    )
    parser.add_argument(
        "--skip-soft-paths",
        action="store_true",
        help="Do not write the local soft_paths file.",
    )
    args = parser.parse_args()

    venv_dir = ROOT / args.venv
    run([args.python, "-m", "venv", str(venv_dir)])

    if not args.skip_install:
        python = venv_python(venv_dir)
        run([str(python), "-m", "pip", "install", "--upgrade", "pip"])
        run([str(python), "-m", "pip", "install", "-r", "requirements-dev.txt"])

    if not args.skip_soft_paths:
        soft_paths_file = Path(args.soft_paths)
        if not soft_paths_file.is_absolute():
            soft_paths_file = ROOT / soft_paths_file
        missing = write_soft_paths(venv_dir, soft_paths_file)
        print(f"Wrote {soft_paths_file}")
        if missing:
            print("Missing external tools: " + ", ".join(missing))

    print("Development environment is ready.")
    print(f"Run: {venv_python(venv_dir)} -m pytest")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
