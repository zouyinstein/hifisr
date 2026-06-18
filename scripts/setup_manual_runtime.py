#!/usr/bin/env python3
"""Install/check the runtime used by the manual HiFiSR workflow."""

from __future__ import annotations

import argparse
import importlib.metadata as metadata
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

REQUIRED_TOOLS = [
    "python",
    "minimap2",
    "samtools",
    "seqkit",
    "blastn",
    "bcftools",
    "bamtools",
    "pigz",
    "flye",
    "simple_draft_asm",
    "gfa_editor_cli",
]

BREW_PACKAGES = [
    "python@3.11",
    "minimap2",
    "samtools",
    "seqkit",
    "blast",
    "bcftools",
    "bamtools",
    "pigz",
    "flye",
]

CONDA_PACKAGES = [
    "python=3.11",
    "pip",
    "minimap2",
    "samtools",
    "seqkit",
    "blast",
    "bcftools",
    "bamtools",
    "pigz",
    "flye",
]

PYTHON_PACKAGES = [
    "biopython==1.87",
    "fastexcel==0.20.2",
    "matplotlib==3.10.9",
    "numpy==2.4.6",
    "openpyxl==3.1.5",
    "pandas==3.0.3",
    "polars==1.40.1",
    "pysam==0.24.0",
    "xlsxwriter==3.2.9",
]

PYTHON_PACKAGE_NAMES = [package.split("==", 1)[0] for package in PYTHON_PACKAGES]

VERSION_ARGS = {
    "python": ["--version"],
    "minimap2": ["--version"],
    "samtools": ["--version"],
    "seqkit": ["version"],
    "blastn": ["-version"],
    "bcftools": ["--version"],
    "bamtools": ["--version"],
    "pigz": ["--version"],
    "flye": ["--version"],
    "simple_draft_asm": ["--version"],
    "gfa_editor_cli": ["--version"],
}


def normalize_path(value: str) -> str:
    path = Path(os.path.expandvars(os.path.expanduser(value)))
    if not path.is_absolute():
        path = Path.cwd() / path
    return str(path)


def run(command: list[str], dry_run: bool = False, cwd: Path | None = None) -> None:
    print("+ " + " ".join(command))
    if not dry_run:
        subprocess.run(command, cwd=cwd, check=True)


def venv_python(venv_dir: Path) -> Path:
    return venv_dir / ("Scripts/python.exe" if os.name == "nt" else "bin/python")


def executable_exists(path: str) -> bool:
    expanded = Path(os.path.expandvars(os.path.expanduser(path)))
    if str(expanded).endswith(".py"):
        return expanded.exists()
    return expanded.exists() and os.access(expanded, os.X_OK)


def first_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def sibling_project(name: str) -> Path:
    return ROOT.parent.parent / name


def detect_python(args: argparse.Namespace) -> str:
    if "python" in args.tool_overrides:
        return args.tool_overrides["python"]
    if args.python:
        return normalize_path(args.python)
    if args.venv and venv_python(args.venv).exists():
        return str(venv_python(args.venv))
    if shutil.which("python3.11"):
        return str(shutil.which("python3.11"))
    return sys.executable


def detect_simple_draft_asm(args: argparse.Namespace) -> str | None:
    if "simple_draft_asm" in args.tool_overrides:
        return args.tool_overrides["simple_draft_asm"]
    if args.simple_draft_asm:
        return normalize_path(args.simple_draft_asm)
    from_path = shutil.which("simple_draft_asm")
    if from_path:
        return from_path
    candidate = first_existing([
        sibling_project("simple_draft_asm") / "target" / "release" / "simple_draft_asm",
        ROOT.parent / "simple_draft_asm" / "target" / "release" / "simple_draft_asm",
    ])
    return str(candidate) if candidate else None


def detect_gfa_editor_cli(args: argparse.Namespace) -> str | None:
    if "gfa_editor_cli" in args.tool_overrides:
        return args.tool_overrides["gfa_editor_cli"]
    if args.gfa_editor_cli:
        return normalize_path(args.gfa_editor_cli)
    from_path = shutil.which("gfa_editor_cli")
    if from_path:
        return from_path
    candidate = first_existing([
        sibling_project("GFA_Editor") / "scripts" / "gfa_editor_cli.py",
        Path("/Applications/GFA_Editor.app/Contents/MacOS/GFA_Editor"),
    ])
    return str(candidate) if candidate else None


def detect_tools(args: argparse.Namespace) -> dict[str, str]:
    tools = {
        "python": detect_python(args),
        "simple_draft_asm": detect_simple_draft_asm(args),
        "gfa_editor_cli": detect_gfa_editor_cli(args),
    }
    command_names = {
        "minimap2": "minimap2",
        "samtools": "samtools",
        "seqkit": "seqkit",
        "blastn": "blastn",
        "bcftools": "bcftools",
        "bamtools": "bamtools",
        "pigz": "pigz",
        "flye": "flye",
    }
    for key, command in command_names.items():
        tools[key] = args.tool_overrides.get(key) or shutil.which(command)
    for key, value in args.tool_overrides.items():
        tools[key] = value
    return {key: value for key, value in tools.items() if value}


def write_soft_paths(path: Path, tools: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fout:
        print("# HiFiSR manual workflow soft_paths", file=fout)
        print("# Format: software_name<TAB>absolute_path_to_executable", file=fout)
        for tool in REQUIRED_TOOLS:
            if tool in tools:
                print(f"{tool}\t{tools[tool]}", file=fout)
        for tool in sorted(set(tools) - set(REQUIRED_TOOLS)):
            print(f"{tool}\t{tools[tool]}", file=fout)


def version_command(tool: str, path: str, python_path: str | None) -> list[str]:
    args = VERSION_ARGS.get(tool, ["--version"])
    if path.endswith(".py"):
        return [python_path or sys.executable, path, *args]
    return [path, *args]


def clean_version_text(text: str) -> str:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return ""
    return " | ".join(lines[:3])


def infer_simple_draft_asm_version(path: str) -> tuple[str, str]:
    exe_path = Path(path)
    for parent in exe_path.parents:
        cargo_toml = parent / "Cargo.toml"
        if not cargo_toml.exists():
            continue
        version = ""
        for line in cargo_toml.read_text().splitlines():
            stripped = line.strip()
            if stripped.startswith("version"):
                version = stripped.split("=", 1)[1].strip().strip('"')
                break
        git_rev = ""
        completed = subprocess.run(
            ["git", "-C", str(parent), "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if completed.returncode == 0:
            git_rev = completed.stdout.strip()
        parts = []
        if version:
            parts.append("Cargo.toml version " + version)
        if git_rev:
            parts.append("git " + git_rev)
        if parts:
            return "; ".join(parts), "inferred"
    return "", "unavailable"


def infer_gfa_editor_version(path: str) -> tuple[str, str]:
    cli_path = Path(path)
    for parent in cli_path.parents:
        backend_main = parent / "backend" / "main.py"
        if not backend_main.exists():
            continue
        version = ""
        for line in backend_main.read_text().splitlines():
            stripped = line.strip()
            if stripped.startswith("version="):
                version = stripped.split("=", 1)[1].strip().strip('",')
                break
        git_rev = ""
        completed = subprocess.run(
            ["git", "-C", str(parent), "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if completed.returncode == 0:
            git_rev = completed.stdout.strip()
        parts = []
        if version:
            parts.append("version " + version)
        if git_rev:
            parts.append("git " + git_rev)
        if parts:
            return "; ".join(parts), "inferred"
    return "", "unavailable"


def tool_version(tool: str, path: str, python_path: str | None) -> tuple[str, str, str]:
    command = version_command(tool, path, python_path)
    command_text = " ".join(command)
    if not executable_exists(path):
        return command_text, "", "path_missing_or_not_executable"
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=10,
        )
    except Exception as exc:
        return command_text, str(exc), "error"

    output = clean_version_text(completed.stdout + "\n" + completed.stderr)
    if completed.returncode == 0 and output:
        return command_text, output, "ok"

    if tool == "simple_draft_asm":
        inferred, status = infer_simple_draft_asm_version(path)
        if inferred:
            return command_text, inferred, status
    if tool == "gfa_editor_cli":
        inferred, status = infer_gfa_editor_version(path)
        if inferred:
            return command_text, inferred, status

    return command_text, output, "unavailable"


def write_soft_versions(path: Path, tools: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    python_path = tools.get("python")
    with path.open("w") as fout:
        print("software\tpath\tversion_command\tversion\tstatus", file=fout)
        for tool in REQUIRED_TOOLS:
            if tool not in tools:
                print(f"{tool}\t.\t.\t.\tmissing_from_soft_paths", file=fout)
                continue
            command, version, status = tool_version(tool, tools[tool], python_path)
            print(
                "\t".join([
                    tool,
                    tools[tool],
                    command,
                    version.replace("\t", " "),
                    status,
                ]),
                file=fout,
            )
        for tool in sorted(set(tools) - set(REQUIRED_TOOLS)):
            command, version, status = tool_version(tool, tools[tool], python_path)
            print(
                "\t".join([
                    tool,
                    tools[tool],
                    command,
                    version.replace("\t", " "),
                    status,
                ]),
                file=fout,
            )


def read_soft_paths(path: Path) -> dict[str, str]:
    tools = {}
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = line.split("\t", 1)
        if len(fields) != 2:
            raise ValueError(f"Invalid soft_paths line: {line!r}")
        tools[fields[0].strip()] = fields[1].strip()
    return tools


def check_soft_paths(path: Path) -> list[str]:
    tools = read_soft_paths(path)
    errors = []
    for tool in REQUIRED_TOOLS:
        value = tools.get(tool)
        if not value:
            errors.append(f"{tool}: missing from {path}")
        elif not executable_exists(value):
            errors.append(f"{tool}: path missing or not executable: {value}")
    return errors


def check_python_packages(python_path: str | None = None) -> list[str]:
    if python_path and Path(python_path).resolve() != Path(sys.executable).resolve():
        packages_repr = repr(PYTHON_PACKAGE_NAMES)
        code = f"""
import importlib.metadata as md

missing = []
for package in {packages_repr}:
    try:
        print(f"{{package}}=={{md.version(package)}}")
    except md.PackageNotFoundError:
        missing.append(package)

raise SystemExit("\\n".join(missing) if missing else 0)
"""
        completed = subprocess.run([python_path, "-c", code], capture_output=True, text=True)
        if completed.stdout:
            print(completed.stdout.rstrip())
        if completed.returncode != 0:
            missing = completed.stderr.strip() or completed.stdout.strip()
            return [f"python packages missing: {missing}"]
        return []

    errors = []
    for package in PYTHON_PACKAGE_NAMES:
        try:
            print(f"{package}=={metadata.version(package)}")
        except metadata.PackageNotFoundError:
            errors.append(f"python package missing: {package}")
    return errors


def install_tools(args: argparse.Namespace) -> None:
    system = platform.system()
    if system == "Darwin":
        run(["brew", "install", *BREW_PACKAGES], dry_run=args.dry_run)
        return
    if system == "Linux":
        manager = shutil.which("mamba") or shutil.which("conda")
        if not manager:
            raise RuntimeError("Linux install needs mamba or conda on PATH.")
        run(
            [
                manager,
                "create",
                "-y",
                "-n",
                args.conda_env,
                "-c",
                "conda-forge",
                "-c",
                "bioconda",
                *CONDA_PACKAGES,
            ],
            dry_run=args.dry_run,
        )
        print(f"Activate with: conda activate {args.conda_env}")
        return
    raise RuntimeError(f"Unsupported platform for automatic install: {system}")


def create_venv(args: argparse.Namespace) -> None:
    python = args.python or shutil.which("python3.11") or sys.executable
    run([python, "-m", "venv", str(args.venv)], dry_run=args.dry_run)


def install_python_packages(args: argparse.Namespace) -> None:
    python = detect_python(args)
    run([python, "-m", "pip", "install", "-U", "pip"], dry_run=args.dry_run)
    run([python, "-m", "pip", "install", *PYTHON_PACKAGES], dry_run=args.dry_run)


def build_simple_draft_asm(args: argparse.Namespace) -> None:
    root = Path(
        normalize_path(args.simple_draft_asm_root)
        if args.simple_draft_asm_root
        else sibling_project("simple_draft_asm")
    )
    if not root.exists():
        raise RuntimeError(f"simple_draft_asm root not found: {root}")
    run(["cargo", "build", "--release"], dry_run=args.dry_run, cwd=root)


def parse_tool_overrides(values: list[str], parser: argparse.ArgumentParser) -> dict[str, str]:
    overrides = {}
    for value in values:
        if "=" not in value:
            parser.error(f"--tool must use name=/absolute/path format: {value!r}")
        name, path = value.split("=", 1)
        name = name.strip()
        path = path.strip()
        if not name or not path:
            parser.error(f"--tool must use name=/absolute/path format: {value!r}")
        overrides[name] = normalize_path(path)
    return overrides


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--install-tools", action="store_true", help="Install OS/Conda command-line tools.")
    parser.add_argument("--create-venv", action="store_true", help="Create a Python virtual environment.")
    parser.add_argument("--install-python-packages", action="store_true", help="Install Python packages for the manual workflow.")
    parser.add_argument("--build-simple-draft-asm", action="store_true", help="Run cargo build --release for simple_draft_asm.")
    parser.add_argument("--write-soft-paths", action="store_true", help="Write the soft_paths file.")
    parser.add_argument("--check", action="store_true", help="Check soft_paths and Python packages.")
    parser.add_argument("--dry-run", action="store_true", help="Print install/build commands without executing them.")
    parser.add_argument("--soft-paths", default="deps/soft_paths_manual.txt", help="soft_paths output/check file.")
    parser.add_argument("--soft-versions", default="", help="soft_versions output file. Defaults to soft_paths directory/soft_versions.txt.")
    parser.add_argument("--venv", type=Path, default=ROOT / ".venv", help="Python virtual environment path.")
    parser.add_argument("--python", default="", help="Python executable to use.")
    parser.add_argument("--conda-env", default="hifisr", help="Conda/Mamba environment name for Linux install.")
    parser.add_argument("--simple-draft-asm-root", default="", help="Path to the simple_draft_asm source checkout.")
    parser.add_argument("--simple-draft-asm", default="", help="Path to the simple_draft_asm executable.")
    parser.add_argument("--gfa-editor-cli", default="", help="Path to GFA_Editor CLI script or executable.")
    parser.add_argument(
        "--tool",
        action="append",
        default=[],
        metavar="NAME=PATH",
        help=(
            "Manually set one software path in soft_paths. Can be repeated, "
            "for example --tool python=/abs/python --tool minimap2=/abs/minimap2."
        ),
    )
    args = parser.parse_args(argv)
    args.tool_overrides = parse_tool_overrides(args.tool, parser)

    if not any([
        args.install_tools,
        args.create_venv,
        args.install_python_packages,
        args.build_simple_draft_asm,
        args.write_soft_paths,
        args.check,
    ]):
        args.write_soft_paths = True
        args.check = True

    args.soft_paths = Path(args.soft_paths)
    if not args.soft_paths.is_absolute():
        args.soft_paths = ROOT / args.soft_paths
    args.soft_versions = Path(args.soft_versions) if args.soft_versions else args.soft_paths.with_name("soft_versions.txt")
    if not args.soft_versions.is_absolute():
        args.soft_versions = ROOT / args.soft_versions
    if not args.venv.is_absolute():
        args.venv = ROOT / args.venv
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    if args.install_tools:
        install_tools(args)
    if args.create_venv:
        create_venv(args)
    if args.install_python_packages:
        install_python_packages(args)
    if args.build_simple_draft_asm:
        build_simple_draft_asm(args)

    if args.write_soft_paths:
        tools = detect_tools(args)
        write_soft_paths(args.soft_paths, tools)
        print(f"Wrote {args.soft_paths}")
        write_soft_versions(args.soft_versions, tools)
        print(f"Wrote {args.soft_versions}")
        missing = [tool for tool in REQUIRED_TOOLS if tool not in tools]
        if missing:
            print("Missing from generated soft_paths: " + ", ".join(missing))

    if args.check:
        errors = []
        errors.extend(check_soft_paths(args.soft_paths))
        python_path = read_soft_paths(args.soft_paths).get("python") if args.soft_paths.exists() else None
        errors.extend(check_python_packages(python_path))
        if errors:
            print("Runtime check failed:", file=sys.stderr)
            for error in errors:
                print("  - " + error, file=sys.stderr)
            return 1
        print("Runtime check passed.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
