# HiFiSR module guide:
# - base: command, file, and soft_paths helpers; import hifisr_functions.base as hfbase
# - reads: read extraction, filtering, sampling, and correction; import hifisr_functions.reads as hfreads
# - references: reference rotation, assembly, polishing, and alignment; import hifisr_functions.references as hfref
# - variants: read-variant calling, grouping, and frequency analysis; import hifisr_functions.variants as hfvar
# - transfer: organelle/nuclear transfer-fragment analysis; import hifisr_functions.transfer as hftrans
# - annotations: annotation tables and feature-level summaries; import hifisr_functions.annotations as hfanno
# - reports: read statistics, plots, Excel tables, and report outputs; import hifisr_functions.reports as hfrps

import os
import json
import re
import shutil
import subprocess
import sys

REQUIRED_WORKFLOW_TOOLS = [
    "python",
    "minimap2",
    "samtools",
    "seqkit",
    "mecat",
    "blastn",
    "bcftools",
    "bamtools",
    "pigz",
    "bandage",
    "hifiasm",
    "flye",
    "canu",
    "simple_draft_asm",
]

PYTHON_IMPORT_NAMES = {
    "biopython": "Bio",
    "pulp": "pulp",
    "xlsxwriter": "xlsxwriter",
    "snakemake": "snakemake.cli",
}

# Function purity marker. "pure" means deterministic from explicit inputs with
# no file, shell, environment, logging, or input-mutation side effects.
FUNCTION_PURITY = {
    "run_checked": "impure",
    "get_cli_output_lines": "impure",
    "get_file_lines": "impure",
    "load_soft_paths": "impure",
    "validate_soft_paths": "impure",
    "require_soft_paths": "impure",
    "parse_python_requirements": "impure",
    "check_python_packages": "impure",
}


def run_checked(commands):
    print("+ " + commands, file=sys.stderr)
    completed = subprocess.run(commands, shell=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {completed.returncode}: {commands}"
        )
    return completed.returncode


def get_cli_output_lines(commands, side_effect=False):
    if side_effect:
        ret = subprocess.call(commands, shell=True)
        return ret

    completed = subprocess.run(commands, capture_output=True, shell=True)
    output_lines = completed.stdout.decode().split("\n")[0:-1]
    return output_lines


def get_file_lines(file):
    with open(file) as fin:
        lines = [line.rstrip("\n") for line in fin.readlines()]
    return lines


def load_soft_paths(soft_paths_file, validate=False, required_tools=None, print_paths=True):
    soft_paths_dict = {}
    errors = []
    for line_number, line in enumerate(get_file_lines(soft_paths_file), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        fields = line.split("\t", 1)
        if len(fields) < 2:
            errors.append(
                f"Invalid soft_paths line {line_number} in {soft_paths_file}: {line!r}"
            )
            continue

        name = fields[0].strip()
        executable = os.path.expandvars(os.path.expanduser(fields[1].strip()))
        if not name:
            errors.append(f"Missing software name on line {line_number} in {soft_paths_file}")
            continue
        if not executable:
            errors.append(f"Missing executable path for {name} on line {line_number}")
            continue
        if name in soft_paths_dict:
            errors.append(f"Duplicate soft_paths entry for {name} on line {line_number}")
            continue
        soft_paths_dict[name] = executable

    if errors:
        raise ValueError(
            "Invalid soft_paths.txt entries:\n"
            + "\n".join("  - " + error for error in errors)
        )

    if validate:
        validate_soft_paths(soft_paths_dict, required_tools, soft_paths_file)

    if print_paths:
        print("", file=sys.stderr)
        print("Loaded paths:", file=sys.stderr)
        for key in soft_paths_dict:
            print(key + " -> " + soft_paths_dict[key], file=sys.stderr)
        print("", file=sys.stderr)
    return soft_paths_dict


def validate_soft_paths(
    soft_paths_dict,
    required_tools=None,
    soft_paths_file="soft_paths.txt",
    check_all=True,
):
    errors = []
    tools_to_check = list(soft_paths_dict) if check_all else []
    for tool in required_tools or []:
        if tool not in tools_to_check:
            tools_to_check.append(tool)

    for tool in tools_to_check:
        executable = soft_paths_dict.get(tool)
        if not executable:
            errors.append(f"{tool}: missing from {soft_paths_file}")
            continue

        if os.path.isabs(executable):
            if not os.path.exists(executable):
                errors.append(f"{tool}: path does not exist: {executable}")
            elif not os.access(executable, os.X_OK):
                errors.append(f"{tool}: path is not executable: {executable}")
        elif shutil.which(executable) is None:
            errors.append(f"{tool}: command not found on PATH: {executable}")

    if errors:
        raise RuntimeError(
            "Invalid soft_paths.txt entries:\n"
            + "\n".join("  - " + error for error in errors)
        )
    return True


def require_soft_paths(soft_paths_dict, required_tools):
    return validate_soft_paths(soft_paths_dict, required_tools, check_all=False)


def parse_python_requirements(requirements_file):
    requirements = []
    requirement_re = re.compile(r"^([A-Za-z0-9_.-]+)\s*(?:==\s*([^;\s]+))?")
    for line_number, line in enumerate(get_file_lines(requirements_file), start=1):
        stripped = line.split("#", 1)[0].strip()
        if not stripped or stripped.startswith("-"):
            continue
        stripped = stripped.split(";", 1)[0].strip()
        match = requirement_re.match(stripped)
        if not match:
            raise ValueError(
                f"Unsupported requirement line {line_number} in {requirements_file}: {line!r}"
            )
        requirements.append({
            "name": match.group(1),
            "version": match.group(2),
        })
    return requirements


def check_python_packages(requirements_file, python_executable=None):
    requirements = parse_python_requirements(requirements_file)
    python_executable = python_executable or sys.executable
    validate_soft_paths({"python": python_executable}, ["python"], "python executable")

    check_code = r"""
    import importlib
    import importlib.metadata as metadata
    import json
    import subprocess
    import sys

    requirements = json.loads(sys.argv[1])
    import_names = json.loads(sys.argv[2])

    missing = []
    version_errors = []
    import_errors = []
    installed = []

    for requirement in requirements:
        name = requirement["name"]
        normalized = name.lower().replace("_", "-")
        try:
            version = metadata.version(name)
        except metadata.PackageNotFoundError:
            missing.append(name)
            continue

        expected_version = requirement.get("version")
        if expected_version and version != expected_version:
            version_errors.append(f"{name}: expected {expected_version}, found {version}")

        import_name = import_names.get(normalized, name.replace("-", "_"))
        try:
            importlib.import_module(import_name)
        except Exception as exc:
            import_errors.append(
                f"{name}: import {import_name} failed: {exc.__class__.__name__}: {exc}"
            )
        installed.append({"name": name, "version": version})

    standard_module_errors = []
    for module_name in ["sqlite3"]:
        try:
            importlib.import_module(module_name)
        except Exception as exc:
            standard_module_errors.append(
                f"{module_name}: import failed: {exc.__class__.__name__}: {exc}"
            )

    snakemake_error = None
    snakemake_required = any(
        requirement["name"].lower() == "snakemake" for requirement in requirements
    )
    if snakemake_required and not any(item["name"].lower() == "snakemake" for item in installed):
        snakemake_error = "snakemake: package is not installed"
    elif snakemake_required:
        completed = subprocess.run(
            [sys.executable, "-m", "snakemake", "--version"],
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            snakemake_error = (
                "snakemake --version failed: "
                + (completed.stderr.strip() or completed.stdout.strip())
            )

    report = {
        "python": sys.executable,
        "packages_checked": len(requirements),
        "installed": installed,
        "missing": missing,
        "version_errors": version_errors,
        "import_errors": import_errors,
        "standard_module_errors": standard_module_errors,
        "snakemake_error": snakemake_error,
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    if missing or version_errors or import_errors or standard_module_errors or snakemake_error:
        sys.exit(1)
    """
    completed = subprocess.run(
        [
            python_executable,
            "-c",
            check_code,
            json.dumps(requirements),
            json.dumps(PYTHON_IMPORT_NAMES),
        ],
        capture_output=True,
        text=True,
    )
    try:
        report = json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError(
            "Could not parse Python dependency check output:\n"
            + completed.stdout
            + completed.stderr
        ) from exc

    if completed.returncode != 0:
        messages = []
        for key in [
            "missing",
            "version_errors",
            "import_errors",
            "standard_module_errors",
        ]:
            messages.extend(report.get(key) or [])
        if report.get("snakemake_error"):
            messages.append(report["snakemake_error"])
        raise RuntimeError(
            "Invalid Python environment:\n"
            + "\n".join("  - " + message for message in messages)
        )
    return report
