import os
import shutil
import subprocess
import sys


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


def load_soft_paths(soft_paths_file):
    soft_paths_dict = {}
    for line_number, line in enumerate(get_file_lines(soft_paths_file), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        fields = line.split("\t")
        if len(fields) < 2:
            raise ValueError(
                f"Invalid soft_paths line {line_number} in {soft_paths_file}: {line!r}"
            )
        soft_paths_dict[fields[0]] = fields[1]

    print("", file=sys.stderr)
    print("Loaded paths:", file=sys.stderr)
    for key in soft_paths_dict:
        print(key + " -> " + soft_paths_dict[key], file=sys.stderr)
    print("", file=sys.stderr)
    return soft_paths_dict


def require_soft_paths(soft_paths_dict, required_tools):
    errors = []
    for tool in required_tools:
        executable = soft_paths_dict.get(tool)
        if not executable:
            errors.append(f"{tool}: missing from soft_paths.txt")
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
