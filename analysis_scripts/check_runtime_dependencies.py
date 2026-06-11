# HiFiSR module guide:
# - base: command, file, and soft_paths helpers; import hifisr_functions.base as hfbase
# - reads: read extraction, filtering, sampling, and correction; import hifisr_functions.reads as hfreads
# - references: reference rotation, assembly, polishing, and alignment; import hifisr_functions.references as hfref
# - variants: read-variant calling, grouping, and frequency analysis; import hifisr_functions.variants as hfvar
# - transfer: organelle/nuclear transfer-fragment analysis; import hifisr_functions.transfer as hftrans
# - annotations: annotation tables and feature-level summaries; import hifisr_functions.annotations as hfanno
# - reports: read statistics, plots, Excel tables, and report outputs; import hifisr_functions.reports as hfrps

import argparse
from pathlib import Path
import sys

import _bootstrap  # noqa: F401
import hifisr_functions.base as hfbase


REPO_ROOT = Path(__file__).resolve().parents[1]


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Check soft_paths executables and Python package dependencies."
    )
    parser.add_argument(
        "soft_paths",
        help="Path to soft_paths.txt or another tab-delimited software path file.",
    )
    parser.add_argument(
        "--requirements",
        default=str(REPO_ROOT / "requirements-dev.txt"),
        help="Python requirements file to validate.",
    )
    parser.add_argument(
        "--python",
        default=None,
        help="Python executable to check. Defaults to the python entry in soft_paths.",
    )
    parser.add_argument(
        "--skip-python-packages",
        action="store_true",
        help="Only validate soft_paths executables.",
    )
    args = parser.parse_args(argv)

    try:
        soft_paths = hfbase.load_soft_paths(
            args.soft_paths,
            validate=True,
            required_tools=hfbase.REQUIRED_WORKFLOW_TOOLS,
        )
        print(f"soft_paths check passed: {args.soft_paths}")

        if not args.skip_python_packages:
            python_executable = args.python or soft_paths.get("python") or sys.executable
            report = hfbase.check_python_packages(args.requirements, python_executable)
            print(f"Python executable check passed: {report['python']}")
            print(f"Python package check passed: {report['packages_checked']} packages")
        return 0
    except Exception as exc:
        print("Dependency check failed:", file=sys.stderr)
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
