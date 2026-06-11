#!/usr/bin/env python3
# HiFiSR module guide:
# - base: command, file, and soft_paths helpers; import hifisr_functions.base as hfbase
# - reads: read extraction, filtering, sampling, and correction; import hifisr_functions.reads as hfreads
# - references: reference rotation, assembly, polishing, and alignment; import hifisr_functions.references as hfref
# - variants: read-variant calling, grouping, and frequency analysis; import hifisr_functions.variants as hfvar
# - transfer: organelle/nuclear transfer-fragment analysis; import hifisr_functions.transfer as hftrans
# - annotations: annotation tables and feature-level summaries; import hifisr_functions.annotations as hfanno
# - reports: read statistics, plots, Excel tables, and report outputs; import hifisr_functions.reports as hfrps

"""Print planned batch inputs from a hifisr sample sheet."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SAMPLE_SHEET = ROOT / "tests" / "fixtures" / "sample_sheet.tsv"


def iter_samples(sample_sheet: Path):
    with sample_sheet.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample = (row.get("sample") or "").strip()
            if not sample:
                continue
            yield row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "sample_sheet",
        nargs="?",
        default=DEFAULT_SAMPLE_SHEET,
        type=Path,
        help="TSV file with hifisr batch sample metadata.",
    )
    args = parser.parse_args()

    sample_sheet = args.sample_sheet
    if not sample_sheet.is_absolute():
        sample_sheet = ROOT / sample_sheet

    for row in iter_samples(sample_sheet):
        fields = [
            f"sample={row.get('sample', '')}",
            f"run_info={row.get('run_info', '')}",
            f"mito_ref={row.get('mito_ref', '')}",
            f"plastid_ref={row.get('plastid_ref', '')}",
            f"reads={row.get('reads', '')}",
        ]
        print("\t".join(fields))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

