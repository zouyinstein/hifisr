#!/usr/bin/env python3
"""Step 3: apply mito manual edits, run run_3, rebuild verified GFA, and report."""

from __future__ import annotations

import argparse
from pathlib import Path

from manual_workflow import (
    add_common_args,
    load_context,
    require_files,
    run_logged,
    sample_path,
    script_path,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    return parser.parse_args()


def print_file_if_exists(path: str) -> None:
    target = Path(path)
    if not target.exists():
        print("[MISSING] " + path)
        return
    print("\n# " + path)
    print(target.read_text().rstrip())


def main() -> int:
    args = parse_args()
    values = load_context(args)
    py = values["python"]

    pos_ref_alt = sample_path(values, "draft_assembly/mito/pos_ref_alt.txt")
    if not args.dry_run:
        require_files([pos_ref_alt], allow_empty=True)

    run_logged(
        values,
        "correct_reference_with_manual_edits.mito.log",
        [
            py,
            script_path(values, "correct_erroneous.py"),
            sample_path(values, "draft_assembly/mito/mito_flye_polish_aligned.fasta"),
            sample_path(values, "draft_assembly/mito/mito_flye_polish_aligned_cor.fasta"),
            pos_ref_alt,
        ],
        dry_run=args.dry_run,
    )
    run_logged(
        values,
        "run_3.mito.log",
        [
            py,
            script_path(values, "get_variants_in_reads.py"),
            values["soft_paths"],
            values["sample"],
            "mito",
            "run_3",
            sample_path(values, "draft_assembly/mito/mito_flye_polish_aligned_cor.fasta"),
            sample_path(values, "reads/mito.fastq.gz"),
            values["threads"],
        ],
        dry_run=args.dry_run,
    )
    run_logged(
        values,
        "verified_gfa.mito.log",
        [
            py,
            script_path(values, "build_verified_gfa_read_support.py"),
            values["soft_paths"],
            values["sample"],
            "mito",
            sample_path(values, "draft_assembly/mito/mito_checked_draft.gfa"),
            sample_path(values, "draft_assembly/mito/mito_flye_polish_aligned_cor.fasta"),
            sample_path(values, "mito/run_3"),
            values["mito_ref"],
            values["threads"],
        ],
        dry_run=args.dry_run,
    )
    run_logged(
        values,
        "verified_gfa.plastid.log",
        [
            py,
            script_path(values, "build_verified_gfa_read_support.py"),
            values["soft_paths"],
            values["sample"],
            "plastid",
            sample_path(values, "draft_assembly/plastid/plastid_checked_draft.gfa"),
            sample_path(values, "draft_assembly/plastid/plastid_flye_polish_aligned.fasta"),
            sample_path(values, "plastid/run_2"),
            values["plastid_ref"],
            values["threads"],
        ],
        dry_run=args.dry_run,
    )
    run_logged(
        values,
        "generate_compact_report_bundle.log",
        [
            py,
            script_path(values, "generate_compact_report_bundle.py"),
            "--base",
            values["sample_dir"],
        ],
        dry_run=args.dry_run,
    )

    report_files = [
        sample_path(values, "reports/index.html"),
        sample_path(values, "reports/project_manifest.json"),
        sample_path(values, "reports/mito/index.html"),
        sample_path(values, "reports/plastid/index.html"),
    ]
    consistency_files = [
        sample_path(values, "draft_assembly/mito/verified_gfa_read_support/verified_mito.auto_repeat_check.tsv"),
        sample_path(values, "draft_assembly/plastid/verified_gfa_read_support/verified_plastid.auto_repeat_check.tsv"),
        sample_path(values, "draft_assembly/mito/verified_gfa_read_support/coordinate_consistency.tsv"),
        sample_path(values, "draft_assembly/plastid/verified_gfa_read_support/coordinate_consistency.tsv"),
    ]
    if not args.dry_run:
        require_files(report_files + consistency_files)
        for path in consistency_files:
            print_file_if_exists(path)

    print("\nStep 3 finished for " + values["sample"])
    print("Reports: " + sample_path(values, "reports/index.html"))
    print("Logs: " + str(Path(values["manual_log_dir"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
