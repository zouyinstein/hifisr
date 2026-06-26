#!/usr/bin/env python3
"""Step 2: polish checked GFA files, align to references, and run run_2."""

from __future__ import annotations

import argparse
from pathlib import Path

from manual_workflow import (
    add_common_args,
    load_context,
    print_next_files,
    require_files,
    run_logged,
    sample_path,
    script_path,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    values = load_context(args)
    py = values["python"]

    checked_mito_gfa = sample_path(values, "draft_assembly/mito/mito_checked_draft.gfa")
    checked_plastid_gfa = sample_path(values, "draft_assembly/plastid/plastid_checked_draft.gfa")
    if not args.dry_run:
        require_files([checked_mito_gfa, checked_plastid_gfa])

    run_logged(
        values,
        "polish_alignment.mito.log",
        [
            py,
            script_path(values, "get_polished_assembly.py"),
            values["soft_paths"],
            values["sample"],
            "mito",
            values["mito_ref"],
            checked_mito_gfa,
            sample_path(values, "reads/mito.fastq.gz"),
            values["threads"],
        ],
        dry_run=args.dry_run,
    )
    run_logged(
        values,
        "polish_alignment.plastid.log",
        [
            py,
            script_path(values, "get_polished_assembly.py"),
            values["soft_paths"],
            values["sample"],
            "plastid",
            values["plastid_ref"],
            checked_plastid_gfa,
            sample_path(values, "reads/plastid.fastq.gz"),
            values["threads"],
        ],
        dry_run=args.dry_run,
    )
    run_logged(
        values,
        "run_2.mito.log",
        [
            py,
            script_path(values, "get_variants_in_reads.py"),
            values["soft_paths"],
            values["sample"],
            "mito",
            "run_2",
            sample_path(values, "draft_assembly/mito/mito_flye_polish_aligned.fasta"),
            sample_path(values, "reads/mito.fastq.gz"),
            values["threads"],
        ],
        dry_run=args.dry_run,
    )
    run_logged(
        values,
        "run_2.plastid.log",
        [
            py,
            script_path(values, "get_variants_in_reads.py"),
            values["soft_paths"],
            values["sample"],
            "plastid",
            "run_2",
            sample_path(values, "draft_assembly/plastid/plastid_flye_polish_aligned.fasta"),
            sample_path(values, "reads/plastid.fastq.gz"),
            values["threads"],
        ],
        dry_run=args.dry_run,
    )

    print_next_files(
        "Manual checkpoint 2: inspect run_2 and provide mito pos_ref_alt.txt before step 3.",
        [
            sample_path(values, "mito/run_2/all_FL_report.txt"),
            sample_path(values, "mito/run_2/variants_anno_combined_depth_frq_filter.xlsx"),
            sample_path(values, "plastid/run_2/all_FL_report.txt"),
            sample_path(values, "plastid/run_2/variants_anno_combined_depth_frq_filter.xlsx"),
            sample_path(values, "draft_assembly/mito/pos_ref_alt.txt"),
        ],
    )
    print("\nStep 2 finished for " + values["sample"])
    print("Logs: " + str(Path(values["manual_log_dir"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
