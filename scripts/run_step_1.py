#!/usr/bin/env python3
"""Step 1: extract organelle reads and build draft GFA files."""

from __future__ import annotations

import argparse
from pathlib import Path

from manual_workflow import (
    add_common_args,
    load_context,
    print_next_files,
    run_logged,
    sample_path,
    script_path,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    parser.add_argument("--mito-modes", default="ms,mh,mx")
    parser.add_argument("--plastid-modes", default="ps,ph")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    values = load_context(args)
    py = values["python"]

    run_logged(
        values,
        "extract_mtpt_reads.log",
        [
            py,
            script_path(values, "get_mtpt_reads.py"),
            values["soft_paths"],
            values["sample"],
            values["mito_ref"],
            values["plastid_ref"],
            values["raw_reads"],
            values["threads"],
            values["mito_read_limit"],
            values["plastid_read_limit"],
        ],
        dry_run=args.dry_run,
    )
    run_logged(
        values,
        "draft_assembly.mito.log",
        [
            py,
            script_path(values, "get_draft_assembly.py"),
            values["soft_paths"],
            values["sample"],
            "mito",
            values["mito_ref"],
            sample_path(values, "reads/mito.fastq.gz"),
            values["threads"],
            args.mito_modes,
            sample_path(values, "reads/mito.fastq.gz"),
        ],
        dry_run=args.dry_run,
    )
    run_logged(
        values,
        "draft_assembly.plastid.log",
        [
            py,
            script_path(values, "get_draft_assembly.py"),
            values["soft_paths"],
            values["sample"],
            "plastid",
            values["plastid_ref"],
            sample_path(values, "reads/plastid.fastq.gz"),
            values["threads"],
            args.plastid_modes,
            sample_path(values, "reads/plastid.fastq.gz"),
        ],
        dry_run=args.dry_run,
    )

    print_next_files(
        "Manual checkpoint 1: provide checked draft GFA files before step 2.",
        [
            sample_path(values, "draft_assembly/mito/mito_checked_draft.gfa"),
            sample_path(values, "draft_assembly/plastid/plastid_checked_draft.gfa"),
        ],
    )
    print("\nStep 1 finished for " + values["sample"])
    print("Logs: " + str(Path(values["manual_log_dir"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
