#!/usr/bin/env python3
"""Build verified GFA read-support outputs and project variants by node."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import sys

import _bootstrap  # noqa: F401
import hifisr_functions.base as hfbase
from hifisr_functions.graph import evidence_projection as hfgproj
from hifisr_functions.graph.verified_gfa import build_verified_gfa


def absolute_path(value: str, label: str) -> str:
    if not os.path.isabs(value):
        raise ValueError(label + " must be an absolute path: " + value)
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Wrapper for build_verified_gfa.py + project_variants_to_verified_gfa.py. "
            "Runs one genome and writes verified GFA read-support outputs."
        )
    )
    parser.add_argument("soft_paths_file")
    parser.add_argument("sample_index")
    parser.add_argument("genome")
    parser.add_argument("raw_gfa")
    parser.add_argument("verified_fasta")
    parser.add_argument("run_dir")
    parser.add_argument("image_reference_fasta")
    parser.add_argument("threads", nargs="?", default="1")
    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Output directory. Defaults to raw_gfa parent / verified_gfa_read_support."
        ),
    )
    parser.add_argument(
        "--variant-table",
        action="append",
        default=None,
        help=(
            "Variant Excel table to project. Defaults to "
            "variants_anno_combined_depth_frq.xlsx and "
            "variants_anno_combined_depth_frq_filter.xlsx under run_dir."
        ),
    )
    parser.add_argument(
        "--merged-gfa-template",
        default=None,
        help="Optional curated merged raw GFA for manual comparison.",
    )
    parser.add_argument(
        "--max-link-gap",
        type=int,
        default=1000,
        help="Maximum verified-coordinate gap for a medium-confidence link candidate.",
    )
    parser.add_argument(
        "--flank-anchor-length",
        type=int,
        default=500,
        help="Anchor length inside each single-copy flank endpoint.",
    )
    parser.add_argument(
        "--repeat-anchor-length",
        type=int,
        default=500,
        help="Anchor length inside each repeat copy core.",
    )
    parser.add_argument(
        "--flank-anchor-offset",
        type=int,
        default=1000,
        help="Distance from repeat-adjacent flank endpoint before placing the anchor.",
    )
    parser.add_argument(
        "--min-anchor-overlap",
        type=int,
        default=200,
        help="Minimum aligned overlap required to count one FL read.",
    )
    parser.add_argument(
        "--gfa-editor-cli",
        default=None,
        help="Optional GFA_Editor CLI path; otherwise soft_paths['gfa_editor_cli'] is used.",
    )
    return parser.parse_args()


def existing_absolute_file(value: str, label: str) -> str:
    path = absolute_path(value, label)
    if not os.path.exists(path):
        raise FileNotFoundError(label + " does not exist: " + path)
    return path


def existing_absolute_dir(value: str, label: str) -> str:
    path = absolute_path(value, label)
    if not os.path.isdir(path):
        raise FileNotFoundError(label + " directory does not exist: " + path)
    return path


def write_one_row_tsv(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(row.keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def project_variants_to_verified_nodes(
    genome: str,
    verified_gfa: Path,
    verified_fasta: Path,
    node_source_map: Path,
    linear_fasta: str,
    run_dir: Path,
    output_dir: Path,
    variant_tables: list[str],
) -> dict[str, str]:
    coordinate_check = hfgproj.compare_fasta_coordinates(
        linear_fasta,
        str(verified_fasta),
    )
    if not coordinate_check["identical_sequence"]:
        raise RuntimeError(
            "linear FASTA and verified FASTA are not coordinate-identical; "
            "refuse direct projection. See coordinate_consistency.tsv inputs."
        )

    intervals = hfgproj.build_node_projection_intervals(
        str(node_source_map),
        str(verified_gfa),
        str(verified_fasta),
    )

    by_nodes_fasta = output_dir / ("verified_" + genome + "_by_nodes.fasta")
    linear_to_node_map = output_dir / ("verified_" + genome + "_linear_to_node_coordinate.tsv")
    hfgproj.write_verified_gfa_by_nodes_fasta(
        by_nodes_fasta,
        str(verified_gfa),
        str(node_source_map),
    )
    hfgproj.write_linear_to_node_coordinate_map(linear_to_node_map, intervals)
    variant_outputs = hfgproj.write_simple_variant_node_outputs(
        variant_tables,
        intervals,
        str(node_source_map),
        output_dir / "variant_by_nodes",
    )

    coordinate_tsv = output_dir / "coordinate_consistency.tsv"
    write_one_row_tsv(coordinate_tsv, coordinate_check)

    outputs = {
        "by_nodes_fasta": str(by_nodes_fasta),
        "linear_to_node_coordinate": str(linear_to_node_map),
        "coordinate_consistency": str(coordinate_tsv),
    }
    outputs.update({key: str(value) for key, value in variant_outputs.items()})
    return outputs


def main() -> int:
    args = parse_args()
    try:
        soft_paths_file = existing_absolute_file(args.soft_paths_file, "soft_paths_file")
        raw_gfa = existing_absolute_file(args.raw_gfa, "raw_gfa")
        verified_fasta = existing_absolute_file(args.verified_fasta, "verified_fasta")
        run_dir = existing_absolute_dir(args.run_dir, "run_dir")
        image_reference_fasta = existing_absolute_file(
            args.image_reference_fasta,
            "image_reference_fasta",
        )
        output_dir = args.output_dir
        if output_dir is None:
            output_dir = str(Path(raw_gfa).parent / "verified_gfa_read_support")
        output_dir = absolute_path(output_dir, "output_dir")
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        merged_gfa_template = args.merged_gfa_template
        if merged_gfa_template is not None:
            merged_gfa_template = existing_absolute_file(
                merged_gfa_template,
                "merged_gfa_template",
            )

        variant_tables = args.variant_table
        if variant_tables is None:
            variant_tables = [
                str(Path(run_dir) / "variants_anno_combined_depth_frq.xlsx"),
                str(Path(run_dir) / "variants_anno_combined_depth_frq_filter.xlsx"),
            ]
        variant_tables = [
            existing_absolute_file(path, "variant_table") for path in variant_tables
        ]

        soft_paths_dict = hfbase.load_soft_paths(soft_paths_file)
        build_outputs = build_verified_gfa(
            genome=args.genome,
            raw_gfa_path=raw_gfa,
            verified_fasta_path=verified_fasta,
            run_dir=run_dir,
            output_dir=output_dir,
            soft_paths_dict=soft_paths_dict,
            threads=args.threads,
            merged_gfa_template=merged_gfa_template,
            max_link_gap=args.max_link_gap,
            flank_anchor_length=args.flank_anchor_length,
            flank_anchor_offset=args.flank_anchor_offset,
            repeat_anchor_length=args.repeat_anchor_length,
            min_anchor_overlap=args.min_anchor_overlap,
            gfa_editor_cli=args.gfa_editor_cli,
            image_reference_fasta=image_reference_fasta,
        )

        output_path = Path(output_dir)
        projection_outputs = project_variants_to_verified_nodes(
            genome=args.genome,
            verified_gfa=output_path / ("verified_" + args.genome + ".gfa"),
            verified_fasta=output_path / ("verified_" + args.genome + ".fasta"),
            node_source_map=output_path / "verified_node_source_map.tsv",
            linear_fasta=verified_fasta,
            run_dir=Path(run_dir),
            output_dir=output_path,
            variant_tables=variant_tables,
        )

        print("sample\t" + args.sample_index)
        print("genome\t" + args.genome)
        print("output_dir\t" + output_dir)
        print("Verified GFA read-support outputs:")
        for key, value in build_outputs.items():
            print(key + "\t" + str(value))
        print("Variant projection outputs:")
        for key, value in projection_outputs.items():
            print(key + "\t" + str(value))
    except Exception as exc:
        print("build_verified_gfa_read_support.py: " + str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
