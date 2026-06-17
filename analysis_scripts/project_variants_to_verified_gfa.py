#!/usr/bin/env python3
"""Split verified FASTA and run variants by verified GFA node coordinates."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hifisr_functions.graph import evidence_projection as hfgproj  # noqa: E402


def parse_args():
    parser = argparse.ArgumentParser(
        description="Project SNV/InDel and SV-like run evidence onto verified GFA nodes."
    )
    parser.add_argument("--verified-gfa", required=True)
    parser.add_argument("--verified-fasta", required=True)
    parser.add_argument("--node-source-map", required=True)
    parser.add_argument("--linear-fasta", required=True)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--genome", default=None)
    parser.add_argument(
        "--variant-table",
        action="append",
        default=None,
        help="Variant Excel table to project. Defaults to the two run_3 variants_anno tables.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    run_dir = Path(args.run_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    variant_tables = args.variant_table
    if variant_tables is None:
        variant_tables = [
            str(run_dir / "variants_anno_combined_depth_frq.xlsx"),
            str(run_dir / "variants_anno_combined_depth_frq_filter.xlsx"),
        ]
    genome = args.genome
    if genome is None:
        stem = Path(args.verified_gfa).stem
        genome = stem.removeprefix("verified_")

    coordinate_check = hfgproj.compare_fasta_coordinates(
        args.linear_fasta,
        args.verified_fasta,
    )
    if not coordinate_check["identical_sequence"]:
        raise RuntimeError(
            "linear FASTA and verified FASTA are not coordinate-identical; "
            "refuse direct projection. See coordinate_consistency.tsv inputs."
        )

    intervals = hfgproj.build_node_projection_intervals(
        args.node_source_map,
        args.verified_gfa,
        args.verified_fasta,
    )

    by_nodes_fasta = output_dir / ("verified_" + genome + "_by_nodes.fasta")
    linear_to_node_map = output_dir / ("verified_" + genome + "_linear_to_node_coordinate.tsv")
    hfgproj.write_verified_gfa_by_nodes_fasta(
        by_nodes_fasta,
        args.verified_gfa,
        args.node_source_map,
    )
    hfgproj.write_linear_to_node_coordinate_map(linear_to_node_map, intervals)
    variant_outputs = hfgproj.write_simple_variant_node_outputs(
        variant_tables,
        intervals,
        args.node_source_map,
        output_dir / "variant_by_nodes",
    )

    coordinate_tsv = output_dir / "coordinate_consistency.tsv"
    pd.DataFrame([coordinate_check]).to_csv(coordinate_tsv, sep="\t", index=False)

    print("by_nodes_fasta\t" + str(by_nodes_fasta))
    print("linear_to_node_coordinate\t" + str(linear_to_node_map))
    print("coordinate_consistency\t" + str(coordinate_tsv))
    for key, value in variant_outputs.items():
        print(key + "\t" + str(value))


if __name__ == "__main__":
    main()
