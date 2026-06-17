# HiFiSR module guide:
# - base: command, file, and soft_paths helpers; import hifisr_functions.base as hfbase
# - reads: read extraction, filtering, sampling, and correction; import hifisr_functions.reads as hfreads
# - references: reference rotation, assembly, polishing, and alignment; import hifisr_functions.references as hfref
# - variants: read-variant calling, grouping, and frequency analysis; import hifisr_functions.variants as hfvar
# - transfer: organelle/nuclear transfer-fragment analysis; import hifisr_functions.transfer as hftrans
# - annotations: annotation tables and feature-level summaries; import hifisr_functions.annotations as hfanno
# - reports: read statistics, plots, Excel tables, and report outputs; import hifisr_functions.reports as hfrps

import argparse
import os
import sys

import _bootstrap  # noqa: F401
import hifisr_functions.base as hfbase
from hifisr_functions.graph.verified_gfa import build_verified_gfa


def absolute_path(value, label):
    if not os.path.isabs(value):
        print(label + " must be an absolute path: " + value, file=sys.stderr)
        sys.exit(1)
    return value


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build a reviewable verified GFA by mapping merged raw GFA nodes "
            "onto a run_3 verified FASTA. No repeat resolution is performed."
        )
    )
    parser.add_argument("soft_paths_file")
    parser.add_argument("genome")
    parser.add_argument("raw_gfa")
    parser.add_argument("verified_fasta")
    parser.add_argument("run_dir")
    parser.add_argument("output_dir")
    parser.add_argument("threads", nargs="?", default="1")
    parser.add_argument(
        "--merged-gfa-template",
        default=None,
        help=(
            "Optional curated merged raw GFA for manual comparison. If omitted, "
            "the script automatically judges whether raw_gfa needs non-repeat "
            "linear compaction and writes merged_{genome}_raw.gfa."
        ),
    )
    parser.add_argument(
        "--max-link-gap",
        type=int,
        default=1000,
        help="Maximum verified-coordinate gap for a link_candidate to be medium confidence.",
    )
    parser.add_argument(
        "--flank-anchor-length",
        type=int,
        default=500,
        help="Anchor length inside each single-copy flank endpoint for read-level repeat path support.",
    )
    parser.add_argument(
        "--repeat-anchor-length",
        type=int,
        default=500,
        help="Anchor length inside each repeat copy core for read-level repeat path support.",
    )
    parser.add_argument(
        "--flank-anchor-offset",
        type=int,
        default=1000,
        help="Distance from the repeat-adjacent flank endpoint before placing the flank anchor.",
    )
    parser.add_argument(
        "--min-anchor-overlap",
        type=int,
        default=200,
        help="Minimum aligned overlap with each flank/repeat anchor required to count one FL read.",
    )
    parser.add_argument(
        "--gfa-editor-cli",
        default=None,
        help=(
            "Optional GFA_Editor gfa_editor_cli.py path. If omitted, the script "
            "uses soft_paths['gfa_editor_cli'] when available, then a local "
            "Codex/GFA_Editor fallback."
        ),
    )
    parser.add_argument(
        "--image-reference-fasta",
        default=None,
        help=(
            "Reference FASTA used for GFA PDF query colouring. Use the initial "
            "rotated workflow reference, not the corrected/verified FASTA, when "
            "reviewing against the original coordinate system."
        ),
    )
    return parser.parse_args()


args = parse_args()
soft_paths_dict = hfbase.load_soft_paths(args.soft_paths_file)

raw_gfa = absolute_path(args.raw_gfa, "raw_gfa")
verified_fasta = absolute_path(args.verified_fasta, "verified_fasta")
run_dir = absolute_path(args.run_dir, "run_dir")
output_dir = absolute_path(args.output_dir, "output_dir")
merged_gfa_template = args.merged_gfa_template
if merged_gfa_template is not None:
    merged_gfa_template = absolute_path(merged_gfa_template, "merged_gfa_template")
image_reference_fasta = args.image_reference_fasta
if image_reference_fasta is not None:
    image_reference_fasta = absolute_path(image_reference_fasta, "image_reference_fasta")

outputs = build_verified_gfa(
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

print("Verified GFA review outputs:")
for key, value in outputs.items():
    print(key + "\t" + str(value))
