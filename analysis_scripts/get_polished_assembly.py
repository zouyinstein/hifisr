# HiFiSR module guide:
# - base: command, file, and soft_paths helpers; import hifisr_functions.base as hfbase
# - reads: read extraction, filtering, sampling, and correction; import hifisr_functions.reads as hfreads
# - references: reference rotation, assembly, polishing, and alignment; import hifisr_functions.references as hfref
# - variants: read-variant calling, grouping, and frequency analysis; import hifisr_functions.variants as hfvar
# - transfer: organelle/nuclear transfer-fragment analysis; import hifisr_functions.transfer as hftrans
# - annotations: annotation tables and feature-level summaries; import hifisr_functions.annotations as hfanno
# - reports: read statistics, plots, Excel tables, and report outputs; import hifisr_functions.reports as hfrps

import sys
import _bootstrap  # noqa: F401
import hifisr_functions.base as hfbase
import hifisr_functions.references as hfref
import hifisr_functions.reports as hfrps
from hifisr_functions.graph.verified_gfa import merge_unambiguous_gfa, parse_gfa, write_gfa
import os
import shutil
import subprocess
from pathlib import Path


def move_backup_info(genome):
    backup_dir = "backup_info"
    keep_prefixes = (genome + "_flye_polish", "simple_draft_asm")
    keep_filenames = {
        "gfa_image_export_protocol.tsv",
        "pos_ref_alt.txt",
        f"{genome}_checked_draft.gfa",
        f"{genome}_checked_draft.unmerged.gfa",
        f"{genome}_checked_draft.fasta",
        f"{genome}_checked_draft.pdf",
    }
    if not os.path.exists(backup_dir):
        os.makedirs(backup_dir)
    for filename in sorted(os.listdir(".")):
        if filename == backup_dir or filename in keep_filenames or filename.startswith(keep_prefixes):
            continue
        backup_path = os.path.join(backup_dir, filename)
        if os.path.isdir(backup_path) and not os.path.islink(backup_path):
            shutil.rmtree(backup_path)
        elif os.path.exists(backup_path):
            os.remove(backup_path)
        shutil.move(filename, backup_path)


def draft_input_kind(path):
    suffixes = "".join(Path(path).suffixes).lower()
    if suffixes.endswith(".gfa"):
        return "gfa"
    if suffixes.endswith(".fa") or suffixes.endswith(".fasta"):
        return "fasta"
    return "unknown"


def checked_gfa_state(path):
    name = Path(path).name.lower()
    if "unmerged" in name:
        return "declared_unmerged"
    if "merged" in name:
        return "declared_merged"
    return "auto_detect"


def write_prepare_report(path, rows):
    with open(path, "wt") as fout:
        print("metric\tvalue", file=fout)
        for key, value in rows:
            print(str(key) + "\t" + str(value), file=fout)


def run_command(command):
    print("+ " + " ".join(command), file=sys.stderr)
    completed = subprocess.run(command, capture_output=True, text=True)
    if completed.returncode != 0:
        raise RuntimeError(
            "Command failed with exit code "
            + str(completed.returncode)
            + ": "
            + " ".join(command)
            + "\n"
            + completed.stderr
        )
    return completed


def command_for_python_script(path, soft_paths_dict):
    path = str(path)
    if not path.endswith(".py"):
        return [path]
    script_path = Path(path)
    local_python = script_path.resolve().parents[1] / ".venv" / "bin" / "python"
    if local_python.exists():
        return [str(local_python), path]
    if not os.access(path, os.X_OK):
        return [soft_paths_dict.get("python", "python"), path]
    return [path]


def prepare_checked_gfa_for_polish(genome, reference_fasta, gfa_path, soft_paths_dict):
    hfbase.require_soft_paths(soft_paths_dict, ["gfa_editor_cli"])
    gfa_editor = command_for_python_script(soft_paths_dict["gfa_editor_cli"], soft_paths_dict)
    input_gfa = Path(gfa_path)
    state = checked_gfa_state(input_gfa)
    copied_input = Path(genome + "_checked_draft.input.gfa")
    merge_input = Path(genome + "_checked_draft.merged_raw.gfa")
    auto_merged = Path(genome + "_checked_draft.auto_merged.gfa")
    auto_resolved = Path(genome + "_checked_draft.auto_repeat_resolved.gfa")
    auto_summary = Path(genome + "_checked_draft.auto_repeat_summary.json")
    auto_history = Path(genome + "_checked_draft.auto_repeat_history.json")
    output_fasta = Path(genome + "_checked_draft.auto_merged.fasta")
    prepare_report = Path(genome + "_checked_draft.gfa_prepare.tsv")

    shutil.copyfile(input_gfa, copied_input)
    raw_gfa = parse_gfa(copied_input)
    candidate_gfa, candidate_mode = merge_unambiguous_gfa(raw_gfa)
    segment_delta = len(raw_gfa.segments) - len(candidate_gfa.segments)
    if state == "declared_merged":
        selected_gfa = raw_gfa
        merge_action = "input_declared_merged_kept"
        merge_mode = "input_declared_merged"
    elif segment_delta > 0:
        selected_gfa = candidate_gfa
        merge_action = "auto_checked_and_merged"
        merge_mode = candidate_mode
    else:
        selected_gfa = raw_gfa
        merge_action = "auto_checked_no_merge_needed"
        merge_mode = "input_already_merged_or_no_linear_compaction"
    write_gfa(selected_gfa, merge_input)

    run_command(gfa_editor + [
        "auto-merge",
        str(merge_input),
        str(auto_merged),
        "--reference-fasta",
        reference_fasta,
        "--resolved-output",
        str(auto_resolved),
        "--summary-json",
        str(auto_summary),
        "--history-json",
        str(auto_history),
    ])
    run_command(gfa_editor + [
        "export",
        str(auto_merged),
        str(output_fasta),
        "--format",
        "fasta",
    ])
    write_prepare_report(prepare_report, [
        ("input_gfa", str(input_gfa)),
        ("input_state", state),
        ("copied_input_gfa", str(copied_input)),
        ("merged_raw_gfa", str(merge_input)),
        ("merge_action", merge_action),
        ("merge_mode", merge_mode),
        ("candidate_merge_mode", candidate_mode),
        ("raw_segment_count", len(raw_gfa.segments)),
        ("merged_segment_count", len(selected_gfa.segments)),
        ("candidate_segment_delta", segment_delta),
        ("reference_fasta", reference_fasta),
        ("auto_merged_gfa", str(auto_merged)),
        ("auto_resolved_gfa", str(auto_resolved)),
        ("auto_summary_json", str(auto_summary)),
        ("auto_history_json", str(auto_history)),
        ("output_fasta", str(output_fasta)),
    ])
    return str(output_fasta)


# Usage: python
# Load the soft paths
soft_paths_file = sys.argv[1]
soft_paths_dict = hfbase.load_soft_paths(soft_paths_file)

# Parse other arguments
sample_index = sys.argv[2] # ATHiFi001
genome = sys.argv[3] # mito or plastid
genome_absolute_path = sys.argv[4] # for bait genome
before_fasta_absolute_path = sys.argv[5] # draft fasta or checked draft GFA
reads_absolute_path = sys.argv[6] # ATHiFi001.fastq
threads = sys.argv[7]

# check absolute paths
if not os.path.isabs(genome_absolute_path):
    print("The path to the genome reference fasta file is not an absolute path.")
    sys.exit(1)
if not os.path.isabs(before_fasta_absolute_path):
    print("The path to the draft fasta/GFA file is not an absolute path.")
    sys.exit(1)
if not os.path.isabs(reads_absolute_path):
    print("The path to the reads fastq file is not an absolute path.")
    sys.exit(1)

# create project folder by sample_index
if not os.path.exists(sample_index):
    os.makedirs(sample_index)
os.chdir(sample_index)
if not os.path.exists("draft_assembly"):
    os.makedirs("draft_assembly")
os.chdir("draft_assembly")
if not os.path.exists(genome):
    os.makedirs(genome)
os.chdir(genome)

input_kind = draft_input_kind(before_fasta_absolute_path)
if input_kind == "gfa":
    before_fasta_absolute_path = prepare_checked_gfa_for_polish(
        genome,
        genome_absolute_path,
        before_fasta_absolute_path,
        soft_paths_dict,
    )
elif input_kind != "fasta":
    print("Unsupported draft input type: " + before_fasta_absolute_path, file=sys.stderr)
    sys.exit(1)

hfref.replace_fasta_id(genome, before_fasta_absolute_path, "input_" + genome + ".fa")
rot_step = hfref.rotate_ref_to_non_repeat_region(genome, "input_" + genome + ".fa", soft_paths_dict, rotation=True)
hfbase.get_cli_output_lines("rm input_" + genome + ".fa")
hfref.flye_polish(genome, genome + "_rotated_" + str(rot_step) + ".fasta", genome + "_flye_polish", reads_absolute_path, soft_paths_dict, "HiFi", threads, correction=True)
# check the SVs by BLASTn
adj_count, rc_or_not, rot_or_not = hfref.aln_to_ref(genome, genome_absolute_path, genome + "_flye_polish.fasta", genome + "_flye_polish_aligned.fasta", "all_sorted_blastn_alignments.txt", soft_paths_dict)
hfrps.convert_blastn_alignments_to_table("all_sorted_blastn_alignments.txt", genome + "_flye_polish_aligned_blastn_alignments.xlsx")
if adj_count == 0:
    print("No change.")
if rc_or_not == 1:
    print("The draft genome has been reverse complemented.")
if rot_or_not == 1:
    print("The draft genome has been rotated.")

move_backup_info(genome)
os.chdir("../../..")
