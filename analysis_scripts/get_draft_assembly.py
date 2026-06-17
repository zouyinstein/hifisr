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
import os


def require_outputs(paths):
    missing = [path for path in paths if not os.path.exists(path)]
    empty = [path for path in paths if os.path.exists(path) and os.path.getsize(path) == 0]
    if missing or empty:
        if missing:
            print("Missing expected draft outputs:", file=sys.stderr)
            for path in missing:
                print("  " + path, file=sys.stderr)
        if empty:
            print("Empty expected draft outputs:", file=sys.stderr)
            for path in empty:
                print("  " + path, file=sys.stderr)
        sys.exit(1)


modes = ["ms", "mh", "mx"]

# Usage: python get_draft_assembly.py soft_paths_file sample_index genome bait_path full_reads_fastq threads [modes] [full_reads_fastq]
# Load the soft paths
soft_paths_file = sys.argv[1]
soft_paths_dict = hfbase.load_soft_paths(soft_paths_file)

# Parse other arguments
sample_index = sys.argv[2] # ATHiFi001
genome = sys.argv[3] # mito or plastid
genome_absolute_path = sys.argv[4] # for bait genome
reads_absolute_path = sys.argv[5] # full mito/plastid reads for simple_draft_asm
threads = sys.argv[6]
if len(sys.argv) > 7:
    modes = [mode.strip() for mode in sys.argv[7].split(",") if mode.strip()]
if genome == "plastid":
    modes = [{"ms": "ps", "mh": "ph"}.get(mode, mode) for mode in modes if mode != "mx"]
modes = list(dict.fromkeys(modes))
full_reads_absolute_path = sys.argv[8] if len(sys.argv) > 8 else reads_absolute_path
simple_draft_asm_modes = [mode for mode in modes if mode in ["ms", "mh", "mx", "pl", "ps", "ph"]]
unsupported_modes = [mode for mode in modes if mode not in ["ms", "mh", "mx", "pl", "ps", "ph"]]
if unsupported_modes:
    print("Unsupported simple_draft_asm mode(s): " + ",".join(unsupported_modes), file=sys.stderr)
    print("Use get_draft_assembly_flye.py for mecat_flye/flye modes.", file=sys.stderr)
    sys.exit(1)
if not modes:
    print("No draft assembly modes were selected.", file=sys.stderr)
    sys.exit(1)

required_tools = ["simple_draft_asm"]
hfbase.require_soft_paths(soft_paths_dict, list(dict.fromkeys(required_tools)))


# check absolute paths
if not os.path.isabs(genome_absolute_path):
    print("The path to the genome reference fasta file is not an absolute path.")
    sys.exit(1)
if not os.path.isabs(full_reads_absolute_path):
    print("The path to the full genome reads fastq file is not an absolute path.")
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

if genome == "plastid":
    genome_size = 150
elif genome == "mito":
    genome_size = 500
else:
    print("Unsupported genome: " + genome, file=sys.stderr)
    sys.exit(1)

gfa_outputs = []

full_reads_link = genome + ".fastq.gz" if full_reads_absolute_path.endswith(".gz") else genome + ".fastq"
command_1 = "rm -f " + genome + ".fastq " + genome + ".fastq.gz && ln -sf " + full_reads_absolute_path + " " + full_reads_link
hfbase.run_checked(command_1)

for mode in simple_draft_asm_modes:
    gfa_outputs.append(hfref.simple_draft_asm(genome, genome_size, mode, full_reads_link, soft_paths_dict, threads))

hfbase.run_checked("rm -f " + full_reads_link)

gfa_outputs = list(dict.fromkeys(gfa_outputs))
image_outputs = hfrps.get_gfa_reference_pdf(genome_absolute_path, soft_paths_dict)
expected_outputs = gfa_outputs + image_outputs
require_outputs(expected_outputs)
os.chdir("../../..")
