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


modes = ["mecat_flye", "flye"]

# Usage: python get_draft_assembly_flye.py soft_paths_file sample_index genome bait_path sampled_reads_fastq threads [modes]
# Load the soft paths
soft_paths_file = sys.argv[1]
soft_paths_dict = hfbase.load_soft_paths(soft_paths_file)

# Parse other arguments
sample_index = sys.argv[2] # ATHiFi001
genome = sys.argv[3] # mito or plastid
genome_absolute_path = sys.argv[4] # for bait genome
reads_absolute_path = sys.argv[5] # sampled reads for mecat + flye and flye modes
threads = sys.argv[6]
if len(sys.argv) > 7:
    modes = [mode.strip() for mode in sys.argv[7].split(",") if mode.strip()]
modes = [{"mecat+flye": "mecat_flye", "mecat": "mecat_flye", "all": "flye"}.get(mode, mode) for mode in modes]
modes = list(dict.fromkeys(modes))
unsupported_modes = [mode for mode in modes if mode not in ["mecat_flye", "flye"]]
if unsupported_modes:
    print("Unsupported Flye draft mode(s): " + ",".join(unsupported_modes), file=sys.stderr)
    print("Use get_draft_assembly.py for simple_draft_asm modes.", file=sys.stderr)
    sys.exit(1)
if not modes:
    print("No Flye draft assembly modes were selected.", file=sys.stderr)
    sys.exit(1)

required_tools = ["bandage", "flye"]
if "mecat_flye" in modes:
    required_tools.append("mecat")
if "flye" in modes:
    required_tools.append("seqkit")
if reads_absolute_path.endswith(".gz"):
    required_tools.append("pigz")
hfbase.require_soft_paths(soft_paths_dict, list(dict.fromkeys(required_tools)))


# check absolute paths
if not os.path.isabs(genome_absolute_path):
    print("The path to the genome reference fasta file is not an absolute path.")
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

if genome == "plastid":
    genome_size = 150
elif genome == "mito":
    genome_size = 500
else:
    print("Unsupported genome: " + genome, file=sys.stderr)
    sys.exit(1)

gfa_outputs = []
reads_link = "reads.fastq"
if reads_absolute_path.endswith(".gz"):
    command_1 = "rm -f reads.fastq reads.fastq.gz && " + soft_paths_dict.get("pigz") + " -dc " + reads_absolute_path + " > reads.fastq"
else:
    command_1 = "rm -f reads.fastq reads.fastq.gz && ln -sf " + reads_absolute_path + " reads.fastq"
hfbase.run_checked(command_1)

if "mecat_flye" in modes:
    hfref.mecat_cns(genome, genome_size, reads_link, soft_paths_dict, threads)
    hfref.flye_assemble("mecat", genome, genome_size, "mecat_" + genome + "_" + str(genome_size) + ".fasta", soft_paths_dict, "HiFi", threads, correction=True)
    for suffix in ["before_rr", "after_rr"]:
        gfa_outputs.append(f"mecat_{genome}_{genome_size}K_{suffix}.gfa")
    hfbase.run_checked("rm -f mecat_" + genome + "_" + str(genome_size) + ".fasta")

if "flye" in modes:
    command_2 = soft_paths_dict.get("seqkit") + " fq2fa " + reads_link + " -o reads.fasta -j " + threads
    hfbase.run_checked(command_2)
    hfref.flye_assemble("all", genome, genome_size, "reads.fasta", soft_paths_dict, "HiFi", threads, correction=False)
    for suffix in ["before_rr", "after_rr"]:
        gfa_outputs.append(f"all_{genome}_{genome_size}K_{suffix}.gfa")
    hfbase.run_checked("rm -f reads.fasta")

hfbase.run_checked("rm -f reads.fastq reads.fastq.gz")
gfa_outputs = list(dict.fromkeys(gfa_outputs))
hfrps.get_gfa_blastn_png(genome_absolute_path, soft_paths_dict)
expected_outputs = gfa_outputs + [gfa_file.replace(".gfa", ".png") for gfa_file in gfa_outputs]
require_outputs(expected_outputs)
os.chdir("../../..")
