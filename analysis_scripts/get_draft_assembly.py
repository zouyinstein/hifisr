import sys
import hifisr_functions.base as hfbase
import hifisr_functions.references as hfref
import hifisr_functions.reports as hfrps
import os


# Usage: python get_draft_assembly.py soft_paths_file genome bait_path reads_fastq threads
# Load the soft paths
soft_paths_file = sys.argv[1]
soft_paths_dict = hfbase.load_soft_paths(soft_paths_file)

# Parse other arguments
sample_index = sys.argv[2] # ATHiFi001
genome = sys.argv[3] # mito or plastid
genome_absolute_path = sys.argv[4] # for bait genome
reads_absolute_path = sys.argv[5] # ATHiFi001.fastq
threads = sys.argv[6]


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

# add sample_index
command_1 = "ln -sf " + reads_absolute_path + " reads.fastq"
ret = hfbase.get_cli_output_lines(command_1, side_effect = True)
if genome == "plastid":
    genome_size = 150
elif genome == "mito":
    genome_size = 500

hfref.mecat_cns(genome, genome_size, "reads.fastq", soft_paths_dict, threads)
hfref.flye_assemble("mecat", genome, genome_size, "mecat_" + genome + "_" + str(genome_size) + ".fasta", "HiFi", threads, correction=True)
# convert reads.fastq to fasta using seqkit
command_2 = soft_paths_dict.get("seqkit") + " fq2fa reads.fastq -o reads.fasta -j " + threads
hfbase.get_cli_output_lines(command_2, side_effect = True)
hfref.flye_assemble("all", genome, genome_size, "reads.fasta", "HiFi", threads, correction=False)
hfrps.get_gfa_blastn_png(genome_absolute_path, soft_paths_dict)
os.chdir("../../..")
