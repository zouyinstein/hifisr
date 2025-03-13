import sys
import hifisr_functions.base as hfbase
import hifisr_functions.references as hfref
import hifisr_functions.reads as hfreads
import hifisr_functions.variants as hfvar
import hifisr_functions.annotations as hfanno
import hifisr_functions.reports as hfrps
from Bio import SeqIO
import os


# Usage: python get_variants_in_reads.py soft_paths.txt ATHiFi001 mito run_1 /mnt/software/bio/hifisr/dev/results/ref/mito_rotated_293434.fasta /mnt/software/bio/hifisr/dev/results/ATHiFi001/reads/sample_reads/sample_4000_mito.fastq 32
# Load the soft paths
soft_paths_file = sys.argv[1] 
soft_paths_dict = hfbase.load_soft_paths(soft_paths_file)

# Parse other arguments
sample_index = sys.argv[2] # ATHiFi001
genome = sys.argv[3] # mito
run_info = sys.argv[4] # run_1
genome_absolute_path = sys.argv[5]
reads_absolute_path = sys.argv[6]
threads = sys.argv[7]

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
if not os.path.exists(genome):
    os.makedirs(genome)
os.chdir(genome)
if not os.path.exists(run_info):
    os.makedirs(run_info)
os.chdir(run_info)
current_dir = os.getcwd()

# run multi-threading blastn, and get all_sorted_blastn_alignments.txt
hfvar.run_multi_threads_blastn(sample_index, genome, run_info, reads_absolute_path, genome_absolute_path, "all_sorted_blastn_alignments.txt", soft_paths_dict, threads)
hfvar.summarize_blastn_results("all_sorted_blastn_alignments.txt")

# plot type_2_rep
type_2_rep_file = current_dir + "/combined_excel/type_2_subtype_rep_NA_summary_anno.xlsx"
IDs_dir = current_dir + "/backup_info/IDs"
if os.path.exists(type_2_rep_file):
    hfrps.plot_bubble_type_2_rep_raw(type_2_rep_file, IDs_dir, genome_absolute_path)

# generate FL.fasta, partial.fasta, variant_cov.fasta
ref_len = hfvar.get_cov_reads(current_dir, IDs_dir, soft_paths_dict, genome_absolute_path) 
hfreads.cal_ID_coverage("FL", genome_absolute_path, "FL.fasta", "HiFi", soft_paths_dict, threads)
hfreads.cal_ID_coverage("partial", genome_absolute_path, "partial.fasta", "HiFi", soft_paths_dict, threads)
hfreads.cal_ID_coverage("variant", genome_absolute_path, "variant_cov.fasta", "HiFi", soft_paths_dict, threads)
hfrps.plot_coverage("FL_cov.txt", "partial_cov.txt", "variant_cov.txt", 1, ref_len, fig_length=12, fig_height=3)

# use variant_cov.fasta for variant calling, and get all_bcftools_calls.txt
hfvar.run_multi_threads_bcftools(sample_index, genome, run_info, "variant_cov.fasta", genome_absolute_path, "HiFi", "all_bcftools_calls.txt", soft_paths_dict, threads)
df = hfvar.snv_or_indel("all_bcftools_calls.txt", "snv_indel_reformat_single.xlsx")

# annotate the variants: types, effects, sequence features
df_anno_raw = hfanno.get_variant_types(genome_absolute_path, "snv_indel_reformat_single.xlsx", "variants_anno_single.xlsx")
hfanno.combine_variant_anno("variants_anno_single.xlsx", "variants_anno_combined.xlsx")

os.chdir("../../..")