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
import hifisr_functions.reads as hfreads
import hifisr_functions.reports as hfrps
import os
import shutil
import shlex


DEFAULT_ANALYSIS_READ_LIMIT = 50000


def parse_read_limit(value, default=DEFAULT_ANALYSIS_READ_LIMIT):
    if value is None or str(value).strip() == "":
        return default
    value = str(value).strip().lower()
    if value in {"none", "all", "no", "false", "0"}:
        return 0
    read_limit = int(value)
    if read_limit < 0:
        raise ValueError("Read limit must be >= 0: " + value)
    return read_limit


def parse_count(value):
    return int(str(value).replace(",", ""))


def move_if_exists(filename, dest_dir):
    if os.path.exists(filename):
        dest_path = os.path.join(dest_dir, filename)
        if os.path.exists(dest_path):
            os.remove(dest_path)
        shutil.move(filename, dest_path)


def prepare_analysis_reads(genome, read_limit, soft_paths_dict, threads, summary_rows):
    backup_info_dir = "backup_info"
    os.makedirs(backup_info_dir, exist_ok=True)

    raw_fastq = genome + "_all.fastq"
    final_fastq = genome + ".fastq"
    final_fastq_gz = final_fastq + ".gz"
    full_stats = genome + "_id_length_qual.txt"
    analysis_prefix = genome + "_analysis"
    analysis_ids = os.path.join(backup_info_dir, genome + "_analysis_ids.txt")

    for filename in [final_fastq, final_fastq_gz]:
        if os.path.exists(filename):
            os.remove(filename)
    for filename in [
        analysis_prefix + "_id_length_qual.txt",
        analysis_prefix + "_length_qual_distribution.pdf",
        analysis_prefix + "_length_qual_2d_distribution.pdf",
    ]:
        if os.path.exists(filename):
            os.remove(filename)
        backup_filename = os.path.join(backup_info_dir, filename)
        if os.path.exists(backup_filename):
            os.remove(backup_filename)
    for filename in [analysis_ids]:
        if os.path.exists(filename):
            os.remove(filename)

    full_stats, raw_read_number, raw_bases = hfrps.get_fastq_stats(genome, raw_fastq, soft_paths_dict, threads)
    raw_read_count = parse_count(raw_read_number)
    if raw_read_count == 0:
        open(final_fastq, "wt").close()
        if os.path.exists(raw_fastq):
            os.remove(raw_fastq)
        summary_rows.append([
            genome,
            "0",
            str(read_limit),
            "no",
            "0",
            "0",
        ])
        return

    hfrps.plot_length_qual(genome, "HiFi", full_stats, raw_read_number, raw_bases)
    should_sample = read_limit > 0 and raw_read_count > read_limit

    if should_sample:
        sampled_stats, final_read_number, final_bases = hfreads.random_sampling(analysis_prefix, full_stats, sample_number=read_limit)
        hfrps.plot_length_qual(analysis_prefix, "HiFi", sampled_stats, final_read_number, final_bases)
        command_1 = "cut -f 1 " + shlex.quote(sampled_stats) + " > " + shlex.quote(analysis_ids)
        command_2 = (
            soft_paths_dict.get("seqkit") + " grep -f " + shlex.quote(analysis_ids)
            + " " + shlex.quote(raw_fastq) + " > " + shlex.quote(final_fastq)
        )
        hfbase.get_cli_output_lines(command_1 + " && " + command_2, side_effect=True)
        os.remove(raw_fastq)
    else:
        os.replace(raw_fastq, final_fastq)
        final_read_number = raw_read_number
        final_bases = raw_bases

    for filename in [
        analysis_prefix + "_id_length_qual.txt",
        analysis_prefix + "_length_qual_distribution.pdf",
        analysis_prefix + "_length_qual_2d_distribution.pdf",
    ]:
        move_if_exists(filename, backup_info_dir)

    summary_rows.append([
        genome,
        str(raw_read_count),
        str(read_limit),
        "yes" if should_sample else "no",
        str(parse_count(final_read_number)),
        str(final_bases),
    ])


# Usage: python get_mtpt_reads.py soft_paths_file sample_index mito_absolute_path plastid_absolute_path reads_absolute_path threads [mito_analysis_read_limit] [plastid_analysis_read_limit]
# Load the soft paths
soft_paths_file = sys.argv[1]
soft_paths_dict = hfbase.load_soft_paths(soft_paths_file)

# Parse other arguments
sample_index = sys.argv[2] # ATHiFi001
mito_absolute_path = sys.argv[3]
plastid_absolute_path = sys.argv[4]
reads_absolute_path = sys.argv[5] # ATHiFi001.fastq
threads = sys.argv[6]
mito_analysis_read_limit = parse_read_limit(sys.argv[7] if len(sys.argv) > 7 else None)
plastid_analysis_read_limit = parse_read_limit(sys.argv[8] if len(sys.argv) > 8 else None)

# check absolute paths
if not os.path.isabs(mito_absolute_path):
    print("The path to the mito reference fasta file is not an absolute path.")
    sys.exit(1)
if not os.path.isabs(plastid_absolute_path):
    print("The path to the plastid reference fasta file is not an absolute path.")
    sys.exit(1)
if not os.path.isabs(reads_absolute_path):
    print("The path to the reads fastq file is not an absolute path.")
    sys.exit(1)

# create project folder by sample_index
if not os.path.exists(sample_index):
    os.makedirs(sample_index)
os.chdir(sample_index)
if not os.path.exists("reads"):
    os.makedirs("reads")
os.chdir("reads")

# split the reads into mito and plastid reads
hfref.replace_fasta_id("mito", mito_absolute_path, "mito.fa")
hfref.replace_fasta_id("plastid", plastid_absolute_path, "plastid.fa")
split_prefix = "split_mtpt"
if reads_absolute_path.endswith(".gz"):
    raw_reads_link = "all.fastq.gz"
else:
    raw_reads_link = "all.fastq"
cleanup_files = [
    "all.fastq", "all.fastq.gz",
    "mito.fastq", "mito.fastq.gz", "mito_all.fastq",
    "plastid.fastq", "plastid.fastq.gz", "plastid_all.fastq",
    split_prefix + "_mito.fastq", split_prefix + "_plastid.fastq",
    sample_index + ".fastq", sample_index + ".fastq.gz",
    sample_index + "_mito.fastq", sample_index + "_mito.fastq.gz",
    sample_index + "_plastid.fastq", sample_index + "_plastid.fastq.gz",
]
hfbase.run_checked("rm -f " + " ".join(cleanup_files))
command_1 = "ln -sf " + reads_absolute_path + " " + raw_reads_link
hfbase.get_cli_output_lines(command_1, side_effect = True)
hfreads.split_mtpt_reads(split_prefix, raw_reads_link, "HiFi", "mito.fa", "plastid.fa", soft_paths_dict, threads)
os.replace(split_prefix + "_mito.fastq", "mito_all.fastq")
os.replace(split_prefix + "_plastid.fastq", "plastid_all.fastq")

id_length_qual_file, total_read_number, total_bases = hfrps.get_fastq_stats("all", raw_reads_link, soft_paths_dict, threads)
hfrps.plot_length_qual("all", "HiFi", id_length_qual_file, total_read_number, total_bases)

sampling_summary_rows = [["genome", "extracted_read_count", "analysis_read_limit", "sampled", "analysis_read_count", "analysis_bases"]]
prepare_analysis_reads("mito", mito_analysis_read_limit, soft_paths_dict, threads, sampling_summary_rows)
prepare_analysis_reads("plastid", plastid_analysis_read_limit, soft_paths_dict, threads, sampling_summary_rows)
with open("backup_info/downstream_read_sampling_summary.tsv", "wt") as fout:
    for row in sampling_summary_rows:
        print("\t".join(row), file=fout)

pigz = soft_paths_dict.get("pigz", "pigz")
command_1 = pigz + " -p " + threads + " -f mito.fastq plastid.fastq"
hfbase.get_cli_output_lines(command_1, side_effect = True)

os.chdir("../..")
