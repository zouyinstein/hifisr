# HiFiSR module guide:
# - base: command, file, and soft_paths helpers; import hifisr_functions.base as hfbase
# - reads: read extraction, filtering, sampling, and correction; import hifisr_functions.reads as hfreads
# - references: reference rotation, assembly, polishing, and alignment; import hifisr_functions.references as hfref
# - variants: read-variant calling, grouping, and frequency analysis; import hifisr_functions.variants as hfvar
# - transfer: organelle/nuclear transfer-fragment analysis; import hifisr_functions.transfer as hftrans
# - annotations: annotation tables and feature-level summaries; import hifisr_functions.annotations as hfanno
# - reports: read statistics, plots, Excel tables, and report outputs; import hifisr_functions.reports as hfrps

import hifisr_functions.base as hfbase
from Bio import SeqIO
import polars as pl
import sys
import os
import pysam

# Function purity marker. "pure" means deterministic from explicit inputs with
# no file, shell, environment, logging, or input-mutation side effects.
FUNCTION_PURITY = {
    "split_mtpt_reads": "impure",
    "split_reads_by_contig": "impure",
    "filt_length_qual": "impure",
    "random_sampling": "impure",
    "replace_reads_id": "impure",
    "cal_ID_coverage": "impure",
    "calc_error_rate": "impure",
    "correct_reads_by_canu": "impure",
}


def split_mtpt_reads(sample_index, sample_fastq_path, sample_platform, mito_fa, plastid_fa, soft_paths_dict, threads):
    platform_dict = {
            "HiFi": "map-hifi",
            "CLR": "map-pb",
            "ONT": "map-ont",
            "ultra-long": "map-ont",
        }
    command_1 = soft_paths_dict.get("seqkit") + " seq -ni " + mito_fa + " > mito_ids.txt"
    command_2 = soft_paths_dict.get("seqkit") + " seq -ni " + plastid_fa + " > plastid_ids.txt"
    command_3 = "cat " + mito_fa + " " + plastid_fa + " > mtpt.fa"
    command_4 = soft_paths_dict.get("minimap2") + " -t " + threads + " -ax " + platform_dict.get(sample_platform) + " mtpt.fa " + sample_fastq_path
    command_5 = soft_paths_dict.get("samtools") + " view -Sb -F 4 -@ " + threads + " - "
    command_6 = soft_paths_dict.get("samtools") + " sort -@ " + threads + " - -o reads.sorted.bam"
    command_7 = soft_paths_dict.get("samtools") + " index reads.sorted.bam"
    command_8 = soft_paths_dict.get("bamtools") + " split -in reads.sorted.bam -reference"
    command_9 = "rm -rf " + sample_index + "_mito.fastq " + sample_index + "_plastid.fastq " + sample_index + "_mito.fastq.gz " + sample_index + "_plastid.fastq.gz"
    command_10 = "cat mito_ids.txt | while read ID; do " + soft_paths_dict.get("samtools") + " fastq reads.sorted.REF_${ID}.bam -@ " + threads + " >> " + sample_index + "_mito.fastq; done"
    command_11 = "cat plastid_ids.txt | while read ID; do " + soft_paths_dict.get("samtools") + " fastq reads.sorted.REF_${ID}.bam -@ " + threads + " >> " + sample_index + "_plastid.fastq; done"
    command_12 = "rm -rf mito_ids.txt plastid_ids.txt mtpt.fa reads.sorted.bam reads.sorted.bam.bai reads.sorted.*.bam" 
    commands = command_1 + " ; " + command_2 + " ; " + command_3 + " ; " + command_4 + " | " + command_5 + " | " + command_6 + " && " + command_7 + " && " + command_8 + " ; " + command_9 + " ; " + command_10 + " ; " + command_11 + " ; " + command_12
    hfbase.get_cli_output_lines(commands, side_effect = True)
    return


def split_reads_by_contig(sample_fastq_path, sample_platform, contigs_fa, soft_paths_dict, threads):
    platform_dict = {
            "HiFi": "map-hifi",
            "CLR": "map-pb",
            "ONT": "map-ont",
            "ultra-long": "map-ont",
        }
    command_1 = soft_paths_dict.get("seqkit") + " seq -ni " + contigs_fa + " > contigs_ids.txt"
    command_2 = soft_paths_dict.get("minimap2") + " -t " + threads + " -ax " + platform_dict.get(sample_platform) + " " + contigs_fa + " " + sample_fastq_path
    command_3 = soft_paths_dict.get("samtools") + " view -Sb -F 4 -@ " + threads + " - "
    command_4 = soft_paths_dict.get("samtools") + " sort -@ " + threads + " - -o reads.sorted.bam"
    command_5 = soft_paths_dict.get("samtools") + " index reads.sorted.bam"
    command_6 = soft_paths_dict.get("bamtools") + " split -in reads.sorted.bam -reference"
    command_7 = "cat contigs_ids.txt | while read ID; do " + soft_paths_dict.get("samtools") + " fastq reads.sorted.REF_${ID}.bam -@ " + threads + " > ${ID}.fastq; done"
    command_8 = "rm -rf contigs_ids.txt reads.sorted.bam reads.sorted.bam.bai reads.sorted.*.bam"
    command_9 = soft_paths_dict.get("seqkit") + " split -i " + contigs_fa
    commands = command_1 + " ; " + command_2 + " | " + command_3 + " | " + command_4 + " && " + command_5 + " && " + command_6 + " ; " + command_7 + " ; " + command_8 + " ; " + command_9
    hfbase.get_cli_output_lines(commands, side_effect = True)
    results =  hfbase.get_cli_output_lines("seqkit stat -T " + contigs_fa + " | tail -n1 | cut -f4", side_effect = False)
    count = results[0].split()[0]
    return count


def filt_length_qual(prefix, id_length_qual_file, filt_length=0, filt_qual=0):
    df = pl.read_csv(id_length_qual_file, separator="\t", has_header=False, new_columns=["id", "length", "qual"])
    df = df.filter(pl.col("id").is_not_null())
    df = df.with_columns([
        pl.col("length").cast(pl.Int64),
        pl.col("qual").cast(pl.Float64)
    ])
    df_filtered = df.filter((pl.col("length") >= filt_length) & (pl.col("qual") >= filt_qual))
    filt_read_number = str(df_filtered.shape[0])
    filt_bases = str(df_filtered["length"].sum())
    id_length_qual_file_filt = "filt_" + prefix + "_id_length_qual.txt"
    df_filtered.write_csv(id_length_qual_file_filt, separator="\t", include_header=False)
    return id_length_qual_file_filt, filt_read_number, filt_bases


def random_sampling(prefix, id_length_qual_file, sample_number):
    df = pl.read_csv(id_length_qual_file, separator="\t", has_header=False, new_columns=["id", "length", "qual"])
    df = df.filter(pl.col("id").is_not_null())
    df = df.with_columns([
        pl.col("length").cast(pl.Int64),
        pl.col("qual").cast(pl.Float64)
    ])
    total_rows = int(df.shape[0])
    if sample_number > total_rows:
        sample_number = total_rows
    df_sampled = df.sample(n=sample_number, with_replacement=False, seed=2025)
    sample_read_number = str(df_sampled.shape[0])
    sample_bases = str(df_sampled["length"].sum())
    id_length_qual_file_sampled = prefix + "_id_length_qual.txt"
    df_sampled.write_csv(id_length_qual_file_sampled, separator="\t", include_header=False)
    return id_length_qual_file_sampled, sample_read_number, sample_bases


def replace_reads_id(reads_file, new_reads_file):
    with open(reads_file, "rt") as fin, open(new_reads_file, "wt") as fout:
        for record in SeqIO.parse(fin, "fasta"):
            record.id = record.id.replace("/", "_")
            record.description = ""
            SeqIO.write(record, fout, "fasta")
    return new_reads_file


def cal_ID_coverage(prefix, ref_fasta, reads_fasta, sample_platform, soft_paths_dict, threads):
    platform_dict = {
            "HiFi": "map-hifi",
            "CLR": "map-pb",
            "ONT": "map-ont",
            "ultra-long": "map-ont",
        }
    command_1 = soft_paths_dict.get("minimap2") + " -t " + threads + " -ax " + platform_dict.get(sample_platform) + " " + ref_fasta + " " + reads_fasta
    command_2 = soft_paths_dict.get("samtools") + " view -Sb -F 0x100 -@ " + threads + " -"
    command_3 = soft_paths_dict.get("samtools") + " sort -@ " + threads + " -o " + prefix + ".sorted.bam "
    command_4 = soft_paths_dict.get("samtools") + " index -@ " + threads + " " + prefix + ".sorted.bam "
    command_5 = soft_paths_dict.get("samtools") + " depth -a -J -@ " + threads + " " + prefix + ".sorted.bam | cut -f2- > " + prefix + "_cov.txt"
    command_6 = "rm " + prefix + ".sorted.bam " + prefix + ".sorted.bam.bai"
    commands = command_1 + " | " + command_2 + " | " + command_3 + " ; " + command_4 + " ; " + command_5 + " ; " + command_6
    hfbase.get_cli_output_lines(commands, side_effect = True)
    return


def calc_error_rate(bam_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    total_bases = 0
    total_errors = 0
    for aln in bam.fetch():
        if aln.is_unmapped:
            continue
        if aln.has_tag("NM"):
            total_errors += aln.get_tag("NM")
            total_bases += aln.query_alignment_length
    bam.close()
    return total_errors / total_bases if total_bases > 0 else 0


def correct_reads_by_canu(sample_index, genome, bait_ref, sample_num, sample_platform, soft_paths_dict, threads, result_dir="canu_corrected_reads"):
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    canu_platform_dict = {
            "HiFi": "-pacbio-hifi",
            "CLR": "-pacbio-raw",
            "ONT": "-nanopore-raw",
            "ultra-long": "-nanopore-raw",
    }
    minimap2_platform_dict = {
            "HiFi": "map-hifi",
            "CLR": "map-pb",
            "ONT": "map-ont",
            "ultra-long": "map-ont",
    }
    canu_memory = "128g"
    canu_outcoverage_dict = {
            "mito": "100",
            "plastid": "200",
    }
    canu_genome_size_dict = {
            "mito": "500k",
            "plastid": "150k",
    }
    command_1 = soft_paths_dict.get("canu") + " -correct -p " + sample_index + "_" + genome + "  -d " + genome + "_canu genomeSize=" + canu_genome_size_dict.get(genome) + " " + canu_platform_dict.get(sample_platform) + " sample_reads/sample_" + str(sample_num) + "_" + genome + ".fastq maxThreads=" + threads + " maxMemory=" + canu_memory + " useGrid=0 corOutCoverage=" + canu_outcoverage_dict.get(genome)
    command_2 = "mv " + genome + "_canu/" + sample_index + "_" + genome + ".correctedReads.fasta.gz " + result_dir + "/" + sample_index + "_" + genome + ".correctedReads.fasta.gz"
    command_3 = "rm -r " + genome + "_canu"
    hfbase.get_cli_output_lines(command_1 + " && " + command_2 + " && " + command_3, side_effect = True)
    # evaluate error rate against the given reference
    command_4 = soft_paths_dict.get("minimap2") + " -t " + threads + " -ax " + minimap2_platform_dict.get(sample_platform) + " " + bait_ref + " sample_reads/sample_" + str(sample_num) + "_" + genome + ".fastq | " + soft_paths_dict.get("samtools") + " view -bS - | " + soft_paths_dict.get("samtools") + " sort -o " + genome + "_raw.bam"
    command_5 = soft_paths_dict.get("samtools") + " index " + genome + "_raw.bam"
    command_6 = soft_paths_dict.get("minimap2") + " -t " + threads + " -ax " + minimap2_platform_dict.get(sample_platform) + " " + bait_ref + " " + result_dir + "/" + sample_index + "_" + genome + ".correctedReads.fasta.gz | " + soft_paths_dict.get("samtools") + " view -bS - | " + soft_paths_dict.get("samtools") + " sort -o " + genome + "_corrected.bam"
    command_7 = soft_paths_dict.get("samtools") + " index " + genome + "_corrected.bam"
    hfbase.get_cli_output_lines(command_4 + " && " + command_5 + " && " + command_6 + " && " + command_7, side_effect = True)
    raw_err = calc_error_rate(genome + "_raw.bam")
    corr_err = calc_error_rate(genome + "_corrected.bam")
    with open(os.path.join(result_dir, sample_index + "_" + genome + "_canu_correction_evaluation.txt"), "wt") as fout:
        fout.write("Canu Correction Evaluation\n")
        fout.write(f"Raw read error rate:       {raw_err:.6%}\n")
        fout.write(f"Corrected read error rate: {corr_err:.6%}\n")
        fout.write(f"Error reduction fold:      {raw_err/corr_err:.2f}x\n")
    command_8 = "rm " + genome + "_raw.bam " + genome + "_raw.bam.bai " + genome + "_corrected.bam " + genome + "_corrected.bam.bai"
    hfbase.get_cli_output_lines(command_8, side_effect = True)
    return
