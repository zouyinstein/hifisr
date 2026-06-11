# hifisr_functions Reference

This page gives a short usage map for the helper modules in
`hifisr_functions/` and records whether each public helper is a pure function.

Purity definition used here:

- `pure`: deterministic from explicit arguments and no file, shell,
  environment, logging, or input-mutation side effects.
- `impure`: reads or writes files, runs shell commands, reads environment or
  filesystem state, prints/logs, uses sampling, mutates an input object, or
  depends on stateful class data.

Only two helpers are currently pure: `references.get_subseq()` and
`variants.get_next_groups()`. Most workflow helpers are intentionally impure
because they wrap command-line tools, parse files, or write workflow outputs.

## Import Aliases

```python
import hifisr_functions.base as hfbase
import hifisr_functions.reads as hfreads
import hifisr_functions.references as hfref
import hifisr_functions.variants as hfvar
import hifisr_functions.transfer as hftrans
import hifisr_functions.annotations as hfanno
import hifisr_functions.reports as hfrps
```

## base

Command execution, file reading, `soft_paths` parsing, and runtime dependency
checks.

| Function | Purity | Typical use | Notes |
| --- | --- | --- | --- |
| `run_checked(commands)` | impure | `hfbase.run_checked("samtools faidx ref.fa")` | Runs a shell command and raises `RuntimeError` on non-zero exit. |
| `get_cli_output_lines(commands, side_effect=False)` | impure | `lines = hfbase.get_cli_output_lines("seqkit stats reads.fq")` | Runs a shell command and returns stdout lines; `side_effect=True` returns the exit status. |
| `get_file_lines(file)` | impure | `lines = hfbase.get_file_lines("ids.txt")` | Reads a text file into stripped lines. |
| `load_soft_paths(soft_paths_file, validate=False, required_tools=None, print_paths=True)` | impure | `soft_paths = hfbase.load_soft_paths("deps/soft_paths.txt", validate=True)` | Parses tab-delimited tool paths and can validate executables. |
| `validate_soft_paths(soft_paths_dict, required_tools=None, soft_paths_file="soft_paths.txt", check_all=True)` | impure | `hfbase.validate_soft_paths(soft_paths, ["python", "minimap2"])` | Checks required tool entries and executable availability. |
| `require_soft_paths(soft_paths_dict, required_tools)` | impure | `hfbase.require_soft_paths(soft_paths, ["samtools"])` | Convenience validator for scripts that need a subset of tools. |
| `parse_python_requirements(requirements_file)` | impure | `reqs = hfbase.parse_python_requirements("requirements-dev.txt")` | Parses package names from a requirements file. |
| `check_python_packages(requirements_file, python_executable=None)` | impure | `report = hfbase.check_python_packages("requirements-dev.txt")` | Verifies Python packages, `sqlite3`, and Snakemake in the selected Python. |

## reads

Read extraction, read filtering, random sampling, ID normalization, read-depth
calculation, and Canu correction.

| Function | Purity | Typical use | Notes |
| --- | --- | --- | --- |
| `split_mtpt_reads(sample_index, sample_fastq_path, sample_platform, mito_fa, plastid_fa, soft_paths_dict, threads)` | impure | `hfreads.split_mtpt_reads("W3-5-2", reads, "hifi", mito, plastid, soft_paths, 8)` | Aligns reads to mito/plastid references and writes organelle-specific FASTQ files. |
| `split_reads_by_contig(sample_fastq_path, sample_platform, contigs_fa, soft_paths_dict, threads)` | impure | `n = hfreads.split_reads_by_contig(reads, "hifi", contigs, soft_paths, 8)` | Splits reads by reference contig after alignment. |
| `filt_length_qual(prefix, id_length_qual_file, filt_length=0, filt_qual=0)` | impure | `hfreads.filt_length_qual("sample", "id_length_qual.txt", 1000, 20)` | Filters a read length/quality table and writes `filt_*` output. |
| `random_sampling(prefix, id_length_qual_file, sample_number)` | impure | `hfreads.random_sampling("sample", "id_length_qual.txt", 4000)` | Samples read IDs with a fixed random seed and writes sample tables. |
| `replace_reads_id(reads_file, new_reads_file)` | impure | `hfreads.replace_reads_id("in.fa", "out.fa")` | Rewrites FASTA read IDs so `/` becomes `_`. |
| `cal_ID_coverage(prefix, ref_fasta, reads_fasta, sample_platform, soft_paths_dict, threads)` | impure | `hfreads.cal_ID_coverage("run", ref, reads, "hifi", soft_paths, 8)` | Runs alignment/depth commands and writes coverage output. |
| `calc_error_rate(bam_file)` | impure | `err = hfreads.calc_error_rate("reads.bam")` | Reads a BAM file and computes a mismatch/error-rate summary. |
| `correct_reads_by_canu(sample_index, genome, bait_ref, sample_num, sample_platform, soft_paths_dict, threads, result_dir="canu_corrected_reads")` | impure | `hfreads.correct_reads_by_canu("W3-5-2", "mito", ref, 4000, "hifi", soft_paths, 8)` | Runs Canu correction and writes corrected-read evaluation outputs. |

## references

Reference sequence editing, reverse-complement generation, reference rotation,
assembly, polishing, and final alignment helpers.

| Function | Purity | Typical use | Notes |
| --- | --- | --- | --- |
| `replace_fasta_id(genome, input_fasta_path, output_fasta_path)` | impure | `hfref.replace_fasta_id("mito", "in.fa", "out.fa")` | Reads FASTA and writes records with normalized IDs. |
| `get_subseq(ref, start, end, flank=0)` | pure | `seq = hfref.get_subseq(seq, 100, 200, flank=20)` | Returns a bounded subsequence from an in-memory sequence. |
| `get_rc(input_fa_path, rc_fa_path)` | impure | `hfref.get_rc("in.fa", "rc.fa")` | Reads FASTA and writes reverse-complement FASTA. |
| `rotate_ref_to_non_repeat_region(genome, genome_fasta_path, soft_paths_dict, rotation=False)` | impure | `hfref.rotate_ref_to_non_repeat_region("mito", ref, soft_paths, True)` | Uses sequence/repeat information to choose or create a rotated reference. |
| `find_continous_zeros(info_list, repeat_pos_array)` | impure | `pos = hfref.find_continous_zeros(info_list, repeat_array)` | Mutates the supplied repeat-position array while finding a non-repeat region. |
| `rotate_fasta(genome_fasta_path, rotated_fasta_path, step)` | impure | `hfref.rotate_fasta("in.fa", "rotated.fa", 1200)` | Writes a rotated FASTA sequence. |
| `mecat_cns(genome, genome_size, reads, soft_paths_dict, threads)` | impure | `hfref.mecat_cns("mito", "500k", reads, soft_paths, 8)` | Runs MECAT consensus workflow outputs. |
| `flye_assemble(prefix, genome, genome_size, reads, soft_paths_dict, sample_platform, threads, correction=True)` | impure | `hfref.flye_assemble("sample", "mito", "500k", reads, soft_paths, "hifi", 8)` | Runs Flye assembly and writes assembly files. |
| `flye_polish(genome, before_fasta_path, after_fasta_prefix, reads, soft_paths_dict, sample_platform, threads, correction=True)` | impure | `hfref.flye_polish("mito", "draft.fa", "polished", reads, soft_paths, "hifi", 8)` | Runs polishing steps and writes polished FASTA files. |
| `aln_to_ref(genome, genome_absolute_path, polished_fasta_absolute_path, final_fasta_name, results_filename, soft_paths_dict)` | impure | `hfref.aln_to_ref("mito", ref, polished, "final.fa", "results", soft_paths)` | Aligns polished sequence to the reference and writes final FASTA/alignment products. |
| `get_flipped_fasta(genome_fasta_path, flipped_fasta_path, flip_start, flip_end)` | impure | `hfref.get_flipped_fasta("in.fa", "flipped.fa", 100, 500)` | Writes a FASTA with one interval reverse-complemented. |

## reports

FASTQ statistics, plotting helpers, graph visualization, and BLAST table
conversion.

| Function | Purity | Typical use | Notes |
| --- | --- | --- | --- |
| `get_fastq_stats(prefix, sample_fastq_path, soft_paths_dict, threads)` | impure | `hfrps.get_fastq_stats("sample", reads, soft_paths, 8)` | Runs read-stat commands and writes summary files. |
| `plot_length_qual(prefix, sample_platform, id_length_qual_file, total_read_number, total_bases)` | impure | `hfrps.plot_length_qual("sample", "hifi", "id_length_qual.txt", n, bases)` | Creates read length/quality plots. |
| `plot_bubble_type_2_rep_raw(table_file, IDs_dir, ref_fasta)` | impure | `hfrps.plot_bubble_type_2_rep_raw("table.xlsx", "IDs", ref)` | Reads variant/report inputs and writes a plot. |
| `plot_coverage(cov_file_1, cov_file_2, cov_file_3, start, end, fig_length=12, fig_height=3)` | impure | `hfrps.plot_coverage(a, b, c, 1, 5000)` | Reads coverage files and produces a coverage plot. |
| `get_gfa_blastn_png(genome_absolute_path, soft_paths_dict)` | impure | `hfrps.get_gfa_blastn_png(ref, soft_paths)` | Runs graph/BLAST visualization tools and writes PNG outputs. |
| `convert_blastn_alignments_to_table(blastn_alignments_file, output_file)` | impure | `hfrps.convert_blastn_alignments_to_table("blastn.txt", "blastn.xlsx")` | Converts BLAST alignment text into a tabular output. |

## annotations

Variant annotation, depth/frequency annotation, and organelle/nuclear transfer
fragment summaries.

| Function | Purity | Typical use | Notes |
| --- | --- | --- | --- |
| `get_variant_types(ref_fasta_file, input_excel, output_excel)` | impure | `hfanno.get_variant_types(ref, "variants.xlsx", "typed.xlsx")` | Reads reference/table inputs and writes variant-type annotations. |
| `combine_variant_anno(input_excel, output_excel)` | impure | `hfanno.combine_variant_anno("typed.xlsx", "combined.xlsx")` | Combines variant annotation sheets. |
| `add_depth_and_frq(input_excel, variant_cov_file, output_excel, filt_excel, engine="openpxl")` | impure | `hfanno.add_depth_and_frq("variants.xlsx", "cov.txt", "anno.xlsx", "filtered.xlsx")` | Adds depth/frequency columns and writes full/filtered Excel files. |
| `run_transfer_blastn(chr_absolute_path, mito_absolute_path, plastid_absolute_path, soft_paths_dict)` | impure | `hfanno.run_transfer_blastn(chr_fa, mito_fa, plastid_fa, soft_paths)` | Runs BLAST for NUMT/NUPT-style transfer-fragment analysis. |
| `merge_fragments_files(chr_name, num_mito, num_plastid)` | impure | `hfanno.merge_fragments_files("Chr1", mito_n, plastid_n)` | Reads per-fragment files and writes merged transfer-fragment tables. |
| `remove_inner_fragments(all_fragments_df)` | impure | `filtered = hfanno.remove_inner_fragments(df)` | Mutates and filters an input DataFrame-like object. |
| `merge_numt_nupt(num_Chr, num_mito, num_plastid)` | impure | `hfanno.merge_numt_nupt(chr_n, mito_n, plastid_n)` | Merges NUMT/NUPT fragment summaries across references. |
| `get_raw_transfer_groups(num_Chr)` | impure | `hfanno.get_raw_transfer_groups(chr_n)` | Reads raw transfer-fragment groups and writes grouped outputs. |

## transfer

Stable `hftrans` import alias for transfer-fragment helpers. These functions are
currently re-exported from `annotations` and have the same purity and behavior.

| Function | Purity | Typical use | Notes |
| --- | --- | --- | --- |
| `run_transfer_blastn(chr_absolute_path, mito_absolute_path, plastid_absolute_path, soft_paths_dict)` | impure | `hftrans.run_transfer_blastn(chr_fa, mito_fa, plastid_fa, soft_paths)` | Delegates to `annotations.run_transfer_blastn()`. |
| `merge_fragments_files(chr_name, num_mito, num_plastid)` | impure | `hftrans.merge_fragments_files("Chr1", mito_n, plastid_n)` | Delegates to `annotations.merge_fragments_files()`. |
| `remove_inner_fragments(all_fragments_df)` | impure | `hftrans.remove_inner_fragments(df)` | Delegates to `annotations.remove_inner_fragments()`. |
| `merge_numt_nupt(num_Chr, num_mito, num_plastid)` | impure | `hftrans.merge_numt_nupt(chr_n, mito_n, plastid_n)` | Delegates to `annotations.merge_numt_nupt()`. |
| `get_raw_transfer_groups(num_Chr)` | impure | `hftrans.get_raw_transfer_groups(chr_n)` | Delegates to `annotations.get_raw_transfer_groups()`. |

## variants

BLAST grouping, read-variant classification, subgroup extraction, coverage-read
selection, BCFtools calling, and SNV/indel summaries.

| Function or class | Purity | Typical use | Notes |
| --- | --- | --- | --- |
| `get_tmp_root()` | impure | `tmp_root = hfvar.get_tmp_root()` | Reads environment variables and filesystem state to choose a temp root. |
| `Index_label_alignments` | impure | `indexer = hfvar.Index_label_alignments(df)` | Stateful helper class for labeling indexed alignments. |
| `run_blastn_sorter_single(read_record, ref_fasta, results_filename, soft_paths_dict, tmp_blastn_dir)` | impure | `hfvar.run_blastn_sorter_single(record, ref, "results", soft_paths, tmp_dir)` | Runs BLAST for one read record and writes temporary results. |
| `run_multi_threads_blastn(sample_index, genome, run_info, reads_filename, genome_absolute_path, results_filename, soft_paths_dict, threads)` | impure | `hfvar.run_multi_threads_blastn("sample", "mito", run_info, reads, ref, "results", soft_paths, 8)` | Runs parallel BLAST read alignment and writes sorted results. |
| `get_type_and_subtype(blastn_info_file, default_num_types, out_dir="read_group_files")` | impure | `hfvar.get_type_and_subtype("blastn.txt", 4)` | Reads BLAST table, assigns read types/subtypes, and writes group files. |
| `check_FL_and_multi(blastn_info_file, default_num_types, out_dir="FL_read_group_files", id_dir="IDs", report_dir="Reports")` | impure | `hfvar.check_FL_and_multi("blastn.txt", 4)` | Filters full-length/multi-alignment read groups and writes reports. |
| `get_next_groups(groups)` | pure | `next_groups = hfvar.get_next_groups(groups)` | Returns next subgroup labels from an in-memory group list. |
| `match_se1_ss2(old_index, num_align, subtype, SE1, SS2, blastn_df, ID_subgroup_dir="ID_subgroup")` | impure | `hfvar.match_se1_ss2(i, n, subtype, se1, ss2, df)` | Matches two-boundary structural patterns and writes subgroup IDs. |
| `match_se1_ss2_se2_ss3(old_index, num_align, subtype, SE1, SS2, SE2, SS3, blastn_df, ID_subgroup_dir="ID_subgroup")` | impure | `hfvar.match_se1_ss2_se2_ss3(i, n, subtype, se1, ss2, se2, ss3, df)` | Matches four-boundary patterns and writes subgroup IDs. |
| `match_se1_ss2_se2_ss3_se3_ss4(old_index, num_align, subtype, SE1, SS2, SE2, SS3, SE3, SS4, blastn_df, ID_subgroup_dir="ID_subgroup")` | impure | `hfvar.match_se1_ss2_se2_ss3_se3_ss4(i, n, subtype, se1, ss2, se2, ss3, se3, ss4, df)` | Matches six-boundary patterns and writes subgroup IDs. |
| `match_se1_ss2_se2_ss3_se3_ss4_se4_ss5(old_index, num_align, subtype, SE1, SS2, SE2, SS3, SE3, SS4, SE4, SS5, blastn_df, ID_subgroup_dir="ID_subgroup")` | impure | `hfvar.match_se1_ss2_se2_ss3_se3_ss4_se4_ss5(i, n, subtype, se1, ss2, se2, ss3, se3, ss4, se4, ss5, df)` | Matches eight-boundary patterns and writes subgroup IDs. |
| `get_subgroups(groups, num_align, input_dir="FL_read_group_files", out_dir="combined_excel", ID_subgroup_dir="ID_subgroup")` | impure | `hfvar.get_subgroups(groups, 2)` | Reads grouped read files and writes subgroup workbooks/ID files. |
| `summarize_blastn_results(blastn_result_file="all_sorted_blastn_alignments.txt")` | impure | `hfvar.summarize_blastn_results("all_sorted_blastn_alignments.txt")` | Reads BLAST result files and writes summary tables. |
| `get_cov_reads(current_dir, IDs_dir, soft_paths_dict, genome_absolute_path)` | impure | `hfvar.get_cov_reads(".", "IDs", soft_paths, ref)` | Extracts coverage reads for downstream variant analysis. |
| `run_bcftools(read_record, ref_fasta, sample_platform, soft_paths_dict, tmp_bcftools_dir)` | impure | `hfvar.run_bcftools(record, ref, "hifi", soft_paths, tmp_dir)` | Runs BCFtools pipeline for one read or read group. |
| `run_multi_threads_bcftools(sample_index, genome, run_info, reads_filename, genome_absolute_path, sample_platform, results_filename, soft_paths_dict, threads)` | impure | `hfvar.run_multi_threads_bcftools("sample", "mito", run_info, reads, ref, "hifi", "results", soft_paths, 8)` | Runs BCFtools in parallel and writes variant tables. |
| `snv_or_indel(input_file, output_file)` | impure | `hfvar.snv_or_indel("variants.tsv", "typed.tsv")` | Reads variant calls and writes SNV/indel labels. |
