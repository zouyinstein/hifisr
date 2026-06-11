# HiFiSR Semi-Automatic Workflow

This Snakemake workflow reruns the `W3-5-2` local test dataset from the
`hifisr` project root. It keeps the two manual decisions explicit:

1. Linearizing GFA files with GFA_Editor and saving edited FASTA files.
2. Reviewing coverage/bubble plots and writing `pos_ref_alt.txt` for corrected-genome verification.

The default config is `workflow/config/w3_5_2.yaml`. On Ubuntu, the workflow is
intended to run with a server-specific `deps/soft_paths.txt`; no macOS-specific
setup rule is used.

The default config sets `project_dir: "{snakefile_dir}"` and
`results_dir: "{project_dir}/results"`, so paths are explicit without hardcoding
a machine-specific absolute project path.

## Required Inputs

Run commands from the `hifisr` project root and prepare:

- HiFi reads: `data/W3-5-2.fastq.gz`
  ([CRR384154.fastq.gz](https://download.cncb.ac.cn/gsa/CRA006060/CRR384154/CRR384154.fastq.gz))
- Mitochondrial reference: `references/Col_mito.fa`
- Plastid reference: `references/Col_plastid.fa`
- Software path file: `deps/soft_paths.txt`

The `soft_paths.txt` file must be tab-delimited and include the executables
used by the workflow, including `python`, `minimap2`, `samtools`, `seqkit`,
`mecat`, `blastn`, `bcftools`, `bamtools`, `pigz`, `bandage`, `hifiasm`,
`flye`, and `canu`. The `python` entry should point to an environment where
Snakemake and the HiFiSR Python dependencies are installed. The workflow does
not require installing the `hifisr` package into the virtual environment;
scripts under `analysis_scripts/` import the local source code from
`hifisr_functions/` via `analysis_scripts/_bootstrap.py`.

If the file is stored elsewhere, override the config:

```bash
python -m snakemake --cores 8 --config soft_paths=/path/to/soft_paths.txt final
```

## Main Targets

Run from the `hifisr` project root:

```bash
cd /path/to/hifisr
python -m pip install snakemake
python -m snakemake --cores 8 references_ready
python -m snakemake --cores 8 reads_ready
python -m snakemake --cores 8 draft_for_manual_edit
python -m snakemake --cores 8 polish_alignment_variant
python -m snakemake --cores 8 polish_alignment_variant_review_inputs
python -m snakemake --cores 8 verify_corrected_genome
python -m snakemake --cores 8 final
```

After reviewing the final outputs, clean large regenerable intermediates:

```bash
python -m snakemake --cores 1 clean
```

If `snakemake` is available on `PATH`, the equivalent command is:

```bash
snakemake clean
```

Use `-n` for a dry run:

```bash
python -m snakemake -n --cores 8 final
```

## Endpoint Resume

Snakemake resumes from file endpoints. For example:

- `references_ready` checks that workflow references are available.
- `reads_ready` stops after mito/plastid read extraction, filtering, and sampling.
- `draft_for_manual_edit` stops after draft GFA and PNG files are ready.
- `polish_alignment_variant` runs polish/alignment plus read-variant analysis for mito and plastid.
- `polish_alignment_variant_review_inputs` verifies the polish/alignment variant
  outputs needed for manual review exist.
- `verify_corrected_genome` runs only the genomes listed in `run3.genomes`.
- `final` returns run_3 output for genomes listed in `run3.genomes`, otherwise run_2 output.

## Manual Breakpoint 1: GFA_Editor

After `draft_for_manual_edit`, inspect/edit the GFA files and save linear FASTA
files at the paths in config:

```yaml
draft_edited:
  mito: "{results_dir}/{sample}/draft_assembly/mito/all_mito_500K_after_rr.edited.fasta"
  plastid: "{results_dir}/{sample}/draft_assembly/plastid/all_plastid_150K_after_rr.edited.fasta"
```

Then continue:

```bash
python -m snakemake --cores 8 polish_alignment_variant
```

## Manual Breakpoint 2: Coverage And Bubble Review

After `polish_alignment_variant` or `polish_alignment_variant_review_inputs`, inspect:

- `{results_dir}/{sample}/{genome}/run_2/coverage_*.png`
- `{results_dir}/{sample}/{genome}/run_2/bubble_type_2_rep_raw.pdf`
- `{results_dir}/{sample}/{genome}/run_2/variants_anno_combined_depth_frq_filter.xlsx`

Write corrections as tab-separated `pos ref alt` rows:

```text
322485	C	A
344025	G	T
```

Save them at the paths in config:

```yaml
pos_ref_alt:
  mito: "{results_dir}/{sample}/draft_assembly/mito/pos_ref_alt.txt"
  plastid: "{results_dir}/{sample}/draft_assembly/plastid/pos_ref_alt.txt"
```

Use an empty `pos_ref_alt.txt` for a genome when no correction is needed. Set
which genomes should run corrected-genome verification:

```yaml
run3:
  name: run_3
  genomes:
    - mito
```

Then continue:

```bash
python -m snakemake --cores 8 verify_corrected_genome
```

## Cleaning Large Intermediates

The `clean` target keeps result-interpretation files, including final and
intermediate Excel tables, coverage plots, bubble plots, edited/corrected FASTA
files, `pos_ref_alt.txt`, small read-statistics files, and logs.

It removes large regenerable FASTQ/FASTA intermediates:

- `results/{sample}/reads/{sample}_{genome}.fastq`
- `results/{sample}/reads/sample_reads/sample_4000_{genome}.fastq`
- `results/{sample}/{genome}/{run}/reads.fasta`
- `results/{sample}/{genome}/{run}/new_reads.fasta`
- `results/{sample}/{genome}/{run}/FL.fasta`
- `results/{sample}/{genome}/{run}/partial.fasta`
- `results/{sample}/{genome}/{run}/variant_cov.fasta`
- `results/.tmp/`
- `results/.matplotlib/`

For the default `W3-5-2` config, `{genome}` expands to `mito` and `plastid`.
The `{run}` paths include `run_2` for both genomes and `run_3` for genomes
listed in `run3.genomes` (`mito` by default).

To add or remove files from the clean list, edit `Snakefile`:

- Add or remove file paths in `clean_intermediate_paths()`.
- Add or remove removable directories in `clean_intermediate_dirs()`.

After changing the clean list, update this README section and verify the target
before deleting files:

```bash
python -m snakemake -n --cores 1 clean
```

After `clean`, the kept files remain under `results/`, but rerunning upstream
analysis steps may require regenerating the removed intermediates.

## Output Layout

All workflow products, logs, caches, and temporary files are written under
`results/`:

```text
results/W3-5-2/reads/
results/W3-5-2/draft_assembly/
results/W3-5-2/mito/run_2/
results/W3-5-2/mito/run_3/
results/W3-5-2/plastid/run_2/
results/W3-5-2/logs/snakemake/
results/.tmp/
results/.cache/
results/.matplotlib/
```

Input reads, reference FASTA files, `deps/soft_paths.txt`, workflow metadata,
and source code remain outside `results/`.

## Reference Rotation

For the current W3-5-2 local rerun, `Col_mito.fa` and `Col_plastid.fa` already
contain the rotated references, so the default config has:

```yaml
references:
  rotate: false
```

For a fresh run from raw selected references, set:

```yaml
references:
  rotate: true
  mito: path/to/raw_mito.fa
  plastid: path/to/raw_plastid.fa
```

The workflow writes stable rotated copies under:

```text
{results_dir}/{sample}/references/mito_rotated.fasta
{results_dir}/{sample}/references/plastid_rotated.fasta
```

It does not overwrite the source reference files.
