# hifisr

HiFi-SR is a Python-based workflow for organellar genome read extraction,
draft assembly, polishing, alignment, and structural-variant review from HiFi
reads.

***We will continuously make upgrades and modifications to hifisr to enhance its functionality and performance.***

## Installation

Detailed environment setup, PyPI mirror installation, and third-party software
installation commands are maintained in [docs/installation.md](docs/installation.md).
The Snakemake workflow runs directly from this source checkout; installing the
`hifisr` package itself is not required.

Function usage notes and pure/impure markers for `hifisr_functions/` are in
[docs/hifisr_functions_reference.md](docs/hifisr_functions_reference.md).

## Dependency Coverage

The current Snakemake workflow has two dependency layers:

- Python packages: `requirements-dev.txt` and `environment.yml` cover the Python
  environment, including Snakemake and the analysis libraries.
- External command-line tools: `soft_paths.txt` points the workflow to the
  server-specific executable paths.

The Conda environment is pinned to Python 3.11; other Python versions may also
work, but they have not been tested for this workflow.

The workflow-facing external tools are `python`, `minimap2`, `samtools`,
`seqkit`, `mecat`, `blastn`, `bcftools`, `bamtools`, `pigz`, `bandage`,
`hifiasm`, `flye`, and `canu`. The macOS example may also include legacy
entries such as `meryl` and `winnowmap`; these are not required by the current
W3-5-2 Snakemake workflow.

## soft_paths Configuration

The default server path file is:

```text
deps/soft_paths.txt
```

The file must be tab-delimited: first column is the software name, second column
is the absolute path to the executable. For Linux server runs, create this file
locally on the server and keep it out of Git.

```text
python	/path/to/hifisr/.venv/bin/python
minimap2	/path/to/minimap2
samtools	/path/to/samtools
seqkit	/path/to/seqkit
mecat	/path/to/mecat.pl
blastn	/path/to/blastn
bcftools	/path/to/bcftools
bamtools	/path/to/bamtools
pigz	/path/to/pigz
bandage	/path/to/Bandage
hifiasm	/path/to/hifiasm
flye	/path/to/flye
canu	/path/to/canu
```

For Conda/Mamba environments, use [deps/soft_paths_conda.txt](deps/soft_paths_conda.txt)
as a template. The `python` entry should match `which python` after activating
the Conda environment:

```text
python	/path/to/miniconda3/envs/hifisr/bin/python
```

`environment.yml` installs the Conda-available third-party tools into the Conda
environment: `minimap2`, `samtools`, `seqkit`, `blast`/`blastn`, `bcftools`,
`bamtools`, `pigz`, `bandage`, `hifiasm`, `flye`, and `canu`. `mecat` still
needs a manual path in `soft_paths.txt`.

The repository also includes a macOS local-development example:

```text
deps/soft_paths_macOS.txt
```

Do not use the macOS file on Ubuntu unless every path has been replaced with
server executables. To use a non-default file, pass it explicitly:

```bash
python -m snakemake --configfile workflow/config/w3_5_2_linux.yaml \
  --config soft_paths=deps/soft_paths.txt \
  --cores 8 final
```

## quick_check

Before running a workflow on a server, check both the `soft_paths` executables
and the Python packages installed in the workflow Python environment.

Third-party tools checked from `soft_paths`: `python`, `minimap2`, `samtools`,
`seqkit`, `mecat`, `blastn`, `bcftools`, `bamtools`, `pigz`, `bandage`,
`hifiasm`, `flye`, and `canu`.

Python packages checked from `requirements-dev.txt`: `biopython`, `pysam`,
`pandas`, `numpy`, `openpyxl`, `xlsxwriter`, `matplotlib`, `polars`,
`fastexcel`, `pulp==2.7.0`, `snakemake`, and `pytest`.

```bash
python -m snakemake --configfile workflow/config/w3_5_2_linux.yaml \
  --config soft_paths=deps/soft_paths.txt \
  --cores 1 check_runtime_dependencies
```

## Local Snakemake Test Workflow

The repository also contains a semi-automatic Snakemake workflow for rerunning the
`W3-5-2` local test dataset. The workflow now lives inside the `hifisr`
project directory. The bundled W3-5-2 example config is
`workflow/config/w3_5_2_macOS.yaml`; copy it for Ubuntu/server runs and replace
the machine-specific paths. Analysis products are written under `results/`.

### Required input files

Run from the `hifisr` project root and prepare the following files before
starting from scratch:

- HiFi reads: `data/W3-5-2.fastq.gz`
  ([CRR384154.fastq.gz](https://download.cncb.ac.cn/gsa/CRA006060/CRR384154/CRR384154.fastq.gz))
- Mitochondrial reference: `references/Col_mito.fa`
- Plastid reference: `references/Col_plastid.fa`
- Software path file: `deps/soft_paths_macOS.txt` for the bundled macOS config,
  or `deps/soft_paths.txt` for Ubuntu/server configs

The `soft_paths.txt` file must be tab-delimited and include the paths for
`python`, `minimap2`, `samtools`, `seqkit`, `mecat`, `blastn`, `bcftools`,
`bamtools`, `pigz`, `bandage`, `hifiasm`, `flye`, and `canu`. The workflow uses
this file to derive executable paths and PATH entries, so on an Ubuntu server
the only machine-specific workflow file should be the corresponding
`soft_paths.txt`. The `python` entry should point to an environment where
Snakemake and the Python package dependencies listed above are installed. It
does not need an installed `hifisr` package because the workflow imports
`hifisr_functions/` from the local source checkout. If the file is stored
elsewhere, override it with `--config soft_paths=/path/to/soft_paths.txt`.

The default location is:

```text
/path/to/hifisr/deps/soft_paths.txt
```

The repository may contain or generate a macOS local-development example named
`deps/soft_paths_macOS.txt`. Do not use that file on Ubuntu unless every path in
it has been replaced with Ubuntu executables. For server runs, create
`deps/soft_paths.txt` yourself or pass a server-specific file with
`--config soft_paths=...`.
For macOS local reruns, either copy the macOS example to `deps/soft_paths.txt`
or pass it explicitly:

```bash
python -m snakemake --config soft_paths=deps/soft_paths_macOS.txt --cores 8 final
```

Example contents:

```text
python	/path/to/hifisr/.venv/bin/python
minimap2	/path/to/minimap2
samtools	/path/to/samtools
seqkit	/path/to/seqkit
mecat	/path/to/mecat.pl
blastn	/path/to/blastn
bcftools	/path/to/bcftools
bamtools	/path/to/bamtools
pigz	/path/to/pigz
bandage	/path/to/Bandage
hifiasm	/path/to/hifiasm
flye	/path/to/flye
canu	/path/to/canu
```

Replace every `/path/to/...` entry with the real Ubuntu executable path. For
example, if Flye is installed in the active virtual environment, find it with
`which flye` and set:

```text
flye	/path/to/hifisr/.venv/bin/flye
```

If `draft_for_manual_edit` fails with `/bin/sh: 1: flye: not found`, the `flye`
line is missing, misspelled, or points to a command that is not on `PATH`.
Update `deps/soft_paths.txt` or pass a server-specific file with
`--config soft_paths=/path/to/soft_paths.txt`.

If `soft_paths.txt` is saved outside `deps/`, either specify it in the sample
config:

```yaml
soft_paths: "/home/user/software_paths/soft_paths.txt"
```

or override it at runtime:

```bash
python -m snakemake --configfile workflow/config/my_sample.yaml \
  --config soft_paths=/home/user/software_paths/soft_paths.txt \
  --cores 8 final
```

### Linux Server Run With Explicit Config

On Linux servers, avoid using the bundled macOS paths. Create a server config
and a server `soft_paths.txt`, then pass both explicitly in every Snakemake
command:

```bash
cd /path/to/hifisr
cp workflow/config/w3_5_2_macOS.yaml workflow/config/w3_5_2_linux.yaml
```

Edit `workflow/config/w3_5_2_linux.yaml` so machine-specific inputs point to
Linux files:

```yaml
reads: data/W3-5-2.fastq.gz
soft_paths: "{project_dir}/deps/soft_paths.txt"

references:
  rotate: false
  mito: "{project_dir}/references/Col_mito.fa"
  plastid: "{project_dir}/references/Col_plastid.fa"
```

Create the Linux software path file:

```bash
mkdir -p deps
nano deps/soft_paths.txt
```

Run with both files specified:

```bash
python -m snakemake \
  --configfile workflow/config/w3_5_2_linux.yaml \
  --config soft_paths=deps/soft_paths.txt \
  --cores 8 \
  references_ready
```

Use the same `--configfile` and `--config soft_paths=...` options for later
targets, for example:

```bash
python -m snakemake --configfile workflow/config/w3_5_2_linux.yaml \
  --config soft_paths=deps/soft_paths.txt --cores 8 reads_ready
python -m snakemake --configfile workflow/config/w3_5_2_linux.yaml \
  --config soft_paths=deps/soft_paths.txt --cores 8 draft_for_manual_edit
python -m snakemake --configfile workflow/config/w3_5_2_linux.yaml \
  --config soft_paths=deps/soft_paths.txt --cores 8 polish_alignment_variant
python -m snakemake --configfile workflow/config/w3_5_2_linux.yaml \
  --config soft_paths=deps/soft_paths.txt --cores 8 verify_corrected_genome
python -m snakemake --configfile workflow/config/w3_5_2_linux.yaml \
  --config soft_paths=deps/soft_paths.txt --cores 8 final
```

The current local test configuration assumes that `Col_mito.fa` and
`Col_plastid.fa` are already in the desired rotated state, so
`references.rotate` is set to `false`.

### Run commands

Run the workflow from the `hifisr` project root. The Python used for Snakemake
should be the same environment listed as `python` in `deps/soft_paths.txt`; the
workflow scripts will import local code from `hifisr_functions/`.

```bash
cd /path/to/hifisr
python -m snakemake --cores 8 references_ready
python -m snakemake --cores 8 reads_ready
python -m snakemake --cores 8 draft_for_manual_edit
```

After `draft_for_manual_edit`, inspect and linearize the GFA files manually,
for example with GFA_Editor, and save the edited FASTA files at:

```text
results/W3-5-2/draft_assembly/mito/all_mito_500K_after_rr.edited.fasta
results/W3-5-2/draft_assembly/plastid/all_plastid_150K_after_rr.edited.fasta
```

Then continue with polish/alignment and variant calling:

```bash
python -m snakemake --cores 8 polish_alignment_variant
```

The `polish_alignment_variant` target generates both mitochondrial and plastid
`run_2` outputs.
Run this target if you want to inspect mitochondrial `run_2` before applying
manual reference corrections.

After `polish_alignment_variant`, review the coverage plots, bubble plots, and
filtered variant tables:

```text
results/W3-5-2/{genome}/run_2/coverage_*.png
results/W3-5-2/{genome}/run_2/bubble_type_2_rep_raw.pdf
results/W3-5-2/{genome}/run_2/variants_anno_combined_depth_frq_filter.xlsx
```

Write manual corrections as tab-delimited `pos ref alt` rows. Use an empty file
when no correction is needed:

```text
results/W3-5-2/draft_assembly/mito/pos_ref_alt.txt
results/W3-5-2/draft_assembly/plastid/pos_ref_alt.txt
```

Then run the corrected-genome verification and final targets:

```bash
python -m snakemake --cores 8 verify_corrected_genome
python -m snakemake --cores 8 final
```

Use `--cores`, not `--core`, when setting the number of CPU cores.
On network filesystems or busy servers, add `--latency-wait 60` or a larger
value if Snakemake reports that a job finished but recently created output
files are not visible yet:

```bash
python -m snakemake --configfile workflow/config/w3_5_2_macOS.yaml \
  --cores 8 --latency-wait 60 draft_for_manual_edit
```

After `final` has completed and the outputs have been reviewed, remove large
regenerable intermediates with:

```bash
python -m snakemake --cores 1 clean
```

If `snakemake` is available on `PATH`, the equivalent command is:

```bash
snakemake clean
```

The `clean` target keeps the result-interpretation files, including final and
intermediate Excel tables, coverage plots, bubble plots, edited/corrected FASTA
files, `pos_ref_alt.txt`, small read-statistics files, and logs.

It deletes only the following regenerable files and directories under
`results/`:

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

Use `-n` for a dry run:

```bash
python -m snakemake -n --cores 8 final
```

### Output layout

Analysis products, Snakemake logs, run-local caches, and temporary files are
written under `results/`, including:

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

The `final` target uses `run_3` for genomes listed in `run3.genomes` and
otherwise uses `run_2`. In the current `W3-5-2` test configuration, `run_3` is
enabled for mitochondria only, so the final target depends on:

```text
results/W3-5-2/mito/run_3/variants_anno_combined_depth_frq_filter.xlsx
results/W3-5-2/plastid/run_2/variants_anno_combined_depth_frq_filter.xlsx
```

Input reads, reference FASTA files, the virtual environment, workflow metadata,
and source code remain outside `results/`.

After `clean`, the kept files remain under `results/`, but rerunning upstream
analysis steps may require regenerating the removed intermediates.

### Running Other Samples

Use one config file per sample. Copy the W3-5-2 example config and edit the
sample-specific paths; do not edit `Snakefile` for routine sample runs:

```bash
cd /path/to/hifisr
cp workflow/config/w3_5_2_macOS.yaml workflow/config/my_sample.yaml
```

For Ubuntu/server runs, set `soft_paths` in the copied config to the server
file, usually `"{project_dir}/deps/soft_paths.txt"`.

Prepare the input files referenced by the new config:

```text
data/MySample.fastq.gz
references/MySample_mito.fa
references/MySample_plastid.fa
deps/soft_paths.txt
```

The config should set the sample name, reads, reference FASTA files, and manual
review paths. A minimal template is:

```yaml
sample: MySample
project_dir: "{snakefile_dir}"
results_dir: "{project_dir}/results"

reads: data/MySample.fastq.gz
soft_paths: "{project_dir}/deps/soft_paths.txt"

genomes:
  - mito
  - plastid

references:
  rotate: false
  mito: "{project_dir}/references/MySample_mito.fa"
  plastid: "{project_dir}/references/MySample_plastid.fa"

draft_edited:
  mito: "{results_dir}/{sample}/draft_assembly/mito/all_mito_500K_after_rr.edited.fasta"
  plastid: "{results_dir}/{sample}/draft_assembly/plastid/all_plastid_150K_after_rr.edited.fasta"

run2:
  name: run_2

pos_ref_alt:
  mito: "{results_dir}/{sample}/draft_assembly/mito/pos_ref_alt.txt"
  plastid: "{results_dir}/{sample}/draft_assembly/plastid/pos_ref_alt.txt"

run3:
  name: run_3
  genomes:
    - mito
```

Run the workflow with `--configfile`:

```bash
python -m snakemake --configfile workflow/config/my_sample.yaml \
  --config soft_paths=deps/soft_paths.txt --cores 8 references_ready
python -m snakemake --configfile workflow/config/my_sample.yaml \
  --config soft_paths=deps/soft_paths.txt --cores 8 reads_ready
python -m snakemake --configfile workflow/config/my_sample.yaml \
  --config soft_paths=deps/soft_paths.txt --cores 8 draft_for_manual_edit
```

After `draft_for_manual_edit`, inspect and linearize the GFA files manually and
save edited FASTA files at the paths listed in `draft_edited`. Then continue:

```bash
python -m snakemake --configfile workflow/config/my_sample.yaml \
  --config soft_paths=deps/soft_paths.txt --cores 8 polish_alignment_variant
```

Review the `run_2` coverage plots, bubble plots, and filtered variant tables.
Write `pos_ref_alt.txt` files at the paths listed in the config; create an empty
file for any genome with no manual correction. Then run:

```bash
python -m snakemake --configfile workflow/config/my_sample.yaml \
  --config soft_paths=deps/soft_paths.txt --cores 8 verify_corrected_genome
python -m snakemake --configfile workflow/config/my_sample.yaml \
  --config soft_paths=deps/soft_paths.txt --cores 8 final
```

Outputs are written under:

```text
results/MySample/
```

If a sample should run only one organellar genome, edit `genomes`. If no
corrected-genome verification is needed, keep `run3.name` but set
`run3.genomes` to an empty list:

```yaml
run3:
  name: run_3
  genomes: []
```

## Exmaples

### Example 1: Wild-type *Arabidopsis thaliana* (Col-0 ecotype)

#### Prepare the reference sequences of mitochondrial and plastid genomes

```bash
# create your working folder
mkdir -p /mnt/software/scripts/results && cd /mnt/software/scripts/results
# create a folder ref to store the references
mkdir ref && cd ref
# Download the references Col_mito.fa, Col_plastid.fa.
# rotate the reference
python adjust_ref_fasta.py /mnt/software/scripts/hifisr/deps/soft_paths.txt mito Col_mito.fa
python adjust_ref_fasta.py /mnt/software/scripts/hifisr/deps/soft_paths.txt plastid Col_plastid.fa
```

#### Prepare the reads (fastq format)

Analysis of an example wild-type *Arabidopsis thaliana* dataset Col-CEN ([ERR6210723](https://www.ncbi.nlm.nih.gov/sra/ERR6210723), 14.6 Gb, [Naish et al., 2021, Science](https://www.science.org/doi/10.1126/science.abi7489))

```bash
# Download the data as Col-CEN.fastq
python get_mtpt_reads.py /mnt/software/scripts/hifisr/deps/soft_paths.txt ATHiFi001 /mnt/software/scripts/results/mito_rotated_293434.fasta /mnt/software/scripts/results/plastid_rotated_61049.fasta /mnt/software/scripts/results/Col-CEN.fastq 32
```

| All                                                                                                                           | mitochondria                                                                                                                            | plastid                                                                                                                               |
| ----------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| ![All](001.植物细胞器基因组组装后校验/版本A/细胞器代码库/hifisr-dev/hifisr/examples/example_1/all_length_qual_distribution.jpg) | ![mitochondria](001.植物细胞器基因组组装后校验/版本A/细胞器代码库/hifisr-dev/hifisr/examples/example_1/mito_length_qual_distribution.jpg) | ![plastid](001.植物细胞器基因组组装后校验/版本A/细胞器代码库/hifisr-dev/hifisr/examples/example_1/plastid_length_qual_distribution.jpg) |

#### Calling and calculating the frequencies of SVs, SNVs and small InDels with plotting

```bash
# filter the reads: remove reads shorter than 10 kb
python filt_read_ids.py /mnt/software/scripts/hifisr/deps/soft_paths.txt ATHiFi001 mito_id_length_qual.txt plastid_id_length_qual.txt 10000 0

# random sampling of 4000 reads for mitochondrial and plastid genomes, respectively
python sample_reads.py /mnt/software/scripts/hifisr/deps/soft_paths.txt ATHiFi001 mito filt_L10K_mito_id_length_qual.txt 4000
python sample_reads.py /mnt/software/scripts/hifisr/deps/soft_paths.txt ATHiFi001 plastid filt_L10K_plastid_id_length_qual.txt 4000

# Estimation of variant frequencies
python get_variants_in_reads.py /mnt/software/scripts/hifisr/deps/soft_paths.txt ATHiFi001 mito run_1 /mnt/software/scripts/results/mito_rotated_293434.fasta /mnt/software/scripts/results/ATHiFi001/reads/sample_reads/sample_4000_mito.fastq 32
python get_variants_in_reads.py /mnt/software/scripts/hifisr/deps/soft_paths.txt ATHiFi001 plastid run_1 /mnt/software/scripts/results/plastid_rotated_61049.fasta /mnt/software/scripts/results/ATHiFi001/reads/sample_reads/sample_4000_plastid.fastq 32
```

| mitochondria                                        | plastid                                                |
| --------------------------------------------------- | ------------------------------------------------------ |
| ![One-rearrangements](mito_bubble_type_2_rep_raw.jpg) | ![One-rearrangements](plastid_bubble_type_2_rep_raw.jpg) |
| ![Coverage](mito_coverage_plot.jpg)                   | ![Coverage](plastid_coverage_plot.jpg)                   |

### Example 2: *Amborella trichopoda* (var. Santa Cruz 75)

#### Prepare the reference sequences of mitochondrial and plastid genomes

```bash
# create your working folder
mkdir -p /mnt/software/scripts/results && cd /mnt/software/scripts/results
# create a folder ref to store the references
mkdir ref && cd ref
# Download the references for mt genome: KF754799.1_mt_5.fasta, KF754800.1_mt_4.fasta, KF754801.1_mt_3.fasta, KF754802.1_mt_2.fasta, KF754803.1_mt_1.fasta.
# Concatenate the files 
cat KF754799.1_mt_5.fasta KF754800.1_mt_4.fasta KF754801.1_mt_3.fasta KF754802.1_mt_2.fasta KF754803.1_mt_1.fasta > mt_all.fasta
# Download the references for pt genome: AJ506156.2_pt.fasta
# rotate the pt reference
python adjust_ref_fasta.py /mnt/software/scripts/hifisr/deps/soft_paths.txt plastid AJ506156.2_pt.fasta
```

#### Prepare the reads (fastq format)

Analysis of the *Amborella trichopoda* dataset (var. Santa Cruz 75) ([SRR28888927](https://www.ncbi.nlm.nih.gov/sra/SRR28888927), 33.1 Gb, [Carey et al., 2025, Nature Plants](https://www.nature.com/articles/s41477-024-01858-x))

```bash
# Download the data as SRR28888927.1.fastq.gz
python get_mtpt_reads.py /mnt/software/scripts/hifisr/deps/soft_paths.txt Amborella_1 /mnt/software/scripts/results/mt_all.fasta /mnt/software/scripts/results/AJ506156.2_pt.fasta /mnt/software/scripts/results/SRR28888927.1.fastq.gz 32
```

| All                                                                                                                           | mitochondria                                                                                                                            | plastid                                                                                                                               |
| ----------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| ![All](001.植物细胞器基因组组装后校验/版本A/细胞器代码库/hifisr-dev/hifisr/examples/example_2/all_length_qual_distribution.jpg) | ![mitochondria](001.植物细胞器基因组组装后校验/版本A/细胞器代码库/hifisr-dev/hifisr/examples/example_2/mito_length_qual_distribution.jpg) | ![plastid](001.植物细胞器基因组组装后校验/版本A/细胞器代码库/hifisr-dev/hifisr/examples/example_2/plastid_length_qual_distribution.jpg) |

#### Get the draft assembly of mitochondrial genome

```bash
# filter the reads: remove reads shorter than 10 kb
python filt_read_ids.py /mnt/software/scripts/hifisr/deps/soft_paths.txt Amborella_1 mito_id_length_qual.txt plastid_id_length_qual.txt 10000 0

# random sampling of 10000 reads for mt genome
python sample_reads.py /mnt/software/scripts/hifisr/deps/soft_paths.txt Amborella_1 mito filt_L10K_mito_id_length_qual.txt 10000

# get the draft assembly
python get_draft_assembly.py /mnt/software/scripts/hifisr/deps/soft_paths.txt Amborella_1 mito mt_all.fasta sample_10000_mito.fastq 20
```

| MECAT2 + metaFlye                      | metaFlye                             |
| -------------------------------------- | ------------------------------------ |
| ![before](mecat_mito_500K_before_rr.jpg) | ![before](all_mito_500K_before_rr.jpg) |
| ![after](mecat_mito_500K_after_rr.jpg)   | ![after](all_mito_500K_after_rr.jpg)   |

#### Choose the assembly (upper right, metaFlye before repeat resolution, and without running MECAT2) for downstream analysis

```bash
# Manually resolve the graph in Bandage, and save the results as draft.fasta
# Split the reads by contig
python get_reads_by_contig.py /mnt/software/scripts/hifisr/deps/soft_paths.txt Amborella_1 mt_draft_contigs draft.fasta sample_10000_mito.fastq 32

# Check the correctness of the contigs by calling and calculating variant frequencies
python get_variants_in_reads.py /mnt/software/scripts/hifisr/deps/soft_paths.txt Amborella_1 mito contig_1 contigs.part_contig_1.fa contig_1.fastq 32
python get_variants_in_reads.py /mnt/software/scripts/hifisr/deps/soft_paths.txt Amborella_1 mito contig_2 contigs.part_contig_2.fa contig_2.fastq 32
python get_variants_in_reads.py /mnt/software/scripts/hifisr/deps/soft_paths.txt Amborella_1 mito contig_3 contigs.part_contig_3.fa contig_3.fastq 32
python get_variants_in_reads.py /mnt/software/scripts/hifisr/deps/soft_paths.txt Amborella_1 mito contig_4 contigs.part_contig_4.fa contig_4.fastq 32
python get_variants_in_reads.py /mnt/software/scripts/hifisr/deps/soft_paths.txt Amborella_1 mito contig_5 contigs.part_contig_5.fa contig_5.fastq 32

# Manually correct the small errors contig_1, contig_2, contig_3 based on the results
```

#### Correct contig_4

```bash
# Re-assemble contig_4 separately
python get_draft_assembly.py /mnt/software/scripts/hifisr/deps/soft_paths.txt Amborella_1 mito contig_4.fastq 32

# Manually resolve the graph in Bandage, and save the results as contig_4_v2.fasta
# Polish the assembly and align it to KF754801.1_mt_3.fasta
python get_polished_assembly.py /mnt/software/scripts/hifisr/deps/soft_paths.txt Amborella_1 mito KF754801.1_mt_3.fasta contig_4_v2.fasta contig_4.fastq 32

# Check the correctness again
python get_variants_in_reads.py /mnt/software/scripts/hifisr/deps/soft_paths.txt Amborella_1 mito contig_4_v2 mito_flye_polish_aligned.fasta contig_4.fastq 32
```

| Coverage in Round 1                     | Coverage in Round 2                     |
| --------------------------------------- | --------------------------------------- |
| ![Round 1](coverage_contig_4_Round_1.jpg) | ![Round 2](coverage_contig_4_Round_2.jpg) |

#### Correct contig_5

```bash
# Re-assemble contig_5 separately
mkdir -p Amborella_1/hifiasm_results && cd Amborella_1/hifiasm_results
/mnt/software/scripts/hifisr/deps/hifiasm/hifiasm -o contig_5_hifiasm -t 32 contig_5.fastq
cd ../..

# Manually resolve the graph in Bandage, and save the results as contig_5_hifiasm_p_ctg.fasta
# Check the structural correctness of contig_5_hifiasm_p_ctg.fasta
python get_variants_in_reads.py /mnt/software/scripts/hifisr/deps/soft_paths.txt Amborella_1 mito contig_5_hifiasm_p_ctg contig_5_hifiasm_p_ctg.fasta contig_5.fastq 32

# Polish the assembly and align it to KF754799.1_mt_5.fasta, using only the fully-aligned reads
python get_polished_assembly.py /mnt/software/scripts/hifisr/deps/soft_paths.txt Amborella_1 mito KF754803.1_mt_1.fasta contig_5_hifiasm_p_ctg.fasta FL.fasta 32

# Check the correctness of the polished assembly
python get_polished_assembly.py /mnt/software/scripts/hifisr/deps/soft_paths.txt Amborella_1 mito contig_5_hifiasm_p_ctg_run_2 mito_flye_polish_aligned.fasta contig_5.fastq 32
# Manually correct the small errors based on the results, and check the correctness again
python get_polished_assembly.py /mnt/software/scripts/hifisr/deps/soft_paths.txt Amborella_1 mito contig_5_hifiasm_p_ctg_run_3 mito_flye_polish_aligned_cor.fasta contig_5.fastq 32
```

| Coverage in Round 1                     | Coverage in Round 2                     |
| --------------------------------------- | --------------------------------------- |
| ![Round 1](coverage_contig_5_Round_1.jpg) | ![Round 2](coverage_contig_5_Round_2.jpg) |

| Coverage in Round 3                     | Coverage in Round 4                     |
| --------------------------------------- | --------------------------------------- |
| ![Round 3](coverage_contig_5_Round_3.jpg) | ![Round 4](coverage_contig_5_Round_4.jpg) |

## Citations

1. Yi Zou, Weidong Zhu, Daniel B. Sloan, Zhiqiang Wu. (2022). Long-read sequencing characterizes mitochondrial and plastid genome variants in *Arabidopsis msh1* mutants. *The Plant journal* *112* (3), 738–755. https://doi.org/10.1111/tpj.15976
2. Yi Zou, Weidong Zhu, Yingke Hou, Daniel B. Sloan, Zhiqiang Wu. (2025). The evolutionary dynamics of organellar pan-genomes in *Arabidopsis thaliana*. *bioRxiv* 2025.01.20.633836; doi: https://doi.org/10.1101/2025.01.20.633836
