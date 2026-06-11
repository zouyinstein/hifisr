# HiFiSR Semi-Automatic Workflow

This Snakemake workflow reruns the `W3-5-2` local test dataset from the
`hifisr` project root. It keeps the two manual decisions explicit:

1. Linearizing GFA files with GFA_Editor and saving edited FASTA files.
2. Reviewing coverage/bubble plots and writing `pos_ref_alt.txt` for corrected-genome verification.

The bundled W3-5-2 example config is `workflow/config/w3_5_2_macOS.yaml`. On
Ubuntu, copy it to a server-specific config and set `soft_paths` to the
server's `deps/soft_paths.txt`; no macOS-specific setup rule is used.

The default config sets `project_dir: "{snakefile_dir}"` and
`results_dir: "{project_dir}/results"`, so paths are explicit without hardcoding
a machine-specific absolute project path.

## Required Inputs

Run commands from the `hifisr` project root and prepare:

- HiFi reads: `data/W3-5-2.fastq.gz`
  ([CRR384154.fastq.gz](https://download.cncb.ac.cn/gsa/CRA006060/CRR384154/CRR384154.fastq.gz))
- Mitochondrial reference: `references/Col_mito.fa`
- Plastid reference: `references/Col_plastid.fa`
- Software path file: `deps/soft_paths_macOS.txt` for the bundled macOS config,
  or `deps/soft_paths.txt` for Ubuntu/server configs

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

On Ubuntu, install the system packages needed to create virtual environments
and provide Python's `sqlite3` module, which Snakemake requires:

```bash
sudo apt update
sudo apt install -y python3-venv sqlite3 libsqlite3-dev
```

Create and activate the Python environment from the `hifisr` project root
before running Snakemake:

```bash
cd /path/to/hifisr
python3 -m venv .venv
source .venv/bin/activate
python -c "import sqlite3; print(sqlite3.sqlite_version)"
python -m pip install --upgrade pip
python -m pip install -r requirements-dev.txt
```

The `source .venv/bin/activate` command makes `python`, `pip`, and
`python -m snakemake` use this environment in the current shell session. Confirm
the active environment and Snakemake installation with:

```bash
which python
python -m snakemake --version
```

If `python -c "import sqlite3"` fails with `No module named '_sqlite3'`, rebuild
the virtual environment with a Python executable that has sqlite support, for
example the system Python that passes the same import check.

On Python 3.10, `pip` may install an older Snakemake release. That release
expects the older PuLP `list_solvers` API, so `requirements-dev.txt` pins
`pulp==2.7.0`. If Snakemake fails with
`AttributeError: module 'pulp' has no attribute 'list_solvers'`, reinstall the
dependencies in the active environment:

```bash
python -m pip install --force-reinstall "pulp==2.7.0"
python -m pip install -r requirements-dev.txt
python -m snakemake --version
```

The default software path file location is:

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

## Linux Server Run With Explicit Config

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

Use `--cores`, not `--core`, when setting the number of CPU cores.
On network filesystems or busy servers, add `--latency-wait 60` or a larger
value if Snakemake reports that a job finished but recently created output
files are not visible yet:

```bash
python -m snakemake --configfile workflow/config/w3_5_2_macOS.yaml \
  --cores 8 --latency-wait 60 draft_for_manual_edit
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

## Running Other Samples

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
