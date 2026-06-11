# HiFiSR Installation Guide

This document contains the detailed environment and third-party software
installation notes. The main `README.md` keeps only the workflow-facing
`soft_paths.txt` configuration.

## Dependency Coverage

The current Snakemake workflow is fully described by two dependency layers:

- Python packages: listed in `requirements-dev.txt` and installed by
  `environment.yml` or `python -m pip install -r requirements-dev.txt`.
- External command-line tools: listed in `deps/soft_paths.txt` or another file
  passed with `--config soft_paths=/path/to/soft_paths.txt`.

The core workflow expects these external tools to be available through
`soft_paths`: `python`, `minimap2`, `samtools`, `seqkit`, `mecat`, `blastn`,
`bcftools`, `bamtools`, `pigz`, `bandage`, `hifiasm`, `flye`, and `canu`.
`meryl` and `winnowmap` may appear in local legacy examples, but they are not
required by the current W3-5-2 Snakemake workflow.

## Option A: Conda/Mamba Environment

Conda or Mamba is recommended on Linux servers because it provides a Python
build with `sqlite3` support without rebuilding system Python.
The provided `environment.yml` uses Python 3.11; other Python versions may also
work, but they have not been tested for this workflow.

```bash
cd /path/to/hifisr
conda env create -f environment.yml
conda activate hifisr
```

If you use Mamba:

```bash
mamba env create -f environment.yml
conda activate hifisr
```

The Conda environment installs the Python runtime plus the third-party tools
that are available from conda-forge/Bioconda:

```text
minimap2
samtools
seqkit
blast      # provides blastn
bcftools
bamtools
pigz
bandage    # executable is usually Bandage
hifiasm
flye
canu
```

`mecat` is not installed by `environment.yml`; install MECAT2 manually and keep
its `mecat.pl` path in `deps/soft_paths.txt`.

Confirm that the active Python and Snakemake are usable:

```bash
which python
python -c "import sqlite3; print(sqlite3.sqlite_version)"
python -m snakemake --version
```

Set the `python` line in the server `deps/soft_paths.txt` to the active Conda
environment Python reported by `which python`:

```text
python	/path/to/miniconda3/envs/hifisr/bin/python
```

For Conda-installed tools, the other `soft_paths` entries should usually point
inside the same environment, for example:

```text
minimap2	/path/to/miniconda3/envs/hifisr/bin/minimap2
blastn	/path/to/miniconda3/envs/hifisr/bin/blastn
bandage	/path/to/miniconda3/envs/hifisr/bin/Bandage
```

Use `which minimap2`, `which blastn`, or `which Bandage` after activating the
environment to confirm the exact executable paths.

After pulling or syncing new code, update the existing Conda environment with:

```bash
conda activate hifisr
python -m pip install -r requirements-dev.txt
```

## Option B: Python venv Environment

Use this path when Conda/Mamba is not available.

### Create and activate a venv environment

On Ubuntu, install the system packages needed to create virtual environments
and provide Python's `sqlite3` module, which Snakemake requires:

```bash
sudo apt update
sudo apt install -y python3-venv sqlite3 libsqlite3-dev
```

Create and activate a virtual environment:

```bash
cd /path/to/hifisr
python3 -m venv .venv
source .venv/bin/activate
python -c "import sqlite3; print(sqlite3.sqlite_version)"
python -m pip install --upgrade pip
python -m pip install -r requirements-dev.txt
```

The `source .venv/bin/activate` command makes `python`, `pip`, and
`python -m snakemake` use this environment in the current shell session.

### Install packages using your local mirror of PyPI

If a local mirror is required, install the Python dependencies with:

```bash
python -m pip install -i https://mirrors.aliyun.com/pypi/simple/ -r requirements-dev.txt
```

For manual package installation without `requirements-dev.txt`, the required
Python packages are:

```bash
python -m pip install -i https://mirrors.aliyun.com/pypi/simple/ \
  biopython pysam pandas numpy openpyxl xlsxwriter matplotlib polars fastexcel \
  "pulp==2.7.0" snakemake pytest
```

If Snakemake fails with `No module named '_sqlite3'`, rebuild the virtual
environment with a Python executable that has sqlite support. If Snakemake fails
with `AttributeError: module 'pulp' has no attribute 'list_solvers'`, reinstall
the pinned PuLP version:

```bash
python -m pip install --force-reinstall "pulp==2.7.0"
python -m pip install -r requirements-dev.txt
python -m snakemake --version
```

## Install Other Dependencies

Create a folder to store third-party software and the runtime path file:

```bash
mkdir -p deps
touch deps/soft_paths.txt
```

Add paths for the required software in `deps/soft_paths.txt`. The file is
tab-delimited:

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

Replace every `/path/to/...` entry with the real executable path on the server.

After editing `deps/soft_paths.txt`, validate both third-party executable paths
and Python packages before running the workflow:

```bash
python analysis_scripts/check_runtime_dependencies.py \
  deps/soft_paths.txt \
  --requirements requirements-dev.txt
```

The same check is available as a Snakemake target:

```bash
python -m snakemake --configfile workflow/config/w3_5_2_linux.yaml \
  --config soft_paths=deps/soft_paths.txt \
  --cores 1 check_runtime_dependencies
```

The repository test suite also protects the Python dependency list. The test
`tests/test_source_checkout.py::test_legacy_python_packages_are_in_runtime_requirements`
checks that the legacy manual-install package set remains listed in
`requirements-dev.txt`: `biopython`, `pysam`, `pandas`, `numpy`, `openpyxl`,
`xlsxwriter`, `matplotlib`, `polars`, and `fastexcel`.

```bash
python -m pytest tests/test_source_checkout.py
```

### Install minimap2

```bash
cd /path/to/hifisr/deps
curl -L https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 | tar -jxvf -
```

### Install samtools

```bash
cd /path/to/hifisr/deps
wget -c https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2
tar -xjf samtools-1.21.tar.bz2
cd samtools-1.21
autoheader
autoconf -Wno-syntax
./configure --prefix=/path/to/hifisr/deps/samtools/samtools-1.21
make -j 20
make install
cd ..
rm -rf samtools-1.21 samtools-1.21.tar.bz2
```

### Install seqkit

```bash
cd /path/to/hifisr/deps
wget -c http://app.shenwei.me/data/seqkit/seqkit_linux_amd64.tar.gz
tar -zxf seqkit_linux_amd64.tar.gz
rm seqkit_linux_amd64.tar.gz
```

### Install mecat

```bash
cd /path/to/hifisr/deps
git clone https://github.com/xiaochuanle/MECAT2.git
cd MECAT2
make -j 20
```

### Install blastn

```bash
cd /path/to/hifisr/deps
wget -c https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz
tar -zxf ncbi-blast-2.16.0+-x64-linux.tar.gz
rm ncbi-blast-2.16.0+-x64-linux.tar.gz
```

### Install bcftools

```bash
cd /path/to/hifisr/deps
wget -c https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2
tar -xjf bcftools-1.21.tar.bz2
cd bcftools-1.21
./configure --prefix=/path/to/hifisr/deps/bcftools/bcftools-1.21
make -j 20
make install
cd ..
rm -rf bcftools-1.21 bcftools-1.21.tar.bz2
```

### Install bamtools

```bash
cd /path/to/hifisr/deps
wget -c https://github.com/pezmaster31/bamtools/archive/refs/tags/v2.5.2.zip
unzip v2.5.2.zip
rm v2.5.2.zip
cd bamtools-2.5.2
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/hifisr/deps/bamtools/bamtools-2.5.2 ..
make -j 20
make install
cd ../..
rm -rf bamtools-2.5.2
```

### Install pigz

```bash
cd /path/to/hifisr/deps
wget -c https://zlib.net/pigz/pigz.tar.gz
tar -zxf pigz.tar.gz
rm pigz.tar.gz
cd pigz
make -j 20
```

### Install Flye

```bash
cd /path/to/hifisr/deps
git clone https://github.com/fenderglass/Flye
cd Flye
python -m pip install .
```

### Install Bandage

```bash
cd /path/to/hifisr/deps
wget -c https://github.com/rrwick/Bandage/releases/download/v0.8.1/Bandage_Ubuntu_static_v0_8_1.zip
unzip Bandage_Ubuntu_static_v0_8_1.zip
rm Bandage_Ubuntu_static_v0_8_1.zip sample_LastGraph
```

### Install hifiasm

```bash
cd /path/to/hifisr/deps
git clone https://github.com/chhylp123/hifiasm
cd hifiasm
make -j 8
```

### Install canu

```bash
cd /path/to/hifisr/deps
wget -c https://github.com/marbl/canu/releases/download/v2.3/canu-2.3.Linux-amd64.tar.xz
tar -xJf canu-2.3.Linux-amd64.tar.xz
rm canu-2.3.Linux-amd64.tar.xz
```

## Source Checkout Layout

The Snakemake workflow runs directly from the source tree and does not require
installing the `hifisr` package into the environment. Keep this layout intact:

```text
hifisr/
  analysis_scripts/
    _bootstrap.py
  hifisr_functions/
  workflow/
  Snakefile
```
