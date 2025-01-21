# hifisr

HiFi-SR is a Python-based pipeline for the detection of plant mitochondrial structural rearrangements based on the mapping of PacBio high-fidelity (HiFi) reads or Circular Consensus Sequencing (ccs) data, to a reference genome (i.e., the hypothetical master cycle DNA).

**Deps**

The pipeline has been tested in WSL2 distribution Ubuntu-20.04. It shall work in other Linux operating system, such as CentOS.

* Download the hifisr repository

```
git clone https://github.com/zouyinstein/hifisr
```

* Create and activate a conda environment.

```bash
conda create -n hifisr python=3.13
conda activate hifisr
```

* Use Anaconda3 to install required packages.

```bash
conda install pigz -c conda-forge
conda install samtools bamtools blast seqkit parafly -c bioconda
# create a soft link to ensure samtools can work
ln -sf ${HOME}/.conda/envs/hifisr/lib/libcrypto.so.1.1 ${HOME}/.conda/envs/hifisr/lib/libcrypto.so.1.0.0  
```

* Install bcftools
  demo
* Install minimap2

```bash
cd hifisr/deps
curl -L https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 | tar -jxvf -
./minimap2-2.24_x64-linux/minimap2
```

* Install MECAT2

```bash
cd hifisr/deps
git clone https://github.com/xiaochuanle/MECAT2.git
cd MECAT2
make -j

export PATH="$PWD/deps/MECAT2/Linux-amd64/bin:$PATH"
```

* Install metaFlye

```bash
cd hifisr/deps
git clone https://github.com/fenderglass/Flye
cd Flye
python setup.py install
```

* Install required Python packages

```bash
pip install biopython pandas openpyxl
```

**Test**

```bash
# Make sure the HiFi-SR repository has been downloaded (git clone https://github.com/zouyinstein/hifisr).
# Make sure the hifisr environment has been activated (conda activate hifisr).
# Make sure the dependent third-party softwares/packages has been installed.
# Change working directory to hifisr

```

**Description of results**

demo

**Example 1**

Analyze of an example wild-type *Arabidopsis thaliana* dataset Col-CEN ([ERR6210723](https://www.ncbi.nlm.nih.gov/sra/ERR6210723), 14.6 Gb, [Naish et al., 2021, Science](https://www.science.org/doi/10.1126/science.abi7489)):
