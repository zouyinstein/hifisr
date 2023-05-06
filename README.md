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
conda create -n hifisr python=3.9
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

* Install Filtlong

```bash
cd hifisr/deps
git clone https://github.com/rrwick/Filtlong.git
cd Filtlong
make -j

export PATH="$PWD/deps/Filtlong/bin:$PATH"
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
cd hifisr
# Add minimap2, filtlong, mecat.pl executables to the system PATH
export PATH="$PWD/deps/minimap2-2.24_x64-linux":$PATH
export PATH="$PWD/deps/Filtlong/bin:$PATH"
export PATH="$PWD/deps/MECAT2/Linux-amd64/bin:$PATH"

# Check and unzip the reference sequences
cd references
ls Col_mito.fa Col_plastid.fa Col_ref.fa.gz
pigz -d -p 8 Col_ref.fa.gz

# Check and unzip the input HiFi reads
cd ../data
pigz -d -p 8 Col.fastq.gz

# Change working directory to test and prepare the file input_files.txt
cd ../test
touch input_files.txt
# Contents of input_files.txt are tab-delimilated columns of sample name, input reads, total genome reference, mt genome reference, and pt  genome reference. The information of multiple samples can be added in different lines.

# Prepare the starting files and directories and a job script work.sh
python ../scripts/start_project.py input_files.txt 8
# Run the job script will start the HiFi-SR pipeline
nohup bash work.sh &
```

**Description of results**

demo

**Example 1**

Analyze of an example wild-type *Arabidopsis thaliana* dataset Col-CEN ([ERR6210723](https://www.ncbi.nlm.nih.gov/sra/ERR6210723), 14.6 Gb, [Naish et al., 2021, Science](https://www.science.org/doi/10.1126/science.abi7489)):

```bash
cd hifisr/pre
mkdir CEN  # CEN is the sample name
cd CEN
pigz -d -k -p 16 ERR6210723.fastq.gz
ln -sf ERR6210723.fastq CEN.fastq
# manually change the FASTA headers of mitochondrial and plastid genome refercences into mito and plastid for easily manipulation
cat Athaliana_447_TAIR10.id_Chr{1,2,3,4,5}.fa refs_cp28673mod.id_mito.fas refs_cp28673mod.id_plastid.fas > CEN_ref.fa
ln -sf refs_cp28673mod.id_mito.fas CEN_mito.fa
ln -sf refs_cp28673mod.id_plastid.fas CEN_plastid.fa
# check conda activate hifisr
# check export PATH=/path/to/hifisr/deps/minimap2-2.24_x64-linux:$PATH
python ../../hifisr.py -s CEN -t 10 -i fastq single > $(date +%s).log 2> $(date +%s).err &

```

**Example 2**

Analyze of an example wild-type *Arabidopsis thaliana* dataset Col-XJTU ([CRR302668](https://ngdc.cncb.ac.cn/gsa/browse/CRA004538/CRR302668), 22.9 Gb, [Wang et al., 2021, GPB](https://www.sciencedirect.com/science/article/pii/S1672022921001741)):

```bash
cd hifisr/pre
mkdir XJTU  # XJTU is the sample name
cd XJTU
pigz -d -k -p 16 CRR302668.fastq.gz
ln -sf CRR302668.fastq XJTU.fastq
# manually change the FASTA headers of mitochondrial and plastid genome refercences into mito and plastid for easily manipulation
cat Athaliana_447_TAIR10.id_Chr{1,2,3,4,5}.fa refs_cp28673mod.id_mito.fas refs_cp28673mod.id_plastid.fas > XJTU_ref.fa
ln -sf refs_cp28673mod.id_mito.fas XJTU_mito.fa
ln -sf refs_cp28673mod.id_plastid.fas XJTU_plastid.fa
# check conda activate hifisr
# check export PATH=/path/to/hifisr/deps/minimap2-2.24_x64-linux:$PATH
python ../../hifisr.py -s XJTU -t 10 -i fastq single > $(date +%s).log 2> $(date +%s).err &
```

**Merge reports of multiple samples**

Merge reports of Col-CEN and Col-XJTU:

```bash
cd hifisr/pre
echo CEN >> merge_1.txt
echo XJTU >> merge_1.txt
# check conda activate hifisr
# check export PATH=/path/to/hifisr/deps/minimap2-2.24_x64-linux:$PATH
python ../hifisr.py -m merge_1 merge > $(date +%s).log 2> $(date +%s).err &
```
