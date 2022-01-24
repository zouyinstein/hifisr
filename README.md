# hifisr

HiFi-SR is a Python-based pipeline for the detection of plant mitochondrial structural rearrangements based on the mapping of HiFi reads, or PacBio Circular Consensus Sequencing (ccs) data, to a reference genome (e.g., hypothetical master cycle DNA).

**Deps**

The pipeline has been tested in WSL2 distribution Ubuntu-20.04. It should work in other Linux operating system, such as CentOS.

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
conda install samtools bamtools blast seqkit -c bioconda
ln -sf /home/zouyi/.conda/envs/hifisr/lib/libcrypto.so.1.1 /home/zouyi/.conda/envs/hifisr/lib/libcrypto.so.1.0.0  # create a soft link
```

* Install minimap2

```bash
cd hifisr/deps
curl -L https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 | tar -jxvf -
./minimap2-2.24_x64-linux/minimap2
export PATH=/mnt/e/02.vol1/03.work/03.HiFi/hifisr/deps/minimap2-2.24_x64-linux:$PATH
```

* Install required Python packages

```bash
pip install biopython
```

**Example**

Analyze of an example wild-type *Arabidopsis thaliana* Col-0 dataset (ERR6210723, 14.6 Gb, Naish et al., 2021, Science): 

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
bash ../../run_hifisr.sh CEN 16 fastq &

```
