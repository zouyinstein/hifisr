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
pip install biopython pandas openpyxl
```

**Test**

```bash
cd hifisr/test
cp -R start Col
cd Col
pigz -d -p 8 *.gz
# check conda activate hifisr
# check export PATH=/mnt/e/02.vol1/03.work/03.HiFi/hifisr/deps/minimap2-2.24_x64-linux:$PATH
python ../../hifisr.py -s Col single &
# OR python ../../hifisr.py -s Col -t 16 -i fastq single &
# clean results except for the start fold and its contents if you want to rerun the test
```


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
# check export PATH=/mnt/e/02.vol1/03.work/03.HiFi/hifisr/deps/minimap2-2.24_x64-linux:$PATH
python ../../hifisr.py -s CEN -t 10 -i fastq single &

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
# check export PATH=/mnt/e/02.vol1/03.work/03.HiFi/hifisr/deps/minimap2-2.24_x64-linux:$PATH
python ../../hifisr.py -s XJTU -t 10 -i fastq single &
```

**Merge reports of multiple samples**

Merge reports of Col-CEN and Col-XJTU:

```bash
cd hifisr/pre
echo CEN >> merge_1.txt
echo XJTU >> merge_1.txt
# check conda activate hifisr
# check export PATH=/mnt/e/02.vol1/03.work/03.HiFi/hifisr/deps/minimap2-2.24_x64-linux:$PATH
python ../../hifisr.py -m merge1 merge &
```

**Parameters**

```bash
./hifisr.py -h
```

usage: hifisr.py [-h] [-s str] [-t int] [-i str] [-n int] [-m str] mode

The Python user interface of the HiFi-SR pipeline

positional arguments:
  mode                  "single" to run hifisr for each sample; "merge" to merge results of multiple samples

optional arguments:
  -h, --help            show this help message and exit
  -s str, --sample_name str
                        In "single" mode, sample name, for example sample_1 (default: )
  -t int, --num_of_thread int
                        In "single" mode, number of threads to run the pipeline (default: 8)
  -i str, --input_type str
                        In "single" mode, input file type fastq or fasta (default: fastq)
  -n int, --type_number int
                        In "single" mode, analyze reads up to four rearrangements in a single read (default: 5)
  -m str, --samples_txt_prefix str
                        In "merge" mode, prefix of txt file containing multiple sample names, for example "merge_1" (default: )
