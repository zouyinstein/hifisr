export PATH="/usr/software/Filtlong/bin:$PATH"
export PATH="/usr/software/MECAT2/Linux-amd64/bin:$PATH"
sample=$1
genome=$2
thread=$3
length="1000"
total="200000000"
if [ $genome = "mito" ]; then genome_size=370; fi
if [ $genome = "plastid" ]; then genome_size=155; fi


filtlong  --min_length ${length} --target_bases ${total} reads/${sample}.sorted.REF_${genome}.fastq > ${sample}.${genome}_filt.fastq

mecat.pl config ${genome}_config.txt
sed s/PROJECT=/PROJECT=mecat_${genome}/ ${genome}_config.txt > new_config && mv new_config ${genome}_config.txt
sed s/RAWREADS=/RAWREADS=${sample}.${genome}_filt.fastq/ ${genome}_config.txt > new_config && mv new_config ${genome}_config.txt
sed s/GENOME_SIZE=/GENOME_SIZE=${genome_size}000/ ${genome}_config.txt > new_config && mv new_config ${genome}_config.txt
sed s/THREADS=/THREADS=${thread}/ ${genome}_config.txt > new_config && mv new_config ${genome}_config.txt
mecat.pl correct ${genome}_config.txt
cp mecat_${genome}/1-consensus/cns_final.fasta mecat_${genome}.fasta

flye --meta --pacbio-hifi mecat_${genome}.fasta --extra-params output_gfa_before_rr=1 --genome-size ${genome_size}K -t 70 -o flye_${genome}_${genome_size}K
cp flye_${genome}_${genome_size}K/20-repeat/graph_before_rr.gfa flye_${genome}_${genome_size}K_before_rr.gfa
cp flye_${genome}_${genome_size}K/assembly_graph.gfa flye_${genome}_${genome_size}K_after_rr.gfa

if [ -e "mecat_${genome}.fasta" ]; then rm -rf mecat_${genome}; fi
if [ -e "flye_${genome}_${genome_size}K_before_rr.gfa" -a -e "flye_${genome}_${genome_size}K_after_rr.gfa" ]; then rm -rf flye_${genome}_${genome_size}K; fi
mkdir draft_assembly
mv flye* draft_assembly
mv mecat* draft_assembly
mv *filt.fastq draft_assembly
mv *config.txt draft_assembly

