sample=$1      # sample name
thread=$2      # number of threads
input_type=$3  # fastq or fasta
type_number=5  # analyze up to 4 rearragements in a single read

# map and separate mitochondrial- and plastid-genome-mapped reads
minimap2 -t ${thread} -ax map-hifi ${sample}_ref.fa ${sample}.${input_type} > ${sample}.sam
samtools view -Sb -@ ${thread} ${sample}.sam > ${sample}.bam
samtools sort ${sample}.bam -@ ${thread} > ${sample}.sorted.bam
samtools index -@ ${thread} ${sample}.sorted.bam
samtools flagstat ${sample}.sorted.bam -@ ${thread} > ${sample}.stat
bamtools split -in ${sample}.sorted.bam -reference
find ./ -name "*.REF_*.bam" > split_bam.list
cat split_bam.list | cut -c3- | while read i; do samtools fastq $i -@ ${thread} -0 ${i/.bam/}.fastq; done
find ./ -name "*.REF_*.fastq" > split_fastq.list_all
cat split_fastq.list_all | cut -c3- | while read i; do seqkit fq2fa $i -j ${thread} -o ${i/.fastq/}.fasta; seqkit stat -T -a -j ${thread} $i > $i.stat; done