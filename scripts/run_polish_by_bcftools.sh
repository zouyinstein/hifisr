
export PATH="/usr/software/bcftools/bcftools-1.17/bin:$PATH"
# export PATH="/usr/software/hifisr/scripts:$PATH"
scripts_path="/usr/software/hifisr/scripts"

reference=$1  #AH-7
reads=$2   # mito
thread=$3    # 20

seqkit fq2fa ${reads} -j ${thread} -o ${reads/.fastq/}.fasta
makeblastdb -dbtype nucl -parse_seqids -in ${reference}
blastn -query ${reads/.fastq/}.fasta -db ${reference} -num_threads ${thread} -outfmt 5 -out blastn_${reads/.fastq/}.xml
python ${scripts_path}/analysis_blastn.py blastn_${reads/.fastq/}.xml ${reference} > blastn_${reads/.fastq/}.txt
python ${scripts_path}/get_type.py blastn_${reads/.fastq/}.txt > blastn_type_count_${reads/.fastq/}.txt
grep "type_1" blastn_type_count_${reads/.fastq/}.txt > blastn_type_1_${reads/.fastq/}.txt
python ${scripts_path}/get_type_1.py blastn_type_1_${reads/.fastq/}.txt > type_1_all_${reads/.fastq/}.txt
## remove type_1 reads that are not in full length
python ${scripts_path}/get_full_length_type_1.py type_1_all_${reads/.fastq/}.txt > type_1_full_${reads/.fastq/}_ids.txt
seqkit grep -f type_1_full_${reads/.fastq/}_ids.txt ${reference} -j ${thread} > type_1_full_${reads/.fastq/}.fasta
minimap2 -t ${thread} -ax map-hifi ${reference} type_1_full_${reads/.fastq/}.fasta > type_1_full_${reads/.fastq/}.sam
samtools view -Sb -@ ${thread} type_1_full_${reads/.fastq/}.sam > type_1_full_${reads/.fastq/}.bam
samtools sort type_1_full_${reads/.fastq/}.bam -@ ${thread} > type_1_full_${reads/.fastq/}.sorted.bam
samtools index -@ ${thread} type_1_full_${reads/.fastq/}.sorted.bam

bcftools mpileup -X pacbio-ccs --threads ${thread} -Ou -f ${reference} type_1_full_${reads/.fastq/}.sorted.bam | bcftools call --threads ${thread} -mv -Ou -o type_1_full_${reads/.fastq/}_bcftools.vcf
grep -v "^#" type_1_full_${reads/.fastq/}_bcftools.vcf | cut -f2,4,5,10 > type_1_full_${reads/.fastq/}_bcftools_small.vcf
python ${scripts_path}/check_snp_indel.py type_1_full_${reads/.fastq/}_bcftools_small.vcf type_1_full_${reads/.fastq/}.sorted.bam ${reference} > ${reads/.fastq/}_bcftools_checked.txt
python ${scripts_path}/replace_by_bcftools.py ${reads/.fastq/}_bcftools_checked.txt ${reference}
