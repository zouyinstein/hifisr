sample=$1
thread=$2
input_type=$3
type_number=5

export PATH=$PATH:/usr/software/minimap2-2.20_x64-linux/
export PATH=$PATH:/usr/software/anaconda3/bin/
export PATH=/usr/software/Python-3.9.9/bin/:$PATH

minimap2 -t ${thread} -ax map-hifi ref.fa ${sample}.${input_type} > ${sample}.sam
samtools view -Sb -@ ${thread} ${sample}.sam > ${sample}.bam
samtools sort ${sample}.bam -@ ${thread} > ${sample}.sorted.bam
samtools index -@ ${thread} ${sample}.sorted.bam
samtools flagstat ${sample}.sorted.bam -@ ${thread} > ${sample}.stat
bamtools split -in ${sample}.sorted.bam -reference
find ./ -name "*.REF_*.bam" > split_bam.list
cat split_bam.list | cut -c3- | while read i; do samtools fastq $i -@ ${thread} -0 ${i/.bam/}.fastq; done
find ./ -name "*.REF_*.fastq" > split_fastq.list_all
cat split_fastq.list_all | cut -c3- | while read i; do seqkit fq2fa $i -j ${thread} -o ${i/.fastq/}.fasta; seqkit stat -T -a -j ${thread} $i > $i.stat; done
makeblastdb -dbtype nucl -parse_seqids -in ${sample}_mito.fa
makeblastdb -dbtype nucl -parse_seqids -in ${sample}_plastid.fa
blastn -query ${sample}.sorted.REF_mito.fasta -db ${sample}_mito.fa -num_threads ${thread} -outfmt 5 -out blastn_${sample}_mito.xml
blastn -query ${sample}.sorted.REF_plastid.fasta -db ${sample}_plastid.fa -num_threads ${thread} -outfmt 5 -out blastn_${sample}_plastid.xml
python ../../scripts/analysis_blastn.py blastn_${sample}_mito.xml > blastn_result_${sample}_mito.txt
python ../../scripts/analysis_blastn.py blastn_${sample}_plastid.xml > blastn_result_${sample}_plastid.txt
python ../../scripts/get_type.py blastn_result_${sample}_mito.txt > blastn_type_count_result_${sample}_mito.txt
seq 1 ${type_number} | while read i; do grep "type_${i}" blastn_type_count_result_${sample}_mito.txt > blastn_type_${i}_result_${sample}_mito.txt; python ../../scripts/get_type_${i}.py blastn_type_${i}_result_${sample}_mito.txt > type_${i}_all_${sample}_mito.txt; done
python scripts/get_type.py blastn_result_${sample}_plastid.txt > blastn_type_count_result_${sample}_plastid.txt
seq 1 ${type_number} | while read i; do grep "type_${i}" blastn_type_count_result_${sample}_plastid.txt > blastn_type_${i}_result_${sample}_plastid.txt; python ../../scripts/get_type_${i}.py blastn_type_${i}_result_${sample}_plastid.txt > type_${i}_all_${sample}_plastid.txt; done
pigz -p ${thread} ${sample}.${input_type}
rm *sam *bam *bai split*
mkdir -p {mito,plastid,reads,pre}
mv type_*_all_${sample}_mito.txt mito
mv type_*_all_${sample}_plastid.txt plastid
mv ${sample}.${input_type}.gz reads
mv ${sample}* reads
mv blastn* pre
cd mito
cut -f3 type_*_all_${sample}_mito.txt | sort -u > mito_reads_type.list
cat mito_reads_type.list | while read i; do mkdir "$i"; grep "$i" type_*_all_${sample}_mito.txt > ${i}/${i}_${sample}_mito.txt; wc -l ${i}/${i}_${sample}_mito.txt; done
cd ../plastid
cut -f3 type_*_all_${sample}_plastid.txt | sort -u > plastid_reads_type.list
cat plastid_reads_type.list | while read i; do mkdir "$i"; grep "$i" type_*_all_${sample}_plastid.txt > ${i}/${i}_${sample}_plastid.txt; wc -l ${i}/${i}_${sample}_plastid.txt; done
cd ..
cd mito
cat mito_reads_type.list | grep "type_2" | while read sub_type; do cd ${sub_type}; python ../../../../scripts/analysis_type_2.py ${sub_type}_${sample}_mito.txt ${sample}; cd ..; done
cd ../plastid
cat plastid_reads_type.list | grep "type_2" | while read sub_type; do cd ${sub_type}; python ../../../../scripts/analysis_type_2.py ${sub_type}_${sample}_plastid.txt ${sample}; cd ..; done
cd ..