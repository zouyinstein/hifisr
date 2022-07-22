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

# realign mitochondrial- and plastid-genome-mapped reads
makeblastdb -dbtype nucl -parse_seqids -in ${sample}_mito.fa
makeblastdb -dbtype nucl -parse_seqids -in ${sample}_plastid.fa
blastn -query ${sample}.sorted.REF_mito.fasta -db ${sample}_mito.fa -num_threads ${thread} -outfmt 5 -out blastn_${sample}_mito.xml
blastn -query ${sample}.sorted.REF_plastid.fasta -db ${sample}_plastid.fa -num_threads ${thread} -outfmt 5 -out blastn_${sample}_plastid.xml

# process raw blastn results
python ../../scripts/analysis_blastn.py blastn_${sample}_mito.xml ${sample}_mito.fa > blastn_result_${sample}_mito.txt
python ../../scripts/analysis_blastn.py blastn_${sample}_plastid.xml ${sample}_plastid.fa > blastn_result_${sample}_plastid.txt
python ../../scripts/get_type.py blastn_result_${sample}_mito.txt > blastn_type_count_result_${sample}_mito.txt
seq 1 ${type_number} | while read i; do grep "type_${i}" blastn_type_count_result_${sample}_mito.txt > blastn_type_${i}_result_${sample}_mito.txt; python ../../scripts/get_type_${i}.py blastn_type_${i}_result_${sample}_mito.txt > type_${i}_all_${sample}_mito.txt; done
python ../../scripts/get_type.py blastn_result_${sample}_plastid.txt > blastn_type_count_result_${sample}_plastid.txt
seq 1 ${type_number} | while read i; do grep "type_${i}" blastn_type_count_result_${sample}_plastid.txt > blastn_type_${i}_result_${sample}_plastid.txt; python ../../scripts/get_type_${i}.py blastn_type_${i}_result_${sample}_plastid.txt > type_${i}_all_${sample}_plastid.txt; done

# tidy intermediate files
pigz -p ${thread} ${sample}.${input_type}
rm ./*sam ./*bam ./*bai ./split*
mkdir -p {mito,plastid,reads,blastn}
mv type_*_all_${sample}_mito.txt mito
mv type_*_all_${sample}_plastid.txt plastid
mv ${sample}.${input_type}.gz reads
mv ${sample}* reads
mv blastn_* blastn

# analyze subtypes
cd mito
cut -f3 type_*_all_${sample}_mito.txt | sort -u > mito_reads_type.list
cat mito_reads_type.list | while read i; do mkdir "$i"; grep "$i" type_*_all_${sample}_mito.txt > ${i}/${i}_${sample}_mito.txt; wc -l ${i}/${i}_${sample}_mito.txt; done
cd ../plastid
cut -f3 type_*_all_${sample}_plastid.txt | sort -u > plastid_reads_type.list
cat plastid_reads_type.list | while read i; do mkdir "$i"; grep "$i" type_*_all_${sample}_plastid.txt > ${i}/${i}_${sample}_plastid.txt; wc -l ${i}/${i}_${sample}_plastid.txt; done
cd ..

cut -f1 mito/type_1_Ref/type_1_Ref_${sample}_mito.txt | cut -d ":" -f2 > mito/type_1_Ref/type_1_Ref_${sample}_mito_ids.list
cut -f1 plastid/type_1_Ref/type_1_Ref_${sample}_plastid.txt | cut -d ":" -f2 > plastid/type_1_Ref/type_1_Ref_${sample}_plastid_ids.list

cd mito
cat mito_reads_type.list | grep "type_2" | while read sub_type; do cd ${sub_type}; python ../../../../scripts/analysis_type_2.py ${sub_type}_${sample}_mito.txt ${sample}; cd ..; done
cd ../plastid
cat plastid_reads_type.list | grep "type_2" | while read sub_type; do cd ${sub_type}; python ../../../../scripts/analysis_type_2.py ${sub_type}_${sample}_plastid.txt ${sample}; cd ..; done
cd ..

cd mito
cat mito_reads_type.list | grep "type_3" | while read sub_type; do cd ${sub_type}; python ../../../../scripts/analysis_type_3.py ${sub_type}_${sample}_mito.txt ${sample}; cd ..; done
cd ../plastid
cat plastid_reads_type.list | grep "type_3" | while read sub_type; do cd ${sub_type}; python ../../../../scripts/analysis_type_3.py ${sub_type}_${sample}_plastid.txt ${sample}; cd ..; done
cd ..
