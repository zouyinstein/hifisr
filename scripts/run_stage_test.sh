sample=$1      # sample name
thread=$2      # number of threads
input_type=$3  # fastq or fasta
type_number=5  # analyze up to 4 rearragements in a single read

# analyze subtypes
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