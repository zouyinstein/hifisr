samples_txt_prefix=$1

mkdir -p ${samples_txt_prefix}/{mito,plastid}
cat ${samples_txt_prefix}.txt | while read sample; do cat ${sample}/mito/mito_reads_type.list >> mito_reads_type_${samples_txt_prefix}.txt; done
cut -f3 mito_reads_type_${samples_txt_prefix}.txt | sort -u > ${samples_txt_prefix}/mito/mito_reads_type.list
rm mito_reads_type_${samples_txt_prefix}.txt
cat ${samples_txt_prefix}/mito/mito_reads_type.list | while read sub_type; do mkdir -p ${samples_txt_prefix}/mito/${sub_type}; done
cat ${samples_txt_prefix}.txt | while read sample; do cat ${sample}/plastid/plastid_reads_type.list >> plastid_reads_type_${samples_txt_prefix}.txt; done
cut -f3 plastid_reads_type_${samples_txt_prefix}.txt | sort -u > ${samples_txt_prefix}/plastid/plastid_reads_type.list
rm plastid_reads_type_${samples_txt_prefix}.txt
cat ${samples_txt_prefix}/plastid/plastid_reads_type.list | while read sub_type; do mkdir -p ${samples_txt_prefix}/plastid/${sub_type}; done

cat ${samples_txt_prefix}.txt | while read sample; do cp ${sample}/mito/type_1_Ref/type_1_Ref_${sample}_mito_ids.list ${samples_txt_prefix}/mito/type_1_Ref/type_1_Ref_${sample}_mito_ids.list; done
cat ${samples_txt_prefix}.txt | while read sample; do cp ${sample}/plastid/type_1_Ref/type_1_Ref_${sample}_plastid_ids.list ${samples_txt_prefix}/plastid/type_1_Ref/type_1_Ref_${sample}_plastid_ids.list; done

cat ${samples_txt_prefix}/mito/mito_reads_type.list | grep "type_2" | while read sub_type; do cat ${samples_txt_prefix}.txt | while read sample; do cat ${sample}/mito/${sub_type}/count_se1_ss2_${sub_type}_${sample}_mito.txt >> ${samples_txt_prefix}/mito/${sub_type}/count_se1_ss2_${sub_type}_${samples_txt_prefix}_mito.txt; done; done
cat ${samples_txt_prefix}/plastid/plastid_reads_type.list | grep "type_2" | while read sub_type; do cat ${samples_txt_prefix}.txt | while read sample; do cat ${sample}/plastid/${sub_type}/count_se1_ss2_${sub_type}_${sample}_plastid.txt >> ${samples_txt_prefix}/plastid/${sub_type}/count_se1_ss2_${sub_type}_${samples_txt_prefix}_plastid.txt; done; done
cat ${samples_txt_prefix}.txt | while read sample; do cat ${samples_txt_prefix}/mito/mito_reads_type.list | grep "type_2" | while read sub_type; do cd ${samples_txt_prefix}/mito/${sub_type}; python ../../../../scripts/merge_type_2.py ${sub_type} ${samples_txt_prefix} mito; cd ../../..; done; done
cat ${samples_txt_prefix}.txt | while read sample; do cat ${samples_txt_prefix}/plastid/plastid_reads_type.list | grep "type_2" | while read sub_type; do cd ${samples_txt_prefix}/plastid/${sub_type}; python ../../../../scripts/merge_type_2.py ${sub_type} ${samples_txt_prefix} plastid; cd ../../..; done; done

cat ${samples_txt_prefix}/mito/mito_reads_type.list | grep "type_3" | while read sub_type; do cat ${samples_txt_prefix}.txt | while read sample; do cat ${sample}/mito/${sub_type}/count_se1_ss2_se2_ss3_${sub_type}_${sample}_mito.txt >> ${samples_txt_prefix}/mito/${sub_type}/count_se1_ss2_se2_ss3_${sub_type}_${samples_txt_prefix}_mito.txt; done; done
cat ${samples_txt_prefix}/plastid/plastid_reads_type.list | grep "type_3" | while read sub_type; do cat ${samples_txt_prefix}.txt | while read sample; do cat ${sample}/plastid/${sub_type}/count_se1_ss2_se2_ss3_${sub_type}_${sample}_plastid.txt >> ${samples_txt_prefix}/plastid/${sub_type}/count_se1_ss2_se2_ss3_${sub_type}_${samples_txt_prefix}_plastid.txt; done; done
cat ${samples_txt_prefix}.txt | while read sample; do cat ${samples_txt_prefix}/mito/mito_reads_type.list | grep "type_3" | while read sub_type; do cd ${samples_txt_prefix}/mito/${sub_type}; python ../../../../scripts/merge_type_3.py ${sub_type} ${samples_txt_prefix} mito; cd ../../..; done; done
cat ${samples_txt_prefix}.txt | while read sample; do cat ${samples_txt_prefix}/plastid/plastid_reads_type.list | grep "type_3" | while read sub_type; do cd ${samples_txt_prefix}/plastid/${sub_type}; python ../../../../scripts/merge_type_3.py ${sub_type} ${samples_txt_prefix} plastid; cd ../../..; done; done