#!/bin/bash
set -e

gene_expre_file=$1 #/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_libsizes_conf2_annotFrame_cap1000_runs_pya0.17.1/TCGA_Ovarian_374/cohort_mutNone/gene_expression_detail_mini.pq
result_file=$2 #/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_libsizes_conf2_annotFrame_cap1000_runs_pya0.17.1/TCGA_Ovarian_374/cohort_mutNone/mini_result
field_num=$(parquet-tools csv ${gene_expre_file}| head -1 |awk -F "," '{print NF; exit}')
for sample in $(seq 2 ${field_num}); do
       	echo "${sample}/${field_num} processed"	
	echo "$(parquet-tools csv ${gene_expre_file}| head -1 |  cut -f${sample} -d ',') $(parquet-tools csv ${gene_expre_file}| tail -n +2 | cut -f${sample} -d ',' | sort  -n | awk 'BEGIN{c=0} {total[c]=$1; c++;} END{print total[int(NR*0.75)]}')" >> ${result_file}
done
