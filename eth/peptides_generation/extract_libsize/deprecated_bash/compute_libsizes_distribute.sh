#!/bin/bash
set -e

gene_expre_file='/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_librarysize_pya.0.17.1_conf2_annot_ref_chrall_cap1000/cohort_mutNone/*tmp*batch*/gene_expression_detail.pq'
#gene_expre_file='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_8f09893_conf2_annotFrame_cap1000_runs_pya0.17.1/TCGA_All_Normals/cohort_mutNone/tmp_out_ref_batch_9*/gene_expression_detail.pq'
#result_file=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_libsizes_conf2_annotFrame_cap1000_runs_pya0.17.1/TCGA_Ovarian_374/cohort_mutNone/mini_result
result_file=/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_librarysize_pya.0.17.1_conf2_annot_ref_chrall_cap1000/cohort_mutNone/GTEX_hg38_library_size_75perce.tsv
echo $gene_expre_file
echo $result_file
field_num=$(parquet-tools csv ${gene_expre_file}| head -1 |awk -F "," '{print NF; exit}')
for sample in $(seq 2 ${field_num}); do
       	echo "${sample}/${field_num} processed"	
	echo "$(parquet-tools csv ${gene_expre_file}| head -1 |  cut -f${sample} -d ',') $(parquet-tools csv ${gene_expre_file}| tail -n +2 | cut -f${sample} -d ',' | sort  -n | awk 'BEGIN{c=0} {total[c]=$1; c++;} END{print total[int(NR*0.75)]}')" >> ${result_file}
done
