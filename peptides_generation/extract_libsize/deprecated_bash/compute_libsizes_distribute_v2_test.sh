#!/bin/bash
set -e

gene_expre_file='/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_librarysize_pya.0.17.1_conf2_annot_ref_chrall_cap1000/cohort_mutNone/*tmp*batch*/gene_expression_detail.pq'
#gene_expre_file='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_8f09893_conf2_annotFrame_cap1000_runs_pya0.17.1/TCGA_All_Normals/cohort_mutNone/tmp_out_ref_batch_9*/gene_expression_detail.pq'
#result_file=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_libsizes_conf2_annotFrame_cap1000_runs_pya0.17.1/TCGA_Ovarian_374/cohort_mutNone/mini_result
result_file=/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_librarysize_pya.0.17.1_conf2_annot_ref_chrall_cap1000/cohort_mutNone/GTEX_hg38_library_size_75perce_v2_test.tsv
sample_file=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/GTEX/GTEx_normal_samples_12-2020_nohead.tsv

echo $gene_expre_file
echo $result_file

while read sample; do
	sample=$(echo $sample | tr -d '.')
       	echo "$(date +%x_%r) ${sample} processed"	
	parquet-tools csv ${gene_expre_file} -c ${sample} > $result_file
	echo "$(date +%x_%r) column extracted with parquet tools"
	tail -n +2 $result_file > ${result_file}_2
	echo "$(date +%x_%r) tail done"
       	sort  -n ${result_file}_2 > ${result_file}_3
	echo "$(date +%x_%r) sort done"
	awk 'BEGIN{c=0} {total[c]=$1; c++;} END{print total[int(NR*0.75)]}' ${result_file}_3  >> ${result_file}_4
	echo "$(date +%x_%r) awk done"
done < ${sample_file}
