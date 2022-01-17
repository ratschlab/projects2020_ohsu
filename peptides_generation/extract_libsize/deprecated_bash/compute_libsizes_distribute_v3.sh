#!/bin/bash
set -e

gene_expre_file='/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_librarysize_pya.0.17.1_conf2_annot_ref_chrall_cap1000/cohort_mutNone/*tmp*batch*/gene_expression_detail.pq'
#gene_expre_file='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_8f09893_conf2_annotFrame_cap1000_runs_pya0.17.1/TCGA_All_Normals/cohort_mutNone/tmp_out_ref_batch_9*/gene_expression_detail.pq'
#result_file=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_libsizes_conf2_annotFrame_cap1000_runs_pya0.17.1/TCGA_Ovarian_374/cohort_mutNone/mini_result
result_file=/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_librarysize_pya.0.17.1_conf2_annot_ref_chrall_cap1000/cohort_mutNone/GTEX_hg38_genes_expr_v3.tsv
sample_file=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/GTEX/GTEx_normal_samples_12-2020_nohead.tsv

echo $gene_expre_file
echo $result_file

echo "${sample} processed"	
time parquet-tools csv ${gene_expre_file} > ${result_file}
