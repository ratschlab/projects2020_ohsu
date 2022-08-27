
run_cohort='gtex'
expression_type='JUNC' #SEGM


if [ ${run_cohort} == 'gtex' ]; 
	dir_links=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/commit_v3_TEST_merged3_372a147_medium_run_conf2_annotFrame_cap0_runs/background_relink
	if [ ${expression_type} == 'JUNC' ] ; then
                file_to_link=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/commit_v3_TEST_merged3_372a147_medium_run_conf2_annotFrame_cap0_runs/selected_genes_rerun1/select_gtex_realign_junct_Expr.txt
        elif  [ ${expression_type} == 'SEGM' ] ; then
                file_to_link=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/commit_v3_TEST_merged3_372a147_medium_run_conf2_annotFrame_cap0_runs/selected_genes_rerun1/select_gtex_realign_segm_Expr.txt

else
	dir_links=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/commit_v3_TEST_merged3_372a147_medium_run_conf2_annotFrame_cap0_runs/TCGA_Breast_1102/cohort_mutNone_relink
	if [ ${expression_type} == 'JUNC' ] ; then 
		file_to_link=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/commit_v3_TEST_merged3_372a147_medium_run_conf2_annotFrame_cap0_runs/selected_genes_rerun1/select_breast_junct_Expr.txt
	elif  [ ${expression_type} == 'SEGM' ] ; then
		file_to_link=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/commit_v3_TEST_merged3_372a147_medium_run_conf2_annotFrame_cap0_runs/selected_genes_rerun1/select_breast_segm_Expr.txt


while read brca; do ln -s $brca $(basename $(dirname $brca) | cut -f5 -d '_')_$(basename $brca)  ;done < ${file_to_link}



