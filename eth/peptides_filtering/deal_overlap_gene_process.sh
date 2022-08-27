### Find genes commun to gtex and TCGA Breast ###


folder_gtex=/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_v3_TEST_merged3_372a147_medium_run_pya.0.17.1_conf2_annot_ref_chrall_cap/cohort_mutNone
folder_breast=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/commit_v3_TEST_merged3_372a147_medium_run_conf2_annotFrame_cap0_runs/TCGA_Breast_1102/cohort_mutNone
outdir=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/commit_v3_TEST_merged3_372a147_medium_run_conf2_annotFrame_cap0_runs/selected_genes_rerun1
mkdir -p ${outdir}
for folder in ${folder_breast}/tmp*; do
	batch=$(basename $folder)
	if [ -f ${folder_breast}/${batch}/ref_graph_kmer_SegmExpr.pq ] && [ -f ${folder_gtex}/${batch}/ref_graph_kmer_SegmExpr.pq ]; then  
      		echo ${batch} >> ${outdir}/genes_commun.txt; echo ${folder_breast}/${batch}/ref_graph_kmer_JuncExpr.pq >> ${outdir}/select_breast_junct_Expr.txt 
		echo ${folder_breast}/${batch}/ref_graph_kmer_SegmExpr.pq >> ${outdir}/select_breast_segm_Expr.txt 
	       	echo ${folder_gtex}/${batch}/ref_graph_kmer_JuncExpr.pq >> ${outdir}/select_gtex_realign_junct_Expr.txt 
	       	echo ${folder_gtex}/${batch}/ref_graph_kmer_SegmExpr.pq >> ${outdir}/select_gtex_realign_segm_Expr.txt 
       	fi
done
