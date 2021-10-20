#!/bin/bash
set -e

suffix_out='Normals'
if [ $suffix_out == 'GTEX2019' ]; then 
	outfolder=/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_v2.5f0752a_pya.0.17.1_conf2_annot_ref_chrall_cap/cohort_mutNone
elif  [ $suffix_out == 'BRCA' ]; then 
	outfolder=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102/cohort_mutNone
elif [ $suffix_out == 'Normals' ]; then
	outfolder=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_All_Normals/cohort_mutNone
fi 

gene_interest_file=./dev_genes
#gene_interest_file=./tmp_genes


run_status='collect_finished' #'collect_finished'

if [ $run_status == 'collect_finished' ] ; then  
	echo "collect processed genes to $(dirname ${outfolder})/ref_gene_folder_eq"
	for folder in ${outfolder}/tmp*/gene_expression_detail.pq ; 
		do echo $(dirname $folder)  $(parquet-tools csv -c gene $folder  | cut -f1 -d ',' | grep ENS) > $(dirname ${outfolder})/ref_gene_folder_eq	
	done
	#TODO add check for folders without gene output 
	else
		file_out=Segm
		if [ ${file_out} == 'Junc' ]; then
		        file_out_full=ref_graph_kmer_JuncExpr.pq
		elif  [ ${file_out} == 'Segm' ]; then
       	 		file_out_full=ref_graph_kmer_SegmExpr.pq
		elif [ ${file_out} == 'annot' ] ; then
        		file_out_full=ref_annot_9mer.pq
		fi
	while read gene; do
		myfile=$(grep -E "${gene}" $(dirname ${outfolder})/ref_gene_folder_eq | sed -z 's/\n/|/g;s/|$/\n/' | cut -f1 -d ' ' | sed "s,$,/${file_out_full},g")
	       	if [ -f "${myfile}" ]; then
			echo $myfile >> ./tmp20211015_input_${file_out}_${suffix_out}
		fi
 	done < ${gene_interest_file} 

fi 
