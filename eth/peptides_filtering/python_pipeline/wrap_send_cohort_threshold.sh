#!/bin/bash

start_gene=$1
end_gene=$2 #60609
hours=120
cpus=$3
mem=60G
run_cohort=gtex

if [ ${run_cohort} == 'gtex' ]; then
	# ---- Cohort dependant Submission Parameters ----
	job_name=back_${start_gene}_${end_gene}
	suffix=py_filter_gtex_large_mem_nooverlap_completionmissing
	launch_script=send_background_threshold.sh
	log_dir=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/lsf
	# ---- Cohort dependant Run Parameters ----
	path_cohort='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/cohort_mutNone'
	whitelist_tag='10-21overlap_order'
	whitelist='/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/GTEX/GTEx_sample_IDs_10-2021_lib_graph_juliannelist'
	path_libsize='/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/ARCHIV_keep_runs/GTEX2019_commit_v3_TEST_merged3_372a147_medium_run_pya.0.17.1_conf2_annot_ref_chrall_cap/expression_counts.libsize.tsv' 
	normalizer_libsize=400000
	sample_pattern='SRR'
elif  [ ${run_cohort} == 'brca' ]; then
	# ---- Cohort dependant Submission Parameters ----
	job_name=brca_${start_gene}_${end_gene}
	suffix=py_filter_brca_large
	launch_script=send_foreground_threshold_b.sh
	log_dir=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/lsf
	# ---- Cohort dependant Run Parameters ----
	path_cohort='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/cohort_mutNone'
	whitelist_tag=''
	whitelist='/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/TCGA_foreground/sample_full_BRCA_1102_format.tsv'
	path_libsize='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/expression_counts.libsize.tsv'
	normalizer_libsize=400000
	sample_pattern='TCGA'
elif  [ ${run_cohort} == 'ov' ]; then
        # ---- Cohort dependant Submission Parameters ----
	job_name=ov_${start_gene}_${end_gene}
	suffix=py_filter_ov_large
	launch_script=send_foreground_threshold_o.sh
	log_dir=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374/lsf
	# ---- Cohort dependant Run Parameters ----
	path_cohort='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374/cohort_mutNone'
	whitelist_tag='full'
	whitelist='/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/TCGA_foreground/sample_full_Ov_378_format.tsv'
	path_libsize='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374/expression_counts.libsize.tsv'
	normalizer_libsize=400000
	sample_pattern='TCGA'

fi 


logging_base=${log_dir}/${suffix}_${cpus}p_${mem}_${start_gene}-${end_gene}
err_file=${logging_base}.err
log_file=${logging_base}.out


# ---- Launch run ----
echo "#!/bin/bash" > ${launch_script}
echo "#SBATCH -o ${log_file}" >> ${launch_script}
echo "#SBATCH -e ${err_file}" >> ${launch_script}
echo "#SBATCH -J ${job_name}" >> ${launch_script}
echo "#SBATCH --time=${hours}:00:00" >> ${launch_script}
echo "#SBATCH --cpus-per-task=${cpus}" >> ${launch_script}
echo "#SBATCH --mem=${mem}" >> ${launch_script}
cmd="python ./cohort_threshold.py --processes ${cpus} --start-id ${start_gene} --end-id ${end_gene} --path-cohort ${path_cohort} --whitelist-tag ${whitelist_tag} --whitelist ${whitelist} --path-libsize ${path_libsize} --normalizer-libsize ${normalizer_libsize} --filters 0.0 1.0 2.0 3.0 5.0 10.0 --sample-pattern ${sample_pattern} --do-normalize " #--do-overwrite"
echo $cmd >>  ${launch_script}

echo "Output to ${log_file}"
echo "Errors to ${err_file}"
sbatch ${launch_script}
squeue -u ${USER} | grep ${job_name}
