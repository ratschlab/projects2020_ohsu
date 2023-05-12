#!/bin/bash

# ----- SLURM Parameters -----
hours=48
cpus=1
mem=150G
tag='order_r_complete_starversion'
suffix="${tag}"
job_name='star_big_mx'

launch_script=send_star_inter_cohorts.sh


# ----- ANALYSIS Parameters -----
sample_type='brca'


path_star='/cluster/work/grlab/projects/GTEx/rna_gencode32_realign/results'
whitelist='/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/GTEX/GTEx_sample_IDs_10-2021_lib_graph_juliannelist'
libsize='/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/ARCHIV_keep_runs/GTEX2019_commit_v3_TEST_merged3_372a147_medium_run_pya.0.17.1_conf2_annot_ref_chrall_cap/expression_counts.libsize.tsv'
normalizer=400000
#filter_thresholds=[0.0, 1.0, 2.0, 3.0, 5.0, 10.0]

if [ "${sample_type}" == 'ov' ]; then
base_cancer='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374'
big_matrix='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374/filtering_intermediate/complete_cancer_candidates_order_r_complete.tsv.gz'
jx_target_list='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374/filtering_samples/filters_22March_order_wany_wAnnot/tmp_all_experiments_jx.txt'
elif [ "${sample_type}" == 'brca' ]; then
base_cancer='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102'
big_matrix='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/filtering_intermediate/complete_cancer_candidates_order_r_complete.tsv.gz'
jx_target_list='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/filtering_samples/filters_22March_order_wany_wAnnot/tmp_all_experiments_jx.txt.gz'
fi

# ----- Logging -----
log_dir=${base_cancer}/lsf
logging_base=${log_dir}/${suffix}_${cpus}p_${mem}
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
cmd="python ./star_intermediate_cohorts.py --path-star ${path_star} --libsize ${libsize} --whitelist ${whitelist} --normalizer ${normalizer} --big-matrix ${big_matrix} --jx-target-list ${jx_target_list} --filter-thresholds 0.0 1.0 2.0 3.0 5.0 10.0"
echo $cmd >>  ${launch_script}

echo "Output to ${log_file}"
echo "Errors to ${err_file}"
sbatch ${launch_script}
squeue -u ${USER} | grep ${job_name}
