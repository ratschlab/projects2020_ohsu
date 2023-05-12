#!/bin/bash

# ----- SLURM Parameters -----
hours=24
cpus=1
mem=50G
tag='filter_param'
suffix="${tag}"
job_name='filtering'

launch_script=send_filtering_from_intermediate.py


# ----- ANALYSIS Parameters -----

sample_type='brca'

if [ "${sample_type}" == 'ov' ]; then
base_cancer='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374'
samples='TCGA-AO-A0JM-01A-21R-A056-07.all TCGA-C8-A12P-01A-11R-A115-07.all TCGA-BH-A18V-01A-11R-A12D-07.all TCGA-A2-A0D2-01A-21R-A034-07.all TCGA-A2-A0SX-01A-12R-A084-07.all'
intermediate_output="${base_cancer}/filtering_intermediate/complete_cancer_candidates_order_r.tsv.gz"
filtering_id='filters_22March_order_wany_wAnnot'
elif [ "${sample_type}" == 'brca' ]; then
base_cancer='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102'
samples='TCGA-25-1319-01A-01R-1565-13.all TCGA-25-1313-01A-01R-1565-13.all TCGA-61-2008-01A-02R-1568-13.all TCGA-24-1431-01A-01R-1566-13.all TCGA-24-2298-01A-01R-1569-13.all'
intermediate_output="${base_cancer}/filtering_intermediate/complete_cancer_candidates_order_r.tsv.gz"
filtering_id='filters_22March_order_wany_wAnnot_REPRODUCE'
fi

Threshold_target='0.0'
Threshold_cancer_cohort='None, 0.0, 2.0'
N_samples_cancer='None, 1, 5'
Threshold_normal_cohort='0.0, 1.0, 3.0, None'
N_samples_normal='1 2 10 None'

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
cmd="python ./star_intermediate_cohorts.py --basedir ${base_cancer} --intermediate-output ${intermediate_output} --filtering-id ${filtering_id} --target-samples ${samples} --Threshold-target ${Threshold_target} --Threshold-cancer-cohort $${Threshold_cancer_cohort} --N-samples-cancer ${N_samples_cancer} --Threshold-normal-cohort ${Threshold_normal_cohort} --N-samples-normal ${N_samples_normal}"
echo $cmd >>  ${launch_script}

echo "Output to ${log_file}"
echo "Errors to ${err_file}"
sbatch ${launch_script}
squeue -u ${USER} | grep ${job_name}
