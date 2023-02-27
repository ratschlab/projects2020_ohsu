#!/bin/bash

hours=24
cpus=30
mem=20G
suffix='pool_interm_filtering_79'
job_name='bix_matrix'

launch_script=send_inter_cohorts.sh
log_dir=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/lsf/

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
cmd="python ./intermediate_cohorts.py"
echo $cmd >>  ${launch_script}

echo "Output to ${log_file}"
echo "Errors to ${err_file}"
sbatch ${launch_script}
squeue -u ${USER} | grep ${job_name}
