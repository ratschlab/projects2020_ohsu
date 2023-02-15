#!/bin/bash

start_gene=0
end_gene=50000
hours=04
cpus=30
mem=20G

job_name=back_${start_gene}_${end_gene}
suffix=py_filter_gtex_large_TEST

log_dir=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/lsf
logging_base=${log_dir}/${suffix}_${cpus}p_${mem}_${start_gene}-${end_gene}
err_file=${logging_base}.err
log_file=${logging_base}.out
launch_script=send_cohort_filter.sh

echo "#!/bin/bash" > ${launch_script}
echo "#SBATCH -o ${log_file}" >> ${launch_script}
echo "#SBATCH -e ${err_file}" >> ${launch_script}
echo "#SBATCH -J ${job_name}" >> ${launch_script}
echo "#SBATCH --time=${hours}:00:00" >> ${launch_script}
echo "#SBATCH --cpus-per-task=${cpus}" >> ${launch_script}
echo "#SBATCH --mem=${mem}" >> ${launch_script}
echo "python ./cohort_filter.py --processes ${cpus} --start-id ${start_gene} --end-id ${end_gene}" >>  ${launch_script}

echo "Output to ${log_file}"
echo "Errors to ${err_file}"
sbatch ${launch_script}
squeue -u ${USER} | grep ${job_name}
