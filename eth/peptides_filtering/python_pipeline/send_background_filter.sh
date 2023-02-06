#!/bin/sh

#log_dir=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/lsf
#suffix=test_run
#log_file=${log_dir}/py_filter_background_${suffix}.err
#res_file=${log_dir}/py_filter_background_${suffix}.out

#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=1G
#SBATCH --job-name="background"
#SBATCH --output=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/lsf/py_filter_background_map_async_rew_dummy.out
#SBATCH --error=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/lsf/py_filter_background_map_async_rew_dummy.err


#Specifies that the job will be requeued after a node


#TODO add the number of processes
python ./background_filter.py --processes 10
