#!/bin/sh

#log_dir=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/lsf
#suffix=test_run
#log_file=${log_dir}/py_filter_background_${suffix}.err
#res_file=${log_dir}/py_filter_background_${suffix}.out

#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=23000
#SBATCH --job-name="background"
#SBATCH --output=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/lsf/py_filter_background_test_run.err
#SBATCH --error=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/lsf/py_filter_background_test_run.err


#Specifies that the job will be requeued after a node


#TODO add the number of processes
python ./background_filter.py
