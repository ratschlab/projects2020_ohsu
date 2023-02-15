#!/bin/bash
#SBATCH -o /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/lsf/py_filter_gtex_large_TEST_30p_20G_0-50000.out
#SBATCH -e /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/lsf/py_filter_gtex_large_TEST_30p_20G_0-50000.err
#SBATCH -J back_0_50000
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=20G
python ./cohort_filter.py --processes 30 --start-id 0 --end-id 50000
