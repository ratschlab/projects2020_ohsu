#!/bin/bash
#SBATCH -o /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/lsf/py_filter_gtex_large_30p_20G_35000-45000.out
#SBATCH -e /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/lsf/py_filter_gtex_large_30p_20G_35000-45000.err
#SBATCH -J back_35000_45000
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=20G
python ./background_filter.py --processes 30 --start-id 35000 --end-id 45000
