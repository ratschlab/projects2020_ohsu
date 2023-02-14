#!/bin/bash
#SBATCH -o /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/lsf/py_filter_gtex_large_30p_20G_30000-50000.out
#SBATCH -e /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/lsf/py_filter_gtex_large_30p_20G_30000-50000.err
#SBATCH -J back_30000_50000
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=20G
python ./background_filter.py --processes 30 --start-id 30000 --end-id 50000
