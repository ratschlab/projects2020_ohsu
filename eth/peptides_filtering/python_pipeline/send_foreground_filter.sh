#!/bin/bash
#SBATCH -o /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/lsf/py_filter_brca_large_30p_20G_0-None.out
#SBATCH -e /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/lsf/py_filter_brca_large_30p_20G_0-None.err
#SBATCH -J brca_0_None
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=20G
python ./foreground_filter.py --processes 30
