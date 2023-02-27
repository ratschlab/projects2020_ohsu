#!/bin/bash
#SBATCH -o /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/lsf//pool_interm_filtering_79_30p_20G.out
#SBATCH -e /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/lsf//pool_interm_filtering_79_30p_20G.err
#SBATCH -J bix_matrix
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=20G
python ./intermediate_cohorts.py
