#!/bin/bash
#SBATCH -o /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/lsf/py_filter_brca_large_30p_20G_0-60608.out
#SBATCH -e /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/lsf/py_filter_brca_large_30p_20G_0-60608.err
#SBATCH -J brca_0_60608
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=20G
python ./cohort_filter.py --processes 30 --start-id 0 --end-id 60608 --path-cohort /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/cohort_mutNone --whitelist-tag --whitelist /cluster/work/grlab/projects/projects2020_OHSU/sample_lists/TCGA_foreground/sample_full_BRCA_1102_format.tsv --path-libsize /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/expression_counts.libsize.tsv --normalizer-libsize 400000 --filters 0.0 1.0 2.0 3.0 5.0 10.0 --sample-pattern TCGA --do-normalize --do-overwrite
