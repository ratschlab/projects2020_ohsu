#!/bin/bash
#SBATCH -o /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/lsf/filter_param_1p_50G.out
#SBATCH -e /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/lsf/filter_param_1p_50G.err
#SBATCH -J filtering
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
python ./filtering_from_intermediate.py --basedir /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102 --intermediate-output /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/filtering_intermediate/complete_cancer_candidates_order_r.tsv.gz --filtering-id TEST --target-samples TCGA-AO-A0JM-01A-21R-A056-07.all TCGA-C8-A12P-01A-11R-A115-07.all TCGA-BH-A18V-01A-11R-A12D-07.all TCGA-A2-A0D2-01A-21R-A034-07.all TCGA-A2-A0SX-01A-12R-A084-07.all --Threshold-target 0.0 --Threshold-cancer-cohort 0.0 2.0 --N-samples-cancer None 1 5 --Threshold-normal-cohort 0.0 1.0 3.0 None --N-samples-normal 1 2 10 None
