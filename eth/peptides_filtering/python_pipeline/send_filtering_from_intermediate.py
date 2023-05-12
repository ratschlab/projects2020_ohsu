#!/bin/bash
#SBATCH -o /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/lsf/filter_param_1p_50G.out
#SBATCH -e /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/lsf/filter_param_1p_50G.err
#SBATCH -J filtering
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
python ./star_intermediate_cohorts.py --basedir /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102 --intermediate-output /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/filtering_intermediate/complete_cancer_candidates_order_r.tsv.gz --filtering-id filters_22March_order_wany_wAnnot_REPRODUCE --target-samples TCGA-25-1319-01A-01R-1565-13.all TCGA-25-1313-01A-01R-1565-13.all TCGA-61-2008-01A-02R-1568-13.all TCGA-24-1431-01A-01R-1566-13.all TCGA-24-2298-01A-01R-1569-13.all --Threshold-target 0.0 --Threshold-cancer-cohort 3424780{Threshold_cancer_cohort} --N-samples-cancer None, 1, 5 --Threshold-normal-cohort 0.0, 1.0, 3.0, None --N-samples-normal 1 2 10 None
