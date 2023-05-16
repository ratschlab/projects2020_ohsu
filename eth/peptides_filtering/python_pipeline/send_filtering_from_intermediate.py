#!/bin/bash
#SBATCH -o /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374/lsf/filter_param_1p_50G.out
#SBATCH -e /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374/lsf/filter_param_1p_50G.err
#SBATCH -J filtering
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
python ./filtering_from_intermediate.py --basedir /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374 --intermediate-output /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374/filtering_intermediate/complete_cancer_candidates_order_r_complete.tsv.gzp --filtering-id filters_15May_order_wany_cp_wAnnot --target-samples TCGA-25-1319-01A-01R-1565-13.all TCGA-25-1313-01A-01R-1565-13.all TCGA-61-2008-01A-02R-1568-13.all TCGA-24-1431-01A-01R-1566-13.all TCGA-24-2298-01A-01R-1569-13.all --Threshold-target 0.0 --Threshold-cancer-cohort None 0.0 2.0 --N-samples-cancer None 1 5 --Threshold-normal-cohort 0.0 1.0 3.0 None --N-samples-normal 1 2 10 None
