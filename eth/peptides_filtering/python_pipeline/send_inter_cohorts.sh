#!/bin/bash
#SBATCH -o /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/lsf/pool_interm_order_r_1p_120G.out
#SBATCH -e /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/lsf/pool_interm_order_r_1p_120G.err
#SBATCH -J bix_matrix
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=120G
python ./intermediate_cohorts.py --base-normal /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref --interm-normal-cohort ref_graph_kmer_normalized_filtered_10-21overlap_order_.gz --metadata kmer coord junctionAnnotated readFrameAnnotated isCrossJunction --normalizer-libsize 400000 --intermediate-output /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/filtering_intermediate/complete_cancer_candidates_order_r.tsv.gz --base-cancer /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102 --cancer-targets TCGAC8A12P01A11RA11507all TCGAAOA0JM01A21RA05607all TCGABHA18V01A11RA12D07all TCGAA2A0D201A21RA03407all TCGAA2A0SX01A12RA08407all --interm-cancer-cohort ref_graph_kmer_normalized_filtered__.gz --path-cancer-libsize /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/expression_counts.libsize.tsv
