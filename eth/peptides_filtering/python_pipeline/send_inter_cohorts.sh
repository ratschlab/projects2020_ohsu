#!/bin/bash
#SBATCH -o /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374/lsf/pool_interm_order_r_1p_120G.out
#SBATCH -e /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374/lsf/pool_interm_order_r_1p_120G.err
#SBATCH -J bix_matrix
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=120G
python ./intermediate_cohorts.py --base-normal /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref --interm-normal-cohort ref_graph_kmer_normalized_filtered_10-21overlap_order_.gz --metadata kmer coord junctionAnnotated readFrameAnnotated isCrossJunction --normalizer-libsize 400000 --intermediate-output /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374/filtering_intermediate/complete_cancer_candidates_order_r.tsv.gz --base-cancer /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374 --cancer-targets TCGA25131901A01R156513all TCGA25131301A01R156513all TCGA61200801A02R156813all TCGA24143101A01R156613all TCGA24229801A01R156913all --interm-cancer-cohort ref_graph_kmer_normalized_filtered_full_.gz --path-cancer-libsize /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374/expression_counts.libsize.tsv
