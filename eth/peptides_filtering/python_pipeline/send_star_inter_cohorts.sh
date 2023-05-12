#!/bin/bash
#SBATCH -o /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374/lsf/order_r_complete_starversion_1p_150G.out
#SBATCH -e /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374/lsf/order_r_complete_starversion_1p_150G.err
#SBATCH -J star_big_mx
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=150G
python ./star_intermediate_cohorts.py --path-star /cluster/work/grlab/projects/GTEx/rna_gencode32_realign/results --libsize /cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/ARCHIV_keep_runs/GTEX2019_commit_v3_TEST_merged3_372a147_medium_run_pya.0.17.1_conf2_annot_ref_chrall_cap/expression_counts.libsize.tsv --whitelist /cluster/work/grlab/projects/projects2020_OHSU/sample_lists/GTEX/GTEx_sample_IDs_10-2021_lib_graph_juliannelist --normalizer 400000 --big-matrix /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374/filtering_intermediate/complete_cancer_candidates_order_r_complete.tsv.gz --jx-target-list /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374/filtering_samples/filters_22March_order_wany_wAnnot/tmp_all_experiments_jx.txt --filter-thresholds 0.0 1.0 2.0 3.0 5.0 10.0
