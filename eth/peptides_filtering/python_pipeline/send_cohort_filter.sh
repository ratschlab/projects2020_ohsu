#!/bin/bash
#SBATCH -o /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/lsf/py_filter_gtex_large_mem_100p_60G_0-60609.out
#SBATCH -e /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/lsf/py_filter_gtex_large_mem_100p_60G_0-60609.err
#SBATCH -J back_0_60609
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=100
#SBATCH --mem=60G
python ./cohort_filter.py --processes 100 --start-id 0 --end-id 60609 --path-cohort /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/cohort_mutNone --whitelist-tag 10-21overlap --whitelist /cluster/work/grlab/projects/projects2020_OHSU/sample_lists/GTEX/GTEx_sample_IDs_10-2021_lib_graph_juliannelist --path-libsize /cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/ARCHIV_keep_runs/GTEX2019_commit_v3_TEST_merged3_372a147_medium_run_pya.0.17.1_conf2_annot_ref_chrall_cap/expression_counts.libsize.tsv --normalizer-libsize 400000 --filters 0.0 1.0 2.0 3.0 5.0 10.0 --sample-pattern SRR --do-normalize
