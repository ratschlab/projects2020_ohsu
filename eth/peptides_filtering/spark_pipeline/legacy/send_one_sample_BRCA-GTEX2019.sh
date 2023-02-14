#!/bin/bash
set -e

mem=20000
time_=120
local_=submit  #run_local  # "run_local"
parallel=8
suffix="test_run"
log_dir=./logs_${suffix}
mkdir -p ${log_dir}

### Immunopepper Run 
echo "WARNING check activation myimmuno3"

cancer_type=BRCA
sample=TCGA-A2-A0SX.all
confidence=spladder_confidence_1
frames=all_frames
path_cancer=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/10c3360_runs_pya0.17.1/${cancer_type}/${sample}/${confidence}/${frames}
path_cancer_libsize=${path_cancer}/expression_counts.libsize.tsv
path_normal=/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_1fc5828_pya.0.17.1_ref_mode2_linked_part_snappy
path_normal_libsize=/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_1fc5828_pya.0.17.1_ref_mode2/libsize_dummy.tsv #expression_counts.libsize.tsv
path_normal_matrix_segm=${path_normal}/SegmExpr_mtx
path_normal_matrix_edge=${path_normal}/JuncExpr_mtx
uniprot=/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/uniprot/9mers_uniprot-human-UP000005640_9606.tsv
whitelist_normals=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/GTEX/GTEx_normal_samples_12-2020_whitelist.tsv
interm_file=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/10c3360_runs_pya0.17.1/BRCA/TCGA-A2-A0SX.all/spladder_confidence_1/all_frames/TCGA-A2-A0SX.all/filter_output/interm_normals_segm-edge_max_expr-in-1-samples-with-0.0-normalized-cts.tsv



expr_high_limit_normal='0'
expr_limit_normal='0' #TODO add 2 here 
expr_n_limit='1' # 10% of samples
kmer='9'
expr_limit_cancer='0'
parallelism='10000'
#scratch_dir=/cluster/work/grlab/projects/tmp_laurie/nobackup/tmp_cancerspecif_filt
out_partitions=2000
scratch_mem=100000

output_dir=${path_cancer}/${sample}/filter_output
mkdir -p $output_dir
#mkdir -p $scratch_dir

# group the 20 conditions as input for one filter background / match the GTEX all samples 
# remove the annotation
# right normalisation scheme 
# later somatic mutations as input  --> maybe soon 
file_interest=ref_sample_9mer.pq.gz
cmd="immunopepper cancerspecif --cores $parallel --mem-per-core $mem --kmer $kmer --expression-fields-c "segment_Expr" "junction_Expr" --path-cancer-libsize $path_cancer_libsize --path-normal-libsize ${path_normal_libsize} --paths-cancer-samples ${path_cancer}/${sample}/${file_interest} --ids-cancer-samples "${sample}" --path-normal-matrix-segm $(readlink -f ${path_normal_matrix_segm}/* | paste -s ) --path-normal-matrix-edge $( readlink -f  ${path_normal_matrix_edge}/* | paste -s ) --output-dir $output_dir --expr-high-limit-normal $expr_high_limit_normal --expr-high-limit-normal $expr_high_limit_normal --expr-limit-normal $expr_limit_normal --n-samples-lim-normal $expr_n_limit --expr-limit-cancer ${expr_limit_cancer} --uniprot ${uniprot} --parallelism ${parallelism} --out-partitions ${out_partitions}" # --path-normal-kmer-list ${interm_file} " #--whitelist ${whitelist_normals}" # --statistical"

			
	if [ "$local_" = "run_local" ] ; then
		echo "running local"
		$cmd
	else
		echo $cmd
		echo $cmd | bsub -n ${parallel} -G ms_raets -J TMP_FULL_all -W ${time_}:00 -R "rusage[mem=${mem}]" -R "span[hosts=1]" -R "rusage[scratch=$scratch_mem]" #-o $logfile #-e ${logfile}.e -o $logfile #-R "span[hosts=1]" -o $logfile
	fi
#done < ../TCGA-BRCA/tmp_samples #TODO real samples
