#!/bin/bash
set -e


### Lsf and Run Parameters
mem=20000
time_=120
local_=run_cluster #"run_local"
parallel=6
#edge_or_segm=edge
suffix="commit83ac08b"
file_tmp='20211015'
#file_tmp=''
echo "WARNING check activation myimmuno3"

### Inputs

base_cancer=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102
path_cancer=${base_cancer}/cohort_mutNone
path_cancer_libsize=${base_cancer}/expression_counts.libsize.tsv

base_normal=/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_v2.5f0752a_pya.0.17.1_conf2_annot_ref_chrall_cap/
path_normal=${base_normal}/cohort_mutNone
path_normal_libsize=${base_normal}/expression_counts.libsize.tsv
whiltelist_normal=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/GTEX/GTEx_sample_IDs_10-2021_lib_graph_juliannelist
libsize_normal=/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_librarysize_pya.0.17.1_conf2_annot_ref_chrall_cap1000/cohort_mutNone/GTEX_hg38_coding_libsize75.tsv

uniprot=/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/uniprot/9mers_uniprot-human-UP000005640_9606.tsv

sample_type=BRCA
if [ ${sample_type} == 'OV' ] ; then 
	whitelist_cancer=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/TCGA_foreground/sample_full_Ov_378.tsv
	target_cancer=''
	libsize_cancer=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_libsizes_conf2_annotFrame_cap1000_runs_pya0.17.1_KEEP/TCGA_Ovarian_374/cohort_mutNone/TCGA_Ovarian_374_coding_libsize75.tsv
else
	whitelist_cancer=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/TCGA_foreground/sample_full_BRCA_1102.tsv
	target_cancer='' #TODO update target 
	libsize_cancer=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_libsizes_conf2_annotFrame_cap1000_runs_pya0.17.1_KEEP/TCGA_Breast_1102/cohort_mutNone/TCGA_Breast_1102_coding_libsize75.tsv
fi

#TODO Update right limit 
cohort_expr_lim_cancer='1'
expr_n_limit_cancer='2'
sample_expr_lim_cancer='2' #already conf2 
kmer='9'
expr_n_limit_normal='3' 
cohort_expr_lim_normal='1'
normalizer_cancer_libsize='400000'
normalizer_normal_libsize=${normalizer_cancer_libsize}

parallelism='10000'
out_partitions=1
scratch_mem=155000
test_='norm'

log_dir=${base_cancer}/lsf
mkdir -p ${log_dir}
### Main 
while read sample; do
	for mutation_canc in ref; do 
		## Organize folders
		
		## Generate instructions
		output_dir=${base_cancer}/filter_${sample}/${suffix}
		mkdir -p ${output_dir}	
		logfile=${log_dir}/${sample}.cancerspec.${suffix}.log
		echo $logfile
		# --path-cancer-matrix-segm PATH_CANCER_MATRIX_SEGM [PATH_CANCER_MATRIX_SEGM ...]
                #        segment expression integrated matrix of kmers * samples for foreground
 # --path-cancer-matrix-edge
 		#--path-normal-kmer-list $(sed -z 's/\n/ /g' ./tmp_input_annot_BRCA)
		cmd="immunopepper cancerspecif --cores $parallel --mem-per-core $mem --kmer $kmer --expression-fields-c "segmentExpr" "junctionExpr" --path-cancer-matrix-edge $(sed -z 's/\n/ /g' ./tmp${file_tmp}_input_Junc_BRCA) --ids-cancer-samples "${sample}" --mut-cancer-samples ${mutation_canc} --whitelist-cancer ${whitelist_cancer} --path-cancer-libsize ${libsize_cancer} --normalizer-cancer-libsize ${normalizer_cancer_libsize} --whitelist-normal ${whiltelist_normal} --path-normal-libsize ${libsize_normal} --normalizer-normal-libsize ${normalizer_normal_libsize} --output-dir $output_dir --cohort-expr-support-cancer ${cohort_expr_lim_cancer} --n-samples-lim-cancer ${expr_n_limit_cancer} --sample-expr-support-cancer ${sample_expr_lim_cancer} --uniprot ${uniprot} --parallelism ${parallelism} --out-partitions ${out_partitions} --path-normal-matrix-segm $(sed -z 's/\n/ /g' ./tmp${file_tmp}_input_Segm_GTEX2019) --path-normal-matrix-edge  $(sed -z 's/\n/ /g' ./tmp${file_tmp}_input_Junc_GTEX2019) --path-normal-kmer-list $(sed -z 's/\n/ /g' ./tmp${file_tmp}_input_annot_BRCA) --cohort-expr-support-norm ${cohort_expr_lim_normal} --n-samples-lim-normal ${expr_n_limit_normal} --scratch-dir 'TMPDIR' > ${base_cancer}/${sample}.${mutation_canc}.run_cancerspecif.${suffix}${test_}.log 2>&1" 
				
		## Run
		if [ "$local_" = "run_local" ] ; then
			echo "running local"
			echo $cmd
			#$cmd
		else
			echo $cmd
			echo $cmd | bsub -n ${parallel} -J Filtipp -W ${time_}:00 -R "rusage[mem=${mem}]" -R "span[hosts=1]" -R "rusage[scratch=$scratch_mem]" #-o $logfile #-e ${logfile}.e -o $logfile #-R "span[hosts=1]" -o $logfile
		fi

	done
done < /cluster/work/grlab/projects/projects2020_OHSU/sample_lists/TCGA_foreground/BRCA_5samples_spladder_full.csv #./tmp_samples
