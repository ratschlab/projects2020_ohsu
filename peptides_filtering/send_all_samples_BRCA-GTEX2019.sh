#!/bin/bash
set -e


### Lsf and Run Parameters
mem=20000
time_=75
local_=run_cluster #"run_local"
parallel=6
#edge_or_segm=edge
suffix="commit_0a02cfe"
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
input_Junc_cancer=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102/cohort_mutNone_relink/Junc_BRCA_19077
input_annot_cancer=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102/cohort_mutNone_relink/annot_BRCA_19077
input_Segm_normal=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102/cohort_mutNone_relink/Segm_GTEX2019_19077
input_Junc_normal=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102/cohort_mutNone_relink/Junc_GTEX2019_19077

normalizer_cancer_libsize='400000'
normalizer_normal_libsize=${normalizer_cancer_libsize}
kmer='9'

parallelism='1000'
out_partitions=1
scratch_mem=155000
test_=''

cohort_expr_lim_cancer='1'
expr_n_limit_cancer='2'
sample_expr_lim_cancer='2' #already conf2 
expr_n_limit_normal='3' 
cohort_expr_lim_normal='1'


log_dir=${base_cancer}/lsf
mkdir -p ${log_dir}
### Main 
##TODO add argument core whitelist; all normal subset and all normals with whitelist
for cohort_expr_lim_cancer in '1' ; do #'none' '5'; do 
	for expr_n_limit_cancer in '1'; do #'2' '10'; do 
		for sample_expr_lim_cancer in '2' ; do #'0' ; do #other type of splicing graph
			for expr_n_limit_normal in '1' '2' '10'; do 
				for cohort_expr_lim_normal in '0' '3' '10'; do  
				   for background in 'gtex'; do 
					while read sample; do
						for mutation_canc in ref; do 
							## Organize folders
							
							## Generate instructions
							output_dir=${base_cancer}/filter_${sample}/${suffix}_${background}
							mkdir -p ${output_dir}	
							logfile=${log_dir}/${sample}.cancerspec.${suffix}.log
							echo $logfile
							
							cmd="immunopepper cancerspecif --cores $parallel --mem-per-core $mem --kmer $kmer --expression-fields-c "segmentExpr" "junctionExpr" --path-cancer-matrix-edge ${input_Junc_cancer} --ids-cancer-samples "${sample}" --mut-cancer-samples ${mutation_canc} --whitelist-cancer ${whitelist_cancer} --path-cancer-libsize ${libsize_cancer} --normalizer-cancer-libsize ${normalizer_cancer_libsize} --whitelist-normal ${whiltelist_normal} --path-normal-libsize ${libsize_normal} --normalizer-normal-libsize ${normalizer_normal_libsize} --output-dir $output_dir --cohort-expr-support-cancer ${cohort_expr_lim_cancer} --n-samples-lim-cancer ${expr_n_limit_cancer} --sample-expr-support-cancer ${sample_expr_lim_cancer} --uniprot ${uniprot} --parallelism ${parallelism} --out-partitions ${out_partitions} --path-normal-matrix-segm ${input_Segm_normal} --path-normal-matrix-edge ${input_Junc_normal} --path-normal-kmer-list ${input_annot_cancer}" #TODO add back scratch for cancer? --scratch-dir 'TMPDIR'" 
							
							if [ ${cohort_expr_lim_cancer} != 'none' ]; then 	
								cmd1="${cmd} --cohort-expr-support-norm ${cohort_expr_lim_normal} --n-samples-lim-normal ${expr_n_limit_normal}"
							else
								cmd1="${cmd}"	
							fi 
							cmd2="${cmd1} > ${base_cancer}/${sample}.${mutation_canc}.run_cancerspecif.${suffix}.Cec${cohort_expr_lim_cancer}.Cn${expr_n_limit_cancer}.Ces.${sample_expr_lim_cancer}.Nec${cohort_expr_lim_normal}.Nn${expr_n_limit_normal}.${test_}.log 2>&1"	
							## Run
							if [ "$local_" = "run_local" ] ; then
								echo "running local"
								echo $cmd2
								#$cmd
							else
								echo $cmd2
								echo $cmd2 | bsub -n ${parallel} -J Filtipp -W ${time_}:00 -R "rusage[mem=${mem}]" -R "span[hosts=1]" -R "rusage[scratch=$scratch_mem]" -o $logfile #-e ${logfile}.e -o $logfile #-R "span[hosts=1]" -o $logfile
							fi

						done
					done < ./tmp_samples #/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/TCGA_foreground/BRCA_5samples_spladder_full.csv #./tmp_samples
				done
			done
		done
	done
done
done
