#!/bin/bash
set -e


### Lsf and Run Parameters
mem=20000
time_=01
local_=run_cluster #"run_local"
parallel=2
#edge_or_segm=edge
suffix="smallGene_155m1_minimalPreprocess"
echo "WARNING check activation myimmuno3"
module load spark
### Inputs
uniprot=/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/uniprot/9mers_uniprot-human-UP000005640_9606.tsv
## Cancer Cohorts 

## Normal Cohorts
base_normal='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref'
sample_back='GTEXcore'
#libsize_normal=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_libsizes_conf2_annotFrame_cap1000_runs_pya0.17.1_KEEP/GTEX2019_copy/GTEX_hg38_TCGA_All_Normals_coding_libsize75.tsv
input_Segm_normal='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/mini_Segm_list.txt'
input_Junc_normal='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/mini_Junc_list.txt'
whitelist_normal='/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/GTEX/GTEx_sample_IDs_10-2021_lib_graph_juliannelist_noBrain_noTestis'
#whitelist_normal=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/GTEX/GTEx_sample_IDs_10-2021_lib_graph_juliannelist_noBrain_noTestis_TMP25perc
tag_normals='Gtexcore'
batch='False'

## Parameters

normalizer_cancer_libsize='400000'
normalizer_normal_libsize=${normalizer_cancer_libsize}
kmer='9'

#TODO adjust parallelism 
parallelism='1000'
out_partitions=1
scratch_mem=1G   #270000 #155000 # 270000 #155000

cohort_expr_lim_cancer='1'
sample_expr_lim_cancer='2' #already conf2 
expr_n_limit_normal='3' 
cohort_expr_lim_normal='1'
expr_n_limit_cancer='none'

log_dir=${base_normal}/lsf
mkdir -p ${log_dir}
### Main 
for mutation_canc in ref; do 
    ## Organize folders
    output_norm=${base_normal}/filtered_backgrounds/${suffix}
    output_dir=${output_norm}
    output_count=${output_norm}/summary_counts
    mkdir -p ${output_dir}	
    mkdir -p ${output_norm}
    mkdir -p ${output_count}
    logfile=${log_dir}/mode.cancerspec.${suffix}.${tag_normals}.log
    ## Extract parameters	
    #cohort_expr_lim_normal=$(echo $expr_nsamples_limit_normal | cut -f1 -d ',')
    #expr_n_limit_normal=$(echo $expr_nsamples_limit_normal | cut -f2 -d ',')
    
    
    ## Cmd
    cmd000="immunopepper cancerspecif --cores $parallel --mem-per-core $mem --kmer $kmer --whitelist-normal ${whitelist_normal} --output-dir $output_dir --parallelism ${parallelism} --out-partitions ${out_partitions} --path-normal-matrix-segm $(cat ${input_Segm_normal}) --path-normal-matrix-edge $(cat ${input_Junc_normal}) --interm-dir-norm ${output_norm} --tag-prefix 'G' --output-count ${output_count}" #TODO add back scratch for cancer? --scratch-dir 'TMPDIR'" #TODO output count remove? #TODO Uniprot add back #--path-normal-libsize ${libsize_normal}  --normalizer-normal-libsize ${normalizer_normal_libsize} --path-normal-kmer-list ${input_annot_cancer} 
    ## None case
    if [ ${cohort_expr_lim_normal} != 'Any' ]; then
        cmd00="${cmd000} --cohort-expr-support-norm ${cohort_expr_lim_normal}"
    else
        cmd00="${cmd000}"
    fi
    if [ ${expr_n_limit_normal} != 'Any' ]; then
        cmd0="${cmd00} --n-samples-lim-normal ${expr_n_limit_normal}"
    else
        cmd0="${cmd00}"
    fi
    
    if [ ${expr_n_limit_cancer} != 'none' ]; then 	
	cmd1="${cmd0} --cohort-expr-support-cancer ${cohort_expr_lim_cancer} --n-samples-lim-cancer ${expr_n_limit_cancer}"
    else 
	cmd1="${cmd0}"
    fi 
    
    ## Batch case
    if [ ${batch} == 'True' ]; then 
	cmd2="${cmd1} --tot-batches ${tot_batches} --batch-id nbtc --tag-normals ${tag_normals}nbtc"
    else
	cmd2="${cmd1} --tag-normals ${tag_normals}"
    fi
		
    output_file=${output_dir}/${mutation_canc}.cancerspecif.${suffix}.Nec${cohort_expr_lim_normal}.Nn${expr_n_limit_normal}_nbtc.log

if [ ! -f "${test_output_exist}/_SUCCESS" ] ; then
	echo "${parallelism} ${scratch_mem}" 	
	echo ${logfile}
	echo ${output_file}
	echo $cmd2
        echo '#!/bin/bash' > ./tmp_file
	echo $cmd2 >> ./tmp_file
	sbatch --job-name=filt${suffix}_GTEX2019 --cpus-per-task=${parallel} --time=${time_}:00:00 --mem=${mem} -e ${logfile} -o ${output_file} ./tmp_file
fi    
done
