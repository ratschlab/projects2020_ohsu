#!/bin/bash

# ----- SLURM Parameters -----
hours=24
cpus=30
mem=20G
tag='missing47'
suffix="pool_interm_${tag}"
job_name='bix_matrix'

launch_script=send_inter_cohorts.sh


# ----- ANALYSIS Parameters -----
sample_type='ov'

base_gtex='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref'
interm_gtex_cohort='ref_graph_kmer_normalized_filtered_10-21overlap_.gz'
metadata=$( echo 'kmer' 'coord' 'junctionAnnotated' 'readFrameAnnotated' 'isCrossJunction')
normalizer_libsize=400000
if [ "${sample_type}" == 'brca' ]; then 
	base_cancer='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374'
	target_samples=$(echo 'TCGAC8A12P01A11RA11507all', 'TCGAAOA0JM01A21RA05607all' 'TCGABHA18V01A11RA12D07all' 'TCGAA2A0D201A21RA03407all' 'TCGAA2A0SX01A12RA08407all')
	interm_cancer_cohort='ref_graph_kmer_normalized_filtered__.gz'
elif [ "${sample_type}" == 'ov' ]; then 
    	base_cancer='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102'
	target_samples=$(echo 'TCGA25131901A01R156513all' 'TCGA25131301A01R156513all' 'TCGA61200801A02R156813all' 'TCGA24143101A01R156613all' 'TCGA24229801A01R156913all')
	interm_cancer_cohort='ref_graph_kmer_normalized_filtered__.gz'
fi
path_cancer_libsize=${base_cancer}/expression_counts.libsize.tsv
intermediate_folder=${base_cancer}/filtering_intermediate
intermediate_output=${intermediate_folder}/complete_cancer_candidates_${tag}.tsv.gz
mkdir -p ${intermediate_folder}

# ----- Logging -----
log_dir=${base_cancer}/lsf
logging_base=${log_dir}/${suffix}_${cpus}p_${mem}
err_file=${logging_base}.err
log_file=${logging_base}.out

# ---- Launch run ----
echo "#!/bin/bash" > ${launch_script}
echo "#SBATCH -o ${log_file}" >> ${launch_script}
echo "#SBATCH -e ${err_file}" >> ${launch_script}
echo "#SBATCH -J ${job_name}" >> ${launch_script}
echo "#SBATCH --time=${hours}:00:00" >> ${launch_script}
echo "#SBATCH --cpus-per-task=${cpus}" >> ${launch_script}
echo "#SBATCH --mem=${mem}" >> ${launch_script}
cmd="python ./intermediate_cohorts.py --base-gtex ${base_gtex} --interm-gtex-cohort ${interm_gtex_cohort} --metadata ${metadata} --normalizer-libsize ${normalizer_libsize} --intermediate-output ${intermediate_output} --base-cancer ${base_cancer} --target-samples ${target_samples} --interm-cancer-cohort ${interm_cancer_cohort} --path-cancer-libsize ${path_cancer_libsize} "
echo $cmd >>  ${launch_script}

echo "Output to ${log_file}"
echo "Errors to ${err_file}"
sbatch ${launch_script}
squeue -u ${USER} | grep ${job_name}
