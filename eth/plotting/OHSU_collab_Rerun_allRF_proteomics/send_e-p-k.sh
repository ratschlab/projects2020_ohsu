#!/bin/bash

script_folder=${PWD}
save_folder='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/data_plots'


for MS_FDR in '_crema' '_crux'; do 
    for MS_strategy in 'pool' 'joint' 'single'; do
    	mkdir -p ${save_folder}/run_${MS_strategy}${MS_FDR}
	cd ${save_folder}/run_${MS_strategy}${MS_FDR}
        sbatch ${script_folder}/run_extract-peptides-to-kmer.sh ${MS_FDR} ${MS_strategy} 
    done
done

cd - 
