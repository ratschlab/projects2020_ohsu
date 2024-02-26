#!/bin/bash

save_folder='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/data_plots'

cd ${save_folder}

for MS_FDR in '_crema' '_crux'; do 
    for MS_strategy in 'pool' 'joint' 'single'; do 
        sbatch run_extract-peptides-to-kmer.sh ${MS_FDR} ${MS_strategy} 
    done
done

cd - 
