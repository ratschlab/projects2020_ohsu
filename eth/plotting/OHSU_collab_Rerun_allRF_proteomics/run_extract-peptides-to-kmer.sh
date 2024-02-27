#!/bin/bash

#SBATCH --job-name=extract_plotting
#SBATCH --output=extract_plotting.out
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --exclude=compute-biomed-10,compute-biomed-01

script_dir=/cluster/home/prelotla/github/projects2020_ohsu/eth/plotting/OHSU_collab_Rerun_allRF_proteomics
proteomicsdir='/cluster/work/grlab/projects/projects2020_OHSU/proteomics_fixMerge_25012024'
fasta_base_OHSU='/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/current'
fasta_base_ETH="'/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/*/filtering_samples/filters_19May_order_5ge_wAnnot_GPstar'"
kmer_files_OHSU='/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/current/kmer_files'
save_folder='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/data_plots'
FDR_limit=0.05

MS_FDR=$1
MS_strategy=$2

cmd="python ${script_dir}/extract-peptides-to-kmer.py --proteomicsdir ${proteomicsdir} --samples-breast TCGA-C8-A12P-01A-11R-A115-07 TCGA-AO-A0JM-01A-21R-A056-07 TCGA-BH-A18V-01A-11R-A12D-07 TCGA-A2-A0D2-01A-21R-A034-07 TCGA-A2-A0SX-01A-12R-A084-07 --samples-ov TCGA-25-1319-01A-01R-1565-13 TCGA-25-1313-01A-01R-1565-13 TCGA-61-2008-01A-02R-1568-13 TCGA-24-1431-01A-01R-1566-13 TCGA-24-2298-01A-01R-1569-13 --fasta-base-OHSU ${fasta_base_OHSU} --fasta-base-ETH ${fasta_base_ETH} --kmer-files-OHSU ${kmer_files_OHSU} --pipelines ETH OHSU --FDR-limit ${FDR_limit} --MS-FDR ${MS_FDR} --MS-strategy ${MS_strategy} --save-folder ${save_folder}"

echo $cmd
rm ${script_dir}/run_example.sh
echo $cmd > ${script_dir}/run_example.sh
$cmd
