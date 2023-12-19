#!/bin/bash
crux_home=/cluster/home/prelotla/util/crux-4.1.Linux.x86_64/bin/crux
output_dir=.
search_output=/cluster/work/grlab/projects/projects2020_OHSU/proteomics/tide_search_joint/TCGA-BH-A18V/TCGA_BH-A18V_A7-A13F_BH-A0E1_117C_W_BI_20130520_H-PM_f10/tide-search.txt
#TODO Pool the partitions for union, pool partitions for ETH, pool partitions for OHSU
sbatch script_confidence.sh  ${crux_home} ${output_dir} ${search_output}
