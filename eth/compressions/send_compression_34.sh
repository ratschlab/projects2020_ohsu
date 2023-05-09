#!/bin/bash
#SBATCH -o /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/lsf/compress_34_1p_4G.out
#SBATCH -e /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/lsf/compress_34_1p_4G.err
#SBATCH -J c_34
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
while read line; do for part in ${line}/*Expr*/* ; do new_name=$(echo $part | sed 's,\.gz,\.lz,g') ; cat $part | gzip -d -c | plzip -6 -o $new_name ; if [ -f $new_name ] ; then echo $part; chmod 440 ${new_name} ; chmod 660 ${part}; rm ${part} ;fi; done; done < files_34.txt
