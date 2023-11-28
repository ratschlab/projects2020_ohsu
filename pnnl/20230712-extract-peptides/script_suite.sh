#!/bin/bash

#SBATCH --job-name=trypsine
#SBATCH --output=trypsine.out
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=1GB

scripts_home=$1
logfile=$2
input=$3

is_zip=$(echo $input | grep -c gz)
if [ $is_zip -eq 1 ]; then
    zcat $input > tmp
    input=tmp
fi
python3 ${scripts_home}/extract-peptides.py $input > ${logfile}/log.txt
python3 ${scripts_home}/filter-peptides.py peptide-extracted.fasta >> ${logfile}/log.txt
python3 ${scripts_home}/plot-pep-cdf.py peptide-extracted-filter.fasta >> ${logfile}/log.txt
rm tmp
