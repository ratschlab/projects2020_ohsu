#!/bin/bash

set -e

#for fname in /cluster/work/grlab/projects/GTEx/rna/results/alignments/*.all.hdf5
for fname in /cluster/work/grlab/projects/GTEx/rna_gencode32_realign/results/alignments/*.all.hdf5
do
    outfile=${fname%.hdf5}.coverage_stats.tsv
    if [ ! -f ${outfile} ]
    then
        echo "python $(pwd)/summarize_coverage_hdf5.py $fname" | bsub -G ms_raets -n 1 -M 1000 -W 1:00 -o /dev/null -J stat_gtex
    fi
done
