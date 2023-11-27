#!/bin/bash

set -e

for fname in /cluster/work/grlab/projects/TCGA/reprocess_2020/rna/results/alignments/*.all.hdf5
do
    outfile=${fname%.hdf5}.coverage_stats.tsv
    if [ ! -f ${outfile} ]
    then
        echo "python $(pwd)/summarize_coverage_hdf5.py $fname" | bsub -G ms_raets -n 1 -M 1000 -W 1:00 -o /dev/null -J stat_tcga
    fi
done
