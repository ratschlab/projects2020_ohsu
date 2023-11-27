#!/bin/bash

set -e

#basedir=/cluster/work/grlab/projects/GTEx/rna/results/expression/simple_counting
basedir=/cluster/work/grlab/projects/GTEx/rna_gencode32_realign/results/expression/simple_counting

find $basedir -name \*.non_alt.tsv > ${basedir}.all_files.txt

cd ~/git/tools/gromics

if [ ! -f ${basedir}.all_counts.hdf5 ]
then
    python -m gromics.counting.collect_counts -f ${basedir}.all_files.txt -o ${basedir}.all_counts.hdf5
fi

python -m gromics.counting.count_hdf2tsv -i ${basedir}.all_counts.hdf5 -o ${basedir}.all_counts.tsv -v
