#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/GTEx/rna_gencode32_realign/results

outfile=${basedir}/alignment.coverage_stats.tsv.gz
rm -f $outfile

for fname in $(find ${basedir}/alignments -name \*.all.coverage_stat.tsv | sort)
do
    tmp=$(tail -n1 $fname | cut -f 2-)
    sample=$(basename $fname | cut -f 1 -d '.')
    echo -e "${sample}\t${tmp}" | gzip >> $outfile
done
