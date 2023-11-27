#!/bin/bash

set -e

#basedir=/cluster/work/grlab/projects/GTEx/rna_gencode32_realign/results
basedir=/cluster/work/grlab/projects/ICGC/alignments/alignments_gtex_2015-03-29_sparse

### collect all coordinates and make unique
echo collecting coordinates
zcat ${basedir}/SRR1068687.aligned.tsv.gz | head -n 1 | cut -f -4 | gzip > ${basedir}.junctions.all_coords.sorted.uniq.tsv.gz
(zcat ${basedir}/*.aligned.tsv.gz | grep -v junction_start | cut -f -4 | gzip > ${basedir}.junctions.all_coords.tsv.gz) && (zcat ${basedir}.junctions.all_coords.tsv.gz | sort -u | gzip >> ${basedir}.junctions.all_coords.sorted.uniq.tsv.gz) && rm ${basedir}.junctions.all_coords.tsv.gz
echo done
