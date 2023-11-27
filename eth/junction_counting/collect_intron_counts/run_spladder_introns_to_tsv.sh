#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/GTEx/rna_gencode32_realign/results/alignments

### filtered case
for fname in $(find ${basedir} -name \*.filt.hdf5)
do
    if [ ! -f ${fname%hdf5}tsv.gz ]
    then
        echo "python $(pwd)/spladder_introns_to_tsv.py ${fname}" | bsub -n 1 -M 1000 -o /dev/null -J hdf2tsv
    fi
done

### unfiltered case
for fname in $(find ${basedir} -name \*.all.hdf5)
do
    if [ ! -f ${fname%hdf5}tsv.gz ]
    then
        echo "python $(pwd)/spladder_introns_to_tsv.py ${fname}" | bsub -n 1 -M 1000 -o /dev/null -J hdf2tsv
    fi
done
