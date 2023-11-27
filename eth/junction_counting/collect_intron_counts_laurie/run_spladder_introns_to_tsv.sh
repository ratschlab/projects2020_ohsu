#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/ICGC/alignments/alignments_gtex_2015-03-29_sparse

### filtered case
for fname in $(find ${basedir} -name \*.aligned.hdf5)
do
    if [ ! -f ${fname%hdf5}tsv.gz ]
    then
        #echo "python $(pwd)/spladder_introns_to_tsv.py ${fname}" | bsub -n 1 -M 1000 -o /dev/null -J hdf2tsv
        sbatch --cpus-per-task 1 --mem=1G --time=4:00:00 --output=/dev/null --job-name hdf2tsv --wrap "python $(pwd)/spladder_introns_to_tsv.py ${fname}"
    fi
done
