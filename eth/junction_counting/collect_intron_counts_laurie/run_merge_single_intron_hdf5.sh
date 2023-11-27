#!/bin/bash

set -e

mem=6000
basedir=/cluster/work/grlab/projects/ICGC/alignments/alignments_gtex_2015-03-29_sparse

outdir=${basedir}_junctions_projected
GLs=""
for i in $(seq 191 230) $(seq 232 249)
do
    GLs="$GLs GL000${i}.1"
done

chrms=""
for chrm in $(seq 1 22) X Y MT NC hs37d5 $GLs
do
    chrms="${chrms},${chrm}"
done
python $(pwd)/merge_single_intron_hdf5.py ${outdir} ${outdir}.all.hdf5 ${chrms:1}
#sbatch --cpus-per-task=1 --mem=${mem} --time=24:00:00 --output=${logfname} --job-name=gtex_prj --wrap "python $(pwd)/collect_single_intron_hdf5.py ${basedir} ${chr} ${basedir}.junctions.all_coords.sorted.uniq.tsv.gz ${outfname}"
