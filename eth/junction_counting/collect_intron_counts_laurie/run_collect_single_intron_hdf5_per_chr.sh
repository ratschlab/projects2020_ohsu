#!/bin/bash

set -e

mem=6000
basedir=/cluster/work/grlab/projects/ICGC/alignments/alignments_gtex_2015-03-29_sparse

outdir=${basedir}_junctions_projected
mkdir -p $outdir
GLs=""
for i in $(seq 191 230) $(seq 232 249)
do
    GLs="$GLs GL000${i}.1"
done

for chr in $(seq 1 22) X Y MT NC hs37d5 $GLs
do
    outfname=${outdir}/all_junctions.projected.${chr}.hdf5
    logfname=${outdir}/all_junctions.projected.${chr}.log
    if [ ! -f ${outfname} ]
    then
        #echo "python $(pwd)/collect_single_intron_hdf5.py ${basedir} chr${chr} ${basedir}.all_junctions.all_coords.sorted.uniq.tsv.gz ${outfname}" | bsub -n 1 -J gtex_prj -M ${mem} -o /dev/null -We 8:00 -R "rusage[mem=${mem}]"
        sbatch --cpus-per-task=1 --mem=${mem} --time=24:00:00 --output=${logfname} --job-name=gtex_prj --wrap "python $(pwd)/collect_single_intron_hdf5.py ${basedir} ${chr} ${basedir}.junctions.all_coords.sorted.uniq.tsv.gz ${outfname}"
        #echo "python $(pwd)/collect_single_intron_hdf5.py ${basedir} ${chr} ${basedir}.junctions.all_coords.sorted.uniq.tsv.gz ${outfname}"
       # exit
    fi
done
