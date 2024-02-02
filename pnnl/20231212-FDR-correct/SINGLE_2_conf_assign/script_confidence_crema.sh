#!/bin/bash

#SBATCH --job-name=assignConfCrema
#SBATCH --output=crema.out
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=50GB

crema_script=$1
input_files=$2
eval_fdr=$3
threshold=$4
outdir=$5

python ${crema_script} --input-files ${input_files} --eval-fdr ${eval_fdr} --threshold ${threshold} --outdir ${outdir}
