#!/bin/bash

#SBATCH --job-name=assignConfCrema
#SBATCH --output=crema.out
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=50GB

crema_script=$1
input_files=$2
outdir=$3

echo "python ${crema_script} --input-files ${input_files} --outdir ${outdir}"
python ${crema_script} --input-files ${input_files} --outdir ${outdir}
