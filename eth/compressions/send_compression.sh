#!/bin/bash
#SBATCH -o /lsf/compress_99_1p_4G.out
#SBATCH -e /lsf/compress_99_1p_4G.err
#SBATCH -J c_${tag}
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
