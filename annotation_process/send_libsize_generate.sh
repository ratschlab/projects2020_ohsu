#!/bin/bash
set -e

mem=20000
time_=24
parallel=1

### Immunopepper Run

cmd='python /cluster/home/prelotla/GitHub/projects2020_ohsu/annotation_process/libsize_generate.py' 

echo $cmd
echo $cmd | bsub -n ${parallel} -J libsize -W ${time_}:00 -R "rusage[mem=${mem}]" -R "span[hosts=1]"  #-o $logfile #-e ${logfile}.e -o $logfile #-R "span[hosts=1]" -o $logfile
