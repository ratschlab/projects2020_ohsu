#!/bin/bash
set -e


time_=120
mem=40000
echo "sh compute_libsizes_distribute_v3.sh "| bsub -J libsizes -W ${time_}:00 -R "rusage[mem=${mem}]" 
