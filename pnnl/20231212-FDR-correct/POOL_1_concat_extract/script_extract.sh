#!/bin/bash
#SBATCH --job-name=FDRpool/union
#SBATCH --output=concat.out
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=1GB

#'''This scripts concatenates the search results across all partitions (except the reference) and extract the relevant peptides (discards the neighbors)'''
search_out_folder=$1

cd ${search_out_folder}
counter=0 
res=tide-search-concat.txt
rm ${search_out_folder}/$res
echo "output is ${search_out_folder}/${res}"

for folder in TCGA*; do 
    if [[  $(echo ${folder} | grep -c fA) -eq 0 ]] && [[ $(echo ${folder} | grep -c POOL) -eq 0  ]]; then 
        counter=$((counter+1)) 
	if [[ ${counter} -eq 1 ]] ;then 
	    head -1 ${folder}/tide-search.txt >> ${res}
	fi 
        tail -n +2 ${folder}/tide-search.txt | grep pepID >> ${res}
    fi
done
