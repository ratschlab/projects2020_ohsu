
dir_links=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/commit_v3_TEST_merged3_372a147_medium_run_conf2_annotFrame_cap0_runs/selected_genes_rerun1

for file_to_link in $(ls ${dir_links}/*select*); do 
	outdir=$( echo ${file_to_link} | sed 's/\.txt//g')
	mkdir ${outdir}
	cd ${outdir}
	while read brca; do ln -s $brca $(basename $(dirname $brca) | cut -f5 -d '_')_$(basename $brca)  ;done < ${file_to_link}
	cd -
done


