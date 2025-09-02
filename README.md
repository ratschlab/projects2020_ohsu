This is the analysis code of the preprint from Prelot and David et. al.
The main analysis scripts are detailed, and key ressource files regarding samples and environements are provided

## PEPTIDE GENERATION STEPS

Generation of foreground peptides:
[projects2020\_ohsu/eth/peptides\_generation/foreground](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptides_generation/foreground)<br />
Generation of GTEX peptides
[projects2020\_ohsu/eth/peptides\_generation/background](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptides_generation/background)

## FILTERING STEPS

Codes in [projects2020\_ohsu/eth/peptides\_filtering/python\_pipeline](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptides_filtering/python_pipeline)<br />
Steps explained in [projects2020\_ohsu/eth/peptides\_filtering/python\_pipeline/README.txt](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptides_filtering/python_pipeline/README.txt)<br />
Note: This is heavily inspired from [projects2020\_ohsu/eth/peptides\_filtering/spark\_pipeline](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptides_filtering/spark_pipeline), but this pipeline has not been used to generate the final results

## FASTA GENERATION 
Files containing the pattern *pool*.fa <br />

Codes in
[projects2020\_ohsu/eth/peptide\_search\_format/notebooks/20230523\_speedup\_kmer\_to\_metadata\_matching.ipynb](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptide_search_format/notebooks/20230523_speedup_kmer_to_metadata_matching.ipynb)<br />
[projects2020\_ohsu/eth/peptide\_search\_format/notebooks/p20230320\_format\_kmers\_to\_peptides\_FASTA-may2023.ipynb](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptide_search_format/notebooks/p20230320_format_kmers_to_peptides_FASTA-may2023.ipynb)

## EXPERIMENT PER PEPTIDES 
Files containing the pattern *experiments\_per\_peptide* <br />

[projects2020\_ohsu/eth/peptide\_search\_format/notebooks/p20230503\_parse\_fasta\_to\_experiments\_per\_peptide.ipynb](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptide_search_format/notebooks/p20230503_parse_fasta_to_experiments_per_peptide.ipynb)<br />
[projects2020\_ohsu/eth/peptide\_search\_format/notebooks/helpers\_map.py](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptide_search_format/notebooks/helpers_map.py)
(The old development version is in
[projects2020\_ohsu/eth/debug/single\_case/p20230726\_extract\_exp\_map\_for\_fasta.ipynb](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/debug/single_case/p20230726_extract_exp_map_for_fasta.ipynb))

## EXPERIMENT MAPS 
Files containing the pattern *experiment\_map* <br />

[projects2020\_ohsu/eth/peptide\_search\_format/notebooks/p20230726\_parse\_experiment\_map.ipynb](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptide_search_format/notebooks/p20230726_parse_experiment_map.ipynb)

## COUNT FILES 
Files containing the pattern *\_samp\_chrt\_norm\_mot\_unip.tsv* <br />

Generated along the
[projects2020\_ohsu/eth/peptides\_filtering/python\_pipeline](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptides_filtering/python_pipeline)<br />
Added the star filtering to the files here
[projects2020\_ohsu/eth/plotting/OHSU\_collab\_Rerun\_allRF/20240204\_count\_files-for-julianne-figures.ipynb](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/plotting/OHSU_collab_Rerun_allRF/20240204_count_files-for-julianne-figures.ipynb)

## PROTEOMICS STEPS

Codes in [projects2020\_ohsu/pnnl](https://github.com/ratschlab/projects2020_ohsu/tree/master/pnnl)<br />
Steps explained in [projects2020\_ohsu/pnnl/README\_proteomics.txt](https://github.com/ratschlab/projects2020_ohsu/tree/master/pnnl/README_proteomics.txt)

## PLOTTING

### Parse the filtering results

[projects2020\_ohsu/eth/plotting/OHSU\_collab\_Rerun\_allRF\_proteomics/p20230315\_extract\_filtering\_results.ipynb](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/plotting/OHSU_collab_Rerun_allRF_proteomics/p20230315_extract_filtering_results.ipynb)

### Plot the filtering results (intersect the filtered OHSU and ETH kmers on a barplot)

[projects2020\_ohsu/eth/plotting/OHSU\_collab\_Rerun\_allRF\_proteomics/20240320\_barplot\_Experiment\_intersection-filtered-kmers-plottingOnly.ipynb](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/plotting/OHSU_collab_Rerun_allRF_proteomics/20240320_barplot_Experiment_intersection-filtered-kmers-plottingOnly.ipynb)

### Parse the proteomics results

[projects2020\_ohsu/eth/plotting/OHSU\_collab\_Rerun\_allRF\_proteomics/send\_e-p-k.sh](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/plotting/OHSU_collab_Rerun_allRF_proteomics/send_e-p-k.sh)<br />
[projects2020\_ohsu/eth/plotting/OHSU\_collab\_Rerun\_allRF\_proteomics/run\_extract-peptides-to-kmer.sh](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/plotting/OHSU_collab_Rerun_allRF_proteomics/run_extract-peptides-to-kmer.sh)<br />
[projects2020\_ohsu/eth/plotting/OHSU\_collab\_Rerun\_allRF\_proteomics/extract-peptides-to-kmer.py](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/plotting/OHSU_collab_Rerun_allRF_proteomics/extract-peptides-to-kmer.py)<br />
[projects2020\_ohsu/eth/plotting/OHSU\_collab\_Rerun\_allRF\_proteomics/run\_example.sh](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/plotting/OHSU_collab_Rerun_allRF_proteomics/run_example.sh)

### Plot the proteomics results (Intersect validated OHSU and ETH kmers on a barplot)

[projects2020\_ohsu/eth/plotting/OHSU\_collab\_Rerun\_allRF\_proteomics/20240116\_barplot\_Experiment\_intersection-peptides-to-kmer-plottingOnly.ipynb](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/plotting/OHSU_collab_Rerun_allRF_proteomics/20240116_barplot_Experiment_intersection-peptides-to-kmer-plottingOnly.ipynb)

### Helpers in

[projects2020\_ohsu/eth/plotting/OHSU\_collab\_Rerun\_allRF\_proteomics](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/plotting/OHSU_collab_Rerun_allRF_proteomics/)

### Perform the intersection of the generated peptides for each pipeline

[projects2020\_ohsu/eth/MY\_Master\_thesis\_rerun\_LP](https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/MY_Master_thesis_rerun_LP)


## RESSOURCES
The samples ID used from the GTEx cohort can be found in [projects2020\_ohsu/ressources/GTEx\_sample\_IDs\_10\-2021.tsv.txt](https://github.com/ratschlab/projects2020_ohsu/tree/master/ressources/GTEx_sample_IDs_10-2021.tsv.txt)
