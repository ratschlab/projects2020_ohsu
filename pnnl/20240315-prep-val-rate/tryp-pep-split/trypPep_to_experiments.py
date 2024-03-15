import os 
from pathlib import Path
import pandas as pd
from collections import defaultdict
import glob
import timeit
import argparse
import numpy as np

from helpers_split import get_pep_ids_df, reader_experiments, experiments_maps, search_result_rows, select_search_result, reconstruct_experiment_FDR


def fasta_to_experiments(list_experiments, base_pipeline_folder, save_folder, create_sample_subfolder, samples):
    exp_all = reader_experiments(list_experiments)

    n_samples_process = 10
    for sample in samples:
        
        fa_path = os.path.join(base_pipeline_folder, sample, 'trypsine_digest/peptide-extracted-filter-unique.fasta')
        print(sample)
        n_samples_process -= 1
        print(n_samples_process)

        print('...read tryptic peptides')
        df_pep_sample = get_pep_ids_df(fa_path)

        print('...process experiment map')
        id_to_pep, id_to_exp, exp_to_id = experiments_maps(exp_all[sample])
        print(len(id_to_exp))

        print('...extract rows IDS corresponding to peptides')
        id_to_SearchRow = search_result_rows(df_pep_sample)
        print(len(id_to_SearchRow))

        print('...extract rows needed to reconstruct experiments')
        select_rows = select_search_result(id_to_exp, id_to_SearchRow)
        print(len(select_rows))


        print('...save experiments')
        reconstruct_experiment_FDR(select_rows, df_pep_sample, save_folder, sample, create_sample_subfolder=create_sample_subfolder)



        if n_samples_process < 1:
            break


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Takes a fasta and splits it between experimental conditions')
    parser.add_argument("--list-experiments", help='file containing the paths to the experiment files per sample')
    parser.add_argument("--base-pipeline-folder", help='path containing the result for the pipeline. Should contain {sample}/trypsine_digest/peptide-extracted-filter-unique.fasta')
    parser.add_argument("--samples", nargs='+',
                        help='list of samples to process')
    parser.add_argument("--save-folder",help='base folder to save results')
    parser.add_argument("--create-sample-subfolder", default=False, 
                        help='if path is provided, creates a folder with the sample name and an additional subfolder when saving')


    args = parser.parse_args()
    print(args)
    fasta_to_experiments(args.list_experiments, args.base_pipeline_folder, args.save_folder,
                       args.create_sample_subfolder, args.samples)
