import os 
from pathlib import Path
import pandas as pd
from collections import defaultdict
import glob
import timeit
import argparse
import numpy as np

from helpers_split import reader_FDR_results, reader_experiments, experiments_maps, search_result_rows, select_search_result, reconstruct_experiment_FDR


def FDR_to_experiments(list_experiments, search_out_folder, save_folder, create_sample_subfolder, sample_search_out_folder):
    exp_all = reader_experiments(list_experiments)
    FDR_res = reader_FDR_results(search_out_folder, sample_search_out_folder)

    n_samples_process = 10
    for sample, FDR_file in FDR_res.items():
        print(sample)
        n_samples_process -= 1
        print(n_samples_process)

        print('...read FDR result')
        df_FDR = pd.read_csv(FDR_file, sep = '\t')
        print(df_FDR.shape)

        print('...extract rows IDS corresponding to peptides')
        id_to_SearchRow = search_result_rows(df_FDR)
        print(len(id_to_SearchRow))

        print('...process experiment map')
        id_to_pep, id_to_exp, exp_to_id = experiments_maps(exp_all[sample])
        print(len(id_to_exp))

        print('...select experiment rows')
        select_rows = select_search_result(id_to_exp, id_to_SearchRow)
        print(len(select_rows))

        print('...save experiments')
        reconstruct_experiment_FDR(select_rows, df_FDR, save_folder, sample, create_sample_subfolder=create_sample_subfolder)


        if n_samples_process < 1:
            break



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Takes results from FDR and splits them between experimental conditions')
    parser.add_argument("--list-experiments", help='file containing the paths to the experiment files per sample')
    parser.add_argument("--search-out-folder",help='path (with wildcards) of the tide search results folder containing the FDR file')
    parser.add_argument("--sample-search-out-folder", default=7, type=int,
                        help='position of the sample name in the path')
    parser.add_argument("--save-folder",help='base folder to save results')
    parser.add_argument("--create-sample-subfolder", default=False, 
                        help='if path is provided, creates a folder with the sample name and an additional subfolder when saving')


    args = parser.parse_args()
    print(args)
    FDR_to_experiments(args.list_experiments, args.search_out_folder, args.save_folder,
                       args.create_sample_subfolder, args.sample_search_out_folder)
