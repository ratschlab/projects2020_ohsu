import os 
from pathlib import Path
import pandas as pd
from collections import defaultdict
import glob
import timeit
import argparse
import numpy as np

from helpers_psm import reader_tide_results, reader_experiments, experiments_maps, search_result_rows, select_search_result, reconstruct_experiment


def psm_to_experiments(list_experiments, search_out_folder, save_folder, create_sample_subfolder, rerank_psm):
    exp_all = reader_experiments(list_experiments)
    search_res = reader_tide_results(search_out_folder)

    n_samples_process = 10 
    for sample, partitions in search_res.items():
        if len(partitions) == 24:
            print(sample)
            n_samples_process -= 1
            print(n_samples_process)

            print('...read search result')
            df_search = pd.concat([pd.read_csv(part, sep = '\t') for part in partitions])
            print(df_search.shape)

            print('...extract rows IDS corresponding to peptides')
            id_to_SearchRow = search_result_rows(df_search)
            print(len(id_to_SearchRow))

            print('...process experiment map')
            id_to_pep, id_to_exp, exp_to_id = experiments_maps(exp_all[sample])
            print(len(id_to_exp))

            print('...select experiment rows')
            select_rows = select_search_result(id_to_exp, id_to_SearchRow)
            print(len(select_rows))
                                      
            print('...save experiments')
            reconstruct_experiment(select_rows, df_search, save_folder, sample, rerank=rerank_psm)


        if n_samples_process < 1:
            break




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Takes results from tide search and splits them between experimental conditions')
    parser.add_argument("--list-experiments", help='file containing the paths to the experiment files per sample')
    parser.add_argument("--search-out-folder",help='path (with wildcards) of the tide search results')
    parser.add_argument("--save-folder",help='base folder to save results')
    parser.add_argument("--create-sample-subfolder", default=False, action='store_true', 
                        help='wheather to create a subfolder with the sample name when saving')
    parser.add_argument("--rerank-psm", default=False, 
                        action='store_true',  
                        help='wheather to apply re-ranking of the psm within condition and partition')
    args = parser.parse_args()
    print(args)
    psm_to_experiments(args.list_experiments, args.search_out_folder, args.save_folder,
                       args.create_sample_subfolder, args.rerank_psm)
