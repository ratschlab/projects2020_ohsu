import os 
from pathlib import Path
import pandas as pd
from collections import defaultdict
import glob
import timeit
import numpy as np


def reader_FDR_results(search_out_folder, sample_search_out_folder):
    search_res = dict()
    for path in glob.glob(search_out_folder):
        sample = path.split('/')[sample_search_out_folder]
        search_res[sample] = path 
    return search_res

def reader_experiments(list_experiments):
    '''Read files with path lists'''
    with open(list_experiments, 'r') as f:
        path_dict = {}
        for i in f.readlines():
            sample_name = os.path.basename(i.strip()).split('_')[1]
            sample_short = '-'.join(sample_name.split('-')[0:3])
            path_dict[sample_short] = i.strip()
    return path_dict


def experiments_maps(path):
    '''Extract experiment maps'''
    df = pd.read_csv(path, sep = '\t')
    id_to_pep = {}
    id_to_exp = {}
    exp_to_id = defaultdict(list)
    for i, row in df.iterrows():
        id_to_pep[row['peptide_id']] = row['peptide_sequence']
        id_to_exp[row['peptide_id']] = row['experiment_ids'].split(';')

    for k, v in id_to_exp.items():
        for ID in v:
            exp_to_id[ID].append(k)
    return id_to_pep, id_to_exp, exp_to_id


def search_result_rows(df_search):
    id_to_row = defaultdict(list)
    for i, idx in enumerate(df_search['protein id']):
        if idx is np.nan:
            print('ERROR: Search not successful on all fractions of sample. Please RERUN')
        for name_ in idx.split(','):
            if 'pepID' not in name_:
                continue
            pep_ix = int(name_.split('-')[1].replace('(1)', ''))
            id_to_row[pep_ix].append(i)
    return id_to_row


def select_search_result(id_to_exp, id_to_SearchRow):
    select_rows = defaultdict(set)
    for pep_idx, exp_list in id_to_exp.items():
        for experiment in exp_list:
            peptide_rows = id_to_SearchRow[pep_idx]
            if peptide_rows:
                select_rows[experiment].update(peptide_rows)
    return select_rows


def reconstruct_experiment_FDR(select_rows_pipeline, df_search, save_folder, sample, create_sample_subfolder):
    '''Selects all the rows from the initial experiment
    Save'''
    df_search_i = df_search.reset_index()
    for experiment_id in select_rows_pipeline:
        print(f'.....{experiment_id}')

        df_experiment = df_search_i.loc[select_rows_pipeline[experiment_id]]
        df_experiment = df_experiment.drop_duplicates()
        

        if create_sample_subfolder:
            path_save = os.path.join(save_folder, sample, create_sample_subfolder) 
        else:
            path_save = os.path.join(save_folder, sample)

        Path(path_save).mkdir(parents=True, exist_ok=True)
        path_save = os.path.join(path_save, f'tsearch-{experiment_id}.txt')
        print(path_save)
        df_experiment.to_csv(path_save, sep = '\t', index=None)
