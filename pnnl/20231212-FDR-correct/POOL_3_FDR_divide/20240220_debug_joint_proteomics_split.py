#!/usr/bin/env python
# coding: utf-8
import pandas as pd
from Bio import SeqIO
import argparse
import os 
from collections import defaultdict
import numpy as np
import timeit
import argparse


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


def select_search_result_simple(id_to_SearchRow, valid_pep_ids):
    select_rows = set()
    for id_ in valid_pep_ids:
        select_rows.update(id_to_SearchRow[id_])
    return select_rows


def subset_pipeline(pipeline_index_mapping, joint_pepID_as_rows, search_res):

    indexes_from_pipeline = set([int(i.split('-')[1]) for i in pipeline_index_mapping['reindex']])
    indexes_joint_search = set(joint_pepID_as_rows.keys())
    subset_pipeline_pepID = indexes_joint_search.intersection(indexes_from_pipeline)

    print(f'Selected {len(subset_pipeline_pepID)} out of {len(indexes_joint_search)} peptide indexes searched')

    subset_pipeline_rows = select_search_result_simple(joint_pepID_as_rows, subset_pipeline_pepID)

    print(f'Selected {len(subset_pipeline_rows)} out of {len(search_res)} peptide rows searched')

    search_res_pipeline = search_res.iloc[list(subset_pipeline_rows), :]
    return search_res_pipeline


def replace_protein_id(search_res_pipeline, pipeline_index_mapping):
    # Generate a dictionary with peptide indexes and rows containing the given index
    ETH_pepID_as_rows = search_result_rows(search_res_pipeline)

    # Dictionary with the indexes from the joint pipeline file (reindexed) as keys
    # and the original indexes from either the ETH or the OHSU pipeline as values

    pipeline_index_mapping_dict = {}
    for joint_index, original_index in zip(pipeline_index_mapping['reindex'], pipeline_index_mapping['fasta_index']): #faster as iterrows
        pipeline_index_mapping_dict[int(joint_index.split('-')[1])] = int(original_index.split('-')[1])

    # Replace the joint IDs by the pipeline-specific ones in the output table restricted to the pipeline
    counter = 0
    indexes_lookup = []
    replacement_value = []
    print(f'Iterating over {len(pipeline_index_mapping_dict)} peptideIDs')
    for id_reindexed, original_index in pipeline_index_mapping_dict.items():
        counter +=1
        if counter % 5000 == 0:
            print(f'...{counter}')
        # Locate the protein IDs to replace    
        ids_to_replace = search_res_pipeline.iloc[ETH_pepID_as_rows[id_reindexed]]['protein id']

        # Generate the strings with the replaced IDs
        ids_back_to_original = [pepID.replace(str(id_reindexed), str(original_index)) 
                                for pepID in ids_to_replace]
        # store
        indexes_lookup.extend(list(ids_to_replace.index))
        replacement_value.extend(ids_back_to_original)


    #Replace the joint IDs by the pipeline-specific ones

    new_ids = pd.DataFrame.from_dict({'index': indexes_lookup, 'protein id new': replacement_value})
    search_res_pipeline = search_res_pipeline.reset_index()
    print('Size before merge', search_res_pipeline.shape)
    search_res_pipeline = search_res_pipeline.merge(new_ids, on = 'index', how = 'inner')
    print('Size after merge', search_res_pipeline.shape)
    search_res_pipeline = search_res_pipeline.rename({'protein id': 'joint id'}, axis = 1)
    search_res_pipeline = search_res_pipeline.rename({'protein id new': 'protein id'}, axis = 1)
    return search_res_pipeline


def tide_pipeline_split(search_res, ETH_index_mapping, OHSU_index_mapping, save_folder):


    # Reading
    search_res = pd.read_csv(search_res, sep = '\t')
    OHSU_index_mapping = pd.read_csv(OHSU_index_mapping, sep = '\t')
    ETH_index_mapping = pd.read_csv(ETH_index_mapping, sep = '\t')

    # Generate a dictionary with peptide indexes and rows containing the given index
    joint_pepID_as_rows = search_result_rows(search_res)

    print('Subsetting OHSU pipeline')
    search_res_OHSU = subset_pipeline(OHSU_index_mapping, joint_pepID_as_rows, search_res)
    print('Subsetting ETH pipeline')
    search_res_ETH = subset_pipeline(ETH_index_mapping, joint_pepID_as_rows, search_res)

    print('Replacing for OHSU pipeline')
    search_res_OHSU = replace_protein_id(search_res_OHSU, OHSU_index_mapping)
    print('Replacing for ETH pipeline')
    search_res_ETH = replace_protein_id(search_res_ETH, ETH_index_mapping)

    #Save
    base_name = 'FDR_file'
    extension = '.tsv'
    search_res_OHSU.to_csv(os.path.join(save_folder, base_name + '_OHSU' + extension), sep = '\t', index = None)
    search_res_ETH.to_csv(os.path.join(save_folder, base_name + '_ETH' + extension), sep = '\t', index = None)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='separates a joint reindexed file into one file per pipeline with the original indexes')
    parser.add_argument("--file-joint", help='tide search file on union of pipelines to separate in two pipelines')
    parser.add_argument("--map-eth-file", help='file for eth containing the mapping table between original ids and shared ids')
    parser.add_argument("--map-ohsu-file",help='file for ohsu containing the mapping table between original ids and shared ids')
    parser.add_argument("--save-folder",help='base folder to save results')
    args = parser.parse_args()
    print(args)
    tide_pipeline_split(args.file_joint, args.map_eth_file, args.map_ohsu_file, args.save_folder)

