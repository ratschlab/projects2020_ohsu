#!/usr/bin/env python
# coding: utf-8
import pandas as pd
from Bio import SeqIO
import argparse
import os 

def get_pep_ids(fa_path):
    Ids = set()
    pep_to_id = {}
    for seq in SeqIO.parse(fa_path,'fasta'):
        Ids.add(seq.seq)
        assert(seq.seq not in pep_to_id)
        pep_to_id[seq.seq] = seq.id
        
    return Ids, pep_to_id
    
def fasta_reindex(fasta_eth, fasta_ohsu, save_folder, path_map_eth, path_map_ohsu):

    # Parse IDS
    print(f'Reading {fasta_eth} \n {fasta_ohsu}')
    pep_eth, pep_to_ID_eth = get_pep_ids(fasta_eth)
    pep_ohsu, pep_to_ID_ohsu = get_pep_ids(fasta_ohsu)

    # Compare
    joined_pep = pep_eth.intersection(pep_ohsu)
    eth_pep = pep_eth.difference(pep_ohsu)
    ohsu_pep = pep_ohsu.difference(pep_eth)
    print('shared peptides', len(pep_ohsu.union(pep_eth)))

    # Reindex
    file_write = []
    shared_to_ohsu = {}
    shared_to_eth = {}
    for idx, pep in enumerate(pep_ohsu.union(pep_eth)):
        id_text = f'pepID-{idx}'
        id_string = '>' + id_text
        file_write.append(id_string)
        file_write.append(str(pep))
        if pep in pep_ohsu:
            shared_to_ohsu[id_text] = pep_to_ID_ohsu[pep]
        if pep in pep_eth:
            shared_to_eth[id_text] = pep_to_ID_eth[pep]

    # Write result
    path_save = os.path.join(save_folder, 'joint-peptide-extracted-filter-unique.fasta')
    with open(path_save, 'w') as f:
        for i, line in enumerate(file_write):
            f.write(line + '\n')
    print(f'Reindexed joint fasta saved in {path_save}')

    # Save OHSU
    path = os.path.join(path_map_ohsu, 'pepID_joint_original.tsv.gz')
    df = pd.DataFrame.from_dict(shared_to_ohsu, orient='index').reset_index()
    df.columns = ['reindex', 'fasta_index']
    df.to_csv(path, sep = '\t', index = None)
    print(f'Map for OHSU save in {path}')

    # Save ETH
    path = os.path.join(path_map_eth, 'pepID_joint_original.tsv.gz')
    df = pd.DataFrame.from_dict(shared_to_eth, orient='index').reset_index()
    df.columns = ['reindex', 'fasta_index']
    df.to_csv(path, sep = '\t', index = None)
    print(f'Map for ETH save in {path}')



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Reindexes the peptides coming from two pipelines while taking into account their overlap')
    parser.add_argument("--file-eth", help='fasta file from eth pipeline with peptide IDs as handles')
    parser.add_argument("--file-ohsu",help='fasta file from ohsu pipeline with peptide IDs as handles')
    parser.add_argument("--map-eth-folder", help='folder for eth where to save the mapping table between original ids and shared ids')
    parser.add_argument("--map-ohsu-folder",help='folder for ohsu where to save the mapping table between original ids and shared ids')
    parser.add_argument("--save-folder",help='base folder to save results')
    args = parser.parse_args()
    print(args)
    fasta_reindex(args.file_eth, args.file_ohsu, args.save_folder, args.map_eth_folder, args.map_ohsu_folder)
