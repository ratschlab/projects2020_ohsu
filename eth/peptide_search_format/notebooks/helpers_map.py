import tarfile
import pandas as pd
import glob
import os
import numpy as np
from collections import defaultdict
import gzip


def path_to_expID(path:str):
    ID = {}
    filter_background = os.path.basename(path).split('_')[-1].replace('.tsv.gz', '')
    filter_foreground = os.path.basename(path).split('_')[-2]
    #print(filter_background, filter_foreground)
    # Extract Values
    ID['filter_background_reads'] = filter_background.split('lim')[-1].split('Across')[0]
    ID['filter_background_samples'] = filter_background.split('Across')[-1]
    ID['filter_background_cohort'] = filter_background.split('Normals')[1].split('Cohort')[0]
    ID['filter_foreground_reads'] = filter_foreground.split('Lim')[-1].split('Across')[0]
    ID['filter_foreground_samples'] = filter_foreground.split('Across')[-1]
    ID['filter_foreground_target'] = (filter_foreground.split('Lim')[1].replace('SampleLim', '').replace('Cohort', ''))
    return ID
    

def expID_to_block(ID_split: dict):
    ID_dict = {'Any': 'A', 'None': 'N', '10': 'X', 'paired': 'P', 'Gtex': 'G', 'Core_GTEx': 'R', 'Full': 'F'}
    motif = 'A' # Any
    for k, v in ID_split.items():
        try:
            assert(v != '10') #replace roman letter
            ID_split[k] = str(int(float(v)))
        except:
            ID_split[k] = ID_dict[v]

    return ID_split['filter_foreground_target'] + \
    ID_split['filter_foreground_reads'] + \
    ID_split['filter_foreground_samples'] + \
    ID_split['filter_background_reads'] + \
    ID_split['filter_background_samples'] + \
    ID_split['filter_background_cohort'] + \
    motif 
    

def preprocess_fasta(fasta):
    '''fasta: str. Path fasta file'''
    print(f'Load {fasta} fasta')
    # Extract the peptides from the 1 target fasta from sample
    peptides_IDS = []
    peptides_sequences = []
    pep_position = []
    between_codons = []
    with gzip.open(fasta, 'rt') as f:
        for line in f.readlines():
            if '>' in line:
                pep_position.append(int(line.split(';')[1].split('-')[-1]))
                between_codons.append(int(line.split(';')[2].split('-')[-1]))
                peptides_IDS.append(line.split(';')[0].split('-')[-1])
            else:
                peptides_sequences.append(line.replace('\n',''))




    # Cut the peptides around the junction 
    peptides_IDS_expand = []
    cuts_peptides = []
    kmer_len = 9 
    for idx in np.arange(len(peptides_IDS)):
        for k in np.arange(kmer_len - between_codons[idx] ): 
            peptides_IDS_expand.append(peptides_IDS[idx])
            cut_kmer = peptides_sequences[idx][pep_position[idx] - k : pep_position[idx] + (kmer_len - k)] 
            cuts_peptides.append(cut_kmer)
    return peptides_IDS, peptides_sequences, peptides_IDS_expand, cuts_peptides

def experiments_preprocess(path_filtered_sample):
    '''path_filtered_sample: str. path for all experiements of sample'''
    # Extract the occurences of the kmers in the experimental files. 


    kmer_to_filesID = defaultdict(list)
    for filtered in path_filtered_sample: # All experiments for given sample
        # GET filtered file ID
        ID_split = path_to_expID(filtered)
        ID_EXPERIMENT = expID_to_block(ID_split)
        # Read filtered file 
        df_filt = pd.read_csv(filtered, sep = '\t')
        for km in set(df_filt['kmer']):
            kmer_to_filesID[km].append(ID_EXPERIMENT)

    kmer_to_filesID_ = []
    for k, v in kmer_to_filesID.items():  # Collapse 
        kmer_to_filesID_.append((k,';'.join(np.unique(v))))
    return kmer_to_filesID_