import pandas as pd
import glob
import os 
import numpy as np 
from collections import defaultdict
from Bio import SeqIO
import tarfile
import gzip



def explode_immunopepper_coord(mx, coord_col = 'coord', sep=':'):
    sep_out = '_'
    coord_mx = mx[coord_col].str.split(sep, expand=True) #7 min

    coord_mx[1] = coord_mx[1].astype('int')
    coord_mx[2] = coord_mx[2].astype('int')

    coord_mx['strand'] = None
    coord_mx.loc[coord_mx[1] < coord_mx[2], 'strand'] = '+'
    coord_mx.loc[coord_mx[1] > coord_mx[2], 'strand'] = '-'

    coord_mx['junction_coordinate1'] = None
    coord_mx['junction_coordinate2'] = None

    
    coord_mx = coord_mx.astype(str) # 7 min
    
    coord_mx['+first'] = coord_mx[1] + sep_out + coord_mx[2]
    coord_mx['+secon'] = coord_mx[3] + sep_out + coord_mx[4]
    coord_mx['-first'] = coord_mx[3] + sep_out + coord_mx[0]
    coord_mx['-secon'] = coord_mx[5] + sep_out + coord_mx[2]
    


    coord_mx.loc[(coord_mx['strand'] == '+'), 'junction_coordinate1'] = coord_mx['+first'] 
    coord_mx.loc[(coord_mx['strand'] == '-'), 'junction_coordinate1'] = coord_mx['-first'] 
    coord_mx.loc[(coord_mx['strand'] == '+') & (coord_mx[4] != 'None') , 'junction_coordinate2'] = coord_mx['+secon']
    coord_mx.loc[(coord_mx['strand'] == '-') & (coord_mx[4] != 'None') , 'junction_coordinate2'] = coord_mx['-secon']
    
    return coord_mx


def search_result_peptides_ids(df_search, col_seq):
    pep_IDs = []
    tryptic_pep = []
    for idx, peptide in zip(df_search['protein id'], df_search[col_seq] ):
        if idx is np.nan:
            print('ERROR: Search not successful on all fractions of sample. Please RERUN')
        for name_ in idx.split(','):
            if 'pepID' not in name_:
                continue
            pep_ix = int(name_.split('-')[1].replace('(1)', ''))
            pep_IDs.append(pep_ix)
            tryptic_pep.append(peptide)
    df = pd.DataFrame({ 'tryptic-pep-passFDR': tryptic_pep, 'pepID': pep_IDs }).drop_duplicates()

    return df


def get_pep_ids(fa_path):
    input_path = fa_path
    if '.gz' in input_path:
        fa_path = gzip.open(fa_path, 'rt')
        
    id_to_pep = {}
    for seq in SeqIO.parse(fa_path,'fasta'):
        pepID = int(seq.id.split(';')[0].split('-')[1])
        if str(seq.seq) in id_to_pep:
            assert(id_to_pep[str(seq.seq)] == pepID)
        id_to_pep[pepID] = str(seq.seq)
    
    if '.gz' in input_path:
        fa_path.close()
    return id_to_pep


def get_pep_coord(fa_path):
    input_path = fa_path
    if '.gz' in input_path:
        fa_path = gzip.open(fa_path, 'rt')
        
    jxs = []
    pepIDs = []
    biexon_peptides = []
    for seq in SeqIO.parse(fa_path,'fasta'):
        header = {}
        for field in seq.id.split(';'):
            last_item = field[-1:]
            fs = field.split('-')
            if len(fs) == 1:
                continue #ill formatting. Not '-' separator
            if last_item == '-':
                header[fs[0]] = fs[1] + last_item #correct for straind '-' being split
            else:
                header[fs[0]] = fs[1]

        jxs.append(header['jx_coord'].replace('|', ';'))
        pepIDs.append(int(header['pepID']))
        biexon_peptides.append(str(seq.seq))
    
    if '.gz' in input_path:
        fa_path.close()

    fasta_coordinates = pd.DataFrame({'jx': jxs, 'pepID': pepIDs, 'bi_exon_pep': biexon_peptides }).drop_duplicates()
    return fasta_coordinates


def tar_reader(tar_arxiv, file_name):
    with tarfile.open(tar_arxiv, "r:*") as tar: #OHSU
        df_ohsu = pd.read_csv(tar.extractfile(file_name), sep="\t")
        df_ohsu.reset_index(inplace = True)
        df_ohsu.columns = ['kmer', 'jx']
    return df_ohsu

def ohsu_tsv_reader(path, file_name):
    df_ohsu = pd.read_csv(os.path.join(path, file_name), sep="\t")
    df_ohsu.reset_index(inplace = True)
    df_ohsu.columns = ['kmer', 'jx']
    return df_ohsu


def kmer_in_bi_exon_peptide(df_kmer):
    '''Filters out kmers which are not in the bi-exon peptides. 
    This can happens is the kmer countains the junction but is not translated in the same reading frame. 
    In this case it needs to be match to the peptide in the correct reading frame'''
    indexes = []
    counter = 0
    for kmer, bi_exon_pep in zip(df_kmer['kmer'], df_kmer['bi_exon_pep'] ):
        if kmer is not np.nan: 
            if kmer in bi_exon_pep:
                indexes.append(counter)
        counter += 1

    df_kmer = df_kmer.iloc[indexes, :]
    return df_kmer


def validated_filtered_kmers(df_filtered, fasta_base_OHSU, kmer_files_OHSU, 
                             fasta_base_ETH, sample, experiment, pipeline, col_seq='unmodified sequence'):
    '''Takes validated tryptic peptides and outputs validated kmers'''
    # Get filtered peptides passing the FDR for the given experiment
    
    if (df_filtered is not None) and (df_filtered.shape[0]):
        df_filtered_summary = search_result_peptides_ids(df_filtered, col_seq)


        # Get Fasta file path 
        # Get kmer files path
        if pipeline == 'OHSU':
            fasta_file = f'{fasta_base_OHSU}/J_{sample}_pool_kmer.fa'
            file_kmer = f'J_{sample}_{experiment[1:]}.tsv' 
            #kmers_experiment = tar_reader(kmer_files_OHSU, file_kmer)
            kmers_experiment = ohsu_tsv_reader(kmer_files_OHSU, file_kmer)
        elif pipeline == 'ETH':
            fasta_base_ETH = fasta_base_ETH.replace('\'','') 
            fasta_file = f'{fasta_base_ETH}/G_{sample}_pool_kmer_25012024.fa.gz'
            file_kmer = os.path.join(fasta_base_ETH, f'G_{sample}_{experiment}.tsv.gz')
            file_kmer = glob.glob(file_kmer)[0]
            kmers_experiment = pd.read_csv(file_kmer, sep = '\t')
            #Extract coordinates
            coord_mx = explode_immunopepper_coord(kmers_experiment, coord_col = 'coord', sep=':')

            kmers_experiment = pd.concat([kmers_experiment, 
                                    coord_mx[['strand', 'junction_coordinate1', 'junction_coordinate2']]], axis = 1)
        else:
            print(f'ERROR {pipeline} not correctly defined')

        # Load fasta file
        fasta_file = glob.glob(fasta_file)[0]
        # Extract ids and peptides from fasta
        id_to_peptide_fasta = get_pep_ids(fasta_file)


        # Get the junction coordinates from the fasta file
        fasta_coordinates = get_pep_coord(fasta_file)


        # Add the junction coordinates to the filtered DF
        df_filtered_jx = fasta_coordinates.merge(df_filtered_summary, on = 'pepID', how = 'right')
        df_filtered_jx = df_filtered_jx.dropna(axis=0)
        test = [pep for pep, bi in zip(df_filtered_jx['tryptic-pep-passFDR'], df_filtered_jx['bi_exon_pep'] ) if pep in bi]
        issues = []
        c=0 
        for pep, bi in zip(df_filtered_jx['tryptic-pep-passFDR'], df_filtered_jx['bi_exon_pep'] ):
            if pep not in bi:
                issues.append(c)
            c+=1
        if df_filtered_jx.shape[0] !=len(test):
            
            assert(df_filtered_jx.shape[0] == len(test))

        # Add the kmers to the filtered DF
        if pipeline == 'OHSU':
            df_filtered_jx_kmers = kmers_experiment.merge(df_filtered_jx, on = 'jx', how = 'right')


        if pipeline == 'ETH':
            first_jx = kmers_experiment.merge(df_filtered_jx, left_on = 'junction_coordinate1', right_on = 'jx', how = 'inner') 
            second_jx = kmers_experiment.merge(df_filtered_jx, left_on = 'junction_coordinate2', right_on = 'jx', how = 'inner') 
            df_filtered_jx_kmers = pd.concat([first_jx, second_jx], axis = 0)

        # Remove kmers in different reading frames than the peptide 
        df_filtered_jx_kmers = kmer_in_bi_exon_peptide(df_filtered_jx_kmers)

        ratio = len(set(df_filtered_jx_kmers['kmer']))/ len(set(kmers_experiment['kmer'])) 
        val_rate = np.round(ratio * 100 , 2)
    else:
        df_filtered_jx_kmers = pd.DataFrame([], columns = ['kmer'])
        val_rate = 0 
    return df_filtered_jx_kmers, val_rate


def reader_assign_conf_pep(path, FDR_threshold, col_seq, col_qval):
    print(f'Reading {path}')
    if os.path.isfile(path):
        df = pd.read_csv(path, sep = '\t')
        tot_peptides = len(df[col_seq].unique())
        print(f'With Shape: {df.shape[0]}')
        print(f'With unique peptides: {tot_peptides}')
        assert('sequence' in df.columns)
        df_filtered = df.loc[df[col_qval] < FDR_threshold]
        print(f'Number of validated psm: {df_filtered.shape}')
        
        return df_filtered
    else:
        return None

    
def compare_OHSU_ETH(samples_store_pep, read_from_disk):
    '''
    Reads a nested dictionary and performs sets operations
    keys: samples
        keys: experiment (e.g. J02112GA)
            keys: pipeline (e.g. ETH, OHSU)
                values: kmer or peptide sets
    Returns a dictionary with peptides or kmers intersections and differences
    '''
    ## Compare peptides
    if read_from_disk:
        compare = {'sample' : [], 
                  'filter_' : [], 
                  'pep_size_ohsu' : [], 
                  'pep_size_eth' : [], 
                  'pep_size_intersection' : [], 
                  'pep_size_ohsu\eth' : [], 
                  'pep_size_eth\ohsu' : []}

        for sample, experiments_ in samples_store_pep.items():
            for experiment, pipelines_ in experiments_.items():
                if ('OHSU' in pipelines_.keys()) and ('ETH' in pipelines_.keys()):
                    compare['sample'].append(sample)
                    compare['filter_'].append(experiment)
                    compare['pep_size_ohsu'].append(len(pipelines_['OHSU']))
                    compare['pep_size_eth'].append(len(pipelines_['ETH']))
                    compare['pep_size_ohsu\eth'].append(len(pipelines_['OHSU'].difference(pipelines_['ETH'])))
                    compare['pep_size_eth\ohsu'].append(len(pipelines_['ETH'].difference(pipelines_['OHSU'])))
                    compare['pep_size_intersection'].append(len(pipelines_['ETH'].intersection(pipelines_['OHSU'])))
        compare = pd.DataFrame(compare)
        print('Data to save', compare.shape)

        return compare
   
    else:
        return None
    
def format_validation_rates(samples_store_rates, read_from_disk):
    if read_from_disk:
        compare = {'sample' : [], 
                  'filter_' : [], 
                  'pipeline': [], 
                   'validation_rate':[]}

        for sample, experiments_ in samples_store_rates.items():
            for experiment, pipelines_ in experiments_.items():
                if ('OHSU' in pipelines_.keys()) and ('ETH' in pipelines_.keys()):
                    for pipeline, rate in pipelines_.items():
                        compare['sample'].append(sample)
                        compare['filter_'].append(experiment)
                        compare['pipeline'].append( pipeline)
                        compare['validation_rate'].append(rate)
        compare = pd.DataFrame(compare)
        print('Data to save', compare.shape)
        return compare
    else:
        return None
