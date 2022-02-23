import click
import pandas as pd
import tarfile
import glob 
import os 
import numpy as np 
import logging
import pickle


def ohsu_to_eth(path):
    cohort = {'NormalCohortcore_GTEx_': 'Gtexcore', 
             'NormalCohortpaired_': 'Matched', 
             'NormalCohortAll_':'Alls'}
    
    key_to_apply = [k for k in cohort if k in path]
    if key_to_apply: 
        key_to_apply = key_to_apply[0]
        path = path.replace(key_to_apply, '')
        path = path.replace('FiltNormalsC','FiltNormals{}C'.format(cohort[key_to_apply]) )

        path = path.replace('J_', 'G_')

        path = path.replace('CohortLim', '.0CohortLim')
        sample = path.split('_')[1]
        #print(path)
        return path, sample
    else:
        return None, None 
        
def get_eth_path(base_folder_ETH, name_eth=None, sample=None):
    path_o = None 
    if name_eth is not None: 
        path_list = os.path.join(base_folder_ETH, 'filter_' + sample, '*', name_eth, 'part*')
        path_list = glob.glob(path_list)
        if path_list:
            path_o = path_list[0]
    return path_o 

def set_stats(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    jaccard = float(intersection) / union
    return intersection, union, intersection/len(list1), intersection/len(list2), jaccard
     
def save_picke(i, sizes_ohsu, sizes_eth, sizes_intersect, 
              sizes_union, percent_ohsu_in_eth, percent_eth_in_ohsus, 
              jaccard):
    name = "/cluster/work/grlab/projects/projects2020_OHSU/plots/heatmaps/stats_single{}.pickle".format(i)
    filehandler = open(name.encode(),"wb")
    res = [sizes_ohsu ,
    sizes_eth ,
    sizes_intersect, 
    sizes_union ,
    percent_ohsu_in_eth, 
    percent_eth_in_ohsus,
    jaccard ]
    pickle.dump(res,filehandler)

@click.command()
@click.option('--io', default=1, help='OHSU pair number')
def main(io):
    tar_file_OHSU = '/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/Aug21_graph_data_current/OHSU_kmer_lists_Nov24.tar.gz'
    
    base_folder_ETH = '/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102'
    
    file_pair = {'eth':[], 'ohsu': []}
    with tarfile.open(tar_file_OHSU, "r:*") as tar:
        file_names_OHSU = tar.getnames()

        for name_ohsu in file_names_OHSU:
            name_eth, sample = ohsu_to_eth(name_ohsu)
            eth_path = get_eth_path(base_folder_ETH, name_eth, sample)
            if (eth_path is not None) and os.path.isfile(eth_path):
                file_pair['eth'].append(eth_path) 
                file_pair['ohsu'].append(name_ohsu)

    len_pairs = len(file_pair['eth'])
    sizes_ohsu = np.zeros(len_pairs)
    sizes_eth = np.zeros(len_pairs)
    sizes_intersect = np.zeros(len_pairs)
    sizes_union = np.zeros(len_pairs)
    percent_ohsu_in_eth = np.zeros(len_pairs)
    percent_eth_in_ohsu = np.zeros(len_pairs)
    jaccard = np.zeros(len_pairs)

    with tarfile.open(tar_file_OHSU, "r:*") as tar:
        name_ohsu = file_pair['ohsu'][io]
        for i_e, eth_path in enumerate(file_pair['eth']):
            if i_e >= io:
                df_ohsu = pd.read_csv(tar.extractfile(name_ohsu), header=0, sep="\t")
                df_eth = pd.read_csv(eth_path, sep="\t", usecols = ['kmer'])

                sizes_ohsu[io], sizes_eth[i_e] = len(df_ohsu['kmer']), len(df_eth['kmer'])

                sizes_intersect[i_e], sizes_union[i_e],\
                percent_ohsu_in_eth[i_e], percent_eth_in_ohsu[i_e],  \
                 jaccard[i_e] =  set_stats(df_ohsu['kmer'], df_eth['kmer'])
        save_picke(io, sizes_ohsu, sizes_eth, sizes_intersect, 
              sizes_union, percent_ohsu_in_eth, percent_eth_in_ohsu, 
              jaccard)



if __name__ == '__main__':
    main()
