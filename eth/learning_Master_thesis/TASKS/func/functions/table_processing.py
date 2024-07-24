import pandas as pd
import tqdm
from config import *

def change_column_names(ohsu_df):
    sample_col = [col for col in ohsu_df if col.startswith('TCGA')]
    new_col = [s.replace('-', '') for s in sample_col]
    #new_col = [s.replace('-', '') + 'all' for s in sample_col]
    d = dict(zip(sample_col, new_col))
    ohsu_df.rename(columns = d, inplace = True)
    return ohsu_df


def preprocess_ohsu(df):
    pd.options.mode.chained_assignment = None #None - ignoring the warning, This line disables a warning that occurs when assigning values to a DataFrame. It prevents the warning from being displayed.
    df = df[df['gene_id'].notna()] #notna() - detect existing (non-missing) values. This line filters out rows where the 'gene_id' column is not null. It removes any rows with missing gene IDs.
    df = df[df['kmer'].notna()] #This line filters out rows where the 'kmer' column is not null. It removes any rows with missing kmer values.
    df.loc[:, 'kmer'] = df['kmer'].apply(lambda s: s.split(';')) #This line splits the values in the 'kmer' column by the ';' delimiter and assigns the resulting list to the 'kmer' column. It splits a string into multiple values.
    df = df.explode('kmer') #This line explodes the 'kmer' column, which means it creates multiple rows for each value in the 'kmer' column. It expands a list into separate rows.
    df.loc[:, 'gene_id'] = df['gene_id'].apply(lambda s: s.split(','))
    df = df.explode('gene_id')
    return df


def ohsu_to_eth_coord(df, col = 'kmer', new_col = 'jx_shifted', sep = ';'):
    tmp_jx = df[col].str.split(sep,  expand = True)
    df[new_col] = tmp_jx[0] + sep + (tmp_jx[1].astype(int) - 1).astype(str) + sep + tmp_jx[2] + sep + tmp_jx[3]
    return df

def get_junction_coordinates(df, coordinates_col, sep=':'):
    df['strand'] = None
    df['junction_coordinate'] = None

    for idx, row in tqdm.tqdm(df.iterrows()):
        kmer_coordinates = [int(x) for x in row[coordinates_col].split(sep) if x !='None']

        if kmer_coordinates[1] < kmer_coordinates[2]: # order strand +

            df.loc[idx, 'strand'] = '+'
            if len(kmer_coordinates) == 4:  # 2 exons
                df.loc[idx, 'junction_coordinate'] = ':'.join([str(x) for x in kmer_coordinates[1:3]])
            elif len(kmer_coordinates) == 6:
                df.loc[idx, 'junction_coordinate'] = ':'.join([str(x) for x in kmer_coordinates[1:5]])
        else: # order strand +
            df.loc[idx, 'strand'] = '-'
            if len(kmer_coordinates) == 4:  # 2 exons
                df.loc[idx, 'junction_coordinate'] = ':'.join([str(x) for x in [kmer_coordinates[3],
                                                                                kmer_coordinates[0]]])
            elif len(kmer_coordinates) == 6:
                df.loc[idx, 'junction_coordinate'] = ':'.join([str(x) for x in [kmer_coordinates[3],
                                                                                kmer_coordinates[0],
                                                                                kmer_coordinates[2],
                                                                                kmer_coordinates[5]
                                                                               ]])
    return df


def fil_definition(filter_foreground,filter_background):
    filfor = []
    filbac = []   
    
    
    # Processing of filters into a given parametr (1, 0 ,4) (0, 2)
    
    for ff in filter_foreground:
        if not OHSU_BRCA_NEW:
            filter_sample=re.findall(PATTERN,re.findall(SAMPLE_PATTERN,ff)[0])[0]
            filter_cohorlim=re.findall(PATTERN,re.findall(COHORTLIM_PATTERN,ff)[0])[0]
            filter_across=re.findall(PATTERN,re.findall(ACROSS_PATTERN,ff)[0])[0]
            filfor.append(f'({filter_sample}, {filter_cohorlim}, {filter_across})')
        else:
            filfor.append(ff)
    for fb in filter_background:
        if not OHSU_BRCA_NEW:
            filter_cohorlim=re.findall(PATTERN,re.findall(COHORTLIM_PATTERN,fb)[0])[0]   
            filter_across=re.findall(PATTERN,re.findall(ACROSS_PATTERN,fb)[0])[0]
            filbac.append(f'({filter_cohorlim}, {filter_across})')
        else:
            filbac.append(fb)
    return (filfor,filbac)


def common_genes(ohsu_df, eth_df, show_log = True):
    g1 = ohsu_df['gene_id'].unique()
    g2 = eth_df['gene_id'].unique()
    gene_names = list(set(g1) & set(g2))
    if show_log:
        print("OHSU number of genes:", len(g1))
        print("ETH number of genes:", len(g2))
        print("Number of common genes:", len(gene_names))
    return gene_names

def filter_df_common_kmers(df1, df2, show_log = True, col = 'kmer'):
    '''
    Filter input datframes, so that they have only records with common pairs gene + kmer
    '''
    if show_log:
        print('Initials sizes: ', len(df1), len(df2))
    gene_names = common_genes(df1, df2, show_log)
    df1 = df1[df1['gene_id'].isin(gene_names)]
    df1_grouped = df1.groupby('gene_id')
    df2 = df2[df2['gene_id'].isin(gene_names)]
    df2_grouped = df2.groupby('gene_id')
    dfs1 = []
    dfs2 = []
    for gene in tqdm.tqdm(gene_names):
        d1 = df1_grouped.get_group(gene)
        d2 = df2_grouped.get_group(gene)
        common_kmer = set(d1[col]) & set(d2[col])
        if len(common_kmer) == 0 and show_log:
            print("There are no common kmers for gene: ", gene)
        else:
            d1 = d1[d1[col].isin(common_kmer)]
            d2 = d2[d2[col].isin(common_kmer)]
            d1['gene_id'] = gene
            dfs1.append(d1)
            d2['gene_id'] = gene
            dfs2.append(d2)
    df1, df2 = pd.concat(dfs1, axis=0, ignore_index=True), pd.concat(dfs2, axis=0, ignore_index=True)
    if show_log:
        print('Result sizes: ', len(df1), len(df2))
    return df1, df2

def filter_df_common_coord(df1, df2):
    return filter_df_common_kmers(df1, df2, False, 'junction_coordinate')

def get_filter(filter):
    f =[]
    for ff in filter:
        if not OHSU_BRCA_NEW:
            filter_sample=re.findall(PATTERN,re.findall(SAMPLE_PATTERN,ff)[0])[0]
            filter_cohorlim=re.findall(PATTERN,re.findall(COHORTLIM_PATTERN,ff)[0])[0]
            filter_across=re.findall(PATTERN,re.findall(ACROSS_PATTERN,ff)[0])[0]
            f.append(f'({filter_sample}, {filter_cohorlim}, {filter_across})')
        else:
            f.append(ff)
    return f

def sorting(data):
    if PLOT_SORT_BY =='x-axis sorted by size of intersection':
    # Pandas sorting by intersection
        data=data.sort_values(by=['sample','size_intersection'], ascending=[False,False])
        print('I')
    
    elif PLOT_SORT_BY =='x-axis sorted by size of intersection and priority':
        data=data.sort_values(['size_intersection', 'priority'],ascending = [False, True])
        print('I+P')
        
    elif PLOT_SORT_BY=='x-axis sorted by size of JP-specific set':
        data=data.sort_values(by=['sample','size_ohsu\eth'], ascending=[False,False])
        print('JP')
        
    elif PLOT_SORT_BY=='x-axis sorted by size of JP-specific set and priority': 
        data=data.sort_values(['size_ohsu\eth', 'priority'],ascending = [False, True])
        print('JP+P')
        
    elif PLOT_SORT_BY=='x-axis sorted by size of GP-specific set':
        data=data.sort_values(by=['sample','size_eth\ohsu'], ascending=[False,False])
        print('GP')
        
    elif PLOT_SORT_BY=='x-axis sorted by size of GP-specific set and priority':
        data=data.sort_values(['size_eth\ohsu', 'priority'],ascending = [False, True])
        print('GP+P')
    else:
        print("Choose right filter")
    return data


def sorting_GP(data):
    if PLOT_GP_SORT_BY =='x-axis sorted by size of intersection':
    # Pandas sorting by intersection
        data=data.sort_values(by=['sample','GP from Inter NF'], ascending=[False,False])
        print('I')
    elif PLOT_GP_SORT_BY=='x-axis sorted by size of GP-specific set':
        data=data.sort_values(['sample', 'GP from GP NF'],ascending = [False, False])
        print('GP')
    else:
        print("Choose right filter")
    return data


def sorting_JP(data):
    if PLOT_JP_SORT_BY =='x-axis sorted by size of intersection':
    # Pandas sorting by intersection
        data = data.sort_values(by=['sample','JP from Inter NF'], ascending=[False,False])
        print('I')
    elif PLOT_JP_SORT_BY=='x-axis sorted by size of JP-specific set':
        data=data.sort_values(by=['sample','JP from JP NF'], ascending=[False,False])
        print('JP')
    else:
        print("Choose right filter")
    return data
