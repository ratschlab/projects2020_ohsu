import pandas as pd 
import h5py 
import numpy as np 
import logging
import os 

def compute_library_size(path_expr, path_gene_list, path_save, percentile, query_no_version = False ):
    logging.info("Path of graph with gene expressions:  {}".format(path_expr))
    h5_gene_expr = h5py.File(open(path_expr, 'rb'))
    coding_list = pd.read_csv(path_gene_list, header = None, names = ['genes'])
    logging.info("Querying coding genes")
    if query_no_version: 
        h5_coding_mx = [list(h5_gene_expr["raw_count"][gene_id, :])
                  for gene_id, gene_name in enumerate(h5_gene_expr['gene_ids'])
                  if gene_name.decode().split('.')[0] in coding_list['genes'].values]
    else:
        h5_coding_mx = [list(h5_gene_expr["raw_count"][gene_id, :])
                  for gene_id, gene_name in enumerate(h5_gene_expr['gene_ids'])
                  if gene_name.decode() in coding_list['genes'].values]
    h5_coding_mx = np.array(h5_coding_mx)
    logging.info("Shape of coding matrix is: {}".format(h5_coding_mx.shape))
    logging.info("Number of coding genes is: {}".format(coding_list.shape))
    sample_names = [strain.decode() for strain in h5_gene_expr['strains']]
    logging.info("Computing library percentiles")
    libsize_count = pd.DataFrame({'sample': sample_names})
    for perc in percentile:
        libsize_count['libsize_{}percent'.format(perc)] = np.round(np.percentile(h5_coding_mx, perc, axis = 0 ), 1)
    logging.info("Saving to {}".format(path_save))
    libsize_count.to_csv(path_save, index = None, sep = '\t')

percentile=[ 50, 60, 70, 75, 80, 90, 95 ]  
outdir = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220'

path_expr = '/cluster/work/grlab/projects/GTEx/rna/results/splicing/spladder/genes_graph_conf2.merge_graphs.count.gene_expression.hdf5' # V30 
path_gene_list = '/cluster/work/grlab/projects/projects2020_OHSU/gene_lists/genes_coding_gencode_v32_inter_v30_wo_version.txt'
tag = ".coding_gencode_v32_inter_v30_wo_version.txt"
path_save = os.path.join(outdir, 'libsize{}.tsv'.format(tag)) 
query_no_version = True 
compute_library_size(path_expr, path_gene_list, path_save, percentile, query_no_version)


path_expr = '/cluster/work/grlab/projects/GTEx/rna/results/splicing/spladder/genes_graph_conf2.merge_graphs.count.gene_expression.hdf5'
path_gene_list = '/cluster/work/grlab/projects/projects2020_OHSU/gene_lists/genes_coding_gencode_v32_inter_v30.txt'
tag = ".coding_gencode_v32_inter_v30.txt"
path_save = path_save = os.path.join(outdir, 'libsize{}.tsv'.format(tag)) 
compute_library_size(path_expr, path_gene_list, path_save, percentile )

path_expr = '/cluster/work/grlab/projects/GTEx/rna/results/splicing/spladder/genes_graph_conf2.merge_graphs.count.gene_expression.hdf5'
path_gene_list = '/cluster/work/grlab/projects/projects2020_OHSU/gene_lists/genes_coding_gencode_v30.txt'
tag = ".coding_gencode_v30.txt"
path_save = path_save = os.path.join(outdir, 'libsize{}.tsv'.format(tag))
compute_library_size(path_expr, path_gene_list, path_save, percentile )

path_expr = '/cluster/work/grlab/projects/GTEx/rna/results/splicing/spladder/genes_graph_conf2.merge_graphs.count.gene_expression.hdf5'
path_gene_list = '/cluster/work/grlab/projects/projects2020_OHSU/gene_lists/genes_coding_gencode_v32.txt'
tag = ".coding_gencode_v32.txt"
path_save = path_save = os.path.join(outdir, 'libsize{}.tsv'.format(tag))
compute_library_size(path_expr, path_gene_list, path_save, percentile )

path_expr = '/cluster/work/grlab/projects/TCGA/PanCancer/rerun2018_alt_splice_GTEx/spladder/genes_graph_conf2.merge_graphs.validated.count.gene_expression.hdf5'
path_gene_list = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/coding_genes_nocap_GTEX2017'
tag = ".coding_gencode_v19.txt"
path_save = path_save = os.path.join(outdir, 'libsize{}.tsv'.format(tag)) 
compute_library_size(path_expr, path_gene_list, path_save, percentile )
