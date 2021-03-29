import pandas as pd 
import h5py 
import numpy as np 
import logging


def compute_library_size(path_expr, path_gene_list, path_save, percentile ):
    logging.info("Path of graph with gene expressions:  {}".format(path_expr))
    h5_gene_expr = h5py.File(open(path_expr, 'rb'))
    coding_list = pd.read_csv(path_gene_list, header = None, names = ['genes'])
    logging.info("Querying coding genes")
    h5_coding_mx = [list(h5_gene_expr["raw_count"][gene_id, :])
                  for gene_id, gene_name in enumerate(h5_gene_expr['gene_ids'])
                  if gene_name.decode() in coding_list['genes'].values]
    h5_coding_mx = np.array(h5_coding_mx)
    logging.info("Shape of coding matrix is: {}".format(h5_coding_mx.shape))
    logging.info("Number of coding genes is: {}".format(coding_list.shape))
    sample_names = [strain.decode() for strain in h5_gene_expr['strains']]
    logging.info("Computing library percentiles")
    libsize_count = pd.DataFrame({'sample': sample_names,
                                 'libsize_75percent': np.percentile(h5_coding_mx, percentile, axis = 0 )})
    logging.info("Saving")
    libsize_count.to_csv(path_save, index = None, sep = '\t')

percentile=75
path_expr = '/cluster/work/grlab/projects/GTEx/rna/results/splicing/spladder/genes_graph_conf2.merge_graphs.count.gene_expression.hdf5'
path_gene_list = '/cluster/work/grlab/projects/projects2020_OHSU/gene_lists/genes_coding_gencode_v32_inter_v30.txt'
path_save = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_1fc5828_pya.0.17.1_ref_mode2'
tag = ".coding_gencode_v32_inter_v30.txt"
path_save = 'expression_counts.libsize{}.tsv'.format(tag)
compute_library_size(path_expr, path_gene_list, path_save, percentile )

path_expr = '/cluster/work/grlab/projects/GTEx/rna/results/splicing/spladder/genes_graph_conf2.merge_graphs.count.gene_expression.hdf5'
path_gene_list = '/cluster/work/grlab/projects/projects2020_OHSU/gene_lists/genes_coding_gencode_v30.txt'
path_save = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_1fc5828_pya.0.17.1_ref_mode2'
tag = ".coding_gencode_v30.txt"
path_save = 'expression_counts.libsize{}.tsv'.format(tag)
compute_library_size(path_expr, path_gene_list, path_save, percentile )

path_expr = '/cluster/work/grlab/projects/TCGA/PanCancer/rerun2018_alt_splice_GTEx/spladder/genes_graph_conf2.merge_graphs.validated.count.gene_expression.hdf5'
path_gene_list = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/coding_genes_nocap_GTEX2017'
path_save = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/GTEX2017_commit_1fc5828_pya.0.17.1_ref'
tag = ".coding_gencode_v19.txt"
path_save = 'expression_counts.libsize{}.tsv'.format(tag)
compute_library_size(path_expr, path_gene_list, path_save, percentile )
