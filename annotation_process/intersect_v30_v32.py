import pandas as pd 
import os 

v32 = '/cluster/work/grlab/projects/projects2020_OHSU//annotation/gencode.v32.annotation.gtf'
genes_30 = '/cluster/work/grlab/projects/projects2020_OHSU/gene_lists/all_genes_gencode_v30.txt'
path_custom_annot = '/cluster/work/grlab/projects/projects2020_OHSU/annotation/gencode.v32_IntersectGenesInV30.gtf'
genes_accepted = pd.read_csv(genes_30, header = None)[0].values


kept = 0 
skipped = 0 
with open(path_custom_annot, 'w') as fp:
    for line in open(v32, 'r'):
            if line[0] == '#':
                continue

            gene = line.strip().split('\t')[8].split(';')[0].split('"')[1]
            if gene not in genes_accepted:
                skipped +=1
                continue

            fp.write(line)
            kept +=1
print("kept {}".format(kept))
print("skipped {}".format(skipped))
