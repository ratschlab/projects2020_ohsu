#!/usr/bin/env python
# coding: utf-8

# In[7]:


import pandas as pd
import os 


# In[8]:


v30 = '/cluster/work/grlab/projects/GTEx/annotation/gencode.v30.annotation.gtf'
v32 = '/cluster/work/grlab/projects/projects2020_OHSU//annotation/gencode.v32.annotation.gtf'
outpath = '/cluster/work/grlab/projects/projects2020_OHSU/gene_lists'


# ## Functions for analysis 

# In[9]:


def append_non_duplic(gene_to_transcript_dict, gene_id, transcript_id):
    if gene_id in gene_to_transcript_dict and transcript_id not in gene_to_transcript_dict[gene_id]:
        gene_to_transcript_dict[gene_id].append(transcript_id)
    else:
        gene_to_transcript_dict[gene_id] = [transcript_id]


# In[10]:


def process_anno(ann_path):
    transcript_to_gene_dict = {}    # transcript -> gene id
    gene_to_transcript_dict = {}    # gene_id -> list of transcripts
    transcript_to_cds_dict = {}     # transcript -> list of CDS exons
    transcript_cds_begin_dict = {}  # transcript -> first exon of the CDS
    gene_cds_begin_dict = {}        # gene -> list of first CDS exons

    file_type = ann_path.split('.')[-1]
    chromesome_set = set()
    # collect information from annotation file
    for line in open(ann_path, 'r'):
        if line[0] == '#':
            continue
        item = line.strip().split('\t')
        chromesome_set.add(item[0])
        feature_type = item[2]
        attribute_item = item[-1]
        attribute_dict = attribute_item_to_dict(attribute_item, file_type, feature_type)
        # store relationship between gene ID and its transcript IDs
        if feature_type in ['transcript', 'mRNA']:
            gene_id = attribute_dict['gene_id']
            transcript_id = attribute_dict['transcript_id']
            if attribute_dict['gene_type'] != 'protein_coding' or attribute_dict['transcript_type']  != 'protein_coding':
                continue
            assert (transcript_id not in transcript_to_gene_dict)
            transcript_to_gene_dict[transcript_id] = gene_id
            if gene_id in gene_to_transcript_dict and transcript_id not in gene_to_transcript_dict[gene_id]:
                gene_to_transcript_dict[gene_id].append(transcript_id)
            else:
                gene_to_transcript_dict[gene_id] = [transcript_id]
       # Todo python is 0-based while gene annotation file(.gtf, .vcf, .maf) is one based
        elif feature_type == "CDS":
            parent_ts = attribute_dict['transcript_id']
            strand_mode = item[6]
            cds_left = int(item[3])-1
            cds_right = int(item[4])
            frameshift = int(item[7])
            if parent_ts in transcript_to_cds_dict:
                transcript_to_cds_dict[parent_ts].append((cds_left, cds_right, frameshift))
            else:
                transcript_to_cds_dict[parent_ts] = [(cds_left, cds_right, frameshift)]
            if strand_mode == "+" :
                cds_start, cds_stop = cds_left, cds_right
            else:
                cds_start, cds_stop = cds_right, cds_left

            # we only consider the start of the whole CoDing Segment
            if parent_ts not in transcript_cds_begin_dict or                leq_strand(cds_start, transcript_cds_begin_dict[parent_ts][0], strand_mode):
                transcript_cds_begin_dict[parent_ts] = (cds_start, cds_stop, item)

    # collect first CDS exons for all transcripts of a gene
    for ts_key in transcript_to_gene_dict:

        target_gene = transcript_to_gene_dict[ts_key]
        if target_gene not in gene_cds_begin_dict:
            gene_cds_begin_dict[target_gene] = []
        if ts_key in transcript_cds_begin_dict:
            gene_cds_begin_dict[target_gene].append(transcript_cds_begin_dict[ts_key])

    # sort list of CDS exons per transcript
    for ts_key in transcript_to_cds_dict:
        transcript_to_cds_dict[ts_key] = sorted(transcript_to_cds_dict[ts_key], key=lambda coordpair: coordpair[0])

    return transcript_to_gene_dict, gene_to_transcript_dict, transcript_to_cds_dict, transcript_cds_begin_dict, gene_cds_begin_dict


# In[11]:


def attribute_item_to_dict(a_item, file_type, feature_type):
    """  From attribute item in annotation file to get corresponding dictionary

    Parameters
    ----------
    a_item: str. attribute item
    file_type: str. Choose from {'gtf', 'gff', 'gff3'}
    feature_type: str. Extract other fields. We only
        consider 'CDS', 'mRNA' and 'transcript'

    Returns
    -------
    gtf_dict: dict. store all the necessary data

    """
    gtf_dict = {}
    if file_type.lower() == 'gtf':
        attribute_list = a_item.split('; ')
        for attribute_pair in attribute_list:
            pair = attribute_pair.split(' ')
            gtf_dict[pair[0]] = pair[1][1:-1]
    elif file_type.lower() == 'gff3':
        attribute_list = a_item.split(';')
        for attribute_pair in attribute_list:
            pair = attribute_pair.split('=')
            gtf_dict[pair[0]] = pair[1]
    elif file_type.lower() == 'gff':
        gff_dict = {}
        attribute_list = a_item.split(';')
        for attribute_pair in attribute_list:
            pair = attribute_pair.split('=')
            gff_dict[pair[0]] = pair[1]  # delete "", currently now work on level 2
        if feature_type == 'CDS':
            gtf_dict['transcript_id'] = gff_dict['Parent']
        elif feature_type in {'mRNA', 'transcript'}:  # mRNA or transcript
            gtf_dict['gene_id'] = gff_dict['geneID']
            gtf_dict['transcript_id'] = gff_dict['ID']
            gtf_dict['gene_type'] = gff_dict['gene_type']
            gtf_dict['transcript_type'] = gff_dict['transcript_type']

    return gtf_dict


# In[12]:


def leq_strand(coord1, coord2, strand):
    if strand == "+":
        return coord1 <= coord2
    else:
        return coord1 >= coord2


# In[13]:


transcript_to_gene_dict30, gene_to_transcript_dict30, transcript_to_cds_dict30, transcript_cds_begin_dict30, gene_cds_begin_dict30 = process_anno(v30)


# In[14]:


transcript_to_gene_dict32, gene_to_transcript_dict32, transcript_to_cds_dict32, transcript_cds_begin_dict32, gene_cds_begin_dict32 = process_anno(v32)


# In[22]:


def compare_(input_a, input_b, label_a, label_b, base_name=False, print_diff=False):
    if base_name:
        input_a_wo_version = [ name.split('.')[0] for name in input_a]
        input_b_wo_version = [ name.split('.')[0] for name in input_b]
        a = set(input_a_wo_version)
        b = set(input_b_wo_version)
    
    else: 
        a = set(input_a)
        b = set(input_b)

    print(label_a)
    print(len(a))
    print('\n')
    print(label_b)
    print(len(b))
    print('\n')
    print("{} - {}".format(label_a, label_b))
    print(len(a.difference(b)))
    if print_diff:
        print(a.difference(b))
    print('\n')
    print("{} - {}".format(label_b, label_a))
    print(len(b.difference(a)))
    if print_diff:
        print(b.difference(a))
    print('\n')
    print("{} and {}".format(label_b, label_a))
    print(len(a.intersection(b)))
    return a.intersection(b)


# ## Compare

# In[34]:


print("transcript_to_gene_dict")
print(len(transcript_to_gene_dict30), len(transcript_to_gene_dict32))
print("gene_to_transcript_dict")
print(len(gene_to_transcript_dict30), len(gene_to_transcript_dict32))
print("transcript_to_cds_dict")
print(len(transcript_to_cds_dict30), len(transcript_to_cds_dict32))
print("transcript_cds_begin_dict")
print(len(transcript_cds_begin_dict30), len(transcript_cds_begin_dict32))
print("gene_cds_begin_dict")
print(len(gene_cds_begin_dict30), len(gene_cds_begin_dict32))


# In[53]:


# input_a = gene_cds_begin_dict30.keys()
# input_b = gene_cds_begin_dict32.keys()
# label_a = "v30"
# label_b = "v32"
# intersect = compare_ (input_a, input_b, label_a, label_b, base_name=True, print_diff=True)


# In[54]:


# input_a = gene_to_transcript_dict30.keys()
# input_b = gene_to_transcript_dict32.keys()
# label_a = "v30"
# label_b = "v32"
# intersect = compare_ (input_a, input_b, label_a, label_b)


# In[55]:


# input_a = gene_to_transcript_dict30.keys()
# input_b = gene_cds_begin_dict30.keys()
# label_a = "gene_to_transcript_dict30"
# label_b = "gene_cds_begin_dict30"
# intersect = compare_ (input_a, input_b, label_a, label_b)


# In[56]:


# input_a = transcript_to_gene_dict30.keys()
# input_b = transcript_cds_begin_dict30.keys()
# label_a = "transcript_to_gene_dict30"
# label_b = "transcript_cds_begin_dict30"
# intersect = compare_ (input_a, input_b, label_a, label_b)


# ## Test the different names 

# ## Save 

# In[40]:


input_a = gene_cds_begin_dict30.keys()
input_b = gene_cds_begin_dict32.keys()
label_a = "v30"
label_b = "v32"
intersect_wo_version = compare_(input_a, input_b, label_a, label_b, base_name=True, print_diff=False)
df = pd.DataFrame(intersect_wo_version)
df.to_csv(os.path.join(outpath, 'genes_coding_gencode_v32_inter_v30_wo_version.txt'),                      index = None, header = False)

gene_to_transcript_dict30_noversion = {}
gene_to_transcript_dict32_noversion = {}

for key, item in gene_to_transcript_dict30.items():
    gene_to_transcript_dict30_noversion[key.split('.')[0]] = item
for key, item in gene_to_transcript_dict32.items():
    gene_to_transcript_dict32_noversion[key.split('.')[0]] = item
gtt30 = [gene_to_transcript_dict30_noversion[gene] for gene in intersect_wo_version]
gtt32 = [gene_to_transcript_dict32_noversion[gene] for gene in intersect_wo_version]
gtt30 = [item for sublist in gtt30 for item in sublist]
gtt32 = [item for sublist in gtt32 for item in sublist]
print(len(gtt30))
print(len(gtt32))
df = pd.DataFrame(gtt32)
df.to_csv(os.path.join(outpath, 'transcripts_withcds_gencode_v32_inter_v30_wo_version.txt'),                      index = None, header = False)

df = pd.DataFrame(gene_cds_begin_dict30.keys())
df.to_csv(os.path.join(outpath, 'genes_coding_gencode_v30.txt'),                      index = None, header = False)


# In[86]:


df = pd.DataFrame(gene_cds_begin_dict30.keys())
df.to_csv(os.path.join(outpath, 'genes_coding_gencode_v30.txt'),                      index = None, header = False)


# In[87]:


df = pd.DataFrame(gene_cds_begin_dict32.keys())
df.to_csv(os.path.join(outpath, 'genes_coding_gencode_v32.txt'),                      index = None, header = False)


# In[89]:


df = pd.DataFrame(transcript_to_cds_dict30.keys())
df.to_csv(os.path.join(outpath, 'transcript_cds_gencode_v30.txt'),                      index = None, header = False)


# In[90]:


df = pd.DataFrame(transcript_to_cds_dict32.keys())
df.to_csv(os.path.join(outpath, 'transcript_cds_gencode_v32.txt'),                      index = None, header = False)


# In[93]:


input_a = gene_cds_begin_dict30.keys()
input_b = gene_cds_begin_dict32.keys()
label_a = "v30"
label_b = "v32"
gene_intersect = compare_ (input_a, input_b, label_a, label_b)

df = pd.DataFrame(gene_intersect)
df.to_csv(os.path.join(outpath, 'genes_coding_gencode_v32_inter_v30.txt'),                      index = None, header = False)

gtt30 = [gene_to_transcript_dict30[gene] for gene in gene_intersect]
gtt32 = [gene_to_transcript_dict32[gene] for gene in gene_intersect]
gtt30 = [item for sublist in gtt30 for item in sublist]
gtt32 = [item for sublist in gtt32 for item in sublist]
print(len(gtt30))
print(len(gtt32))
df = pd.DataFrame(gtt32)
df.to_csv(os.path.join(outpath, 'transcripts_withcds_gencode_v32_inter_v30.txt'),                      index = None, header = False)


# In[ ]:


intersect for annotation cds and reading frame different 


# In[63]:


# Delta genes due to version taken into account for immunopepper

print("intersection without version tag")
input_a = gene_cds_begin_dict30.keys()
input_b = gene_cds_begin_dict32.keys()
label_a = "v30"
label_b = "v32"
intersect_wo_version = compare_(input_a, input_b, label_a, label_b, base_name=True, print_diff=False)
df = pd.DataFrame(intersect_wo_version)

print("intersection with version tag")
input_a = gene_cds_begin_dict30.keys()
input_b = gene_cds_begin_dict32.keys()
label_a = "v30"
label_b = "v32"
intersect_with_version = compare_(input_a, input_b, label_a, label_b, base_name=False, print_diff=False)
df = pd.DataFrame(intersect_with_version)

core_intersect_with_version  = {}
for gene in intersect_with_version:
    core_intersect_with_version[gene] = gene.split('.')[0] 
    core_intersect_with_version_g = set(core_intersect_with_version.values())
print('missing genes due to version naming')
missing_genes = intersect_wo_version.difference(core_intersect_with_version_g)
print(len(missing_genes))
df = pd.DataFrame(missing_genes)
df.to_csv(os.path.join(outpath, 'missing_genes_untagged_from_gencode_v32_inter_v30_with_version.txt'),                      index = None, header = False)

print('missing genes due to version naming, mapped back to v32')
missing_genes_with_32_tag = []
for gene in gene_cds_begin_dict32.keys():
    if gene.split('.')[0] in missing_genes:
        missing_genes_with_32_tag.append(gene)
print(len(missing_genes_with_32_tag))
df = pd.DataFrame(missing_genes_with_32_tag)
df.to_csv(os.path.join(outpath, 'missing_genes_taggedv32_from_gencode_v32_inter_v30_with_version.txt'),                      index = None, header = False)


print('missing genes due to version naming, mapped back to v30')
missing_genes_with_30_tag = []
for gene in gene_cds_begin_dict30.keys():
    if gene.split('.')[0] in missing_genes:
        missing_genes_with_30_tag.append(gene)
print(len(missing_genes_with_30_tag))
df = pd.DataFrame(missing_genes_with_30_tag)
df.to_csv(os.path.join(outpath, 'missing_genes_taggedv30_from_gencode_v32_inter_v30_with_version.txt'),                      index = None, header = False)


# ## Prototype get the intersection annotations

# In[57]:


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

        
        
        


# In[ ]:


fp.close()


# In[ ]:





