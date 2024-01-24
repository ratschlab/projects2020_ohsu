import os 
import pandas as pd
from collections import defaultdict
import numpy as np 
import glob
import gzip
import pickle

from immunopepper.preprocess import attribute_item_to_dict
from immunopepper.preprocess import leq_strand


def preprocess_annot_custom(ann_path):
    # Partial copy From Immunopepper preprocess.py 
    transcript_to_gene_dict = {}    # transcript -> gene id


    gene_to_transcript_dict = {}    # gene_id -> list of transcripts
    gene_cds_begin_dict = {}        # gene -> list of first CDS exons

    transcript_to_cds_dict = {}     # transcript -> list of CDS exons
    transcript_cds_begin_dict = {}  # transcript -> first exon of the CDS
    transcript_to_strand = {}

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
            transcript_to_strand[parent_ts] = strand_mode
            if parent_ts in transcript_to_cds_dict:
                transcript_to_cds_dict[parent_ts].append((cds_left, cds_right, frameshift))
            else:
                transcript_to_cds_dict[parent_ts] = [(cds_left, cds_right, frameshift)]
            if strand_mode == "+" :
                cds_start, cds_stop = cds_left, cds_right
            else:
                cds_start, cds_stop = cds_right, cds_left

            # we only consider the start of the whole CoDing Segment
            if parent_ts not in transcript_cds_begin_dict or \
               leq_strand(cds_start, transcript_cds_begin_dict[parent_ts][0], strand_mode):
                transcript_cds_begin_dict[parent_ts] = (cds_start, cds_stop, item)

    # collect first CDS exons for all transcripts of a gene
    for ts_key in transcript_to_gene_dict:
        target_gene = transcript_to_gene_dict[ts_key]
        if target_gene not in gene_cds_begin_dict:
            gene_cds_begin_dict[target_gene] = []
        if ts_key in transcript_cds_begin_dict:
            gene_cds_begin_dict[target_gene].append(transcript_cds_begin_dict[ts_key])


    return transcript_to_gene_dict, gene_to_transcript_dict, \
    gene_cds_begin_dict, transcript_to_cds_dict, transcript_cds_begin_dict, transcript_to_strand
    

def my_CDS_collection(transcript_to_gene_dict, gene_cds_begin_dict, transcript_to_cds_dict, 
                      transcript_cds_begin_dict, transcript_to_strand):
    # Custom collection of CDS 
    transcript_cds_begin_dict_bis = {}
    transcript_cds_end_dict_bis = {}

    gene_cds_begin_dict_bis = defaultdict(list)
    gene_cds_end_dict_bis = defaultdict(list)

    # will be in reading order 
    for ts_key in transcript_to_cds_dict:
        if transcript_to_strand[ts_key] == '+': # '+'
            transcript_cds_begin_dict_bis[ts_key] = (transcript_to_cds_dict[ts_key][0][0],
                                                     transcript_to_cds_dict[ts_key][0][1], '+')
            transcript_cds_end_dict_bis[ts_key] = (transcript_to_cds_dict[ts_key][-1][0],
                                                     transcript_to_cds_dict[ts_key][-1][1], '+')


        else: 
            transcript_cds_begin_dict_bis[ts_key] = (transcript_to_cds_dict[ts_key][0][1],
                                                     transcript_to_cds_dict[ts_key][0][0], '-')
            transcript_cds_end_dict_bis[ts_key] = (transcript_to_cds_dict[ts_key][-1][1],
                                                     transcript_to_cds_dict[ts_key][-1][0], '-')

        assert(transcript_cds_begin_dict_bis[ts_key][0] == transcript_cds_begin_dict[ts_key][0])
        assert(transcript_cds_begin_dict_bis[ts_key][1] == transcript_cds_begin_dict[ts_key][1])

    # collect first, last CDS exons for all transcripts of a gene
    for ts_key in transcript_to_gene_dict:
        target_gene = transcript_to_gene_dict[ts_key]
        gene_cds_begin_dict_bis[target_gene].append(transcript_cds_begin_dict_bis[ts_key])
        gene_cds_end_dict_bis[target_gene].append(transcript_cds_end_dict_bis[ts_key])
    
    return gene_cds_begin_dict_bis, gene_cds_end_dict_bis
     


def extract_end_starts(pep_orig_coord, strand):
    ''' Get peptide end and start coordinates'''
    if strand == '+': # Do - strand 
        pep_start = int(pep_orig_coord[0])
        pep_end = int(pep_orig_coord[-1])
    else: 
        pep_start = int(pep_orig_coord[1])
        pep_end = int(pep_orig_coord[-2])
    return pep_start, pep_end


def get_include_flag(start_cds, end_cds, pep_start, pep_end, has_stop_codon ):
    '''Use end and start coordinates for 3' 5' include flag'''
    if pep_start in start_cds: # We will always miss things that are new in the graph 
        pep_5include = 1
    else: 
        pep_5include = 0 
    if (pep_end in end_cds) or (has_stop_codon) == '1':
        pep_3include = 1
    else: 
        pep_3include = 0 
    return pep_5include, pep_3include


def get_nt_len_with_aa_shift(pep_modi_coord):
    '''Get nt length of each exon involved -> jx_list, shift_list'''
    tot_len = 0 
    shift = 0 
    jx_list = []
    jx_list_ori = []
    shift_list = []
    for pair in np.arange(0, len(pep_modi_coord), 2):
        cds = int(pep_modi_coord[pair + 1]) - int(pep_modi_coord[pair])  # 0 based, open right 
        jx_list_ori.append(cds)
        cds += shift 
        shift = cds % 3
        jx_list.append(cds - shift)
        shift_list.append(shift)
        
    return jx_list, shift_list, jx_list_ori


def get_aaPos_betweenFlag(shift_list, jx_list):
    '''Get aa position of the junction
    the junction coordinate jx_pos is the 0-based position in the peptide 
    of the amino acid that either overlaps the junction (if the junction is 
    in the middle of a codon), or is immediately before it if the junction 
    occurs between codons'''
    if shift_list[0]: # junction is inside an amino acid
        aa_junction_pos0 = int((jx_list[0] / 3)) # because 0 based
        between_codons0 = 0 
    else: # junction is between amino acids 
        aa_junction_pos0 = int((jx_list[0] / 3) - 1)  # because 0 based
        between_codons0 = 1
        
    if len(shift_list) > 2: #third exon 
        if shift_list[1]: # junction is inside an amino acid
            aa_junction_pos1 = int((jx_list[1] / 3)) # because 0 based
            between_codons1 = 0 
        else: # junction is between amino acids 
            aa_junction_pos1 = int((jx_list[1] / 3) - 1)  # because 0 based
            between_codons1 = 1 
        aa_junction_pos1_from_start = aa_junction_pos1 + aa_junction_pos0 + 1 
    else:
        aa_junction_pos1 = None
        between_codons1 = None
        aa_junction_pos1_from_start = None
    
    return aa_junction_pos0, between_codons0, aa_junction_pos1, between_codons1, \
           aa_junction_pos1_from_start


def get_genomic_coordinates(pep_modi_coord, strand):
    '''We have in + case: exon1_start, exon 1_stop, exon2_start, exon2_stop, exon3_start, exon3_stop
     In the - case: exon1_stop, exon 1_start, exon 2_stop, exon2_start, exon3_stop, exon3_start'''
    genome_junction_pos1 = None
    if strand == '+':
        genome_junction_pos0 = '{}_{}'.format(pep_modi_coord[1], pep_modi_coord[2])
        if len(pep_modi_coord) > 4:
            genome_junction_pos1 = '{}_{}'.format(pep_modi_coord[3], pep_modi_coord[4])
    else:
        genome_junction_pos0 = '{}_{}'.format(pep_modi_coord[0], pep_modi_coord[3])
        if len(pep_modi_coord) > 4:
            genome_junction_pos1 = '{}_{}'.format(pep_modi_coord[2], pep_modi_coord[5])
    return genome_junction_pos0, genome_junction_pos1


def split_coord(pep_coord):
    pep_coord = pep_coord.split(';')
    pep_coord = [coord for coord in pep_coord if (coord != 'None') and (coord != 'nan')]
    return pep_coord


def write_peptide_to_experiment(filepointer, pep_idx=None, pep_seq=None,\
                                idx=None, header=False):
    if header:
        header_exp = 'peptide_id\tpeptide_sequence\texperiment_ids\n'
        filepointer.write(header_exp)
    elif pep_idx is not None:
        exp_line = '{}\t{}\t{}\n'.format(pep_idx, 
                                         pep_seq,
                                         idx)
        filepointer.write(exp_line)

            
def write_fasta(write_, filepointer, pep_seq, pep_idx, aa_junction_pos, 
                aa_junction_pos1_from_start, between_codons, between_codons1,
                pep_5include, pep_3include, pep_gene, 
                genome_junction_pos, genome_junction_pos1, 
                kmer, jx_pep1, jx_pep2, readFrameAnnotated, \
                kmer_coord, kmer_type, strand, do_write=True):
        
    if write_:
        pep_handle1 = (f'>pepID-{pep_idx};jx_pos-{aa_junction_pos};between_codons-{between_codons}'
                       f';includes_5\'-{pep_5include};includes_3\'-{pep_3include};gene-{pep_gene};'
                       f'jx_coord-{genome_junction_pos};kmer-{kmer};readFrameAnnotated-{readFrameAnnotated};'
                       f'kmer_coord-{kmer_coord};origin-{kmer_type};strand{strand}')

        pep_handle2 = (f'>pepID-{pep_idx};jx_pos-{aa_junction_pos1_from_start};between_codons-{between_codons1}'
                       f';includes_5\'-{pep_5include};includes_3\'-{pep_3include};gene-{pep_gene};'
                       f'jx_coord-{genome_junction_pos1};kmer-{kmer};readFrameAnnotated-{readFrameAnnotated};'
                       f'kmer_coord-{kmer_coord};origin-{kmer_type};strand{strand}')

        if kmer in jx_pep1:
            pep_idx+=1
            filepointer.write(pep_handle1 + '\n')
            filepointer.write(pep_seq + '\n')
        elif kmer in jx_pep2:
            pep_idx+=1
            filepointer.write(pep_handle2 + '\n')
            filepointer.write(pep_seq + '\n')
        #### 01/2024 ADDITION
        elif kmer in pep_seq:
            filepointer.write(pep_handle1 + '\n')
            filepointer.write(pep_seq + '\n')
            filepointer.write(pep_handle2 + '\n')
            filepointer.write(pep_seq + '\n')
        #### END 01/2024 ADDDITION
        return pep_idx


def cut_peptides(pep_seq, jx_list, between_codons, between_codons1, aa_junction_pos, 
                 aa_junction_pos1, aa_junction_pos1_from_start, 
                 print_ = False):
    peptide_cut = []
    aa_junction_pos_shift = aa_junction_pos + 1 
    exon1, aa_jx1, exon2, aa_jx2, exon3 = '', '', '', '', ''

    if len(jx_list) == 2:
        if between_codons:
            exon1 = pep_seq[:aa_junction_pos_shift]
            exon2 = pep_seq[aa_junction_pos_shift:]
        else:
            exon1 = pep_seq[:aa_junction_pos]
            aa_jx1 = pep_seq[aa_junction_pos:aa_junction_pos_shift]
            exon2 = pep_seq[aa_junction_pos_shift:]
    elif len(jx_list) == 3:
        aa_junction2_pos_shift =  aa_junction_pos1 + aa_junction_pos + 1
        aa_junction2_pos_sshift = aa_junction_pos1 + aa_junction_pos + 2
        assert(aa_junction_pos1_from_start == aa_junction2_pos_shift)
        if between_codons and between_codons1: 
            exon1 = pep_seq[:aa_junction_pos_shift]
            exon2 = pep_seq[aa_junction_pos_shift:
                          aa_junction2_pos_sshift]
            exon3 = pep_seq[aa_junction2_pos_sshift:]
        if (not between_codons) and between_codons1: 
            exon1 = pep_seq[:aa_junction_pos ]
            aa_jx1 = pep_seq[aa_junction_pos:aa_junction_pos_shift]
            exon2 = pep_seq[aa_junction_pos_shift:
                          aa_junction2_pos_sshift]
            exon3 = pep_seq[aa_junction2_pos_sshift:]
        if (between_codons) and (not between_codons1): 
            exon1 = pep_seq[:aa_junction_pos_shift]
            exon2 = pep_seq[aa_junction_pos_shift :
                          aa_junction2_pos_shift]
            aa_jx2 = pep_seq[aa_junction2_pos_shift:
                         aa_junction2_pos_sshift]
            exon3 = pep_seq[aa_junction2_pos_sshift:]
        if (not between_codons) and (not between_codons1): 
            exon1 = pep_seq[:aa_junction_pos ]
            aa_jx1 = pep_seq[aa_junction_pos:aa_junction_pos_shift]
            exon2 = pep_seq[aa_junction_pos_shift :
                          aa_junction2_pos_shift]
            aa_jx2 = pep_seq[aa_junction2_pos_shift:
                         aa_junction2_pos_sshift]
            exon3 = pep_seq[aa_junction2_pos_sshift:]
    if print_:
        print(f'exon1:{exon1}, aa_containing_jx1:{aa_jx1}, exon2:{exon2}, aa_containing_jx2:{aa_jx2}, exon3:{exon3}')
        print(f'junction positions jx1: {aa_junction_pos}, jx2:{aa_junction_pos1_from_start}')
        print(f'is junction between a codon jx1: {between_codons}, jx2: {between_codons1}')
        print('\n')
    return exon1 + aa_jx1 + exon2, exon2 + aa_jx2 + exon3 


def print_stats(print_, kmer, pep_seq, strand, pep_orig_coord, pep_modi_coord, jx_list,
               jx_list_ori, genome_junction_pos, genome_junction_pos1,
               aa_junction_pos, aa_junction_pos1, between_codons, between_codons1):
    
    if print_:
        p_ori_coord = ';'.join(pep_orig_coord)
        p_modif_coord = ';'.join(pep_modi_coord)
        print(f'INSTANCE: \n kmer {kmer}/ sequence {pep_seq}/ strand {strand} / \n original coordinates {p_ori_coord} / \n modif coordinates {p_modif_coord} /  \n junction list origin {jx_list_ori}/ junction list {jx_list} / \n junction coordinates 1 {genome_junction_pos} / junction coordinates 2 {genome_junction_pos1}')
        #print(aa_junction_pos, between_codons, aa_junction_pos1, between_codons1)
        print('peptide length', len(pep_seq))


def readlines_custom(file_=None, pandas_df=None):
    if file_ is not None:
        with gzip.open(file_, 'rb') as fp: 
            lines = fp.readlines()
            for j, line in enumerate(lines):
                line = line.decode().replace('\n', '').split('\t')
                if '\\t' in line[0]:
                    splitted_p = line[0].split('\\t') #correction temp
                    kmer, coord, peptide = splitted_p
                    line = [kmer, coord, peptide] + line[1:]
                else:
                    line = [kmer, coord] + line #Use previous kmer instance
                lines[j] = line
        return lines
    elif pandas_df is not None:
        lines = pandas_df.iterrows()
        return lines