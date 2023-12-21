import os 
import pandas as pd
from collections import defaultdict
import numpy as np 
import glob


def extract_peptide_fields_pq(pep):
    '''0 peptide
    1 id
    2 readFrame
    3 readFrameAnnotated
    4 geneName
    5 geneChr
    6 geneStrand
    7 mutationMode
    8 junctionAnnotated
    9 hasStopCodon
    10 isInJunctionList
    11 isIsolated
    12 variantComb
    13 variantSegExpr
    14 modifiedExonsCoord
    15 originalExonsCoord
    16 vertexIdx
    ?junctionExpr # sample mode only
    ?segmentExpr #sample mode only
    17 kmerType 
    18 minorIntron'''
    
    pep_seq = pep['peptide']
    pep_gene = pep['geneName']
    pep_orig_coordS = pep['originalExonsCoord'].split('/')
    pep_modif_coordS = pep['modifiedExonsCoord'].split('/')
    strand = pep['geneStrand']
    has_stop_codon = pep['hasStopCodon']
    readFrameAnnotated = pep['readFrameAnnotated']
    junctionAnnotated = pep['junctionAnnotated'].split(';')
    if len(junctionAnnotated) == 1:
        junctionAnnotated.append(np.nan)
    kmer_type = pep['kmerType'].replace('-','').replace('\n','')
    is_isolated = pep['isIsolated']
    junctionExpr = pep['junctionExpr']
    minorIntron = pep['minorIntron']

    return pep, pep_seq, pep_gene, pep_orig_coordS, pep_modif_coordS, \
           strand, has_stop_codon, readFrameAnnotated, junctionAnnotated, \
           kmer_type, is_isolated, junctionExpr, minorIntron


def extract_end_starts(pep_orig_coord, strand):
    ''' Get peptide end and start coordinates'''
    if strand == '+': # Do - strand 
        pep_start = np.int(pep_orig_coord[0])
        pep_end = np.int(pep_orig_coord[-1])
    else: 
        pep_start = np.int(pep_orig_coord[1])
        pep_end = np.int(pep_orig_coord[-2])
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

def preprocess_line(line):
    assert(False)
    line = line.replace('3-exons_9-mer ', '3-exons_9-mer@').replace('2-exons ', '2-exons@')
    kmer = line.split(',')[0]
    peptides = ','.join(line.split(',')[1:])
    peptides = peptides.split('@')
    return line, kmer, peptides

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


def write_fasta_option_MI(sp, pep_seq, pep_idx, aa_junction_pos, 
                aa_junction_pos1_from_start, between_codons, between_codons1,
                pep_5include, pep_3include, pep_gene, gene_strand,
                genome_junction_pos, genome_junction_pos1, 
                kmer, jx_pep1, jx_pep2, readFrameAnnotated, \
                junctionAnnotated, kmer_type, minorIntron, junctionExpr, sample, 
                jx_to_expression, do_write=True):
    
        
    expr_val = junctionExpr.split(';')
    
    pep_handle1 = '>pepID-{};jx_pos-{};between_codons-{};includes_5\'-{};includes_3\'-{};gene-{};jx_coord-{};kmer-{};readFrameAnnotated-{};junctionAnnotated-{};origin-{}'.format(
    pep_idx, aa_junction_pos, between_codons, pep_5include, 
    pep_3include, pep_gene, genome_junction_pos, kmer, readFrameAnnotated, 
    junctionAnnotated[0], kmer_type)
    expr1 = float(expr_val[0])
    
    if aa_junction_pos1_from_start:
        pep_handle2 = '>pepID-{};jx_pos-{};between_codons-{};includes_5\'-{};includes_3\'-{};gene-{};jx_coord-{};kmer-{};readFrameAnnotated-{};junctionAnnotated-{};origin-{}'.format(
        pep_idx, aa_junction_pos1_from_start, between_codons1, pep_5include, 
        pep_3include, pep_gene, genome_junction_pos1, kmer, readFrameAnnotated, 
        junctionAnnotated[1], kmer_type)
        expr2 = float(expr_val[1])
    else:
        pep_handle2 = pep_handle1
        expr2 = expr1
        genome_junction_pos1 = genome_junction_pos
        
    # Whether to write the first, the second junction or both 
    if (minorIntron == 1) or (minorIntron == 3) :

        if gene_strand == '+':
            #write_first 
            if expr1 > 0: 
                pep_idx+=1
                sp.write(pep_handle1 + '\n')
                sp.write(pep_seq + '\n')
                jx_to_expression[0].append(genome_junction_pos)
                jx_to_expression[3].append(expr1)
                jx_to_expression[1].append(pep_seq)
                jx_to_expression[2].append(sample)
        else:
            #write_second 
            if expr2 > 0: 
                pep_idx+=1
                sp.write(pep_handle2 + '\n')
                sp.write(pep_seq + '\n')
                jx_to_expression[0].append(genome_junction_pos1)
                jx_to_expression[3].append(expr2)
                jx_to_expression[1].append(pep_seq)
                jx_to_expression[2].append(sample)
    if (minorIntron == 2) or (minorIntron == 3) :

        if gene_strand == '+':
            #write_second 
            if expr2 > 0:
                pep_idx+=1
                sp.write(pep_handle2 + '\n')
                sp.write(pep_seq + '\n')
                jx_to_expression[0].append(genome_junction_pos1)
                jx_to_expression[3].append(expr2)
                jx_to_expression[1].append(pep_seq)
                jx_to_expression[2].append(sample)
        else:
            #write_first 
            if expr1 > 0: 
                pep_idx+=1
                sp.write(pep_handle1 + '\n')
                sp.write(pep_seq + '\n')
                jx_to_expression[0].append(genome_junction_pos)
                jx_to_expression[3].append(expr1)
                jx_to_expression[1].append(pep_seq)
                jx_to_expression[2].append(sample)
    return pep_idx, jx_to_expression



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
