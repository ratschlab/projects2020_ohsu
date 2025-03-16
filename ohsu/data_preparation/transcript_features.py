from collections import deque
from copy import deepcopy
from math import floor
import os
import subprocess as sp
import sys
import re
sys.path.append(
    os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
)
from classes import transcript_states as ts
from classes import insert_junction as ij

# Note: * == stop
_TRANSLATION_DICT = {
    "GCT": "A", "GCG": "A", "GCA": "A", "GCC": "A", "GGT": "G", "GGC": "G",
    "GGA": "G", "GGG": "G", "ATT": "I", "ATC": "I", "ATA": "I", "CTT": "L",
    "CTC": "L", "CTA": "L", "CTG": "L", "TTA": "L", "TTG": "L", "CCT": "P",
    "CCC": "P", "CCA": "P", "CCG": "P", "GTT": "V", "GTC": "V", "GTA": "V",
    "GTG": "V", "TTT": "F", "TTC": "F", "TGG": "W", "TAT": "Y", "TAC": "Y",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R", "CAT": "H", "CAC": "H",
    "AAA": "K", "AAG": "K", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "AGT": "S", "AGC": "S", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "TGT": "C", "TGC": "C", "ATG": "M", "AAT": "N", "AAC": "N", "CAA": "Q",
    "CAG": "Q", "TAA": "*", "TAG": "*", "TGA": "*"
}
_SELENO_TRANSLATION_DICT = deepcopy(_TRANSLATION_DICT)
_SELENO_TRANSLATION_DICT['TGA'] = "U"
_MITOCHONDRIAL_TRANSLATION_DICT = {
    "GCT": "A", "GCG": "A", "GCA": "A", "GCC": "A", "GGT": "G", "GGC": "G",
    "GGA": "G", "GGG": "G", "ATT": "I", "ATC": "I", "ATA": "M", "CTT": "L",
    "CTC": "L", "CTA": "L", "CTG": "L", "TTA": "L", "TTG": "L", "CCT": "P",
    "CCC": "P", "CCA": "P", "CCG": "P", "GTT": "V", "GTC": "V", "GTA": "V",
    "GTG": "V", "TTT": "F", "TTC": "F", "TGG": "W", "TAT": "Y", "TAC": "Y",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGA": "*", "AGG": "*", "CAT": "H", "CAC": "H",
    "AAA": "K", "AAG": "K", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "AGT": "S", "AGC": "S", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "TGT": "C", "TGC": "C", "ATG": "M", "AAT": "N", "AAC": "N", "CAA": "Q",
    "CAG": "Q", "TAA": "*", "TAG": "*", "TGA": "W"
}
_SUPPORTED_NUCLEOTIDES = ['G', 'C', 'A', 'T']
_REV_COMP = str.maketrans("ATCG", "TAGC")
_MOCK_LEFT = '1'
_MOCK_RIGHT = '2'
_STRAND_MAPPER = {'+': 'PLUS', '-': 'MINUS'}
_TARGET_LENGTH = 153


def jx_modified_tx_name(orig_tx_id, jx):
    chrom, left, right, strand = jx.split(';')
    new_jx = '.'.join([chrom.upper(), left, right, _STRAND_MAPPER[strand]])
    return orig_tx_id + '.MOD.' + new_jx


def get_trans_type(gtf_info):
    return re.sub(r".*transcript_type \"(\S+?)\"[;].*", r"\1", gtf_info)


def get_trans_id(gtf_info):
    return re.sub(r".*transcript_id \"(\S+?)\"[;].*", r"\1", gtf_info)


def get_gene_name(gtf_info):
    return re.sub(r".*gene_name \"(\S+?)\"[;].*", r"\1", gtf_info)


def get_gene_id(gtf_info):
    return re.sub(r".*gene_id \"(\S+?)\"[;].*", r"\1", gtf_info)


def cache_reference_fasta_output(curr_range, samtools_path, reference_fasta):
    try:
        output = sp.check_output([
            '{}'.format(samtools_path), 'faidx',
            '{}'.format(reference_fasta), '{}'.format(curr_range)
        ])
    except sp.CalledProcessError:
        output = None
    else:
        Transcript.FASTA_MAP[curr_range] = output
    return output


def get_dna_seq(chrom, range_list, reference_fasta, samtools_path, strand):
    """Given an input junction and reference genome, transcribes RNA sequence.

    Input:
    chrom: junction chromosome (string)
    left: junction global left coordinate, 0-based, closed (string)
    right: junction global right coordinate, 0-based, closed (string)
    reference_fasta: reference genome fasta file, previously sorted and
        indexed by samtools (path to fasta file, string)
    samtools_path: to be called from subprocess to collect the sequence (path
        to samtools executable, string)
    protein_length: desired minimum length of returned protein sequence (int)
    print_output: whether to print DNA sequence to std out (bool)
    dna_len: to override desired minimum protein length with an exact length
        for desired DNA sequence returned (int)
    split_seq: to return two sequences, from the left and right sides of the
        junction, separately instead of a joined RNA tx sequence (bool)
    strand: input required if returning split sequences, in the direction of
        the tx instead of global chromosome direction (string)

    Returns a left-to-right nucleotide sequence on either side of the aberrant
    junction sufficiently long to generate the desired protein sequence.
    """
    full_seq = ''
    for left, right in range_list:
        curr_range = chrom + ':' + str(left) + '-' + str(right)
        output = Transcript.FASTA_MAP.get(
            curr_range, cache_reference_fasta_output(
                curr_range, samtools_path, reference_fasta
            )
        )
        if output is None:
            return ''
        seq = ''.join(output.decode("utf-8").splitlines()[1:])
        if strand == '-':
            full_seq = seq + full_seq
        else:
            full_seq += seq
    return full_seq


def codon_to_aa(codon, translation_dict):
    """Checks a codon for acceptable characters, and returns an amino acid.

    Input:
        codon: three-nucleotide amino acid sequence (string)

    Returns the appropriate amino acid, stop signal, or not-supported signal.
    """
    for char in codon:
        if char == '-':
            return '-'
        if char not in _SUPPORTED_NUCLEOTIDES:
            return 'X'
    return translation_dict[codon]


def seq_to_protein(dna_sequence, strand, minimum_length, translation_dict):
    """Translates an input RNA sequence to (an) amino acid sequence(s).

    Input:
        dna_sequence: string consisting of A, T, C, and G only (string)
        strand: '+' or '-' (string)
        pep_length: target maximum length of protein to return; also
            determines the minimum acceptable length of protein (int)
        print_output: whether to print translated sequence to stdout (boolean)

    Determines the minimum acceptable length (0.75 of target length) of the
    returned protein.  Translates each codon and appends to the protein in the
    appropriate reading frame.

    Returns a list of three strings; each is a protein sequence for one of
    three possible reading frames, if not eliminated due to stop codons.
    """
    stop_seq = [False, False, False]
    amino_acid_seq = ['' for _ in range(3)]
    orig_seqs = ['' for _ in range(3)]
    if strand == '-':
        dna_sequence = dna_sequence.translate(_REV_COMP)[::-1]
    for k, char in enumerate(dna_sequence[:-2]):
        window = dna_sequence[k:k + 3]
        frame = k % 3
        amino_acid = codon_to_aa(window, translation_dict)
        orig_seqs[frame] += amino_acid
        if amino_acid == '*':
            stop_seq[frame] = True
        else:
            if not stop_seq[frame]:
                amino_acid_seq[frame] += amino_acid
    amino_acid_seq = [
        seq if len(seq) >= minimum_length else '' for seq in amino_acid_seq
    ]
    return amino_acid_seq, stop_seq


class ExonNode(object):
    """Single transcript feature (exon)"""
    def __init__(self, left, right, strand, chrom, gtf_line_list):
        self.left = int(left)
        self.right = int(right)
        self.strand = strand
        self.chrom = chrom
        self.gtf_line_list = gtf_line_list
        self.has_start_codon = False
        self.start_codon_add = 2
        self.stop_codon_add = 2
        self.has_stop_codon = False
        self.cds_left = 0
        self.cds_right = 0
        self.has_exon_cds = False
        self.cds_rf = 0
        self.last_exon = False

    @property
    def cds_length(self):
        if self.has_exon_cds:
            left = max(self.left, self.cds_left)
            if self.cds_right:
                right = self.cds_right
            else:
                right = self.right
            return 1 + (right - left)
        else:
            return 0

    @property
    def start_codon_left(self):
        if self.has_start_codon:
            if self.strand == '-':
                return self.cds_right - self.start_codon_add
            else:
                return self.cds_left
        else:
            return 0

    @property
    def stop_codon_left(self):
        if self.has_stop_codon:
            if self.strand == '-':
                return self.cds_left - 1 - self.stop_codon_add
            else:
                return self.cds_right + 1
        else:
            return 0

    def to_string(self):
        start = ''
        stop = ''
        cds_left = str(self.cds_left)
        cds_right = str(self.cds_right)
        if self.has_start_codon:
            start = 'start codon'
        if self.has_stop_codon:
            stop = 'stop codon'
        cds_string = ', '.join([cds_left, cds_right, start, stop])
        le = '; last exon' if self.last_exon else ''
        return '{}:{}-{}, {}; cds: {}{}'.format(
            self.chrom, self.left, self.right, self.strand, cds_string, le
        )

    def add_cds(self, left, right, first_cds):
        if first_cds:
            if self.strand == '-':
                self.cds_right = right
            else:
                self.cds_left = left
        else:
            if left != self.left:
                self.cds_left = left
                if left != 0:
                    if self.strand == '-':
                        self.last_exon = True
            if right != self.right:
                self.cds_right = right
                if right != 0:
                    if self.strand == '+':
                        self.last_exon = True
        self.has_exon_cds = True
        return None

    def add_start_codon(self, start_codon_left_pos, right_pos):
        self.has_start_codon = True
        self.start_codon_add = right_pos - start_codon_left_pos
        if self.strand == '-':
            if not self.cds_right:
                self.cds_right = right_pos
        else:
            if not self.cds_left:
                self.cds_left = start_codon_left_pos
        self.has_exon_cds = True
        return None

    def add_stop_codon(self, stop_codon_left_pos, right_pos):
        prev_exon_last = False
        self.has_stop_codon = True
        self.stop_codon_add = right_pos - stop_codon_left_pos
        if self.strand == '-':
            if not self.cds_left:
                if right_pos != self.right:
                    self.cds_left = right_pos + 1
                    self.has_exon_cds = True
                    self.last_exon = True
                else:
                    prev_exon_last = True
        else:
            if not self.cds_right:
                if stop_codon_left_pos != self.left:
                    self.cds_right = stop_codon_left_pos - 1
                    self.has_exon_cds = True
                    self.last_exon = True
                else:
                    prev_exon_last = True
        return prev_exon_last


class Exons(deque):
    def __init__(self, chrom, strand, seleno):
        super(Exons, self).__init__()
        self.strand = strand
        self.chrom = chrom
        self.mock_node = ExonNode(
            _MOCK_LEFT, _MOCK_RIGHT, self.strand, self.chrom, []
        )
        self.reverse = strand == '-'
        if seleno:
            self.codon_table = _SELENO_TRANSLATION_DICT
        elif self.chrom in ["chrM", "M", "chrMT", "MT"]:
            self.codon_table = _MITOCHONDRIAL_TRANSLATION_DICT
        else:
            self.codon_table = _TRANSLATION_DICT
        self.junction_fiveprime = None
        self.has_tx_cds = False
        self.hanging = True
        self.annotated_end = True

    @property
    def has_stop_codon(self):
        for exon in self:
            if exon.has_stop_codon:
                return True
        return False

    @property
    def has_start_codon(self):
        for exon in self:
            if exon.has_start_codon:
                return True
        return False

    @property
    def translate(self):
        start = False
        stop = False
        for exon in self:
            if exon.has_start_codon:
                start = True
            if exon.has_stop_codon:
                stop = True
        return self.has_tx_cds or (start and stop)

    @property
    def cds_min(self):
        if self.reverse:
            for exon in reversed(self):
                if exon.has_stop_codon or exon.cds_left:
                    return max(exon.cds_left-3, exon.left)
            if self.has_tx_cds:
                return self[-1].left
            else:
                return 0
        else:
            for exon in self:
                if exon.has_start_codon or exon.cds_left:
                    # previous +3, which misses junction immediately after
                    # start codon
                    return exon.cds_left + 2
            return 0

    @property
    def cds_max(self):
        if self.reverse:
            for exon in self:
                if exon.has_start_codon or exon.cds_right:
                    # previously -3, which misses junction immediately after
                    # start codon)
                    return exon.cds_right - 2
            return 0
        else:
            for exon in reversed(self):
                if exon.has_stop_codon or exon.cds_right:
                    return exon.cds_right
            if self.has_tx_cds:
                return self[-1].right
            else:
                return 0

    def add_exon_node(self, left, right, gtf_line_list):
        feature = ExonNode(left, right, self.strand, self.chrom, gtf_line_list)
        return self.append(feature)

    def add_cds_to_exon_node(self, left, right, reading_frame):
        if not self.has_tx_cds:
            self.has_tx_cds = True
            first_cds = True
            if self[-1].left != int(left) or self[-1].right != int(right):
                self.hanging = False
            if self.reverse:
                cds_left = int(left)
                cds_right = int(right) - reading_frame
            else:
                cds_left = int(left) + reading_frame
                cds_right = int(right)
        else:
            first_cds = False
            cds_left = int(left)
            cds_right = int(right)
        self[-1].add_cds(cds_left, cds_right, first_cds)
        return None

    def add_start_codon_to_exon_node(self, left, right):
        self.hanging = False
        self[-1].add_start_codon(int(left), int(right))
        return None

    def add_stop_codon_to_exon_node(self, left, right):
        for i, ex in enumerate(reversed(self), 1):
            if not (ex.left <= int(left) <= ex.right):
                continue
            prev_exon_has_final_cds = ex.add_stop_codon(int(left), int(right))
            if prev_exon_has_final_cds:
                prev_exon = len(self) - (i + 1)
                self[prev_exon].last_exon = True
                if self.reverse:
                    if not self[prev_exon].cds_left:
                        self[prev_exon].cds_left = self[prev_exon].left
                else:
                    if not self[prev_exon].cds_right:
                        self[prev_exon].cds_right = self[prev_exon].right
            break
        return None

    def insert_exon_from_existing(self, left, right, index, model_line_list,
                                  new_cds_coords, new_exon_stop_codon,
                                  new_exon_last):
        new_line = deepcopy(model_line_list)
        new_line[3] = str(left)
        new_line[4] = str(right)
        try:
            new_line[8]['exon_number'] += '+'
        except KeyError:
            pass
        new_line[8]['exon_id'] = new_line[8]['exon_id'] + '_added'
        new_exon = ExonNode(left, right, self.strand, self.chrom, new_line)
        if new_cds_coords:
            new_exon.add_cds(new_cds_coords[0], new_cds_coords[1], False)
        new_exon.has_stop_codon = new_exon_stop_codon
        new_exon.last_exon = new_exon_last
        return self.insert(index, new_exon)

    def insert_junction(self, junction, verbose):
        chrom, left, right, strand = junction.split(';')
        if chrom != self.chrom:
            if verbose:
                print(chrom, self.chrom)
                print("jx not added; doesn't match transcript chromosome")
            return None
        if strand != self.strand:
            if verbose:
                print("jx not added; strand doesn't match transcript strand")
            return None
        self.insert_intron_with_cds_checking(int(left), int(right), verbose)
        if self.reverse:
            self.junction_fiveprime = int(right)
        else:
            self.junction_fiveprime = int(left)
        return None

    def remove_all_cds(self):
        # Any transcript on which this is called will not be printed into the
        # gtf, since CDS or start/stop codons are required for translation.
        for exon in self:
            exon.has_start_codon = False
            exon.has_stop_codon = False
            exon.cds_left = 0
            exon.cds_right = 0
        self.has_tx_cds = False
        return None

    def insert_intron_with_cds_checking(self, intron_left, intron_right,
                                        verbose):
        if self.reverse:
            tx_handler = ij.setup_reverse(intron_left, intron_right)
        else:
            tx_handler = ij.setup_forward(intron_left, intron_right)
        init_cds_min = self.cds_min
        init_cds_max = self.cds_max
        if not (init_cds_min <= tx_handler.intron_5prime <= init_cds_max):
            if verbose:
                print("junction 5' end does not fall within CDS; skipping.")
                print("cds min, intron 5', cds_max")
                print(init_cds_min, tx_handler.intron_5prime, init_cds_max)
                print(
                    'current exon left and right:', self[0].right, self[0].left
                )
            self.remove_all_cds()
            return None
        search_jx_3prime_end = False
        orig_stop_codon = self.has_stop_codon
        for i, curr_exon in enumerate(self):
            # Set up logic and locations for + vs - strand transcripts
            tx_handler.current_exon = curr_exon
            if tx_handler.jx_starts_before_exon and not search_jx_3prime_end:
                # New junction does not start in an exon; don't include exon
                # in new gtf.
                if verbose:
                    print('junction starts before exon; skipping transcript')
                self.remove_all_cds()
                return None
            if tx_handler.jx_starts_in_current_exon:
                if verbose:
                    print('junction starts in current exon')
                # Upstream jx end is in the current exon; replace the 3' end
                # of the exon with the 5' end of the jx.
                orig_cds_left = curr_exon.cds_left
                orig_cds_right = curr_exon.cds_right
                tx_handler.reset_exon_3prime()
                # Stop codon is moved to the following exon.
                curr_exon.has_stop_codon = False
                curr_exon.last_exon = False
                if tx_handler.jx_ends_before_exon_end:
                    # The full jx is within this exon; create a second exon
                    # from the 3' end of the jx - 3' end of the exon.
                    new_stop_cdn = curr_exon.has_stop_codon
                    new_exon_last = curr_exon.last_exon
                    new_left, new_right = tx_handler.new_exon_coords()
                    if init_cds_min <= tx_handler.intron_3prime <= init_cds_max:
                        # Junction is entirely contained within CDS. Old stop
                        # codon and old CDS downstream end should be retained.
                        new_cds_coords = tx_handler.new_cds_coords(
                            orig_cds_left, orig_cds_right
                        )
                    else:
                        # Junction cuts out the end of the CDS and the stop
                        # codon if there is one.
                        results = self.add_novel_cds_end(
                            new_left, new_right, i, inplace=False
                        )
                        new_stop_cdn, new_cds_coords = results
                    self.insert_exon_from_existing(
                        new_left, new_right, i+1, curr_exon.gtf_line_list,
                        new_cds_coords, new_stop_cdn, new_exon_last
                    )
                    break
                search_jx_3prime_end = True
            elif tx_handler.jx_5prime_end_annotated:
                if verbose:
                    print("junction 5' end is annotated")
                # The 5' end of the junction matches the annotated 3' end of
                # the current exon; look for downstream jx end.
                search_jx_3prime_end = True
                curr_exon.last_exon = False
                curr_exon.has_stop_codon = False
                tx_handler.wipe_cds_5prime()
            elif search_jx_3prime_end:
                threeprime_state = self.get_threeprime_state(
                    tx_handler.jx_3prime_end_annotated,
                    tx_handler.jx_ends_before_exon_end
                )
                three_prime_end_found = self.process_threeprime_state(
                    threeprime_state, curr_exon, i, tx_handler.intron_3prime,
                    init_cds_min, init_cds_max
                )
                if three_prime_end_found:
                    search_jx_3prime_end = False
                    break
        if search_jx_3prime_end:
            self.add_mini_exon_to_end_of_transcript(
                intron_left, intron_right, orig_stop_codon
            )
        self.remove_mock_nodes()
        return None

    def process_threeprime_state(self, threeprime_state, curr_exon, i,
                                 intron_3prime, init_cds_min, init_cds_max):
        if threeprime_state == ts.ThreePrimeState.annotated:
            # the 3' end of the junction matches the annotated 5' end
            # of the current exon; both ends accounted for.
            return True
        elif threeprime_state == ts.ThreePrimeState.in_exon:
            # Current exon 5' end becomes junction 3' coordinate;
            # assumes exon is simply extended on 5' end.
            if self.reverse:
                curr_exon.right = intron_3prime
                curr_exon.gtf_line_list[4] = str(intron_3prime)
            else:
                curr_exon.left = intron_3prime
                curr_exon.gtf_line_list[3] = str(intron_3prime)

            if not init_cds_min < intron_3prime < init_cds_max:
                # Junction cuts out end of CDS and stop codon.
                self.add_novel_cds_end(
                    curr_exon.left, curr_exon.right, i, inplace=True
                )
            return True
        else:
            # If junction 5' end does not occur before right of exon:
            # Junction 3' end not yet reached; current exon skipped;
            # delete current exon node.
            self[i] = self.mock_node
            return False

    def get_threeprime_state(self, jx_3prime_end_annotated,
                             jx_ends_before_exon_end):
        if jx_3prime_end_annotated:
            return ts.ThreePrimeState.annotated
        elif jx_ends_before_exon_end:
            # Current exon 5' end becomes junction 3' coordinate;
            # assumes exon is simply extended on 5' end.
            return ts.ThreePrimeState.in_exon
        else:
            # If junction 5' end does not occur before right of exon:
            # Junction 3' end not yet reached; current exon skipped;
            # delete current exon node.
            return ts.ThreePrimeState.beyond_exon

    def remove_mock_nodes(self):
        while self.mock_node in self:
            self.remove(self.mock_node)
        return None

    def add_novel_cds_end(self, new_exon_left, new_exon_right, curr_idx,
                          inplace=False):
        """
        :param new_exon_left:
        :param new_exon_right:
        :param curr_idx:
        :param inplace:

        Junction cuts out the end of the CDS and the stop codon if there is
        one.
        - CDS should be continued if present; if not, set stop codon 33
          bases downstream.
        - Stop codon is removed, so cds left_pos and right values should match
          new exon left_pos and right.
        :return:
        """
        new_exon_has_stop_codon = self[curr_idx].has_stop_codon
        new_cds_coords = []
        exon_len = abs(new_exon_right - new_exon_left)
        if exon_len < 31 and len(self) > curr_idx + 1:
            # End of exon is too short for full 9mer
            # translation; continue translation into next exon
            if self.reverse:
                self[curr_idx + 1].cds_left = max(
                    self[curr_idx + 1].left + 4,
                    self[curr_idx + 1].right - (_TARGET_LENGTH - exon_len)
                )
            else:
                self[curr_idx + 1].cds_right = min(
                    self[curr_idx + 1].right - 4,
                    self[curr_idx + 1].left + (_TARGET_LENGTH - exon_len)
                )
            self[curr_idx + 1].has_stop_codon = self[curr_idx].has_stop_codon
            if inplace:
                self[curr_idx].cds_left = 0
                self[curr_idx].cds_right = 0
                self[curr_idx].has_stop_codon = False
                self[curr_idx].last_exon = False
            else:
                new_exon_has_stop_codon = False
        else:
            # Add arbitrary end to translation at end of exon
            if self.reverse:
                if inplace:
                    self[curr_idx].cds_left = new_exon_left + 4
                else:
                    new_cds_coords = [new_exon_left + 4, new_exon_right]
            else:
                if inplace:
                    self[curr_idx].cds_right = new_exon_right - 4
                else:
                    new_cds_coords = [new_exon_left, new_exon_right - 4]
        return new_exon_has_stop_codon, new_cds_coords

    def add_mini_exon_to_end_of_transcript(self, intron_left, intron_right,
                                           orig_stop_codon):
        # Add tiny mock exon at new junction's 3' end for translation
        # new_exon_length = 33
        new_exon_length = 153
        if self.reverse:
            new_left = intron_left - new_exon_length
            new_right = intron_left
            cds_left = new_left + 4
            cds_right = new_right
            stop_left = new_left + 1
        else:
            new_left = intron_right
            new_right = intron_right + new_exon_length
            cds_left = new_left
            cds_right = new_right - 4
            stop_left = new_right - 3
        new_line_list = []
        for exon in reversed(self):
            if exon.gtf_line_list:
                new_line_list = deepcopy(exon.gtf_line_list)
                break
        new_line_list[3] = str(new_left)
        new_line_list[4] = str(new_right)
        new_line_list[8]['exon_number'] += '+'
        new_line_list[8]['exon_id'] += '_novel'
        self.add_exon_node(new_left, new_right, new_line_list)
        self[-1].last_exon = True
        # So that the "include 3'" flag isn't wrongly raised in translation
        self.annotated_end = False
        if self.has_tx_cds:
            self.add_cds_to_exon_node(cds_left, cds_right, 0)
        if orig_stop_codon:
            self.add_stop_codon_to_exon_node(stop_left, stop_left + 2)
        return None

    def print(self):
        for item in self:
            print(item.to_string())
        return None

    def copy_exons(self):
        new_exons = Exons(self.chrom, self.strand)
        for exon in self:
            new_exons.add_exon_node(
                exon.left, exon.right, exon.gtf_line_list
            )
        return new_exons

    def sequence_ranges(self, verbose):
        start_pos = 0
        cds_count_to_jx = 0
        sequence_ranges = []
        continue_cds_count = True
        for i, exon in enumerate(self):
            if start_pos and not (exon.cds_left or exon.cds_right):
                sequence_ranges.append([exon.left, exon.right])
                if continue_cds_count:
                    cds_count_to_jx += exon.cds_length
            # Junctions have been inserted into transcripts; it's impossible
            # to have the entire CDS range in one exon.
            if exon.last_exon:
                # Ensure translation continues past annotated stop codon
                # coordinate in case of an alternate reading frame
                if verbose:
                    print('in exon.last_exon')
                if self.reverse:
                    left = min(exon.left, exon.right - _TARGET_LENGTH)
                    right = exon.right
                else:
                    left = exon.left
                    right = max(exon.right, exon.left + _TARGET_LENGTH)
                sequence_ranges.append([left, right])
                break
            if exon.cds_left:
                if self.reverse:
                    # Ensure translation continues past annotated stop codon
                    # coordinate in case of an alternate reading frame
                    sequence_ranges.append([
                        min(exon.cds_left, exon.right-_TARGET_LENGTH),
                        exon.right
                    ])
                    break
                else:
                    sequence_ranges.append([exon.cds_left, exon.right])
                    start_pos = exon.cds_left
                    cds_count_to_jx += exon.cds_length
            if exon.cds_right:
                if verbose:
                    print('in exon.cds_right')
                if self.reverse:
                    sequence_ranges.append([exon.left, exon.cds_right])
                    start_pos = exon.cds_right
                    cds_count_to_jx += exon.cds_length
                else:
                    # Ensure translation continues past annotated stop codon
                    # coordinate in case of an alternate reading frame
                    sequence_ranges.append([
                        exon.left,
                        max(exon.cds_right, exon.left+_TARGET_LENGTH)
                    ])
                    if verbose:
                        print(
                            'exon cds right/target lenghts are {}/{}'.format(
                                exon.cds_right, exon.left + _TARGET_LENGTH
                            )
                        )
                    break
            if self.junction_fiveprime in {exon.left, exon.right}:
                continue_cds_count = False
        return start_pos, cds_count_to_jx, sequence_ranges

    def trans_range_biexons(self, cds_count_to_jx, nucleotide_ct,
                            sequence_ranges, verbose):
        current_length = 0
        intron_rf = cds_count_to_jx % 3
        if verbose:
            print('initial rf is {}'.format(intron_rf))
        if intron_rf == 0:
            post_required = nucleotide_ct - 1
        elif intron_rf == 2:
            post_required = nucleotide_ct
        else:
            post_required = nucleotide_ct + 1
        post_found = 0
        nucs_to_transl = max(
            0, cds_count_to_jx - ((nucleotide_ct - 1) + intron_rf)
        )
        if verbose:
            print(
                'initial rf is {}; post required {}; nucleotides to '
                'translation {}'.format(
                    intron_rf, post_required, nucs_to_transl
                )
            )
        trans_start_count = nucs_to_transl
        translation_ranges = []
        in_trans_region = False
        more_exons_needed = True
        jx_fiveprime_reached = False
        jx_pos = 0
        for [l_pos, r_pos] in sequence_ranges:
            seq_len = (r_pos - l_pos) + 1
            if not in_trans_region:
                if seq_len <= nucs_to_transl:
                    nucs_to_transl -= seq_len
                    continue
                if self.junction_fiveprime in {l_pos, r_pos}:
                    translation_ranges.append([l_pos, r_pos])
                    current_length += seq_len
                    trans_start_count -= nucs_to_transl
                    jx_fiveprime_reached = True
                    for l, r in translation_ranges:
                        jx_pos += (r - l) + 1
                else:
                    if self.reverse:
                        trans_start = r_pos - nucs_to_transl
                        translation_ranges.append([l_pos, trans_start])
                    else:
                        trans_start = l_pos + nucs_to_transl
                        translation_ranges.append([trans_start, r_pos])
                in_trans_region = True
            elif more_exons_needed or not jx_fiveprime_reached:
                translation_ranges.append([l_pos, r_pos])
                current_length += seq_len
                if self.junction_fiveprime in {l_pos, r_pos}:
                    jx_fiveprime_reached = True
                    for l, r in translation_ranges:
                        jx_pos += (r - l) + 1
                else:
                    post_found += seq_len
                more_exons_needed = (
                    (current_length < _TARGET_LENGTH)
                    or
                    (post_found < post_required)
                )
                if not more_exons_needed:
                    break
        correct_rf = (3 - trans_start_count) % 3
        if verbose:
            print('updated rf is {} (translation start count {})'.format(
                correct_rf, trans_start_count
            ))
        fivepr_end = int(sequence_ranges[0][0] == translation_ranges[0][0])
        threepr_end = int(sequence_ranges[-1][1] == translation_ranges[-1][1])
        return correct_rf, translation_ranges, jx_pos, fivepr_end, threepr_end

    def translate_jx(self, kmer_len, ref_fasta, samtools_path, verbose):
        results = self.get_transl_ranges(kmer_len, verbose)
        reading_frame, transl_ranges, jx_pos, fivepr_end, threepr_end = results
        dna_seq = get_dna_seq(
            self.chrom, transl_ranges, ref_fasta, samtools_path, self.strand
        )
        min_peptide_length = floor(jx_pos / 3) + 1
        min_peptide_length -= reading_frame
        pep_seqs, stopcodons = seq_to_protein(
            dna_seq, self.strand, min_peptide_length, self.codon_table
        )
        inframe_peptide = pep_seqs[reading_frame]
        if stopcodons[reading_frame]:
            # In case of novel stop codon/frame shift junction
            threepr_end = 1
        if verbose:
            print(pep_seqs)
            print('reading frame is:', reading_frame)
            print(
                "overlaps 5' end: {} and 3' end: {}".format(
                    bool(fivepr_end), bool(threepr_end)
                )
            )
            print('predicted peptide is:', inframe_peptide)
        flags = [jx_pos, fivepr_end, threepr_end]
        return [inframe_peptide, pep_seqs], flags, reading_frame

    def get_transl_ranges(self, kmer_len, verbose):
        # # For bi-exon or 9-mer translation around the junction only:
        # nucleotide_ct = (kmer_len * 3) - 2
        # For 50 amino acids on either side of the junction:
        nucleotide_ct = _TARGET_LENGTH
        start_pos, cds_count_to_jx, sequence_ranges = self.sequence_ranges(verbose)
        if verbose:
            print('sequence ranges are:')
            print(sequence_ranges)
        results = self.trans_range_biexons(
            cds_count_to_jx, nucleotide_ct, sequence_ranges, verbose
        )
        reading_frame, transl_ranges, jx_pos, fivepr_end, threepr_end = results
        if not self.annotated_end:
            threepr_end = 0
        if verbose:
            print('translation ranges are:')
            print(transl_ranges)
        return reading_frame, transl_ranges, jx_pos, fivepr_end, threepr_end


class TranscriptInfo(deque):
    def __init__(self):
        super(TranscriptInfo, self).__init__()

    def add_info(self, line_info_list):
        self.append(line_info_list)
        return None


class Transcript(object):
    FASTA_MAP = {}
    TRYPSIN_SITES = {'exon_size': [], 'cleavage_sites': []}

    def __init__(self, name, chrom, strand, seleno):
        self.name = name
        self.chrom = chrom
        self.strand = strand
        self.exons = Exons(chrom, strand, seleno)
        self.tx_info = TranscriptInfo()
        self.print_cds = False
        self.seleno = seleno

    @property
    def translate(self):
        return self.exons.translate

    @property
    def reverse(self):
        return self.strand == '-'

    @property
    def has_cds(self):
        return self.exons.has_tx_cds

    def add_exon(self, left, right, line):
        self.exons.add_exon_node(left, right, self.parse_gtf_line(line))
        return None

    def add_feature_to_exon(self, feature_type, left, right, reading_frame):
        if feature_type == 'CDS' and reading_frame != '.':
            self.print_cds = True
            self.exons.add_cds_to_exon_node(left, right, int(reading_frame))
        elif feature_type == 'start_codon':
            self.exons.add_start_codon_to_exon_node(left, right)
        elif feature_type == 'stop_codon':
            self.exons.add_stop_codon_to_exon_node(left, right)
        return None

    def add_non_exon_info(self, line):
        self.tx_info.add_info(self.parse_gtf_line(line))
        return None

    def parse_gtf_line(self, line):
        if type(line) == list:
            return line
        else:
            items = line.split('\t')
            vals_dict = {}
            for attr in items[8].split(';')[:-1]:
                attr, _, val = attr.strip().partition(' ')
                vals_dict[attr] = val.strip('"')
            items[8] = vals_dict
            return items

    def print(self):
        print('transcript has cds: {}'.format(self.exons.has_tx_cds))
        print('transcript is hanging: {}'.format(self.exons.hanging))
        for item in self.exons:
            print('exon has cds: {}'.format(item.has_exon_cds))
            print(item.cds_left, item.left, item.cds_right, item.right)
            print(item.to_string())
        return None

    def copy_transcript(self, new_name):
        new_tx = Transcript(new_name, self.chrom, self.strand, self.seleno)
        for exon in self.exons:
            new_line_list = deepcopy(exon.gtf_line_list)
            new_line_list[8]['transcript_id'] = new_name
            new_tx.add_exon(exon.left, exon.right, new_line_list)
            new_tx.exons.add_cds_to_exon_node(
                str(exon.cds_left), str(exon.cds_right), exon.cds_rf
            )
            new_tx.print_cds = self.print_cds
            if exon.has_start_codon:
                new_tx.add_feature_to_exon(
                    'start_codon', str(exon.start_codon_left),
                    str(exon.start_codon_left + exon.start_codon_add),
                    exon.cds_rf
                )
            if exon.has_stop_codon:
                new_tx.add_feature_to_exon(
                    'stop_codon', str(exon.stop_codon_left),
                    str(exon.stop_codon_left + exon.stop_codon_add),
                    exon.cds_rf
                )
            if exon.last_exon:
                new_tx.exons[-1].last_exon = True
                if self.reverse:
                    if new_tx.exons[-1].cds_left == 0:
                        new_tx.exons[-1].cds_left = new_tx.exons[-1].left
                else:
                    if new_tx.exons[-1].cds_right == 0:
                        new_tx.exons[-1].cds_right = new_tx.exons[-1].right
            new_tx.exons[-1].has_exon_cds = exon.has_exon_cds
        for line_list in self.tx_info:
            new_line_list = deepcopy(line_list)
            if 'transcript_id' in new_line_list[8]:
                new_line_list[8]['transcript_id'] = new_name
            new_tx.add_non_exon_info(new_line_list)
        new_tx.exons.hanging = self.exons.hanging
        return new_tx

    def insert_junction(self, junction, verbose=False):
        if self.translate:
            if verbose:
                print('inserting junction')
                print(junction, self.chrom)
            self.exons.insert_junction(junction, verbose)
        return None

    def translate_jx_peptides(self, ref_fasta, samtools_path, kmer_len=9,
                              verbose=False, incl_hang=False,
                              hanging_only=False, kmers=True):
        meets_translation_criteria = (
            self.exons.has_tx_cds and
            ((hanging_only and self.exons.hanging)
             or incl_hang or not self.exons.hanging)
        )
        if meets_translation_criteria:
            peptide_list, flags, reading_frame = self.exons.translate_jx(
                kmer_len, ref_fasta, samtools_path, verbose
            )
            jx_pos = flags[0] - reading_frame
            if jx_pos % 3 == 0:
                jx_between_codons = 1
            else:
                jx_between_codons = 0
            jx_pos = int(floor(jx_pos / 3))
            if jx_between_codons:
                # Fixes issue found Apr 2023
                jx_pos = jx_pos - 1
            flags[0] = jx_pos
            flags.append(jx_between_codons)

            if len(peptide_list[0]) <= jx_pos:
                if verbose:
                    print(
                        'junction position {} past end of peptide ({}), '
                        'removing:'.format(jx_pos, len(peptide_list[0]))
                    )
                peptide_list[0] = ''

            left_side = peptide_list[0][:jx_pos+1]
            right_side = peptide_list[0][jx_pos+1:]
            for exon_peptide in [left_side, right_side]:
                Transcript.TRYPSIN_SITES['exon_size'].append(len(exon_peptide))
                clv_site_num = exon_peptide.count('K')
                clv_site_num += exon_peptide.count('R')
                clv_site_num -= exon_peptide.count('KP')
                clv_site_num -= exon_peptide.count('RP')
                Transcript.TRYPSIN_SITES['cleavage_sites'].append(clv_site_num)

            if kmers:
                inframe_peptide = peptide_list[0]
                start_loc = max(0, jx_pos - (kmer_len - 1))
                end_loc = jx_pos + (kmer_len - 1)
                if not jx_between_codons:
                    end_loc += 1
                inframe_peptide = inframe_peptide[start_loc:end_loc]
                if len(inframe_peptide) < kmer_len:
                    inframe_peptide = ''
                peptide_list.append(inframe_peptide)
            else:
                peptide_list.append('')
        else:
            if verbose:
                print(
                    'not translating; has CDS {}, including hanging {}, exons '
                    'are hanging {}'.format(
                        self.exons.has_tx_cds, incl_hang, self.exons.hanging
                    )
                )
            peptide_list = ['', '', '']
            flags = [0, 0, 0, 0]
        return peptide_list, flags

    def get_translation_ranges(self, kmer_len=9):
        if self.exons.has_tx_cds:
            return self.exons.translate_jx(kmer_len)
        else:
            return None, None
