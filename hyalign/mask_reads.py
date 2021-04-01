"""
Copyright 2021 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Hyalign

This file is part of Hyalign. Hyalign is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Hyalign is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Hyalign.
If not, see <http://www.gnu.org/licenses/>.
"""

import re

from .alignment import get_multi_alignment_read_names
from .log import log, section_header, explanation
from .misc import load_fasta


def mask_read_sequences(read_pair_names, alignments, target):
    section_header('Masking read sequences')
    explanation('Before continuing, Hyalign masks out read bases which do not match the reference '
                'at any location. I.e. if a read base appears to be a mistake at each of the '
                'read\'s alignment locations, it is assumed to be a read error and ignored from '
                'this point onward.')

    target_seqs = dict(load_fasta(target))
    multi_alignment_read_names = get_multi_alignment_read_names(read_pair_names, alignments)

    total_masked_bases = 0
    for name in multi_alignment_read_names:
        read_alignments = alignments[name]
        mask_positions = get_mask_positions(read_alignments, target_seqs)
        if mask_positions:
            for a in read_alignments:
                create_masked_read_seq(a, mask_positions)
        total_masked_bases += len(mask_positions)

    log(f'Masked {total_masked_bases:,} bases in {len(multi_alignment_read_names):,} reads')
    log()


def get_mask_positions(read_alignments, target_seqs):
    """
    Returns a list of positions in the read to mask. Positions are given in terms of the forward
    strand.
    """
    mask_positions = None
    for a in read_alignments:
        ref_seq = target_seqs[a.ref_name][a.ref_start:a.ref_end]
        aligned_read_seq, aligned_ref_seq, diffs, read_mismatches = \
            align_seqs_with_cigar(a.read_seq, ref_seq, a.cigar)

        # Sanity check: if there are too many mismatches, something has gone terribly wrong!
        assert len(read_mismatches) < len(a.read_seq) // 2

        if a.is_on_reverse_strand():
            read_mismatches = flip_positions(read_mismatches, len(a.read_seq))
        if mask_positions is None:
            mask_positions = read_mismatches
        else:
            mask_positions &= read_mismatches
    return mask_positions


def create_masked_read_seq(alignment, mask_positions):
    """
    This function creates a masked version of the read sequence in the alignment, where all mask
    positions are replaced with 'N' bases.
    """
    alignment.masked_read_seq = alignment.read_seq[:]
    if alignment.is_on_reverse_strand():
        mask_positions = flip_positions(mask_positions, len(alignment.read_seq))
    for i in mask_positions:
        alignment.masked_read_seq = \
            alignment.masked_read_seq[:i] + 'N' + alignment.masked_read_seq[i + 1:]


def flip_positions(positions, seq_length):
    """
    This function takes a set of positions (0-based) and flips them to the reverse strand positions.
    """
    flipped_positions = set()
    for p in positions:
        flipped_positions.add(seq_length - p - 1)
    return flipped_positions


def align_seqs_with_cigar(seq_1, seq_2, cigar):
    expanded_cigar = get_expanded_cigar(cigar)
    i, j = 0, 0
    aligned_seq_1, aligned_seq_2, differences = [], [], []
    seq_1_mismatches = set()
    for c in expanded_cigar:
        if c == 'M':
            b_1 = seq_1[i]
            b_2 = seq_2[j]
            if b_1 == b_2:
                diff = ' '
            else:
                diff = '*'
                seq_1_mismatches.add(i)
            i += 1
            j += 1
        elif c == 'I':
            b_1 = seq_1[i]
            b_2 = '-'
            diff = '*'
            i += 1
        elif c == 'D':
            b_1 = '-'
            b_2 = seq_2[j]
            diff = '*'
            j += 1
        else:
            assert False
        aligned_seq_1.append(b_1)
        aligned_seq_2.append(b_2)
        differences.append(diff)
    assert i == len(seq_1)
    assert j == len(seq_2)
    aligned_seq_1 = ''.join(aligned_seq_1)
    aligned_seq_2 = ''.join(aligned_seq_2)
    differences = ''.join(differences)
    assert aligned_seq_1.replace('-', '') == seq_1
    assert aligned_seq_2.replace('-', '') == seq_2
    return aligned_seq_1, aligned_seq_2, differences, seq_1_mismatches


def get_expanded_cigar(cigar):
    expanded_cigar = []
    cigar_parts = re.findall(r'\d+[MIDNSHP=X]', cigar)
    for p in cigar_parts:
        size = int(p[:-1])
        letter = p[-1]
        assert letter == 'M' or letter == 'D' or letter == 'I'  # no clips in end-to-end alignment
        expanded_cigar.append(letter * size)
    return ''.join(expanded_cigar)
