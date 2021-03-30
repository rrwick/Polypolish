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

from .log import log, section_header, explanation
from .misc import load_fasta


def mark_read_sequences(read_names, alignments, target):
    section_header('Masking read sequences')
    explanation('Before continuing, Hyalign masks out read bases which do not match the reference '
                'at any location. I.e. if a read base appears to be a mistake at each of the '
                'read\'s alignment locations, it is assumed to be a read error.')

    target_seqs = dict(load_fasta(target))

    multi_alignment_read_names = []
    for name in sorted(read_names):
        name_1, name_2 = name + '/1', name + '/2'
        alignments_1, alignments_2 = alignments[name_1], alignments[name_2]
        if len(alignments_1) > 1:
            multi_alignment_read_names.append(name_1)
        if len(alignments_2) > 1:
            multi_alignment_read_names.append(name_2)

    for name in multi_alignment_read_names:
        read_alignments = alignments[name]
        read_mismatches_all_alignments = None
        for a in read_alignments:
            ref_seq = target_seqs[a.ref_name][a.ref_start:a.ref_end]

            aligned_read_seq, aligned_ref_seq, diffs, read_mismatches, _ = \
                align_seqs(a.read_seq, ref_seq, a.cigar)
            if a.is_on_reverse_strand():
                read_mismatches = flip_positions(read_mismatches, len(a.read_seq))
            if read_mismatches_all_alignments is None:
                read_mismatches_all_alignments = read_mismatches
            else:
                read_mismatches_all_alignments &= read_mismatches

            print()  # TEMP
            print(a)  # TEMP
            print(a.cigar, a.is_on_forward_strand())  # TEMP
            print(aligned_read_seq)  # TEMP
            print(aligned_ref_seq)  # TEMP
            print(diffs)  # TEMP
        print(read_mismatches_all_alignments)  # TEMP

        for a in read_alignments:
            create_masked_read_seq(a, read_mismatches_all_alignments)
            print(a.masked_read_seq)  # TEMP


        print()  # TEMP
        print()  # TEMP
        print()  # TEMP


def align_seqs(seq_1, seq_2, cigar):
    expanded_cigar = get_expanded_cigar(cigar)
    i, j = 0, 0
    aligned_seq_1, aligned_seq_2, differences = [], [], []
    seq_1_mismatchs, seq_1_mismatchs = set(), set()
    for c in expanded_cigar:
        if c == 'M':
            b_1 = seq_1[i]
            b_2 = seq_2[j]
            if b_1 == b_2:
                diff = ' '
            else:
                diff = '*'
                seq_1_mismatchs.add(i)
                seq_1_mismatchs.add(j)
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
    return aligned_seq_1, aligned_seq_2, differences, seq_1_mismatchs, seq_1_mismatchs


def get_expanded_cigar(cigar):
    expanded_cigar = []
    cigar_parts = re.findall(r'\d+[MIDNSHP=X]', cigar)
    for p in cigar_parts:
        size = int(p[:-1])
        letter = p[-1]
        assert letter == 'M' or letter == 'D' or letter == 'I'  # no clips in end-to-end alignment
        expanded_cigar.append(letter * size)
    return ''.join(expanded_cigar)


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
