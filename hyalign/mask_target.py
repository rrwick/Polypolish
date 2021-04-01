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

import collections
import re

from .alignment import get_multi_alignment_read_names
from .log import log, section_header, explanation
from .misc import load_fasta


def mask_target_sequences(read_pair_names, alignments, target):
    section_header('Masking target sequences')
    explanation('Hyalign now masks out any bases of the target sequence which lack significant '
                'ambiguity. E.g. if a target base is "A" and all read alignments agree with that '
                'base, then that position is not informative and will be masked. But if the read '
                'alignments show ambiguity (some say "A" and some say "T"), then that position '
                'is informative and will not be masked.')

    target_seqs = load_fasta(target)
    repetitive_regions = find_repetitive_regions(target_seqs, read_pair_names, alignments)
    for target_name, target_seq in target_seqs:
        mask_one_target_sequence(target_name, target_seq, alignments, read_pair_names,
                                 repetitive_regions[target_name])


def find_repetitive_regions(target_seqs, read_pair_names, alignments):
    repetitive_regions = {name: set() for name, _ in target_seqs}
    multi_alignment_read_names = get_multi_alignment_read_names(read_pair_names, alignments)

    # TODO: find repetitive regions of the target by looking for all positions covered by
    #       multi-alignment reads.

    return repetitive_regions


def mask_one_target_sequence(target_name, target_seq, alignments, read_pair_names,
                             repetitive_regions):
    pileup = collections.defaultdict(list)
    read_names = [n + '/1' for n in read_pair_names] + [n + '/2' for n in read_pair_names]
    for name in read_names:
        read_alignments = alignments[name]
        for a in read_alignments:
            if a.ref_name != target_name:
                continue
            if not alignment_overlaps_repeats(a, repetitive_regions):
                continue
            ref_seq = target_seq[a.ref_start:a.ref_end]
            if a.masked_read_seq is None:
                read_seq = a.read_seq
            else:
                read_seq = a.masked_read_seq
            aligned_bases = get_read_bases_for_each_target_base(read_seq, ref_seq, a.cigar)
            for i, bases in enumerate(aligned_bases):
                ref_pos = a.ref_start + i
                if 'N' not in bases:
                    pileup[ref_pos].append(bases)

    for ref_pos in sorted(pileup.keys()):
        print(ref_pos, ' '.join(pileup[ref_pos]))


def alignment_overlaps_repeats(a, repetitive_regions):
    # TODO: return True if this alignment's ref range overlaps any of the positions in the
    #       repetitive_regions set.
    return True  # TEMP


def get_read_bases_for_each_target_base(read_seq, ref_seq, cigar):
    expanded_cigar = get_expanded_cigar(cigar)
    i, j = 0, 0
    read_bases = ['' for _ in range(len(ref_seq))]
    for c in expanded_cigar:
        if c == 'M':
            read_bases[j] += read_seq[i]
            i += 1
            j += 1
        elif c == 'I':
            read_bases[j] += read_seq[i]
            i += 1
        elif c == 'D':
            read_bases[j] += '-'
            j += 1
        else:
            assert False

    assert i == len(read_seq)
    assert j == len(ref_seq)
    return read_bases


def get_expanded_cigar(cigar):
    expanded_cigar = []
    cigar_parts = re.findall(r'\d+[MIDNSHP=X]', cigar)
    for p in cigar_parts:
        size = int(p[:-1])
        letter = p[-1]
        assert letter == 'M' or letter == 'D' or letter == 'I'  # no clips in end-to-end alignment
        expanded_cigar.append(letter * size)
    return ''.join(expanded_cigar)
