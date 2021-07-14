"""
Copyright 2021 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Polypolish

This file is part of Polypolish. Polypolish is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Polypolish is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Polypolish.
If not, see <http://www.gnu.org/licenses/>.
"""

import collections
import re

from .log import log, section_header, quit_with_error
from .misc import reverse_complement


def load_alignments(sam_filenames, max_errors):
    section_header('Loading alignments')
    alignments = collections.defaultdict(list)
    for i, sam_filename in enumerate(sam_filenames):
        read_name_suffix = f'_{i+1}'
        load_alignments_one_file(sam_filename, alignments, read_name_suffix)
    log()
    add_secondary_read_seqs(alignments)
    filter_alignments(alignments, max_errors)
    return alignments


def load_alignments_one_file(sam_filename, alignments, read_name_suffix):
    log(f'{sam_filename}: ', end='')
    alignment_count = 0
    read_names = set()
    with open(sam_filename, 'rt') as sam_file:
        for line in sam_file:
            if line.startswith('@'):  # header line
                continue
            alignment = Alignment(line)
            if alignment.is_aligned():
                alignment.read_name += read_name_suffix
                alignments[alignment.read_name].append(alignment)
                alignment_count += 1
                read_names.add(alignment.read_name)
    log(f'{alignment_count:,} alignments from {len(read_names):,} reads')
    return alignments


def filter_alignments(alignments, max_errors):
    log('Filtering for high-quality end-to-end alignments:')
    keep_count, discard_count = 0, 0
    for name in list(alignments.keys()):
        good_alignments = []
        for a in alignments[name]:
            if a.starts_and_ends_with_match() and a.mismatches <= max_errors:
                good_alignments.append(a)
                keep_count += 1
            else:
                discard_count += 1
            alignments[name] = good_alignments
    log(f'  {keep_count:,} alignment{"" if keep_count == 1 else "s"} kept')
    log(f'  {discard_count:,} alignment{"" if keep_count == 1 else "s"} discarded')
    log()


def add_secondary_read_seqs(alignments):
    """
    Secondary alignments don't have read sequences or qualities, but we will need those later, so
    we add them in now.
    """
    for name in list(alignments.keys()):
        if len(alignments[name]) == 0:
            continue
        seq, qual = None, None
        for a in alignments[name]:
            read_seq = a.read_seq
            read_qual = a.read_qual
            if a.is_secondary():
                assert read_seq == '*' and read_qual == '*'
            else:
                assert read_seq != '*' and read_qual != '*'
                if a.is_on_forward_strand():
                    seq = read_seq
                    qual = read_qual
                else:
                    seq = reverse_complement(read_seq)
                    qual = read_qual[::-1]
                break
        assert seq is not None and qual is not None
        for a in alignments[name]:
            if a.is_secondary():
                if a.is_on_forward_strand():
                    a.read_seq = seq
                    a.read_qual = qual
                else:
                    a.read_seq = reverse_complement(seq)
                    a.read_qual = qual[::-1]
        for a in alignments[name]:
            assert a.read_seq != '*' and a.read_qual != '*'


class Alignment(object):

    def __init__(self, sam_line):
        parts = sam_line.strip().split('\t')
        if len(parts) < 11:
            quit_with_error('\nError: alignment file does not seem to be in SAM format')

        self.read_name = parts[0]
        self.sam_flags = int(parts[1])
        self.ref_name = parts[2]
        self.ref_start = int(parts[3]) - 1  # switch from SAM's 1-based to Python's 0-based
        self.cigar = parts[5]
        self.read_seq = parts[9]
        self.read_qual = parts[10]

        self.ref_end = get_ref_end(self.ref_start, self.cigar)

        self.mismatches = -1
        for p in parts[11:]:
            if p.startswith('NM:i:'):
                self.mismatches = int(p[5:])

    def __repr__(self):
        return f'{self.read_name}:{self.ref_name}:{self.ref_start}-{self.ref_end}:{self.mismatches}'

    def is_aligned(self):
        return not self.has_flag(4)

    def is_secondary(self):
        return self.has_flag(256)

    def is_on_forward_strand(self):
        return not self.has_flag(16)

    def has_flag(self, flag):
        return bool(self.sam_flags & flag)

    def starts_and_ends_with_match(self):
        cigar_parts = re.findall(r'\d+[MIDNSHP=X]', self.cigar)
        first_part, last_part = cigar_parts[0], cigar_parts[-1]
        return first_part[-1] == 'M' and last_part[-1] == 'M'

    def get_read_bases_for_each_target_base(self, ref_seq):
        expanded_cigar = get_expanded_cigar(self.cigar)
        i = 0
        read_bases = []
        for c in expanded_cigar:
            if c == 'M':
                read_bases.append(self.read_seq[i])
                i += 1
            elif c == 'I':
                read_bases[-1] += self.read_seq[i]
                i += 1
            elif c == 'D':
                read_bases.append('-')
            else:
                assert False
        assert i == len(self.read_seq)
        assert len(read_bases) == len(ref_seq)
        return read_bases


def get_ref_end(ref_start, cigar):
    ref_end = ref_start
    cigar_parts = re.findall(r'\d+[MIDNSHP=X]', cigar)
    for p in cigar_parts:
        size = int(p[:-1])
        letter = p[-1]
        if letter == 'M' or letter == 'D' or letter == 'N' or letter == '=' or letter == 'X':
            ref_end += size
    return ref_end


def get_expanded_cigar(cigar):
    expanded_cigar = []
    cigar_parts = re.findall(r'\d+[MIDNSHP=X]', cigar)
    for p in cigar_parts:
        size = int(p[:-1])
        letter = p[-1]
        assert letter == 'M' or letter == 'D' or letter == 'I'  # no clips in end-to-end alignment
        expanded_cigar.append(letter * size)
    return ''.join(expanded_cigar)
