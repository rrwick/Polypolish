"""
Copyright 2021 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Poligner

This file is part of Poligner. Poligner is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Poligner is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Poligner.
If not, see <http://www.gnu.org/licenses/>.
"""

import collections
import pathlib
import re
import tempfile

from .log import log, section_header, explanation, quit_with_error
from .misc import run_command, iterate_fastq, reverse_complement


def align_reads(target, short1, short2, threads, max_errors):
    section_header('Aligning short reads to target sequence')
    explanation('Poligner uses minimap2 to align the short reads to the target sequence. The '
                'alignment is done in an unpaired manner, and all end-to-end alignments are kept '
                'for each read.')

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)
        log('Creating temp directory for working files:')
        log(f'  {temp_dir}')
        log()
        read_filename, read_count, read_pair_names = \
            combine_reads_into_one_file(short1, short2,temp_dir)
        alignments = align_with_minimap2(read_filename, read_pair_names, target, threads)
        add_secondary_read_seqs(alignments, read_pair_names)
        filter_alignments(alignments, read_pair_names, max_errors)
        print_alignment_info(alignments, read_count, read_pair_names)
        return alignments, read_pair_names, read_count


def combine_reads_into_one_file(short1, short2, temp_dir):
    log(f'Combining paired reads into one file:')
    read_filename = str(temp_dir / 'reads.fastq')
    first_count, second_count, total_count = 0, 0, 0
    first_names, second_names = set(), set()
    with open(read_filename, 'wt') as r:
        for name, header, sequence, qualities in iterate_fastq(short1):
            if not name.endswith('/1'):
                quit_with_error('Error: read names in first file must end with "/1"')
            r.write(f'{header}\n{sequence}\n+\n{qualities}\n')
            short_name = name[:-2]
            if short_name in first_names:
                quit_with_error(f'Error: duplicate read name {name}')
            first_names.add(short_name)
            first_count += 1
            total_count += 1
        for name, header, sequence, qualities in iterate_fastq(short2):
            if not name.endswith('/2'):
                quit_with_error('Error: read names in second file must end with "/2"')
            r.write(f'{header}\n{sequence}\n+\n{qualities}\n')
            short_name = name[:-2]
            if short_name in second_names:
                quit_with_error(f'Error: duplicate read name {name}')
            second_names.add(short_name)
            second_count += 1
            total_count += 1
    if first_count != second_count:
        quit_with_error('Error: unequal number of reads in each file')
    if first_names != second_names:
        quit_with_error('Error: read names in first and second files do not match')
    log(f'  {total_count:,} reads')
    log(f'  {read_filename}')
    log()
    return read_filename, total_count, first_names


def align_with_minimap2(reads, read_pair_names, target, threads):
    log(f'Aligning reads with minimap2:')
    command = ['minimap2', '-ax', 'sr', '--secondary=yes', '-p', '0.1', '-N', '1000000',
               '-t', str(threads), target, reads]
    log('  ' + ' '.join(command))
    stdout, stderr, return_code = run_command(command)
    alignments = {}
    alignment_count = 0
    for name in read_pair_names:
        alignments[name + '/1'] = []
        alignments[name + '/2'] = []
    for line in stdout.splitlines():
        if line.startswith('@'):
            continue
        alignment = Alignment(line)
        if alignment.is_aligned():
            alignments[alignment.read_name].append(alignment)
            alignment_count += 1
    log(f'  {alignment_count:,} total alignments')
    log()
    return alignments


def add_secondary_read_seqs(alignments, read_pair_names):
    """
    Secondary alignments don't have read sequences or qualities, but we will need those later, so
    we add them in now.
    """
    read_names = [n + '/1' for n in read_pair_names] + [n + '/2' for n in read_pair_names]
    for name in read_names:
        if len(alignments[name]) == 0:
            continue
        seq, qual = None, None
        for a in alignments[name]:
            if a.is_secondary():
                assert a.read_seq == '*' and a.read_qual == '*'
            else:
                assert a.read_seq != '*' and a.read_qual != '*'
                if a.is_on_forward_strand():
                    seq = a.read_seq
                    qual = a.read_qual
                else:
                    seq = reverse_complement(a.read_seq)
                    qual = a.read_qual[::-1]
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


def filter_alignments(alignments, read_pair_names, max_errors):
    log('Filtering for high quality end-to-end alignments... ', end='')
    read_names = [n + '/1' for n in read_pair_names] + [n + '/2' for n in read_pair_names]
    for name in read_names:
        alignments[name] = [a for a in alignments[name]
                            if a.starts_and_ends_with_match() and a.nm_tag <= max_errors]
    log('done')
    log()


def print_alignment_info(alignments, read_count, read_pair_names):
    log('Summary:')
    alignment_count, zero_count, single_count, multi_count = 0, 0, 0, 0
    incomplete_pair_count, unique_pair_count, multi_pair_count = 0, 0, 0

    for name in read_pair_names:
        name_1 = name + '/1'
        name_2 = name + '/2'
        count_1 = len(alignments[name_1])
        count_2 = len(alignments[name_2])
        alignment_count += count_1
        alignment_count += count_2
        if count_1 == 0:
            zero_count += 1
        elif count_1 == 1:
            single_count += 1
        else:
            multi_count += 1
        if count_2 == 0:
            zero_count += 1
        elif count_2 == 1:
            single_count += 1
        else:
            multi_count += 1
        if count_1 == 0 or count_2 == 0:
            incomplete_pair_count += 1
        elif count_1 == 1 and count_2 == 1:
            unique_pair_count += 1
        else:
            assert count_1 > 1 or count_2 > 1
            multi_pair_count += 1

    log(f'  {alignment_count:,} alignments')
    number_size = len(f'{alignment_count:,}')
    log()
    log(f'  {f"{zero_count:,}".rjust(number_size)} reads have no alignments')
    log(f'  {f"{single_count:,}".rjust(number_size)} reads have one alignment')
    log(f'  {f"{multi_count:,}".rjust(number_size)} reads have multiple alignments')
    log()
    log(f'  {f"{incomplete_pair_count:,}".rjust(number_size)} read pairs are incomplete')
    log(f'  {f"{unique_pair_count:,}".rjust(number_size)} read pairs are uniquely aligned')
    log(f'  {f"{multi_pair_count:,}".rjust(number_size)} read pairs have multiple combinations')
    log()
    assert read_count == zero_count + single_count + multi_count
    assert read_count // 2 == incomplete_pair_count + unique_pair_count + multi_pair_count


def get_multi_alignment_read_names(read_pair_names, alignments):
    """
    Returns a list of the names of individual reads (i.e. not read pairs) with more than one
    alignment.
    """
    multi_alignment_read_names = []
    for name in sorted(read_pair_names):
        name_1, name_2 = name + '/1', name + '/2'
        alignments_1, alignments_2 = alignments[name_1], alignments[name_2]
        if len(alignments_1) > 1:
            multi_alignment_read_names.append(name_1)
        if len(alignments_2) > 1:
            multi_alignment_read_names.append(name_2)
    return multi_alignment_read_names


class Alignment(object):

    def __init__(self, sam_line):
        self.sam_line = sam_line.strip()
        parts = self.sam_line.split('\t')
        if len(parts) < 11:
            quit_with_error('\nError: alignment file does not seem to be in SAM format')

        self.sam_line = sam_line
        self.read_name = parts[0]
        self.flags = int(parts[1])
        self.ref_name = parts[2]
        self.ref_start = int(parts[3]) - 1  # switch from SAM's 1-based to Python's 0-based
        self.map_q = int(parts[4])
        self.cigar = parts[5]
        self.ref_end = get_ref_end(self.ref_start, self.cigar)
        self.read_seq = parts[9]
        self.read_qual = parts[10]
        self.masked_read_seq = None

        self.nm_tag = None
        for p in parts:
            if p.startswith('NM:i:'):
                self.nm_tag = int(p[5:])
        assert self.nm_tag is not None

        # The following pieces of information are left empty for the time being. They will be
        # filled in later as needed for multi-alignment reads.
        self.aligned_read_seq = None
        self.aligned_ref_seq = None
        self.diffs = None
        self.read_error_positions = None
        self.ref_error_positions = None
        self.read_positions_to_ref_positions = None
        self.ref_positions_to_read_positions = None

    def __repr__(self):
        return f'{self.read_name}:{self.ref_name}:{self.ref_start}-{self.ref_end}'

    def has_flag(self, flag: int):
        return bool(self.flags & flag)

    def is_aligned(self):
        return not self.has_flag(4)

    def is_secondary(self):
        return self.has_flag(256)

    def is_on_reverse_strand(self):
        return self.has_flag(16)

    def is_on_forward_strand(self):
        return not self.has_flag(16)

    def has_no_indels(self):
        return self.cigar.endswith('M') and self.cigar[:-1].isdigit()

    def starts_and_ends_with_match(self):
        cigar_parts = re.findall(r'\d+[MIDNSHP=X]', self.cigar)
        first_part, last_part = cigar_parts[0], cigar_parts[-1]
        return first_part[-1] == 'M' and last_part[-1] == 'M'

    def add_detailed_alignment_info(self, ref_seq):
        """
        This function adds a lot of extra info to the alignment: aligned versions of the sequences,
        positions of errors in both sequences and the relationship between the positions in the
        sequences.

        All read positions are stored in terms of the reference-aligned strand. E.g. if a read
        aligned to the reverse strand of the reference, then a position of 0 according to this
        function actually corresponds to the end of the original read sequence.

        All reference positions are stored in terms of the entire reference sequence, not just the
        aligned part. E.g. if an alignment starts at position 431 of the reference, then the first
        position of the alignment according to this function is 431 (not 0).
        """
        self.aligned_read_seq, self.aligned_ref_seq, self.diffs = [], [], []
        self.read_error_positions, self.ref_error_positions = set(), set()
        self.read_positions_to_ref_positions = collections.defaultdict(set)
        self.ref_positions_to_read_positions = collections.defaultdict(set)

        expanded_cigar = get_expanded_cigar(self.cigar)
        i, j = 0, 0
        for c in expanded_cigar:
            if c == 'M':
                b_1 = self.read_seq[i]
                b_2 = ref_seq[j]
                if b_1 == b_2:
                    diff = ' '
                else:
                    diff = '*'
                    self.read_error_positions.add(i)
                    self.ref_error_positions.add(j + self.ref_start)
                self.read_positions_to_ref_positions[i].add(j + self.ref_start)
                self.ref_positions_to_read_positions[j + self.ref_start].add(i)
                i += 1
                j += 1
            elif c == 'I':
                b_1 = self.read_seq[i]
                b_2 = '-'
                self.read_error_positions.add(i)
                self.ref_error_positions.add(j + self.ref_start - 1)
                self.ref_error_positions.add(j + self.ref_start)
                self.read_positions_to_ref_positions[i].add(j + self.ref_start - 1)
                self.read_positions_to_ref_positions[i].add(j + self.ref_start)
                diff = '*'
                i += 1
            elif c == 'D':
                b_1 = '-'
                b_2 = ref_seq[j]
                self.read_error_positions.add(i - 1)
                self.read_error_positions.add(i)
                self.ref_error_positions.add(j + self.ref_start)
                self.ref_positions_to_read_positions[j + self.ref_start].add(i - 1)
                self.ref_positions_to_read_positions[j + self.ref_start].add(i)
                diff = '*'
                j += 1
            else:
                assert False
            self.aligned_read_seq.append(b_1)
            self.aligned_ref_seq.append(b_2)
            self.diffs.append(diff)
        assert i == len(self.read_seq)
        assert j == len(ref_seq)
        self.aligned_read_seq = ''.join(self.aligned_read_seq)
        self.aligned_ref_seq = ''.join(self.aligned_ref_seq)
        self.diffs = ''.join(self.diffs)
        assert self.aligned_read_seq.replace('-', '') == self.read_seq
        assert self.aligned_ref_seq.replace('-', '') == ref_seq

        # Convert from dictionary of sets to dictionary of lists.
        self.read_positions_to_ref_positions = \
            {read_pos: sorted(ref_pos)
             for read_pos, ref_pos in self.read_positions_to_ref_positions.items()}
        self.ref_positions_to_read_positions = \
            {ref_pos: sorted(read_pos)
             for ref_pos, read_pos in self.ref_positions_to_read_positions.items()}

        # Sanity check: if there are too many mismatches, something has gone terribly wrong!
        assert len(self.read_error_positions) < len(self.read_seq) // 2

    def print_detailed_alignment_info(self):
        """
        For debug purposes.
        """
        log(self)
        log(self.aligned_read_seq)
        log(self.aligned_ref_seq)
        log(self.diffs)


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


def flip_positions(positions, seq_length):
    """
    This function takes a set of positions (0-based) and flips them to the reverse strand positions.
    """
    flipped_positions = set()
    for p in positions:
        flipped_positions.add(seq_length - p - 1)
    return flipped_positions
