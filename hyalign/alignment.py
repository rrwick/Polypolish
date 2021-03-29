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
import sys

from .log import log, section_header, explanation, quit_with_error
from .misc import run_command, iterate_fastq


def align_reads(target, short1, short2, temp_dir, threads):
    section_header('Aligning short reads to target sequence')
    explanation('Hyalign uses Bowtie2 to align the short reads to the target sequence. The '
                'alignment is done in an unpaired end-to-end manner, and all alignments are kept '
                'for each read.')
    index = build_bowtie2_index(target, temp_dir)
    read_filename, read_count = combine_reads_into_one_file(short1, short2, temp_dir)
    alignments = align_with_bowtie2(read_filename, index, threads)
    print_alignment_info(alignments, read_count)
    return alignments


def build_bowtie2_index(target, temp_dir):
    log(f'Building a Bowtie2 index for {target}:')
    index = str(temp_dir / 'index')
    stdout, stderr, return_code = run_command(['bowtie2-build', target, index])
    if return_code != 0:
        log('\n' + stderr)
        quit_with_error(f'Error: bowtie2-build command failed on {target}')
    for f in sorted(temp_dir.glob('*')):
        log(f'  {f}')
    log()
    return index


def combine_reads_into_one_file(short1, short2, temp_dir):
    log(f'Combining paired reads into one file:')
    read_filename = str(temp_dir / 'reads.fastq')
    read_count = 0
    with open(read_filename, 'wt') as r:
        for _, header, sequence, qualities in iterate_fastq(short1):
            r.write(f'{header}\n{sequence}\n+\n{qualities}\n')
            read_count += 1
        for _, header, sequence, qualities in iterate_fastq(short2):
            r.write(f'{header}\n{sequence}\n+\n{qualities}\n')
            read_count += 1
    log(f'  {read_count:,} reads')
    log(f'  {read_filename}')
    log()
    return read_filename, read_count


def align_with_bowtie2(reads, index, threads):
    log(f'Aligning reads with Bowtie2:')
    command = ['bowtie2', '-U', reads, '-x', index, '-a', '--threads', str(threads), '--end-to-end']
    log(' '.join(command))
    stdout, stderr, return_code = run_command(command)
    alignments = collections.defaultdict(list)
    for line in stdout.splitlines():
        if line.startswith('@'):
            continue
        alignment = Alignment(line)
        alignments[alignment.read_name].append(alignment)
    log()
    return alignments


def print_alignment_info(alignments, read_count):
    sam_count, alignment_count, zero_count, single_count, multi_count = 0, 0, 0, 0, 0
    for read, read_alignments in alignments.items():
        sam_count += len(read_alignments)
        for a in read_alignments:
            if not a.has_flag(4):  # read is aligned
                alignment_count += 1
        if len(read_alignments) == 1:
            if read_alignments[0].is_aligned():
                single_count += 1
            else:
                zero_count += 1
        else:
            multi_count += 1
    log(f'  {sam_count:,} SAM lines')
    number_size = len(f'{sam_count:,}')
    log(f'  {f"{alignment_count:,}".rjust(number_size)} alignments')
    log(f'  {f"{zero_count:,}".rjust(number_size)} reads have no alignments')
    log(f'  {f"{single_count:,}".rjust(number_size)} reads have one alignment')
    log(f'  {f"{multi_count:,}".rjust(number_size)} reads have multiple alignments')
    log()
    assert zero_count + single_count + multi_count == read_count


class Alignment(object):

    def __init__(self, sam_line):
        self.sam_line = sam_line.strip()
        parts = self.sam_line.split('\t')
        if len(parts) < 11:
            sys.exit('\nError: alignment file does not seem to be in SAM format')
        self.sam_line = sam_line
        self.read_name = parts[0]
        self.flags = int(parts[1])
        self.ref_name = parts[2]
        self.ref_start = int(parts[3]) - 1  # switch from SAM's 1-based to Python's 0-based
        self.map_q = int(parts[4])
        self.cigar = parts[5]
        self.ref_end = get_ref_end(self.ref_start, self.cigar)

    def has_flag(self, flag: int):
        return bool(self.flags & flag)

    def is_aligned(self):
        return not self.has_flag(4)

    def is_on_reverse_strand(self):
        return self.has_flag(16)

    def is_on_forward_strand(self):
        return not self.has_flag(16)

    def has_no_indels(self):
        return self.cigar.endswith('M') and self.cigar[:-1].isdigit()


def get_ref_end(ref_start, cigar):
    ref_end = ref_start
    cigar_parts = re.findall(r'\d+[MIDNSHP=X]', cigar)
    for p in cigar_parts:
        size = int(p[:-1])
        letter = p[-1]
        if letter == 'M' or letter == 'D' or letter == 'N' or letter == '=' or letter == 'X':
            ref_end += size
    return ref_end
