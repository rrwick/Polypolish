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

import pathlib
import re
import tempfile

from .log import log, section_header, explanation, quit_with_error
from .misc import run_command, iterate_fastq, reverse_complement


def align_reads(target, short1, short2, threads, max_errors, debug):
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
            combine_reads_into_one_file(short1, short2, temp_dir)
        alignments, header_lines, unaligned = \
            align_with_minimap2(read_filename, read_pair_names, target, threads, debug)
        add_secondary_read_seqs(alignments, read_pair_names)
        filter_alignments(alignments, read_pair_names, max_errors, unaligned)
        print_alignment_info(alignments, read_count, read_pair_names)
        return alignments, read_pair_names, read_count, header_lines, unaligned


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
    return read_filename, total_count, sorted(first_names)


def align_with_minimap2(reads, read_pair_names, target, threads, debug):
    """
    Runs minimap2 to align the reads. We use minimap2 using mostly the -x sr option (short-read
    preset), but with the following changes:
    * I set `--secondary=yes -p0.1 -N1000000` to include all alignments
    * I do NOT use --sr (short-read alignment heuristics), which resulted in some goofy alignments
      in my tests
    """
    log(f'Aligning reads with minimap2:')
    command = ['minimap2', '-a', '-k21', '-w11', '--frag=yes', '-A2', '-B8', '-O12,32', '-E2,1',
               '-r50', '-f1000,5000', '-n2', '-m20', '-s40', '-g200', '-2K50m', '--heap-sort=yes',
               '--secondary=yes', '-p0.1', '-N1000000', '-t', str(threads), target, reads]
    log('  ' + ' '.join(command))
    stdout, stderr, return_code = run_command(command)
    if debug:
        log()
        log(stderr)
        log()
    # TODO: check return code here
    alignments, unaligned = {}, {}
    header_lines = []
    alignment_count = 0
    for name in read_pair_names:
        alignments[name + '/1'] = []
        alignments[name + '/2'] = []
    for line in stdout.splitlines():
        if line.startswith('@'):
            header_lines.append(line)
            continue
        alignment = Alignment(line)
        if alignment.is_aligned():
            alignments[alignment.read_name].append(alignment)
            alignment_count += 1
        else:
            unaligned[alignment.read_name] = alignment

    log(f'  {alignment_count:,} total alignments')
    log()
    return alignments, header_lines, unaligned


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
                    a.set_read_seq(seq)
                    a.set_read_qual(qual)
                else:
                    a.set_read_seq(reverse_complement(seq))
                    a.set_read_qual(qual[::-1])
        for a in alignments[name]:
            assert a.read_seq != '*' and a.read_qual != '*'


def filter_alignments(alignments, read_pair_names, max_errors, unaligned):
    log('Filtering for high quality end-to-end alignments... ', end='')
    read_names = [n + '/1' for n in read_pair_names] + [n + '/2' for n in read_pair_names]
    for name in read_names:
        good_alignments, bad_alignments = [], []
        for a in alignments[name]:
            if a.starts_and_ends_with_match() and a.nm_tag <= max_errors:
                good_alignments.append(a)
            else:
                bad_alignments.append(a)
            alignments[name] = good_alignments
            if not good_alignments:
                unaligned[name] = bad_alignments[0]
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
    for name in read_pair_names:
        name_1, name_2 = name + '/1', name + '/2'
        alignments_1, alignments_2 = alignments[name_1], alignments[name_2]
        if len(alignments_1) > 1:
            multi_alignment_read_names.append(name_1)
        if len(alignments_2) > 1:
            multi_alignment_read_names.append(name_2)
    return multi_alignment_read_names


class Alignment(object):

    def __init__(self, sam_line):
        sam_line = sam_line.strip()
        self.parts = sam_line.split('\t')
        if len(self.parts) < 11:
            quit_with_error('\nError: alignment file does not seem to be in SAM format')

        self.read_name = self.parts[0]
        self.flags = int(self.parts[1])
        self.ref_name = self.parts[2]
        self.ref_start = int(self.parts[3]) - 1  # switch from SAM's 1-based to Python's 0-based
        self.cigar = self.parts[5]
        self.ref_end = get_ref_end(self.ref_start, self.cigar)
        self.read_seq = self.parts[9]
        self.read_qual = self.parts[10]

        self.nm_tag = None
        for p in self.parts:
            if p.startswith('NM:i:'):
                self.nm_tag = int(p[5:])

        # We now change a few parts to make the SAM be in a paired-read format.
        self.parts[0] = self.parts[0][:-2]  # remove '/1' or '/2' from the read name
        assert self.parts[6] == '*'
        self.parts[6] = '='

    def __repr__(self):
        return f'{self.read_name}:{self.ref_name}:{self.ref_start}-{self.ref_end}:{self.nm_tag}'

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

    def set_flags(self, flags):
        self.flags = flags
        self.parts[1] = str(flags)

    def set_ref_name(self, ref_name):
        self.ref_name = ref_name
        self.parts[2] = ref_name

    def set_ref_start(self, ref_start):
        self.ref_start = ref_start
        self.parts[3] = str(ref_start + 1)

    def set_mapq(self, mapq):
        self.parts[4] = str(mapq)

    def set_cigar(self, cigar):
        self.cigar = cigar
        self.parts[5] = cigar

    def set_rnext(self, rnext):
        self.parts[6] = rnext

    def set_pnext(self, pnext):
        self.parts[7] = str(pnext + 1)

    def set_tlen(self, tlen):
        self.parts[8] = str(tlen)

    def set_read_seq(self, seq):
        self.read_seq = seq
        self.parts[9] = seq

    def set_read_qual(self, qual):
        self.read_qual = qual
        self.parts[10] = qual

    def make_unaligned(self):
        self.flags = 5  # read paired (1) and read unmapped (4)
        self.set_mapq(0)
        self.set_cigar('*')

    def get_sam_line(self):
        return '\t'.join(self.parts)

    def has_no_indels(self):
        return self.cigar.endswith('M') and self.cigar[:-1].isdigit()

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


def verify_no_multi_alignments(alignments, read_pair_names):
    for name in read_pair_names:
        name_1, name_2 = name + '/1', name + '/2'
        alignments_1, alignments_2 = alignments[name_1], alignments[name_2]
        assert len(alignments_1) <= 1 and len(alignments_2) <= 1


def fix_sam_pairing(alignments, unaligned, read_pair_names):
    for name in read_pair_names:
        name_1, name_2 = name + '/1', name + '/2'
        alignments_1, alignments_2 = alignments[name_1], alignments[name_2]

        # If neither read in the pair aligned:
        if not alignments_1 and not alignments_2:
            unaligned[name_1].set_ref_name('*')
            unaligned[name_2].set_ref_name('*')
            unaligned[name_1].set_ref_start(-1)
            unaligned[name_2].set_ref_start(-1)
            unaligned[name_1].set_rnext('*')
            unaligned[name_2].set_rnext('*')
            unaligned[name_1].set_pnext(-1)
            unaligned[name_2].set_pnext(-1)
            unaligned[name_1].set_tlen(0)
            unaligned[name_2].set_tlen(0)

        # If only the first read in the pair aligned:
        elif alignments_1 and not alignments_2:
            unaligned[name_2].set_ref_name(alignments_1[0].ref_name)
            unaligned[name_2].set_ref_start(alignments_1[0].ref_start)
            unaligned[name_2].set_pnext(alignments_1[0].ref_start)
            alignments_1[0].set_pnext(alignments_1[0].ref_start)

        # If only the second read in the pair aligned:
        elif not alignments_1 and alignments_2:
            unaligned[name_1].set_ref_name(alignments_2[0].ref_name)
            unaligned[name_1].set_ref_start(alignments_2[0].ref_start)
            unaligned[name_1].set_pnext(alignments_2[0].ref_start)
            alignments_2[0].set_pnext(alignments_2[0].ref_start)

        # If both reads aligned:
        elif alignments_1 and alignments_2:
            a_1, a_2 = alignments_1[0], alignments_2[0]
            a_1.set_pnext(a_2.ref_start)
            a_2.set_pnext(a_1.ref_start)
            template_min = min(a_1.ref_start, a_2.ref_start)
            template_max = max(a_1.ref_end, a_2.ref_end)
            tlen = template_max - template_min
            if a_1.ref_start < a_2.ref_start:
                a_1.set_tlen(tlen)
                a_2.set_tlen(-tlen)
            else:
                a_1.set_tlen(-tlen)
                a_2.set_tlen(tlen)

        else:
            assert False


def output_alignments_to_stdout(alignments, read_pair_names, header_lines, unaligned):
    for line in header_lines:
        print(line)
    for name in read_pair_names:
        name_1, name_2 = name + '/1', name + '/2'
        alignments_1, alignments_2 = alignments[name_1], alignments[name_2]
        for a in alignments_1:
            print(a.get_sam_line())
        if not alignments_1:
            print(unaligned[name_1].get_sam_line())
        for a in alignments_2:
            print(a.get_sam_line())
        if not alignments_2:
            print(unaligned[name_2].get_sam_line())
