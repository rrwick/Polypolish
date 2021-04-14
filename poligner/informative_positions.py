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
import math
import statistics

from .alignment import get_multi_alignment_read_names, get_expanded_cigar, print_alignment_info, \
    get_read_bases_for_each_target_base
from .log import log, section_header, explanation
from .misc import load_fasta


def find_informative_positions(read_pair_names, alignments, target):
    section_header('Finding informative target positions')
    explanation('Poligner now finds bases in the target sequence which have significant ambiguity. '
                'E.g. if a target base is "A" and all read alignments agree with that base, then '
                'that position is not informative. But if the read alignments show ambiguity '
                '(e.g. some say "A" and some say "T"), then that position is informative.')

    target_seqs = load_fasta(target)
    repetitive_regions = find_repetitive_regions(target_seqs, read_pair_names, alignments)
    informative_positions = {}
    for target_name, target_seq in target_seqs:
        informative_positions[target_name] = \
            find_positions(target_name, target_seq, alignments, read_pair_names,
                           repetitive_regions[target_name])
    return informative_positions


def find_repetitive_regions(target_seqs, read_pair_names, alignments):
    repetitive_regions = {name: set() for name, _ in target_seqs}
    multi_alignment_read_names = get_multi_alignment_read_names(read_pair_names, alignments)
    for name in multi_alignment_read_names:
        read_alignments = alignments[name]
        for a in read_alignments:
            for i in range(a.ref_start, a.ref_end):
                repetitive_regions[a.ref_name].add(i)
    return repetitive_regions


def find_positions(target_name, target_seq, alignments, read_pair_names, repetitive_regions):
    log(f'Analysing {target_name}:')
    log(f'  {len(target_seq):,} bp total')
    repeat_length = len(repetitive_regions)
    repeat_percent = 100.0 * repeat_length / len(target_seq)
    log(f'  {repeat_length:,} bp repetitive ({repeat_percent:.2f}%)')

    pileup = get_pileup(read_pair_names, alignments, target_name, target_seq)
    non_repeat_depth = get_mean_non_repeat_depth(pileup, repetitive_regions)
    non_ref_base_threshold = int(math.ceil(non_repeat_depth / 10))
    log(f'  Mean non-repeat depth: {non_repeat_depth:.1f}')
    log(f'  Minimum ambiguity threshold: {non_ref_base_threshold}')

    informative_positions = set()
    log('  Finding informative positions: ', end='')
    for i in sorted(pileup.keys()):
        if i in repetitive_regions:
            # pileup_str = f'{i}  {" ".join(sorted(pileup[i]))}'  # DEBUGGING
            counts = [pileup[i].count(b) for b in set(pileup[i])]
            if len(counts) == 1:
                second_best_count = 0
            else:
                second_best_count = sorted(counts)[-2]
            if second_best_count >= non_ref_base_threshold:
                informative_positions.add(i)
            #     log(bold_yellow(pileup_str))  # DEBUGGING
            # else:  # DEBUGGING
            #     log(pileup_str)  # DEBUGGING
    log(f'{len(informative_positions):,} positions found')
    log()
    return informative_positions


def get_pileup(read_pair_names, alignments, target_name, target_seq):
    pileup = collections.defaultdict(list)
    read_names = [n + '/1' for n in read_pair_names] + [n + '/2' for n in read_pair_names]
    for name in read_names:
        read_alignments = alignments[name]
        for a in read_alignments:
            if a.ref_name != target_name:
                continue
            ref_seq = target_seq[a.ref_start:a.ref_end]
            aligned_bases = a.get_read_bases_for_each_target_base(ref_seq)
            for i, bases in enumerate(aligned_bases):
                ref_pos = a.ref_start + i
                if 'N' not in bases:
                    pileup[ref_pos].append(bases)
    return pileup


def get_mean_non_repeat_depth(pileup, repetitive_regions):
    non_repeat_depths = []
    for ref_pos in sorted(pileup.keys()):
        if ref_pos not in repetitive_regions:
            non_repeat_depths.append(len(pileup[ref_pos]))
    mean_depth = statistics.mean(non_repeat_depths)
    return mean_depth


def alignment_overlaps_repeats(a, repetitive_regions):
    for i in range(a.ref_start, a.ref_end):
        if i in repetitive_regions:
            return True
    return False


def select_alignments_using_informative_positions(alignments, informative_positions,
                                                  read_pair_names, read_count):
    section_header('Selecting alignments informative target positions')
    explanation('Poligner now filters alignments using informative positions in '
                'the target sequence.')
    multi_alignment_read_names = get_multi_alignment_read_names(read_pair_names, alignments)
    for name in multi_alignment_read_names:
        best_alignments = choose_best_alignments_one_read(alignments[name], informative_positions)
        alignments[name] = best_alignments
    print_alignment_info(alignments, read_count, read_pair_names)


def choose_best_alignments_one_read(read_alignments, informative_target_positions):
    """
    This function chooses the best alignment for a multi-alignment read. Informative positions are
    used to mask the read sequence, if available. If a tie occurs, then a best position is chosen
    at random.
    """
    informative_read_positions = set()
    read_length = None
    for a in read_alignments:
        # a.print_detailed_alignment_info()  # DEBUGGING
        for i in range(a.ref_start, a.ref_end):
            if i in informative_target_positions[a.ref_name]:
                for read_pos in a.ref_positions_to_read_positions[i]:
                    informative_read_positions.add(read_pos)
        #             log(f'{i}  {read_pos}')  # DEBUGGING
        # log()  # DEBUGGING
        if read_length is None:
            read_length = len(a.read_seq)
        else:
            assert read_length == len(a.read_seq)
    # log(f'{informative_read_positions}')  # DEBUGGING

    best_error_count = None
    best_alignments = []
    for a in read_alignments:
        error_count = 0
        for i in a.read_error_positions:
            if i in informative_read_positions:
                error_count += 1
        if best_error_count is None or error_count < best_error_count:
            best_error_count = error_count
            best_alignments = [a]
        elif error_count == best_error_count:
            best_alignments.append(a)
    # log(f'{len(best_alignments)} best alignments\n\n\n\n\n')  # DEBUGGING

    return best_alignments
