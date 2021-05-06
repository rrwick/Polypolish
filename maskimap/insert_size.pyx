"""
Copyright 2021 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Maskimap

This file is part of Maskimap. Maskimap is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Maskimap is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Maskimap.
If not, see <http://www.gnu.org/licenses/>.
"""

import random

from .alignment import print_alignment_info
from .log import log, section_header, explanation
from .misc import get_percentile_sorted
from . import settings


def get_insert_size_distribution(alignments):
    section_header('Finding insert size distribution')
    explanation('Cleanly aligned read pairs (one alignment for the first read and one alignment '
                'for the second read on opposite strands) are used to determine the insert size '
                'distribution for the read set.')
    unique_pairs = get_uniquely_aligned_pairs(alignments)
    proper_pairs = get_properly_aligned_pairs(unique_pairs)
    log(f'{len(proper_pairs):,} read pairs are cleanly aligned (no indels and good orientation)')
    log()
    insert_sizes = []
    for alignment_1, alignment_2 in proper_pairs:
        insert_size = get_insert_size(alignment_1, alignment_2)
        if insert_size <= settings.MAX_ALLOWED_INSERT_SIZE:
            insert_sizes.append(insert_size)
    return get_distribution(insert_sizes)


def get_uniquely_aligned_pairs(alignments):
    unique_pairs = []
    for name_1, alignments_1 in alignments.items():
        if name_1.endswith('/1') and len(alignments_1) == 1:
            name_2 = name_1[:-2] + '/2'
            alignments_2 = alignments[name_2]
            if len(alignments_2) == 1:
                unique_pairs.append((alignments_1[0], alignments_2[0]))
    return unique_pairs


def get_properly_aligned_pairs(unique_pairs):
    proper_pairs = []
    for alignment_1, alignment_2 in unique_pairs:
        if alignment_1.has_no_indels() and alignment_2.has_no_indels():
            if alignment_1.is_on_forward_strand() and alignment_2.is_on_reverse_strand():
                proper_pairs.append((alignment_1, alignment_2))
            elif alignment_1.is_on_reverse_strand() and alignment_2.is_on_forward_strand():
                proper_pairs.append((alignment_1, alignment_2))
    return proper_pairs


def get_insert_size(alignment_1, alignment_2):
    starts_ends = sorted([alignment_1.ref_start, alignment_1.ref_end,
                          alignment_2.ref_start, alignment_2.ref_end])
    return starts_ends[-1] - starts_ends[0]


def get_distribution(insert_sizes):
    insert_sizes = sorted(insert_sizes)
    distribution = {}
    log('Percentile   Insert size')
    for p in [0.001, 0.01, 0.1, 1, 10, 50, 90, 99, 99.9, 99.99, 99.999]:
        distribution[p] = get_percentile_sorted(insert_sizes, p)
        log(f'    {p:6.3f}        {distribution[p]:6}')
    log()
    return distribution


def select_alignments_using_insert_size(alignments, distribution, read_pair_names, read_count):
    """
    This function modifies the alignments dictionary, removing alignments that are not part of a
    pair with a good insert size.
    """
    section_header('Selecting alignments using insert size')
    explanation('Maskimap now filters alignments using insert size. Whenever there is a read pair '
                'with multiple possible combinations, Maskimap will discard any alignments that '
                'do not appear to be part of a combination with a valid insert size.')

    for name in read_pair_names:
        name_1, name_2 = name + '/1', name + '/2'
        alignments_1, alignments_2 = alignments[name_1], alignments[name_2]
        count_1, count_2 = len(alignments_1), len(alignments_2)
        if count_1 >= 1 and count_2 >= 1 and count_1+count_2 >= 3:
            all_scores = []
            for a_1 in alignments_1:
                for a_2 in alignments_2:
                    all_scores.append(score_insert_size(get_insert_size(a_1, a_2), distribution))
            min_score = max(all_scores) - 2
            good_1 = select_good_alignments(alignments_1, alignments_2, distribution, min_score)
            good_2 = select_good_alignments(alignments_2, alignments_1, distribution, min_score)
            alignments[name_1] = good_1
            alignments[name_2] = good_2

    print_alignment_info(alignments, read_count, read_pair_names)


def select_good_alignments(alignments_1, alignments_2, distribution, min_score):
    """
    This function looks at all pairwise combinations of alignments from the two groups and returns
    a list of alignments from the first group which seem to be in good pairs.
    """
    good_alignments = []
    for a_1 in alignments_1:
        for a_2 in alignments_2:
            insert_size = get_insert_size(a_1, a_2)
            score = score_insert_size(insert_size, distribution)
            if score >= min_score:
                good_alignments.append(a_1)
                break
    return good_alignments


def score_insert_size(insert_size, distribution):
    """
    Returns a score for the insert size, with higher values being better. In this context, 'better'
    means closer to a typical insert size based on the empirical distribution.
    """
    if distribution[10.0] <= insert_size <= distribution[90.0]:
        return 5
    if distribution[1.0] <= insert_size <= distribution[99.0]:
        return 4
    if distribution[0.1] <= insert_size <= distribution[99.9]:
        return 3
    if distribution[0.01] <= insert_size <= distribution[99.99]:
        return 2
    if distribution[0.001] <= insert_size <= distribution[99.999]:
        return 1
    return 0


def final_alignment_selection(alignments, distribution, read_pair_names, read_count):

    section_header('Final alignment selection')
    explanation('After the last round of selection, all reads with multiple alignments have a '
                'tie between equally good alignments. Maskimap now makes a final selection, '
                'keeping only one alignment for each multi-alignment read. This final selection '
                'is made using insert size, if possible. Otherwise the kept alignment is chosen '
                'at random.')

    insert_decision_count, random_decision_count = 0, 0
    for name in read_pair_names:
        name_1, name_2 = name + '/1', name + '/2'
        alignments_1, alignments_2 = alignments[name_1], alignments[name_2]
        count_1, count_2 = len(alignments_1), len(alignments_2)

        if count_1 <= 1 and count_2 <= 1:
            pass

        elif count_1 > 1 and count_2 == 0:
            alignments[name_1] = [random.choice(alignments[name_1])]

        elif count_1 == 0 and count_2 > 1:
            alignments[name_2] = [random.choice(alignments[name_2])]

        elif count_1 >= 1 and count_2 >= 1 and count_1+count_2 >= 3:
            all_scores = []
            for a_1 in alignments_1:
                for a_2 in alignments_2:
                    insert_size = get_insert_size(a_1, a_2)
                    all_scores.append(score_insert_size(insert_size, distribution))
            max_score = max(all_scores)
            good_pairs = []
            for a_1 in alignments_1:
                for a_2 in alignments_2:
                    insert_size = get_insert_size(a_1, a_2)
                    if score_insert_size(insert_size, distribution) == max_score:
                        good_pairs.append((a_1, a_2))
            assert len(good_pairs) >= 1
            if len(good_pairs) == 1:
                insert_decision_count += 1
            else:
                random_decision_count += 1
            a_1, a_2 = random.choice(good_pairs)
            alignments[name_1] = [a_1]
            alignments[name_2] = [a_2]

        else:
            assert False
    log(f'Ties broken with insert size:   {insert_decision_count:,}')
    log(f'Ties broken with random choice: {random_decision_count:,}')
    log()
    print_alignment_info(alignments, read_count, read_pair_names)


def set_sam_flags(alignments, unaligned, read_pair_names, insert_size_distribution):
    for a in unaligned.values():
        a.make_unaligned()
    read_names = [n + '/1' for n in read_pair_names] + [n + '/2' for n in read_pair_names]
    for name in read_names:
        if not alignments[name]:
            a = unaligned[name]
        else:
            assert len(alignments[name]) == 1
            a = alignments[name][0]
        new_flags = 1                                             # 1 = read paired
        if not a.is_aligned():
            new_flags += 4                                        # 4 = read unmapped
        if name.endswith('/1'):
            pair_alignments = alignments[name[:-2] + '/2']
            new_flags += 64                                       # 64 = first in pair
        elif name.endswith('/2'):
            pair_alignments = alignments[name[:-2] + '/1']
            new_flags += 128                                      # 128 = second in pair
        else:
            assert False
        if a.is_on_reverse_strand():
            new_flags += 16                                       # 16 = read reverse strand
        if not pair_alignments:
            new_flags += 8                                        # 8 = mate unmapped
        else:
            pair_a = pair_alignments[0]
            insert_size = get_insert_size(a, pair_a)
            if score_insert_size(insert_size, insert_size_distribution) >= 3:
                new_flags += 2                                    # 2 = read mapped in proper pair
            if pair_a.is_on_reverse_strand():
                new_flags += 32                                   # 32 = mate reverse strand
        a.set_flags(new_flags)
