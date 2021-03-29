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

from .log import log, section_header, explanation
from .misc import get_percentile_sorted
from . import settings


def get_insert_size_distribution(alignments):
    section_header('Finding insert size distribution')
    explanation('Uniquely aligned read pairs (one alignment for the first read and one alignment '
                'for the second read) are used to determine the insert size distribution.')
    unique_pairs = get_uniquely_aligned_pairs(alignments)
    log(f'{len(unique_pairs):,} pairs are uniquely aligned')
    proper_pairs = get_properly_aligned_pairs(unique_pairs)
    log(f'{len(proper_pairs):,} pairs are cleanly aligned (no indels and good orientation)')
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
        if name_1.endswith('/1') and len(alignments_1) == 1 and alignments_1[0].is_aligned():
            name_2 = name_1[:-2] + '/2'
            alignments_2 = alignments[name_2]
            if len(alignments_2) == 1 and alignments_2[0].is_aligned():
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
    log('Percentile   Insert')
    for p in [0.001, 0.01, 0.1, 1, 10, 50, 90, 99, 99.9, 99.99, 99.999]:
        distribution[p] = get_percentile_sorted(insert_sizes, p)
        log(f'    {p:6.3f}   {distribution[p]:6}')
    return distribution
