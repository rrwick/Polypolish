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

from .log import log, section_header, explanation
from .misc import load_fasta


def mask_target_sequences(read_pair_names, alignments, target_filename, kmer_size):
    section_header('Masking target sequences')
    explanation('Poligner masks out any target bases which appear to be in error. It does this '
                'tallying up all of the k-mers in the aligned reads and finding positions in the '
                'repetitive parts of the target sequence which have unusually low k-mer depth.')

    target_seqs = load_fasta(target_filename)
    mask_positions = {name: set() for name, _ in target_seqs}
    for target_name, target_seq in target_seqs:
        mask_one_target(target_name, target_seq, mask_positions, alignments, kmer_size)
    quit()  # TEMP
    return mask_positions


def mask_one_target(target_name, target_seq, mask_positions, alignments, kmer_size):
    log(f'Masking {target_name}')
    log(f'  {len(target_seq):,} bp total')
    kmer_counts = count_all_read_kmers(alignments, kmer_size, target_name)
    mean_depth = sum(kmer_counts.values()) / len(target_seq)
    log(f'  Mean k-mer depth:      {mean_depth:.1f}x')
    threshold = math.ceil(mean_depth / 10)
    log(f'  Threshold k-mer depth: {threshold:.1f}x')
    depths_by_pos = get_depths_by_pos(target_seq, kmer_counts, kmer_size)
    potential_problem_sites = get_potential_problem_sites(depths_by_pos, kmer_size, threshold,
                                                          target_seq)
    print(potential_problem_sites)  # TEMP

    problem_sites = []
    target_seq_list = [b for b in target_seq]
    for i in potential_problem_sites:
        if improvement_possible(target_seq_list, i, kmer_size, kmer_counts, threshold):
            problem_sites.append(i)

    print(problem_sites)  # TEMP


def count_all_read_kmers(alignments, kmer_size, target_name):
    """
    Builds a dictionary:
    * key =   k-mer sequence
    * value = the number of times that k-mer occurs in reads which aligned to this target sequence
    """
    kmer_counts = collections.defaultdict(int)
    for _, alignments in alignments.items():
        if not alignments:
            continue
        if not any(a.ref_name == target_name for a in alignments):
            continue
        a = alignments[0]  # just process one alignment per read
        for i in range(len(a.read_seq) - kmer_size + 1):
            kmer = a.read_seq[i:i+kmer_size]
            kmer_counts[kmer] += 1
    return kmer_counts


def get_depths_by_pos(target_seq, kmer_counts, kmer_size):
    """
    Builds a dictionary:
    * key =   target sequence position
    * value = a list of k-mer counts for each k-mer overlapping that position
    """
    depths_by_pos = {i: [] for i in range(len(target_seq))}
    for i in range(len(target_seq) - kmer_size + 1):
        kmer = target_seq[i:i+kmer_size]
        depth = kmer_counts[kmer]
        for j in range(i, i+kmer_size):
            depths_by_pos[j].append(depth)
    return depths_by_pos


def get_potential_problem_sites(depths_by_pos, kmer_size, threshold, target_seq):
    """
    Returns a list of positions in the target sequence which look like they may be problems. They
    qualify if any of the following is true:
    * all of the k-mer counts at that position fail to meet the threshold
    * all of the k-mer counts at that position except the first fail to meet the threshold
    * most of the k-mer counts at that position fail to meet the threshold and it doesn't look any
      than its neighbouring positions
    """
    potential_problem_sites = []
    for i, depths in depths_by_pos.items():
        problem = False
        if len(depths) == kmer_size:
            if all(d < threshold for d in depths[1:]):
                problem = True
                potential_problem_sites.append(i)
            elif 0 < i < len(target_seq) - 1:
                this_fail_count = sum(1 if d < threshold else 0 for d in depths)
                if this_fail_count > kmer_size / 2:
                    prev_fail_count = sum(1 if d < threshold else 0 for d in depths_by_pos[i-1])
                    next_fail_count = sum(1 if d < threshold else 0 for d in depths_by_pos[i+1])
                    if this_fail_count >= prev_fail_count and this_fail_count >= next_fail_count:
                        potential_problem_sites.append(i)
                        problem = True
        print(i, target_seq[i], depths, problem)  # TEMP
    return potential_problem_sites


def improvement_possible(target_seq_list, i, kmer_size, kmer_counts, threshold):
    """
    This function tries to improve on the target sequence at the given position by mutating the
    base at that position in all possible ways and seeing if any of the mutations give better k-mer
    counts.
    """
    start_pos = max(0, i - kmer_size)
    end_pos = min(i + kmer_size + 1, len(target_seq_list))
    preceding_seq = ''.join(target_seq_list[start_pos:i])
    following_seq = ''.join(target_seq_list[i+1:end_pos])
    unmutated_base = target_seq_list[i]
    unmutated_failures = count_kmer_failures(preceding_seq + unmutated_base + following_seq,
                                             kmer_size, kmer_counts, threshold)
    print(i, unmutated_base, unmutated_failures)  # TEMP

    best_base = unmutated_base
    best_failures = unmutated_failures
    for mutated_base in get_possible_mutations(unmutated_base):
        mutated_failures = count_kmer_failures(preceding_seq + mutated_base + following_seq,
                                               kmer_size, kmer_counts, threshold)
        if mutated_failures < best_failures:
            best_base = mutated_base
            best_failures = mutated_failures
    print(i, best_base, best_failures)  # TEMP
    print()  # TEMP


def count_kmer_failures(seq, kmer_size, kmer_counts, threshold):
    total_count, fail_count = 0, 0
    for i in range(len(seq) - kmer_size + 1):
        kmer = seq[i:i + kmer_size]
        total_count += 1
        if kmer_counts[kmer] < threshold:
            fail_count += 1
    return fail_count / total_count


def mean_kmer_count(seq, kmer_size, kmer_counts):
    counts = []
    for i in range(len(seq) - kmer_size + 1):
        kmer = seq[i:i + kmer_size]
        counts.append(kmer_counts[kmer])
    return statistics.mean(counts)


def get_possible_mutations(base):
    last_base = base[-1]
    deletion = base[:-1]
    substitutions = [deletion + b for b in {'A', 'C', 'G', 'T'} - {last_base}]
    insertions = [base + ins for ins in ['A', 'C', 'G', 'T']]
    return sorted(substitutions + [deletion] + insertions)
