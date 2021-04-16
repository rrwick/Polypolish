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

from .alignment import get_multi_alignment_read_names, print_alignment_info
from .log import log, section_header, explanation


def mask_target_sequences(read_pair_names, alignments, target_seqs, debug):
    section_header('Masking target sequences')
    explanation('Poligner masks out any target bases which appear to be in error. It does this '
                'tallying up all of the k-mers in the aligned reads and finding positions in the '
                'repetitive parts of the target sequence which have unusually low k-mer depth.')

    mask_positions = {}
    for target_name, target_seq in target_seqs:
        mask_positions_one_target = \
            get_mask_positions(target_name, target_seq, read_pair_names, alignments, debug)
        mask_positions[target_name] = mask_positions_one_target
    return mask_positions


def get_mask_positions(target_name, target_seq, read_pair_names, alignments, debug):
    log(f'Masking {target_name}')
    log(f'  {len(target_seq):,} bp total')
    depths_by_pos = get_depths_by_pos(target_name, target_seq, alignments)
    mask_positions = set()

    read_names = [n + '/1' for n in read_pair_names] + [n + '/2' for n in read_pair_names]
    pileup = get_pileup(read_names, alignments, target_name, target_seq)
    if debug:
        log()

    for i in range(len(target_seq)):
        read_bases = pileup[i]
        depth = depths_by_pos[i]
        ref_base = target_seq[i]
        match_count = sum(1 if b == ref_base else 0 for b in read_bases)
        if depth > 0.0:
            match_fraction = match_count / depth
            bad_position = match_fraction < 0.5
            if bad_position:
                mask_positions.add(i)
            if debug:
                result = 'FAIL' if bad_position else 'pass'
                debug_str = f'{i}  {ref_base}  {depth}  {match_fraction:.2f}  {result}  ' \
                    f'{" ".join(read_bases)}'
                log(debug_str)
    mask_count = len(mask_positions)
    mask_percent = 100.0 * mask_count / len(target_seq)
    if debug:
        log('\nMasked positions:')
        log(', '.join([str(i) for i in sorted(mask_positions)]))
        log()
    log(f'  {mask_count:,} positions masked ({mask_percent:.3f}% of total positions)')
    log()
    return mask_positions


def get_depths_by_pos(target_name, target_seq, alignments):
    depths_by_pos = {i: 0.0 for i in range(len(target_seq))}
    for _, read_alignments in alignments.items():
        if not read_alignments:
            continue
        depth_contribution = 1.0 / len(read_alignments)
        for a in read_alignments:
            if a.ref_name == target_name:
                for i in range(a.ref_start, a.ref_end):
                    depths_by_pos[i] += depth_contribution
    return depths_by_pos


def get_pileup(read_names, alignments, target_name, target_seq):
    pileup = collections.defaultdict(list)
    for name in read_names:
        read_alignments = alignments[name]
        for a in read_alignments:
            if a.ref_name != target_name:
                continue
            ref_seq = target_seq[a.ref_start:a.ref_end]
            aligned_bases = a.get_read_bases_for_each_target_base(ref_seq)
            for i, bases in enumerate(aligned_bases):
                ref_pos = a.ref_start + i
                if a.ref_positions_to_read_positions is not None:
                    corresponding_read_bases = a.ref_positions_to_read_positions[ref_pos]
                else:
                    corresponding_read_bases = set()
                if not any(b in a.masked_read_positions for b in corresponding_read_bases):
                    pileup[ref_pos].append(bases)
    return pileup


def select_best_alignments(alignments, mask_positions, read_pair_names, read_count, target_seqs):
    section_header('Selecting best alignments')
    explanation('Poligner now chooses the best alignment(s) for each read, ignoring the masked '
                'positions of the target sequence.')
    target_seqs = dict(target_seqs)
    multi_alignment_read_names = get_multi_alignment_read_names(read_pair_names, alignments)
    for name in multi_alignment_read_names:
        best_alignments = \
            choose_best_alignments_one_read(alignments[name], mask_positions, target_seqs)
        alignments[name] = best_alignments
    print_alignment_info(alignments, read_count, read_pair_names)


def choose_best_alignments_one_read(read_alignments, mask_positions, target_seqs):
    """
    This function chooses the best alignment(s) for a multi-alignment read. Masked positions of the
    target sequence are ignored. If a tie occurs, then all best positions are returned.
    """
    best_error_count = None
    best_alignments = []
    for a in read_alignments:
        ref_seq = target_seqs[a.ref_name][a.ref_start:a.ref_end]
        ref_range = list(range(a.ref_start, a.ref_end))
        aligned_bases = a.get_read_bases_for_each_target_base(ref_seq)
        assert len(ref_seq) == len(ref_range) == len(aligned_bases)
        error_count = 0
        for i, ref_base, read_base in zip(ref_range, ref_seq, aligned_bases):
            if ref_base != read_base and i not in mask_positions[a.ref_name]:
                error_count += 1
        if best_error_count is None or error_count < best_error_count:  # new best
            best_error_count = error_count
            best_alignments = [a]
        elif error_count == best_error_count:  # tie for best
            best_alignments.append(a)
    return best_alignments
