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
import statistics

from .alignment import get_multi_alignment_read_names, print_alignment_info
from .log import log, section_header, explanation


def polish_target_sequences(alignments, assembly_seqs, debug):
    section_header('Masking assembly sequences')
    explanation('Polypolish masks out any target bases which appear to be in error. It does this '
                'by generating a pileup using all of the alignments and looking for positions '
                'where the matching pileup bases (read bases which match the reference) total '
                'less than half of the expected depth.')

    mask_positions = {}
    for target_name, target_seq in assembly_seqs:
        polish_target_sequence(target_name, target_seq, alignments, debug)


def polish_target_sequence(target_name, target_seq, alignments, debug):
    log(f'Polishing {target_name}')
    log(f'  {len(target_seq):,} bp total')

    changed_positions = set()
    pileup, depths_by_pos = get_pileup(alignments, target_name, target_seq)
    mean_depth = statistics.mean(depths_by_pos.values())
    log(f'  mean read depth: {mean_depth:.1f}x')
    uncovered = sum(1 if d == 0.0 else 0 for d in depths_by_pos.values())
    coverage = 100.0 * (len(target_seq) - uncovered) / len(target_seq)
    have = 'has' if uncovered == 1 else 'have'
    log(f'  {uncovered:,} bp {have} a depth of zero ({coverage:.3f}% coverage)')

    if debug:
        log()
        log('Column_1=ref_pos  Column_2=ref_base  Column_3=depth  Column_4=match_fraction  '
            'Column_5=read_pileup')

    new_bases = []
    for i in range(len(target_seq)):
        read_base_counts = pileup[i]
        depth = depths_by_pos[i]
        ref_base = target_seq[i]

        if depth > 0.0:
            target_count = 0.5 * depth
            valid_bases = [b for b, c in read_base_counts.items() if c >= target_count]
        else:
            valid_bases = []

        if len(valid_bases) == 0:
            new_base = ref_base
            status = 'no_reads'
        elif len(valid_bases) == 1:
            new_base = valid_bases[0]
            if new_base == ref_base:
                status = 'kept'
            else:
                changed_positions.add(i)
                status = 'changed'
        else:  # multiple valid bases
            new_base = ref_base
            status = 'multiple'

        if debug:
            count_str = ','.join([f'{b}x{c}' for b, c in read_base_counts.items()])
            debug_str = f'  {i}  {ref_base}  {depth:.1f}  {count_str}  {new_base}  {status}'
            log(debug_str)

        new_bases.append(new_base)

    changed_count = len(changed_positions)
    changed_percent = 100.0 * changed_count / len(target_seq)
    estimated_accuracy = 100.0 - changed_percent

    positions = 'position' if changed_count == 1 else 'positions'
    log(f'  {changed_count:,} {positions} changed ({changed_percent:.3f}% of total positions)')
    log(f'  estimated assembly sequence accuracy: {estimated_accuracy:.3f}%')
    log()

    new_sequence = ''.join(new_bases)
    new_sequence = new_sequence.replace('-', '')

    print(f'>{target_name}_polished')
    print(new_sequence)


def print_match_fraction_distribution(match_fraction_distribution):
    step_size = 0.01
    inverse_step = int(round(1.0 / step_size))
    log('Distribution of match fractions over the target seq:')
    histogram, more_than_one = collections.defaultdict(int), 0
    for fraction, count in match_fraction_distribution.items():
        if fraction > 1.0:
            more_than_one += count
        else:
            rounded_frac = int(round(fraction * inverse_step))
            histogram[rounded_frac] += count
    for fraction in sorted(histogram.keys()):
        count = histogram[fraction]
        fraction /= inverse_step
        log(f'   {fraction:.2f}: {count} bp')
    log(f'  >1.00: {more_than_one} bp')
    log()


def get_pileup(alignments, target_name, target_seq):
    pileup = {i: {} for i in range(len(target_seq))}
    depths_by_pos = {i: 0.0 for i in range(len(target_seq))}
    for _, read_alignments in alignments.items():
        for a in read_alignments:
            if a.ref_name != target_name:
                continue
            ref_seq = target_seq[a.ref_start:a.ref_end]
            aligned_bases = a.get_read_bases_for_each_target_base(ref_seq)
            depth_contribution = 1.0 / len(read_alignments)

            # Alignments that end on in a homopolymer can cause trouble, as they can align cleanly
            # (without an indel) even when an indel is needed.
            #
            # For example, an alignment should look like this:
            #   read: ... T G A G T A C AG G G G G A A G T
            #   ref:  ... T G A G T A C A  G G G G A A G T C C A G T ...
            #
            # But if the read ends in the homopolymer, it could look like this:
            #   read: ... T G A G T A C A G G
            #   ref:  ... T G A G T A C A G G G G A A G T C C A G T ...
            #
            # Which results in a clean alignment on the 'A' that should be 'AG'. To avoid this, we
            # trim off the last couple unique bases of the alignment, so the example becomes:
            #   read: ... T G A G T A C
            #   ref:  ... T G A G T A C A G G G G A A G T C C A G T ...

            last_base = aligned_bases[-1]
            while aligned_bases[-1] == last_base:
                aligned_bases.pop()
            aligned_bases.pop()

            for i, bases in enumerate(aligned_bases):
                if bases in pileup[a.ref_start + i]:
                    pileup[a.ref_start + i][bases] += 1
                else:
                    pileup[a.ref_start + i][bases] = 1
                depths_by_pos[a.ref_start + i] += depth_contribution
    return pileup, depths_by_pos


def select_best_alignments(alignments, mask_positions, read_pair_names, read_count, target_seqs):
    section_header('Selecting best alignments')
    explanation('Polypolish now chooses the best alignment(s) for each read, ignoring the masked '
                'positions of the target sequence. I.e. each read\'s alignments are ranked from '
                'fewest-errors to most-errors in unmasked positions of the reference, and only '
                'the best alignments are kept.')
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
