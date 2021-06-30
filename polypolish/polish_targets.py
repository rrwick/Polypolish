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

import statistics

from .log import log, section_header, explanation
from . import settings


def polish_target_sequences(alignments, assembly_seqs, debug, min_depth, min_fraction):
    section_header('Polishing assembly sequences')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. '
                'Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu '
                'fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in '
                'culpa qui officia deserunt mollit anim id est laborum.')
    new_lengths = []
    for target_name, target_seq in assembly_seqs:
        new_name, new_length = polish_target_sequence(target_name, target_seq, alignments, debug,
                                                      min_depth, min_fraction)
        new_lengths.append((new_name, new_length))
    return new_lengths


def polish_target_sequence(target_name, target_seq, alignments, debug, min_depth, min_fraction):
    log(f'Polishing {target_name} ({len(target_seq):,} bp):')

    changed_positions = set()
    log('  Building read pileup:')
    pileup, depths_by_pos = get_pileup(alignments, target_name, target_seq)
    mean_depth = statistics.mean(depths_by_pos.values())
    log(f'    mean read depth: {mean_depth:.1f}x')
    uncovered = sum(1 if d == 0.0 else 0 for d in depths_by_pos.values())
    coverage = 100.0 * (len(target_seq) - uncovered) / len(target_seq)
    have = 'has' if uncovered == 1 else 'have'
    log(f'    {uncovered:,} bp {have} a depth of zero ({coverage:.4f}% coverage)')

    log('  Repairing assembly sequence:')
    if debug:
        log()
        log('Column_1=ref_pos  Column_2=ref_base  Column_3=depth  Column_4=match_fraction  '
            'Column_5=read_pileup')

    new_bases = []
    for i in range(len(target_seq)):
        if i % settings.POLISH_PROGRESS_INTERVAL == 0 and not debug:
            log(f'\r    {i:,} / {len(target_seq):,} bases', end='')

        read_base_counts = pileup[i]
        depth = depths_by_pos[i]
        ref_base = target_seq[i]

        if depth > 0.0:
            target_count = max(min_depth, min_fraction * depth)
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

    if not debug:
        log(f'\r    {len(target_seq):,} / {len(target_seq):,} bases')
    changed_count = len(changed_positions)
    changed_percent = 100.0 * changed_count / len(target_seq)
    estimated_accuracy = 100.0 - changed_percent

    positions = 'position' if changed_count == 1 else 'positions'
    log(f'    {changed_count:,} {positions} changed ({changed_percent:.4f}% of total positions)')
    log(f'    estimated pre-polishing sequence accuracy: {estimated_accuracy:.4f}%')
    log()

    new_name = f'{target_name}_polypolish'
    new_sequence = ''.join(new_bases)
    new_sequence = new_sequence.replace('-', '')

    print(f'>{new_name}')
    print(new_sequence)
    return new_name, len(new_sequence)


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
            try:
                while aligned_bases[-1] == last_base:
                    aligned_bases.pop()
                aligned_bases.pop()
            except IndexError:  # popped all the bases off
                pass

            for i, bases in enumerate(aligned_bases):
                if bases in pileup[a.ref_start + i]:
                    pileup[a.ref_start + i][bases] += 1
                else:
                    pileup[a.ref_start + i][bases] = 1
                depths_by_pos[a.ref_start + i] += depth_contribution

    return pileup, depths_by_pos
