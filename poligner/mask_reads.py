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

from .alignment import get_multi_alignment_read_names, flip_positions
from .log import log, section_header, explanation
from .misc import load_fasta


def mask_read_sequences(read_pair_names, alignments, target):
    section_header('Masking read sequences')
    explanation('Before continuing, Poligner masks out read bases which do not match the reference '
                'at any location. I.e. if a read base appears to be a mistake at each of the '
                'read\'s alignment locations, it is assumed to be a read error and ignored from '
                'this point onward.')

    target_seqs = dict(load_fasta(target))
    multi_alignment_read_names = get_multi_alignment_read_names(read_pair_names, alignments)

    for name in multi_alignment_read_names:
        for a in alignments[name]:
            ref_seq = target_seqs[a.ref_name][a.ref_start:a.ref_end]
            a.add_detailed_alignment_info(ref_seq)

    total_masked_bases = 0
    for name in multi_alignment_read_names:
        read_alignments = alignments[name]
        mask_positions = get_mask_positions(read_alignments)
        if mask_positions:
            for a in read_alignments:
                a.masked_read_positions = mask_positions
        total_masked_bases += len(mask_positions)

    log(f'Masked {total_masked_bases:,} bases in {len(multi_alignment_read_names):,} reads')
    log()


def get_mask_positions(read_alignments):
    """
    Returns a list of positions in the read to mask. Positions are given in terms of the forward
    strand.
    """
    mask_positions = None
    for a in read_alignments:
        if a.is_on_reverse_strand():
            read_errors = flip_positions(a.read_error_positions, len(a.read_seq))
        else:
            read_errors = a.read_error_positions

        if mask_positions is None:
            mask_positions = set(read_errors)
        else:
            mask_positions &= set(read_errors)
    return mask_positions
