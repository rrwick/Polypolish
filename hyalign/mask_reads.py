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

from .alignment import get_multi_alignment_read_names, flip_positions
from .log import log, section_header, explanation
from .misc import load_fasta


def mask_read_sequences(read_pair_names, alignments, target):
    section_header('Masking read sequences')
    explanation('Before continuing, Hyalign masks out read bases which do not match the reference '
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
                create_masked_read_seq(a, mask_positions)
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
            mask_positions = read_errors
        else:
            mask_positions &= read_errors
    return mask_positions


def create_masked_read_seq(alignment, mask_positions):
    """
    This function creates a masked version of the read sequence in the alignment, where all mask
    positions are replaced with 'N' bases.
    """
    alignment.masked_read_seq = alignment.read_seq[:]
    if alignment.is_on_reverse_strand():
        mask_positions = flip_positions(mask_positions, len(alignment.read_seq))
    for i in mask_positions:
        alignment.masked_read_seq = \
            alignment.masked_read_seq[:i] + 'N' + alignment.masked_read_seq[i + 1:]