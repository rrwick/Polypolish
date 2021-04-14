"""
This module contains some tests for Poligner. To run them, execute `pytest` from the root Poligner
directory.

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

import poligner.alignment


def test_flags_1():
    flags = 345
    a = poligner.alignment.Alignment(f'read\t{flags}\tref\t0\t0\t4M\t*\t0\t0\tACGT\tAAAA\tNM:i:0')
    assert a.has_flag(1)
    assert not a.has_flag(2)
    assert not a.has_flag(4)
    assert a.has_flag(8)
    assert a.has_flag(16)
    assert not a.has_flag(32)
    assert a.has_flag(64)
    assert not a.has_flag(128)
    assert a.has_flag(256)
    assert not a.has_flag(512)
    assert not a.has_flag(1024)
    assert not a.has_flag(2048)


def test_flags_2():
    flags = 1044
    a = poligner.alignment.Alignment(f'read\t{flags}\tref\t0\t0\t4M\t*\t0\t0\tACGT\tAAAA\tNM:i:0')
    assert not a.has_flag(1)
    assert not a.has_flag(2)
    assert a.has_flag(4)
    assert not a.has_flag(8)
    assert a.has_flag(16)
    assert not a.has_flag(32)
    assert not a.has_flag(64)
    assert not a.has_flag(128)
    assert not a.has_flag(256)
    assert not a.has_flag(512)
    assert a.has_flag(1024)
    assert not a.has_flag(2048)


def test_ref_end_1():
    assert poligner.alignment.get_ref_end(100, '100M') == 200


def test_ref_end_2():
    assert poligner.alignment.get_ref_end(1000, '10M3D8M4I10M') == 1031


def prep_alignment():
    """
    cigar: MMMMMMMMMMMIMMDDMMMMMMMMMMIIIMMMMMMMMMDMMMMMMMMMDMMMMMM

    read:  CTCTATGACGACGA--AACGTCGCTCTGTACGAGCGAC-TATAGCGTT-AAAATA
    ref:   CTCTATGACGA-GACGAACGTCGCTA---ACGAGCGACCTATAGCGTTTAAAATA
    diff:             *  **         ****         *         *

    read   00000000001111  1111112222222222333333 333344444 444455
    pos:   01234567890123  4567890123456789012345 678901234 567890

    ref    33344444444 44555555555566   66666666777777777788888888
    pos:   78901234567 89012345678901   23456789012345678901234567
    """
    ref_start = 37
    cigar = '11M1I2M2D10M3I9M1D9M1D6M'
    read_seq = 'CTCTATGACGACGAAACGTCGCTCTGTACGAGCGACTATAGCGTTAAAATA'
    ref_seq = 'CTCTATGACGAGACGAACGTCGCTAACGAGCGACCTATAGCGTTTAAAATA'
    read_qual = 'A' * len(read_seq)

    a = poligner.alignment.Alignment(f'read\t0\tref\t{ref_start+1}\t0\t{cigar}\t*\t0\t0\t'
                                     f'{read_seq}\t{read_qual}\tNM:i:0')
    a.add_detailed_alignment_info(ref_seq)
    return a


def test_prep_for_detailed_alignment_info_1():
    a = prep_alignment()
    assert a.aligned_read_seq == 'CTCTATGACGACGA--AACGTCGCTCTGTACGAGCGAC-TATAGCGTT-AAAATA'
    assert a.aligned_ref_seq == 'CTCTATGACGA-GACGAACGTCGCTA---ACGAGCGACCTATAGCGTTTAAAATA'
    assert a.diffs == '           *  **         ****         *         *      '


def test_prep_for_detailed_alignment_info_2():
    a = prep_alignment()
    assert a.read_error_positions == [11, 13, 13, 23, 24, 25, 26, 35, 44]
    assert a.ref_error_positions == [47, 50, 51, 61, 61, 61, 61, 71, 81]


def test_prep_for_detailed_alignment_info_3():
    a = prep_alignment()
    assert a.read_positions_to_ref_positions[5] == [42]
    assert a.read_positions_to_ref_positions[19] == [57]
    assert a.read_positions_to_ref_positions[31] == [66]
    assert a.read_positions_to_ref_positions[48] == [85]

    assert a.ref_positions_to_read_positions[42] == [5]
    assert a.ref_positions_to_read_positions[57] == [19]
    assert a.ref_positions_to_read_positions[66] == [31]
    assert a.ref_positions_to_read_positions[85] == [48]


def test_prep_for_detailed_alignment_info_4():
    a = prep_alignment()
    assert a.read_positions_to_ref_positions[10] == [47]
    assert a.read_positions_to_ref_positions[11] == [47]
    assert a.read_positions_to_ref_positions[12] == [48]

    assert a.ref_positions_to_read_positions[47] == [10, 11]
    assert a.ref_positions_to_read_positions[48] == [12]

    assert a.read_positions_to_ref_positions[22] == [60]
    assert a.read_positions_to_ref_positions[23] == [61]
    assert a.read_positions_to_ref_positions[24] == [61]
    assert a.read_positions_to_ref_positions[25] == [61]
    assert a.read_positions_to_ref_positions[26] == [61]
    assert a.read_positions_to_ref_positions[27] == [62]

    assert a.ref_positions_to_read_positions[60] == [22]
    assert a.ref_positions_to_read_positions[61] == [23, 24, 25, 26]
    assert a.ref_positions_to_read_positions[62] == [27]


def test_prep_for_detailed_alignment_info_5():
    a = prep_alignment()
    assert a.ref_positions_to_read_positions[49] == [13]
    assert a.ref_positions_to_read_positions[50] == [13]
    assert a.ref_positions_to_read_positions[51] == [13]
    assert a.ref_positions_to_read_positions[52] == [14]

    assert a.read_positions_to_ref_positions[13] == [49, 50, 51]
    assert a.read_positions_to_ref_positions[14] == [52]

    assert a.ref_positions_to_read_positions[70] == [35]
    assert a.ref_positions_to_read_positions[71] == [35]
    assert a.ref_positions_to_read_positions[72] == [36]

    assert a.read_positions_to_ref_positions[35] == [70, 71]
    assert a.read_positions_to_ref_positions[36] == [72]

    assert a.ref_positions_to_read_positions[80] == [44]
    assert a.ref_positions_to_read_positions[81] == [44]
    assert a.ref_positions_to_read_positions[82] == [45]

    assert a.read_positions_to_ref_positions[44] == [80, 81]
    assert a.read_positions_to_ref_positions[45] == [82]


def test_get_read_bases_for_each_target_base():
    a = prep_alignment()
    ref_seq = 'CTCTATGACGAGACGAACGTCGCTAACGAGCGACCTATAGCGTTTAAAATA'
    read_bases = a.get_read_bases_for_each_target_base(ref_seq)
    assert read_bases == ['C', 'T', 'C', 'T', 'A', 'T', 'G', 'A', 'C', 'G', 'AC', 'G', 'A', '-',
                          '-', 'A', 'A', 'C', 'G', 'T', 'C', 'G', 'C', 'T', 'CTGT', 'A', 'C', 'G',
                          'A', 'G', 'C', 'G', 'A', 'C', '-', 'T', 'A', 'T', 'A', 'G', 'C', 'G',
                          'T', 'T', '-', 'A', 'A', 'A', 'A', 'T', 'A']


def test_flip_positions_1():
    """
    ***
    0123456789
    9876543210
    """
    assert poligner.alignment.flip_positions([0, 1, 2], 10) == [7, 8, 9]


def test_flip_positions_2():
    """
     *   ** *
    0123456789
    9876543210
    """
    assert poligner.alignment.flip_positions([1, 5, 6, 8], 10) == [1, 3, 4, 8]


def test_flip_positions_3():
    """
     *   ** *
    012345678
    876543210
    """
    assert poligner.alignment.flip_positions([1, 5, 6, 8], 9) == [0, 2, 3, 7]
