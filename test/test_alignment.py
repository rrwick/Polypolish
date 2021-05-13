"""
This module contains some tests for Maskimap. To run them, execute `pytest` from the root Maskimap
directory.

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

import maskimap.alignment


def test_flags_1():
    flags = 345
    a = maskimap.alignment.Alignment(f'read\t{flags}\tref\t0\t0\t4M\t*\t0\t0\tACGT\tAAAA\tNM:i:0')
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
    a = maskimap.alignment.Alignment(f'read\t{flags}\tref\t0\t0\t4M\t*\t0\t0\tACGT\tAAAA\tNM:i:0')
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
    assert maskimap.alignment.get_ref_end(100, '100M') == 200


def test_ref_end_2():
    assert maskimap.alignment.get_ref_end(1000, '10M3D8M4I10M') == 1031


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
    read_qual = 'A' * len(read_seq)

    a = maskimap.alignment.Alignment(f'read\t0\tref\t{ref_start+1}\t0\t{cigar}\t*\t0\t0\t'
                                     f'{read_seq}\t{read_qual}\tNM:i:0')
    return a


def test_get_read_bases_for_each_target_base():
    a = prep_alignment()
    ref_seq = 'CTCTATGACGAGACGAACGTCGCTAACGAGCGACCTATAGCGTTTAAAATA'
    read_bases = a.get_read_bases_for_each_target_base(ref_seq)
    assert read_bases == ['C', 'T', 'C', 'T', 'A', 'T', 'G', 'A', 'C', 'G', 'AC', 'G', 'A', '-',
                          '-', 'A', 'A', 'C', 'G', 'T', 'C', 'G', 'C', 'T', 'CTGT', 'A', 'C', 'G',
                          'A', 'G', 'C', 'G', 'A', 'C', '-', 'T', 'A', 'T', 'A', 'G', 'C', 'G',
                          'T', 'T', '-', 'A', 'A', 'A', 'A', 'T', 'A']


def test_get_mapq():
    assert maskimap.alignment.get_mapq(0) == 0
    assert maskimap.alignment.get_mapq(1) == 60
    assert maskimap.alignment.get_mapq(2) == 3
    assert maskimap.alignment.get_mapq(3) == 2
    assert maskimap.alignment.get_mapq(4) == 1
    assert maskimap.alignment.get_mapq(5) == 1
    assert maskimap.alignment.get_mapq(10) == 1
    assert maskimap.alignment.get_mapq(20) == 1
    assert maskimap.alignment.get_mapq(1000) == 1
