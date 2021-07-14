"""
This module contains some tests for Polypolish. To run them, execute `pytest` from the root
Polypolish directory.

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

import polypolish.alignment


def test_get_expanded_cigar():
    assert polypolish.alignment.get_expanded_cigar('5M2I3M1D4M') == 'MMMMMIIMMMDMMMM'


def test_get_ref_end():
    assert polypolish.alignment.get_ref_end(50, '10M') == 60
    assert polypolish.alignment.get_ref_end(100, '10M2I10M') == 120
    assert polypolish.alignment.get_ref_end(100, '10M2D10M') == 122


def test_add_secondary_read_seqs():
    # Both reads on + strand.
    a_1 = polypolish.alignment.Alignment('read\t0\tref\t10\t\t10M\t\t\t\tACGTACGTCG\t0123456789')
    a_2 = polypolish.alignment.Alignment('read\t256\tref\t10\t\t10M\t\t\t\t*\t*')
    alignments = {'read': [a_1, a_2]}
    polypolish.alignment.add_secondary_read_seqs(alignments)
    assert a_2.read_seq == 'ACGTACGTCG'

    # Both reads on - strand.
    a_1 = polypolish.alignment.Alignment('read\t16\tref\t10\t\t10M\t\t\t\tACGTACGTCG\t0123456789')
    a_2 = polypolish.alignment.Alignment('read\t272\tref\t10\t\t10M\t\t\t\t*\t*')
    alignments = {'read': [a_1, a_2]}
    polypolish.alignment.add_secondary_read_seqs(alignments)
    assert a_2.read_seq == 'ACGTACGTCG'

    # Primary read on + strand, secondary read on - strand.
    a_1 = polypolish.alignment.Alignment('read\t0\tref\t10\t\t10M\t\t\t\tACGTACGTCG\t0123456789')
    a_2 = polypolish.alignment.Alignment('read\t272\tref\t10\t\t10M\t\t\t\t*\t*')
    alignments = {'read': [a_1, a_2]}
    polypolish.alignment.add_secondary_read_seqs(alignments)
    assert a_2.read_seq == 'CGACGTACGT'

    # Primary read on - strand, secondary read on + strand.
    a_1 = polypolish.alignment.Alignment('read\t16\tref\t10\t\t10M\t\t\t\tACGTACGTCG\t0123456789')
    a_2 = polypolish.alignment.Alignment('read\t256\tref\t10\t\t10M\t\t\t\t*\t*')
    alignments = {'read': [a_1, a_2]}
    polypolish.alignment.add_secondary_read_seqs(alignments)
    assert a_2.read_seq == 'CGACGTACGT'


def test_starts_and_ends_with_match():
    a = polypolish.alignment.Alignment('read\t0\tref\t10\t\t10M\t\t\t\tACGTACGTCG\t0123456789')
    assert a.starts_and_ends_with_match()

    a = polypolish.alignment.Alignment('read\t0\tref\t10\t\t2S8M\t\t\t\tACGTACGTCG\t0123456789')
    assert not a.starts_and_ends_with_match()

    a = polypolish.alignment.Alignment('read\t0\tref\t10\t\t8M2S\t\t\t\tACGTACGTCG\t0123456789')
    assert not a.starts_and_ends_with_match()
