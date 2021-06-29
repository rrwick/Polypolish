"""
This module contains some tests for Repeatish. To run them, execute `pytest` from the root Repeatish
directory.

Copyright 2021 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Repeatish

This file is part of Repeatish. Repeatish is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Repeatish is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Repeatish.
If not, see <http://www.gnu.org/licenses/>.
"""

import repeatish.mask_targets
import test.test_alignment


def test_get_pileup():
    a = test.test_alignment.prep_alignment()
    ref_seq = 'ATCGCAATTGTAGAAGGACCTAGGAAGCAAAAGTTTC' \
              'CTCTATGACGAGACGAACGTCGCTAACGAGCGACCTATAGCGTTTAAAATA' \
              'GGATTATCGAACACCGGTAG'
    alignments = {'read': [a]}
    pileup, _ = repeatish.mask_targets.get_pileup(alignments, 'ref', ref_seq)
    for i in range(37):
        assert pileup[i] == {}
    assert pileup[37] == {'C': 1}
    assert pileup[38] == {'T': 1}
    assert pileup[39] == {'C': 1}
    assert pileup[40] == {'T': 1}
    assert pileup[41] == {'A': 1}
    assert pileup[42] == {'T': 1}
    assert pileup[43] == {'G': 1}
    assert pileup[44] == {'A': 1}
    assert pileup[45] == {'C': 1}
    assert pileup[46] == {'G': 1}
    assert pileup[47] == {'AC': 1}
    assert pileup[48] == {'G': 1}
    assert pileup[49] == {'A': 1}
    assert pileup[50] == {'-': 1}
    assert pileup[51] == {'-': 1}
    assert pileup[52] == {'A': 1}
    assert pileup[53] == {'A': 1}
    assert pileup[54] == {'C': 1}
    assert pileup[55] == {'G': 1}
    assert pileup[56] == {'T': 1}
    assert pileup[57] == {'C': 1}
    assert pileup[58] == {'G': 1}
    assert pileup[59] == {'C': 1}
    assert pileup[60] == {'T': 1}
    assert pileup[61] == {'CTGT': 1}
    assert pileup[62] == {'A': 1}
    assert pileup[63] == {'C': 1}
    assert pileup[64] == {'G': 1}
    assert pileup[65] == {'A': 1}
    assert pileup[66] == {'G': 1}
    assert pileup[67] == {'C': 1}
    assert pileup[68] == {'G': 1}
    assert pileup[69] == {'A': 1}
    assert pileup[70] == {'C': 1}
    assert pileup[71] == {'-': 1}
    assert pileup[72] == {'T': 1}
    assert pileup[73] == {'A': 1}
    assert pileup[74] == {'T': 1}
    assert pileup[75] == {'A': 1}
    assert pileup[76] == {'G': 1}
    assert pileup[77] == {'C': 1}
    assert pileup[78] == {'G': 1}
    assert pileup[79] == {'T': 1}
    assert pileup[80] == {'T': 1}
    assert pileup[81] == {'-': 1}
    assert pileup[82] == {'A': 1}
    assert pileup[83] == {'A': 1}
    assert pileup[84] == {'A': 1}
    assert pileup[85] == {'A': 1}
    # Last couple bases are trimmed off.
    for i in range(86, 108):
        assert pileup[i] == {}
