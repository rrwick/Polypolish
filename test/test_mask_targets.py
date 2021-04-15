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

import poligner.mask_targets


def test_get_possible_mutations_1():
    assert poligner.mask_targets.get_possible_mutations('A') == \
        sorted(['C', 'G', 'T', '', 'AA', 'AC', 'AG', 'AT'])
    assert poligner.mask_targets.get_possible_mutations('C') == \
        sorted(['A', 'G', 'T', '', 'CA', 'CC', 'CG', 'CT'])
    assert poligner.mask_targets.get_possible_mutations('G') == \
        sorted(['A', 'C', 'T', '', 'GA', 'GC', 'GG', 'GT'])
    assert poligner.mask_targets.get_possible_mutations('T') == \
        sorted(['A', 'C', 'G', '', 'TA', 'TC', 'TG', 'TT'])


def test_get_possible_mutations_2():
    assert poligner.mask_targets.get_possible_mutations('AA') == \
        sorted(['AC', 'AG', 'AT', 'A', 'AAA', 'AAC', 'AAG', 'AAT'])
    assert poligner.mask_targets.get_possible_mutations('AG') == \
        sorted(['AA', 'AC', 'AT', 'A', 'AGA', 'AGC', 'AGG', 'AGT'])
