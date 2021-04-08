"""
This module contains some tests for Hyalign. To run them, execute `pytest` from the root Hyalign
directory.

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

import hyalign.mask_reads


def test_flip_positions_1():
    """
    ***
    0123456789
    9876543210
    """
    assert hyalign.mask_reads.flip_positions({0, 1, 2}, 10) == {7, 8, 9}


def test_flip_positions_2():
    """
     *   ** *
    0123456789
    9876543210
    """
    assert hyalign.mask_reads.flip_positions({1, 5, 6, 8}, 10) == {1, 3, 4, 8}


def test_flip_positions_3():
    """
     *   ** *
    012345678
    876543210
    """
    assert hyalign.mask_reads.flip_positions({1, 5, 6, 8}, 9) == {0, 2, 3, 7}
