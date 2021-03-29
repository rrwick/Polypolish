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

import hyalign.alignment


def test_flags_1():
    flags = 345
    a = hyalign.alignment.Alignment(f'read\t{flags}\tref\t0\t0\t4M\t*\t0\t0\tACGT\tAAAA')
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
    a = hyalign.alignment.Alignment(f'read\t{flags}\tref\t0\t0\t4M\t*\t0\t0\tACGT\tAAAA')
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
    assert hyalign.alignment.get_ref_end(100, '100M') == 200


def test_ref_end_2():
    assert hyalign.alignment.get_ref_end(1000, '10M3D8M4I10M') == 1031
