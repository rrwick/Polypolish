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

import gzip
import pytest
import sys
import unittest.mock

import polypolish.misc


def test_get_compression_type_1():
    assert polypolish.misc.get_compression_type('test/test_misc/test.txt') == 'plain'


def test_get_compression_type_2():
    assert polypolish.misc.get_compression_type('test/test_misc/test.gz') == 'gz'


def test_get_compression_type_3():
    with pytest.raises(SystemExit) as e:
        polypolish.misc.get_compression_type('test/test_misc/test.bz2')
    assert 'cannot use bzip2' in str(e.value)


def test_get_compression_type_4():
    with pytest.raises(SystemExit) as e:
        polypolish.misc.get_compression_type('test/test_misc/test.zip')
    assert 'cannot use zip' in str(e.value)


def test_get_open_func_1():
    assert polypolish.misc.get_open_func('test/test_misc/test.txt') == open


def test_get_open_func_2():
    assert polypolish.misc.get_open_func('test/test_misc/test.gz') == gzip.open


def test_load_fasta_1():
    seqs = polypolish.misc.load_fasta('test/test_misc/test.fasta')
    assert len(seqs) == 2
    assert seqs[0][0] == 'A'
    assert seqs[0][1].startswith('TTGCCTGTAGTCGGGACC')
    assert seqs[1][0] == 'B'
    assert seqs[1][1].startswith('ATTCTCAGAATGGCGTAG')


def test_load_fasta_2():
    seqs = polypolish.misc.load_fasta('test/test_misc/test.fasta', include_full_header=True)
    assert len(seqs) == 2
    assert seqs[0][0] == 'A'
    assert seqs[0][1] == 'A info'
    assert seqs[0][2].startswith('TTGCCTGTAGTCGGGACC')
    assert seqs[1][0] == 'B'
    assert seqs[1][1] == 'B stuff'
    assert seqs[1][2].startswith('ATTCTCAGAATGGCGTAG')


def test_load_fasta_3():
    seqs = polypolish.misc.load_fasta('test/test_misc/test.fasta.gz')
    assert len(seqs) == 2
    assert seqs[0][0] == 'A'
    assert seqs[0][1].startswith('TTGCCTGTAGTCGGGACC')
    assert seqs[1][0] == 'B'
    assert seqs[1][1].startswith('ATTCTCAGAATGGCGTAG')


def test_load_fasta_4():
    seqs = polypolish.misc.load_fasta('test/test_misc/bad_1.fasta')
    assert len(seqs) == 2
    assert seqs[0][0] == 'A'
    assert seqs[0][1].startswith('TTGCCTGTAGTCGGGACC')
    assert seqs[1][0] == 'B'
    assert seqs[1][1].startswith('ATTCTCAGAATGGCGTAG')


def test_load_fasta_5():
    seqs = polypolish.misc.load_fasta('test/test_misc/lowercase.fasta')
    assert len(seqs) == 2
    assert seqs[0][0] == 'A'
    assert seqs[0][1].startswith('TTGCCTGTAGTCGGGACC')
    assert seqs[1][0] == 'B'
    assert seqs[1][1].startswith('ATTCTCAGAATGGCGTAG')


def test_reverse_complement_1():
    assert polypolish.misc.reverse_complement('GGGGaaaaaaaatttatatat') == 'atatataaattttttttCCCC'


def test_reverse_complement_2():
    assert polypolish.misc.reverse_complement('atatataaattttttttCCCC') == 'GGGGaaaaaaaatttatatat'


def test_reverse_complement_3():
    assert polypolish.misc.reverse_complement('ACGT123') == 'NNNACGT'


def test_check_python_version_1():
    with unittest.mock.patch.object(sys, 'version_info') as v_info:
        v_info.major = 3
        v_info.minor = 6
        polypolish.misc.check_python_version()


def test_check_python_version_2():
    with unittest.mock.patch.object(sys, 'version_info') as v_info:
        v_info.major = 3
        v_info.minor = 8
        polypolish.misc.check_python_version()


def test_check_python_version_3():
    with pytest.raises(SystemExit) as e:
        with unittest.mock.patch.object(sys, 'version_info') as v_info:
            v_info.major = 3
            v_info.minor = 5
            polypolish.misc.check_python_version()
    assert 'requires Python 3.6 or later' in str(e.value)


def test_check_python_version_4():
    with pytest.raises(SystemExit) as e:
        with unittest.mock.patch.object(sys, 'version_info') as v_info:
            v_info.major = 2
            v_info.minor = 7
            polypolish.misc.check_python_version()
    assert 'requires Python 3.6 or later' in str(e.value)


def test_get_ascii_art():
    assert "| |                     | |(_)     | |" in polypolish.misc.get_ascii_art()
