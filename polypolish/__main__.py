#!/usr/bin/env python3
"""
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

import argparse
import random
import sys

from .alignment import align_reads
from .help_formatter import MyParser, MyHelpFormatter
from .insert_size import get_insert_size_distribution, select_alignments_using_insert_size
from .log import log, bold, section_header, explanation
from .polish_targets import polish_target_sequences
from .misc import get_default_thread_count, get_ascii_art, load_fasta
from .software import check_python, check_minimap2
from .version import __version__


def main():
    args = parse_args()
    starting_message(args)
    random.seed(args.seed)
    assembly_seqs = load_fasta(args.assembly)

    alignments, read_pair_names, read_count = \
        align_reads(args.assembly, args.short1, args.short2, args.threads, args.max_errors,
                    args.debug)

    insert_size_distribution = get_insert_size_distribution(alignments)
    select_alignments_using_insert_size(alignments, insert_size_distribution, read_pair_names,
                                        read_count)

    polish_target_sequences(alignments, assembly_seqs, args.debug)


def parse_args():
    description = 'R|' + get_ascii_art() + '\n' + \
                  bold('Polypolish: a tool for polishing a long-read assembly - repeats included')
    parser = MyParser(description=description, formatter_class=MyHelpFormatter, add_help=False)

    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument('-a', '--assembly', type=str, required=True,
                               help='Assembly to polish (FASTA format)')
    required_args.add_argument('-1', '--short1', type=str, required=True,
                               help='Input short reads, first in pair (FASTQ format)')
    required_args.add_argument('-2', '--short2', type=str, required=True,
                               help='Input short reads, second in pair (FASTQ format)')

    setting_args = parser.add_argument_group('Settings')
    setting_args.add_argument('-m', '--max_errors', type=int, default=10,
                              help='Ignore alignments with more than this number of mismatches '
                                   'and indels')
    setting_args.add_argument('-s', '--seed', type=int, default=0,
                              help='Seed for random number generator')
    setting_args.add_argument('-t', '--threads', type=int, default=get_default_thread_count(),
                              help='Number of threads')

    setting_args.add_argument('--debug', action='store_true',
                              help='Output lots of extra information (for debugging purposes)')

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version', version='Polypolish v' + __version__,
                           help="Show program's version number and exit")

    # If no arguments were used, print the base-level help which lists possible commands.
    if len(sys.argv) == 1:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    return parser.parse_args()


def starting_message(args):
    section_header('Starting Polypolish')
    explanation('Polypolish is a paired-end read aligner which aims to produce high-quality '
                'alignments for polishing. Specifically, it takes as input paired-end short reads '
                'and a haploid genome assembly to be polished (e.g. a long-read assembly). '
                'Instead of simply aligning each read to its best location (as most aligners do), '
                'it aims to align reads to where they will be most useful for polishing.')
    log(f'Polypolish version: v{__version__}')
    log(f'Using {args.threads} threads')
    log()
    log('Input short reads:')
    log(f'  {args.short1}')
    log(f'  {args.short2}')
    log()
    log('Input assembly:')
    log(f'  {args.assembly}')
    log()
    check_requirements()
    log()


def check_requirements():
    log('Checking requirements:')
    check_python()
    check_minimap2()


if __name__ == '__main__':
    main()