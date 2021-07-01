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
import datetime
import sys

from .alignment import load_alignments
from .help_formatter import MyParser, MyHelpFormatter
from .log import log, bold, section_header, explanation
from .polish_targets import polish_target_sequences
from .misc import get_ascii_art, load_fasta
from .version import __version__


def main():
    args = parse_args()
    start_time = starting_message(args)
    assembly_seqs = load_assembly(args.assembly)
    alignments = load_alignments(args.sam1, args.sam2, args.max_errors)
    new_lengths = polish_target_sequences(alignments, assembly_seqs, args.debug, args.min_depth,
                                          args.min_fraction)
    finished_message(start_time, new_lengths)


def parse_args():
    description = 'R|' + get_ascii_art() + '\n' + \
                  bold('Polypolish: a tool for polishing a long-read assembly - repeats included')
    parser = MyParser(description=description, formatter_class=MyHelpFormatter, add_help=False)

    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument('-a', '--assembly', type=str, required=True,
                               help='Assembly to polish (FASTA format)')
    required_args.add_argument('-1', '--sam1', type=str, required=True,
                               help='Input short read alignments, unpaired or first in pair '
                                    '(SAM format)')
    required_args.add_argument('-2', '--sam2', type=str, required=False,
                               help='Input short read alignments, second in pair (SAM format)')

    setting_args = parser.add_argument_group('Settings')
    setting_args.add_argument('-m', '--max_errors', type=int, default=10,
                              help='Ignore alignments with more than this number of mismatches '
                                   'and indels')
    setting_args.add_argument('-d', '--min_depth', type=int, default=5,
                              help='A base must occur at least this many times in the pileup to '
                                   'be considered valid')
    setting_args.add_argument('-f', '--min_fraction', type=float, default=0.5,
                              help='A base must make up at least this fraction of the pileup to '
                                   'be considered valid')
    setting_args.add_argument('--debug', type=str,
                              help='Optional file in which to store per-base information for '
                                   'debugging purposes')

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
    explanation('Polypolish is a tool for polishing genome assemblies with short reads. Unlike '
                'other tools in this category, Polypolish uses SAM files where each read has been '
                'aligned to all possible locations (not just a single best location). This allows '
                'it to repair errors in repeat regions that other alignment-based polishers '
                'cannot fix.')

    log(f'Polypolish version: v{__version__}')
    log()
    log('Input assembly:')
    log(f'  {args.assembly}')
    log()
    log('Input short-read alignments:')
    log(f'  {args.sam1}')
    if args.sam2 is not None:
        log(f'  {args.sam2}')
    log()
    log('Settings:')
    log(f'  --max_errors {args.max_errors}')
    log(f'  --min_depth {args.min_depth}')
    log(f'  --min_fraction {args.min_fraction}')
    if args.debug is None:
        log(f'  not logging debugging information')
    else:
        log(f'  --debug {args.debug}')
    log()
    check_python()
    return datetime.datetime.now()


def finished_message(start_time, new_lengths, debug):
    section_header('Finished!')
    log('Polished sequence (to stdout):')
    for new_name, new_length in new_lengths:
        log(f'  {new_name} ({new_length:,} bp)')
    log()
    if debug is not None:
        log(f'Per-base debugging info written to {debug}')
    time_to_run = datetime.datetime.now() - start_time
    log(f'Time to run: {time_to_run}')
    log()


def load_assembly(assembly_filename):
    section_header('Loading assembly')
    assembly_seqs = load_fasta(assembly_filename)
    count = len(assembly_seqs)
    if count == 0:
        sys.exit(f'Error: no sequences in {assembly_filename}')
    log(f'{count} sequence{"" if count == 1 else "s"} in {assembly_filename}:')
    for name, seq in assembly_seqs:
        log(f'  {name} ({len(seq):,} bp)')
    log()
    return assembly_seqs


def check_python():
    major, minor = sys.version_info.major, sys.version_info.minor
    good_version = (major >= 3 and minor >= 6)
    if not good_version:
        sys.exit('Error: Polypolish requires Python 3.6 or later')


if __name__ == '__main__':
    main()
