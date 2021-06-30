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
import sys

from .alignment import load_alignments
from .help_formatter import MyParser, MyHelpFormatter
from .log import log, bold, section_header, explanation
from .polish_targets import polish_target_sequences
from .misc import get_default_thread_count, get_ascii_art, load_fasta
from .software import check_python, check_minimap2
from .version import __version__


def main():
    args = parse_args()
    starting_message(args)
    assembly_seqs = load_fasta(args.assembly)
    alignments = load_alignments(args.sam1, args.sam2, args.max_errors)
    polish_target_sequences(alignments, assembly_seqs, args.debug)


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
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. '
                'Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu '
                'fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in '
                'culpa qui officia deserunt mollit anim id est laborum.')

    log(f'Polypolish version: v{__version__}')
    log(f'Using {args.threads} threads')
    log()
    log('Input short-read alignments:')
    log(f'  {args.sam1}')
    if args.sam2 is not None:
        log(f'  {args.sam2}')
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
