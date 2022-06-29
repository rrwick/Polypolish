#!/usr/bin/env python3
"""
This script produces a human-readable file showing the changes that Polypolish has made to an
assembly. It can assist in spotting troublesome regions of the genome.

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
import edlib
import gzip
import mappy
import os
import re
import pathlib
import shutil
import subprocess
import sys
import textwrap

__version__ = '0.5.0'


def main():
    args = parse_args()
    check_inputs(args)
    start_time = starting_message(args)
    before, after = load_assemblies(args.before, args.after)
    align_sequences(before, after, args.padding, args.merge, args.aligner)
    finished_message(start_time)


def parse_args():
    description = 'R|' + get_ascii_art() + '\n' + \
                  'Polypolish human readable v' + __version__ + '\n' + \
                  'github.com/rrwick/Polypolish'
    parser = MyParser(description=description, formatter_class=MyHelpFormatter, add_help=False)

    input_args = parser.add_argument_group('Inputs')
    input_args.add_argument('before', type=str,
                            help='Assembly before polishing')
    input_args.add_argument('after', type=str,
                            help='Assembly after polishing')

    setting_args = parser.add_argument_group('Settings')
    setting_args.add_argument('--padding', type=int, default=15,
                              help='Bases of additional sequence to show before/after each change')
    setting_args.add_argument('--merge', type=int, default=30,
                              help='Changes this close are merged together in the output')
    setting_args.add_argument('--aligner', type=str, choices=['mappy', 'edlib'], default='mappy',
                              help='Aligner library: mappy has affine-gap, edlib is more robust')

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version',
                           version='Polypolish human readable v' + __version__,
                           help="Show program's version number and exit")

    # If no arguments were used, print the base-level help which lists possible commands.
    if len(sys.argv) == 1:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    return parser.parse_args()


def check_inputs(args):
    check_python_version()
    if not pathlib.Path(args.before).is_file():
        quit_with_error(f'Error: {args.before} is not a file')
    if not pathlib.Path(args.after).is_file():
        quit_with_error(f'Error: {args.after} is not a file')
    if args.padding < 0 or args.padding > 1000:
        quit_with_error('Error: the value of --padding must be >= 0 and <= 1000')


def starting_message(args):
    section_header('Polypolish human readable')
    explanation('This script produces a human-readable file showing all of the changes made to '
                'an assembly by Polypolish. This can help to spot troublesome regions of the '
                'assembly with a large number of changes.')
    log('Pre-polishing assembly:')
    log(f'  {args.before}')
    log()
    log('Post-polishing assembly:')
    log(f'  {args.after}')
    log()
    log('Settings:')
    log(f'  --padding {args.padding}')
    log(f'  --merge {args.merge}')
    log(f'  --aligner {args.aligner}')
    log()
    return datetime.datetime.now()


def load_assemblies(before_filename, after_filename):
    section_header('Loading assemblies')
    log(before_filename)
    before = load_fasta(before_filename)
    for name, seq in before:
        log(f'  {name}: {len(seq):,} bp')
    log()
    log(after_filename)
    after = load_fasta(after_filename)
    for name, seq in after:
        log(f'  {name}: {len(seq):,} bp')
    log()
    if len(before) != len(after):
        quit_with_error('Error: the before and after assemblies need to contain the same number '
                        'of sequences')
    for b, a in zip(before, after):
        before_name, before_seq = b
        before_seq_len = len(before_seq)
        after_name, after_seq = a
        after_seq_len = len(after_seq)
        if before_seq_len == 0 or after_seq_len == 0:
            quit_with_error('Error: zero-length sequences are not allowed')
        ratio = before_seq_len / after_seq_len
        if ratio < 0.9 or ratio > 1.11111111:
            quit_with_error(f'Error: {before_name} and {after_name} are too different in length '
                            f'- are the files in the same order?')
    return before, after


def align_sequences(before, after, padding, merge, aligner):
    section_header('Aligning sequences')
    longest_label = get_longest_label(before, after)
    for b, a in zip(before, after):
        before_name, before_seq = b
        after_name, after_seq = a
        output_differences(before_name, before_seq, after_name, after_seq, padding, merge,
                           longest_label, aligner)


def get_longest_label(before, after):
    longest_name, longest_seq = 0, 0
    for name, seq in before + after:
        longest_name = max(longest_name, len(name))
        longest_seq = max(longest_seq, len(str(len(seq))))
    return longest_name + 2*longest_seq + 3


def output_differences(before_name, before_seq, after_name, after_seq, padding, merge,
                       longest_label, aligner):
    log(f'Aligning {before_name} to {after_name}:')
    before_aligned, after_aligned, differences, before_pos, after_pos, diff_pos = \
        get_aligned_seqs(before_seq, after_seq, aligner)
    log(f'  {len(diff_pos):,} differences')

    aligned_len = len(before_aligned)
    diff_ranges = make_diff_ranges(diff_pos, padding, merge, aligned_len)

    for start, end in diff_ranges:
        # Convert positions in alignment to positions in unaligned sequences:
        before_start, before_end = before_pos[start], before_pos[end]
        after_start, after_end = after_pos[start], after_pos[end]

        # Sanity check:
        assert before_aligned[start:end].replace('-', '') == before_seq[before_start:before_end]
        assert after_aligned[start:end].replace('-', '') == after_seq[after_start:after_end]

        # Add 1 to starts to convert from 0-based exclusive ranges to 1-based inclusive ranges.
        before_label = f'{before_name} {before_start+1}-{before_end}:'
        after_label = f'{after_name} {after_start+1}-{after_end}:'
        assert len(before_label) <= longest_label
        assert len(after_label) <= longest_label
        before_label = before_label.rjust(longest_label)
        after_label = after_label.rjust(longest_label)

        print(f'{before_label}', before_aligned[start:end])
        print(f'{after_label}', after_aligned[start:end])
        print(' ' * longest_label, differences[start:end])
        print()
    log()


def make_diff_ranges(diff_pos, padding, merge, aligned_len):
    diff_ranges = []
    last_diff_pos = None
    for p in diff_pos:
        start = max(0, p-padding)
        end = min(aligned_len, p+padding+1)
        if not last_diff_pos:  # this is the first diff
            diff_ranges.append((start, end))
        elif p - last_diff_pos <= merge:   # this diff is close to the previous diff
            prev_start = diff_ranges[-1][0]
            diff_ranges.pop()
            diff_ranges.append((prev_start, end))
        else:   # this diff is far from the previous diff
            diff_ranges.append((start, end))
        last_diff_pos = p
    return diff_ranges


def get_aligned_seqs(before_seq, after_seq, aligner):
    cigar = get_cigar(before_seq, after_seq, aligner)
    expanded_cigar = get_expanded_cigar(cigar)
    i, j = 0, 0
    before_aligned, after_aligned, differences = [], [], []
    before_positions, after_positions, diff_positions = [], [], []
    for c in expanded_cigar:
        before_positions.append(i)
        after_positions.append(j)
        if c == 'M':
            b_1 = before_seq[i]
            b_2 = after_seq[j]
            if b_1 == b_2:
                diff = ' '
            else:
                diff = '*'
                diff_positions.append(len(differences))
            i += 1
            j += 1
        elif c == '=':
            b_1 = before_seq[i]
            b_2 = after_seq[j]
            diff = ' '
            i += 1
            j += 1
            assert b_1 == b_2
        elif c == 'X':
            b_1 = before_seq[i]
            b_2 = after_seq[j]
            diff = '*'
            diff_positions.append(len(differences))
            i += 1
            j += 1
            assert b_1 != b_2
        elif c == 'I':
            b_1 = before_seq[i]
            b_2 = '-'
            diff = '*'
            diff_positions.append(len(differences))
            i += 1
        elif c == 'D':
            b_1 = '-'
            b_2 = after_seq[j]
            diff = '*'
            diff_positions.append(len(differences))
            j += 1
        else:
            assert False
        before_aligned.append(b_1)
        after_aligned.append(b_2)
        differences.append(diff)
    assert i == len(before_seq)
    assert j == len(after_seq)
    before_aligned = ''.join(before_aligned)
    after_aligned = ''.join(after_aligned)
    differences = ''.join(differences)
    for p in diff_positions:
        assert differences[p] == '*'
    assert before_aligned.replace('-', '') == before_seq
    assert after_aligned.replace('-', '') == after_seq
    return before_aligned, after_aligned, differences, \
        before_positions, after_positions, diff_positions


def get_cigar(before_seq, after_seq, aligner):
    if aligner == 'mappy':
        return get_cigar_with_mappy(before_seq, after_seq)
    elif aligner == 'edlib':
        return get_cigar_with_edlib(before_seq, after_seq)
    else:
        assert False


def get_cigar_with_mappy(before_seq, after_seq):
    a = mappy.Aligner(seq=after_seq, preset='map-ont')
    for result in a.map(before_seq):
        full_length_query = (result.q_st == 0 and result.q_en == len(before_seq))
        full_length_ref = (result.r_st == 0 and result.r_en == len(after_seq))
        if full_length_query and full_length_ref:
            return result.cigar_str
    quit_with_error('Error: mappy alignment failed, try using --aligner edlib')


def get_cigar_with_edlib(before_seq, after_seq):
    result = edlib.align(before_seq, after_seq, mode='NW', task='path')
    return result['cigar']


def get_expanded_cigar(cigar):
    expanded_cigar = []
    cigar_parts = re.findall(r'\d+[IDX=M]', cigar)
    for p in cigar_parts:
        size = int(p[:-1])
        letter = p[-1]
        expanded_cigar.append(letter * size)
    return ''.join(expanded_cigar)


def finished_message(start_time):
    section_header('Finished!')
    time_to_run = datetime.datetime.now() - start_time
    log(f'Time to run: {time_to_run}')
    log()


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)
    unknown_file = open(str(filename), 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('\nError: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('\nError: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


def load_fasta(fasta_filename):
    fasta_seqs = []
    with get_open_func(fasta_filename)(fasta_filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs.append((name.split()[0], ''.join(sequence)))
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line.upper())
        if name:
            fasta_seqs.append((name.split()[0], ''.join(sequence)))
    return fasta_seqs


class MyParser(argparse.ArgumentParser):
    """
    This subclass of ArgumentParser changes the error messages, such that if a command is run with
    no other arguments, it will display the help text. If there is a different error, it will give
    the normal response (usage and error).
    """
    def error(self, message):
        if len(sys.argv) == 2:  # if a command was given but nothing else
            self.print_help(file=sys.stderr)
            sys.exit(2)
        else:
            super().error(message)


class MyHelpFormatter(argparse.HelpFormatter):
    """
    This is a custom formatter class for argparse. It allows for some custom formatting,
    in particular for the help texts with multiple options (like bridging mode and verbosity level).
    http://stackoverflow.com/questions/3853722
    """
    def __init__(self, prog):
        terminal_width = shutil.get_terminal_size().columns
        os.environ['COLUMNS'] = str(terminal_width)
        max_help_position = min(max(24, terminal_width // 3), 40)
        self.colours = get_colours_from_tput()
        super().__init__(prog, max_help_position=max_help_position)

    def _get_help_string(self, action):
        """
        Override this function to add default values, but only when 'default' is not already in the
        help text.
        """
        help_text = action.help
        if action.default != argparse.SUPPRESS and action.default is not None:
            if 'default' not in help_text.lower():
                help_text += ' (default: {})'.format(action.default)
            elif 'default: DEFAULT' in help_text:
                help_text = help_text.replace('default: DEFAULT',
                                              'default: {}'.format(action.default))
        return help_text

    def start_section(self, heading):
        """
        Override this method to add bold underlining to section headers.
        """
        if self.colours > 1:
            heading = BOLD + heading + END_FORMATTING
        super().start_section(heading)

    def _split_lines(self, text, width):
        """
        Override this method to add special behaviour for help texts that start with:
          'R|' - loop text one option per line
        """
        if text.startswith('R|'):
            text_lines = text[2:].splitlines()
            wrapped_text_lines = []
            for line in text_lines:
                if len(line) <= width:
                    wrapped_text_lines.append(line)
                else:
                    wrap_column = 2
                    line_parts = line.split(', ')
                    join = ','
                    current_line = line_parts[0]
                    for part in line_parts[1:]:
                        if len(current_line) + len(join) + 1 + len(part) <= width:
                            current_line += join + ' ' + part
                        else:
                            wrapped_text_lines.append(current_line + join)
                            current_line = ' ' * wrap_column + part
                    wrapped_text_lines.append(current_line)
            return wrapped_text_lines
        else:
            return argparse.HelpFormatter._split_lines(self, text, width)

    def _fill_text(self, text, width, indent):
        if text.startswith('R|'):
            return ''.join(indent + line for line in text[2:].splitlines(keepends=True))
        else:
            return argparse.HelpFormatter._fill_text(self, text, width, indent)

    def _format_action(self, action):
        """
        Override this method to make help descriptions dim.
        """
        # determine the required width and the entry label
        help_position = min(self._action_max_length + 2,
                            self._max_help_position)
        help_width = self._width - help_position
        action_width = help_position - self._current_indent - 2
        action_header = self._format_action_invocation(action)

        # ho nelp; start on same line and add a final newline
        if not action.help:
            tup = self._current_indent, '', action_header
            action_header = '%*s%s\n' % tup
            indent_first = 0

        # short action name; start on the same line and pad two spaces
        elif len(action_header) <= action_width:
            tup = self._current_indent, '', action_width, action_header
            action_header = '%*s%-*s  ' % tup
            indent_first = 0

        # long action name; start on the next line
        else:
            tup = self._current_indent, '', action_header
            action_header = '%*s%s\n' % tup
            indent_first = help_position

        # collect the pieces of the action help
        parts = [action_header]

        # if there was help for the action, add lines of help text
        if action.help:
            help_text = self._expand_help(action)
            help_lines = self._split_lines(help_text, help_width)
            first_line = help_lines[0]
            if self.colours > 8:
                first_line = DIM + first_line + END_FORMATTING
            parts.append('%*s%s\n' % (indent_first, '', first_line))
            for line in help_lines[1:]:
                if self.colours > 8:
                    line = DIM + line + END_FORMATTING
                parts.append('%*s%s\n' % (help_position, '', line))

        # or add a newline if the description doesn't end with one
        elif not action_header.endswith('\n'):
            parts.append('\n')

        # if there are any sub-actions, add their help as well
        for subaction in self._iter_indented_subactions(action):
            parts.append(self._format_action(subaction))

        # return a single string
        return self._join_parts(parts)


def get_colours_from_tput():
    try:
        return int(subprocess.check_output(['tput', 'colors']).decode().strip())
    except (ValueError, subprocess.CalledProcessError, FileNotFoundError, AttributeError):
        return 1


def get_ascii_art():
    ascii_art = (r"  _____        _                       _  _       _     " + '\n' +
                 r" |  __ \      | |                     | |(_)     | |    " + '\n' +
                 r" | |__) |___  | | _   _  _ __    ___  | | _  ___ | |__  " + '\n' +
                 r" |  ___// _ \ | || | | || '_ \  / _ \ | || |/ __|| '_ \ " + '\n' +
                 r" | |   | (_) || || |_| || |_) || (_) || || |\__ \| | | |" + '\n' +
                 r" |_|    \___/ |_| \__, || .__/  \___/ |_||_||___/|_| |_|" + '\n' +
                 r"                   __/ || |                             " + '\n' +
                 r"                  |___/ |_|                             " + '\n')
    return ascii_art


def log(message='', end='\n'):
    print(message, file=sys.stderr, flush=True, end=end)


def section_header(text):
    log()
    time = get_timestamp()
    time_str = dim('(' + time + ')')
    header = bold_yellow_underline(text)
    print(header + ' ' + time_str, file=sys.stderr, flush=True)


END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
YELLOW = '\033[93m'
DIM = '\033[2m'


def bold_yellow_underline(text):
    return YELLOW + BOLD + UNDERLINE + text + END_FORMATTING


def dim(text):
    return DIM + text + END_FORMATTING


def explanation(text, indent_size=4):
    text = ' ' * indent_size + text
    terminal_width, _ = get_terminal_size_stderr()
    for line in textwrap.wrap(text, width=terminal_width - 1):
        log(dim(line))
    log()


def quit_with_error(text):
    terminal_width, _ = get_terminal_size_stderr()
    log()
    for line in textwrap.wrap(text, width=terminal_width - 1):
        log(line)
    log()
    sys.exit()


def get_terminal_size_stderr(fallback=(80, 24)):
    """
    Unlike shutil.get_terminal_size, which looks at stdout, this looks at stderr.
    """
    try:
        size = os.get_terminal_size(sys.__stderr__.fileno())
    except (AttributeError, ValueError, OSError):
        size = os.terminal_size(fallback)
    return size


def get_timestamp():
    return '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())


def check_python_version():
    if sys.version_info.major < 3 or sys.version_info.minor < 6:
        sys.exit('\nError: polypolish_human_readable.py requires Python 3.6 or later')


if __name__ == '__main__':
    main()


# Unit tests for Pytest
# =====================

def test_get_longest_label():
    before = [('seq1', 'ACGT'), ('seq2', 'ACGTACGTACGTACGT')]
    after = [('seq1_polished', 'ACGT'), ('seq2_polished', 'ACGT')]
    assert get_longest_label(before, after) == 20


def test_get_expanded_cigar():
    assert get_expanded_cigar('5=') == '====='
    assert get_expanded_cigar('3=2X4=1I6=3D3=') == '===XX====I======DDD==='


def test_make_diff_ranges():
    diff_pos = [100, 110]
    padding = 10
    merge = 20
    aligned_len = 200
    diff_ranges = make_diff_ranges(diff_pos, padding, merge, aligned_len)
    assert diff_ranges == [(90, 121)]

    diff_pos = [100, 120]
    padding = 10
    merge = 20
    aligned_len = 200
    diff_ranges = make_diff_ranges(diff_pos, padding, merge, aligned_len)
    assert diff_ranges == [(90, 131)]

    diff_pos = [100, 121]
    padding = 10
    merge = 20
    aligned_len = 200
    diff_ranges = make_diff_ranges(diff_pos, padding, merge, aligned_len)
    assert diff_ranges == [(90, 111), (111, 132)]

    diff_pos = [100, 150]
    padding = 10
    merge = 20
    aligned_len = 200
    diff_ranges = make_diff_ranges(diff_pos, padding, merge, aligned_len)
    assert diff_ranges == [(90, 111), (140, 161)]

    diff_pos = [2, 195]
    padding = 10
    merge = 20
    aligned_len = 200
    diff_ranges = make_diff_ranges(diff_pos, padding, merge, aligned_len)
    assert diff_ranges == [(0, 13), (185, 200)]
