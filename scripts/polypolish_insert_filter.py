#!/usr/bin/env python3


import argparse
import collections
import datetime
import os
import pathlib
import re
import shutil
import subprocess
import sys
import textwrap

__version__ = '0.5.0'


def main():
    args = parse_args()
    check_inputs(args)
    start_time = starting_message(args)

    alignments, before_count = load_alignments(args.in1, args.in2)
    distribution = get_insert_size_distribution(alignments, args.orientation)

    after_count = 0  # TEMP
    finished_message(start_time, before_count, after_count)


def parse_args():
    description = 'R|' + get_ascii_art() + '\n' + \
                  'Polypolish insert size alignment filter\n' + \
                  'github.com/rrwick/Polypolish'
    parser = MyParser(description=description, formatter_class=MyHelpFormatter, add_help=False)

    input_args = parser.add_argument_group('Inputs')
    input_args.add_argument('--in1', type=str, required=True,
                            help='Input SAM file (first read in pairs)')
    input_args.add_argument('--in2', type=str, required=True,
                            help='Input SAM file (first second in pairs)')

    output_args = parser.add_argument_group('Outputs')
    output_args.add_argument('--out1', type=str, required=True,
                             help='Output SAM file (first read in pairs)')
    output_args.add_argument('--out2', type=str, required=True,
                             help='Output SAM file (first second in pairs)')

    setting_args = parser.add_argument_group('Settings')
    setting_args.add_argument('--orientation', choices=['FR', 'RF', 'FF', 'RR', 'auto'],
                              default='auto',
                              help='Expected pair orientation (default: determine automatically)')

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version',
                           version='Polypolish insert size alignment filter v' + __version__,
                           help="Show program's version number and exit")

    # If no arguments were used, print the base-level help which lists possible commands.
    if len(sys.argv) == 1:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    return parser.parse_args()


def check_inputs(args):
    pass
    # TODO: make sure that none of the required inputs are the same


def starting_message(args):
    section_header('Polypolish insert size alignment filter')
    explanation('This script is a pre-processing filter that can be run on SAM alignments before '
                'they are used in Polypolish. It looks at each read pair and discards alignments '
                'that do not seem to be part of a concordant pair. This serves to reduce the '
                'number of alignments (especially near the edges of repeats), making ')

    log('Input alignments:')
    log(f'  {args.in1}')
    log(f'  {args.in2}')
    log()
    log('Output alignments:')
    log(f'  {args.out1}')
    log(f'  {args.out2}')
    log()
    log('Settings:')
    log(f'  --orientation {args.orientation}')
    log()
    check_python_version()
    return datetime.datetime.now()


def finished_message(start_time, before_count, after_count):
    section_header('Finished!')
    log(f'Alignments before filtering: {before_count}')
    log(f'Alignments after filtering:  {after_count}')
    time_to_run = datetime.datetime.now() - start_time
    log(f'Time to run: {time_to_run}')
    log()


def load_alignments(sam_1, sam_2):
    section_header('Loading alignments')
    alignments = collections.defaultdict(list)
    load_alignments_one_file(sam_1, alignments, '_1')
    load_alignments_one_file(sam_2, alignments, '_2')
    log()
    count = sum(len(a) for a in alignments.values())
    return alignments, count


def load_alignments_one_file(sam_filename, alignments, read_name_suffix):
    log(f'{sam_filename}: ', end='')
    alignment_count = 0
    read_names = set()
    with open(sam_filename, 'rt') as sam_file:
        for line in sam_file:
            if line.startswith('@'):  # header line
                continue
            alignment = Alignment(line)
            if alignment.is_aligned():
                alignment.read_name += read_name_suffix
                alignments[alignment.read_name].append(alignment)
                alignment_count += 1
                read_names.add(alignment.read_name)
    log(f'{alignment_count:,} alignments from {len(read_names):,} reads')
    return alignments


def get_insert_size_distribution(alignments, correct_orientation):
    section_header('Finding insert size distribution')
    explanation('Read pairs with exactly one alignment per read are used to determine the '
                'orientation and insert size distribution for the read set.')
    insert_sizes = {'FR': [], 'RF': [], 'FF': [], 'RR': []}
    for name_1, alignments_1 in alignments.items():
        if not name_1.endswith('_1'):
            continue
        name_2 = name_1[:-2] + '_2'
        alignments_2 = alignments[name_2]
        if len(alignments_1) == 1 and len(alignments_2) == 1:
            a_1, a_2 = alignments_1[0], alignments_2[0]
            if a_1.ref_name == a_2.ref_name:
                orientation = get_orientation(a_1, a_2)
                insert_size = get_insert_size(a_1, a_2)
                insert_sizes[orientation].append(insert_size)
    if correct_orientation == 'auto':
        correct_orientation = auto_determine_orientation(insert_sizes)

    return distribution


def get_orientation(alignment_1, alignment_2):
    strand_1 = 'F' if alignment_1.is_on_forward_strand() else 'R'
    strand_2 = 'F' if alignment_2.is_on_forward_strand() else 'R'
    if strand_1 != strand_2:  # on different strands
        if alignment_1.ref_start < alignment_2.ref_start:
            return strand_1 + strand_2
        else:
            return strand_2 + strand_1
    elif strand_1 == 'F':  # both on forward strand
        if alignment_1.ref_start < alignment_2.ref_start:
            return 'FF'
        else:
            return 'RR'
    elif strand_1 == 'R':  # both on reverse strand
        if alignment_2.ref_start < alignment_1.ref_start:
            return 'FF'
        else:
            return 'RR'
    else:
        assert False


def auto_determine_orientation(insert_sizes):
    log('Automatically determining read orientation:')
    max_count = max(len(i) for i in insert_sizes.values())
    orientations = []
    for orientation in ['FR', 'RF', 'FF', 'RR']:
        count = len(insert_sizes[orientation])
        log(f'  {orientation}: {count} pairs')
        if count == max_count:
            orientations.append(orientation)
    if len(orientations) == 1:
        correct_orientation = orientations[0]
        log(f'  Correct orientation: {correct_orientation}')
        return correct_orientation
    quit_with_error('\nError: could not automatically determine read pair orientation')


def get_insert_size(alignment_1, alignment_2):
    insert_start = min(alignment_1.ref_start, alignment_1.ref_end,
                       alignment_2.ref_start, alignment_2.ref_end)
    insert_end = max(alignment_1.ref_start, alignment_1.ref_end,
                     alignment_2.ref_start, alignment_2.ref_end)
    return insert_end - insert_start


class Alignment(object):

    def __init__(self, sam_line):
        parts = sam_line.strip().split('\t')
        if len(parts) < 11:
            quit_with_error('\nError: alignment file does not seem to be in SAM format')

        self.read_name = parts[0]
        self.sam_flags = int(parts[1])
        self.ref_name = parts[2]
        self.ref_start = int(parts[3]) - 1  # switch from SAM's 1-based to Python's 0-based
        cigar = parts[5]
        self.ref_end = get_ref_end(self.ref_start, cigar)

    def __repr__(self):
        return f'{self.read_name}:{self.ref_name}:{self.ref_start}-{self.ref_end}:{self.mismatches}'

    def is_aligned(self):
        return not self.has_flag(4)

    def is_secondary(self):
        return self.has_flag(256)

    def is_on_forward_strand(self):
        return not self.has_flag(16)

    def has_flag(self, flag):
        return bool(self.sam_flags & flag)

    def starts_and_ends_with_match(self):
        cigar_parts = re.findall(r'\d+[MIDNSHP=X]', self.cigar)
        first_part, last_part = cigar_parts[0], cigar_parts[-1]
        return first_part[-1] == 'M' and last_part[-1] == 'M'


def get_ref_end(ref_start, cigar):
    ref_end = ref_start
    cigar_parts = re.findall(r'\d+[MIDNSHP=X]', cigar)
    for p in cigar_parts:
        size = int(p[:-1])
        letter = p[-1]
        if letter == 'M' or letter == 'D' or letter == 'N' or letter == '=' or letter == 'X':
            ref_end += size
    return ref_end


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
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
YELLOW = '\033[93m'
DIM = '\033[2m'


def bold(text):
    return BOLD + text + END_FORMATTING


def bold_yellow(text):
    return YELLOW + BOLD + text + END_FORMATTING


def bold_yellow_underline(text):
    return YELLOW + BOLD + UNDERLINE + text + END_FORMATTING


def dim(text):
    return DIM + text + END_FORMATTING


def green(text):
    return GREEN + text + END_FORMATTING


def red(text):
    return RED + text + END_FORMATTING


def bold_red(text):
    return RED + BOLD + text + END_FORMATTING


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
        sys.exit('\nError: polypolish_insert_filter.py requires Python 3.6 or later')


if __name__ == '__main__':
    main()


# Unit tests with Pytest
# ======================

def test_get_orientation_1():
    # 1------>
    #            <------2
    pos_1, pos_2 = 100000, 200000
    strand_1, strand_2 = 0, 16
    a_1 = Alignment(f'read_1\t{strand_1}\tref\t{pos_1}\t60\t150M\t*\t0\t0\tACTG\tKKKK')
    a_2 = Alignment(f'read_2\t{strand_2}\tref\t{pos_2}\t60\t150M\t*\t0\t0\tACTG\tKKKK')
    assert get_orientation(a_1, a_2) == 'FR'


def test_get_orientation_2():
    # 2------>
    #            <------1
    pos_1, pos_2 = 200000, 100000
    strand_1, strand_2 = 16, 0
    a_1 = Alignment(f'read_1\t{strand_1}\tref\t{pos_1}\t60\t150M\t*\t0\t0\tACTG\tKKKK')
    a_2 = Alignment(f'read_2\t{strand_2}\tref\t{pos_2}\t60\t150M\t*\t0\t0\tACTG\tKKKK')
    assert get_orientation(a_1, a_2) == 'FR'


def test_get_orientation_3():
    #            1------>
    # <------2
    pos_1, pos_2 = 200000, 100000
    strand_1, strand_2 = 0, 16
    a_1 = Alignment(f'read_1\t{strand_1}\tref\t{pos_1}\t60\t150M\t*\t0\t0\tACTG\tKKKK')
    a_2 = Alignment(f'read_2\t{strand_2}\tref\t{pos_2}\t60\t150M\t*\t0\t0\tACTG\tKKKK')
    assert get_orientation(a_1, a_2) == 'RF'


def test_get_orientation_4():
    # <------1
    #            2------>
    pos_1, pos_2 = 100000, 200000
    strand_1, strand_2 = 16, 0
    a_1 = Alignment(f'read_1\t{strand_1}\tref\t{pos_1}\t60\t150M\t*\t0\t0\tACTG\tKKKK')
    a_2 = Alignment(f'read_2\t{strand_2}\tref\t{pos_2}\t60\t150M\t*\t0\t0\tACTG\tKKKK')
    assert get_orientation(a_1, a_2) == 'RF'


def test_get_orientation_5():
    # 1------>   2------>
    pos_1, pos_2 = 100000, 200000
    strand_1, strand_2 = 0, 0
    a_1 = Alignment(f'read_1\t{strand_1}\tref\t{pos_1}\t60\t150M\t*\t0\t0\tACTG\tKKKK')
    a_2 = Alignment(f'read_2\t{strand_2}\tref\t{pos_2}\t60\t150M\t*\t0\t0\tACTG\tKKKK')
    assert get_orientation(a_1, a_2) == 'FF'


def test_get_orientation_6():
    # <------2   <------1
    pos_1, pos_2 = 200000, 100000
    strand_1, strand_2 = 16, 16
    a_1 = Alignment(f'read_1\t{strand_1}\tref\t{pos_1}\t60\t150M\t*\t0\t0\tACTG\tKKKK')
    a_2 = Alignment(f'read_2\t{strand_2}\tref\t{pos_2}\t60\t150M\t*\t0\t0\tACTG\tKKKK')
    assert get_orientation(a_1, a_2) == 'FF'


def test_get_orientation_7():
    # 2------>   1------>
    pos_1, pos_2 = 200000, 100000
    strand_1, strand_2 = 0, 0
    a_1 = Alignment(f'read_1\t{strand_1}\tref\t{pos_1}\t60\t150M\t*\t0\t0\tACTG\tKKKK')
    a_2 = Alignment(f'read_2\t{strand_2}\tref\t{pos_2}\t60\t150M\t*\t0\t0\tACTG\tKKKK')
    assert get_orientation(a_1, a_2) == 'RR'


def test_get_orientation_8():
    # <------1   <------2
    pos_1, pos_2 = 100000, 200000
    strand_1, strand_2 = 16, 16
    a_1 = Alignment(f'read_1\t{strand_1}\tref\t{pos_1}\t60\t150M\t*\t0\t0\tACTG\tKKKK')
    a_2 = Alignment(f'read_2\t{strand_2}\tref\t{pos_2}\t60\t150M\t*\t0\t0\tACTG\tKKKK')
    assert get_orientation(a_1, a_2) == 'RR'


def test_auto_determine_orientation_1():
    insert_sizes = {'FR': [100, 100, 100], 'RF': [200], 'FF': [300], 'RR': [400]}
    assert auto_determine_orientation(insert_sizes) == 'FR'


def test_auto_determine_orientation_2():
    insert_sizes = {'FR': [100], 'RF': [200, 200, 200], 'FF': [300], 'RR': [400]}
    assert auto_determine_orientation(insert_sizes) == 'RF'


def test_auto_determine_orientation_3():
    insert_sizes = {'FR': [100], 'RF': [200], 'FF': [300, 300, 300], 'RR': [400]}
    assert auto_determine_orientation(insert_sizes) == 'FF'


def test_auto_determine_orientation_4():
    insert_sizes = {'FR': [100], 'RF': [200], 'FF': [300], 'RR': [400, 400, 400]}
    assert auto_determine_orientation(insert_sizes) == 'RR'
