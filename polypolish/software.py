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

import subprocess
import sys

from .log import log


def check_python():
    major, minor, micro = sys.version_info.major, sys.version_info.minor, sys.version_info.micro
    python_version = f'v{major}.{minor}.{micro}'
    good_version = (major >= 3 and minor >= 6)
    if good_version:
        log(f'  Python:   {python_version}')
    else:
        log(f'  Python:   {python_version}')
        sys.exit('\nError: Polypolish requires Python 3.6 or later')


def check_minimap2():
    try:
        output = subprocess.check_output(['minimap2', '--version'], stderr=subprocess.STDOUT)
    except FileNotFoundError:
        sys.exit('\nError: unable to find minimap2 - make sure that minimap2 is installed and '
                 'available on the path, then try again.')
    except subprocess.CalledProcessError:
        sys.exit('\nError: unable to determine minimap2 version - make sure that minimap2 is '
                 'correctly installed, then try again.')
    output = output.decode().strip()
    log(f'  minimap2: v{output}')
