#!/usr/bin/env python3
"""
This is the Polypolish installation script. Assuming you're in the same directory, it can be run
like this: `python3 setup.py install`, or (probably better) like this: `pip3 install .`

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

from setuptools import setup
from Cython.Build import cythonize
from distutils.extension import Extension


def readme():
    with open('README.md') as f:
        return f.read()


# Get the program version from another file.
__version__ = '0.0.0'
exec(open('polypolish/version.py').read())


extensions = [Extension(name='Polypolish', sources=['polypolish/*.pyx'], extra_compile_args=['-w'])]

setup(name='Polypolish',
      version=__version__,
      description='Polypolish: a hybrid short-read aligner',
      long_description=readme(),
      long_description_content_type='text/markdown',
      url='https://github.com/rrwick/Polypolish',
      author='Ryan Wick',
      author_email='rrwick@gmail.com',
      license='GPLv3',
      packages=['polypolish'],
      install_requires=['Cython', 'pytest'],
      entry_points={"console_scripts": ['polypolish = polypolish.__main__:main']},
      include_package_data=True,
      ext_modules=cythonize(extensions, compiler_directives={'language_level': "3"}),
      zip_safe=False,
      python_requires='>=3.6')
