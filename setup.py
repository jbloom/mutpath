"""Setup script for ``mutpath``.

This script uses ``distutils``, the standard python mechanism for installing
packages. To build and install the package, use the following
commands::

    python setup.py build
    python setup.py install

If the user does not have permissions to write to the install directory,
the last command may need to be replaced by::

    sudo python setup.py install

In order for plotting to be enabled, ``pylab`` and ``matplotlib`` must be
installed and available. If they are not available, this script prints a
warning indicating that fact.

Written by Jesse Bloom.
"""


import sys
import os
from distutils.core import setup
from distutils.core import Extension
from distutils.core import Command



# list of C extensions
ctree = Extension('mutpath.ctree', sources=['src/ctree.c'])


# main setup command
setup(
    name = 'mutpath', 
    version = '0.1', 
    author = 'Jesse D. Bloom', 
    author_email = 'jbloom@fhcrc.org', 
    url = 'https://github.com/jbloom/mutpath', 
    description = 'Builds mutational paths through sequence space using BEAST.',
    classifiers = [
        'Intended Audience :: Science/Research',
        'License :: Free for non-commercial use',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    platforms = 'Tested on Mac OS X.',
    packages = ['mutpath'],
    package_dir = {'mutpath':'src'},
    ext_modules = [ctree],
    scripts = [
            'scripts/mutpath_parse_to_beastxml.py',
            'scripts/mutpath_compact_trees.py',
            'scripts/mutpath_get_paths.py',
            'scripts/mutpath_annotate_tree.py',
            'scripts/mutpath_make_digraph.py',
            ],
)
