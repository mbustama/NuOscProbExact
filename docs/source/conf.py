# -*- coding: utf-8 -*-
r"""Sphinx configuration for the NuOscProbExact documentation.

Build the HTML documentation with::

    pip install -e ".[docs]"
    sphinx-build -b html docs/source docs/build

The library modules live in ``src/`` as top-level modules rather than
inside a package, so that directory is put on ``sys.path`` below and
``autodoc`` imports them by their bare names, exactly as user code does.
"""

import os
import sys

sys.path.insert(0, os.path.abspath('../../src'))

project = 'NuOscProbExact'
copyright = '2019-2026, Mauricio Bustamante'
author = 'Mauricio Bustamante'
release = '1.1.0'
version = '1.1'

extensions = [
    'sphinx.ext.autodoc',       # Pull the API reference from the docstrings
    'sphinx.ext.mathjax',       # Render the LaTeX math in the docstrings
    'sphinx.ext.viewcode',      # Link API entries to highlighted source
    'sphinx.ext.intersphinx',
    'numpydoc',                 # Parse the numpydoc-style docstrings
    'myst_parser',              # Lets changelog.rst .. include:: CHANGELOG.md
]

# numpydoc would otherwise try to document the members of every class it
# meets; the library exposes only module-level routines.
numpydoc_show_class_members = False

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable', None),
}

master_doc = 'index'
templates_path = ['_templates']
exclude_patterns = ['_build']

html_theme = 'alabaster'
html_static_path = []
