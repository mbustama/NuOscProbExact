# -*- coding: utf-8 -*-
r"""Sphinx configuration for the NuOscProbExact documentation.

Build the HTML documentation with::

    pip install -r docs/requirements.txt
    cd docs && make html

The result lands in ``docs/build/html``.

The library modules live in ``src/`` as top-level modules rather than
inside a package, so that directory is put on ``sys.path`` below and
``autodoc`` imports them by their bare names, exactly as user code does.
"""

import os
import sys

sys.path.insert(0, os.path.abspath('../../src'))

# -- Project information -----------------------------------------------------

project = 'NuOscProbExact'
copyright = '2019-2026, Mauricio Bustamante'
author = 'Mauricio Bustamante'
release = '1.3.0'
version = '1.3'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',       # API reference, generated from the docstrings
    'sphinx.ext.mathjax',       # Renders the LaTeX math in the docstrings
    'sphinx.ext.viewcode',      # Links API entries to highlighted source
    'sphinx.ext.intersphinx',   # Cross-links to the Python and NumPy docs
    'numpydoc',                 # Parses the numpydoc-style docstrings
    'sphinx_copybutton',        # Copy-to-clipboard button on code blocks
    'sphinxcontrib.bibtex',     # References page (refs.bib)
    'myst_parser',              # Lets changelog.rst .. include:: CHANGELOG.md
]

# numpydoc would otherwise try to document the members of every class it
# meets; the library exposes only module-level routines and constants.
numpydoc_show_class_members = False

# The docstrings are numpydoc, not Google style.
napoleon_google_docstring = False
napoleon_numpy_docstring = True

bibtex_bibfiles = ['refs.bib']

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable', None),
    'scipy': ('https://docs.scipy.org/doc/scipy', None),
}

# Keep the source order of the routines, which groups them the way the
# method is derived, rather than sorting them alphabetically.
autodoc_member_order = 'bysource'

master_doc = 'index'
templates_path = ['_templates']
exclude_patterns = ['_build']

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = []
html_title = 'NuOscProbExact %s' % version
