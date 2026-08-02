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
import pathlib
import re
import sys

sys.path.insert(0, os.path.abspath('../../src'))

# -- Project information -----------------------------------------------------

project = 'NuOscProbExact'
copyright = '2019-2026, Mauricio Bustamante'
author = 'Mauricio Bustamante'
# The version is not written here.  It lives in pyproject.toml, which is
# what `python -m build` reads and what publish.yml stamps for a TestPyPI
# rehearsal, so that file is the only one entitled to state it.  Keeping a
# second copy here meant two hand-edits per release and no check that they
# agreed -- and the short form below is not a fact at all but a derivation,
# which is the kind of thing that has no business being typed.
#
# Read from the file rather than from package metadata so that the
# documentation builds from a source tree that has not been installed,
# which is what `cd docs && make html` does; metadata is the fallback for a
# build that runs from an installed distribution instead.
def _project_version():
    """Returns the version, from pyproject.toml or installed metadata."""
    pyproject = pathlib.Path(__file__).resolve().parents[2]/'pyproject.toml'
    try:
        found = re.search(r'^version = "([^"]+)"',
                          pyproject.read_text(), re.M)
    except OSError:
        found = None
    if found is not None:
        return found.group(1)

    from importlib.metadata import version as installed_version

    return installed_version('nuoscprobexact')


release = _project_version()
version = '.'.join(release.split('.')[:2])

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
    'jupyter_sphinx',           # Runs the narrative examples at build time
]

# numpydoc would otherwise try to document the members of every class it
# meets; the library exposes only module-level routines and constants.
numpydoc_show_class_members = False

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
html_static_path = ['_static']
html_title = 'NuOscProbExact %s' % version

# Shown at the top of the sidebar, in place of the project name.  The theme
# scales it to the sidebar width, so the source image being larger than it
# needs to be costs nothing but a few hundred kilobytes in the repository.
html_logo = '_static/nuoscprobexact_logo.png'

# `logo_only = False` keeps `html_title` --- and so the version --- visible
# under the logo, rather than letting the image stand alone.  The option names
# here are checked against sphinx_rtd_theme's theme.conf: `display_version`
# looks right and is not one of them, and an unsupported option is a warning,
# which this project's -W build treats as an error.
html_theme_options = {
    'logo_only': False,
}
