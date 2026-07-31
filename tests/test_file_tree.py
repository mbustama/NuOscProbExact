# -*- coding: utf-8 -*-
r"""Tests that the file tree stays accurate and stays in one piece.

The repository layout is documented twice --- in ``README.md`` and in
``docs/source/installation.rst`` --- because the two are read in
different places.  Two copies drift, so these tests check that they
agree with each other and that both agree with what is actually on
disk.
"""

import os
import re
import subprocess

import pytest


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
README = os.path.join(ROOT, 'README.md')
INSTALL_RST = os.path.join(ROOT, 'docs', 'source', 'installation.rst')

# Entries that are directories holding generated output, so the tree
# lists the directory but not the files inside it.
GENERATED_DIRS = {'fig/'}


def extract_tree(path):
    r"""Returns the file-tree block of ``path`` as a list of lines."""
    with open(path) as handle:
        text = handle.read()
    lines = text.split('\n')
    start = next(i for i, line in enumerate(lines)
                 if line.strip() == 'NuOscProbExact/')
    tree = [lines[start].strip()]
    for line in lines[start+1:]:
        stripped = line.strip()
        if not stripped or not re.match(r'^[─-╿ ]*[─-╿]',
                                        stripped):
            break
        tree.append(stripped)
    return tree


def tree_names(tree):
    r"""Returns the bare file and directory names listed in ``tree``."""
    names = set()
    for line in tree[1:]:
        # Strip the box-drawing prefix and any trailing comment
        entry = re.sub(r'^[─-╿\s]+', '', line)
        entry = entry.split('#')[0].strip()
        if entry:
            names.add(entry)
    return names


def test_readme_and_docs_trees_are_identical():
    r"""The tree in README.md and in installation.rst do not drift."""
    readme_tree = extract_tree(README)
    docs_tree = extract_tree(INSTALL_RST)
    assert readme_tree == docs_tree, (
        'the file tree in README.md and docs/source/installation.rst have '
        'diverged; update both')


def test_tree_is_not_empty():
    r"""The extraction actually found a tree, rather than nothing."""
    assert len(extract_tree(README)) > 20


def test_every_tracked_file_appears_in_the_tree():
    r"""No tracked file is missing from the documented layout."""
    tracked = subprocess.check_output(
        ['git', 'ls-files'], cwd=ROOT).decode().split()
    listed = tree_names(extract_tree(README))

    missing = []
    for path in tracked:
        name = os.path.basename(path)
        directory = os.path.dirname(path)
        # A file is covered either by its own entry, or by an entry for a
        # directory whose contents the tree deliberately summarises.
        if name in listed:
            continue
        if any(d.rstrip('/') == directory.split('/')[0]
               for d in GENERATED_DIRS):
            continue
        missing.append(path)

    assert not missing, 'not listed in the file tree: %s' % ', '.join(missing)


def test_every_tree_entry_exists_on_disk():
    r"""The tree does not list files that are not there."""
    listed = tree_names(extract_tree(README))
    names_on_disk = set()
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames
                       if d not in {'.git', '__pycache__', 'build',
                                    '_build'}
                       and not d.endswith('.egg-info')]
        for name in filenames + dirnames:
            names_on_disk.add(name)
            names_on_disk.add(name + '/')

    stale = [entry for entry in listed
             if entry not in names_on_disk
             and entry.rstrip('/') not in names_on_disk]
    assert not stale, 'listed in the file tree but absent: %s' % ', '.join(
        sorted(stale))
