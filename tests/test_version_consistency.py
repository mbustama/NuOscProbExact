# -*- coding: utf-8 -*-
r"""Checks that the version agrees everywhere it is stated or implied.

The version is written in exactly one place, ``pyproject.toml``, because
that is the file :mod:`build` reads and the one ``publish.yml`` stamps
for a TestPyPI rehearsal.  Everything else either derives it or refers
to it, and this module is what keeps "refers to it" honest.

Three things cannot be derived and so are checked instead.

The **changelog heading**, because a release is a claim about what
changed, not a string: a version bumped in ``pyproject.toml`` with no
entry beside it ships a number nobody can look up.

The **documentation's copy**, which is derived in ``conf.py`` rather
than typed --- this asserts the derivation still works, since a change
to the version line's spelling would silently fall back to installed
metadata and could then report a *different* version than the tree
being documented.

The ``versionadded`` and ``versionchanged`` **directives**, which name
releases across all ten modules.  A directive naming a version that the
changelog does not record sends a reader looking for an entry that does
not exist, and the audit that produced this file found the reverse case
too: a release whose behaviour changed and said so nowhere.
"""

import pathlib
import re
import runpy

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
PYPROJECT = ROOT/'pyproject.toml'
CHANGELOG = ROOT/'CHANGELOG.md'
CONF = ROOT/'docs'/'source'/'conf.py'
SRC = ROOT/'src'


def declared_version():
    r"""Returns the version from ``pyproject.toml``, the only source."""
    found = re.search(r'^version = "([^"]+)"', PYPROJECT.read_text(), re.M)
    assert found, 'pyproject.toml has no [project] version line'

    return found.group(1)


def changelog_versions():
    r"""Returns every version with a changelog heading, newest first."""
    return re.findall(r'^## \[(\d+\.\d+\.\d+)\]', CHANGELOG.read_text(), re.M)


def test_the_changelog_leads_with_the_declared_version():
    r"""A bump without an entry ships a number nobody can look up."""
    declared = declared_version()
    newest = changelog_versions()[0]

    assert newest == declared, (
        'pyproject.toml declares %s but the newest CHANGELOG.md heading is '
        '%s; a release needs both' % (declared, newest))


def test_the_changelog_is_ordered_and_free_of_duplicates():
    r"""Headings descend, and no version appears twice.

    Editing this file by anchoring on the previous version heading has
    deleted it before now, which shows up here as a missing version
    rather than as a merge conflict.
    """
    versions = changelog_versions()
    keys = [tuple(int(part) for part in v.split('.')) for v in versions]

    assert len(set(versions)) == len(versions), (
        'duplicate CHANGELOG.md headings: %s'
        % sorted({v for v in versions if versions.count(v) > 1}))
    assert keys == sorted(keys, reverse=True), (
        'CHANGELOG.md headings are not in descending order: %s' % versions)
    assert versions[-1] == '1.0.0', (
        'the oldest CHANGELOG.md heading is %s, not 1.0.0; an entry has '
        'been lost from the end' % versions[-1])


def test_the_documentation_derives_the_version_it_reports():
    r"""``conf.py`` computes it, and computes it correctly.

    Asserted by executing ``conf.py`` rather than by reading it, because
    the failure this guards against is the derivation silently falling
    back to installed metadata and reporting a version other than the
    one in the tree being documented.
    """
    declared = declared_version()
    namespace = runpy.run_path(str(CONF))

    assert namespace['release'] == declared
    assert namespace['version'] == '.'.join(declared.split('.')[:2])

    # And it is derived, not typed: neither value is a literal here
    text = CONF.read_text()
    assert "release = '%s'" % declared not in text
    assert re.search(r"^version = '\d", text, re.M) is None


def test_every_version_directive_names_a_released_version():
    r"""``versionadded`` and ``versionchanged`` point at real entries."""
    released = set(changelog_versions())

    unknown = {}
    for module in sorted(SRC.glob('*.py')):
        for version in re.findall(
                r'\.\. version(?:added|changed):: (\d+\.\d+\.\d+)',
                module.read_text()):
            if version not in released:
                unknown.setdefault(version, set()).add(module.name)

    assert not unknown, (
        'these versions are named in a directive but have no CHANGELOG.md '
        'entry: %s' % {v: sorted(m) for v, m in unknown.items()})


def test_no_directive_announces_a_version_that_has_not_shipped():
    r"""Nothing claims to be newer than the declared version.

    A ``versionadded:: 1.12.0`` written while ``pyproject.toml`` still
    says 1.11.0 documents a promise rather than a fact, and reads as
    released once the docs are built.
    """
    declared = tuple(int(part) for part in declared_version().split('.'))

    ahead = {}
    for module in sorted(SRC.glob('*.py')):
        for version in re.findall(
                r'\.\. version(?:added|changed):: (\d+\.\d+\.\d+)',
                module.read_text()):
            if tuple(int(p) for p in version.split('.')) > declared:
                ahead.setdefault(version, set()).add(module.name)

    assert not ahead, (
        'these directives name a version newer than the declared %s: %s'
        % (declared_version(), {v: sorted(m) for v, m in ahead.items()}))


@pytest.mark.parametrize('spelling', ['__version__', 'VERSION'])
def test_the_version_is_not_restated_in_the_library(spelling):
    r"""``src/`` states the version nowhere, so it cannot drift.

    A module-level ``__version__`` is the usual place for a fourth copy
    to appear.  Callers who want it should ask the installed
    distribution, which cannot disagree with what was installed.
    """
    for module in sorted(SRC.glob('*.py')):
        text = module.read_text()
        assert not re.search(r'^%s\s*=' % spelling, text, re.M), (
            '%s defines %s; the version belongs in pyproject.toml alone'
            % (module.name, spelling))
