# -*- coding: utf-8 -*-
r"""Checks that the performance figures agree everywhere they are quoted.

The same measurements are stated in five places --- ``README.md``,
``docs/source/index.rst``, ``docs/source/quickstart.rst``,
``docs/source/methodology.rst`` and :mod:`fastkernels` --- because each
has to stand on its own: the README is also the PyPI long description,
and a reader of the methodology page should not have to visit the
landing page to find out what the code costs.

Repeating a measurement is fine.  Repeating it *without a check* is what
went wrong before: 1.8.4 set out to reconcile these numbers and updated
four of the five files, leaving `methodology.rst` contradicting itself
because it states the scalar timing twice; and `index.rst` claimed a
speedup span of "25 to 60 times" three lines above its own table, which
gives 21x to 99x.  Both were found by hand, twice, a release apart.

So the copies stay, and these tests make them agree by construction.  A
figure changed in one place and not the others fails here, naming both
the file that disagrees and the value it should carry.
"""

import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

README = os.path.join(ROOT, 'README.md')
INDEX_RST = os.path.join(ROOT, 'docs', 'source', 'index.rst')
QUICKSTART_RST = os.path.join(ROOT, 'docs', 'source', 'quickstart.rst')
METHODOLOGY_RST = os.path.join(ROOT, 'docs', 'source', 'methodology.rst')
FASTKERNELS = os.path.join(ROOT, 'src', 'fastkernels.py')

# The scalar cost of one probability, in microseconds.  Re-measured in
# 1.8.6; see the changelog entry for the method.  Stated in words on the
# documentation pages and in digits in the README and the module, so both
# spellings are searched for.
SCALAR_3NU_US = 8
SCALAR_2NU_US = 1

WORDS = {1: 'one', 2: 'two', 8: 'eight', 13: 'thirteen', 16: 'sixteen'}

# The 2000-point benchmark table, shared by README.md and index.rst:
# (scan, loop, arrays, arrays with Numba), in milliseconds.  `None` is the
# two-flavor row, where the compiled kernel is deliberately not used.
#
# Held as the strings the documents actually print, not as floats, so that
# a trailing zero is part of what is checked: both tables write "0.20 ms",
# and a float would compare equal to a document that had drifted to "0.2".
BENCHMARKS = [
    ('three flavors, versus baseline', '38', '1.8', '0.31'),
    ('three flavors, versus energy', '34', '1.5', '0.20'),
    ('oscillogram, 100 x 100', '197', '5.3', '0.85'),
    ('two flavors, versus baseline', '6.9', '0.07', None),
]


def read(path):
    r"""Returns the text of ``path``."""
    with open(path) as handle:
        return handle.read()


def _numbers_in(text, pattern):
    r"""Returns every float captured by ``pattern`` in ``text``."""
    return [float(match) for match in re.findall(pattern, text)]


def test_scalar_timings_agree_across_every_document():
    r"""One probability costs the same number wherever it is quoted."""
    digits_3nu = re.compile(r'\b%d\s*(?:us|µs|microsecond)' % SCALAR_3NU_US)
    words_3nu = re.compile(r'\b%s microsecond' % WORDS[SCALAR_3NU_US])

    for path in (README, INDEX_RST, QUICKSTART_RST, METHODOLOGY_RST,
                 FASTKERNELS):
        text = read(path)
        found = digits_3nu.search(text) or words_3nu.search(text)
        assert found, (
            '%s does not state the three-flavor scalar timing as %d '
            'microseconds; every document that quotes it must agree'
            % (os.path.relpath(path, ROOT), SCALAR_3NU_US))


def test_no_document_still_carries_a_superseded_scalar_timing():
    r"""The figures replaced in 1.8.6 do not survive anywhere."""
    superseded = [13, 16]
    for path in (README, INDEX_RST, QUICKSTART_RST, METHODOLOGY_RST,
                 FASTKERNELS):
        text = read(path)
        for value in superseded:
            for spelling in (r'\b%d\s*(?:us|µs|microsecond)' % value,
                             r'\b%s microsecond' % WORDS[value]):
                assert not re.search(spelling, text), (
                    '%s still quotes %d microseconds for one probability; '
                    'the measured value is %d'
                    % (os.path.relpath(path, ROOT), value, SCALAR_3NU_US))


def test_benchmark_table_agrees_between_readme_and_index():
    r"""The 2000-point table is the same in both documents."""
    readme = read(README)
    index = read(INDEX_RST)

    for scan, loop, arrays, numba in BENCHMARKS:
        for value in (loop, arrays, numba):
            if value is None:
                continue
            # Written as "38 ms" in the README table and in the rst
            # list-table alike, so one pattern serves both.
            pattern = r'%s ms\b' % re.escape(value)
            assert re.search(pattern, readme), (
                'README.md is missing "%s ms", from the "%s" row of the '
                'benchmark table' % (value, scan))
            assert re.search(pattern, index), (
                'docs/source/index.rst is missing "%s ms", from the "%s" '
                'row of the benchmark table' % (value, scan))


def test_quoted_speedup_span_matches_the_table():
    r"""The "20 to 90 times" claim is what the table actually shows."""
    ratios = [float(loop)/float(arrays) for _, loop, arrays, _ in BENCHMARKS]
    low, high = min(ratios), max(ratios)

    # The prose rounds outward to the nearest ten, which is the honest way
    # to quote a span whose ends are 21 and 99.
    assert 20 <= low <= 30, 'lowest speedup in the table is %.0fx' % low
    assert 90 <= high <= 100, 'highest speedup in the table is %.0fx' % high

    for path in (README, INDEX_RST):
        text = read(path)
        assert re.search(r'20\s*(?:to|–|-)\s*90', text), (
            '%s should quote the array speedup as 20 to 90 times, which is '
            'what its own benchmark table gives (%.0fx to %.0fx)'
            % (os.path.relpath(path, ROOT), low, high))


def test_two_flavor_row_is_marked_as_not_using_the_backend():
    r"""The blank cell is explained, not left as an apparent omission."""
    for path in (README, INDEX_RST):
        text = read(path)
        assert 'not used' in text, (
            '%s should say the compiled backend is "not used" for the '
            'two-flavor row rather than leaving it blank'
            % os.path.relpath(path, ROOT))
