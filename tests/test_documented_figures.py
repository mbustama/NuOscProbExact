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

The second half of this module guards the *cross-check* numbers in the
same way, and for a sharper reason: those have drifted three times.  The
frozen nuSQuIDS set grew from seven configurations to eleven, and on each
occasion prose describing it was left behind --- a case count in two test
modules, an agreement range quoted before the stiffest case existed, and
a summary in the notebook.  Every one was found by reading rather than by
the suite, because a test that checks numbers agree says nothing about
the sentences beside them.  These derive both the count and the range
from the data itself.
"""

import json
import os
import re

import numpy as np

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


# The compiled kernel's speedups against the NumPy path, stated in
# `fastkernels`' own module docstring and again in the methodology page's
# table.  Held as the "~19x" spelling both documents use.
#
# This used to hold the two four-flavor rows alone, and the other four
# drifted apart underneath it: the audit before 1.11.0 was released found
# "20 000 energies, three flavors" quoted as ~9x in the module and ~13x on
# the page, and the two-flavor row as ~1.5x against ~1.4x, while the fifth
# row was a different scan entirely in each.  A guard that covers part of a
# table is an invitation to the rest of it.  `fastkernels` is the source
# these are measured in, so it is the one the page is checked against.
SPEEDUP_KERNEL = [('200 000 energies, four flavors', '19'),
                  ('20 000 energies, four flavors', '18'),
                  ('200 000 energies, three flavors', '15'),
                  ('20 000 energies, three flavors', '9'),
                  ('100 x 100 oscillogram', '3.5'),
                  ('200 000 baselines, two flavors', '1.5')]

# The span quoted for the backend as a whole, in the two documents that
# summarise it in one sentence rather than tabulating it.
SPEEDUP_SPAN = '1.5x to 20x'


def read(path):
    r"""Returns the text of ``path``."""
    with open(path) as handle:
        return handle.read()


def read_flowed(path):
    r"""Returns the text of ``path`` with runs of whitespace collapsed.

    Prose in docstrings and rst is hard-wrapped, so a phrase like "eleven
    configurations" is routinely split across a line break and a search
    for it as written finds nothing.  Collapsing whitespace first is the
    difference between a check and a check that always passes --- the
    first version of the case-count guard below matched nothing for
    exactly this reason.
    """
    return re.sub(r'\s+', ' ', read(path))


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


def test_kernel_speedups_agree_between_module_and_methodology():
    r"""Every row says the same thing in both places.

    Added with the four-flavor kernel in 1.10.0, because the figure is
    stated twice the moment it is stated at all --- which is the shape of
    every drift this module already guards against.  Widened from the two
    four-flavor rows to the whole table once the other four were found to
    have drifted apart in exactly that way.
    """
    for stack, speedup in SPEEDUP_KERNEL:
        rows = stack.split(',')[0]
        for path in (FASTKERNELS, METHODOLOGY_RST):
            text = read_flowed(path)
            assert '~%sx' % speedup in text, (
                '%s does not quote ~%sx, the measured speedup for the "%s" '
                'row of the kernel benchmark'
                % (os.path.relpath(path, ROOT), speedup, stack))
            assert rows in text, (
                '%s quotes kernel speedups but no longer names the "%s" '
                'stack they were measured on'
                % (os.path.relpath(path, ROOT), stack))


def test_the_quoted_speedup_span_covers_the_tabulated_rows():
    r"""The one-sentence summary brackets the tabulated figures.

    The README and the quickstart compress the whole table into a single
    range.  A row outside either end of that range is exactly the kind of
    contradiction that has been found by hand here before.
    """
    speedups = [float(speedup) for _, speedup in SPEEDUP_KERNEL]
    low, high = min(speedups), max(speedups)
    quoted_low, quoted_high = (float(value) for value in
                               re.search(r'([\d.]+)x to ([\d.]+)x',
                                         SPEEDUP_SPAN).groups())
    assert quoted_low <= low, (
        'the documents quote "%s" but the table reaches down to ~%gx'
        % (SPEEDUP_SPAN, low))
    assert quoted_high >= high, (
        'the documents quote "%s" but the table reaches ~%gx'
        % (SPEEDUP_SPAN, high))

    for path in (README, QUICKSTART_RST):
        text = read_flowed(path)
        assert SPEEDUP_SPAN in text, (
            '%s should quote the backend as worth "%s"; the tabulated '
            'rows span ~%gx to ~%gx' % (os.path.relpath(path, ROOT),
                                        SPEEDUP_SPAN, low, high))


def test_the_methodology_speedups_are_the_landing_page_ratios():
    r"""One set of measurements, quoted two ways, kept consistent.

    ``index.rst`` and ``README.md`` give the four scans as absolute
    timings; ``methodology.rst`` gives the same four as ratios.  They
    were written independently and had drifted --- 30x where the timings
    give 21x, and 70x where they give 99x --- so the ratios are now
    derived here from the timings rather than trusted.
    """
    text = read_flowed(METHODOLOGY_RST)

    for scan, loop, arrays, _ in BENCHMARKS:
        ratio = round(float(loop)/float(arrays))
        assert '~%dx' % ratio in text, (
            'methodology.rst should quote ~%dx for the "%s" row, which is '
            'what the %s ms and %s ms in the benchmark table give'
            % (ratio, scan, loop, arrays))


def test_the_readme_states_its_own_benchmark_ratios_correctly():
    r"""The README divides its own table in the margin, so check it.

    ``README.md`` annotates each row of the benchmark table with the
    ratio of the two timings on that same line --- "0.07 ms (~99x)".
    Nothing checked those, and one of them was wrong: the two-flavor row
    read ~93x where 6.9/0.07 is 99, contradicting both the numbers beside
    it and the ratio quoted on the methodology page.  Three of the four
    were right, which is how it survived three passes by eye.
    """
    text = read_flowed(README)

    for scan, loop, arrays, numba in BENCHMARKS:
        ratio = round(float(loop)/float(arrays))
        assert '(~%d×)' % ratio in text, (
            'README.md should annotate the "%s" row with (~%d×), which '
            'is what the %s ms and %s ms on that same line give'
            % (scan, ratio, loop, arrays))

        if numba is not None:
            with_numba = round(float(loop)/float(numba), -1)
            assert '~%d×' % with_numba in text, (
                'README.md should annotate the "%s" row with ~%d× for '
                'the compiled path, which is what the %s ms and %s ms on '
                'that same line give' % (scan, with_numba, loop, numba))


def test_two_flavor_row_is_marked_as_not_using_the_backend():
    r"""The blank cell is explained, not left as an apparent omission."""
    for path in (README, INDEX_RST):
        text = read(path)
        assert 'not used' in text, (
            '%s should say the compiled backend is "not used" for the '
            'two-flavor row rather than leaving it blank'
            % os.path.relpath(path, ROOT))


# ---------------------------------------------------------------------------
# The cross-check numbers, which are described in prose in several places
# ---------------------------------------------------------------------------

NUSQUIDS_JSON = os.path.join(ROOT, 'tests', 'nusquids_reference.json')
MAKE_NOTEBOOKS = os.path.join(ROOT, 'notebooks', 'make_notebooks.py')
NUSQUIDS_TEST = os.path.join(ROOT, 'tests', 'test_nusquids_comparison.py')
EIGENVALUE_TEST = os.path.join(ROOT, 'tests', 'test_matter_eigenvalues.py')

COUNT_WORDS = {7: 'seven', 8: 'eight', 9: 'nine', 10: 'ten', 11: 'eleven',
               12: 'twelve'}

# Files whose prose says how many configurations are frozen.
COUNT_FILES = [MAKE_NOTEBOOKS, NUSQUIDS_TEST, EIGENVALUE_TEST]


def frozen_cases():
    r"""Returns the frozen nuSQuIDS reference cases."""
    with open(NUSQUIDS_JSON) as handle:
        return json.load(handle)['cases']


def measured_agreement():
    r"""Returns ``(best, worst)`` agreement over every frozen case.

    Recomputed here rather than read from anywhere, so that the range the
    documents quote is checked against what the comparison actually
    produces today.
    """
    import test_nusquids_comparison as comparison

    deviations = []
    for case in comparison.CASES:
        density = case['density_g_cm3']
        potentials = ((None, None) if density is None
                      else comparison.matched_potentials(density))
        ours = comparison.our_probabilities(
            case, comparison.matched_km(), potentials)
        theirs = np.asarray(case['probabilities'])
        deviations.append(np.max(np.abs(ours - theirs)))

    return min(deviations), max(deviations)


def test_the_frozen_case_count_is_described_correctly():
    r"""Prose saying how many cases are frozen agrees with the data.

    This is the claim that has gone stale most often: the set grew from
    seven to eleven and three separate files went on saying seven.
    """
    count = len(frozen_cases())
    correct = COUNT_WORDS[count]

    for path in COUNT_FILES:
        text = read_flowed(path)
        for number, word in COUNT_WORDS.items():
            if number == count:
                continue
            stale = re.search(r'\b%s configurations\b' % word, text)
            assert not stale, (
                '%s says "%s configurations"; there are %d, so it should '
                'say "%s"' % (os.path.relpath(path, ROOT), word, count,
                              correct))


def test_the_quoted_agreement_range_matches_the_measured_one():
    r"""The range the documents quote brackets what is measured.

    Checked as decades rather than digits: the point is that a document
    claiming agreement "to 1e-15" when the worst case is 3e-10 is wrong
    by five orders of magnitude, not that the mantissa has drifted.
    """
    best, worst = measured_agreement()

    assert best < 1.0e-14, (
        'the best frozen case now agrees only to %.2e; the documents '
        'describe round-off agreement' % best)
    assert worst < 1.0e-9, (
        'the worst frozen case now agrees only to %.2e, which is outside '
        'what oscprob4nu.POLISH_ROOTS documents' % worst)

    # The superseded upper bound, written before the antineutrino case
    # existed.  It appeared in two files and was wrong in both.
    for path in (MAKE_NOTEBOOKS, NUSQUIDS_TEST, INDEX_RST):
        text = read_flowed(path)
        assert 'between 1e-11 and 1e-15' not in text, (
            '%s quotes a superseded agreement range; the worst case is '
            'now %.2e' % (os.path.relpath(path, ROOT), worst))


def test_every_frozen_case_is_reachable_by_the_notebook():
    r"""The notebook does not silently skip a case it should cover.

    ``ours()`` in the notebook handles antineutrinos and a mass-ordering
    override.  A case added to the JSON with some further variation would
    be computed with the wrong Hamiltonian and quietly widen the range,
    so the keys the notebook understands are pinned here.
    """
    understood = {'name', 'n_flavors', 'energy_gev', 'baseline_km',
                  'density_g_cm3', 'antineutrino', 'dm31', 'probabilities'}

    for case in frozen_cases():
        unexpected = set(case) - understood
        assert not unexpected, (
            'the frozen case %r carries %s, which neither the notebook nor '
            'the comparison test knows how to apply'
            % (case['name'], sorted(unexpected)))


# ---------------------------------------------------------------------------
# Code the documentation shows but nothing runs
# ---------------------------------------------------------------------------

QUICKSTART_RST_PATH = QUICKSTART_RST
INDEX_RST_PATH = INDEX_RST


def _code_blocks(path):
    r"""Returns the ``.. code-block:: python`` bodies in an rst file.

    These are *not* the ``jupyter-execute`` blocks, which the
    documentation build runs and which therefore cannot rot silently.
    A ``code-block`` is inert: Sphinx renders it and never executes it.
    """
    text = read(path)
    blocks = []
    for match in re.finditer(r'\.\. code-block:: python\n\n((?:   +.*\n|\n)+)',
                             text):
        blocks.append('\n'.join(
            line[3:] if line.startswith('   ') else line
            for line in match.group(1).split('\n')))

    return blocks


def test_the_landing_page_getting_started_example_runs():
    r"""The first code a reader sees is executed by nothing else.

    ``index.rst`` shows a complete worked example under *Getting
    started* as a ``code-block``, which Sphinx renders without running.
    Every other snippet in the documentation is either a
    ``jupyter-execute`` block, run at build time, or a deliberate
    fragment; this one is a whole program and is checked here.
    """
    import subprocess
    import sys

    complete = [b for b in _code_blocks(INDEX_RST_PATH)
                if 'probabilities_3nu' in b]
    assert complete, 'index.rst no longer shows a getting-started example'

    for block in complete:
        result = subprocess.run(
            [sys.executable, '-c',
             'import sys; sys.path.insert(0, %r)\n' % os.path.join(ROOT, 'src')
             + block],
            capture_output=True, text=True)
        assert result.returncode == 0, (
            'the getting-started example in index.rst does not run:\n%s'
            % result.stderr[-800:])


def test_the_documented_switches_are_spelled_correctly():
    r"""``CHECK_HERMITICITY`` and ``USE_NUMBA`` snippets actually work.

    Both are shown as two-line ``code-block`` snippets that set a module
    attribute.  A typo in either --- a renamed switch, a wrong module ---
    renders perfectly and does nothing, which is the worst way for a
    documented escape hatch to fail.
    """
    import importlib

    for block in _code_blocks(QUICKSTART_RST_PATH):
        found = re.search(r'import (\w+)\n\s*(\w+)\.(\w+) = (\w+)', block)
        if not found:
            continue
        module_name, target, attribute, value = found.groups()
        assert module_name == target, block
        module = importlib.import_module(module_name)
        assert hasattr(module, attribute), (
            'quickstart.rst sets %s.%s, which does not exist'
            % (module_name, attribute))
        assert isinstance(getattr(module, attribute), bool), (
            '%s.%s is not a switch' % (module_name, attribute))
