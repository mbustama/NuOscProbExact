# -*- coding: utf-8 -*-
r"""Checks the paper against the repository it describes.

``resources/paper/main.tex`` restates things the repository already
knows: how many notebooks there are, what each one is called, and how
many tests ship.  Every such restatement is a copy, and copies drift ---
this repository has already carried "eight PDFs" when fourteen were
written, a test table summing to one total while stating another, and a
notebook count that stayed at eighteen through two additions.  Those
were all found by reading, one release later.

The paper is not built by continuous integration and no reviewer reruns
it, so a stale table there survives longer than a stale sentence
anywhere else.  These tests derive the numbers and the names from the
repository and fail the build when the paper falls behind.

The blurbs in the paper's table are *not* checked: they are written for
a journal reader, deliberately plainer than the ones in
``docs/source/notebooks.rst``, and pinning them here would make two
voices into one.  What is checked is everything a reader could use to
find a notebook --- its number, its position in the reading order, and
its title.
"""

import ast
import glob
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

PAPER = os.path.join(ROOT, 'resources', 'paper', 'main.tex')
MAKE_NOTEBOOKS = os.path.join(ROOT, 'notebooks', 'make_notebooks.py')

# Words for counts the notebook collection plausibly reaches, as
# `test_documented_figures.py` holds them for the other documents.
NOTEBOOK_COUNT_WORDS = {
    15: 'fifteen', 16: 'sixteen', 17: 'seventeen', 18: 'eighteen',
    19: 'nineteen', 20: 'twenty', 21: 'twenty-one', 22: 'twenty-two',
    23: 'twenty-three', 24: 'twenty-four', 25: 'twenty-five',
}


def read(path):
    r"""Returns the text of ``path``."""
    with open(path) as handle:
        return handle.read()


def reading_order():
    r"""Returns ``READING_ORDER`` from the notebook generator.

    That list is the one place the notebooks are named, titled and
    ordered; the generator asserts it against the notebooks it builds,
    and ``test_file_tree.py`` checks the documentation page against it.
    The paper's table is a third reader of the same list, so it is
    checked against the list rather than against the page.
    """
    block = re.search(r'READING_ORDER = (\[.*?\n\])\n', read(MAKE_NOTEBOOKS),
                      re.S)
    assert block, 'READING_ORDER not found in notebooks/make_notebooks.py'

    return ast.literal_eval(block.group(1))


def tabular_labelled(label):
    r"""Returns the body of the tabular whose caption carries ``label``.

    Found by looking *backwards* from the label rather than forwards
    from a ``\begin{tabular}``: the caption follows the tabular it
    belongs to, and a forward search --- the first thing written here
    --- matched the first tabular in the whole document and then ran on
    to the label, so it read the wrong table and reported it as empty.
    """
    text = read(PAPER)
    caption = text.find(r'\label{%s}' % label)
    assert caption > 0, 'no table is labelled %s' % label

    end = text.rfind(r'\end{tabular}', 0, caption)
    start = text.rfind(r'\begin{tabular}', 0, end)
    assert 0 < start < end, 'no tabular precedes the caption of %s' % label

    return text[start:end]


def notebook_table_rows():
    r"""Returns ``(number, title)`` for each row of the paper's table.

    The rows are the numbered ones of the ``tab:notebooks`` tabular; the
    group headings between them are ``\multicolumn`` rows and carry no
    number, so this skips them without having to know what they say.
    """
    rows = []
    for line in tabular_labelled('tab:notebooks').split(r'\\'):
        cells = [cell.strip() for cell in line.split('&')]
        if len(cells) == 3 and cells[0].isdigit():
            rows.append((int(cells[0]), cells[1]))

    return rows


def latex_to_plain(title):
    r"""Strips the LaTeX a title picks up when it is set in the paper.

    ``SU(n)`` is written ``SU($n$)`` in the paper and plainly in the
    generator, and a title with a tilde or a monospaced word would
    differ the same way.  Comparing the stripped forms keeps the two
    files free to spell mathematics as each medium wants while still
    failing on a title that has actually changed.
    """
    title = re.sub(r'\\(?:tt|emph|textit|mathrm)\s*\{([^{}]*)\}', r'\1', title)

    return title.replace('$', '').replace('~', ' ').strip()


def test_the_paper_lists_every_notebook_in_reading_order():
    r"""The paper's table is the notebooks, all of them, in order.

    A notebook added at the end is the easy case and this catches it.
    The case worth having a test for is one inserted in the middle,
    which renumbers everything after it: the table would still have the
    right length and the right titles, and every number below the
    insertion would point at the wrong notebook.
    """
    order = reading_order()
    rows = notebook_table_rows()

    assert len(rows) == len(order), (
        'the paper lists %d notebooks; %d ship' % (len(rows), len(order)))

    for (number, title), (name, wanted, _) in zip(rows, order):
        expected = int(name.split('_')[0])
        assert number == expected, (
            'the paper numbers %s as %d; it is notebook %d'
            % (name, number, expected))
        assert latex_to_plain(title) == wanted, (
            'the paper calls notebook %d "%s"; it is titled "%s"'
            % (number, latex_to_plain(title), wanted))


def test_the_paper_lists_notebooks_that_exist():
    r"""Nothing in the table is a notebook that was renamed away.

    ``READING_ORDER`` is checked against the built notebooks by the
    generator, but only when the generator runs.  This closes the loop
    from the paper to the directory without going through it.
    """
    shipped = sorted(os.path.basename(path) for path
                     in glob.glob(os.path.join(ROOT, 'notebooks', '*.ipynb')))

    assert len(notebook_table_rows()) == len(shipped), (
        'the paper lists %d notebooks; notebooks/ holds %d'
        % (len(notebook_table_rows()), len(shipped)))

    for number, _ in notebook_table_rows():
        assert any(name.startswith('%02d_' % number) for name in shipped), (
            'the paper lists notebook %d, which is not in notebooks/'
            % number)


def test_the_paper_states_the_notebook_count_correctly():
    r"""Prose in the paper saying how many notebooks there are is right.

    The count appears in the body and again in the appendix, and the
    table is no help to either: a reader meets the sentence first.  The
    equivalent guard for the README, the generator and the release notes
    is in ``test_documented_figures.py``; the paper is not covered there
    because it spells the count with a noun in the way
    ("twenty executable notebooks") that the guard's pattern does not
    admit.
    """
    count = len(glob.glob(os.path.join(ROOT, 'notebooks', '*.ipynb')))
    assert count in NOTEBOOK_COUNT_WORDS, (
        'no spelling on file for %d notebooks; extend NOTEBOOK_COUNT_WORDS'
        % count)

    text = re.sub(r'\s+', ' ', read(PAPER))
    assert re.search(r'\b%s (?:\w+ )?notebooks\b' % NOTEBOOK_COUNT_WORDS[count],
                     text, re.IGNORECASE), (
        'the paper never states the notebook count; it should say "%s"'
        % NOTEBOOK_COUNT_WORDS[count])

    for number, word in NOTEBOOK_COUNT_WORDS.items():
        if number == count:
            continue
        # One optional adjective between the count and the noun, which
        # is how the paper writes it, and bounded on both sides so that
        # "twenty" does not match inside "twenty-two".
        stale = re.search(r'\b%s (?:\w+ )?notebooks\b' % word, text,
                          re.IGNORECASE)
        assert not stale, (
            'the paper says "%s"; there are %d notebooks'
            % (stale.group(0), count))


def test_the_test_table_sums_to_the_total_it_states():
    r"""Table C.1 adds up, and covers every test file.

    It is a hand-written breakdown of a number that changes whenever a
    test is added, in a document nothing rebuilds.  Both halves have
    gone wrong before in this repository's documents: rows that no
    longer summed to the stated total, and a total that no longer
    matched the suite.  The row count is checked too, since the table is
    one row per test module.
    """
    rows, total = [], None
    for line in tabular_labelled('tab:tests').split(r'\\'):
        cells = [cell.strip() for cell in line.split('&')]
        if len(cells) != 2 or not cells[1].isdigit():
            continue
        if cells[0].replace(r'\hline', '').strip().lower() == 'total':
            total = int(cells[1])
        else:
            rows.append(int(cells[1]))

    assert total is not None, 'the test table states no total'
    assert sum(rows) == total, (
        'the rows of the test table sum to %d, against a stated total of %d'
        % (sum(rows), total))

    modules = glob.glob(os.path.join(ROOT, 'tests', 'test_*.py'))
    assert len(rows) == len(modules), (
        'the test table has %d rows; there are %d test modules'
        % (len(rows), len(modules)))
