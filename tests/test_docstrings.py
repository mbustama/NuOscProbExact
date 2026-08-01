# -*- coding: utf-8 -*-
r"""Runs the examples embedded in the module docstrings.

Every ``Examples`` block in the library is a ``.. jupyter-execute::``
directive, which Sphinx runs while building the documentation, so the
output on the API page is produced by the code rather than pasted beside
it and cannot drift.

That leaves one gap this file closes.  The documentation is built by a
single job on a single interpreter; the library is supported on five.
An example that works on 3.12 and not on 3.9 would reach the published
page unnoticed, and a broken example would be found only at documentation
build time --- late, and in a job whose failure reads as a documentation
problem rather than a code one.  So the blocks are extracted and executed
here as well, on every supported Python.

They are executed, not compared against an expected value: there is no
longer an expected value to compare with, which is the point of having
Sphinx compute it.  What is checked is that every example still runs.
"""

import doctest

import pytest

import earth
import fastkernels
import globaldefs
import hamiltonians2nu
import hamiltonians3nu
import oscprob2nu
import oscprob3nu
import slabs


MODULES = [globaldefs, oscprob2nu, oscprob3nu, hamiltonians2nu,
           hamiltonians3nu, fastkernels, slabs, earth]


def jupyter_execute_blocks(docstring):
    r"""Returns the source of each ``.. jupyter-execute::`` block.

    Parameters
    ----------
    docstring : str
        A docstring, already dedented by `inspect.getdoc`.

    Returns
    -------
    list of str
        One entry per block, ready to be executed.
    """
    blocks, lines, i = [], docstring.split('\n'), 0
    while i < len(lines):
        if lines[i].strip() == '.. jupyter-execute::':
            i += 1
            body = []
            while i < len(lines) and (not lines[i].strip()
                                      or lines[i].startswith('    ')):
                body.append(lines[i][4:] if lines[i].startswith('    ')
                            else '')
                i += 1
            blocks.append('\n'.join(body).strip('\n'))
        else:
            i += 1
    return blocks


@pytest.mark.parametrize('module', MODULES, ids=[m.__name__ for m in MODULES])
def test_docstring_examples_run(module):
    r"""Every example in the module executes without raising."""
    import inspect

    failures = []
    for name, routine in inspect.getmembers(module, inspect.isfunction):
        if routine.__module__ != module.__name__:
            continue
        for source in jupyter_execute_blocks(inspect.getdoc(routine) or ''):
            try:
                exec(compile(source, '<%s.%s>' % (module.__name__, name),
                             'exec'), {})
            except Exception as error:                # noqa: BLE001
                failures.append('%s: %s: %s'
                                % (name, type(error).__name__, error))

    assert not failures, 'examples failed in %s:\n  %s' % (
        module.__name__, '\n  '.join(failures))


@pytest.mark.parametrize('module', MODULES, ids=[m.__name__ for m in MODULES])
def test_no_doctests_are_left_behind(module):
    r"""No ``>>>`` example survives in the docstrings.

    A stray doctest would render as a prompt with a pasted result on the
    API page --- the thing converting these blocks was meant to remove ---
    and would not be run by the test above, which only knows about
    ``jupyter-execute`` blocks.
    """
    results = doctest.testmod(module, verbose=False, report=False)
    assert results.attempted == 0, \
        '%s still has %d doctest example(s); convert them to ' \
        'jupyter-execute' % (module.__name__, results.attempted)


@pytest.mark.parametrize('module', MODULES, ids=[m.__name__ for m in MODULES])
def test_every_public_routine_has_a_docstring(module):
    r"""Every public routine carries a docstring, as Sphinx requires."""
    import inspect
    missing = [name for name, obj in inspect.getmembers(module,
                                                        inspect.isfunction)
               if obj.__module__ == module.__name__
               and not name.startswith('_')
               and not (obj.__doc__ or '').strip()]
    assert not missing, 'missing docstrings: %s' % ', '.join(missing)


@pytest.mark.parametrize('module', MODULES, ids=[m.__name__ for m in MODULES])
def test_docstrings_have_no_unresolved_placeholders(module):
    r"""No placeholder text such as 'arXiv:1904.XXXXX' survives."""
    import inspect
    offenders = []
    for name, obj in [(module.__name__, module)] + \
            inspect.getmembers(module, inspect.isfunction):
        doc = inspect.getdoc(obj) or ''
        if 'XXXX' in doc or 'TODO' in doc:
            offenders.append(name)
    assert not offenders, 'placeholder text in: %s' % ', '.join(offenders)
