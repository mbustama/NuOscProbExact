# -*- coding: utf-8 -*-
r"""Runs the examples embedded in the module docstrings as doctests.

Every ``Example`` block in the library is executable, so that the
numbers quoted in the documentation --- which Sphinx renders verbatim
--- cannot drift away from what the code actually returns.
"""

import doctest

import pytest

import fastkernels
import globaldefs
import hamiltonians2nu
import hamiltonians3nu
import oscprob2nu
import oscprob3nu


MODULES = [globaldefs, oscprob2nu, oscprob3nu, hamiltonians2nu,
           hamiltonians3nu, fastkernels]


@pytest.mark.parametrize('module', MODULES, ids=[m.__name__ for m in MODULES])
def test_docstring_examples(module):
    r"""Every doctest in the module runs and produces the documented
    output."""
    results = doctest.testmod(module, verbose=False, report=True,
                              optionflags=doctest.NORMALIZE_WHITESPACE)
    assert results.failed == 0, \
        '%d of %d doctests failed in %s' % (results.failed, results.attempted,
                                            module.__name__)


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
