# -*- coding: utf-8 -*-
r"""Tests that the type annotations exist and agree with the docstrings.

Every routine carries annotations on its parameters and its return
value.  Sphinx renders both the annotation and the docstring type, so
the two disagreeing is worse than either being absent: the reader gets
two different answers to the same question.
"""

import inspect
import re
import typing

import pytest

import earth
import fastkernels
import globaldefs
import hamiltonians2nu
import hamiltonians3nu
import hamiltonians4nu
import oscprob2nu
import oscprob3nu
import oscprob4nu
import slabs


MODULES = [oscprob2nu, oscprob3nu, oscprob4nu, hamiltonians2nu,
           hamiltonians3nu, hamiltonians4nu, fastkernels, slabs, earth]


def routines(module, private=False):
    r"""Returns the routines defined in `module`."""
    return [(name, fn) for name, fn in inspect.getmembers(module,
                                                          inspect.isfunction)
            if fn.__module__ == module.__name__
            and (private or not name.startswith('_'))]


def documented_types(func):
    r"""Returns the ``name -> type`` map from the Parameters section."""
    doc = inspect.getdoc(func) or ''
    match = re.search(r'Parameters\n-+\n(.*?)(\n\n[A-Z]|\Z)', doc, re.S)
    if not match:
        return {}
    return dict(re.findall(r'^(\w+) : (.+)$', match.group(1), re.M))


@pytest.mark.parametrize('module', MODULES, ids=[m.__name__ for m in MODULES])
def test_every_parameter_is_annotated(module):
    r"""No parameter of any routine is left unannotated."""
    missing = []
    for name, fn in routines(module, private=True):
        hints = typing.get_type_hints(fn)
        for pname in inspect.signature(fn).parameters:
            if pname not in hints:
                missing.append('%s(%s)' % (name, pname))
    assert not missing, 'unannotated parameters: %s' % ', '.join(missing)


@pytest.mark.parametrize('module', MODULES, ids=[m.__name__ for m in MODULES])
def test_every_routine_annotates_its_return(module):
    r"""Every routine declares what it returns."""
    missing = [name for name, fn in routines(module, private=True)
               if 'return' not in typing.get_type_hints(fn)]
    assert not missing, 'no return annotation: %s' % ', '.join(missing)


@pytest.mark.parametrize('module', MODULES, ids=[m.__name__ for m in MODULES])
def test_every_parameter_is_documented(module):
    r"""Every annotated parameter of a public routine is documented."""
    missing = []
    for name, fn in routines(module):
        documented = documented_types(fn)
        for pname in inspect.signature(fn).parameters:
            if pname not in documented:
                missing.append('%s(%s)' % (name, pname))
    assert not missing, 'undocumented parameters: %s' % ', '.join(missing)


@pytest.mark.parametrize('module', MODULES, ids=[m.__name__ for m in MODULES])
def test_docstring_types_agree_with_the_annotations(module):
    r"""A parameter whose annotation accepts an array says so in its
    docstring, and vice versa.

    This is the mismatch that actually misleads: an annotation that
    admits ``np.ndarray`` while the prose promises a scalar, or the
    reverse.
    """
    problems = []
    for name, fn in routines(module):
        hints = typing.get_type_hints(fn)
        documented = documented_types(fn)
        for pname, annotation in hints.items():
            if pname == 'return' or pname not in documented:
                continue
            text = str(annotation)
            prose = documented[pname]
            takes_array = 'ndarray' in text or 'list' in text
            says_array = 'array' in prose or 'list' in prose
            if takes_array != says_array:
                problems.append(
                    '%s(%s): annotation %r vs docstring %r'
                    % (name, pname, text, prose))
    assert not problems, '; '.join(problems)


@pytest.mark.parametrize('module', MODULES, ids=[m.__name__ for m in MODULES])
def test_annotations_resolve(module):
    r"""Every annotation is a real, importable type.

    ``typing.get_type_hints`` raises if an annotation names something
    that does not exist, which catches a typo that a bare string
    annotation would otherwise hide until a reader tripped over it.
    """
    for name, fn in routines(module, private=True):
        typing.get_type_hints(fn)


def test_globaldefs_defines_no_routines():
    r"""globaldefs holds constants only, so it has nothing to annotate."""
    assert routines(globaldefs, private=True) == []
