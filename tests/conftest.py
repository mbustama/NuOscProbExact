# -*- coding: utf-8 -*-
r"""Shared pytest configuration and fixtures for the NuOscProbExact tests.

Puts ``src/`` on ``sys.path`` so that the core modules can be imported
without installing the package, and provides the random Hermitian
matrices and reference Hamiltonians used across several test modules.
"""

import os
import sys

import numpy as np
import pytest

SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   'src')
if SRC not in sys.path:
    sys.path.insert(0, SRC)


# Absolute tolerance used whenever a quantity is expected to agree with an
# independent computation to within floating-point round-off.
ATOL = 1.0e-12


def random_hermitian(rng, n):
    r"""Returns a random ``n``-by-``n`` complex Hermitian matrix."""
    a = rng.normal(size=(n, n)) + 1.0j*rng.normal(size=(n, n))
    return (a + a.conj().T)/2.0


def traceless(h):
    r"""Returns the traceless part of the square matrix ``h``."""
    h = np.asarray(h, dtype=complex)
    n = h.shape[0]
    return h - np.trace(h)/n*np.eye(n)


def as_nested_list(h):
    r"""Returns the matrix ``h`` as the nested list the modules expect."""
    return [list(row) for row in np.asarray(h)]


@pytest.fixture(scope='session')
def rng():
    r"""Returns a seeded random-number generator, for reproducibility."""
    return np.random.default_rng(20190420)


@pytest.fixture(scope='session')
def hermitian_2nu(rng):
    r"""Returns 100 random 2x2 Hermitian matrices with baselines."""
    return [(random_hermitian(rng, 2), rng.uniform(0.1, 5.0))
            for _ in range(100)]


@pytest.fixture(scope='session')
def hermitian_3nu(rng):
    r"""Returns 100 random 3x3 Hermitian matrices with baselines."""
    return [(random_hermitian(rng, 3), rng.uniform(0.1, 5.0))
            for _ in range(100)]


@pytest.fixture(params=['numpy', 'numba'])
def backend(request, monkeypatch):
    r"""Runs the test once per available backend.

    The ``numba`` case is skipped when Numba is absent; the ``numpy``
    case always runs, with the compiled path forced off.

    This lives here rather than beside the tests of :mod:`fastkernels`
    because the four-flavor tests need it too: with a compiled kernel in
    play, every batched assertion in :mod:`oscprob4nu`'s own suite is
    otherwise made about whichever backend happens to be installed.
    """
    import fastkernels

    if request.param == 'numba':
        if not fastkernels.HAVE_NUMBA:
            pytest.skip('Numba is not installed')
        monkeypatch.setattr(fastkernels, 'USE_NUMBA', True)
    else:
        monkeypatch.setattr(fastkernels, 'USE_NUMBA', False)
    return request.param


@pytest.fixture
def kernel_spy(monkeypatch):
    r"""Counts calls to the compiled kernels, by name.

    A test that means to exercise the Numba backend can assert that it
    actually did.  Without that, raising a dispatch threshold silently
    turns such a test into a comparison of the NumPy path with itself,
    which is how it passes for the wrong reason.
    """
    import fastkernels

    # The dict handed back is the one the counters write to, not a copy
    counts = _CountingDict()

    for name in ('probabilities_2nu_kernel', 'probabilities_3nu_kernel',
                 'probabilities_4nu_kernel',
                 'evolution_operator_3nu_kernel',
                 'slab_product_3nu_kernel'):
        original = getattr(fastkernels, name, None)
        if original is None:
            continue

        def counted(*args, __name=name, __original=original, **kwargs):
            counts[__name] = counts[__name] + 1
            return __original(*args, **kwargs)

        monkeypatch.setattr(fastkernels, name, counted)

    return counts


class _CountingDict(dict):
    r"""A dict that reports zero for names never called."""

    def __missing__(self, key):
        return 0
