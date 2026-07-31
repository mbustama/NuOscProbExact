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
