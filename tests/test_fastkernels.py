# -*- coding: utf-8 -*-
r"""Tests of the optional Numba backend and of the fallback path.

:mod:`fastkernels` compiles the batched expansions when Numba is
installed.  The point of these tests is that the two paths are
interchangeable: whichever is in use, the answers must be the same to
round-off, and turning the compiled one off must change nothing but the
speed.

Every equivalence test here runs with ``fastkernels.USE_NUMBA`` forced
off as well as on, so the NumPy path is exercised even on a machine
where Numba is installed --- otherwise it would rot unnoticed.
"""

import numpy as np
import pytest

import fastkernels
import oscprob2nu
import oscprob3nu

from conftest import ATOL, as_nested_list, random_hermitian


needs_numba = pytest.mark.skipif(not fastkernels.HAVE_NUMBA,
                                 reason='Numba is not installed')


@pytest.fixture(params=['numpy', 'numba'])
def backend(request, monkeypatch):
    r"""Runs the test once per available backend.

    The ``numba`` case is skipped when Numba is absent; the ``numpy``
    case always runs, with the compiled path forced off.
    """
    if request.param == 'numba':
        if not fastkernels.HAVE_NUMBA:
            pytest.skip('Numba is not installed')
        monkeypatch.setattr(fastkernels, 'USE_NUMBA', True)
    else:
        monkeypatch.setattr(fastkernels, 'USE_NUMBA', False)
    return request.param


# --------------------------------------------------------------------------
# The switch itself
# --------------------------------------------------------------------------

def test_available_follows_the_switch(monkeypatch):
    r"""available() is the conjunction of the import and the switch."""
    monkeypatch.setattr(fastkernels, 'USE_NUMBA', False)
    assert fastkernels.available() is False
    monkeypatch.setattr(fastkernels, 'USE_NUMBA', True)
    assert fastkernels.available() is fastkernels.HAVE_NUMBA


def test_library_imports_without_numba(monkeypatch):
    r"""With Numba reported absent, the library still works.

    This is the state of an installation without the optional extra.
    """
    monkeypatch.setattr(fastkernels, 'HAVE_NUMBA', False)
    assert fastkernels.available() is False

    h_stack = np.stack([np.diag([1.0, 2.0, -3.0]).astype(complex)]*40)
    prob = oscprob3nu.probabilities_3nu(h_stack, 1.5)
    assert prob.shape == (40, 9)
    assert np.allclose(prob.reshape(-1, 3, 3).sum(axis=-1), 1.0, atol=ATOL)


# --------------------------------------------------------------------------
# The two paths agree
# --------------------------------------------------------------------------

@needs_numba
@pytest.mark.parametrize('n', [1, 5, 17, 64, 300, 1000])
def test_3nu_kernel_matches_numpy_path(rng, monkeypatch, n):
    r"""Compiled and NumPy paths agree for stacks of every size.

    The sizes straddle both thresholds: `oscprob3nu.SMALL_BATCH`, below
    which the NumPy path defers to the scalar one, and
    `fastkernels.PARALLEL_THRESHOLD`, above which the kernel spreads
    over cores.
    """
    h_stack = np.stack([random_hermitian(rng, 3) for _ in range(n)])
    l_stack = rng.uniform(0.1, 50.0, n)

    monkeypatch.setattr(fastkernels, 'USE_NUMBA', True)
    with_numba = oscprob3nu.probabilities_3nu(h_stack, l_stack)
    monkeypatch.setattr(fastkernels, 'USE_NUMBA', False)
    without = oscprob3nu.probabilities_3nu(h_stack, l_stack)

    assert with_numba.shape == without.shape == (n, 9)
    assert np.allclose(with_numba, without, atol=1.0e-12)


@needs_numba
@pytest.mark.parametrize('n', [1, 17, 300])
def test_2nu_kernel_matches_numpy_path(rng, monkeypatch, n):
    r"""The same, for two flavors."""
    h_stack = np.stack([random_hermitian(rng, 2) for _ in range(n)])
    l_stack = rng.uniform(0.1, 50.0, n)

    monkeypatch.setattr(fastkernels, 'USE_NUMBA', True)
    with_numba = oscprob2nu.probabilities_2nu(h_stack, l_stack)
    monkeypatch.setattr(fastkernels, 'USE_NUMBA', False)
    without = oscprob2nu.probabilities_2nu(h_stack, l_stack)

    assert np.allclose(with_numba, without, atol=ATOL)


@needs_numba
def test_kernel_matches_the_scalar_path(rng, monkeypatch):
    r"""The compiled kernel agrees with the scalar routine element by
    element, which is the reference both batched paths answer to."""
    monkeypatch.setattr(fastkernels, 'USE_NUMBA', True)
    h_stack = np.stack([random_hermitian(rng, 3) for _ in range(40)])
    l_stack = rng.uniform(0.1, 50.0, 40)

    batched = oscprob3nu.probabilities_3nu(h_stack, l_stack)
    scalar = np.array([oscprob3nu.probabilities_3nu(as_nested_list(h), l)
                       for h, l in zip(h_stack, l_stack)])
    assert np.allclose(batched, scalar, atol=1.0e-12)


# --------------------------------------------------------------------------
# Everything the NumPy path guarantees, the kernel must too
# --------------------------------------------------------------------------

def test_broadcasting_patterns(backend, rng):
    r"""All three broadcasting patterns work on whichever backend."""
    h_stack = np.stack([random_hermitian(rng, 3) for _ in range(7)])
    l_stack = np.linspace(0.5, 4.0, 5)

    assert oscprob3nu.probabilities_3nu(h_stack, 2.0).shape == (7, 9)
    assert oscprob3nu.probabilities_3nu(
        as_nested_list(h_stack[0]), l_stack).shape == (5, 9)
    grid = oscprob3nu.probabilities_3nu(h_stack[:, None, :, :],
                                        l_stack[None, :])
    assert grid.shape == (7, 5, 9)

    for i in range(7):
        for j in range(5):
            expected = oscprob3nu.probabilities_3nu(
                as_nested_list(h_stack[i]), l_stack[j])
            assert np.allclose(grid[i, j], expected, atol=1.0e-12)


DEGENERATE = [np.eye(3, dtype=complex),
              np.zeros((3, 3), dtype=complex),
              np.diag([1.0, 1.0, -2.0]).astype(complex),
              np.diag([2.0, -1.0, -1.0]).astype(complex),
              1.0e6*np.eye(3, dtype=complex)]


def test_degenerate_hamiltonians(backend, rng):
    r"""Degenerate spectra survive whichever backend evaluates them.

    The compiled kernel carries its own copy of the degeneracy branch,
    so this is not covered by the NumPy path's tests.
    """
    ordinary = [random_hermitian(rng, 3) for _ in range(20)]
    h_stack = np.stack(DEGENERATE + ordinary)
    l_stack = np.full(len(h_stack), 3.0)

    prob = oscprob3nu.probabilities_3nu(h_stack, l_stack)
    assert np.all(np.isfinite(prob))
    assert np.allclose(prob.reshape(-1, 3, 3).sum(axis=-1), 1.0, atol=1.0e-12)

    scalar = np.array([oscprob3nu.probabilities_3nu(as_nested_list(h), 3.0)
                       for h in h_stack])
    assert np.allclose(prob, scalar, atol=1.0e-12)


def test_probabilities_are_normalized(backend, rng):
    r"""Normalization and positivity hold on whichever backend."""
    h_stack = np.stack([random_hermitian(rng, 3) for _ in range(400)])
    l_stack = rng.uniform(0.1, 500.0, 400)

    prob = oscprob3nu.probabilities_3nu(h_stack, l_stack).reshape(400, 3, 3)
    assert np.allclose(prob.sum(axis=-1), 1.0, atol=1.0e-12)
    assert np.allclose(prob.sum(axis=-2), 1.0, atol=1.0e-12)
    assert np.all(prob >= -ATOL) and np.all(prob <= 1.0+ATOL)


def test_empty_and_single_element_stacks(backend):
    r"""Degenerate stack sizes behave on whichever backend."""
    empty = oscprob3nu.probabilities_3nu(np.zeros((0, 3, 3), dtype=complex),
                                         np.zeros(0))
    assert empty.shape == (0, 9)

    one = oscprob3nu.probabilities_3nu(np.eye(3, dtype=complex)[None], 1.0)
    assert one.shape == (1, 9)
    assert np.allclose(one.reshape(3, 3), np.eye(3), atol=ATOL)


def test_zero_baseline(backend, rng):
    r"""L = 0 gives the identity on whichever backend."""
    h_stack = np.stack([random_hermitian(rng, 3) for _ in range(30)])
    prob = oscprob3nu.probabilities_3nu(h_stack, np.zeros(30))
    assert np.allclose(prob.reshape(30, 3, 3), np.eye(3), atol=ATOL)


def test_small_batch_matches_large_batch(backend, rng):
    r"""The small-stack shortcut returns what the array path would.

    `oscprob3nu.SMALL_BATCH` sends short stacks through the scalar
    routine; padding the same Hamiltonians out past the threshold must
    not change the answers for the ones in common.
    """
    h_stack = np.stack([random_hermitian(rng, 3) for _ in range(64)])
    l_stack = rng.uniform(0.1, 20.0, 64)

    short = oscprob3nu.probabilities_3nu(h_stack[:5], l_stack[:5])
    long = oscprob3nu.probabilities_3nu(h_stack, l_stack)
    assert np.allclose(short, long[:5], atol=1.0e-12)


@needs_numba
def test_kernels_are_cached_on_disk():
    r"""The kernels are declared with cache=True.

    Without it every interpreter start pays seconds of compilation,
    which would make the optional backend a poor trade.
    """
    from numba.core.caching import NullCache
    for name in ('_run_3nu_serial', '_run_3nu_parallel',
                 '_run_2nu_serial', '_run_2nu_parallel'):
        dispatcher = getattr(fastkernels, name)
        assert not isinstance(dispatcher._cache, NullCache), name
