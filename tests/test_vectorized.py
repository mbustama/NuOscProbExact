# -*- coding: utf-8 -*-
r"""Tests of the vectorised evaluation path.

`probabilities_2nu`, `probabilities_3nu` and the two evolution
operators accept a stack of Hamiltonians, a stack of baselines, or both,
and evaluate the whole stack at once.  These tests check that the fast
path agrees with the scalar one element by element --- including for the
degenerate spectra that the scalar path handles with a separate,
non-vectorisable branch --- and that broadcasting behaves as documented.
"""

import numpy as np
import pytest
from scipy.linalg import expm

import oscprob2nu
import oscprob3nu

from conftest import ATOL, as_nested_list, random_hermitian, traceless


def stack_of_hermitians(rng, n, size):
    r"""Returns ``n`` random Hermitian matrices of the given size."""
    return np.stack([random_hermitian(rng, size) for _ in range(n)])


# --------------------------------------------------------------------------
# Agreement with the scalar path
# --------------------------------------------------------------------------

def test_3nu_batch_matches_scalar_elementwise(rng):
    r"""A stack of Hamiltonians and baselines gives the scalar answer."""
    h_stack = stack_of_hermitians(rng, 40, 3)
    l_stack = rng.uniform(0.1, 5.0, 40)

    batched = oscprob3nu.probabilities_3nu(h_stack, l_stack)
    scalar = np.array([oscprob3nu.probabilities_3nu(as_nested_list(h), l)
                       for h, l in zip(h_stack, l_stack)])

    assert batched.shape == (40, 9)
    assert np.allclose(batched, scalar, atol=ATOL)


def test_2nu_batch_matches_scalar_elementwise(rng):
    r"""The same, for two flavors."""
    h_stack = stack_of_hermitians(rng, 40, 2)
    l_stack = rng.uniform(0.1, 5.0, 40)

    batched = oscprob2nu.probabilities_2nu(h_stack, l_stack)
    scalar = np.array([oscprob2nu.probabilities_2nu(as_nested_list(h), l)
                       for h, l in zip(h_stack, l_stack)])

    assert batched.shape == (40, 4)
    assert np.allclose(batched, scalar, atol=ATOL)


def test_3nu_batched_evolution_operator_matches_scalar(rng):
    r"""The batched evolution operator agrees with the scalar one."""
    h_stack = stack_of_hermitians(rng, 20, 3)
    l_stack = rng.uniform(0.1, 5.0, 20)

    batched = oscprob3nu.evolution_operator_3nu(h_stack, l_stack)
    scalar = np.array([oscprob3nu.evolution_operator_3nu(as_nested_list(h), l)
                       for h, l in zip(h_stack, l_stack)])

    assert batched.shape == (20, 3, 3)
    assert np.allclose(batched, scalar, atol=ATOL)


def test_2nu_batched_evolution_operator_matches_scalar(rng):
    r"""The same, for two flavors."""
    h_stack = stack_of_hermitians(rng, 20, 2)
    l_stack = rng.uniform(0.1, 5.0, 20)

    batched = oscprob2nu.evolution_operator_2nu(h_stack, l_stack)
    scalar = np.array([oscprob2nu.evolution_operator_2nu(as_nested_list(h), l)
                       for h, l in zip(h_stack, l_stack)])

    assert batched.shape == (20, 2, 2)
    assert np.allclose(batched, scalar, atol=ATOL)


def test_batched_evolution_operator_matches_matrix_exponential(rng):
    r"""The batched operator is still exp(-i H_0 L), to round-off."""
    h_stack = stack_of_hermitians(rng, 30, 3)
    l_stack = rng.uniform(0.1, 5.0, 30)

    batched = oscprob3nu.evolution_operator_3nu(h_stack, l_stack)
    for u, h, l in zip(batched, h_stack, l_stack):
        assert np.allclose(u, expm(-1.0j*traceless(h)*l), atol=ATOL)


# --------------------------------------------------------------------------
# Broadcasting
# --------------------------------------------------------------------------

def test_one_hamiltonian_many_baselines(rng):
    r"""A single Hamiltonian broadcasts against an array of baselines.

    This is the scan-vs-baseline case.
    """
    h = random_hermitian(rng, 3)
    l_stack = np.linspace(0.1, 20.0, 50)

    batched = oscprob3nu.probabilities_3nu(as_nested_list(h), l_stack)
    scalar = np.array([oscprob3nu.probabilities_3nu(as_nested_list(h), l)
                       for l in l_stack])

    assert batched.shape == (50, 9)
    assert np.allclose(batched, scalar, atol=ATOL)


def test_many_hamiltonians_one_baseline(rng):
    r"""An array of Hamiltonians broadcasts against a single baseline.

    This is the scan-vs-energy case, where only the Hamiltonian changes.
    """
    h_stack = stack_of_hermitians(rng, 50, 3)

    batched = oscprob3nu.probabilities_3nu(h_stack, 2.5)
    scalar = np.array([oscprob3nu.probabilities_3nu(as_nested_list(h), 2.5)
                       for h in h_stack])

    assert batched.shape == (50, 9)
    assert np.allclose(batched, scalar, atol=ATOL)


def test_two_dimensional_grid(rng):
    r"""Energy and baseline axes broadcast into a grid, as in an
    oscillogram."""
    h_stack = stack_of_hermitians(rng, 7, 3)
    l_stack = np.linspace(0.5, 4.0, 5)

    grid = oscprob3nu.probabilities_3nu(h_stack[:, None, :, :],
                                        l_stack[None, :])
    assert grid.shape == (7, 5, 9)

    for i in range(7):
        for j in range(5):
            expected = oscprob3nu.probabilities_3nu(
                as_nested_list(h_stack[i]), l_stack[j])
            assert np.allclose(grid[i, j], expected, atol=ATOL)


def test_scalar_arguments_keep_the_original_return_type(rng):
    r"""Scalar input still returns a plain tuple and nested list, so
    existing code that unpacks them is unaffected."""
    h = as_nested_list(random_hermitian(rng, 3))

    prob = oscprob3nu.probabilities_3nu(h, 1.0)
    assert isinstance(prob, tuple)
    assert len(prob) == 9
    assert all(isinstance(p, float) for p in prob)

    u = oscprob3nu.evolution_operator_3nu(h, 1.0)
    assert isinstance(u, list)
    assert isinstance(u[0], list)

    prob2 = oscprob2nu.probabilities_2nu(as_nested_list(random_hermitian(
        rng, 2)), 1.0)
    assert isinstance(prob2, tuple)
    assert len(prob2) == 4


# --------------------------------------------------------------------------
# Physical properties, on the batched path
# --------------------------------------------------------------------------

def test_batched_probabilities_are_normalized(rng):
    r"""Sum_b P(nu_a -> nu_b) = 1 for every element of the stack."""
    h_stack = stack_of_hermitians(rng, 50, 3)
    l_stack = rng.uniform(0.1, 50.0, 50)

    prob = oscprob3nu.probabilities_3nu(h_stack, l_stack).reshape(50, 3, 3)
    assert np.allclose(prob.sum(axis=-1), 1.0, atol=ATOL)
    assert np.allclose(prob.sum(axis=-2), 1.0, atol=ATOL)
    assert np.all(prob >= -ATOL)
    assert np.all(prob <= 1.0+ATOL)


def test_batched_evolution_operator_is_unitary(rng):
    r"""U^dagger U = 1 for every element of the stack."""
    h_stack = stack_of_hermitians(rng, 50, 3)
    l_stack = rng.uniform(0.1, 50.0, 50)

    u = oscprob3nu.evolution_operator_3nu(h_stack, l_stack)
    product = np.einsum('...ji,...jk->...ik', u.conj(), u)
    assert np.allclose(product, np.eye(3), atol=ATOL)


# --------------------------------------------------------------------------
# Degenerate spectra, which the vectorised path cannot branch on
# --------------------------------------------------------------------------

DEGENERATE = [
    np.eye(3, dtype=complex),
    np.zeros((3, 3), dtype=complex),
    np.diag([1.0, 1.0, -2.0]).astype(complex),
    np.diag([2.0, -1.0, -1.0]).astype(complex),
    1.0e6*np.eye(3, dtype=complex),
]


def test_degenerate_elements_inside_a_batch(rng):
    r"""Degenerate Hamiltonians mixed into a stack of ordinary ones are
    still evaluated exactly.

    The vectorised expression divides by zero for these, so they are
    recomputed individually; this checks that the fallback fires and
    lands on the right elements.
    """
    ordinary = [random_hermitian(rng, 3) for _ in range(5)]
    h_stack = np.stack(DEGENERATE + ordinary)
    l_stack = np.full(len(h_stack), 3.0)

    batched = oscprob3nu.probabilities_3nu(h_stack, l_stack)
    assert np.all(np.isfinite(batched))

    scalar = np.array([oscprob3nu.probabilities_3nu(as_nested_list(h), 3.0)
                       for h in h_stack])
    assert np.allclose(batched, scalar, atol=1.0e-10)

    prob = batched.reshape(-1, 3, 3)
    assert np.allclose(prob.sum(axis=-1), 1.0, atol=1.0e-10)


@pytest.mark.parametrize('h_matrix', DEGENERATE,
                         ids=['identity', 'zero', 'double-low', 'double-high',
                              'large-identity'])
def test_degenerate_batch_matches_matrix_exponential(h_matrix):
    r"""Each degenerate case is exp(-i H_0 L) on the batched path too."""
    l_stack = np.array([0.0, 1.0, 7.0, 1.0e3])
    u = oscprob3nu.evolution_operator_3nu(
        np.broadcast_to(h_matrix, (4, 3, 3)), l_stack)
    for operator, l in zip(u, l_stack):
        assert np.allclose(operator, expm(-1.0j*traceless(h_matrix)*l),
                           atol=1.0e-10)


def test_all_degenerate_batch(rng):
    r"""A stack containing nothing but degenerate Hamiltonians works."""
    h_stack = np.stack(DEGENERATE)
    prob = oscprob3nu.probabilities_3nu(h_stack, 2.0)
    assert np.all(np.isfinite(prob))
    assert np.allclose(prob.reshape(-1, 3, 3).sum(axis=-1), 1.0, atol=1.0e-10)


def test_2nu_degenerate_batch():
    r"""The two-flavor path handles |h| = 0 inside a stack."""
    h_stack = np.stack([np.eye(2, dtype=complex),
                        np.zeros((2, 2), dtype=complex),
                        np.array([[1.0, 0.5], [0.5, -1.0]], dtype=complex)])
    prob = oscprob2nu.probabilities_2nu(h_stack, 2.0)
    assert np.all(np.isfinite(prob))
    assert np.allclose(prob.reshape(-1, 2, 2).sum(axis=-1), 1.0, atol=ATOL)
    # The first two have no traceless part, so no flavor transitions
    assert np.allclose(prob[0], [1.0, 0.0, 0.0, 1.0], atol=ATOL)
    assert np.allclose(prob[1], [1.0, 0.0, 0.0, 1.0], atol=ATOL)


# --------------------------------------------------------------------------
# Shapes and edge cases
# --------------------------------------------------------------------------

@pytest.mark.parametrize('n', [1, 2, 17])
def test_batch_of_any_length(rng, n):
    r"""Stacks of length one and of odd lengths behave."""
    h_stack = stack_of_hermitians(rng, n, 3)
    prob = oscprob3nu.probabilities_3nu(h_stack, np.full(n, 1.5))
    assert prob.shape == (n, 9)


def test_empty_batch():
    r"""An empty stack returns an empty result rather than raising."""
    prob = oscprob3nu.probabilities_3nu(np.zeros((0, 3, 3), dtype=complex),
                                        np.zeros(0))
    assert prob.shape == (0, 9)


def test_zero_baseline_in_a_batch(rng):
    r"""L = 0 inside a stack gives the identity."""
    h_stack = stack_of_hermitians(rng, 10, 3)
    prob = oscprob3nu.probabilities_3nu(h_stack, np.zeros(10))
    assert np.allclose(prob.reshape(10, 3, 3), np.eye(3), atol=ATOL)


# --------------------------------------------------------------------------
# Broadcasting patterns for two flavors, and shape errors for both
# --------------------------------------------------------------------------

def test_2nu_one_hamiltonian_many_baselines(rng):
    r"""A single 2nu Hamiltonian broadcasts against an array of
    baselines.

    The three-flavor equivalent was covered from the start; this one was
    not, and it is the shape that the plotting modules use most.
    """
    h = random_hermitian(rng, 2)
    l_stack = np.linspace(0.1, 20.0, 50)

    batched = oscprob2nu.probabilities_2nu(as_nested_list(h), l_stack)
    scalar = np.array([oscprob2nu.probabilities_2nu(as_nested_list(h), l)
                       for l in l_stack])

    assert batched.shape == (50, 4)
    assert np.allclose(batched, scalar, atol=ATOL)


def test_2nu_many_hamiltonians_one_baseline(rng):
    r"""An array of 2nu Hamiltonians broadcasts against one baseline."""
    h_stack = stack_of_hermitians(rng, 40, 2)

    batched = oscprob2nu.probabilities_2nu(h_stack, 2.5)
    scalar = np.array([oscprob2nu.probabilities_2nu(as_nested_list(h), 2.5)
                       for h in h_stack])

    assert batched.shape == (40, 4)
    assert np.allclose(batched, scalar, atol=ATOL)


def test_2nu_two_dimensional_grid(rng):
    r"""Two flavors broadcast into a grid as well."""
    h_stack = stack_of_hermitians(rng, 6, 2)
    l_stack = np.linspace(0.5, 4.0, 5)

    grid = oscprob2nu.probabilities_2nu(h_stack[:, None, :, :],
                                        l_stack[None, :])
    assert grid.shape == (6, 5, 4)
    for i in range(6):
        for j in range(5):
            expected = oscprob2nu.probabilities_2nu(
                as_nested_list(h_stack[i]), l_stack[j])
            assert np.allclose(grid[i, j], expected, atol=ATOL)


def test_2nu_batched_evolution_operator_broadcasts(rng):
    r"""The 2nu evolution operator broadcasts the same way."""
    h = random_hermitian(rng, 2)
    l_stack = np.linspace(0.1, 8.0, 25)

    batched = oscprob2nu.evolution_operator_2nu(as_nested_list(h), l_stack)
    scalar = np.array([oscprob2nu.evolution_operator_2nu(as_nested_list(h), l)
                       for l in l_stack])

    assert batched.shape == (25, 2, 2)
    assert np.allclose(batched, scalar, atol=ATOL)


@pytest.mark.parametrize('module,size', [(oscprob2nu, 2), (oscprob3nu, 3)])
def test_incompatible_shapes_raise(module, size, rng):
    r"""Hamiltonians and baselines that do not broadcast are rejected,
    rather than silently producing a wrong shape."""
    h_stack = stack_of_hermitians(rng, 5, size)
    probabilities = getattr(module, 'probabilities_%dnu' % size)
    with pytest.raises(ValueError):
        probabilities(h_stack, np.zeros(4))


@pytest.mark.parametrize('module,size', [(oscprob2nu, 2), (oscprob3nu, 3)])
def test_three_dimensional_batch(module, size, rng):
    r"""Stacks with more than one leading axis work, and agree
    elementwise with the scalar path."""
    h_stack = stack_of_hermitians(rng, 12, size).reshape(3, 4, size, size)
    l_stack = rng.uniform(0.1, 5.0, (3, 4))
    probabilities = getattr(module, 'probabilities_%dnu' % size)

    batched = probabilities(h_stack, l_stack)
    assert batched.shape == (3, 4, size*size)
    for i in range(3):
        for j in range(4):
            expected = probabilities(as_nested_list(h_stack[i, j]),
                                     l_stack[i, j])
            assert np.allclose(batched[i, j], expected, atol=ATOL)
