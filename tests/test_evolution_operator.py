# -*- coding: utf-8 -*-
r"""Tests of the two- and three-neutrino time-evolution operators.

The evolution operator is the heart of the SU(2) and SU(3) expansions,
so it is checked against an independent computation of
:math:`e^{-i H_0 L}` (``scipy.linalg.expm``, with :math:`H_0` the
traceless part of the Hamiltonian), and against the defining properties
of a time-evolution operator: unitarity, :math:`U(0) = \mathbb{1}`, and
the group property.
"""

import numpy as np
import pytest
from scipy.linalg import expm

import oscprob2nu
import oscprob3nu

from conftest import ATOL, as_nested_list, traceless


def test_evolution_operator_2nu_matches_matrix_exponential(hermitian_2nu):
    r"""U_2(L) equals exp(-i H_0 L) for random Hermitian Hamiltonians."""
    for h_matrix, l in hermitian_2nu:
        u = np.array(oscprob2nu.evolution_operator_2nu(
            as_nested_list(h_matrix), l))
        assert np.allclose(u, expm(-1.0j*traceless(h_matrix)*l), atol=ATOL)


def test_evolution_operator_3nu_matches_matrix_exponential(hermitian_3nu):
    r"""U_3(L) equals exp(-i H_0 L) for random Hermitian Hamiltonians."""
    for h_matrix, l in hermitian_3nu:
        u = np.array(oscprob3nu.evolution_operator_3nu(
            as_nested_list(h_matrix), l))
        assert np.allclose(u, expm(-1.0j*traceless(h_matrix)*l), atol=ATOL)


@pytest.mark.parametrize('n,module,routine', [
    (2, oscprob2nu, 'evolution_operator_2nu'),
    (3, oscprob3nu, 'evolution_operator_3nu'),
])
def test_evolution_operator_is_unitary(n, module, routine, request):
    r"""U^dagger U = 1, for both the 2nu and the 3nu operators."""
    cases = request.getfixturevalue(
        'hermitian_2nu' if n == 2 else 'hermitian_3nu')
    evolve = getattr(module, routine)
    for h_matrix, l in cases:
        u = np.array(evolve(as_nested_list(h_matrix), l))
        assert np.allclose(u.conj().T @ u, np.eye(n), atol=ATOL)


@pytest.mark.parametrize('n,module,routine', [
    (2, oscprob2nu, 'evolution_operator_2nu'),
    (3, oscprob3nu, 'evolution_operator_3nu'),
])
def test_evolution_operator_at_zero_baseline_is_the_identity(
        n, module, routine, request):
    r"""U(L = 0) = 1."""
    cases = request.getfixturevalue(
        'hermitian_2nu' if n == 2 else 'hermitian_3nu')
    evolve = getattr(module, routine)
    for h_matrix, _ in cases[:20]:
        u = np.array(evolve(as_nested_list(h_matrix), 0.0))
        assert np.allclose(u, np.eye(n), atol=ATOL)


@pytest.mark.parametrize('n,module,routine', [
    (2, oscprob2nu, 'evolution_operator_2nu'),
    (3, oscprob3nu, 'evolution_operator_3nu'),
])
def test_evolution_operator_group_property(n, module, routine, request):
    r"""U(L_1 + L_2) = U(L_2) U(L_1), since H is time-independent."""
    cases = request.getfixturevalue(
        'hermitian_2nu' if n == 2 else 'hermitian_3nu')
    evolve = getattr(module, routine)
    for h_matrix, l in cases[:20]:
        h_list = as_nested_list(h_matrix)
        u_total = np.array(evolve(h_list, l))
        u_half = np.array(evolve(h_list, l/2.0))
        assert np.allclose(u_half @ u_half, u_total, atol=ATOL)


def test_u_coefficients_rebuild_the_2nu_operator(hermitian_2nu):
    r"""U_2 = u_0*1 + i*u_k*sigma_k reproduces the returned matrix."""
    sigma = [np.array([[0, 1], [1, 0]], dtype=complex),
             np.array([[0, -1j], [1j, 0]], dtype=complex),
             np.array([[1, 0], [0, -1]], dtype=complex)]
    for h_matrix, l in hermitian_2nu[:20]:
        h_list = as_nested_list(h_matrix)
        u_coeffs = oscprob2nu.evolution_operator_2nu_u_coefficients(h_list, l)
        rebuilt = u_coeffs[0]*np.eye(2) + \
            1.0j*sum(c*s for c, s in zip(u_coeffs[1:], sigma))
        assert np.allclose(
            rebuilt, np.array(oscprob2nu.evolution_operator_2nu(h_list, l)),
            atol=ATOL)


def test_u_coefficients_rebuild_the_3nu_operator(hermitian_3nu):
    r"""U_3 = u_0*1 + i*u_k*lambda_k reproduces the returned matrix."""
    from test_su3_algebra import GELL_MANN
    for h_matrix, l in hermitian_3nu[:20]:
        h_list = as_nested_list(h_matrix)
        u_coeffs = oscprob3nu.evolution_operator_3nu_u_coefficients(h_list, l)
        rebuilt = u_coeffs[0]*np.eye(3) + \
            1.0j*sum(c*lam for c, lam in zip(u_coeffs[1:], GELL_MANN))
        assert np.allclose(
            rebuilt, np.array(oscprob3nu.evolution_operator_3nu(h_list, l)),
            atol=ATOL)


def test_modulus_matches_numpy_norm(rng):
    r"""modulus() returns the Euclidean norm of the coefficient vector."""
    for _ in range(20):
        h_coeffs = list(rng.normal(size=3) + 1.0j*rng.normal(size=3))
        assert oscprob2nu.modulus(h_coeffs) == \
            pytest.approx(float(np.linalg.norm(h_coeffs)), abs=ATOL)
