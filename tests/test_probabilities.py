# -*- coding: utf-8 -*-
r"""Tests of the two- and three-neutrino oscillation probabilities.

These tests pin down the properties that any set of oscillation
probabilities must satisfy --- they lie in [0, 1] and they sum to one
over the final flavor --- and check the closed-form probabilities
against the modulus-squared entries of the corresponding evolution
operator.

The two-neutrino case is checked with *complex* off-diagonal
Hamiltonian entries as well as real ones: with a real off-diagonal
entry the coefficient :math:`h_2` vanishes, which hides any error in
its contribution to :math:`P_{e\mu}`.
"""

import numpy as np
import pytest

import oscprob2nu
import oscprob3nu

from conftest import ATOL, as_nested_list, random_hermitian


def test_probabilities_2nu_match_the_evolution_operator(hermitian_2nu):
    r"""P_ab equals |U_ba|^2, including for complex Hamiltonians."""
    for h_matrix, l in hermitian_2nu:
        h_list = as_nested_list(h_matrix)
        p_ee, p_em, p_me, p_mm = oscprob2nu.probabilities_2nu(h_list, l)
        u = np.array(oscprob2nu.evolution_operator_2nu(h_list, l))
        assert p_ee == pytest.approx(abs(u[0][0])**2.0, abs=ATOL)
        assert p_em == pytest.approx(abs(u[1][0])**2.0, abs=ATOL)
        assert p_me == pytest.approx(abs(u[0][1])**2.0, abs=ATOL)
        assert p_mm == pytest.approx(abs(u[1][1])**2.0, abs=ATOL)


def test_probabilities_2nu_need_the_h2_coefficient(rng):
    r"""A purely imaginary off-diagonal entry still drives oscillations.

    With H_12 purely imaginary, h_1 = 0 and h_2 != 0.  Dropping the
    h_2 contribution would make P_em vanish identically here.
    """
    h_list = [[1.0+0.0j, 0.0+2.0j], [0.0-2.0j, 3.0+0.0j]]
    p_em = oscprob2nu.probabilities_2nu(h_list, 1.0)[1]
    assert p_em > 0.1


def test_probabilities_3nu_match_the_evolution_operator(hermitian_3nu):
    r"""P_ab equals |U_ba|^2 for the three-neutrino probabilities."""
    for h_matrix, l in hermitian_3nu:
        h_list = as_nested_list(h_matrix)
        prob = np.array(oscprob3nu.probabilities_3nu(h_list, l)).reshape(3, 3)
        u = np.array(oscprob3nu.evolution_operator_3nu(h_list, l))
        assert np.allclose(prob, (np.abs(u)**2.0).T, atol=ATOL)


@pytest.mark.parametrize('n', [2, 3])
def test_probabilities_sum_to_one_over_final_flavor(n, request):
    r"""Sum_b P(nu_a -> nu_b) = 1, from unitarity of U."""
    cases = request.getfixturevalue(
        'hermitian_2nu' if n == 2 else 'hermitian_3nu')
    for h_matrix, l in cases:
        h_list = as_nested_list(h_matrix)
        if n == 2:
            prob = np.array(oscprob2nu.probabilities_2nu(h_list, l))
        else:
            prob = np.array(oscprob3nu.probabilities_3nu(h_list, l))
        prob = prob.reshape(n, n)
        assert np.allclose(prob.sum(axis=1), np.ones(n), atol=ATOL)
        assert np.allclose(prob.sum(axis=0), np.ones(n), atol=ATOL)


@pytest.mark.parametrize('n', [2, 3])
def test_probabilities_lie_in_the_unit_interval(n, request):
    r"""Every probability lies in [0, 1]."""
    cases = request.getfixturevalue(
        'hermitian_2nu' if n == 2 else 'hermitian_3nu')
    for h_matrix, l in cases:
        h_list = as_nested_list(h_matrix)
        if n == 2:
            prob = np.array(oscprob2nu.probabilities_2nu(h_list, l))
        else:
            prob = np.array(oscprob3nu.probabilities_3nu(h_list, l))
        assert np.all(prob >= -ATOL)
        assert np.all(prob <= 1.0 + ATOL)


@pytest.mark.parametrize('n', [2, 3])
def test_probabilities_at_zero_baseline_are_the_identity(n, request):
    r"""P(nu_a -> nu_b, L = 0) = delta_ab."""
    cases = request.getfixturevalue(
        'hermitian_2nu' if n == 2 else 'hermitian_3nu')
    for h_matrix, _ in cases[:20]:
        h_list = as_nested_list(h_matrix)
        if n == 2:
            prob = np.array(oscprob2nu.probabilities_2nu(h_list, 0.0))
        else:
            prob = np.array(oscprob3nu.probabilities_3nu(h_list, 0.0))
        assert np.allclose(prob.reshape(n, n), np.eye(n), atol=ATOL)


def test_probabilities_are_symmetric_for_real_hamiltonians(rng):
    r"""P_ab = P_ba when the Hamiltonian is real, i.e. when CP is
    conserved."""
    for _ in range(30):
        a = rng.normal(size=(3, 3))
        h_matrix = ((a + a.T)/2.0).astype(complex)
        l = rng.uniform(0.1, 5.0)
        prob = np.array(oscprob3nu.probabilities_3nu(
            as_nested_list(h_matrix), l)).reshape(3, 3)
        assert np.allclose(prob, prob.T, atol=ATOL)


def test_2nu_probability_is_a_special_case_of_the_3nu_probability(rng):
    r"""A 3nu Hamiltonian that is block diagonal, with the third flavor
    decoupled, reproduces the 2nu probabilities in the first block."""
    for _ in range(20):
        block = random_hermitian(rng, 2)
        l = rng.uniform(0.1, 5.0)
        h3 = np.zeros((3, 3), dtype=complex)
        h3[:2, :2] = block
        prob_3nu = np.array(oscprob3nu.probabilities_3nu(
            as_nested_list(h3), l)).reshape(3, 3)
        prob_2nu = np.array(oscprob2nu.probabilities_2nu(
            as_nested_list(block), l)).reshape(2, 2)
        assert np.allclose(prob_3nu[:2, :2], prob_2nu, atol=1.0e-10)
        assert prob_3nu[2, 2] == pytest.approx(1.0, abs=1.0e-10)
