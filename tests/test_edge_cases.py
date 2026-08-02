# -*- coding: utf-8 -*-
r"""Tests of numerically awkward inputs.

Degenerate Hamiltonians make the SU(3) invariant :math:`|h|^2` vanish,
which the latent-root formula divides by; near-degenerate ones can push
the argument of the arc cosine marginally outside :math:`[-1, 1]`
through round-off.  Neither may produce NaN or a non-unitary operator.
"""

import importlib

import numpy as np
import pytest

import hamiltonians2nu
import hamiltonians3nu
import oscprob2nu
import oscprob3nu

from conftest import ATOL, as_nested_list, random_hermitian


DEGENERATE_3NU = [
    ('proportional to the identity', np.eye(3, dtype=complex)),
    ('zero', np.zeros((3, 3), dtype=complex)),
    ('large multiple of the identity', 1.0e6*np.eye(3, dtype=complex)),
    ('two degenerate eigenvalues', np.diag([1.0, 1.0, -2.0]).astype(complex)),
    ('nearly degenerate', np.diag([1.0, 1.0 + 1.0e-14, -2.0]).astype(complex)),
]


@pytest.mark.parametrize('label,h_matrix', DEGENERATE_3NU,
                         ids=[c[0] for c in DEGENERATE_3NU])
@pytest.mark.parametrize('l', [0.0, 1.0, 1.0e3])
def test_3nu_evolution_operator_is_finite_for_degenerate_hamiltonians(
        label, h_matrix, l):
    r"""U_3 stays finite and unitary for degenerate Hamiltonians."""
    u = np.array(oscprob3nu.evolution_operator_3nu(
        as_nested_list(h_matrix), l))
    assert np.all(np.isfinite(u)), 'non-finite U for a %s Hamiltonian' % label
    assert np.allclose(u.conj().T @ u, np.eye(3), atol=1.0e-10)


@pytest.mark.parametrize('label,h_matrix', DEGENERATE_3NU,
                         ids=[c[0] for c in DEGENERATE_3NU])
def test_3nu_probabilities_are_finite_for_degenerate_hamiltonians(
        label, h_matrix):
    r"""The 3nu probabilities stay finite and normalized."""
    prob = np.array(oscprob3nu.probabilities_3nu(
        as_nested_list(h_matrix), 1.0)).reshape(3, 3)
    assert np.all(np.isfinite(prob))
    assert np.allclose(prob.sum(axis=1), np.ones(3), atol=1.0e-10)


def test_identity_hamiltonian_gives_no_oscillations():
    r"""A Hamiltonian proportional to the identity has no traceless
    part, so no flavor transitions occur."""
    for scale in [0.0, 1.0, 1.0e6]:
        prob = np.array(oscprob3nu.probabilities_3nu(
            as_nested_list(scale*np.eye(3, dtype=complex)), 10.0)
        ).reshape(3, 3)
        assert np.allclose(prob, np.eye(3), atol=1.0e-10)


@pytest.mark.parametrize('h_matrix', [
    np.eye(2, dtype=complex),
    np.zeros((2, 2), dtype=complex),
    1.0e6*np.eye(2, dtype=complex),
])
def test_2nu_evolution_operator_is_finite_for_degenerate_hamiltonians(
        h_matrix):
    r"""U_2 stays finite and unitary when |h| vanishes."""
    u = np.array(oscprob2nu.evolution_operator_2nu(
        as_nested_list(h_matrix), 1.0))
    assert np.all(np.isfinite(u))
    assert np.allclose(u.conj().T @ u, np.eye(2), atol=1.0e-10)


def test_2nu_probabilities_are_finite_for_degenerate_hamiltonians():
    r"""The 2nu probabilities stay finite when |h| vanishes."""
    prob = np.array(oscprob2nu.probabilities_2nu(
        as_nested_list(np.eye(2, dtype=complex)), 1.0))
    assert np.all(np.isfinite(prob))
    assert np.allclose(prob.reshape(2, 2), np.eye(2), atol=1.0e-10)


def test_psi_roots_are_real_for_hermitian_hamiltonians(rng):
    r"""The latent roots of a Hermitian Hamiltonian are real.

    Round-off can push the arc-cosine argument marginally outside
    [-1, 1]; the implementation must clip it rather than return complex
    roots, which would spoil unitarity.
    """
    from conftest import random_hermitian
    for _ in range(200):
        h_matrix = random_hermitian(rng, 3)
        h_coeffs = oscprob3nu.hamiltonian_3nu_coefficients(
            as_nested_list(h_matrix))
        h2, h3 = oscprob3nu.su3_invariants(h_coeffs)
        psi = np.array(oscprob3nu.psi_roots(h2, h3))
        assert np.allclose(np.imag(psi), 0.0, atol=1.0e-9)


@pytest.mark.parametrize('h2,h3', [
    (1.0, 2.0/3.0/np.sqrt(3.0)*3.0),   # arc-cosine argument at -1
    (1.0, -2.0/3.0/np.sqrt(3.0)*3.0),  # arc-cosine argument at +1
    (1.0, 10.0),                       # argument far outside [-1, 1]
    (1.0, -10.0),
])
def test_psi_roots_clip_the_arc_cosine_argument(h2, h3):
    r"""psi_roots returns real, finite roots even when its argument is
    pushed outside the domain of the arc cosine."""
    psi = np.array(oscprob3nu.psi_roots(h2, h3))
    assert np.all(np.isfinite(psi))
    assert np.allclose(np.imag(psi), 0.0, atol=ATOL)


def test_large_baseline_keeps_unitarity(rng):
    r"""Unitarity survives many oscillation cycles."""
    from conftest import random_hermitian
    for _ in range(20):
        h_matrix = random_hermitian(rng, 3)
        u = np.array(oscprob3nu.evolution_operator_3nu(
            as_nested_list(h_matrix), 1.0e6))
        assert np.allclose(u.conj().T @ u, np.eye(3), atol=1.0e-8)


# --------------------------------------------------------------------------
# Input the expansions cannot make sense of, which used to pass silently
# --------------------------------------------------------------------------

NOT_HERMITIAN = {
    2: [[1.0+0.j, 2.0+0.j], [0.0+0.j, 1.0+0.j]],
    3: [[1.0+0.j, 2.0+0.j, 0.j], [0.j, 1.0+0.j, 0.j], [0.j, 0.j, 1.0+0.j]],
    4: [[1.0+0.j, 2.0+0.j, 0.j, 0.j], [0.j, 1.0+0.j, 0.j, 0.j],
        [0.j, 0.j, 1.0+0.j, 0.j], [0.j, 0.j, 0.j, 1.0+0.j]],
}


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
@pytest.mark.parametrize('routine', ['probabilities', 'evolution_operator'])
def test_a_non_hermitian_hamiltonian_is_refused(n_flavors, routine):
    r"""A non-Hermitian Hamiltonian raises rather than returning numbers.

    This is the failure that had to be caught by checking, because
    nothing downstream reveals it: the expansion returns probabilities
    that still **sum to one**, so the sanity check a caller would
    actually apply cannot tell that the answer is meaningless.  Before
    the check existed, the matrix below returned 1.000000 for the first
    row and no warning of any kind.
    """
    module = importlib.import_module('oscprob%dnu' % n_flavors)
    call = getattr(module, '%s_%dnu' % (routine, n_flavors))

    with pytest.raises(ValueError, match='not Hermitian'):
        call(NOT_HERMITIAN[n_flavors], 1.0)


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
def test_a_non_hermitian_stack_is_refused(n_flavors):
    r"""The batched path checks too, not only the scalar one."""
    module = importlib.import_module('oscprob%dnu' % n_flavors)
    call = getattr(module, 'probabilities_%dnu' % n_flavors)
    stack = np.stack([np.array(NOT_HERMITIAN[n_flavors], dtype=complex)]*50)

    with pytest.raises(ValueError, match='not Hermitian'):
        call(stack, np.full(50, 1.0))


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
def test_the_hermiticity_check_can_be_switched_off(n_flavors):
    r"""The switch is the escape hatch for scans that cannot afford it.

    Checking costs a pass over the stack, which is the same order as the
    evaluation; `CHECK_HERMITICITY` exists so that a caller whose
    Hamiltonians are Hermitian by construction can decline it.
    """
    module = importlib.import_module('oscprob%dnu' % n_flavors)
    call = getattr(module, 'probabilities_%dnu' % n_flavors)

    original = module.CHECK_HERMITICITY
    try:
        module.CHECK_HERMITICITY = False
        prob = call(NOT_HERMITIAN[n_flavors], 1.0)
    finally:
        module.CHECK_HERMITICITY = original

    assert len(prob) == n_flavors*n_flavors


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
def test_hamiltonians_the_library_builds_pass_the_check(n_flavors, rng):
    r"""The tolerance is loose enough for matrices assembled in floating
    point, which is the whole difficulty of checking this at all."""
    module = importlib.import_module('oscprob%dnu' % n_flavors)
    call = getattr(module, 'probabilities_%dnu' % n_flavors)

    for _ in range(20):
        h_matrix = random_hermitian(rng, n_flavors)
        prob = np.asarray(call(as_nested_list(h_matrix), 1.5))
        assert np.isclose(prob[:n_flavors].sum(), 1.0, atol=ATOL)


SINE_OUT_OF_RANGE = [
    ('hamiltonians2nu', 'mixing_matrix_2nu', (1.5,)),
    ('hamiltonians3nu', 'pmns_mixing_matrix', (1.5, 0.1, 0.1, 0.0)),
    ('hamiltonians4nu', 'mixing_matrix_4nu',
     (1.5, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0)),
]


@pytest.mark.parametrize('module_name,routine,args', SINE_OUT_OF_RANGE)
def test_a_sine_outside_its_range_is_refused_the_same_way(module_name,
                                                          routine, args):
    r"""All three flavor counts refuse it, and say the same thing.

    They did not.  Two and three flavors took the cosine with
    :mod:`math` and raised ``math domain error``, which names neither
    the parameter nor the value; four flavors took it with
    :func:`numpy.sqrt`, which returns ``nan`` and let it run silently
    into the probabilities.  Same invalid input, three behaviours.
    """
    module = importlib.import_module(module_name)

    with pytest.raises(ValueError, match='sine of an angle'):
        getattr(module, routine)(*args)


def test_a_sine_of_exactly_one_is_still_allowed():
    r"""The boundary is inclusive: maximal mixing is a legal angle."""
    assert np.isclose(hamiltonians2nu.mixing_matrix_2nu(1.0)[0][0], 0.0)
    assert np.isclose(hamiltonians3nu.pmns_mixing_matrix(1.0, 0.1, 0.1,
                                                         0.0)[0][0].real, 0.0)
