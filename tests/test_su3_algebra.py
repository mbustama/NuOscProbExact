# -*- coding: utf-8 -*-
r"""Tests of the SU(3) algebra: the d tensor, the star product, and the
two SU(3) invariants.

The d tensor is hard-coded in :mod:`oscprob3nu`; these tests check it
against its definition, :math:`d_{ijk} = \frac{1}{4}
\mathrm{Tr}(\{\lambda_i,\lambda_j\}\lambda_k)`.
"""

import numpy as np
import pytest

import oscprob3nu

from conftest import ATOL, as_nested_list, random_hermitian


SQRT3 = np.sqrt(3.0)

# The eight Gell-Mann matrices, written out independently of the module
# under test.
GELL_MANN = [
    np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex),
    np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex),
    np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex),
    np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex),
    np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], dtype=complex),
    np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex),
    np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=complex),
    np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]], dtype=complex)/SQRT3,
]


def d_reference(i, j, k):
    r"""Returns d_ijk computed from its definition."""
    lam_i, lam_j, lam_k = GELL_MANN[i], GELL_MANN[j], GELL_MANN[k]
    anticommutator = lam_i @ lam_j + lam_j @ lam_i
    return 0.25*np.trace(anticommutator @ lam_k).real


@pytest.mark.parametrize('i', range(8))
def test_tensor_d_matches_definition(i):
    r"""Every d_ijk agrees with (1/4)Tr({lambda_i,lambda_j}lambda_k)."""
    for j in range(8):
        for k in range(8):
            assert oscprob3nu.tensor_d(i, j, k) == \
                pytest.approx(d_reference(i, j, k), abs=ATOL)


def test_tensor_d_is_totally_symmetric():
    r"""d_ijk is invariant under any permutation of its indices."""
    d = np.array([[[oscprob3nu.tensor_d(i, j, k) for k in range(8)]
                   for j in range(8)] for i in range(8)])
    assert np.allclose(d, d.transpose(1, 0, 2), atol=ATOL)
    assert np.allclose(d, d.transpose(0, 2, 1), atol=ATOL)
    assert np.allclose(d, d.transpose(2, 1, 0), atol=ATOL)


def test_tensor_d_is_traceless_in_first_two_indices():
    r"""d_iik = 0, a standard property of the SU(3) d tensor."""
    for k in range(8):
        assert sum(oscprob3nu.tensor_d(i, i, k) for i in range(8)) == \
            pytest.approx(0.0, abs=ATOL)


@pytest.mark.parametrize('bad', [(8, 0, 0), (-1, 0, 0), (99, 3, 4)])
def test_tensor_d_rejects_out_of_range_indices(bad):
    r"""Out-of-range indices raise instead of silently returning None."""
    with pytest.raises((IndexError, ValueError)):
        oscprob3nu.tensor_d(*bad)


def test_star_product_matches_explicit_sum(rng):
    r"""(h*h)_i equals d_ijk h_j h_k summed explicitly."""
    for _ in range(20):
        h = rng.normal(size=8)
        for i in range(8):
            expected = sum(oscprob3nu.tensor_d(i, j, k)*h[j]*h[k]
                           for j in range(8) for k in range(8))
            assert oscprob3nu.star(i, h) == pytest.approx(expected, abs=ATOL)


def test_su3_invariants_against_matrix_traces(rng):
    r"""|h|^2 and <h> reproduce Tr(H^2)/2 and Tr(H^3)/2 of the traceless
    Hamiltonian, which is how the two SU(3) invariants are defined."""
    for _ in range(50):
        h_matrix = random_hermitian(rng, 3)
        h_matrix -= np.trace(h_matrix)/3.0*np.eye(3)
        h_coeffs = oscprob3nu.hamiltonian_3nu_coefficients(
            as_nested_list(h_matrix))
        h2, h3 = oscprob3nu.su3_invariants(h_coeffs)
        assert np.real(h2) == \
            pytest.approx(np.trace(h_matrix @ h_matrix).real/2.0, abs=ATOL)
        assert np.real(h3) == pytest.approx(
            np.trace(h_matrix @ h_matrix @ h_matrix).real/2.0, abs=ATOL)


def test_su3_invariant_h3_equals_h_dot_star(rng):
    r"""<h> = h_i (h*h)_i, the contraction the implementation may exploit."""
    for _ in range(20):
        h = rng.normal(size=8)
        _, h3 = oscprob3nu.su3_invariants(h)
        star = np.array([oscprob3nu.star(i, h) for i in range(8)])
        assert np.real(h3) == pytest.approx(float(np.dot(h, star)), abs=ATOL)


def test_hamiltonian_3nu_coefficients_reconstruct_the_matrix(rng):
    r"""H_traceless = h_k lambda_k recovers the original Hamiltonian."""
    for _ in range(50):
        h_matrix = random_hermitian(rng, 3)
        h_matrix -= np.trace(h_matrix)/3.0*np.eye(3)
        h_coeffs = oscprob3nu.hamiltonian_3nu_coefficients(
            as_nested_list(h_matrix))
        rebuilt = sum(np.real(c)*lam
                      for c, lam in zip(h_coeffs, GELL_MANN))
        assert np.allclose(rebuilt, h_matrix, atol=ATOL)


def test_psi_roots_are_eigenvalues_of_minus_the_hamiltonian(rng):
    r"""The latent roots psi are the eigenvalues of -H_traceless."""
    for _ in range(50):
        h_matrix = random_hermitian(rng, 3)
        h_matrix -= np.trace(h_matrix)/3.0*np.eye(3)
        h_coeffs = oscprob3nu.hamiltonian_3nu_coefficients(
            as_nested_list(h_matrix))
        h2, h3 = oscprob3nu.su3_invariants(h_coeffs)
        psi = np.sort(np.real(oscprob3nu.psi_roots(h2, h3)))
        expected = np.sort(-np.linalg.eigvalsh(h_matrix))
        assert np.allclose(psi, expected, atol=1.0e-10)


def test_star_all_matches_the_tensor_contraction(rng):
    r"""The hand-written sparse expansion agrees with the d tensor.

    ``_star_all`` writes out :math:`(h \star h)_i = d_{ijk} h_j h_k`
    term by term, because contracting the dense table costs more NumPy
    dispatch than the arithmetic is worth for eight numbers.  That
    expansion has to stay in step with `tensor_d`, which is what this
    checks.
    """
    for _ in range(200):
        h = list(rng.normal(size=8))
        explicit = oscprob3nu._star_all(h)
        reference = [oscprob3nu.star(i, h) for i in range(8)]
        assert np.allclose(explicit, reference, atol=ATOL)


def test_star_all_handles_arrays_and_lists(rng):
    r"""_star_all accepts either a list or an array of coefficients."""
    h = rng.normal(size=8)
    assert np.allclose(oscprob3nu._star_all(list(h)),
                       oscprob3nu._star_all(h), atol=ATOL)
