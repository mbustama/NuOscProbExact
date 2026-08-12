# -*- coding: utf-8 -*-
r"""Tests of the sample two- and three-neutrino Hamiltonians.

These tests fix the *sign convention* of the vacuum Hamiltonians, which
matters as soon as a matter potential is added: an overall sign flip of
the vacuum term is invisible in vacuum but reverses the sign of the
matter potential relative to it, turning neutrinos into antineutrinos.
"""

import math

import numpy as np
import pytest

import globaldefs as gd
import hamiltonians2nu
import hamiltonians3nu
import hamiltonians4nu
from globaldefs import (D21_NO_BF, D31_NO_BF, DCP_NO_BF, EPS_2, EPS_3, LAMBDA,
                        S12_NO_BF, S13_NO_BF, S23_NO_BF, VCC_EARTH_CRUST)

from conftest import ATOL


def textbook_hamiltonian_2nu_vacuum(sth, dm2):
    r"""Returns the standard traceless 2nu vacuum Hamiltonian, times E.

    H = (Dm2/4) [[-cos(2t), sin(2t)], [sin(2t), cos(2t)]], the convention
    in which the mass eigenstate with the *larger* mass-squared value is
    the second one, so that a positive matter potential added to the ee
    entry describes neutrinos (not antineutrinos).
    """
    th = np.arcsin(sth)
    c2th, s2th = np.cos(2.0*th), np.sin(2.0*th)
    return dm2/4.0*np.array([[-c2th, s2th], [s2th, c2th]])


def test_mixing_matrix_2nu_is_orthogonal():
    r"""The 2nu rotation matrix is orthogonal with unit determinant."""
    for sth in [0.0, 0.2, 0.5, 0.8, 1.0]:
        r = np.array(hamiltonians2nu.mixing_matrix_2nu(sth))
        assert np.allclose(r @ r.T, np.eye(2), atol=ATOL)
        assert np.linalg.det(r) == pytest.approx(1.0, abs=ATOL)


def test_pmns_mixing_matrix_is_unitary():
    r"""The 3x3 PMNS matrix is unitary for a range of parameters."""
    for dcp in [0.0, 0.7, np.pi, 1.4*np.pi]:
        u = np.array(hamiltonians3nu.pmns_mixing_matrix(
            S12_NO_BF, S23_NO_BF, S13_NO_BF, dcp))
        assert np.allclose(u.conj().T @ u, np.eye(3), atol=ATOL)


def test_hamiltonian_2nu_vacuum_uses_the_standard_sign_convention():
    r"""The 2nu vacuum Hamiltonian matches the textbook convention."""
    for sth, dm2 in [(S12_NO_BF, D21_NO_BF), (S23_NO_BF, D31_NO_BF),
                     (0.3, 1.0e-3)]:
        h = np.array(
            hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(
                sth, dm2))
        assert np.allclose(h, textbook_hamiltonian_2nu_vacuum(sth, dm2),
                           atol=ATOL*abs(dm2))


def test_hamiltonian_3nu_vacuum_matches_u_m2_udagger():
    r"""H_vac = U.diag(0, D21, D31).U^dagger / 2."""
    u = np.array(hamiltonians3nu.pmns_mixing_matrix(
        S12_NO_BF, S23_NO_BF, S13_NO_BF, DCP_NO_BF))
    m2 = np.diag([0.0, D21_NO_BF, D31_NO_BF])
    h = np.array(hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
        S12_NO_BF, S23_NO_BF, S13_NO_BF, DCP_NO_BF, D21_NO_BF, D31_NO_BF))
    assert np.allclose(h, u @ m2 @ u.conj().T/2.0, atol=ATOL*D31_NO_BF)


@pytest.mark.parametrize('compute_matrix_multiplication', [False, True])
def test_vacuum_hamiltonian_branches_agree(compute_matrix_multiplication):
    r"""The closed-form and the live matrix-multiplication branches of
    both vacuum Hamiltonians give the same matrix."""
    h2_closed = np.array(
        hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(
            S12_NO_BF, D21_NO_BF, compute_matrix_multiplication=False))
    h2_matmul = np.array(
        hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(
            S12_NO_BF, D21_NO_BF, compute_matrix_multiplication=True))
    assert np.allclose(h2_closed, h2_matmul, atol=ATOL*D21_NO_BF)

    h3_closed = np.array(
        hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
            S12_NO_BF, S23_NO_BF, S13_NO_BF, DCP_NO_BF, D21_NO_BF, D31_NO_BF,
            compute_matrix_multiplication=False))
    h3_matmul = np.array(
        hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
            S12_NO_BF, S23_NO_BF, S13_NO_BF, DCP_NO_BF, D21_NO_BF, D31_NO_BF,
            compute_matrix_multiplication=True))
    assert np.allclose(h3_closed, h3_matmul, atol=ATOL*D31_NO_BF)


def test_matter_resonance_occurs_at_the_neutrino_resonance_energy():
    r"""The 2nu matter Hamiltonian resonates at E = Dm2 cos(2t)/(2 V).

    A sign error in the vacuum Hamiltonian would move the resonance to
    negative energies, i.e. to antineutrinos.
    """
    sth, dm2 = S12_NO_BF, D21_NO_BF
    th = np.arcsin(sth)
    energy_res = dm2*np.cos(2.0*th)/(2.0*VCC_EARTH_CRUST)
    assert energy_res > 0.0

    h_vacuum_energy_indep = \
        hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(sth, dm2)
    # At resonance the effective mixing is maximal: the diagonal entries
    # of the matter Hamiltonian become degenerate.
    h_res = np.array(hamiltonians2nu.hamiltonian_2nu_matter(
        h_vacuum_energy_indep, energy_res, VCC_EARTH_CRUST))
    assert np.real(h_res[0][0]) == pytest.approx(np.real(h_res[1][1]),
                                                 rel=1.0e-9)


def test_matter_hamiltonian_adds_the_potential_with_the_right_sign():
    r"""H_matter - H_vacuum/E = diag(V, 0), i.e. the potential enters
    the ee entry of the *standard* vacuum Hamiltonian."""
    energy = 1.0e9
    h_vacuum_energy_indep = \
        hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(
            S12_NO_BF, D21_NO_BF)
    h_matter = np.array(hamiltonians2nu.hamiltonian_2nu_matter(
        h_vacuum_energy_indep, energy, VCC_EARTH_CRUST))
    expected = textbook_hamiltonian_2nu_vacuum(S12_NO_BF, D21_NO_BF)/energy \
        + np.array([[VCC_EARTH_CRUST, 0.0], [0.0, 0.0]])
    assert np.allclose(h_matter, expected, atol=ATOL*VCC_EARTH_CRUST)


@pytest.mark.parametrize('sxi', [0.0, 0.3, 0.6, 0.9])
def test_hamiltonian_2nu_liv_is_finite_and_hermitian(sxi):
    r"""The 2nu LIV Hamiltonian is finite and Hermitian for any sin(xi).

    A cos(xi) computed as sqrt(1 - 2 sin(xi)) rather than
    sqrt(1 - sin(xi)^2) returns NaN for sin(xi) > 1/2.
    """
    h_vacuum_energy_indep = \
        hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(
            S23_NO_BF, D31_NO_BF)
    h = np.array(hamiltonians2nu.hamiltonian_2nu_liv(
        h_vacuum_energy_indep, 1.0e9, sxi, 1.0e-9, 2.0e-9, LAMBDA),
        dtype=complex)
    assert np.all(np.isfinite(h))
    assert np.allclose(h, h.conj().T, atol=ATOL)


@pytest.mark.parametrize('sxi', [0.0, 0.3, 0.6, 0.9])
def test_liv_term_has_the_prescribed_eigenvalues(sxi):
    r"""The LIV block of the 2nu Hamiltonian has eigenvalues b1, b2,
    scaled by E/Lambda: it is a rotation of diag(b1, b2)."""
    b1, b2, energy = 1.0e-9, 2.0e-9, 1.0e9
    zero_vacuum = [[0.0, 0.0], [0.0, 0.0]]
    h = np.array(hamiltonians2nu.hamiltonian_2nu_liv(
        zero_vacuum, energy, sxi, b1, b2, LAMBDA), dtype=complex)
    expected = np.sort(np.array([b1, b2])*energy/LAMBDA)
    assert np.allclose(np.sort(np.linalg.eigvalsh(h)), expected, atol=ATOL)


def test_liv_term_3nu_has_the_prescribed_eigenvalues():
    r"""The 3nu LIV term is a unitary rotation of diag(b1, b2, b3)."""
    b1, b2, b3, energy = 1.0e-9, 1.5e-9, 2.0e-9, 1.0e9
    zero_vacuum = np.zeros((3, 3), dtype=complex)
    h = np.array(hamiltonians3nu.hamiltonian_3nu_liv(
        zero_vacuum, energy, 0.3, 0.4, 0.5, 0.7, b1, b2, b3, LAMBDA),
        dtype=complex)
    expected = np.sort(np.array([b1, b2, b3])*energy/LAMBDA)
    assert np.allclose(np.sort(np.linalg.eigvalsh(h)), expected, atol=ATOL)


def test_nsi_hamiltonian_2nu_keeps_complex_epsilon():
    r"""A complex eps_em survives into the 2nu NSI Hamiltonian.

    If the vacuum Hamiltonian is stored as a real array, adding a
    complex NSI parameter silently discards its imaginary part.
    """
    eps_em = -0.06 + 0.03j
    h_vacuum_energy_indep = \
        hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(
            S23_NO_BF, D31_NO_BF)
    h = np.array(hamiltonians2nu.hamiltonian_2nu_nsi(
        h_vacuum_energy_indep, 1.0e9, VCC_EARTH_CRUST,
        [0.06, eps_em, 1.2]), dtype=complex)
    assert np.imag(h[0][1]) == \
        pytest.approx(VCC_EARTH_CRUST*eps_em.imag, abs=ATOL)
    assert np.allclose(h, h.conj().T, atol=ATOL)


def test_nsi_hamiltonian_3nu_keeps_complex_epsilon():
    r"""A complex eps_em, eps_et, eps_mt survive into the 3nu NSI
    Hamiltonian, which stays Hermitian."""
    eps = [0.06, -0.06+0.03j, 0.01-0.02j, 1.2, 0.005+0.001j, 0.0]
    h_vacuum_energy_indep = \
        hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
            S12_NO_BF, S23_NO_BF, S13_NO_BF, DCP_NO_BF, D21_NO_BF, D31_NO_BF)
    h = np.array(hamiltonians3nu.hamiltonian_3nu_nsi(
        h_vacuum_energy_indep, 1.0e9, VCC_EARTH_CRUST, eps), dtype=complex)
    assert np.imag(h[0][1]) == \
        pytest.approx(VCC_EARTH_CRUST*eps[1].imag, abs=ATOL)
    assert np.allclose(h, h.conj().T, atol=ATOL)


def test_all_sample_hamiltonians_are_hermitian():
    r"""Every sample Hamiltonian is Hermitian, so that U stays unitary."""
    energy = 1.0e9
    h2_vac = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(
        S23_NO_BF, D31_NO_BF)
    h3_vac = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
        S12_NO_BF, S23_NO_BF, S13_NO_BF, DCP_NO_BF, D21_NO_BF, D31_NO_BF)
    candidates = [
        hamiltonians2nu.hamiltonian_2nu_matter(h2_vac, energy,
                                               VCC_EARTH_CRUST),
        hamiltonians2nu.hamiltonian_2nu_nsi(h2_vac, energy, VCC_EARTH_CRUST,
                                            EPS_2),
        hamiltonians2nu.hamiltonian_2nu_liv(h2_vac, energy, 0.3, 1.0e-9,
                                            2.0e-9, LAMBDA),
        hamiltonians3nu.hamiltonian_3nu_matter(h3_vac, energy,
                                               VCC_EARTH_CRUST),
        hamiltonians3nu.hamiltonian_3nu_nsi(h3_vac, energy, VCC_EARTH_CRUST,
                                            EPS_3),
        hamiltonians3nu.hamiltonian_3nu_liv(h3_vac, energy, 0.3, 0.4, 0.5,
                                            0.7, 1.0e-9, 1.5e-9, 2.0e-9,
                                            LAMBDA),
    ]
    for h in candidates:
        a = np.array(h, dtype=complex)
        assert np.allclose(a, a.conj().T, atol=ATOL)


def test_hamiltonian_builders_do_not_mutate_their_input():
    r"""Building a matter/NSI/LIV Hamiltonian leaves the vacuum
    Hamiltonian passed in unchanged."""
    h2_vac = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(
        S23_NO_BF, D31_NO_BF)
    before = np.array(h2_vac, dtype=complex).copy()
    hamiltonians2nu.hamiltonian_2nu_matter(h2_vac, 1.0e9, VCC_EARTH_CRUST)
    hamiltonians2nu.hamiltonian_2nu_nsi(h2_vac, 1.0e9, VCC_EARTH_CRUST, EPS_2)
    hamiltonians2nu.hamiltonian_2nu_liv(h2_vac, 1.0e9, 0.3, 1.0e-9, 2.0e-9,
                                        LAMBDA)
    assert np.allclose(np.array(h2_vac, dtype=complex), before, atol=ATOL)


###############################################################################
# Angle conventions
###############################################################################

ANGLE_CASES = [
    (hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent,
     (gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF), (gd.DCP_NO_BF,),
     (gd.D21_NO_BF, gd.D31_NO_BF)),
    (hamiltonians3nu.pmns_mixing_matrix,
     (gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF), (gd.DCP_NO_BF,), ()),
    (hamiltonians4nu.mixing_matrix_4nu,
     (gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, 0.3, 0.3, 0.1),
     (gd.DCP_NO_BF,), ()),
    (hamiltonians2nu.mixing_matrix_2nu, (gd.S12_NO_BF,), (), ()),
]


@pytest.mark.parametrize('routine, sines, phases, rest', ANGLE_CASES)
def test_the_four_angle_conventions_describe_the_same_thing(
        routine, sines, phases, rest):
    r"""One parameter set, four spellings, one answer.

    This is the whole point of the keyword: the published numbers are
    not sines, and a reader who types them under the default gets a
    legal sine of a different angle rather than an error.
    """
    reference = np.asarray(routine(*sines, *phases, *rest))
    radians = tuple(math.asin(s) for s in sines)

    from_sin2 = routine(*(s*s for s in sines), *phases, *rest, angles='sin2')
    from_rad = routine(*radians, *phases, *rest, angles='rad')
    from_deg = routine(*(math.degrees(t) for t in radians),
                       *(math.degrees(p) for p in phases), *rest, angles='deg')

    assert np.array_equal(np.asarray(from_sin2), reference)
    assert np.array_equal(np.asarray(from_rad), reference)
    assert np.allclose(np.asarray(from_deg), reference, rtol=0.0, atol=1.0e-15)


def test_the_default_convention_is_unchanged():
    r"""Every frozen reference and every figure rests on this."""
    explicit = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, gd.DCP_NO_BF,
        gd.D21_NO_BF, gd.D31_NO_BF, angles='sin')
    implied = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, gd.DCP_NO_BF,
        gd.D21_NO_BF, gd.D31_NO_BF)

    assert np.array_equal(np.asarray(explicit), np.asarray(implied))


def test_the_published_mixing_parameters_can_be_typed_as_printed():
    r"""NuFit quotes sin^2; `globaldefs` carries the sines of the same.

    Passing the published number under the default is the error this
    keyword exists to make unnecessary, and it is silent: 0.310 is a
    perfectly good sine, of a different angle.
    """
    from_published = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
        0.310, 0.582, 0.02240, gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF,
        angles='sin2')
    from_constants = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, gd.DCP_NO_BF,
        gd.D21_NO_BF, gd.D31_NO_BF)

    assert np.allclose(np.asarray(from_published), np.asarray(from_constants),
                       rtol=0.0, atol=1.0e-12)
    # And the silent error it replaces really is silent
    wrong = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
        0.310, 0.582, 0.02240, gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF)
    assert not np.allclose(np.asarray(wrong), np.asarray(from_constants),
                           rtol=0.0, atol=1.0e-12)


@pytest.mark.parametrize('module', [hamiltonians2nu, hamiltonians3nu,
                                    hamiltonians4nu])
def test_every_convention_is_refused_outside_its_range(module):
    r"""And an unknown name is refused by name."""
    with pytest.raises(ValueError, match='must be one of'):
        module.hamiltonian_2nu_vacuum_energy_independent(
            0.5, gd.D21_NO_BF, angles='sine') if module is hamiltonians2nu \
            else module._sine_from(0.5, 'sine', 's12', 'test')
    with pytest.raises(ValueError, match=r"angles='sin2'"):
        module._sine_from(1.4, 'sin2', 's12', 'test')
    with pytest.raises(ValueError, match=r"angles='sin2'"):
        module._sine_from(-0.1, 'sin2', 's12', 'test')


@pytest.mark.parametrize('module', [hamiltonians2nu, hamiltonians3nu,
                                    hamiltonians4nu])
def test_a_phase_follows_the_angles_only_in_radians_or_degrees(module):
    r"""A phase has no sine, so `sin` and `sin2` leave it alone."""
    assert module._phase_from(2.5, 'sin') == 2.5
    assert module._phase_from(2.5, 'sin2') == 2.5
    assert module._phase_from(2.5, 'rad') == 2.5
    assert module._phase_from(180.0, 'deg') == pytest.approx(math.pi)
