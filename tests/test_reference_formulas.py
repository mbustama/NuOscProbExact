# -*- coding: utf-8 -*-
r"""Tests that the standard oscillation formulas agree with the exact
SU(2) and SU(3) expansions.

The routines ``probabilities_*_std`` implement the textbook oscillation
formulas and exist to validate the exact computation.  They are only
useful if they agree with it to round-off; a rounded numerical prefactor
in the oscillation phase --- 1.27 instead of 1.266933... --- shifts the
phase by 0.24%, which near an oscillation minimum shows up as a tens-of-
percent discrepancy in the probability itself.

All routines in the library take the neutrino energy in eV and the
baseline in eV^{-1}.
"""

import numpy as np
import pytest

import hamiltonians2nu
import hamiltonians3nu
import oscprob2nu
import oscprob3nu
from globaldefs import (CONV_KM_TO_INV_EV, D21_NO_BF, D31_NO_BF, DCP_NO_BF,
                        S12_NO_BF, S13_NO_BF, S23_NO_BF, VCC_EARTH_CRUST)


# Baselines [km] and energies [GeV] spanning several oscillation cycles,
# including points close to an oscillation minimum, where a small phase
# error is most visible.
CASES = [(1.0e3, 1.0), (300.0, 2.0), (1.3e3, 3.0), (500.0, 0.7),
         (7.0e3, 5.0), (100.0, 0.4)]

TOL = 1.0e-12


@pytest.mark.parametrize('l_km,energy_gev', CASES)
def test_2nu_vacuum_std_matches_the_exact_computation(l_km, energy_gev):
    r"""The standard 2nu vacuum formula reproduces the SU(2) result."""
    energy = energy_gev*1.0e9
    baseline = l_km*CONV_KM_TO_INV_EV

    h_vacuum = np.multiply(
        1.0/energy,
        hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(
            S23_NO_BF, D31_NO_BF))
    exact = oscprob2nu.probabilities_2nu(h_vacuum, baseline)
    std = hamiltonians2nu.probabilities_2nu_vacuum_std(
        S23_NO_BF, D31_NO_BF, energy, baseline)

    assert np.allclose(np.array(exact), np.array(std), atol=TOL)


@pytest.mark.parametrize('l_km,energy_gev', CASES)
def test_2nu_matter_std_matches_the_exact_computation(l_km, energy_gev):
    r"""The standard 2nu matter formula reproduces the SU(2) result.

    This also pins the sign of the matter potential relative to the
    vacuum Hamiltonian: the standard formula describes neutrinos.
    """
    energy = energy_gev*1.0e9
    baseline = l_km*CONV_KM_TO_INV_EV

    h_vacuum_energy_indep = \
        hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(
            S23_NO_BF, D31_NO_BF)
    h_matter = hamiltonians2nu.hamiltonian_2nu_matter(
        h_vacuum_energy_indep, energy, VCC_EARTH_CRUST)
    exact = oscprob2nu.probabilities_2nu(h_matter, baseline)
    std = hamiltonians2nu.probabilities_2nu_matter_std(
        S23_NO_BF, D31_NO_BF, VCC_EARTH_CRUST, energy, baseline)

    assert np.allclose(np.array(exact), np.array(std), atol=1.0e-10)


@pytest.mark.parametrize('l_km,energy_gev', CASES)
def test_3nu_vacuum_std_matches_the_exact_computation(l_km, energy_gev):
    r"""The standard 3nu vacuum formula reproduces the SU(3) result."""
    energy = energy_gev*1.0e9
    baseline = l_km*CONV_KM_TO_INV_EV

    h_vacuum = np.multiply(
        1.0/energy,
        hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
            S12_NO_BF, S23_NO_BF, S13_NO_BF, DCP_NO_BF, D21_NO_BF, D31_NO_BF))
    exact = oscprob3nu.probabilities_3nu(h_vacuum, baseline)

    pmns = hamiltonians3nu.pmns_mixing_matrix(S12_NO_BF, S23_NO_BF, S13_NO_BF,
                                              DCP_NO_BF)
    std = hamiltonians3nu.probabilities_3nu_vacuum_std(
        pmns, D21_NO_BF, D31_NO_BF, energy, baseline)

    assert np.allclose(np.array(exact), np.array(std), atol=TOL)


@pytest.mark.parametrize('l_km,energy_gev', CASES)
def test_std_probabilities_are_normalized(l_km, energy_gev):
    r"""The standard formulas return probabilities that sum to one."""
    energy = energy_gev*1.0e9
    baseline = l_km*CONV_KM_TO_INV_EV

    std_2nu = np.array(hamiltonians2nu.probabilities_2nu_vacuum_std(
        S23_NO_BF, D31_NO_BF, energy, baseline)).reshape(2, 2)
    assert np.allclose(std_2nu.sum(axis=1), np.ones(2), atol=TOL)

    pmns = hamiltonians3nu.pmns_mixing_matrix(S12_NO_BF, S23_NO_BF, S13_NO_BF,
                                              DCP_NO_BF)
    std_3nu = np.array(hamiltonians3nu.probabilities_3nu_vacuum_std(
        pmns, D21_NO_BF, D31_NO_BF, energy, baseline)).reshape(3, 3)
    assert np.allclose(std_3nu.sum(axis=1), np.ones(3), atol=1.0e-10)


def test_delta_and_jarlskog_helpers():
    r"""The Kronecker delta and the quartic product J behave as
    advertised."""
    assert hamiltonians3nu.delta(0, 0) == 1
    assert hamiltonians3nu.delta(0, 1) == 0

    pmns = hamiltonians3nu.pmns_mixing_matrix(S12_NO_BF, S23_NO_BF, S13_NO_BF,
                                              DCP_NO_BF)
    u = np.array(pmns)
    for alpha, beta, k, j in [(0, 1, 1, 0), (1, 2, 2, 0), (0, 2, 2, 1)]:
        expected = u[alpha][k].conjugate()*u[beta][k] \
            * u[alpha][j]*u[beta][j].conjugate()
        assert hamiltonians3nu.J(pmns, alpha, beta, k, j) == \
            pytest.approx(expected, abs=TOL)
        # J with k == j is real and positive
        assert np.imag(hamiltonians3nu.J(pmns, alpha, beta, k, k)) == \
            pytest.approx(0.0, abs=TOL)
