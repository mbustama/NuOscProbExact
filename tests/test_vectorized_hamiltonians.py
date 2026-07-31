# -*- coding: utf-8 -*-
r"""Tests of the sample Hamiltonians built for an array of energies.

The matter, NSI and LIV Hamiltonians all take a neutrino energy.  Given
an array of energies they return one Hamiltonian per energy, stacked
along a leading axis, which is the shape that
:func:`oscprob3nu.probabilities_3nu` consumes --- so an energy scan is
two vectorised calls rather than a Python loop over both.
"""

import numpy as np
import pytest

import hamiltonians2nu
import hamiltonians3nu
import oscprob2nu
import oscprob3nu
from globaldefs import (CONV_KM_TO_INV_EV, D21_NO_BF, D31_NO_BF, DCP_NO_BF,
                        EPS_2, EPS_3, LAMBDA, S12_NO_BF, S13_NO_BF, S23_NO_BF,
                        VCC_EARTH_CRUST)

from conftest import ATOL


ENERGIES = np.logspace(-1.0, 1.0, 13)*1.0e9      # [eV]


@pytest.fixture(scope='module')
def h_vacuum_2nu():
    r"""Returns the energy-independent 2nu vacuum Hamiltonian."""
    return hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(
        S23_NO_BF, D31_NO_BF)


@pytest.fixture(scope='module')
def h_vacuum_3nu():
    r"""Returns the energy-independent 3nu vacuum Hamiltonian."""
    return hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
        S12_NO_BF, S23_NO_BF, S13_NO_BF, DCP_NO_BF, D21_NO_BF, D31_NO_BF)


def builders_3nu(h_vacuum):
    r"""Returns the three-flavor builders, each bound to its parameters."""
    return [
        ('matter', lambda e: hamiltonians3nu.hamiltonian_3nu_matter(
            h_vacuum, e, VCC_EARTH_CRUST)),
        ('nsi', lambda e: hamiltonians3nu.hamiltonian_3nu_nsi(
            h_vacuum, e, VCC_EARTH_CRUST, EPS_3)),
        ('nsi complex', lambda e: hamiltonians3nu.hamiltonian_3nu_nsi(
            h_vacuum, e, VCC_EARTH_CRUST,
            [0.06, -0.06+0.03j, 0.01-0.02j, 1.2, 0.005+0.001j, 0.0])),
        ('liv', lambda e: hamiltonians3nu.hamiltonian_3nu_liv(
            h_vacuum, e, 0.3, 0.4, 0.5, 0.7, 1.0e-9, 1.5e-9, 2.0e-9,
            LAMBDA)),
    ]


def builders_2nu(h_vacuum):
    r"""Returns the two-flavor builders, each bound to its parameters."""
    return [
        ('matter', lambda e: hamiltonians2nu.hamiltonian_2nu_matter(
            h_vacuum, e, VCC_EARTH_CRUST)),
        ('nsi', lambda e: hamiltonians2nu.hamiltonian_2nu_nsi(
            h_vacuum, e, VCC_EARTH_CRUST, EPS_2)),
        ('nsi complex', lambda e: hamiltonians2nu.hamiltonian_2nu_nsi(
            h_vacuum, e, VCC_EARTH_CRUST, [0.06, -0.06+0.03j, 1.2])),
        ('liv', lambda e: hamiltonians2nu.hamiltonian_2nu_liv(
            h_vacuum, e, 0.3, 1.0e-9, 2.0e-9, LAMBDA)),
    ]


def test_3nu_builders_match_the_scalar_loop(h_vacuum_3nu):
    r"""One call with an array of energies equals a loop over them."""
    for name, build in builders_3nu(h_vacuum_3nu):
        batched = build(ENERGIES)
        scalar = np.stack([build(e) for e in ENERGIES])
        assert batched.shape == (len(ENERGIES), 3, 3), name
        assert np.array_equal(batched, scalar), name


def test_2nu_builders_match_the_scalar_loop(h_vacuum_2nu):
    r"""The same, for two flavors."""
    for name, build in builders_2nu(h_vacuum_2nu):
        batched = build(ENERGIES)
        scalar = np.stack([build(e) for e in ENERGIES])
        assert batched.shape == (len(ENERGIES), 2, 2), name
        assert np.array_equal(batched, scalar), name


def test_scalar_energy_still_returns_a_single_matrix(h_vacuum_3nu,
                                                     h_vacuum_2nu):
    r"""A scalar energy returns one matrix, not a stack of length one."""
    for name, build in builders_3nu(h_vacuum_3nu):
        assert build(1.0e9).shape == (3, 3), name
    for name, build in builders_2nu(h_vacuum_2nu):
        assert build(1.0e9).shape == (2, 2), name


def test_batched_hamiltonians_are_hermitian(h_vacuum_3nu, h_vacuum_2nu):
    r"""Every Hamiltonian in the stack is Hermitian."""
    for name, build in builders_3nu(h_vacuum_3nu) + \
            builders_2nu(h_vacuum_2nu):
        stack = build(ENERGIES)
        assert np.allclose(stack, np.conj(np.swapaxes(stack, -1, -2)),
                           atol=ATOL), name


def test_builders_do_not_mutate_the_vacuum_hamiltonian(h_vacuum_3nu):
    r"""Building a stack leaves the vacuum Hamiltonian passed in alone."""
    before = np.array(h_vacuum_3nu, dtype=complex).copy()
    for _, build in builders_3nu(h_vacuum_3nu):
        build(ENERGIES)
    assert np.array_equal(np.asarray(h_vacuum_3nu, dtype=complex), before)


def test_end_to_end_energy_scan(h_vacuum_3nu):
    r"""The whole scan --- build the stack, then evaluate it --- agrees
    with the equivalent double loop.

    This is the pattern the plotting modules use.
    """
    baseline = 1.3e3*CONV_KM_TO_INV_EV

    for name, build in builders_3nu(h_vacuum_3nu):
        batched = oscprob3nu.probabilities_3nu(build(ENERGIES), baseline)
        scalar = np.array([oscprob3nu.probabilities_3nu(build(e), baseline)
                           for e in ENERGIES])
        assert batched.shape == (len(ENERGIES), 9), name
        assert np.allclose(batched, scalar, atol=ATOL), name


def test_end_to_end_energy_scan_2nu(h_vacuum_2nu):
    r"""The same, for two flavors."""
    baseline = 1.3e3*CONV_KM_TO_INV_EV

    for name, build in builders_2nu(h_vacuum_2nu):
        batched = oscprob2nu.probabilities_2nu(build(ENERGIES), baseline)
        scalar = np.array([oscprob2nu.probabilities_2nu(build(e), baseline)
                           for e in ENERGIES])
        assert batched.shape == (len(ENERGIES), 4), name
        assert np.allclose(batched, scalar, atol=ATOL), name


def test_energy_and_baseline_grid(h_vacuum_3nu):
    r"""An energy stack and a baseline array combine into an
    oscillogram."""
    baselines = np.linspace(1.0, 1.3e4, 11)*CONV_KM_TO_INV_EV
    h_stack = hamiltonians3nu.hamiltonian_3nu_matter(
        h_vacuum_3nu, ENERGIES, VCC_EARTH_CRUST)

    grid = oscprob3nu.probabilities_3nu(h_stack[:, None, :, :],
                                        baselines[None, :])
    assert grid.shape == (len(ENERGIES), len(baselines), 9)

    for i in (0, 5, len(ENERGIES)-1):
        for j in (0, len(baselines)-1):
            expected = oscprob3nu.probabilities_3nu(h_stack[i], baselines[j])
            assert np.allclose(grid[i, j], expected, atol=ATOL)


def test_array_valued_matter_potential(h_vacuum_3nu):
    r"""The matter potential may vary alongside the energy, e.g. for a
    scan across a density profile."""
    potentials = np.linspace(0.5, 2.0, len(ENERGIES))*VCC_EARTH_CRUST

    batched = hamiltonians3nu.hamiltonian_3nu_matter(
        h_vacuum_3nu, ENERGIES, potentials)
    scalar = np.stack([hamiltonians3nu.hamiltonian_3nu_matter(
        h_vacuum_3nu, e, v) for e, v in zip(ENERGIES, potentials)])

    assert batched.shape == (len(ENERGIES), 3, 3)
    assert np.array_equal(batched, scalar)
