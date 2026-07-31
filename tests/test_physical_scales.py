# -*- coding: utf-8 -*-
r"""Tests at the scales the library is actually used at.

Everything else that compares the evaluation paths does so on random
Hermitian matrices with entries of order one and baselines of order ten.
Real Hamiltonians are nothing like that: their entries are around
10\ :sup:`-13` eV, baselines around 10\ :sup:`13` eV\ :sup:`-1`, and the
oscillation phase is the product of two numbers twenty-six orders of
magnitude apart.  Agreement in one regime does not imply agreement in
the other, so these tests use the bundled Hamiltonians at NuFit
parameters over the energies and baselines an experiment would.

They run on whichever backend is available, since the point is that the
compiled kernels and the NumPy path agree *here*, not only on tidy
inputs.
"""

import numpy as np
import pytest

import fastkernels
import hamiltonians2nu
import hamiltonians3nu
import oscprob2nu
import oscprob3nu
from globaldefs import (CONV_KM_TO_INV_EV, D21_NO_BF, D31_NO_BF, DCP_NO_BF,
                        EPS_2, EPS_3, LAMBDA, S12_NO_BF, S13_NO_BF, S23_NO_BF,
                        VCC_EARTH_CRUST)

from conftest import ATOL


# Energies from 10 MeV to 100 GeV, and baselines from a short-baseline
# reactor experiment to the diameter of the Earth and beyond
ENERGIES = np.logspace(-2.0, 2.0, 200)*1.0e9                    # [eV]
BASELINES = np.array([1.0, 300.0, 1.3e3, 1.3e4, 1.0e5])         # [km]


@pytest.fixture(scope='module')
def h_vacuum_3nu():
    r"""Returns the energy-independent 3nu vacuum Hamiltonian."""
    return hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
        S12_NO_BF, S23_NO_BF, S13_NO_BF, DCP_NO_BF, D21_NO_BF, D31_NO_BF)


@pytest.fixture(scope='module')
def h_vacuum_2nu():
    r"""Returns the energy-independent 2nu vacuum Hamiltonian."""
    return hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(
        S23_NO_BF, D31_NO_BF)


@pytest.fixture(params=['numpy', 'numba'])
def backend(request, monkeypatch):
    r"""Runs the test once per available backend."""
    if request.param == 'numba':
        if not fastkernels.HAVE_NUMBA:
            pytest.skip('Numba is not installed')
        monkeypatch.setattr(fastkernels, 'USE_NUMBA', True)
        # The two-flavor threshold keeps the kernel away from stacks this
        # size; lower it so the compiled path is what is being compared
        monkeypatch.setattr(fastkernels, 'MIN_BATCH', {2: 1, 3: 1})
    else:
        monkeypatch.setattr(fastkernels, 'USE_NUMBA', False)
    return request.param


def scenarios_3nu(h_vacuum):
    r"""Returns the bundled three-flavor Hamiltonians, per energy."""
    return [
        ('vacuum', lambda e: np.asarray(h_vacuum)/np.asarray(e)[..., None,
                                                               None]),
        ('matter', lambda e: hamiltonians3nu.hamiltonian_3nu_matter(
            h_vacuum, e, VCC_EARTH_CRUST)),
        ('nsi', lambda e: hamiltonians3nu.hamiltonian_3nu_nsi(
            h_vacuum, e, VCC_EARTH_CRUST, EPS_3)),
        ('liv', lambda e: hamiltonians3nu.hamiltonian_3nu_liv(
            h_vacuum, e, 0.3, 0.4, 0.5, 0.7, 1.0e-9, 1.5e-9, 2.0e-9,
            LAMBDA)),
    ]


@pytest.mark.parametrize('l_km', BASELINES)
def test_3nu_energy_scan_matches_the_scalar_path(backend, h_vacuum_3nu, l_km):
    r"""A realistic energy scan agrees with the scalar routine.

    The scalar path is the reference: it is the one checked against
    ``scipy.linalg.expm`` and against the standard oscillation formulas.
    """
    baseline = l_km*CONV_KM_TO_INV_EV
    for name, build in scenarios_3nu(h_vacuum_3nu):
        h_stack = build(ENERGIES)
        batched = oscprob3nu.probabilities_3nu(h_stack, baseline)
        scalar = np.array([oscprob3nu.probabilities_3nu(h.tolist(), baseline)
                           for h in h_stack])
        assert np.allclose(batched, scalar, atol=1.0e-12), name


@pytest.mark.parametrize('l_km', BASELINES)
def test_3nu_energy_scan_is_normalized(backend, h_vacuum_3nu, l_km):
    r"""Unitarity survives phases of order 10^5 radians."""
    baseline = l_km*CONV_KM_TO_INV_EV
    for name, build in scenarios_3nu(h_vacuum_3nu):
        prob = oscprob3nu.probabilities_3nu(build(ENERGIES),
                                            baseline).reshape(-1, 3, 3)
        assert np.allclose(prob.sum(axis=-1), 1.0, atol=1.0e-12), name
        assert np.allclose(prob.sum(axis=-2), 1.0, atol=1.0e-12), name
        assert np.all(prob >= -ATOL) and np.all(prob <= 1.0+ATOL), name


def test_3nu_baseline_scan_at_physical_scales(backend, h_vacuum_3nu):
    r"""A realistic baseline scan, one Hamiltonian over many distances."""
    h_matter = hamiltonians3nu.hamiltonian_3nu_matter(
        h_vacuum_3nu, 1.0e9, VCC_EARTH_CRUST)
    baselines = np.logspace(0.0, 5.0, 200)*CONV_KM_TO_INV_EV

    batched = oscprob3nu.probabilities_3nu(h_matter, baselines)
    scalar = np.array([oscprob3nu.probabilities_3nu(h_matter.tolist(), l)
                       for l in baselines])
    assert np.allclose(batched, scalar, atol=1.0e-12)


def test_2nu_energy_scan_at_physical_scales(backend, h_vacuum_2nu):
    r"""The two-flavor equivalent, in vacuum, matter, NSI and LIV."""
    baseline = 1.3e3*CONV_KM_TO_INV_EV
    builders = [
        ('matter', lambda e: hamiltonians2nu.hamiltonian_2nu_matter(
            h_vacuum_2nu, e, VCC_EARTH_CRUST)),
        ('nsi', lambda e: hamiltonians2nu.hamiltonian_2nu_nsi(
            h_vacuum_2nu, e, VCC_EARTH_CRUST, EPS_2)),
        ('liv', lambda e: hamiltonians2nu.hamiltonian_2nu_liv(
            h_vacuum_2nu, e, 0.3, 1.0e-9, 2.0e-9, LAMBDA)),
    ]
    for name, build in builders:
        h_stack = build(ENERGIES)
        batched = oscprob2nu.probabilities_2nu(h_stack, baseline)
        scalar = np.array([oscprob2nu.probabilities_2nu(h.tolist(), baseline)
                           for h in h_stack])
        assert np.allclose(batched, scalar, atol=ATOL), name


def test_oscillogram_at_physical_scales(backend, h_vacuum_3nu):
    r"""An energy-by-baseline grid, the shape an oscillogram takes."""
    energies = np.logspace(-1.0, 1.0, 40)*1.0e9
    baselines = np.logspace(1.0, 4.0, 30)*CONV_KM_TO_INV_EV
    h_stack = hamiltonians3nu.hamiltonian_3nu_matter(
        h_vacuum_3nu, energies, VCC_EARTH_CRUST)

    grid = oscprob3nu.probabilities_3nu(h_stack[:, None, :, :],
                                        baselines[None, :])
    assert grid.shape == (40, 30, 9)
    assert np.allclose(grid.reshape(-1, 3, 3).sum(axis=-1), 1.0, atol=1.0e-12)

    for i in (0, 20, 39):
        for j in (0, 29):
            expected = oscprob3nu.probabilities_3nu(h_stack[i].tolist(),
                                                    baselines[j])
            assert np.allclose(grid[i, j], expected, atol=1.0e-12)


def test_very_long_baselines(backend, h_vacuum_3nu):
    r"""Phases of order 10^7 radians still give sane probabilities.

    Astrophysical baselines wind the phase through millions of cycles.
    Nothing about the closed form breaks there, but the accumulated
    round-off in the argument of the exponential is worth pinning.
    """
    h_matter = hamiltonians3nu.hamiltonian_3nu_matter(
        h_vacuum_3nu, 1.0e9, VCC_EARTH_CRUST)
    baselines = np.array([1.0e5, 1.0e6, 1.0e7])*CONV_KM_TO_INV_EV

    prob = oscprob3nu.probabilities_3nu(h_matter, baselines).reshape(-1, 3, 3)
    assert np.all(np.isfinite(prob))
    assert np.allclose(prob.sum(axis=-1), 1.0, atol=1.0e-11)
    assert np.all(prob >= -ATOL) and np.all(prob <= 1.0+ATOL)


def test_near_degenerate_at_scale(backend):
    r"""A near-degenerate spectrum, split by one part in 10^14.

    This sits just on the ordinary side of the degeneracy threshold, so
    it takes the general expression with denominators that are nearly
    zero --- the numerically worst case the expansion has.
    """
    h_stack = np.stack([np.diag([1.0, 1.0 + 1.0e-14, -2.0]).astype(complex)]
                       * 40)
    prob = oscprob3nu.probabilities_3nu(h_stack, 3.0)

    assert np.all(np.isfinite(prob))
    assert np.allclose(prob.reshape(-1, 3, 3).sum(axis=-1), 1.0, atol=1.0e-10)

    scalar = np.array(oscprob3nu.probabilities_3nu(h_stack[0].tolist(), 3.0))
    assert np.allclose(prob[0], scalar, atol=1.0e-10)
