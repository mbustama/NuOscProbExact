# -*- coding: utf-8 -*-
r"""Checks the probabilities against frozen nuSQuIDS reference values.

Everything else in this suite is, in the end, self-referential about
*conventions*.  ``scipy.linalg.expm`` confirms that :math:`U` is the
matrix exponential of the Hamiltonian --- but of the Hamiltonian this
library built.  The standard oscillation formulas are transcribed from
the same papers the builders were.  If a mixing-matrix ordering, a sign,
or a unit conversion were wrong in a way that is consistent across the
package, none of those tests would notice.

`nuSQuIDS <https://github.com/arguelles/nuSQuIDS>`_ would.  It is an
independent C++ code that evolves the density matrix numerically, and it
was configured here from mixing *angles* and a density in g/cm^3 rather
than from a Hamiltonian, so agreement exercises
:mod:`hamiltonians3nu` and :mod:`hamiltonians4nu` as much as the
expansions themselves.

The reference values are frozen into ``nusquids_reference.json`` rather
than computed here.  nuSQuIDS ships manylinux wheels and would install
in CI, but a comparison that runs in CI is a dependency on somebody
else's release cadence, and this repository has already written down why
it declines those --- see the ruff rule pinning in ``pyproject.toml``.
``nusquids_reference.py``, beside the data, regenerates it deliberately.

Two conventions are matched before comparing, and matching them is the
point rather than a fudge.  Both are differences in definition:

* ``globaldefs.CONV_KM_TO_INV_EV`` carries six significant figures where
  nuSQuIDS uses :math:`1\,\mathrm{km}/\hbar c` in full.  The relative
  difference is 1.4e-7, which the accumulated phase turns into 5e-7 in
  three-flavor probabilities and 1e-4 in four-flavor ones.
* The electron density is one gram over the mean nucleon mass here, and
  :math:`\rho N_A Y_e` in nuSQuIDS --- 0.82% apart, the mass defect.

Unmatched, the two codes appear to disagree at 1e-4.  Matched, they
agree to between 1e-11 and 1e-15, which is the real answer and the
reason both offsets are pinned by tests below rather than left as
folklore.
"""

import json
import os

import numpy as np
import pytest

import globaldefs as gd
import hamiltonians3nu
import hamiltonians4nu
import oscprob3nu
import oscprob4nu


HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE = os.path.join(HERE, 'nusquids_reference.json')

# Matched against nuSQuIDS in vacuum and in matter respectively.  The
# four-flavor cases are looser because a stiff 3+1 spectrum is, and by
# exactly the amount `oscprob4nu.POLISH_ROOTS` documents.
TOLERANCE = {3: 1.0e-11, 4: 1.0e-10}


def reference_data():
    r"""Returns the frozen nuSQuIDS reference values."""
    with open(REFERENCE) as handle:
        return json.load(handle)


DATA = reference_data()
CASES = DATA['cases']
IDS = [case['name'] for case in CASES]


def matched_km():
    r"""Returns nuSQuIDS' km-to-inverse-eV factor."""
    return DATA['constants']['km']


def matched_potentials(density):
    r"""Returns ``(VCC, VNC)`` in nuSQuIDS' convention, in eV.

    Electrons per cm^3 as :math:`\rho N_A Y_e`, converted to eV^3 by the
    length unit, rather than one gram over the mean nucleon mass.
    """
    constants = DATA['constants']
    fraction = DATA['electron_fraction']

    electrons = density*constants['avogadro']*fraction/constants['cm']**3
    neutrons = density*constants['avogadro']*(1.0-fraction)/constants['cm']**3

    return (np.sqrt(2.0)*constants['GF']*electrons,
            -constants['GF']*neutrons/np.sqrt(2.0))


def our_probabilities(case, km, potentials):
    r"""Returns this library's probabilities for one reference case."""
    parameters = DATA['parameters']
    sine = np.sin
    energy = case['energy_gev']*1.0e9
    baseline = case['baseline_km']*km
    density = case['density_g_cm3']

    if case['n_flavors'] == 3:
        vacuum = np.asarray(
            hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
                sine(parameters['th12']), sine(parameters['th23']),
                sine(parameters['th13']), parameters['dcp'],
                parameters['dm21'], parameters['dm31']))
        hamiltonian = (vacuum/energy if density is None
                       else hamiltonians3nu.hamiltonian_3nu_matter(
                           vacuum, energy, potentials[0]))
        return np.asarray(oscprob3nu.probabilities_3nu(
            hamiltonian, baseline)).reshape(3, 3)

    vacuum = hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(
        sine(parameters['th12']), sine(parameters['th23']),
        sine(parameters['th13']), sine(parameters['th14']),
        sine(parameters['th24']), sine(parameters['th34']),
        parameters['dcp'], parameters['dm21'], parameters['dm31'],
        parameters['dm41'])
    hamiltonian = (vacuum/energy if density is None
                   else hamiltonians4nu.hamiltonian_4nu_matter(
                       vacuum, energy, potentials[0], potentials[1]))

    return np.asarray(oscprob4nu.probabilities_4nu(
        hamiltonian, baseline)).reshape(4, 4)


@pytest.mark.parametrize('case', CASES, ids=IDS)
def test_agrees_with_nusquids(case):
    r"""Every frozen case agrees once the conventions are matched."""
    density = case['density_g_cm3']
    potentials = ((None, None) if density is None
                  else matched_potentials(density))

    ours = our_probabilities(case, matched_km(), potentials)
    theirs = np.asarray(case['probabilities'])

    deviation = np.max(np.abs(ours - theirs))
    assert deviation < TOLERANCE[case['n_flavors']], (
        '%s: NuOscProbExact and nuSQuIDS differ by %.2e' %
        (case['name'], deviation))


@pytest.mark.parametrize('case', CASES, ids=IDS)
def test_reference_probabilities_are_normalized(case):
    r"""The frozen values are themselves a sane probability table."""
    probabilities = np.asarray(case['probabilities'])

    assert probabilities.shape == (case['n_flavors'], case['n_flavors'])
    assert np.all(probabilities >= -1.0e-12)
    assert np.max(np.abs(probabilities.sum(axis=1) - 1.0)) < 1.0e-9


def test_the_length_unit_is_the_only_vacuum_difference():
    r"""Pins the six-figure ``CONV_KM_TO_INV_EV`` and what it costs.

    Not a defect --- the constant is documented as it stands, and 1e-7 is
    far below any measurable effect --- but it is the whole of the
    apparent three-flavor disagreement, and it grows with the phase, so
    it is worth having written down and guarded.
    """
    exact = matched_km()

    assert abs(gd.CONV_KM_TO_INV_EV/exact - 1.0) < 3.0e-7
    assert abs(gd.CONV_KM_TO_INV_EV/exact - 1.0) > 1.0e-8

    case = next(c for c in CASES if c['name'].startswith('3nu vacuum, 1 GeV'))
    theirs = np.asarray(case['probabilities'])

    matched = np.max(np.abs(
        our_probabilities(case, exact, (None, None)) - theirs))
    unmatched = np.max(np.abs(
        our_probabilities(case, gd.CONV_KM_TO_INV_EV, (None, None)) - theirs))

    assert matched < 1.0e-13
    assert unmatched > 100.0*matched


def test_the_matter_potential_convention_differs_by_the_mass_defect():
    r"""Pins the 0.82% between one gram over m_N and rho N_A Y_e."""
    ours = gd.VCC_EARTH_CRUST
    theirs, _ = matched_potentials(gd.DENSITY_MATTER_CRUST_G_PER_CM3)

    ratio = theirs/ours
    assert 1.005 < ratio < 1.012, (
        'the two electron-density conventions are %.3f%% apart, which is '
        'not the nuclear mass defect this test expects' % (100*(ratio-1)))
