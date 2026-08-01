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

# The four-flavor tolerance is looser because a stiff 3+1 spectrum is, by
# exactly the amount `oscprob4nu.POLISH_ROOTS` documents.  It is set at
# 1e-9 rather than tighter because of the antineutrino case, which is the
# hardest here: with the potential reversed the eigenvalue cluster tightens
# further, and the closed form lands at 3e-10.
#
# That is *our* limit rather than a disagreement, and
# `test_the_four_flavor_residual_is_ours_not_a_convention` below proves it
# instead of asserting it: against `scipy.linalg.expm` on the same
# Hamiltonian we are 3.0e-10 out and nuSQuIDS is 3.4e-12, so the whole
# residual is on this side and none of it is a mismatched convention.
TOLERANCE = {3: 1.0e-11, 4: 1.0e-9}


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
    r"""Returns this library's probabilities for one reference case.

    Antineutrinos need **both** changes, which is the point of including
    them: conjugate the vacuum Hamiltonian *and* reverse the sign of the
    potentials.  Doing only one is silently plausible, and is the single
    most common way to get a matter calculation wrong.
    """
    parameters = DATA['parameters']
    sine = np.sin
    energy = case['energy_gev']*1.0e9
    baseline = case['baseline_km']*km
    density = case['density_g_cm3']
    antineutrino = case.get('antineutrino', False)
    dm31 = case.get('dm31', parameters['dm31'])

    sign = -1.0 if antineutrino else 1.0
    VCC = None if density is None else sign*potentials[0]
    VNC = None if density is None else sign*potentials[1]

    if case['n_flavors'] == 3:
        vacuum = np.asarray(
            hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
                sine(parameters['th12']), sine(parameters['th23']),
                sine(parameters['th13']), parameters['dcp'],
                parameters['dm21'], dm31))
        if antineutrino:
            vacuum = np.conj(vacuum)
        hamiltonian = (vacuum/energy if density is None
                       else hamiltonians3nu.hamiltonian_3nu_matter(
                           vacuum, energy, VCC))
        return np.asarray(oscprob3nu.probabilities_3nu(
            hamiltonian, baseline)).reshape(3, 3)

    vacuum = hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(
        sine(parameters['th12']), sine(parameters['th23']),
        sine(parameters['th13']), sine(parameters['th14']),
        sine(parameters['th24']), sine(parameters['th34']),
        parameters['dcp'], parameters['dm21'], dm31,
        parameters['dm41'])
    if antineutrino:
        vacuum = np.conj(vacuum)
    hamiltonian = (vacuum/energy if density is None
                   else hamiltonians4nu.hamiltonian_4nu_matter(
                       vacuum, energy, VCC, VNC))

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


@pytest.mark.parametrize(
    'case', [c for c in CASES if c['n_flavors'] == 4], ids=[
        c['name'] for c in CASES if c['n_flavors'] == 4])
def test_the_four_flavor_residual_is_ours_not_a_convention(case):
    r"""Attributes the four-flavor residual, rather than tolerating it.

    The four-flavor cases agree with nuSQuIDS less well than the
    three-flavor ones, and a loose tolerance alone would leave it
    ambiguous whether that is our stiff-spectrum accuracy or a
    convention we failed to match.

    It is the former, and the way to show it is a third party.  Against
    ``scipy.linalg.expm`` on the very Hamiltonian we hand to the
    expansion, our error and our disagreement with nuSQuIDS are the same
    size --- so nuSQuIDS is simply where ``expm`` is, and there is no
    convention gap hiding underneath.  If a mismatched convention ever
    did creep in, the two would come apart and this fails.
    """
    scipy_linalg = pytest.importorskip('scipy.linalg')

    density = case['density_g_cm3']
    potentials = ((None, None) if density is None
                  else matched_potentials(density))
    ours = our_probabilities(case, matched_km(), potentials)
    theirs = np.asarray(case['probabilities'])

    parameters = DATA['parameters']
    sine = np.sin
    energy = case['energy_gev']*1.0e9
    baseline = case['baseline_km']*matched_km()
    sign = -1.0 if case.get('antineutrino', False) else 1.0

    vacuum = hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(
        sine(parameters['th12']), sine(parameters['th23']),
        sine(parameters['th13']), sine(parameters['th14']),
        sine(parameters['th24']), sine(parameters['th34']),
        parameters['dcp'], parameters['dm21'],
        case.get('dm31', parameters['dm31']), parameters['dm41'])
    if case.get('antineutrino', False):
        vacuum = np.conj(vacuum)

    hamiltonian = (vacuum/energy if density is None
                   else hamiltonians4nu.hamiltonian_4nu_matter(
                       vacuum, energy, sign*potentials[0],
                       sign*potentials[1]))
    hamiltonian = np.asarray(hamiltonian)
    traceless = hamiltonian - np.trace(hamiltonian).real/4.0*np.eye(4)
    exponential = np.abs(
        scipy_linalg.expm(-1.j*traceless*baseline).T)**2

    against_expm = np.max(np.abs(ours - exponential))
    against_nusquids = np.max(np.abs(ours - theirs))

    # The two agree in size, so the residual is entirely ours.  The
    # window is generous because both are tiny; what would fail here is a
    # convention gap, which would push them apart by orders of magnitude.
    assert against_nusquids < 10.0*against_expm + 1.0e-15


def test_the_matter_potential_convention_differs_by_the_mass_defect():
    r"""Pins the 0.82% between one gram over m_N and rho N_A Y_e."""
    ours = gd.VCC_EARTH_CRUST
    theirs, _ = matched_potentials(gd.DENSITY_MATTER_CRUST_G_PER_CM3)

    ratio = theirs/ours
    assert 1.005 < ratio < 1.012, (
        'the two electron-density conventions are %.3f%% apart, which is '
        'not the nuclear mass defect this test expects' % (100*(ratio-1)))
