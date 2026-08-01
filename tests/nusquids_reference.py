# -*- coding: utf-8 -*-
r"""Regenerates the frozen nuSQuIDS reference probabilities.

This is **not** a test and pytest does not collect it.  It is the script
that produced ``tests/nusquids_reference.json``, kept beside the data so
that the numbers can be reproduced and extended rather than trusted.

`nuSQuIDS <https://github.com/arguelles/nuSQuIDS>`_ is an independent
code: it evolves the neutrino density matrix numerically, in C++, and
shares no lineage with the closed-form SU(N) expansions here.  Agreement
with it therefore tests something the internal suite cannot --- the
*conventions*, not only the algebra.  ``scipy.linalg.expm`` checks that
:math:`U` is the matrix exponential of the Hamiltonian we built; nothing
internal checks that we built the right Hamiltonian.

Run it, from the repository root, in an environment that has nuSQuIDS::

    python -m venv /tmp/nsq && /tmp/nsq/bin/pip install nusquids
    /tmp/nsq/bin/python tests/nusquids_reference.py > tests/nusquids_reference.json

nuSQuIDS ships manylinux wheels, so that needs no compiler.  It is
deliberately kept out of the test extras: the comparison is frozen into
the JSON rather than run in CI, so that a release of somebody else's
code cannot turn this repository red.

Two conventions have to be matched, and they are the whole reason this
file records constants alongside probabilities.  Both are differences in
*definition*, not errors in either code:

* **The length unit.**  ``globaldefs.CONV_KM_TO_INV_EV`` is ``5.06773e9``,
  six significant figures; nuSQuIDS uses ``5.0677307162e9``, which is
  :math:`1\,\mathrm{km}/\hbar c`.  The relative difference is
  :math:`1.4\times10^{-7}`, which the accumulated phase amplifies --- to
  :math:`5\times10^{-7}` in the probabilities at three flavors, and to
  :math:`10^{-4}` at four, where an eV-scale sterile drives the phase to
  some :math:`10^4` radian.
* **The matter potential.**  We take the electron density as one gram
  divided by the mean nucleon mass; nuSQuIDS takes
  :math:`\rho N_A Y_e`.  The two differ by 0.82%, which is the nuclear
  mass defect.

Matched on both, the two codes agree to between :math:`10^{-11}` and
:math:`10^{-15}` on every case here.
"""

import json
import sys

import numpy as np

import nuSQuIDS as nsq                                        # noqa: E402


UNITS = nsq.Const()

# NuFit 4.0, normal ordering, as `globaldefs` carries it.  Angles, not
# sines: nuSQuIDS takes the angle, this library takes its sine.
TH12 = np.arcsin(np.sqrt(0.310))
TH13 = np.arcsin(np.sqrt(2.240e-2))
TH23 = np.arcsin(np.sqrt(0.582))
DCP = 217.0/180.0*np.pi
DM21 = 7.39e-5
DM31 = 2.525e-3

# The 3+1 additions, matching notebook 16.
DM41 = 1.0
TH14 = np.arcsin(np.sqrt(0.10))
TH24 = np.arcsin(np.sqrt(0.10))
TH34 = 0.0

CASES = [
    ('3nu vacuum, 1 GeV, 1300 km', 3, 1.0, 1300.0, None),
    ('3nu vacuum, 0.6 GeV, 295 km', 3, 0.6, 295.0, None),
    ('3nu matter rho=3, 1 GeV, 1300 km', 3, 1.0, 1300.0, 3.0),
    ('3nu matter rho=3, 5 GeV, 1300 km', 3, 5.0, 1300.0, 3.0),
    ('4nu vacuum, 1 GeV, 1300 km', 4, 1.0, 1300.0, None),
    ('4nu vacuum, 0.5 GeV, 0.6 km', 4, 0.5, 0.6, None),
    ('4nu matter rho=3, 1 GeV, 1300 km', 4, 1.0, 1300.0, 3.0),
]

ELECTRON_FRACTION = 0.5


def configured(n_flavors, body, track):
    r"""Returns a nuSQuIDS instance with the mixing parameters set."""
    solver = nsq.nuSQUIDS(n_flavors, nsq.NeutrinoType.neutrino)
    solver.Set_MixingAngle(0, 1, TH12)
    solver.Set_MixingAngle(0, 2, TH13)
    solver.Set_MixingAngle(1, 2, TH23)
    solver.Set_SquareMassDifference(1, DM21)
    solver.Set_SquareMassDifference(2, DM31)
    solver.Set_CPPhase(0, 2, DCP)

    if n_flavors == 4:
        solver.Set_MixingAngle(0, 3, TH14)
        solver.Set_MixingAngle(1, 3, TH24)
        solver.Set_MixingAngle(2, 3, TH34)
        solver.Set_SquareMassDifference(3, DM41)

    # Far tighter than needed.  Loosening these to 1e-5 changes nothing,
    # which is how we know the residual is convention rather than solver
    # tolerance.
    solver.Set_rel_error(1.0e-12)
    solver.Set_abs_error(1.0e-12)
    solver.Set_Body(body)
    solver.Set_Track(track)

    return solver


def probabilities(n_flavors, energy_gev, baseline_km, density):
    r"""Returns ``P[alpha][beta]``, the initial flavor varying slowest."""
    if density is None:
        body = nsq.Vacuum()
        track = nsq.Vacuum.Track(baseline_km*UNITS.km)
    else:
        body = nsq.ConstantDensity(density, ELECTRON_FRACTION)
        track = nsq.ConstantDensity.Track(baseline_km*UNITS.km)

    rows = []
    for alpha in range(n_flavors):
        solver = configured(n_flavors, body, track)
        solver.Set_E(energy_gev*UNITS.GeV)
        initial = np.zeros(n_flavors)
        initial[alpha] = 1.0
        solver.Set_initial_state(initial, nsq.Basis.flavor)
        solver.EvolveState()
        rows.append([solver.EvalFlavor(beta) for beta in range(n_flavors)])

    return rows


def main():
    r"""Writes the reference data as JSON on standard output."""
    reference = {
        'generated_by': 'tests/nusquids_reference.py',
        'nusquids_version': '1.13.3',
        'numpy_version': np.__version__,
        'python_version': sys.version.split()[0],
        'electron_fraction': ELECTRON_FRACTION,
        'constants': {
            'km': UNITS.km,
            'GeV': UNITS.GeV,
            'GF': UNITS.GF,
            'cm': UNITS.cm,
            'avogadro': 6.02214076e23,
        },
        'parameters': {
            'th12': TH12, 'th13': TH13, 'th23': TH23, 'dcp': DCP,
            'dm21': DM21, 'dm31': DM31,
            'th14': TH14, 'th24': TH24, 'th34': TH34, 'dm41': DM41,
        },
        'cases': [],
    }

    for name, n_flavors, energy, baseline, density in CASES:
        reference['cases'].append({
            'name': name,
            'n_flavors': n_flavors,
            'energy_gev': energy,
            'baseline_km': baseline,
            'density_g_cm3': density,
            'probabilities': probabilities(n_flavors, energy, baseline,
                                           density),
        })
        print('generated: %s' % name, file=sys.stderr)

    print(json.dumps(reference, indent=1))


if __name__ == '__main__':
    main()
