# -*- coding: utf-8 -*-
r"""Regenerates the frozen nuSQuIDS energy scan used by the paper's figures.

Companion to :mod:`tests.nusquids_reference`, which freezes single points
for the test suite.  This one freezes a *curve*: :math:`P_{\nu_\mu \to
\nu_e}` against energy at one baseline, so that a figure can show an
independently written code beside the closed form without the notebook
having to build or import nuSQuIDS.

Run it, from the repository root, in an environment that has nuSQuIDS::

    python tests/nusquids_scan.py > tests/nusquids_scan.json

One convention has to be matched, and it is the same one documented in
:mod:`tests.nusquids_reference`.  This library takes the electron density
as :math:`Y_e \rho / \bar{m}`, with :math:`\bar{m}` the mean free-nucleon
mass; nuSQuIDS takes :math:`\rho N_A Y_e`.  The two differ by the nuclear
mass defect.  Handing nuSQuIDS the rescaled density below makes both
codes propagate with the same :math:`V_{\rm CC}`, so that what the figure
compares is the oscillation physics rather than a choice of nuclear mass.
"""

import json
import sys

import numpy as np

import nuSQuIDS as nsq

sys.path.insert(0, 'src')
import globaldefs as gd                                       # noqa: E402

UNITS = nsq.Const()
AVOGADRO = 6.02214076e23

# NuFit 4.0, normal ordering, as `globaldefs` carries it.  Angles, not
# sines: nuSQuIDS takes the angle, this library takes its sine.
TH12 = np.arcsin(np.sqrt(0.310))
TH13 = np.arcsin(np.sqrt(2.240e-2))
TH23 = np.arcsin(np.sqrt(0.582))
DCP = 217.0/180.0*np.pi
DM21 = 7.39e-5
DM31 = 2.525e-3

BASELINE_KM = 1300.0
DENSITY_G_CM3 = 3.0
N_POINTS = 150
E_MIN_GEV, E_MAX_GEV = 0.6, 20.0

# Electrons per gram: ours is one over the mean nucleon mass, theirs is
# Avogadro's number.  The ratio is the mass defect, and is the factor by
# which the density handed to nuSQuIDS must be scaled.
_ELECTRONS_PER_GRAM = gd.CONV_G_TO_EV/((gd.MASS_PROTON+gd.MASS_NEUTRON)/2.0)
DENSITY_SCALE = _ELECTRONS_PER_GRAM/AVOGADRO
DENSITY_FOR_NUSQUIDS = DENSITY_G_CM3*DENSITY_SCALE


def p_mue(energy_gev, density_g_cm3):
    r"""Returns :math:`P_{\nu_\mu \to \nu_e}` from nuSQuIDS."""
    solver = nsq.nuSQUIDS(3, nsq.NeutrinoType.neutrino)
    solver.Set_MixingAngle(0, 1, TH12)
    solver.Set_MixingAngle(0, 2, TH13)
    solver.Set_MixingAngle(1, 2, TH23)
    solver.Set_SquareMassDifference(1, DM21)
    solver.Set_SquareMassDifference(2, DM31)
    solver.Set_CPPhase(0, 2, DCP)
    solver.Set_rel_error(1.0e-12)
    solver.Set_abs_error(1.0e-12)
    solver.Set_Body(nsq.ConstantDensity(density_g_cm3, 0.5))
    solver.Set_Track(nsq.ConstantDensity.Track(BASELINE_KM*UNITS.km))
    solver.Set_E(energy_gev*UNITS.GeV)
    solver.Set_initial_state(np.array([0.0, 1.0, 0.0]), nsq.Basis.flavor)
    solver.EvolveState()
    return solver.EvalFlavor(0)


def main():
    energies = np.logspace(np.log10(E_MIN_GEV), np.log10(E_MAX_GEV),
                           N_POINTS)
    print(json.dumps({
        'generated_by': 'tests/nusquids_scan.py',
        'nusquids_version': '1.13.3',
        'channel': 'P(nu_mu -> nu_e)',
        'baseline_km': BASELINE_KM,
        'density_g_cm3': DENSITY_G_CM3,
        'density_handed_to_nusquids': DENSITY_FOR_NUSQUIDS,
        'density_scale': DENSITY_SCALE,
        'energy_gev': energies.tolist(),
        'probability': [p_mue(e, DENSITY_FOR_NUSQUIDS) for e in energies],
    }, indent=1))


if __name__ == '__main__':
    main()
