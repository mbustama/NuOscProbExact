"""Extend the frozen nuSQuIDS timing to the sizes the other codes reach.

Same problem and same methodology as the rows already in
tests/timing_other_codes.json: P(nu_mu -> nu_e) at L = 1300 km through
constant matter of 3 g/cm^3, three flavors, multiple-energy mode, rel and
abs error 1e-12, timed through the Python interface including building the
solver and reading the results back, best of three.
"""

import json
import sys
import time

import numpy as np

sys.path.insert(0, '/home/mbustamante/Research/NuOscProb/NuOscProbExact/src')

import globaldefs as gd
import nuSQuIDS as nsq

U = nsq.Const()
AVOGADRO = 6.02214076e23
TH12 = np.arcsin(np.sqrt(0.310))
TH13 = np.arcsin(np.sqrt(2.240e-2))
TH23 = np.arcsin(np.sqrt(0.582))
DCP = 217.0/180.0*np.pi
_EPG = gd.CONV_G_TO_EV/((gd.MASS_PROTON+gd.MASS_NEUTRON)/2.0)
RHO = 3.0*(_EPG/AVOGADRO)


def scan(n):
    energies = np.logspace(-1.0, 1.5, n)*U.GeV
    if n == 1:
        s = nsq.nuSQUIDS(3, nsq.NeutrinoType.neutrino)
        s.Set_E(energies[0])
        state = np.array([0.0, 1.0, 0.0])
    else:
        s = nsq.nuSQUIDS(energies, 3, nsq.NeutrinoType.neutrino, False)
        state = np.tile(np.array([0.0, 1.0, 0.0]), (n, 1))
    s.Set_MixingAngle(0, 1, TH12)
    s.Set_MixingAngle(0, 2, TH13)
    s.Set_MixingAngle(1, 2, TH23)
    s.Set_SquareMassDifference(1, 7.39e-5)
    s.Set_SquareMassDifference(2, 2.525e-3)
    s.Set_CPPhase(0, 2, DCP)
    s.Set_rel_error(1.0e-12)
    s.Set_abs_error(1.0e-12)
    s.Set_Body(nsq.ConstantDensity(RHO, 0.5))
    s.Set_Track(nsq.ConstantDensity.Track(1300.0*U.km))
    s.Set_initial_state(state, nsq.Basis.flavor)
    s.EvolveState()
    if n == 1:
        return [s.EvalFlavor(0)]
    return [s.EvalFlavor(0, e) for e in energies]


def best_of(n, repeat=3):
    best = np.inf
    for _ in range(repeat):
        t0 = time.perf_counter()
        scan(n)
        best = min(best, time.perf_counter()-t0)
    return best


if __name__ == '__main__':
    out = {}
    for n in (10000, 30000):
        out[n] = best_of(n)
        print('%8d  %.6g s   (%.3f us/probability)'
              % (n, out[n], out[n]/n*1e6), flush=True)
    print(json.dumps({str(k): v for k, v in out.items()}))
