"""The exact reference for the constant-density comparison figure.

That figure had been plotting every other code's residual against
NuOscProbExact, which makes this library its own referee.  It is refereed
here the same way the speed-accuracy plane is: by exponentiating the same
Hamiltonian in fifty-digit arithmetic.

Emits tests/const_density_scan.json, on the same 150-energy grid as
tests/nusquids_scan.json and tests/nufast_scan.json, and folds in the
GLoBES and Prob3++ scans produced by the two drivers beside it.
"""

import json
import os
import sys

import numpy as np
import mpmath as mp

sys.path.insert(0, 'src')

import globaldefs as gd
import earth
import hamiltonians3nu

BASELINE_KM = 1300.0
DENSITY = 3.0
E_GEV = np.logspace(np.log10(0.6), np.log10(20.0), 150)

H_VAC = np.asarray(hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
    gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF,
    gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF))


def mp_reference(energy_ev, dps=50):
    r"""P(nu_mu -> nu_e) from a `dps`-digit matrix exponential."""
    mp.mp.dps = dps
    h = np.asarray(hamiltonians3nu.hamiltonian_3nu_matter(
        H_VAC, energy_ev, earth.matter_potential(DENSITY)))
    m = mp.matrix(3, 3)
    for i in range(3):
        for j in range(3):
            m[i, j] = mp.mpc(complex(h[i, j]))
    u = mp.expm(-1j*m*(mp.mpf(BASELINE_KM)*mp.mpf(gd.CONV_KM_TO_INV_EV)))
    return float(abs(u[0, 1])**2.0)


def read(path):
    return [float(line.split()[1]) for line in open(path) if line.strip()]


if __name__ == '__main__':
    # Where globes_scan.c and prob3_scan.cpp leave their tables.  This was
    # an absolute path into one session's scratch directory, which made the
    # dataset unreproducible by anyone else; it now defaults beside this
    # script and is overridden with NUOSC_CONST_SCAN_DIR.
    S = os.environ.get(
        'NUOSC_CONST_SCAN_DIR',
        os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scans')) + os.sep
    for name in ('globes_scan.txt', 'prob3_scan.txt'):
        if not os.path.exists(S + name):
            raise SystemExit(
                '%s not found.  Build and run globes_scan.c and '
                'prob3_scan.cpp as the README beside this file describes, '
                'then point NUOSC_CONST_SCAN_DIR at their output.'
                % (S + name))
    payload = {
        'generated_by': 'tests/external_drivers/const_scan.py',
        'note': (
            'The exact reference for the constant-density comparison '
            'figure, and the two codes that had not been on it.  '
            'P(nu_mu -> nu_e) at L = 1300 km through 3 g/cm^3, on the same '
            '150-energy grid as nusquids_scan.json and nufast_scan.json.  '
            'The reference is a 50-digit mpmath matrix exponential of the '
            'same Hamiltonian, so that no code in that figure is its own '
            'referee.  GLoBES and Prob3++ are at their standard settings, '
            'each with the density and baseline conventions matched as '
            'documented in the README beside this file.'),
        'channel': 'P(nu_mu -> nu_e)',
        'baseline_km': BASELINE_KM,
        'density_g_cm3': DENSITY,
        'mpmath_dps': 50,
        'energy_gev': E_GEV.tolist(),
        'probability_reference': [mp_reference(e*1.0e9) for e in E_GEV],
        'probability_globes': read(S+'globes_scan.txt'),
        'probability_prob3pp': read(S+'prob3_scan.txt'),
    }
    print(json.dumps(payload, indent=1))
