"""Generate and execute the NuOscProbExact example notebooks.

Run from the repository root::

    python notebooks/make_notebooks.py

The notebooks are generated from this one script rather than edited by hand,
so that all fifteen share a setup cell, a plotting style and a voice, and so
that regenerating them after an API change is a single command.

**Edit this file, not the notebooks.**  Anything written directly into a
``.ipynb`` here is lost the next time this runs.

Generation and execution are deliberately one step.  Generating rewrites every
notebook, discarding its stored outputs; executing puts them back.  Doing only
the first leaves fifteen notebooks that still open and still run but show a
reader nothing, which is a silent failure --- it happened twice during the
work that produced them, caught both times by the guard in
``.github/workflows/lint.yml`` rather than locally.  The same guard runs at the
end of this script so that it cannot happen again without being noticed here
first.
"""

import pathlib
import time

import nbformat as nbf

OUT = pathlib.Path('notebooks')

SETUP = '''\
import sys
import os

# Works whether or not the package is installed: these notebooks sit one
# level below the repository root, so src/ is one directory up.
sys.path.insert(0, os.path.abspath(os.path.join('..', 'src')))

import numpy as np
import matplotlib.pyplot as plt

import globaldefs as gd
import oscprob2nu
import oscprob3nu
import oscprob4nu
import hamiltonians2nu
import hamiltonians3nu
import hamiltonians4nu

%matplotlib inline
plt.rcParams.update({'figure.figsize': (7.2, 4.2), 'figure.dpi': 90,
                     'axes.grid': True, 'grid.alpha': 0.3,
                     'font.size': 11, 'legend.frameon': False,
                     # Axes end exactly at the data: no padding beyond
                     # x_min/x_max or y_min/y_max.
                     'axes.xmargin': 0.0, 'axes.ymargin': 0.0})

# The three-flavor vacuum Hamiltonian at the NuFit best fit, normal ordering.
# It is energy-independent; dividing by the energy is what the builders do.
H_VAC_3NU = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
    gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, gd.DCP_NO_BF,
    gd.D21_NO_BF, gd.D31_NO_BF)

KM = gd.CONV_KM_TO_INV_EV      # multiply a length in km to get eV^-1
GEV = 1.0e9                    # multiply an energy in GeV to get eV
'''


def md(text):
    return nbf.v4.new_markdown_cell(text.rstrip())


def code(text):
    return nbf.v4.new_code_cell(text.rstrip())


def notebook(title, intro, cells):
    nb = nbf.v4.new_notebook()
    nb.cells = [md('# %s\n\n%s' % (title, intro)), code(SETUP)] + cells
    nb.metadata = {
        'kernelspec': {'display_name': 'Python 3', 'language': 'python',
                       'name': 'python3'},
        'language_info': {'name': 'python', 'pygments_lexer': 'ipython3'},
    }
    return nb


books = {}

# ---------------------------------------------------------------- 01 basics
books['01_basics.ipynb'] = notebook(
    'Basics',
    'The shortest path to a probability, and the conventions everything else '
    'in these notebooks assumes.\n\n'
    '**NuOscProbExact** computes oscillation probabilities *exactly*, with no '
    'approximation beyond floating-point round-off, for any time-independent '
    'Hamiltonian. The method expands the Hamiltonian and the evolution '
    'operator in the SU(2), SU(3) and SU(4) bases, which gives closed forms '
    'rather than a numerical integration.',
    [
        md('## Units\n\n'
           'The library is unit-agnostic: it asks only that the Hamiltonian '
           'and the baseline be given in reciprocal units, so that $HL$ is '
           'dimensionless. Everything below uses **eV** for energies and '
           '**eV$^{-1}$** for baselines, with `globaldefs` supplying the '
           'conversions.'),
        code('print('
             '"1 km          = %.4e eV^-1" % KM)\n'
             'print("1 GeV         = %.0e eV" % GEV)\n'
             'print("crust V_CC    = %.4e eV" % gd.VCC_EARTH_CRUST)'),
        md('## One three-flavor probability\n\n'
           'A 1 GeV neutrino travelling 1300 km in vacuum — roughly the DUNE '
           'baseline. The nine probabilities come back with the *initial* '
           'flavor varying slowest.'),
        code('energy = 1.0*GEV\n'
             'baseline = 1300.0*KM\n\n'
             '# The vacuum Hamiltonian at this energy\n'
             'H = np.asarray(H_VAC_3NU)/energy\n\n'
             'prob = oscprob3nu.probabilities_3nu(H, baseline)\n\n'
             'labels = ["ee", "emu", "etau", "mue", "mumu", "mutau",\n'
             '          "taue", "taumu", "tautau"]\n'
             'for name, value in zip(labels, prob):\n'
             '    print("P_%-7s = %.6f" % (name, value))'),
        md('Each initial flavor conserves probability, which is the first '
           'thing to check of any oscillation code:'),
        code('prob = np.array(prob)\n'
             'for start, flavor in zip((0, 3, 6), ("e", "mu", "tau")):\n'
             '    print("sum over final flavors, from nu_%-4s = %.15f"\n'
             '          % (flavor, prob[start:start+3].sum()))'),
        md('## Two flavors\n\n'
           'The same interface, with four probabilities instead of nine.'),
        code('# The builders take sin(theta), not theta itself.\n'
             's12 = np.sqrt(0.310)\n'
             'H2_vac = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_'
             'independent(\n'
             '    s12, gd.D21_NO_BF)\n\n'
             'Pee, Pem, Pme, Pmm = oscprob2nu.probabilities_2nu(\n'
             '    np.asarray(H2_vac)/energy, baseline)\n'
             'print("Pee = %.6f   Pem = %.6f" % (Pee, Pem))'),
        md('## Any Hermitian matrix\n\n'
           'Nothing above is special. The core routines take an arbitrary '
           'Hermitian matrix, which is what makes the library useful for '
           'non-standard scenarios: new interactions, Lorentz-invariance '
           'violation and sterile-like perturbations are all just entries in '
           'that matrix.'),
        code('H_arbitrary = np.array([[0.0, 1.0+2.0j, 0.5],\n'
             '                        [1.0-2.0j, 1.0, 0.0],\n'
             '                        [0.5, 0.0, -1.0]], dtype=complex)\n\n'
             'print(oscprob3nu.probabilities_3nu(H_arbitrary, 1.0))'),
        md('## Pass arrays, do not loop\n\n'
           'The single most useful thing to know about performance here: the '
           'routines broadcast. A stack of Hamiltonians, an array of '
           'baselines, or both, are evaluated in one pass — tens of times '
           'faster than the equivalent Python loop, and the results are '
           'identical.'),
        code('baselines = np.linspace(1.0, 3000.0, 5)*KM\n'
             'prob_scan = oscprob3nu.probabilities_3nu(H, baselines)\n\n'
             'print("shape:", np.shape(prob_scan))\n'
             'print("P_ee along the scan: ", np.round(prob_scan[:, 0], 4))\n'
             'print("all nine at point 0:", np.round(prob_scan[0], 4))'),
    ])

# ---------------------------------------------------------------- 02 vacuum
books['02_vacuum_oscillations.ipynb'] = notebook(
    'Oscillations in vacuum',
    'The probabilities against baseline and against energy, in vacuum — the '
    'figures that `run_testsuite.py` writes to `fig/`, produced here inline '
    'instead.',
    [
        md('## Three flavors, against baseline\n\n'
           'Fixed energy, varying distance. The fast oscillation is the '
           'atmospheric frequency $\\Delta m^2_{31}$; the slow envelope is '
           'the solar one.'),
        code('energy = 1.0*GEV\n'
             'L_km = np.logspace(0.0, 4.5, 3000)\n'
             'H = np.asarray(H_VAC_3NU)/energy\n\n'
             'p = oscprob3nu.probabilities_3nu(H, L_km*KM)\n\n'
             'fig, ax = plt.subplots()\n'
             'ax.semilogx(L_km, p[:, 0], label=r"$P_{ee}$")\n'
             'ax.semilogx(L_km, p[:, 1], label=r"$P_{e\\mu}$")\n'
             'ax.semilogx(L_km, p[:, 2], label=r"$P_{e\\tau}$")\n'
             'ax.set_xlabel("Baseline [km]")\n'
             'ax.set_ylabel("Probability")\n'
             'ax.set_title("Three-flavor vacuum oscillations, "\n'
             '             "E = 1 GeV")\n'
             'ax.set_ylim(0.0, 1.0)\n'
             'ax.legend(loc="center left")\n'
             'plt.show()'),
        md('## Three flavors, against energy\n\n'
           'Fixed baseline, varying energy. The first oscillation maximum is '
           'what a long-baseline experiment is designed to sit on.'),
        code('L = 1300.0*KM\n'
             'E_gev = np.logspace(-1.0, 1.5, 3000)\n'
             'H_stack = np.asarray(H_VAC_3NU)/(E_gev[:, None, None]*GEV)\n\n'
             'p = oscprob3nu.probabilities_3nu(H_stack, L)\n\n'
             'fig, ax = plt.subplots()\n'
             'ax.semilogx(E_gev, p[:, 0], label=r"$P_{ee}$")\n'
             'ax.semilogx(E_gev, p[:, 1], label=r"$P_{e\\mu}$")\n'
             'ax.semilogx(E_gev, p[:, 2], label=r"$P_{e\\tau}$")\n'
             'ax.set_xlabel("Energy [GeV]")\n'
             'ax.set_ylabel("Probability")\n'
             'ax.set_title("Three-flavor vacuum oscillations, L = 1300 km")\n'
             'ax.set_ylim(0.0, 1.0)\n'
             'ax.legend(loc="upper left")\n'
             'plt.show()'),
        md('Note the shape of the call above: the Hamiltonian was built as a '
           'stack of 3000 matrices and passed with a *scalar* baseline. The '
           'library broadcasts the two against each other.'),
        md('## Two flavors\n\n'
           'The two-flavor case has the closed form everyone knows, '
           '$P_{e\\mu} = \\sin^2 2\\theta \\sin^2(\\Delta m^2 L / 4E)$. The '
           'exact SU(2) result reproduces it, which is one of the checks in '
           'the regression suite.'),
        code('# `sth` is sin(theta), not the angle: passing arcsin(s) here\n'
             '# silently changes the mixing angle and the comparison below\n'
             '# then disagrees by about 5%.\n'
             's12 = np.sqrt(0.310)\n'
             'H2 = np.asarray(hamiltonians2nu.'
             'hamiltonian_2nu_vacuum_energy_independent(\n'
             '    s12, gd.D21_NO_BF))/(1.0*GEV)\n\n'
             'L_km = np.logspace(0.0, 5.0, 3000)\n'
             'p2 = oscprob2nu.probabilities_2nu(H2, L_km*KM)\n'
             'Pem = p2[:, 1]\n\n'
             'sin2_2theta = 4.0*s12**2.0*(1.0-s12**2.0)\n'
             'analytic = (sin2_2theta\n'
             '            * np.sin(gd.D21_NO_BF*L_km*KM/(4.0*GEV))**2.0)\n\n'
             'fig, ax = plt.subplots()\n'
             'ax.semilogx(L_km, Pem, label="exact SU(2)")\n'
             'ax.semilogx(L_km, analytic, "--", lw=1.2,\n'
             '            label=r"$\\sin^2 2\\theta \\sin^2(\\Delta m^2 '
             'L/4E)$")\n'
             'ax.set_xlabel("Baseline [km]")\n'
             'ax.set_ylabel(r"$P_{e\\mu}$")\n'
             'ax.set_title("Two-flavor vacuum oscillations, E = 1 GeV")\n'
             'ax.legend(loc="upper left")\n'
             'plt.show()\n\n'
             'print("max |exact - textbook| = %.2e" '
             '% np.max(np.abs(Pem-analytic)))'),
    ])

# ------------------------------------------------------- 03 matter, NSI, LIV
books['03_matter_nsi_liv.ipynb'] = notebook(
    'Matter, NSI, and Lorentz-invariance violation',
    'Constant-density matter, and two kinds of new physics. Each is a '
    'different Hamiltonian passed to the same routine — which is the point '
    'the library is built around.',
    [
        md('## Matter of constant density\n\n'
           'Adding the charged-current potential $V_{CC}$ to the $ee$ entry '
           'produces the MSW resonance: a peak in $P_{e\\mu}$ at the energy '
           'where the matter term cancels the vacuum one.'),
        code('E_gev = np.logspace(-1.0, 2.0, 2000)\n'
             'L = 1300.0*KM\n\n'
             'H_vac = np.asarray(H_VAC_3NU)/(E_gev[:, None, None]*GEV)\n'
             'H_mat = hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '    H_VAC_3NU, E_gev*GEV, gd.VCC_EARTH_CRUST)\n\n'
             'p_vac = oscprob3nu.probabilities_3nu(H_vac, L)\n'
             'p_mat = oscprob3nu.probabilities_3nu(H_mat, L)\n\n'
             'fig, ax = plt.subplots()\n'
             'ax.semilogx(E_gev, p_vac[:, 1], label="vacuum")\n'
             'ax.semilogx(E_gev, p_mat[:, 1], label="matter, 3 g/cm$^3$")\n'
             'ax.set_xlabel("Energy [GeV]")\n'
             'ax.set_ylabel(r"$P_{e\\mu}$")\n'
             'ax.set_title("Matter effects at L = 1300 km")\n'
             'ax.legend()\n'
             'plt.show()'),
        md('The builders take an **array** of energies and return one '
           'Hamiltonian per energy, stacked. That is why the matter curve '
           'above needed no loop.'),
        md('## Non-standard interactions\n\n'
           'NSI add an effective potential with its own flavor structure, '
           'parametrised by the $\\epsilon$ matrix.'),
        code('H_nsi = hamiltonians3nu.hamiltonian_3nu_nsi(\n'
             '    H_VAC_3NU, E_gev*GEV, gd.VCC_EARTH_CRUST, gd.EPS_3)\n'
             'p_nsi = oscprob3nu.probabilities_3nu(H_nsi, L)\n\n'
             'fig, ax = plt.subplots()\n'
             'ax.semilogx(E_gev, p_mat[:, 1], label="standard matter")\n'
             'ax.semilogx(E_gev, p_nsi[:, 1], label="with NSI")\n'
             'ax.set_xlabel("Energy [GeV]")\n'
             'ax.set_ylabel(r"$P_{e\\mu}$")\n'
             'ax.set_title("Non-standard interactions, L = 1300 km")\n'
             'ax.legend()\n'
             'plt.show()'),
        md('## Lorentz-invariance violation\n\n'
           'A CPT-odd LIV background contributes a term that grows with '
           'energy rather than falling, so it shows up at the high end.'),
        code('H_liv = hamiltonians3nu.hamiltonian_3nu_liv(\n'
             '    H_VAC_3NU, E_gev*GEV, gd.SXI12, gd.SXI23, gd.SXI13,\n'
             '    gd.DXICP, gd.B1, gd.B2, gd.B3, gd.LAMBDA)\n'
             'p_liv = oscprob3nu.probabilities_3nu(H_liv, L)\n\n'
             'fig, ax = plt.subplots()\n'
             'ax.semilogx(E_gev, p_vac[:, 1], label="vacuum")\n'
             'ax.semilogx(E_gev, p_liv[:, 1], label="LIV background")\n'
             'ax.set_xlabel("Energy [GeV]")\n'
             'ax.set_ylabel(r"$P_{e\\mu}$")\n'
             'ax.set_title("Lorentz-invariance violation, L = 1300 km")\n'
             'ax.legend()\n'
             'plt.show()'),
    ])

# ----------------------------------------------------------- 04 oscillogram
books['04_oscillogram.ipynb'] = notebook(
    'Oscillograms',
    'A probability over a two-dimensional grid of energy and baseline. This '
    'is the case broadcasting was built for: the whole map is one call, with '
    'no Python loop anywhere.',
    [
        md('## Building the grid\n\n'
           'The trick is the shape of the two arguments. Hamiltonians get a '
           'trailing pair of matrix axes and a leading energy axis; baselines '
           'get a leading axis of their own. Indexing with `None` lines them '
           'up so they broadcast into a grid.'),
        code('n_e, n_l = 240, 240\n'
             'E_gev = np.logspace(-1.0, 1.5, n_e)\n'
             'L_km = np.linspace(50.0, 12000.0, n_l)\n\n'
             'H = hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '    H_VAC_3NU, E_gev*GEV, gd.VCC_EARTH_CRUST)\n\n'
             '# (n_e, 1, 3, 3) against (1, n_l) -> an (n_e, n_l) grid\n'
             'p = oscprob3nu.probabilities_3nu(H[:, None, :, :],\n'
             '                                 (L_km*KM)[None, :])\n'
             'print("returned shape:", np.shape(p))'),
        md('## $P_{\\mu e}$ over energy and baseline'),
        code('fig, ax = plt.subplots(figsize=(7.6, 4.8))\n'
             'mesh = ax.pcolormesh(L_km, E_gev, p[:, :, 3], shading="auto",\n'
             '                     cmap="viridis", vmin=0.0)\n'
             'ax.set_yscale("log")\n'
             'ax.set_xlabel("Baseline [km]")\n'
             'ax.set_ylabel("Energy [GeV]")\n'
             'ax.set_title(r"$P_{\\mu e}$ in matter of constant density")\n'
             'ax.grid(False)\n'
             'fig.colorbar(mesh, ax=ax, label=r"$P_{\\mu e}$")\n'
             'plt.show()'),
        md('## Survival, on the same grid\n\n'
           'The muon-neutrino survival probability shows the oscillation '
           'valleys that atmospheric experiments measure.'),
        code('fig, ax = plt.subplots(figsize=(7.6, 4.8))\n'
             'mesh = ax.pcolormesh(L_km, E_gev, p[:, :, 4], shading="auto",\n'
             '                     cmap="magma", vmin=0.0, vmax=1.0)\n'
             'ax.set_yscale("log")\n'
             'ax.set_xlabel("Baseline [km]")\n'
             'ax.set_ylabel("Energy [GeV]")\n'
             'ax.set_title(r"$P_{\\mu\\mu}$ in matter of constant density")\n'
             'ax.grid(False)\n'
             'fig.colorbar(mesh, ax=ax, label=r"$P_{\\mu\\mu}$")\n'
             'plt.show()'),
        md('Both maps above are a single call each on a '
           '240 x 240 grid — 57 600 probabilities. Written as a Python loop '
           'this is tens of times slower; see the Performance section of the '
           'README.'),
    ])

# --------------------------------------------------------- 05 biprobability
books['05_biprobability.ipynb'] = notebook(
    'Bi-probability plots',
    'The classic way to see CP violation: plot $P(\\nu_\\mu \\to \\nu_e)$ '
    'against $P(\\bar\\nu_\\mu \\to \\bar\\nu_e)$ as the CP phase runs '
    'through $2\\pi$. The locus is an ellipse, and it collapses to a point '
    'on the diagonal when there is nothing to see.',
    [
        md('## Neutrinos and antineutrinos\n\n'
           'The antineutrino Hamiltonian is the complex conjugate of the '
           'neutrino one with the matter potential reversed. Building both '
           'and scanning $\\delta_{CP}$ gives the ellipse.'),
        code('def biprobability(energy_gev, baseline_km, vcc, n_phase=400):\n'
             '    """Returns P(numu->nue) and its antineutrino partner over '
             'delta_CP."""\n'
             '    dcp = np.linspace(0.0, 2.0*np.pi, n_phase)\n'
             '    p_nu, p_nubar = [], []\n'
             '    for d in dcp:\n'
             '        h_vac = hamiltonians3nu.'
             'hamiltonian_3nu_vacuum_energy_independent(\n'
             '            gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, d,\n'
             '            gd.D21_NO_BF, gd.D31_NO_BF)\n'
             '        h_nu = hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '            h_vac, energy_gev*GEV, vcc)\n'
             '        # Antineutrinos: conjugate the Hamiltonian, flip V_CC\n'
             '        h_bar = hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '            np.conj(h_vac), energy_gev*GEV, -vcc)\n'
             '        p_nu.append(oscprob3nu.probabilities_3nu(\n'
             '            h_nu, baseline_km*KM)[3])\n'
             '        p_nubar.append(oscprob3nu.probabilities_3nu(\n'
             '            h_bar, baseline_km*KM)[3])\n'
             '    return np.array(p_nu), np.array(p_nubar), dcp'),
        md('## In vacuum\n\n'
           'With no matter the ellipse is centred on the diagonal: CP '
           'violation moves you along it, but neutrinos and antineutrinos '
           'stay symmetric.'),
        code('fig, ax = plt.subplots(figsize=(5.4, 5.2))\n'
             'for e_gev in (0.6, 1.0, 2.0):\n'
             '    p, pbar, dcp = biprobability(e_gev, 1300.0, 0.0)\n'
             '    ax.plot(p, pbar, label="E = %.1f GeV" % e_gev)\n'
             'lims = ax.get_xlim()\n'
             'ax.plot(lims, lims, "k:", lw=1, label="CP symmetric")\n'
             'ax.set_xlabel(r"$P(\\nu_\\mu \\to \\nu_e)$")\n'
             'ax.set_ylabel(r"$P(\\bar\\nu_\\mu \\to \\bar\\nu_e)$")\n'
             'ax.set_title("Bi-probability in vacuum, L = 1300 km")\n'
             'ax.legend()\n'
             'ax.set_aspect("equal", adjustable="datalim")\n'
             'plt.show()'),
        md('## In matter\n\n'
           'Matter breaks the symmetry between neutrinos and antineutrinos on '
           'its own, and pushes the ellipses off the diagonal. Disentangling '
           'that from genuine CP violation is the central difficulty of a '
           'long-baseline measurement.'),
        code('fig, ax = plt.subplots(figsize=(5.4, 5.2))\n'
             'for e_gev in (0.6, 1.0, 2.0):\n'
             '    p, pbar, dcp = biprobability(e_gev, 1300.0,\n'
             '                                 gd.VCC_EARTH_CRUST)\n'
             '    ax.plot(p, pbar, label="E = %.1f GeV" % e_gev)\n'
             'lims = ax.get_xlim()\n'
             'ax.plot(lims, lims, "k:", lw=1, label="CP symmetric")\n'
             'ax.set_xlabel(r"$P(\\nu_\\mu \\to \\nu_e)$")\n'
             'ax.set_ylabel(r"$P(\\bar\\nu_\\mu \\to \\bar\\nu_e)$")\n'
             'ax.set_title("Bi-probability in matter, L = 1300 km")\n'
             'ax.legend()\n'
             'ax.set_aspect("equal", adjustable="datalim")\n'
             'plt.show()'),
        md('## Where on the ellipse\n\n'
           'Marking a few values of $\\delta_{CP}$ shows which part of the '
           'curve corresponds to which phase.'),
        code('p, pbar, dcp = biprobability(1.0, 1300.0, gd.VCC_EARTH_CRUST)\n\n'
             'fig, ax = plt.subplots(figsize=(5.4, 5.2))\n'
             'ax.plot(p, pbar, color="0.4")\n'
             'for target, label in ((0.0, r"$\\delta_{CP} = 0$"),\n'
             '                      (np.pi/2, r"$\\pi/2$"),\n'
             '                      (np.pi, r"$\\pi$"),\n'
             '                      (3*np.pi/2, r"$3\\pi/2$")):\n'
             '    k = int(np.argmin(np.abs(dcp-target)))\n'
             '    ax.plot(p[k], pbar[k], "o")\n'
             '    ax.annotate(label, (p[k], pbar[k]),\n'
             '                textcoords="offset points", xytext=(8, 4))\n'
             'ax.set_xlabel(r"$P(\\nu_\\mu \\to \\nu_e)$")\n'
             'ax.set_ylabel(r"$P(\\bar\\nu_\\mu \\to \\bar\\nu_e)$")\n'
             'ax.set_title("E = 1 GeV, L = 1300 km, matter")\n'
             'plt.show()'),
    ])

# ---------------------------------------------------------- 06 Earth / PREM
books['06_earth_and_prem.ipynb'] = notebook(
    'The Earth: PREM, chords, and slabs',
    'A neutrino crossing the Earth does not see a constant Hamiltonian, so '
    'the exact expansions do not apply to its whole trajectory. They apply to '
    'any piece over which the density is taken constant — which is what the '
    '`earth` and `slabs` modules build.',
    [
        code('import earth\nimport slabs'),
        md('## The density profile\n\n'
           'PREM is a piecewise-polynomial fit to seismological data. Note '
           'the jumps between shells, and that the density is *not* constant '
           'within a shell.'),
        code('r = np.linspace(0.0, gd.EARTH_RADIUS, 4000)\n'
             'rho = earth.density_prem(r)\n\n'
             'fig, ax = plt.subplots()\n'
             'ax.plot(r, rho)\n'
             'for b in earth.PREM_BOUNDARIES:\n'
             '    ax.axvline(b, color="0.8", lw=0.8, zorder=0)\n'
             'ax.set_xlabel("Radius [km]")\n'
             'ax.set_ylabel(r"Density [g cm$^{-3}$]")\n'
             'ax.set_title("Preliminary Reference Earth Model")\n'
             'plt.show()'),
        md('A strong check on the profile: integrating it over the volume of '
           'the Earth must give the Earth\'s mass.'),
        code('integrand = 4.0*np.pi*(r*1.e5)**2.0*rho\n'
             'x = r*1.e5\n'
             'mass_kg = np.sum(0.5*(integrand[1:]+integrand[:-1])'
             '*np.diff(x))*1.e-3\n'
             'print("PREM integrates to  %.4e kg" % mass_kg)\n'
             'print("the Earth\'s mass is 5.972e+24 kg")'),
        md('## The chord, and where it crosses shells\n\n'
           'A neutrino arriving from below travels a chord set by '
           '$\\cos\\theta_z$. Along it the radius falls to a minimum and '
           'rises again, so the chord crosses most shells **twice**.'),
        code('fig, ax = plt.subplots()\n'
             'for costhz in (-0.2, -0.5, -0.8, -1.0):\n'
             '    d = earth.distance_traveled_inside_earth(costhz)\n'
             '    l = np.linspace(0.0, d, 800)\n'
             '    ax.plot(l, earth.density_prem(\n'
             '        earth.earth_radial_distance_from_depth(costhz, l)),\n'
             '        label=r"$\\cos\\theta_z = %.1f$" % costhz)\n'
             'ax.set_xlabel("Distance along the chord [km]")\n'
             'ax.set_ylabel(r"Density [g cm$^{-3}$]")\n'
             'ax.set_title("Density seen along a trajectory")\n'
             'ax.legend()\n'
             'plt.show()'),
        code('costhz = -1.0\n'
             'edges = earth.prem_layer_edges_along_chord(costhz)\n'
             'print("PREM shells:            10")\n'
             'print("boundary crossings:     %d" % len(edges))\n'
             'print("chord segments:         %d" % (len(edges)+1))\n'
             'print("-> a chord enters and leaves each shell it reaches")'),
        md('## From chord to slabs\n\n'
           'Two different things decide where the slabs are cut. **Between** '
           'shells the density jumps, so the chord is split at every boundary '
           'crossing. **Within** a shell it varies smoothly, so each segment '
           'is divided further into `n_slabs_per_segment` equal sub-slabs '
           'with the density taken at the midpoint.'),
        code('widths, densities = earth.earth_slabs(-0.8, '
             'n_slabs_per_segment=6)\n\n'
             'fig, ax = plt.subplots()\n'
             'ax.step(np.concatenate(([0.0], np.cumsum(widths))),\n'
             '        np.concatenate((densities[:1], densities)),\n'
             '        where="pre", label="slabs")\n'
             'l = np.linspace(0.0, widths.sum(), 1500)\n'
             'ax.plot(l, earth.density_prem(\n'
             '    earth.earth_radial_distance_from_depth(-0.8, l)),\n'
             '    "k--", lw=1, label="PREM")\n'
             'ax.set_xlabel("Distance along the chord [km]")\n'
             'ax.set_ylabel(r"Density [g cm$^{-3}$]")\n'
             'ax.set_title(r"Slab approximation, $\\cos\\theta_z = -0.8$, "\n'
             '             "6 sub-slabs per segment")\n'
             'ax.legend()\n'
             'plt.show()\n\n'
             'print("%d slabs, total %.1f km" % (len(widths), widths.sum()))'),
        md('## Convergence\n\n'
           'The slicing is the only approximation in the whole calculation, '
           'so it is worth watching it converge rather than trusting it. '
           'Midpoint sampling is second order: past about 32 sub-slabs per '
           'segment, each doubling cuts the error by roughly four.'),
        code('ref = np.array(earth.probabilities_3nu_earth(\n'
             '    H_VAC_3NU, 10.0*GEV, -0.8, n_slabs_per_segment=512))\n\n'
             'ns = [1, 2, 4, 8, 16, 32, 64, 128]\n'
             'err = [np.max(np.abs(np.array(earth.probabilities_3nu_earth(\n'
             '           H_VAC_3NU, 10.0*GEV, -0.8, n_slabs_per_segment=n))\n'
             '       - ref)) for n in ns]\n\n'
             'fig, ax = plt.subplots()\n'
             'ax.loglog(ns, err, "o-", label="measured")\n'
             'ax.loglog(ns, err[0]*(np.array(ns)/ns[0])**-2.0, "k--", lw=1,\n'
             '          label=r"second order, $n^{-2}$")\n'
             'ax.set_xlabel("Sub-slabs per segment")\n'
             'ax.set_ylabel("Max deviation from the converged result")\n'
             'ax.set_title("Convergence of the slab approximation")\n'
             'ax.legend()\n'
             'plt.show()'),
    ])

# -------------------------------------------------- 07 Earth probabilities
books['07_earth_probabilities.ipynb'] = notebook(
    'Probabilities through the Earth',
    'Putting PREM and the slab machinery together: oscillations for '
    'neutrinos that have crossed the Earth, and between two points on its '
    'surface.',
    [
        code('import earth'),
        md('## Against zenith angle\n\n'
           'At fixed energy, sweeping $\\cos\\theta_z$ from $0$ (horizontal) '
           'to $-1$ (straight through the centre) traces out the matter '
           'effect as the neutrino passes through progressively denser '
           'material.'),
        code('costhz = np.linspace(-0.999, -0.02, 220)\n\n'
             'for e_gev, style in ((5.0, "-"), (10.0, "--")):\n'
             '    p = np.array([earth.probabilities_3nu_earth(\n'
             '        H_VAC_3NU, e_gev*GEV, c, n_slabs_per_segment=6)\n'
             '        for c in costhz])\n'
             '    plt.plot(costhz, p[:, 4], style,\n'
             '             label="E = %.0f GeV" % e_gev)\n\n'
             'plt.xlabel(r"$\\cos\\theta_z$")\n'
             'plt.ylabel(r"$P_{\\mu\\mu}$")\n'
             'plt.title("Muon-neutrino survival through the Earth")\n'
             'plt.legend()\n'
             'plt.show()'),
        md('The core-mantle boundary leaves a visible feature: near '
           '$\\cos\\theta_z \\approx -0.83$ the trajectory starts to clip the '
           'outer core, and the density it crosses jumps.'),
        md('## An Earth oscillogram\n\n'
           'Energy against zenith angle — the map an atmospheric-neutrino '
           'experiment measures. Each point is a full slab calculation, so '
           'this one does loop; the grid is kept modest for that reason.'),
        code('n_e, n_c = 60, 60\n'
             'E_gev = np.logspace(0.0, 2.0, n_e)\n'
             'cz = np.linspace(-0.999, -0.05, n_c)\n\n'
             'grid = np.empty((n_e, n_c))\n'
             'for i, e in enumerate(E_gev):\n'
             '    for j, c in enumerate(cz):\n'
             '        grid[i, j] = earth.probabilities_3nu_earth(\n'
             '            H_VAC_3NU, e*GEV, c, n_slabs_per_segment=4)[4]\n\n'
             'fig, ax = plt.subplots(figsize=(7.6, 4.8))\n'
             'mesh = ax.pcolormesh(cz, E_gev, grid, shading="auto",\n'
             '                     cmap="magma", vmin=0.0, vmax=1.0)\n'
             'ax.set_yscale("log")\n'
             'ax.set_xlabel(r"$\\cos\\theta_z$")\n'
             'ax.set_ylabel("Energy [GeV]")\n'
             'ax.set_title(r"$P_{\\mu\\mu}$ through the Earth (PREM)")\n'
             'ax.grid(False)\n'
             'fig.colorbar(mesh, ax=ax, label=r"$P_{\\mu\\mu}$")\n'
             'plt.show()'),
        md('## Between two places\n\n'
           'The library knows a set of named locations — the same set as the '
           'sibling Magnus package — and will work out the chord between any '
           'two of them.'),
        code('print("known locations:")\n'
             'print("  " + ", ".join(sorted(earth.LOC_COORDS_DMS)))'),
        code('pairs = [("cern", "gran_sasso"),\n'
             '         ("fermilab", "homestake"),\n'
             '         ("tokai", "kamioka"),\n'
             '         ("cern", "kamioka"),\n'
             '         ("south_pole", "north_pole")]\n\n'
             'print("%-26s %10s %12s" % ("baseline", "chord [km]", '
             '"cos(theta_z)"))\n'
             'for a, b in pairs:\n'
             '    la, lo = earth.coordinates_of_named_location(a)\n'
             '    lb, ob = earth.coordinates_of_named_location(b)\n'
             '    d = earth.chord_length_inside_earth(la, lo, lb, ob)\n'
             '    c = earth.costhz_between_points_on_surface(la, lo, lb, ob)\n'
             '    print("%-26s %10.1f %12.4f" % (a+" -> "+b, d, c))'),
        md('These are a check on the geometry, not just a convenience. The '
           'chords come out at the baselines the experiments actually quote: '
           'CERN to Gran Sasso is the 730 km of CNGS, Tokai to Kamioka the '
           '295 km of T2K, and pole to pole is exactly one Earth diameter. '
           'Fermilab to Homestake gives 1285 km against DUNE\'s quoted '
           '1300 km, the difference being that DUNE quotes the distance to '
           'the detector hall rather than the surface chord.'),
        md('And the probability along one of those trajectories, against '
           'energy:'),
        code('E_gev = np.logspace(-0.5, 1.5, 150)\n\n'
             'fig, ax = plt.subplots()\n'
             'for a, b in (("cern", "gran_sasso"), ("fermilab", '
             '"homestake")):\n'
             '    p = np.array([earth.probabilities_3nu_between_locations(\n'
             '        H_VAC_3NU, e*GEV, a, b, n_slabs_per_segment=6)\n'
             '        for e in E_gev])\n'
             '    ax.semilogx(E_gev, p[:, 3], label="%s to %s" % (a, b))\n\n'
             'ax.set_xlabel("Energy [GeV]")\n'
             'ax.set_ylabel(r"$P_{\\mu e}$")\n'
             'ax.set_title("Appearance along real baselines, through PREM")\n'
             'ax.legend()\n'
             'plt.show()'),
    ])

# -------------------------------------------------- 08 unusual profiles
books['08_unusual_density_profiles.ipynb'] = notebook(
    'Unusual density profiles',
    'The Earth is one profile among many. `slabs` takes any sequence of '
    'widths and Hamiltonians, so a profile can be built by hand — and some '
    'of the interesting ones are not remotely Earth-like.\n\n'
    'The point of this notebook is that the *shape* of the profile matters, '
    'not only its average. A castle wall and a uniform slab of the same mean '
    'density give different probabilities, sometimes very different.',
    [
        code('import earth\nimport slabs'),
        code('def probabilities_profile(widths_km, densities, energy_ev,\n'
             '                          h_vac=None):\n'
             '    """Nine probabilities through an arbitrary matter profile."""\n'
             '    h_vac = H_VAC_3NU if h_vac is None else h_vac\n'
             '    vcc = earth.matter_potential(np.asarray(densities, '
             'dtype=float))\n'
             '    H = hamiltonians3nu.hamiltonian_3nu_matter(h_vac, energy_ev,'
             ' vcc)\n'
             '    return np.array(slabs.probabilities_3nu_slabs(\n'
             '        H, np.asarray(widths_km, dtype=float)*KM))\n\n\n'
             'def show_profile(ax, widths_km, densities, **kw):\n'
             '    """Draws a profile as the step function it is."""\n'
             '    edges = np.concatenate(([0.0], np.cumsum(widths_km)))\n'
             '    ax.step(edges, np.concatenate((densities[:1], densities)),\n'
             '            where="pre", **kw)\n'
             '    ax.set_xlim(0.0, edges[-1])'),
        md('## Three profiles, one mean density\n\n'
           'A **castle wall** alternates between two densities. A **serrated** '
           'profile ramps up in steps and drops back. A **random** wall keeps '
           'the same two levels but shuffles them. All three below have the '
           'same total length and the same mean density as the uniform case, '
           'so any difference in the probabilities is due to the arrangement '
           'alone.'),
        code('total_km = 6000.0\n'
             'n_slab = 24\n'
             'rho_lo, rho_hi = 2.0, 8.0\n'
             'mean_rho = 0.5*(rho_lo+rho_hi)\n\n'
             'widths = np.full(n_slab, total_km/n_slab)\n\n'
             'castle = np.where(np.arange(n_slab) % 2 == 0, rho_lo, rho_hi)\n'
             'serrated = np.tile(np.linspace(rho_lo, rho_hi, 6), 4)\n'
             'rng = np.random.default_rng(20260801)\n'
             '# A *permutation* of the castle wall, not a fresh random draw: that\n'
             '# guarantees the identical multiset of slabs, and so exactly the\n'
             '# same mean density, which is the whole premise of the comparison.\n'
             'random_wall = rng.permutation(castle)\n'
             'uniform = np.full(n_slab, mean_rho)\n\n'
             'profiles = [("castle wall", castle),\n'
             '            ("serrated", serrated),\n'
             '            ("random wall", random_wall),\n'
             '            ("uniform", uniform)]\n\n'
             'for name, rho in profiles:\n'
             '    print("%-12s mean = %.3f g/cm^3" % (name, rho.mean()))'),
        code('fig, axes = plt.subplots(4, 1, figsize=(7.2, 7.6), '
             'sharex=True)\n'
             'for ax, (name, rho) in zip(axes, profiles):\n'
             '    show_profile(ax, widths, rho, color="C0")\n'
             '    ax.axhline(mean_rho, color="0.6", ls=":", lw=1)\n'
             '    ax.set_ylabel(name, fontsize=9)\n'
             '    ax.set_ylim(0.0, rho_hi*1.0)\n'
             'axes[-1].set_xlabel("Distance [km]")\n'
             'fig.suptitle("Four profiles with the same mean density")\n'
             'fig.tight_layout()\n'
             'plt.show()'),
        md('## The probabilities they produce\n\n'
           'Scanning energy through each profile. Where the curves separate, '
           'the arrangement of the matter is doing the work — the mean '
           'density is identical throughout.'),
        code('E_gev = np.logspace(-0.7, 1.7, 400)\n\n'
             'fig, ax = plt.subplots()\n'
             'for name, rho in profiles:\n'
             '    p = np.array([probabilities_profile(widths, rho, e*GEV)[3]\n'
             '                  for e in E_gev])\n'
             '    ax.semilogx(E_gev, p, label=name,\n'
             '                lw=1.6 if name == "uniform" else 1.2,\n'
             '                ls="--" if name == "uniform" else "-")\n'
             'ax.set_xlim(E_gev[0], E_gev[-1])\n'
             'ax.set_xlabel("Energy [GeV]")\n'
             'ax.set_ylabel(r"$P_{\\mu e}$")\n'
             'ax.set_title("Same mean density, different arrangement")\n'
             'ax.legend()\n'
             'plt.show()'),
        md('## Parametric enhancement\n\n'
           'A castle wall is not an arbitrary choice. When the width of one '
           'period is tuned against the oscillation length, successive layers '
           'add coherently and the transition probability can be driven far '
           'above what either density gives on its own — a parametric '
           'resonance. Sweeping the number of periods at fixed total length '
           'shows it building up.'),
        code('E_fixed = 3.0*GEV\n'
             'periods = np.arange(1, 41)\n\n'
             'p_castle, p_uniform = [], []\n'
             'for n_period in periods:\n'
             '    n = 2*n_period\n'
             '    w = np.full(n, total_km/n)\n'
             '    rho = np.where(np.arange(n) % 2 == 0, rho_lo, rho_hi)\n'
             '    p_castle.append(probabilities_profile(w, rho, E_fixed)[3])\n'
             '    p_uniform.append(probabilities_profile(\n'
             '        w, np.full(n, mean_rho), E_fixed)[3])\n\n'
             'fig, ax = plt.subplots()\n'
             'ax.plot(periods, p_castle, "o-", ms=3, label="castle wall")\n'
             'ax.plot(periods, p_uniform, "--", label="uniform, same mean")\n'
             'ax.set_xlim(periods[0], periods[-1])\n'
             'ax.set_xlabel("Number of periods over the same 6000 km")\n'
             'ax.set_ylabel(r"$P_{\\mu e}$")\n'
             'ax.set_title("Castle-wall profile at E = 3 GeV")\n'
             'ax.legend()\n'
             'plt.show()'),
        md('## Sensitivity to the arrangement\n\n'
           'Twenty random walls, all with the same two densities and the same '
           'mean. The spread between them is a direct measure of how much the '
           'ordering matters at this energy.'),
        code('E_fixed = 3.0*GEV\n'
             'rng = np.random.default_rng(7)\n\n'
             'curves = []\n'
             'E_gev = np.logspace(-0.3, 1.3, 200)\n'
             'for _ in range(20):\n'
             '    rho = rng.permutation(castle)\n'
             '    curves.append([probabilities_profile(widths, rho, e*GEV)[3]\n'
             '                   for e in E_gev])\n'
             'curves = np.array(curves)\n\n'
             'fig, ax = plt.subplots()\n'
             'ax.semilogx(E_gev, curves.T, color="C0", alpha=0.25, lw=0.9)\n'
             'ax.semilogx(E_gev, [probabilities_profile(\n'
             '    widths, uniform, e*GEV)[3] for e in E_gev],\n'
             '    "k--", lw=1.4, label="uniform, same mean")\n'
             'ax.set_xlim(E_gev[0], E_gev[-1])\n'
             'ax.set_xlabel("Energy [GeV]")\n'
             'ax.set_ylabel(r"$P_{\\mu e}$")\n'
             'ax.set_title("Twenty random castle walls, identical mean '
             'density")\n'
             'ax.legend()\n'
             'plt.show()'),
        md('Every curve above is exact. Nothing here is a perturbative or '
           'adiabatic approximation: each slab is solved in closed form and '
           'the operators are multiplied, so the only thing separating these '
           'results from the continuum is how the profile itself was '
           'defined.'),
    ])

# ------------------------------------------------------------ 09 performance
books['09_performance.ipynb'] = notebook(
    'Performance',
    'Two things make this library fast, and they are worth knowing in order '
    'of importance: **pass arrays instead of looping**, which costs nothing '
    'and wins the most, and **install the optional compiled backend** if the '
    'scans are large.\n\n'
    'Every number below is measured when the notebook runs, so it reflects '
    'the machine that produced this copy rather than a figure carried over '
    'from a README.',
    [
        code('import time\n\n'
             'import fastkernels\n\n\n'
             'def best_of(func, repeat=5):\n'
             '    """Returns the fastest of `repeat` runs, in seconds.\n\n'
             '    The minimum, not the mean: timing noise is one-sided, so\n'
             '    the fastest run is the one least polluted by whatever else\n'
             '    the machine was doing.\n'
             '    """\n'
             '    best = float("inf")\n'
             '    for _ in range(repeat):\n'
             '        t0 = time.perf_counter()\n'
             '        func()\n'
             '        best = min(best, time.perf_counter()-t0)\n'
             '    return best\n\n\n'
             'print("Numba available:", fastkernels.HAVE_NUMBA)'),
        md('## Looping versus broadcasting\n\n'
           'The same energy scan, done both ways. The batched call is one '
           'invocation; the loop is one per point.'),
        code('n_points = 2000\n'
             'E_gev = np.logspace(-1.0, 1.5, n_points)\n'
             'L = 1300.0*KM\n\n'
             'H_stack = hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '    H_VAC_3NU, E_gev*GEV, gd.VCC_EARTH_CRUST)\n\n\n'
             'def batched():\n'
             '    return oscprob3nu.probabilities_3nu(H_stack, L)\n\n\n'
             'def looped():\n'
             '    return [oscprob3nu.probabilities_3nu(h, L) for h in '
             'H_stack]\n\n\n'
             't_batch = best_of(batched)\n'
             't_loop = best_of(looped, repeat=3)\n\n'
             'print("loop over %d points : %8.2f ms" % (n_points, '
             't_loop*1e3))\n'
             'print("one batched call   : %8.2f ms" % (t_batch*1e3))\n'
             'print("speedup            : %8.1fx" % (t_loop/t_batch))'),
        md('The two agree exactly — this is the same arithmetic, organised '
           'differently:'),
        code('print("max |batched - looped| = %.2e"\n'
             '      % np.max(np.abs(np.array(batched()) - '
             'np.array(looped()))))'),
        md('## How the gain scales\n\n'
           'The advantage grows with the size of the scan, because the fixed '
           'cost of a batched call is amortised over more points. Below a '
           'handful of elements the library detects this and takes the scalar '
           'path instead — that is what `oscprob3nu.SMALL_BATCH` is for.'),
        code('sizes = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000]\n'
             'rows = []\n'
             'for n in sizes:\n'
             '    E = np.logspace(-1.0, 1.5, n)\n'
             '    Hs = hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '        H_VAC_3NU, E*GEV, gd.VCC_EARTH_CRUST)\n'
             '    tb = best_of(lambda: oscprob3nu.probabilities_3nu(Hs, L))\n'
             '    tl = best_of(lambda: [oscprob3nu.probabilities_3nu(h, L)\n'
             '                          for h in Hs], repeat=3)\n'
             '    rows.append((n, tl, tb))\n\n'
             'rows = np.array(rows)\n'
             'print("%8s %12s %12s %9s" % ("N", "loop [ms]", "batched [ms]",\n'
             '                            "speedup"))\n'
             'for n, tl, tb in rows:\n'
             '    print("%8d %12.3f %12.3f %8.1fx" % (n, tl*1e3, tb*1e3, '
             'tl/tb))'),
        code('fig, ax = plt.subplots()\n'
             'ax.loglog(rows[:, 0], rows[:, 1]*1e3, "o-", label="Python '
             'loop")\n'
             'ax.loglog(rows[:, 0], rows[:, 2]*1e3, "s-", label="batched")\n'
             'ax.set_xlim(rows[0, 0], rows[-1, 0])\n'
             'ax.set_xlabel("Points in the scan")\n'
             'ax.set_ylabel("Time [ms]")\n'
             'ax.set_title("Cost of an energy scan, three flavors")\n'
             'ax.legend()\n'
             'plt.show()'),
        md('## The compiled backend\n\n'
           'With `numba` installed, the batched paths hand large stacks to a '
           'compiled kernel. It is switchable at runtime, which is how the '
           'regression suite checks the two paths against each other — and is '
           'what makes the comparison below possible in one process.'),
        code('if not fastkernels.HAVE_NUMBA:\n'
             '    print("numba is not installed; skipping this comparison.")\n'
             'else:\n'
             '    n = 200000\n'
             '    E = np.logspace(-1.0, 1.5, n)\n'
             '    Hs = hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '        H_VAC_3NU, E*GEV, gd.VCC_EARTH_CRUST)\n\n'
             '    # Warm the cache: the first call compiles.\n'
             '    oscprob3nu.probabilities_3nu(Hs[:100], L)\n\n'
             '    fastkernels.USE_NUMBA = True\n'
             '    t_kernel = best_of(lambda: '
             'oscprob3nu.probabilities_3nu(Hs, L))\n'
             '    fastkernels.USE_NUMBA = False\n'
             '    t_numpy = best_of(lambda: '
             'oscprob3nu.probabilities_3nu(Hs, L))\n'
             '    fastkernels.USE_NUMBA = True\n\n'
             '    print("%d energies, three flavors" % n)\n'
             '    print("  NumPy path     : %8.1f ms" % (t_numpy*1e3))\n'
             '    print("  compiled kernel: %8.1f ms" % (t_kernel*1e3))\n'
             '    print("  speedup        : %8.1fx" % (t_numpy/t_kernel))'),
        md('The two paths agree to round-off, which is the property that '
           'makes the backend safe to install and forget about:'),
        code('if fastkernels.HAVE_NUMBA:\n'
             '    Hs = hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '        H_VAC_3NU, np.logspace(-1, 1.5, 5000)*GEV,\n'
             '        gd.VCC_EARTH_CRUST)\n'
             '    fastkernels.USE_NUMBA = True\n'
             '    a = np.array(oscprob3nu.probabilities_3nu(Hs, L))\n'
             '    fastkernels.USE_NUMBA = False\n'
             '    b = np.array(oscprob3nu.probabilities_3nu(Hs, L))\n'
             '    fastkernels.USE_NUMBA = True\n'
             '    print("max |kernel - numpy| = %.2e" % np.max(np.abs(a-b)))'),
        md('## Where the backend does *not* help\n\n'
           'The two-flavor expansion reduces to a square root and a sine per '
           'element, which NumPy already does about as well as compiled code '
           'can. The library therefore declines the kernel below a measured '
           'threshold rather than using it everywhere — a backend that is '
           'sometimes slower than the path it replaces is worse than no '
           'backend.'),
        code('print("thresholds at which the compiled kernel is used:")\n'
             'for n_flavors, threshold in sorted('
             'fastkernels.MIN_BATCH.items()):\n'
             '    print("  %d flavors: stacks of %d or more"\n'
             '          % (n_flavors, threshold))\n'
             'print()\n'
             'print("below these sizes the NumPy path is used, whether or "\n'
             '      "not numba is installed")'),
        md('## What to take away\n\n'
           '1. Replace loops with array arguments. This is the large win, it '
           'needs no extra dependency, and the results are identical.\n'
           '2. If the scans are large and three-flavor, `pip install '
           '"nuoscprobexact[fast]"` on top of that.\n'
           '3. Do not bother for a handful of points: the library already '
           'takes the quicker path there on its own.'),
    ])

# ------------------------------------------------------- 10 the paper figures
books['10_paper_figures.ipynb'] = notebook(
    "The paper's figures",
    'The two figures in [arXiv:1904.12391](https://arxiv.org/abs/1904.12391), '
    'reproduced here. They were previously drawn by two plotting modules '
    'that this notebook replaces: three-flavor and two-flavor '
    'probabilities against energy at the DUNE baseline, for four scenarios at '
    'once.',
    [
        md('## Three flavors\n\n'
           'Vacuum, constant-density matter, matter with non-standard '
           'interactions, and a CPT-odd Lorentz invariance-violating '
           'background, at $L = 1300$ km. All four use the default parameters '
           'in `globaldefs`, so the figure is reproducible from the '
           'repository alone.'),
        code('L = 1300.0*KM\n'
             'E_gev = np.logspace(np.log10(0.5), np.log10(30.0), 400)\n'
             'E = E_gev*GEV\n\n'
             'p = {\n'
             '    "Vacuum": oscprob3nu.probabilities_3nu(\n'
             '        np.asarray(H_VAC_3NU)/E[:, None, None], L),\n'
             '    "Matter": oscprob3nu.probabilities_3nu(\n'
             '        hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '            H_VAC_3NU, E, gd.VCC_EARTH_CRUST), L),\n'
             '    "NSI": oscprob3nu.probabilities_3nu(\n'
             '        hamiltonians3nu.hamiltonian_3nu_nsi(\n'
             '            H_VAC_3NU, E, gd.VCC_EARTH_CRUST, gd.EPS_3), L),\n'
             '    "CPT-odd LIV": oscprob3nu.probabilities_3nu(\n'
             '        hamiltonians3nu.hamiltonian_3nu_liv(\n'
             '            H_VAC_3NU, E, gd.SXI12, gd.SXI23, gd.SXI13,\n'
             '            gd.DXICP, gd.B1, gd.B2, gd.B3, gd.LAMBDA), L),\n'
             '}\n\n'
             'styles = ["-", "--", ":", "-."]\n'
             'panels = [(0, r"$P_{\\nu_e \\to \\nu_e}$"),\n'
             '          (3, r"$P_{\\nu_\\mu \\to \\nu_e}$"),\n'
             '          (4, r"$P_{\\nu_\\mu \\to \\nu_\\mu}$")]\n\n'
             'fig, axes = plt.subplots(3, 1, figsize=(7.2, 9.6), '
             'sharex=True)\n'
             'for ax, (k, ylabel) in zip(axes, panels):\n'
             '    for (name, prob), ls in zip(p.items(), styles):\n'
             '        ax.semilogx(E_gev, prob[:, k], ls, lw=1.8, '
             'label=name)\n'
             '    ax.set_ylabel(ylabel)\n'
             '    ax.set_xlim(E_gev[0], E_gev[-1])\n'
             '    ax.set_ylim(0.0, 1.0)\n'
             'axes[0].legend(loc="lower left", ncol=2)\n'
             'axes[-1].set_xlabel("Neutrino energy [GeV]")\n'
             'fig.tight_layout()\n'
             'plt.show()'),
        md('## Two flavors\n\n'
           'The same four scenarios in the two-flavor system, driven by '
           '$\\Delta m^2_{31}$ and $\\theta_{23}$, over a slightly wider '
           'energy range.\n\n'
           'Note the first argument: the builders take '
           '$\\sin\\theta$, not $\\theta$.'),
        code('E_gev = np.logspace(np.log10(0.5), np.log10(40.0), 400)\n'
             'E = E_gev*GEV\n\n'
             'H2_VAC = hamiltonians2nu.'
             'hamiltonian_2nu_vacuum_energy_independent(\n'
             '    gd.S23_NO_BF, gd.D31_NO_BF)\n\n'
             'p2 = {\n'
             '    "Vacuum": oscprob2nu.probabilities_2nu(\n'
             '        np.asarray(H2_VAC)/E[:, None, None], L),\n'
             '    "Matter": oscprob2nu.probabilities_2nu(\n'
             '        hamiltonians2nu.hamiltonian_2nu_matter(\n'
             '            H2_VAC, E, gd.VCC_EARTH_CRUST), L),\n'
             '    "NSI": oscprob2nu.probabilities_2nu(\n'
             '        hamiltonians2nu.hamiltonian_2nu_nsi(\n'
             '            H2_VAC, E, gd.VCC_EARTH_CRUST, gd.EPS_2), L),\n'
             '    "CPT-odd LIV": oscprob2nu.probabilities_2nu(\n'
             '        hamiltonians2nu.hamiltonian_2nu_liv(\n'
             '            H2_VAC, E, gd.SXI12, gd.B1, gd.B3, gd.LAMBDA), L),\n'
             '}\n\n'
             'panels2 = [(0, r"$P_{\\nu_e \\to \\nu_e}$"),\n'
             '           (2, r"$P_{\\nu_\\mu \\to \\nu_e}$"),\n'
             '           (3, r"$P_{\\nu_\\mu \\to \\nu_\\mu}$")]\n\n'
             'fig, axes = plt.subplots(3, 1, figsize=(7.2, 9.6), '
             'sharex=True)\n'
             'for ax, (k, ylabel) in zip(axes, panels2):\n'
             '    for (name, prob), ls in zip(p2.items(), styles):\n'
             '        ax.semilogx(E_gev, prob[:, k], ls, lw=1.8, '
             'label=name)\n'
             '    ax.set_ylabel(ylabel)\n'
             '    ax.set_xlim(E_gev[0], E_gev[-1])\n'
             '    ax.set_ylim(0.0, 1.0)\n'
             'axes[0].legend(loc="lower left", ncol=2)\n'
             'axes[-1].set_xlabel("Neutrino energy [GeV]")\n'
             'fig.tight_layout()\n'
             'plt.show()'),
    ])

# --------------------------------------------- 11 exact vs the approximations
books['11_exact_vs_approximations.ipynb'] = notebook(
    'Exact versus the textbook approximations',
    'The selling point of this library is that it is exact. This notebook is '
    'the case for caring: where the familiar approximations agree with the '
    'exact answer, and where they do not.',
    [
        md('## The two-flavor vacuum formula is not an approximation\n\n'
           'Worth establishing first, because it sets the baseline. In vacuum '
           'with two flavors, $P_{e\\mu} = \\sin^2 2\\theta \\, '
           '\\sin^2(\\Delta m^2 L/4E)$ is *exact*, and the SU(2) expansion '
           'reproduces it to round-off.'),
        code('s12 = np.sqrt(0.310)\n'
             'H2 = np.asarray(hamiltonians2nu.'
             'hamiltonian_2nu_vacuum_energy_independent(\n'
             '    s12, gd.D21_NO_BF))/(1.0*GEV)\n'
             'L_km = np.logspace(0.0, 5.0, 4000)\n\n'
             'exact = oscprob2nu.probabilities_2nu(H2, L_km*KM)[:, 1]\n'
             'textbook = (4.0*s12**2.0*(1.0-s12**2.0)\n'
             '            * np.sin(gd.D21_NO_BF*L_km*KM/(4.0*GEV))**2.0)\n\n'
             'print("max |exact - textbook| = %.2e" % np.max(np.abs(\n'
             '    exact-textbook)))'),
        md('## Three flavors in vacuum: still exact\n\n'
           '`hamiltonians3nu.probabilities_3nu_vacuum_std` evaluates the '
           'standard sum over mass eigenstates. In vacuum this is also exact, '
           'and the two agree — which is what makes it a useful cross-check '
           'of the SU(3) machinery rather than a competitor to it.'),
        code('U = hamiltonians3nu.pmns_mixing_matrix(\n'
             '    gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, gd.DCP_NO_BF)\n\n'
             'E_gev = np.logspace(-1.0, 1.5, 300)\n'
             'L = 1300.0*KM\n\n'
             'exact = oscprob3nu.probabilities_3nu(\n'
             '    np.asarray(H_VAC_3NU)/(E_gev[:, None, None]*GEV), L)\n'
             'std = np.array([hamiltonians3nu.probabilities_3nu_vacuum_std(\n'
             '    U, gd.D21_NO_BF, gd.D31_NO_BF, e*GEV, L) for e in E_gev])\n\n'
             'print("max |exact - standard formula| = %.2e"\n'
             '      % np.max(np.abs(exact-std)))'),
        md('## Where it does break: constant density for a varying profile\n\n'
           'The approximation that actually costs you is replacing a varying '
           'density by a constant. Below, the exact PREM calculation against '
           'the common shortcut of using a single average density along the '
           'chord.'),
        code('import earth\n\n'
             'costhz = -0.9\n'
             'widths, densities = earth.earth_slabs(costhz, '
             'n_slabs_per_segment=8)\n'
             'mean_density = np.average(densities, weights=widths)\n'
             'total_km = widths.sum()\n\n'
             'print("chord           : %.0f km" % total_km)\n'
             'print("mean density    : %.3f g/cm^3" % mean_density)\n'
             'print("range along it  : %.3f to %.3f g/cm^3"\n'
             '      % (densities.min(), densities.max()))'),
        code('E_gev = np.logspace(0.0, 2.0, 160)\n\n'
             'p_prem, p_const = [], []\n'
             'for e in E_gev:\n'
             '    p_prem.append(earth.probabilities_3nu_earth(\n'
             '        H_VAC_3NU, e*GEV, costhz, n_slabs_per_segment=8)[4])\n'
             '    H = hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '        H_VAC_3NU, e*GEV, earth.matter_potential(mean_density))\n'
             '    p_const.append(oscprob3nu.probabilities_3nu(\n'
             '        H, total_km*KM)[4])\n'
             'p_prem, p_const = np.array(p_prem), np.array(p_const)\n\n'
             'fig, ax = plt.subplots()\n'
             'ax.semilogx(E_gev, p_prem, label="PREM, exact")\n'
             'ax.semilogx(E_gev, p_const, "--",\n'
             '            label="one average density")\n'
             'ax.set_xlim(E_gev[0], E_gev[-1])\n'
             'ax.set_xlabel("Energy [GeV]")\n'
             'ax.set_ylabel(r"$P_{\\mu\\mu}$")\n'
             'ax.set_title(r"$\\cos\\theta_z = -0.9$, through the Earth")\n'
             'ax.legend()\n'
             'plt.show()\n\n'
             'print("max difference : %.3f" % np.max(np.abs(p_prem-p_const)))\n'
             'print("mean difference: %.3f" % np.mean(np.abs('
             'p_prem-p_const)))'),
        md('## And the two-flavor reduction of a three-flavor problem\n\n'
           'Appearance is often estimated with a two-flavor formula driven by '
           'the atmospheric parameters. Against the full three-flavor result '
           'that neglects the solar term and CP violation entirely.'),
        code('E_gev = np.logspace(-0.5, 1.2, 400)\n'
             'L = 1300.0*KM\n\n'
             'p3 = oscprob3nu.probabilities_3nu(\n'
             '    np.asarray(H_VAC_3NU)/(E_gev[:, None, None]*GEV), L)[:, 3]\n\n'
             '# The usual two-flavor estimate for nu_mu -> nu_e\n'
             'sin2_2th13 = 4.0*gd.S13_NO_BF**2.0*(1.0-gd.S13_NO_BF**2.0)\n'
             's23sq = gd.S23_NO_BF**2.0\n'
             'p2 = (s23sq*sin2_2th13\n'
             '      * np.sin(gd.D31_NO_BF*L/(4.0*E_gev*GEV))**2.0)\n\n'
             'fig, ax = plt.subplots()\n'
             'ax.semilogx(E_gev, p3, label="exact, three flavors")\n'
             'ax.semilogx(E_gev, p2, "--",\n'
             '            label=r"two-flavor estimate")\n'
             'ax.set_xlim(E_gev[0], E_gev[-1])\n'
             'ax.set_xlabel("Energy [GeV]")\n'
             'ax.set_ylabel(r"$P_{\\mu e}$")\n'
             'ax.set_title("Appearance at L = 1300 km, in vacuum")\n'
             'ax.legend()\n'
             'plt.show()\n\n'
             'print("max difference near the first maximum: %.4f"\n'
             '      % np.max(np.abs(p3-p2)[E_gev > 1.0]))'),
        md('The gap is the solar term and the CP-violating interference — '
           'precisely the physics a long-baseline experiment is built to '
           'measure. An approximation that discards it is not a small error '
           'in the quantity of interest; it is the quantity of interest.'),
    ])

# ------------------------------------------------ 12 mass ordering and octant
books['12_ordering_and_octant.ipynb'] = notebook(
    'Mass ordering and the $\\theta_{23}$ octant',
    'Two open questions, and how they show up in the probabilities. '
    '`globaldefs` carries the NuFit best fit for both orderings, so both are '
    'available without typing any numbers.',
    [
        code('import earth\n\n\n'
             'def h_vacuum(ordering="NO", s23=None, dcp=None):\n'
             '    """Vacuum Hamiltonian for either ordering, with overrides."""\n'
             '    if ordering == "NO":\n'
             '        s12, s23d, s13 = gd.S12_NO_BF, gd.S23_NO_BF, '
             'gd.S13_NO_BF\n'
             '        dcpd, d21, d31 = gd.DCP_NO_BF, gd.D21_NO_BF, '
             'gd.D31_NO_BF\n'
             '    else:\n'
             '        s12, s23d, s13 = gd.S12_IO_BF, gd.S23_IO_BF, '
             'gd.S13_IO_BF\n'
             '        dcpd, d21, d31 = gd.DCP_IO_BF, gd.D21_IO_BF, '
             'gd.D31_IO_BF\n'
             '    return hamiltonians3nu.'
             'hamiltonian_3nu_vacuum_energy_independent(\n'
             '        s12, s23 if s23 is not None else s23d, s13,\n'
             '        dcp if dcp is not None else dcpd, d21, d31)\n\n\n'
             'print("normal   : Dm31 = %+.3e eV^2" % gd.D31_NO_BF)\n'
             'print("inverted : Dm31 = %+.3e eV^2" % gd.D31_IO_BF)'),
        md('## In vacuum the ordering barely shows\n\n'
           'The sign of $\\Delta m^2_{31}$ enters the vacuum probabilities '
           'only through the small solar and CP-violating terms.'),
        code('E_gev = np.logspace(-0.5, 1.2, 500)\n'
             'L = 1300.0*KM\n\n'
             'fig, ax = plt.subplots()\n'
             'for ordering, ls in (("NO", "-"), ("IO", "--")):\n'
             '    p = oscprob3nu.probabilities_3nu(\n'
             '        np.asarray(h_vacuum(ordering))'
             '/(E_gev[:, None, None]*GEV), L)\n'
             '    ax.semilogx(E_gev, p[:, 3], ls, label=ordering)\n'
             'ax.set_xlim(E_gev[0], E_gev[-1])\n'
             'ax.set_xlabel("Energy [GeV]")\n'
             'ax.set_ylabel(r"$P_{\\mu e}$")\n'
             'ax.set_title("Vacuum, L = 1300 km")\n'
             'ax.legend()\n'
             'plt.show()'),
        md('## Matter is what separates them\n\n'
           'The matter potential enters with a definite sign, so it enhances '
           'the resonance for one ordering and suppresses it for the other. '
           'This is the whole basis of the experimental programme.'),
        code('fig, ax = plt.subplots()\n'
             'for ordering, ls in (("NO", "-"), ("IO", "--")):\n'
             '    H = hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '        h_vacuum(ordering), E_gev*GEV, gd.VCC_EARTH_CRUST)\n'
             '    ax.semilogx(E_gev, oscprob3nu.probabilities_3nu(H, L)'
             '[:, 3],\n'
             '                ls, label=ordering)\n'
             'ax.set_xlim(E_gev[0], E_gev[-1])\n'
             'ax.set_xlabel("Energy [GeV]")\n'
             'ax.set_ylabel(r"$P_{\\mu e}$")\n'
             'ax.set_title("In matter, L = 1300 km")\n'
             'ax.legend()\n'
             'plt.show()'),
        md('## Through the Earth the separation is larger still\n\n'
           'A trajectory through the mantle sees much more matter than a '
           '1300 km chord through the crust.'),
        code('costhz = np.linspace(-0.999, -0.05, 140)\n'
             'E_fixed = 8.0*GEV\n\n'
             'fig, ax = plt.subplots()\n'
             'for ordering, ls in (("NO", "-"), ("IO", "--")):\n'
             '    hv = h_vacuum(ordering)\n'
             '    p = [earth.probabilities_3nu_earth(\n'
             '        hv, E_fixed, c, n_slabs_per_segment=5)[3] '
             'for c in costhz]\n'
             '    ax.plot(costhz, p, ls, label=ordering)\n'
             'ax.set_xlim(costhz[0], costhz[-1])\n'
             'ax.set_xlabel(r"$\\cos\\theta_z$")\n'
             'ax.set_ylabel(r"$P_{\\mu e}$")\n'
             'ax.set_title("Through the Earth at E = 8 GeV")\n'
             'ax.legend()\n'
             'plt.show()'),
        md('## The octant\n\n'
           'Disappearance depends on $\\theta_{23}$ mainly through '
           '$\\sin^2 2\\theta_{23}$, which is symmetric about '
           '$\\theta_{23} = \\pi/4$. Two values either side of maximal are '
           'therefore nearly indistinguishable in $P_{\\mu\\mu}$ — the '
           'octant degeneracy.'),
        code('s23_lower = np.sqrt(0.45)      # below maximal\n'
             's23_upper = np.sqrt(0.55)      # above maximal\n\n'
             'fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.0), sharex=True)\n'
             'for s23, ls, lab in ((s23_lower, "-", r"$\\sin^2\\theta_{23} = '
             '0.45$"),\n'
             '                     (s23_upper, "--", r"$\\sin^2\\theta_{23} = '
             '0.55$")):\n'
             '    H = hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '        h_vacuum("NO", s23=s23), E_gev*GEV, gd.VCC_EARTH_CRUST)\n'
             '    p = oscprob3nu.probabilities_3nu(H, L)\n'
             '    axes[0].semilogx(E_gev, p[:, 4], ls, label=lab)\n'
             '    axes[1].semilogx(E_gev, p[:, 3], ls, label=lab)\n'
             'axes[0].set_ylabel(r"$P_{\\mu\\mu}$  (disappearance)")\n'
             'axes[1].set_ylabel(r"$P_{\\mu e}$  (appearance)")\n'
             'axes[1].set_xlabel("Energy [GeV]")\n'
             'for ax in axes:\n'
             '    ax.set_xlim(E_gev[0], E_gev[-1])\n'
             'axes[0].legend()\n'
             'fig.tight_layout()\n'
             'plt.show()'),
        md('Disappearance can hardly tell the two apart; appearance can, '
           'because it carries a factor of $\\sin^2\\theta_{23}$ rather than '
           '$\\sin^2 2\\theta_{23}$. That asymmetry is why resolving the '
           'octant needs the appearance channel.'),
    ])

# ---------------------------------------------------------- 13 antineutrinos
books['13_antineutrinos.ipynb'] = notebook(
    'Antineutrinos, done properly',
    'Getting antineutrinos right is the most common way to get a matter '
    'calculation wrong. Two things change together, and changing only one is '
    'silently plausible.',
    [
        md('## The rule\n\n'
           'For antineutrinos the Hamiltonian is the **complex conjugate** of '
           'the neutrino one *and* the matter potential changes sign:\n\n'
           '$$ H_{\\bar\\nu} = H_\\nu^* \\big|_{V \\to -V} $$\n\n'
           'The conjugation flips the sign of $\\delta_{CP}$; the potential '
           'flip is because matter contains electrons and not positrons. In '
           'practice: conjugate the *vacuum* Hamiltonian, then add the '
           'potential with a minus sign.'),
        code('def h_antineutrino_matter(h_vac, energy, vcc):\n'
             '    """Antineutrino Hamiltonian in matter."""\n'
             '    return hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '        np.conj(h_vac), energy, -vcc)\n\n\n'
             'E_gev = np.logspace(-0.5, 1.2, 500)\n'
             'E = E_gev*GEV\n'
             'L = 1300.0*KM\n\n'
             'p_nu = oscprob3nu.probabilities_3nu(\n'
             '    hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '        H_VAC_3NU, E, gd.VCC_EARTH_CRUST), L)\n'
             'p_bar = oscprob3nu.probabilities_3nu(\n'
             '    h_antineutrino_matter(H_VAC_3NU, E, gd.VCC_EARTH_CRUST), '
             'L)'),
        code('fig, ax = plt.subplots()\n'
             'ax.semilogx(E_gev, p_nu[:, 3], label=r"$\\nu_\\mu \\to '
             '\\nu_e$")\n'
             'ax.semilogx(E_gev, p_bar[:, 3], "--",\n'
             '            label=r"$\\bar\\nu_\\mu \\to \\bar\\nu_e$")\n'
             'ax.set_xlim(E_gev[0], E_gev[-1])\n'
             'ax.set_xlabel("Energy [GeV]")\n'
             'ax.set_ylabel("Appearance probability")\n'
             'ax.set_title("Matter separates neutrinos from antineutrinos")\n'
             'ax.legend()\n'
             'plt.show()'),
        md('## Getting it half right\n\n'
           'The two common mistakes: flipping only the potential, or '
           'conjugating only the Hamiltonian. Both produce a smooth, '
           'plausible curve that is simply wrong.'),
        code('only_potential = oscprob3nu.probabilities_3nu(\n'
             '    hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '        H_VAC_3NU, E, -gd.VCC_EARTH_CRUST), L)\n'
             'only_conjugate = oscprob3nu.probabilities_3nu(\n'
             '    hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '        np.conj(H_VAC_3NU), E, gd.VCC_EARTH_CRUST), L)\n\n'
             'fig, ax = plt.subplots()\n'
             'ax.semilogx(E_gev, p_bar[:, 3], lw=2.0, label="correct")\n'
             'ax.semilogx(E_gev, only_potential[:, 3], "--",\n'
             '            label="potential flipped only")\n'
             'ax.semilogx(E_gev, only_conjugate[:, 3], ":",\n'
             '            label="conjugated only")\n'
             'ax.set_xlim(E_gev[0], E_gev[-1])\n'
             'ax.set_xlabel("Energy [GeV]")\n'
             'ax.set_ylabel(r"$P_{\\bar\\mu\\bar e}$")\n'
             'ax.set_title("Two ways to be wrong")\n'
             'ax.legend()\n'
             'plt.show()\n\n'
             'print("max error, potential only : %.4f"\n'
             '      % np.max(np.abs(p_bar[:, 3]-only_potential[:, 3])))\n'
             'print("max error, conjugate only : %.4f"\n'
             '      % np.max(np.abs(p_bar[:, 3]-only_conjugate[:, 3])))'),
        md('## Checks that must hold\n\n'
           'Two limits pin the treatment down. In **vacuum**, CPT requires '
           '$P(\\nu_\\alpha \\to \\nu_\\beta) = '
           'P(\\bar\\nu_\\beta \\to \\bar\\nu_\\alpha)$ — the transpose of '
           'the probability matrix. And with **$\\delta_{CP} = 0$ and no '
           'matter**, neutrinos and antineutrinos are identical.'),
        code('# CPT in vacuum: P_ab(nu) = P_ba(nubar)\n'
             'e_one = 2.0*GEV\n'
             'p_v = np.array(oscprob3nu.probabilities_3nu(\n'
             '    np.asarray(H_VAC_3NU)/e_one, L)).reshape(3, 3)\n'
             'p_v_bar = np.array(oscprob3nu.probabilities_3nu(\n'
             '    np.conj(np.asarray(H_VAC_3NU))/e_one, L)).reshape(3, 3)\n\n'
             'print("max |P(nu) - P(nubar)^T| in vacuum = %.2e"\n'
             '      % np.max(np.abs(p_v - p_v_bar.T)))'),
        code('# delta_CP = 0 and no matter: the two are the same\n'
             'h_nocp = hamiltonians3nu.'
             'hamiltonian_3nu_vacuum_energy_independent(\n'
             '    gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, 0.0,\n'
             '    gd.D21_NO_BF, gd.D31_NO_BF)\n\n'
             'a = np.array(oscprob3nu.probabilities_3nu(\n'
             '    np.asarray(h_nocp)/e_one, L))\n'
             'b = np.array(oscprob3nu.probabilities_3nu(\n'
             '    np.conj(np.asarray(h_nocp))/e_one, L))\n'
             'print("max |nu - nubar| with dCP = 0, vacuum = %.2e"\n'
             '      % np.max(np.abs(a-b)))'),
        md('## Through the Earth\n\n'
           'The same rule applies to a slabbed calculation: conjugate the '
           'vacuum Hamiltonian and pass the negative potential. Because '
           '`earth` builds its Hamiltonians from the vacuum one it is given, '
           'the antineutrino case is one conjugation away — but the sign of '
           'the potential has to be handled where the potential is built, so '
           'the profile is assembled by hand here.'),
        code('import earth\nimport slabs\n\n\n'
             'def probabilities_earth_antinu(h_vac, energy, costhz,\n'
             '                               n_slabs_per_segment=6):\n'
             '    """Antineutrino probabilities across the Earth."""\n'
             '    widths_km, densities = earth.earth_slabs(\n'
             '        costhz, n_slabs_per_segment)\n'
             '    vcc = -earth.matter_potential(densities)     # antineutrinos\n'
             '    H = hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '        np.conj(h_vac), energy, vcc)             # antineutrinos\n'
             '    return np.array(slabs.probabilities_3nu_slabs(\n'
             '        H, widths_km*KM))\n\n\n'
             'costhz = np.linspace(-0.999, -0.05, 140)\n'
             'E_fixed = 8.0*GEV\n\n'
             'p_nu = [earth.probabilities_3nu_earth(\n'
             '    H_VAC_3NU, E_fixed, c, n_slabs_per_segment=6)[3]\n'
             '    for c in costhz]\n'
             'p_bar = [probabilities_earth_antinu(H_VAC_3NU, E_fixed, c)[3]\n'
             '         for c in costhz]\n\n'
             'fig, ax = plt.subplots()\n'
             'ax.plot(costhz, p_nu, label="neutrinos")\n'
             'ax.plot(costhz, p_bar, "--", label="antineutrinos")\n'
             'ax.set_xlim(costhz[0], costhz[-1])\n'
             'ax.set_xlabel(r"$\\cos\\theta_z$")\n'
             'ax.set_ylabel(r"$P_{\\mu e}$")\n'
             'ax.set_title("Through the Earth at E = 8 GeV")\n'
             'ax.legend()\n'
             'plt.show()'),
    ])

# -------------------------------------------------------- 14 solar neutrinos
books['14_solar_and_adiabatic_msw.ipynb'] = notebook(
    'Solar neutrinos, the MSW resonance, and the limits of slabs',
    'The adiabatic MSW resonance in the Sun is the textbook application of '
    'matter effects. This notebook does it with slabs, gets the right answer, '
    'and then explains why you should probably not do it this way.\n\n'
    '**The short version:** the method works and is validated here against '
    'the analytic adiabatic result, but the cost of a realistic solar '
    'calculation makes a Magnus-type method the better tool. The last section '
    'gives the numbers.',
    [
        code('import slabs\n'
             'import earth   # matter_potential is not Earth-specific'),
        md('## A solar density profile\n\n'
           'The density falls roughly exponentially through the radiative '
           'zone: $\\rho(r) \\simeq \\rho_0 e^{-r/r_0}$ with '
           '$\\rho_0 \\approx 150$ g/cm$^3$ and '
           '$r_0 \\approx R_\\odot/10.54$.'),
        code('R_SUN_KM = 6.957e5\n'
             'RHO_CENTRE = 150.0                 # [g cm^-3]\n'
             'R_SCALE = R_SUN_KM/10.54\n\n\n'
             'def density_sun(r_km):\n'
             '    """Exponential fit to the solar density profile."""\n'
             '    return RHO_CENTRE*np.exp(-np.asarray(r_km, dtype=float)'
             '/R_SCALE)\n\n\n'
             'H2_VAC_SOLAR = hamiltonians2nu.'
             'hamiltonian_2nu_vacuum_energy_independent(\n'
             '    gd.S12_NO_BF, gd.D21_NO_BF)\n\n'
             'r = np.linspace(0.0, R_SUN_KM, 2000)\n'
             'fig, ax = plt.subplots()\n'
             'ax.semilogy(r/R_SUN_KM, density_sun(r))\n'
             'ax.set_xlim(0.0, 1.0)\n'
             'ax.set_xlabel(r"$r/R_\\odot$")\n'
             'ax.set_ylabel(r"Density [g cm$^{-3}$]")\n'
             'ax.set_title("Solar density, exponential fit")\n'
             'plt.show()'),
        md('## Where the resonance sits\n\n'
           'The MSW resonance is where the matter term matches the vacuum '
           'one, $2EV = \\Delta m^2 \\cos 2\\theta$. A 10 MeV neutrino is '
           'produced *above* resonance and crosses it on the way out, which '
           'is the condition for MSW conversion.'),
        code('E_mev = 10.0\n'
             'E = E_mev*1.0e6\n'
             'cos2th = 1.0 - 2.0*gd.S12_NO_BF**2.0\n\n'
             'v_res = gd.D21_NO_BF*cos2th/(2.0*E)\n'
             'v_of_r = earth.matter_potential(density_sun(r))\n'
             'k = int(np.argmin(np.abs(v_of_r - v_res)))\n\n'
             'print("V at the centre    : %.3e eV" % v_of_r[0])\n'
             'print("resonant potential : %.3e eV" % v_res)\n'
             'print("crossed at r       : %.3f R_sun" % (r[k]/R_SUN_KM))\n\n'
             'fig, ax = plt.subplots()\n'
             'ax.semilogy(r/R_SUN_KM, v_of_r, label="matter potential")\n'
             'ax.axhline(v_res, color="C1", ls="--",\n'
             '           label="resonance, E = 10 MeV")\n'
             'ax.axvline(r[k]/R_SUN_KM, color="0.6", lw=1)\n'
             'ax.set_xlim(0.0, 1.0)\n'
             'ax.set_xlabel(r"$r/R_\\odot$")\n'
             'ax.set_ylabel(r"$V_{CC}$ [eV]")\n'
             'ax.set_title("Where the MSW resonance sits")\n'
             'ax.legend()\n'
             'plt.show()'),
        md('## Slabbing the Sun\n\n'
           'A slab is exact only if the Hamiltonian is constant across it, so '
           'the slabs must be short compared with the oscillation length. '
           'That sets the count.'),
        code('def probabilities_sun(n_slabs, energy):\n'
             '    """Two-flavor survival through the solar profile."""\n'
             '    edges = np.linspace(0.0, R_SUN_KM, n_slabs+1)\n'
             '    widths = np.diff(edges)\n'
             '    mid = 0.5*(edges[:-1]+edges[1:])\n'
             '    vcc = earth.matter_potential(density_sun(mid))\n'
             '    H = hamiltonians2nu.hamiltonian_2nu_matter(\n'
             '        H2_VAC_SOLAR, energy, vcc)\n'
             '    return np.array(slabs.probabilities_2nu_slabs(H, '
             'widths*KM))\n\n\n'
             'L_osc_km = 4.0*np.pi*E/gd.D21_NO_BF/KM\n'
             'print("oscillation length at 10 MeV : %.1f km" % L_osc_km)\n'
             'print("path to cross                : %.3e km" % R_SUN_KM)\n'
             'print("slabs at ~10 per oscillation : %.1e"\n'
             '      % (10.0*R_SUN_KM/L_osc_km))\n\n'
             'for n in (100, 300, 1000, 3000, 10000, 30000):\n'
             '    print("%7d slabs -> P_ee = %.6f" % (n, '
             'probabilities_sun(n, E)[0]))'),
        md('It converges, at around ten thousand slabs — exactly the estimate '
           'above. So far so good.\n\n'
           'And then the trap.'),
        md('## The converged number is meaningless on its own\n\n'
           'That $P_{ee}$ is the survival probability for one energy, at one '
           'production point, evaluated at one radius. Move the energy by a '
           'tenth of a percent and it changes completely.'),
        code('for e in np.linspace(0.999e7, 1.001e7, 7):\n'
             '    print("E = %.4f MeV -> P_ee = %.4f"\n'
             '          % (e/1.0e6, probabilities_sun(8000, e)[0]))'),
        md('The probability is oscillating far faster than any detector '
           'resolves, and faster than the production region is localised. '
           'The physical quantity is the **average**, and nothing about a '
           'single evaluation approximates it.'),
        md('## Averaging recovers the textbook result\n\n'
           'Averaged over a narrow energy band, the slab calculation should '
           'reproduce the analytic adiabatic MSW prediction\n\n'
           '$$ \\langle P_{ee} \\rangle = \\tfrac{1}{2}\\left(1 + '
           '\\cos 2\\theta_m^0 \\cos 2\\theta\\right), $$\n\n'
           'with $\\theta_m^0$ the mixing angle at the production point. '
           'This is the check that the whole machinery is right.'),
        code('x = 2.0*E*v_of_r[0]/gd.D21_NO_BF\n'
             'cos2th_m = (cos2th-x)/np.sqrt((cos2th-x)**2.0 + '
             '(1.0-cos2th**2.0))\n'
             'p_adiabatic = 0.5*(1.0 + cos2th_m*cos2th)\n\n'
             'energies = np.linspace(0.99e7, 1.01e7, 200)\n'
             'p_avg = np.mean([probabilities_sun(4000, e)[0] '
             'for e in energies])\n\n'
             'print("adiabatic MSW prediction : %.4f" % p_adiabatic)\n'
             'print("slabs, averaged          : %.4f" % p_avg)\n'
             'print("difference               : %.4f"\n'
             '      % abs(p_avg-p_adiabatic))\n'
             'print()\n'
             'print("for contrast, no MSW at all (vacuum-averaged): %.4f"\n'
             '      % (1.0 - 0.5*(1.0-cos2th**2.0)))'),
        md('It agrees to a fraction of a percent, and is clearly distinct '
           'from the vacuum-averaged value — the conversion is real and the '
           'slab calculation captures it.'),
        md('## Why you should still not do it this way\n\n'
           'Count what that took. One averaged point needed **4000 slabs x '
           '200 energies = 800 000 slab evaluations**, and that is the '
           'cheapest possible version of the problem: two flavors, one '
           'production radius, one narrow energy band, an analytic density '
           'profile.\n\n'
           'A real solar calculation needs, on top of that:\n\n'
           '* three flavors rather than two;\n'
           '* an average over the production region, which spans a range of '
           'radii for each reaction and differs between the pp, $^7$Be and '
           '$^8$B components;\n'
           '* a spectrum rather than a single energy;\n'
           '* the day-night effect, which adds an Earth crossing on top.\n\n'
           'Each of those is a multiplier. The cost grows as the path length '
           'divided by the oscillation length — and the oscillation length '
           'shrinks with energy, so the problem gets worse exactly where the '
           'interesting physics is.\n\n'
           'None of this is a defect in the method. Every slab is still '
           'exact. The difficulty is that a piecewise-constant approximation '
           'is the wrong shape for a profile that varies smoothly over many '
           'oscillation lengths: it needs its step size set by the '
           '*oscillation*, not by the *variation of the Hamiltonian*, and '
           'those differ here by four orders of magnitude.\n\n'
           '### Use a Magnus expansion instead\n\n'
           'A Magnus-type method integrates the varying Hamiltonian directly '
           'rather than freezing it in pieces, so its step size is set by how '
           'fast the Hamiltonian changes, and it converges at high order in '
           'that step. For a smooth profile crossed over many oscillations '
           'this is the difference between $10^4$ steps and a few dozen.\n\n'
           'The sibling package '
           '[Magnus](https://github.com/mbustama/Magnus) does exactly this, '
           'and shares this library\'s conventions and its set of named '
           'locations.\n\n'
           '**Rule of thumb.** Use NuOscProbExact where the Hamiltonian is '
           'constant or genuinely piecewise-constant: a beam through the '
           'crust, a trajectory through the Earth, a castle wall. Use a '
           'Magnus-type method where it varies smoothly on the scale of an '
           'oscillation length, as it does in the Sun.'),
    ])

books['15_numerical_edge_cases.ipynb'] = notebook(
    'Numerical edge cases',
    'The expansions divide by quantities that can vanish. This notebook '
    'walks the cases where a naive implementation returns NaN, and shows what '
    'this one returns instead — every one of them was a real defect fixed in '
    'the 1.1.0 audit, and every one has a regression test.',
    [
        md('## A Hamiltonian proportional to the identity\n\n'
           'If $H \\propto \\mathbb{1}$ the traceless part vanishes, so '
           '$|h| = 0$ and the SU(3) expansion divides by zero. The correct '
           'answer is that nothing oscillates: $U = \\mathbb{1}$.'),
        code('H = np.eye(3, dtype=complex)*3.7\n'
             'prob = np.array(oscprob3nu.probabilities_3nu(H, 5.0))\n'
             'print("P = ", np.round(prob, 12))\n'
             'print("no NaN:", not np.any(np.isnan(prob)))\n'
             'print("identity:", np.allclose(prob.reshape(3, 3), '
             'np.eye(3)))'),
        md('## A doubly degenerate spectrum\n\n'
           '`diag(1, 1, -2)` has two equal latent roots. The Lagrange '
           'interpolation used for the coefficients divides by '
           '$3\\psi_m^2 - |h|^2$, which vanishes exactly at a repeated root; '
           'the confluent limit is taken instead.'),
        code('from scipy.linalg import expm\n\n'
             'H = np.diag([1.0, 1.0, -2.0]).astype(complex)\n'
             'L = 0.9\n\n'
             'U = np.array(oscprob3nu.evolution_operator_3nu(H, L))\n'
             'U_ref = expm(-1.0j*(H - np.trace(H)/3.0*np.eye(3))*L)\n\n'
             'print("max |U - expm| = %.2e" % np.max(np.abs(U-U_ref)))\n'
             'print("unitary to     = %.2e"\n'
             '      % np.max(np.abs(U.conj().T @ U - np.eye(3))))'),
        md('## Approaching degeneracy\n\n'
           'The dangerous region is not the degenerate point itself but its '
           'neighbourhood, where the denominator is small but not zero. '
           'Sweeping a splitting down to $10^{-15}$ against an independent matrix '
           'exponential shows what actually happens: the result stays finite '
           'and unitary throughout, but accuracy is *not* round-off across '
           'the whole range. In the crossover band, where the denominator is '
           'small but the degenerate branch has not yet been taken, the error '
           'peaks around $10^{-8}$ before dropping back. That is worth knowing '
           'and is reported below rather than glossed.'),
        code('errors, splittings = [], np.logspace(-1.0, -15.0, 60)\n'
             'for eps in splittings:\n'
             '    H = np.diag([1.0, 1.0+eps, -2.0-eps]).astype(complex)\n'
             '    U = np.array(oscprob3nu.evolution_operator_3nu(H, 0.9))\n'
             '    U_ref = expm(-1.0j*(H-np.trace(H)/3.0*np.eye(3))*0.9)\n'
             '    errors.append(np.max(np.abs(U-U_ref)))\n\n'
             'fig, ax = plt.subplots()\n'
             'ax.loglog(splittings, errors, "o-", ms=3)\n'
             'ax.set_xlim(splittings[-1], splittings[0])\n'
             'ax.set_xlabel("Splitting between the two degenerate roots")\n'
             'ax.set_ylabel(r"max $|U - e^{-iH_0L}|$")\n'
             'ax.set_title("Stability through the degenerate limit")\n'
             'plt.show()\n\n'
             'print("worst error across the sweep : %.2e" % max(errors))\n'
             'print("error at the tightest splitting: %.2e" % errors[-1])\n'
             'worst = splittings[int(np.argmax(errors))]\n'
             'print("worst at a splitting of        : %.1e" % worst)'),
        md('## Two flavors: the same limit\n\n'
           'The SU(2) expansion carries $\\sin(|h|L)/|h|$, which needs its '
           'limit $\\to L$ at $|h| = 0$.'),
        code('H2 = np.eye(2, dtype=complex)*2.5\n'
             'p = np.array(oscprob2nu.probabilities_2nu(H2, 3.0))\n'
             'print("P =", np.round(p, 12))\n'
             'print("no NaN:", not np.any(np.isnan(p)))'),
        md('## Index validation\n\n'
           '`tensor_d` used to fall off the end of an `elif` chain and return '
           '`None` for an out-of-range index, so the failure surfaced far '
           'from its cause as a `TypeError` in whatever consumed it. It now '
           'raises `IndexError` at the point of the mistake.'),
        code('print("d_007 = %.6f" % oscprob3nu.tensor_d(0, 0, 7))\n'
             'try:\n'
             '    oscprob3nu.tensor_d(0, 1, 99)\n'
             'except IndexError as error:\n'
             '    print("out of range ->", type(error).__name__, error)'),
        md('## Very long and very short baselines\n\n'
           'Nothing special happens at either extreme, but it is worth '
           'confirming: at $L = 0$ the operator is the identity, and at very '
           'large $HL$ the probabilities stay bounded and normalised rather '
           'than accumulating round-off.'),
        code('rng = np.random.default_rng(1)\n'
             'a = rng.normal(size=(3, 3)) + 1.0j*rng.normal(size=(3, 3))\n'
             'H = (a + a.conj().T)/2.0\n\n'
             'U0 = np.array(oscprob3nu.evolution_operator_3nu(H, 0.0))\n'
             'print("U(0) = identity:", np.allclose(U0, np.eye(3)))\n\n'
             'for L in (1.0e3, 1.0e6, 1.0e9):\n'
             '    p = np.array(oscprob3nu.probabilities_3nu(H, L))\n'
             '    print("L = %.0e : sum over final flavors = %.12f, "\n'
             '          "min = %.2e" % (L, p[:3].sum(), p.min()))'),
        md('Every case above is covered by `tests/test_edge_cases.py`, which '
           'is what keeps them from regressing. They are collected here '
           'because a user meeting a NaN wants to know whether it is their '
           'Hamiltonian or the library — and the answer, for all of these, is '
           'that it is neither.'),
    ])

# ------------------------------------------------------- 16 four neutrinos
books['16_four_neutrinos.ipynb'] = notebook(
    'Four neutrinos, and a sterile state',
    'Everything so far has been two or three flavors. `oscprob4nu` carries '
    'the same closed-form treatment to **four**, through the SU(4) algebra — '
    'and four is where the road ends, for a reason worth knowing.\n\n'
    'The payoff is 3+1 sterile scenarios. A 3+1 system is often called '
    '"leaky", because probability drains out of the three active flavors; but '
    'that is a statement about the $3\\times3$ block, not about the physics. '
    'Over all four states the evolution is closed and unitary, which is '
    'exactly what the method assumes.',
    [
        code('# A 3+1 scenario: the NuFit best fit for the active sector,\n'
             '# plus an eV-scale sterile with sin^2(th14) = sin^2(th24) = 0.1\n'
             'DM41 = 1.0                       # [eV^2]\n'
             'S14 = S24 = np.sqrt(0.10)\n'
             'S34 = 0.0\n\n'
             'H_VAC_4NU = '
             'hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(\n'
             '    gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF,\n'
             '    S14, S24, S34,\n'
             '    gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF, DM41)\n\n'
             'print("H_vac shape:", H_VAC_4NU.shape)'),
        md('## One four-flavor probability\n\n'
           'Sixteen probabilities come back, with the *initial* flavor '
           'varying slowest — the same ordering as at three flavors. The '
           'fourth state is the sterile one.'),
        code('energy = 1.0*GEV\n'
             'baseline = 1300.0*KM\n\n'
             'H = hamiltonians4nu.hamiltonian_4nu_matter(\n'
             '    H_VAC_4NU, energy, gd.VCC_EARTH_CRUST, gd.VNC_EARTH_CRUST)\n'
             'prob = np.array(oscprob4nu.probabilities_4nu(H, baseline))\n\n'
             'flavors = ["e", "mu", "tau", "s"]\n'
             'table = prob.reshape(4, 4)\n'
             'print("P(alpha -> beta), rows = initial flavor\\n")\n'
             'print("%-8s" % "" + "".join("%10s" % ("-> " + f) '
             'for f in flavors))\n'
             'for i, a in enumerate(flavors):\n'
             '    print("%-8s" % ("nu_" + a)\n'
             '          + "".join("%10.5f" % v for v in table[i]))\n'
             'print("\\nrows sum to:", np.round(table.sum(axis=1), 12))'),
        md('## The matter potential has a sterile entry\n\n'
           'This is the four-flavor trap, and it is the analogue of getting '
           'antineutrinos wrong: it is **invisible in vacuum** and it moves '
           'the resonance.\n\n'
           'Active neutrinos feel the charged-current potential $V_{CC}$ (the '
           'electron flavor) and the flavor-universal neutral-current '
           'potential $V_{NC}$ (all three). A sterile state feels neither. '
           'Since a term proportional to the identity is only a global phase, '
           'subtracting $V_{NC}\\mathbb{1}$ from all four costs nothing and '
           'leaves\n\n'
           '$$A_4 = \\mathrm{diag}(V_{CC},\\, 0,\\, 0,\\, -V_{NC})$$\n\n'
           'so the sterile entry is *not* zero.'),
        code('added = H - H_VAC_4NU/energy\n'
             'print("diagonal of the matter term [eV]:")\n'
             'for f, v in zip(flavors, np.diag(added).real):\n'
             '    print("  nu_%-4s %+.6e" % (f, v))\n'
             'print("\\n-V_NC / V_CC = %.4f"\n'
             '      % (-gd.VNC_EARTH_CRUST/gd.VCC_EARTH_CRUST))'),
        md('## Sterile oscillations against energy\n\n'
           'An eV-scale splitting oscillates *fast*: the phase goes as '
           '$\\Delta m^2_{41} L / 4E$, so at a long-baseline energy and '
           'distance it has already averaged out and only its mean survives. '
           'That is precisely why short-baseline experiments are the probe — '
           'at $L/E \\sim 1$ km/GeV the oscillation is resolved rather than '
           'averaged.\n\n'
           'So this is a short baseline, 600 m, which is where the structure '
           'actually lives.'),
        code('L_sbl = 0.6*KM                      # 600 m\n'
             'E_gev = np.linspace(0.02, 3.0, 1500)\n\n'
             'H_stack = hamiltonians4nu.hamiltonian_4nu_matter(\n'
             '    H_VAC_4NU, E_gev*GEV, gd.VCC_EARTH_CRUST,\n'
             '    gd.VNC_EARTH_CRUST)\n'
             'p = oscprob4nu.probabilities_4nu(H_stack, L_sbl)\n\n'
             'fig, ax = plt.subplots()\n'
             'ax.plot(E_gev, p[:, 5], label=r"$P_{\\mu\\mu}$")\n'
             'ax.plot(E_gev, p[:, 7], label=r"$P_{\\mu s}$")\n'
             'ax.plot(E_gev, p[:, 0], label=r"$P_{ee}$")\n'
             'ax.set_xlim(E_gev[0], E_gev[-1])\n'
             'ax.set_ylim(0.0, 1.05)\n'
             'ax.set_xlabel("Energy [GeV]")\n'
             'ax.set_ylabel("Probability")\n'
             'ax.set_title(r"3+1 at $L = 600$ m, "\n'
             '             r"$\\Delta m^2_{41} = 1$ eV$^2$, "\n'
             '             r"$\\sin^2\\theta_{14} = \\sin^2\\theta_{24} '
             '= 0.1$")\n'
             'ax.legend(ncol=3, loc="lower right")\n'
             'plt.show()'),
        md('The active flavors are barely oscillating at 600 m through the '
           'standard splittings — those need hundreds of kilometres — so '
           'everything visible here is the sterile state, and $P_{\\mu\\mu}$ '
           'and $P_{ee}$ dip exactly where $P_{\\mu s}$ rises. At '
           '$L = 1300$ km the same curves would be a solid averaged band, '
           'which is the point: the two regimes probe different things.'),
        md('## The sterile matter resonance through the Earth\n\n'
           'The reason four flavors is worth having in a code like this. An '
           'eV-scale sterile driven through Earth matter develops a resonance '
           'in the TeV range for antineutrinos — the signature IceCube '
           'searches for. Here it is as a map in energy and zenith angle, '
           'built the same way as the three-flavor oscillogram: broadcast the '
           'two axes and let one call fill the grid.'),
        code('import earth\n\n'
             'n_e, n_c = 150, 150\n'
             'E_tev = np.logspace(-0.5, 1.5, n_e)          # 0.3 - 30 TeV\n'
             'costhz = np.linspace(-1.0, -0.05, n_c)\n\n'
             '# Antineutrinos: conjugate the vacuum term and flip both\n'
             '# potentials.  Constant mantle density is enough to show the\n'
             '# resonance; `earth` does the full PREM profile.\n'
             'rho_mantle = 4.5                              # [g cm^-3]\n'
             'vcc = earth.matter_potential(rho_mantle)\n'
             'H_bar = hamiltonians4nu.hamiltonian_4nu_matter(\n'
             '    np.conj(H_VAC_4NU), E_tev[:, None]*1.0e3*GEV,\n'
             '    -vcc, -gd.VNC_EARTH_CRUST*rho_mantle/3.0)\n\n'
             '# Path length through the Earth at each zenith angle\n'
             'L_grid = np.array([earth.distance_traveled_inside_earth(c)\n'
             '                   for c in costhz])*KM\n\n'
             'grid = oscprob4nu.probabilities_4nu(H_bar, L_grid[None, :])\n'
             'p_mumu = grid[:, :, 5]\n\n'
             'fig, ax = plt.subplots()\n'
             'mesh = ax.pcolormesh(costhz, E_tev, p_mumu,\n'
             '                     shading="auto", cmap="viridis",\n'
             '                     vmin=0.0, vmax=1.0)\n'
             'ax.set_yscale("log")\n'
             'ax.set_xlabel(r"$\\cos\\theta_z$")\n'
             'ax.set_ylabel("Energy [TeV]")\n'
             'ax.set_title(r"$P_{\\bar\\mu\\bar\\mu}$, 3+1 through the Earth")\n'
             'ax.grid(False)\n'
             'fig.colorbar(mesh, ax=ax, label=r"$P_{\\bar\\mu\\bar\\mu}$")\n'
             'plt.show()'),
        md('The dark band is the sterile matter resonance: a region of '
           'energy and zenith angle where muon antineutrinos are strongly '
           'depleted into the sterile state. Its position depends on '
           '$\\Delta m^2_{41}$, which is what makes it a search channel.'),
        md('## What is actually new at four flavors\n\n'
           'Three things, and each follows from SU(4) having rank three '
           'where SU(3) has rank two.'),
        code('I2, I3, I4 = oscprob4nu.su4_invariants(H)\n'
             'print("three invariants, not two:")\n'
             'print("  I2 = %+.6e" % I2)\n'
             'print("  I3 = %+.6e" % I3)\n'
             'print("  I4 = %+.6e" % I4)\n\n'
             'psi = oscprob4nu.psi_roots_4nu(I2, I3, I4)\n'
             'print("\\nthe characteristic equation is a quartic; its roots:")\n'
             'print(" ", np.array(psi))\n'
             'print("  vs numpy.linalg.eigvalsh:")\n'
             'H_traceless = H - np.trace(H).real/4.0*np.eye(4)\n'
             'print(" ", np.sort(np.linalg.eigvalsh(H_traceless)))'),
        md('And the star-product tower grows a rung. At three flavors '
           '$(h \\star h) \\star h = \\frac{1}{3}|h|^2 h$, which is what lets '
           '`oscprob3nu` stop after two terms. That identity is a '
           'Cayley–Hamilton accident of $n = 3$ and is simply **false** at '
           '$n = 4$:'),
        code('lam = oscprob4nu.LAMBDA_SU4\n'
             'rng = np.random.default_rng(7)\n'
             'h = rng.standard_normal(15)\n'
             'Ht = np.einsum("a,aij->ij", h, lam)\n'
             'i2, _, _ = oscprob4nu.su4_invariants(Ht)\n'
             'tower = (0.5*np.einsum("aij,ji->a", lam, Ht @ Ht @ Ht).real\n'
             '         - 0.5*i2*h)\n'
             'deviation = (np.max(np.abs(tower - i2*h/3.0))\n'
             '             / np.max(np.abs(tower)))\n'
             'print("relative deviation from the SU(3) identity: %.0f%%"\n'
             '      % (100*deviation))\n'
             'print("-> ((h*h)*h)_a is independent data at n = 4")'),
        md('## Accuracy, honestly\n\n'
           'Before the numbers: **none of this is near a measurable '
           'effect.** Oscillation probabilities meet data at the per-cent '
           'level, and both figures below sit several orders of magnitude '
           'beneath anything an experiment resolves. What follows matters '
           'for the library\'s claim to be *exact*, for error accumulating '
           'when `slabs` and `earth` multiply operators across many layers, '
           'and for keeping the regression suite tight enough to catch a '
           'real mistake — not for whether these probabilities are good '
           'enough to use.\n\n'
           'With that said: for a generic Hamiltonian the four-flavor '
           'expansion is exact to round-off, like its two- and three-flavor '
           'siblings. A **stiff** spectrum is different, and 3+1 with an '
           'eV-scale sterile is stiff — the eigenvalues span four orders of '
           'magnitude with three of them clustered, and forming the three '
           'invariants in double precision compresses the $4\\times4$ matrix '
           'into three numbers and loses what separates the cluster. That is '
           'ordinary ill-conditioning of polynomial roots against their '
           'coefficients, so no better root-finder helps.\n\n'
           'The library refines the roots against the *matrix* instead — one '
           'Newton step on $\\chi(\\psi) = \\det(\\psi\\mathbb{1} - \\tilde '
           'H)$, which reads the Hamiltonian entries at full precision. '
           '`POLISH_ROOTS` controls it, and here is what it buys:'),
        code('def reference(Hs, L):\n'
             '    Ht = Hs - np.einsum("...ii->...", Hs).real[..., None, None]'
             '/4.0*np.eye(4)\n'
             '    w, V = np.linalg.eigh(Ht)\n'
             '    U = np.einsum("...ik,...k,...jk->...ij", V,\n'
             '                  np.exp(-1j*w*L), V.conj())\n'
             '    return np.abs(np.swapaxes(U, -1, -2))**2\n\n\n'
             'P_ref = reference(H_stack, baseline).reshape(-1, 16)\n'
             'for flag in (True, False):\n'
             '    oscprob4nu.POLISH_ROOTS = flag\n'
             '    P = oscprob4nu.probabilities_4nu(H_stack, baseline)\n'
             '    print("POLISH_ROOTS = %-5s : max |P - P_eigh| = %.2e"\n'
             '          % (flag, np.max(np.abs(P - P_ref))))\n'
             'oscprob4nu.POLISH_ROOTS = True'),
        md('Three alternatives were measured against `mpmath` at fifty '
           'decimal digits before settling on this one, and two of them '
           'lose in instructive ways.\n\n'
           '| Strategy for the roots | Relative error | Cost, 200k points | '
           'Closed form? |\n'
           '|---|---|---|---|\n'
           '| Closed form alone | 8.3e-11 | 0.17 s | yes |\n'
           '| **Closed form + one Newton step** | **1.1e-16** | 0.41 s | '
           'yes |\n'
           '| `numpy.linalg.eigvalsh` | 7.4e-16 | 0.17 s | no |\n'
           '| Closed form in `numpy.longdouble` | 4.5e-11 | 0.43 s | yes '
           '|\n\n'
           'The Newton step is about **seven times more accurate than '
           'LAPACK**, because `eigvalsh` reduces the matrix by similarity '
           'transforms that each carry a backward error $\\sim\\epsilon\\|H\\|$, '
           'while this converges onto the root of $\\det(\\psi\\mathbb{1} - '
           '\\tilde H)$ for the matrix it was handed.\n\n'
           'Extended precision buys under one digit rather than the three '
           'its extra mantissa suggests — the cluster amplifies coefficient '
           'error — and is slower, since `float128` is not '
           'hardware-vectorised. It is also silently `float64` on Apple '
           'Silicon and Windows, so it would appear to work while doing '
           'nothing.'),
        md('## A decoupled sterile state reproduces three flavors\n\n'
           'The strongest check available, and the one worth running if you '
           'ever doubt a four-flavor number: switch the three sterile angles '
           'off and the active block must be exactly what `oscprob3nu` gives '
           '— a different module, written against a different algebra.'),
        code('H_decoupled = '
             'hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(\n'
             '    gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, 0.0, 0.0, 0.0,\n'
             '    gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF, DM41)\n\n'
             'four = np.array(oscprob4nu.probabilities_4nu(\n'
             '    H_decoupled/energy, baseline)).reshape(4, 4)\n'
             'three = np.array(oscprob3nu.probabilities_3nu(\n'
             '    np.asarray(H_VAC_3NU)/energy, baseline)).reshape(3, 3)\n\n'
             'print("max |4nu active block - 3nu| = %.2e"\n'
             '      % np.max(np.abs(four[:3, :3] - three)))\n'
             'print("leakage into the sterile state = %.2e"\n'
             '      % np.max(np.abs(four[:3, 3])))'),
        md('## Why it stops at four\n\n'
           'Not for want of trying. The whole construction rests on writing '
           'the eigenvalues of the traceless Hamiltonian in closed form — '
           'that is, solving its characteristic polynomial **in radicals**. '
           'That polynomial has degree $n$.\n\n'
           'Quadratics, cubics and quartics are solvable in radicals. The '
           'general quintic is not: that is the Abel–Ruffini theorem, and '
           'Galois theory says why — $S_5$ is not a soluble group, while '
           '$S_2$, $S_3$ and $S_4$ are. At $n = 5$ there is simply no formula '
           'to write down, and the shortfall is a theorem rather than a gap '
           'in anyone\'s algebra.\n\n'
           'What does *not* stop is the philosophy. Nothing above the '
           'eigenvalues needs radicals — the interpolation over the roots, '
           'the fact that eigenvectors are never required, the probability '
           'construction — so feeding numerically computed eigenvalues into '
           'the same machinery degrades gracefully at any $n$. It would just '
           'no longer be a closed form, which is this library\'s reason to '
           'exist.\n\n'
           'So four flavors is both the natural stopping point and a useful '
           'one: it is exactly what 3+1 needs.'),
    ])

# ------------------------------------------------ 17 against other codes
books['17_against_other_codes.ipynb'] = notebook(
    'Against other codes',
    'Every check in the regression suite is, in the end, self-referential '
    'about **conventions**. `scipy.linalg.expm` confirms that $U$ is the '
    'matrix exponential of the Hamiltonian — but of the Hamiltonian this '
    'library built. The standard oscillation formulas were transcribed from '
    'the same papers the builders were. A mixing-matrix ordering, a sign, or '
    'a unit conversion that was wrong *consistently* across the package '
    'would pass all of it.\n\n'
    'This notebook compares against two things that would not: an '
    'independent external code, and a published closed form derived a '
    'different way.',
    [
        md('## Two comparators\n\n'
           '**[nuSQuIDS](https://github.com/arguelles/nuSQuIDS)** evolves the '
           'neutrino density matrix numerically, in C++. It shares no lineage '
           'with the closed-form SU($n$) expansions here, and it is '
           'configured from mixing *angles* and a density in g/cm³ rather '
           'than from a Hamiltonian — so agreement exercises '
           '`hamiltonians3nu` and `hamiltonians4nu` as much as the '
           'expansions.\n\n'
           'Its values are **frozen** into `tests/nusquids_reference.json` '
           'rather than computed here. nuSQuIDS ships manylinux wheels and '
           'would install in CI, but a comparison that runs in CI is a '
           'dependency on somebody else\'s release cadence. '
           '`tests/nusquids_reference.py` regenerates it deliberately.\n\n'
           '**Zaglauer–Schwarzer** (Z. Phys. C 40, 273 (1988)) give the exact '
           'eigenvalues of the matter Hamiltonian as roots of an explicit '
           'cubic in the vacuum parameters — never forming the matrix. That '
           'runs live, below, at any configuration and with no dependency.'),
        code('import json\n\n'
             'with open(os.path.join("..", "tests",\n'
             '                       "nusquids_reference.json")) as handle:\n'
             '    REF = json.load(handle)\n\n'
             'print("nuSQuIDS version :", REF["nusquids_version"])\n'
             'print("cases            :", len(REF["cases"]))\n'
             'for case in REF["cases"]:\n'
             '    print("   ", case["name"])'),
        md('## The convention dictionary\n\n'
           'Comparing two codes is mostly the work of matching what they '
           'mean, not what they compute. There are exactly two differences '
           'here, and **both are definitional** — neither is an error in '
           'either code.'),
        code('C = REF["constants"]\n\n'
             '# 1. The length unit\n'
             'print("CONV_KM_TO_INV_EV")\n'
             'print("  ours      : %.10e" % gd.CONV_KM_TO_INV_EV)\n'
             'print("  nuSQuIDS  : %.10e   (= 1 km / hbar c)" % C["km"])\n'
             'print("  relative  : %.2e"\n'
             '      % abs(gd.CONV_KM_TO_INV_EV/C["km"] - 1.0))\n\n'
             '# 2. The electron density behind V_CC\n'
             'rho, ye = gd.DENSITY_MATTER_CRUST_G_PER_CM3, '
             'REF["electron_fraction"]\n'
             'ne_theirs = rho*C["avogadro"]*ye/C["cm"]**3\n'
             'vcc_theirs = np.sqrt(2.0)*C["GF"]*ne_theirs\n'
             'print("\\nV_CC at rho = %.1f g/cm^3" % rho)\n'
             'print("  ours      : %.6e   (1 g / mean nucleon mass)"\n'
             '      % gd.VCC_EARTH_CRUST)\n'
             'print("  nuSQuIDS  : %.6e   (rho N_A Y_e)" % vcc_theirs)\n'
             'print("  ratio     : %.5f  -> %.2f%%, the nuclear mass defect"\n'
             '      % (vcc_theirs/gd.VCC_EARTH_CRUST,\n'
             '         100*(vcc_theirs/gd.VCC_EARTH_CRUST - 1.0)))'),
        md('The length unit looks harmless at $1.4\\times10^{-7}$, but it '
           'multiplies the **accumulated phase**. At 1300 km and 1 GeV that '
           'phase is some $2.5\\times10^3$ radian; with an eV-scale sterile '
           'it is $10^4$. So a seventh-decimal-place constant becomes a '
           'fourth-decimal-place probability.\n\n'
           '**One asymmetry to be honest about.** Both adjustments are made '
           'on *our* side, to match nuSQuIDS — so a sceptic could say the '
           'agreement below was fitted rather than found. Two things answer '
           'that. The length unit is not a free choice: $1\\,$km$/\\hbar c$ '
           'from CODATA is $5.0677307162\\times10^9$, so nuSQuIDS is right '
           'and our constant is simply rounded, which is checkable without '
           'reference to either code. And two constants cannot fit the '
           'result: what follows is eleven cases carrying between nine and '
           'sixteen probabilities each, agreeing to round-off. Two knobs do '
           'not buy that.'),
        md('## The comparison\n\n'
           'Every frozen case, computed both with our constants and with '
           'nuSQuIDS\' — so the effect of matching is visible rather than '
           'asserted.'),
        code('def ours(case, km, vcc, vnc):\n'
             '    """This library\'s probabilities for one frozen case.\n\n'
             '    Antineutrinos need *both* changes: conjugate the vacuum\n'
             '    Hamiltonian and reverse the sign of the potentials.\n'
             '    """\n'
             '    p = REF["parameters"]\n'
             '    s = np.sin\n'
             '    E = case["energy_gev"]*GEV\n'
             '    L = case["baseline_km"]*km\n'
             '    rho = case["density_g_cm3"]\n'
             '    anti = case.get("antineutrino", False)\n'
             '    dm31 = case.get("dm31", p["dm31"])\n'
             '    sign = -1.0 if anti else 1.0\n'
             '    if case["n_flavors"] == 3:\n'
             '        H0 = np.asarray(hamiltonians3nu.'
             'hamiltonian_3nu_vacuum_energy_independent(\n'
             '            s(p["th12"]), s(p["th23"]), s(p["th13"]),\n'
             '            p["dcp"], p["dm21"], dm31))\n'
             '        if anti:\n'
             '            H0 = np.conj(H0)\n'
             '        H = (H0/E if rho is None else\n'
             '             hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '                 H0, E, sign*vcc))\n'
             '        return np.array(oscprob3nu.probabilities_3nu(\n'
             '            H, L)).reshape(3, 3)\n'
             '    H0 = hamiltonians4nu.'
             'hamiltonian_4nu_vacuum_energy_independent(\n'
             '        s(p["th12"]), s(p["th23"]), s(p["th13"]),\n'
             '        s(p["th14"]), s(p["th24"]), s(p["th34"]),\n'
             '        p["dcp"], p["dm21"], dm31, p["dm41"])\n'
             '    if anti:\n'
             '        H0 = np.conj(H0)\n'
             '    H = (H0/E if rho is None else\n'
             '         hamiltonians4nu.hamiltonian_4nu_matter(\n'
             '             H0, E, sign*vcc, sign*vnc))\n'
             '    return np.array(oscprob4nu.probabilities_4nu(\n'
             '        H, L)).reshape(4, 4)\n\n\n'
             'def potentials(rho, matched):\n'
             '    """(V_CC, V_NC) in either convention."""\n'
             '    if rho is None:\n'
             '        return None, None\n'
             '    if not matched:\n'
             '        scale = rho/gd.DENSITY_MATTER_CRUST_G_PER_CM3\n'
             '        return (gd.VCC_EARTH_CRUST*scale,\n'
             '                gd.VNC_EARTH_CRUST*scale)\n'
             '    ye = REF["electron_fraction"]\n'
             '    ne = rho*C["avogadro"]*ye/C["cm"]**3\n'
             '    nn = rho*C["avogadro"]*(1.0 - ye)/C["cm"]**3\n'
             '    return (np.sqrt(2.0)*C["GF"]*ne,\n'
             '            -C["GF"]*nn/np.sqrt(2.0))\n\n\n'
             'names, unmatched, matched_ = [], [], []\n'
             'for case in REF["cases"]:\n'
             '    theirs = np.array(case["probabilities"])\n'
             '    rho = case["density_g_cm3"]\n'
             '    a = np.max(np.abs(ours(case, gd.CONV_KM_TO_INV_EV,\n'
             '                           *potentials(rho, False)) - theirs))\n'
             '    b = np.max(np.abs(ours(case, C["km"],\n'
             '                           *potentials(rho, True)) - theirs))\n'
             '    names.append(case["name"])\n'
             '    unmatched.append(a)\n'
             '    matched_.append(b)\n\n'
             'print("%-36s %11s %11s" % ("case", "our conv.", "matched"))\n'
             'for n, a, b in zip(names, unmatched, matched_):\n'
             '    print("%-36s %11.2e %11.2e" % (n, a, b))'),
        code('fig, ax = plt.subplots(figsize=(7.6, 4.4))\n'
             'y = np.arange(len(names))\n'
             'ax.barh(y - 0.2, unmatched, height=0.38,\n'
             '        label="our constants")\n'
             'ax.barh(y + 0.2, matched_, height=0.38,\n'
             '        label="conventions matched")\n'
             'ax.set_yticks(y)\n'
             'ax.set_yticklabels([n.replace(", ", "\\n", 1) for n in names],\n'
             '                   fontsize=8)\n'
             'ax.set_xscale("log")\n'
             'ax.set_xlim(1e-16, 1e-3)\n'
             'ax.set_xlabel(r"max $|P_{\\rm NuOscProbExact} - '
             'P_{\\rm nuSQuIDS}|$")\n'
             'ax.set_title("Agreement with an independent code")\n'
             'ax.invert_yaxis()\n'
             'ax.legend(loc="lower right")\n'
             'ax.grid(True, axis="x", alpha=0.3)\n'
             'fig.tight_layout()\n'
             'plt.show()'),
        md('Matched, the two codes agree between $10^{-10}$ and '
           '$10^{-15}$ — round-off, for two implementations that share '
           'nothing. Everything to the left of that in the "our constants" '
           'bars is convention, not physics.\n\n'
           'The **four-flavor rows** deserve their own note. `oscprob4nu` is '
           'new, and until this comparison it was validated only against '
           '`scipy.linalg.expm` and against its own three-flavor limit. It '
           'is now checked by an outside code, across a range that tracks '
           'how stiff the spectrum is: $3.9\\times10^{-16}$ at a short '
           'baseline where it is benign, $2.4\\times10^{-11}$ in matter, and '
           '$3.0\\times10^{-10}$ for antineutrinos in matter, which is the '
           'stiffest case here. That last one is the accuracy '
           '`POLISH_ROOTS` documents, arrived at independently; it is '
           'attributed at the end of this notebook rather than left '
           'hanging.'),
        md('## Zaglauer–Schwarzer, live\n\n'
           'The frozen data covers eleven configurations. This covers '
           '**all** of them, with no dependency.\n\n'
           'Zaglauer and Schwarzer give the eigenvalues of\n\n'
           '$$M^2_{\\rm eff} = U M^2 U^\\dagger + '
           '\\mathrm{diag}(A, 0, 0), \\qquad A = 2EV_{CC}$$\n\n'
           'as roots of an explicit cubic in the vacuum parameters and $A$ — '
           'the matrix is never formed. Comparing those against the spectrum '
           'of the Hamiltonian we *do* form tests the construction: the PMNS '
           'parametrization, where the potential is added, and the relative '
           'normalization of the two terms.\n\n'
           'The $|U_{ei}|^2$ below come straight from the angles, not from '
           '`pmns_mixing_matrix`, so a mistake in our mixing matrix cannot '
           'cancel against itself.'),
        code('def zaglauer_schwarzer(s12, s13, dm21, dm31, A):\n'
             '    """Exact eigenvalues of M^2_eff, from the vacuum '
             'parameters."""\n'
             '    c13sq = 1.0 - s13*s13\n'
             '    ue1, ue2, ue3 = (1.0 - s12*s12)*c13sq, s12*s12*c13sq, '
             's13*s13\n'
             '    x = dm21 + dm31 + A\n'
             '    y = dm21*dm31 + A*(dm21*(1.0 - ue2) + dm31*(1.0 - ue3))\n'
             '    z = A*dm21*dm31*ue1\n'
             '    p = y - x*x/3.0\n'
             '    q = -2.0*x**3/27.0 + x*y/3.0 - z\n'
             '    t = 2.0*np.sqrt(-p/3.0)\n'
             '    th = np.arccos(np.clip(3.0*q/(p*t), -1.0, 1.0))\n'
             '    return np.sort(t*np.cos((th + 2.0*np.pi*np.arange(3))/3.0)\n'
             '                   + x/3.0)\n\n\n'
             'VCC = gd.VCC_EARTH_CRUST\n'
             'E_gev = np.logspace(-0.5, 2.0, 400)\n\n'
             'tracks, deviations = [], []\n'
             'for e in E_gev:\n'
             '    E = e*GEV\n'
             '    H = np.asarray(hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '        H_VAC_3NU, E, VCC))\n'
             '    mine = np.sort(np.linalg.eigvalsh(2.0*E*H))\n'
             '    theirs = zaglauer_schwarzer(gd.S12_NO_BF, gd.S13_NO_BF,\n'
             '                                gd.D21_NO_BF, gd.D31_NO_BF,\n'
             '                                2.0*E*VCC)\n'
             '    tracks.append(mine)\n'
             '    deviations.append(np.max(np.abs(mine - theirs))\n'
             '                      / np.max(np.abs(mine)))\n'
             'tracks = np.array(tracks)\n'
             'deviations = np.array(deviations)\n\n'
             'print("worst relative deviation over the scan: %.2e"\n'
             '      % deviations.max())'),
        code('fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.2),\n'
             '                         sharex=True,\n'
             '                         gridspec_kw={"height_ratios": '
             '[2, 1]})\n\n'
             'for i in range(3):\n'
             '    axes[0].loglog(E_gev, np.abs(tracks[:, i]),\n'
             '                   label=r"$|\\lambda_%d|$" % (i + 1))\n'
             'c2th13 = 1.0 - 2.0*gd.S13_NO_BF**2\n'
             'E_res = gd.D31_NO_BF*c2th13/(2.0*VCC)/GEV\n'
             'axes[0].axvline(E_res, color="k", ls=":", lw=1)\n'
             'axes[0].text(E_res*1.1, 3e-3, "MSW resonance", fontsize=9)\n'
             'axes[0].set_ylabel(r"$|\\lambda_i|$  [eV$^2$]")\n'
             'axes[0].set_title("Matter eigenvalues, and their agreement "\n'
             '                  "with Zaglauer-Schwarzer")\n'
             'axes[0].legend(ncol=3)\n\n'
             'axes[1].loglog(E_gev, np.maximum(deviations, 1e-18), "-")\n'
             'axes[1].axhline(1e-15, color="k", ls="--", lw=1)\n'
             'axes[1].text(0.4, 1.6e-15, "round-off", fontsize=9)\n'
             'axes[1].axvline(E_res, color="k", ls=":", lw=1)\n'
             'axes[1].set_ylim(1e-18, 1e-10)\n'
             'axes[1].set_xlim(E_gev[0], E_gev[-1])\n'
             'axes[1].set_xlabel("Energy [GeV]")\n'
             'axes[1].set_ylabel("relative deviation")\n'
             'fig.tight_layout()\n'
             'plt.show()'),
        md('The agreement holds through the resonance, which is where the '
           'two upper eigenvalues come closest and any error in the '
           'construction would be most visible.\n\n'
           'One transcription note, because it is the kind of mistake that '
           'usually produces a plausible wrong curve rather than an obvious '
           'one: the depressed cubic\'s constant term is '
           '$-2x^3/27 + xy/3 - z$. Written first with the signs of the first '
           'two terms flipped, it disagreed with the spectrum by 50% — large '
           'enough to be unmissable, which was luck. The corrected form is '
           'identical to the resolvent cubic already inside `oscprob4nu`, '
           'which is where it should have been checked against from the '
           'start.'),
        md('## Which one should you use?\n\n'
           'The comparison above is about agreement, but the question a '
           'reader actually has is which tool to reach for. The two are not '
           'competitors so much as different scopes that happen to overlap.\n\n'
           '| | NuOscProbExact | nuSQuIDS |\n'
           '|---|---|---|\n'
           '| Constant or piecewise-constant density | ✔ closed form | ✔ |\n'
           '| Arbitrary Hermitian Hamiltonian | ✔ | ✔ |\n'
           '| The evolution operator, returned directly | ✔ | works in the '
           'density-matrix picture |\n'
           '| Two, three or four flavors | ✔ | ✔ |\n'
           '| More than four flavors | — (see notebook 16) | ✔ |\n'
           '| Smoothly varying profiles | approximate by slabs | ✔ '
           'integrated |\n'
           '| Neutrino–nucleon interactions, attenuation | — | ✔ |\n'
           '| Tau regeneration | — | ✔ |\n'
           '| Non-unitary / open-system evolution | — | ✔ |\n\n'
           'On the overlap, the closed form is **twenty to sixty times '
           'quicker** on a large energy scan — sixty with the optional Numba '
           'extra, twenty on the default NumPy-only install — because it '
           'evaluates an expression where nuSQuIDS integrates a density '
           'matrix.\n\n'
           'That is not a defect in nuSQuIDS, and it is not really a '
           'benchmark result — it is what the two are *for*. The integration '
           'buys generality this library does not have: interactions, '
           'attenuation, tau regeneration, open-system evolution, more than '
           'four flavors, and profiles that vary continuously. If a '
           'calculation needs any of those, the cost of the general machinery '
           'is simply what it costs, and no closed form is going to help. If '
           'it needs none of them, the closed form is exact, has no tolerance '
           'to set, and is much faster.\n\n'
           'The same logic points the other way at the boundary of *this* '
           'library: for a profile that varies appreciably over an '
           'oscillation length — the Sun, adiabatic MSW — slabbing stops '
           'being practical and a Magnus-type method is the right tool. '
           'Notebook 14 shows exactly where that wall is.'),
        md('## What this covers, and what it does not\n\n'
           'Worth being precise, since "validated against an independent '
           'code" is easy to over-claim.\n\n'
           '| Checked | By |\n'
           '|---|---|\n'
           '| $U$ is $e^{-iHL}$ | `scipy.linalg.expm`, in the suite |\n'
           '| The Hamiltonian we build is the right one | nuSQuIDS, and '
           'Zaglauer–Schwarzer for the matter spectrum |\n'
           '| Mixing-matrix conventions, $\\theta_{23}$, $\\delta_{CP}$ | '
           'nuSQuIDS, and the vacuum reference formulas |\n'
           '| **The antineutrino rule** — conjugate *and* flip | nuSQuIDS, '
           'three and four flavors, vacuum and matter |\n'
           '| **The mass ordering sign** | nuSQuIDS, inverted-ordering case |\n'
           '| Four-flavor probabilities | nuSQuIDS, at 3e-10 to 4e-16 |\n'
           '| Arbitrary matter configurations | Zaglauer–Schwarzer, live |\n\n'
           'The antineutrino and ordering rows were added after the fact, and '
           'the omission is worth admitting: the first version of this '
           'comparison covered only neutrinos in the normal ordering. Those '
           'are precisely the two sign conventions this library has a history '
           'of getting wrong — the 1.1.0 release fixed exactly such a bug — '
           'so a comparison meant to catch convention errors was, at first, '
           'blind to the ones worth catching. The three-flavor ones now '
           'agree at $10^{-13}$ to $10^{-15}$; the four-flavor '
           'antineutrino case is the stiffest configuration here and is '
           'discussed below.\n\n'
           'What is still **not** covered: Zaglauer–Schwarzer eigenvalues are '
           'independent of $\\theta_{23}$ and $\\delta_{CP}$, since only the '
           'electron row of $U$ enters the matter term — and eigenvalues are '
           'not probabilities. The frozen nuSQuIDS set is eleven '
           'configurations, not a parameter space: $\\delta_{CP}$, '
           '$\\theta_{23}$ and the sterile angles are each held at one value. '
           'And no external code here exercises `slabs` or `earth`, whose '
           'layered composition is checked internally against `expm` '
           'instead.\n\n'
           'One residual is worth naming rather than hiding. The four-flavor '
           '*antineutrino* case agrees least well, at $3\\times10^{-10}$. '
           'That is not a convention gap: against `scipy.linalg.expm` on the '
           'same Hamiltonian, our error is $3.0\\times10^{-10}$ and '
           'nuSQuIDS\' is $3.4\\times10^{-12}$ — so the whole residual is '
           'ours, and it is the stiff-spectrum limit that `POLISH_ROOTS` '
           'documents, reached independently. A test in the suite pins that '
           'attribution, so a convention error could not later hide behind '
           'the same loose tolerance.'),
    ])

# Figures lifted out of the executed notebooks for the README gallery and the
# documentation's recipes page.  Extracting them rather than drawing them
# again is what keeps the three in step: there is one piece of code behind
# each picture, and it is the one a reader can go and run.
#
#     (notebook, index of the figure within it) -> file name
GALLERY = {
    ('02_vacuum_oscillations.ipynb', 1): 'gallery_vacuum.png',
    ('03_matter_nsi_liv.ipynb', 0): 'gallery_matter.png',
    ('04_oscillogram.ipynb', 0): 'gallery_oscillogram.png',
    ('05_biprobability.ipynb', 1): 'gallery_biprobability.png',
    ('06_earth_and_prem.ipynb', 0): 'gallery_prem.png',
    ('07_earth_probabilities.ipynb', 1): 'gallery_earth.png',
    ('08_unusual_density_profiles.ipynb', 1): 'gallery_profiles.png',
    ('12_ordering_and_octant.ipynb', 2): 'gallery_ordering.png',
    ('16_four_neutrinos.ipynb', 0): 'gallery_sterile.png',
    ('16_four_neutrinos.ipynb', 1): 'gallery_sterile_earth.png',
}

GALLERY_DIR = pathlib.Path('img') / 'gallery'


def extract_gallery():
    """Writes the gallery figures out of the executed notebooks."""
    import base64

    GALLERY_DIR.mkdir(parents=True, exist_ok=True)
    written = 0
    for (notebook, index), filename in sorted(GALLERY.items()):
        nb = nbf.read(OUT / notebook, as_version=4)
        images = [output['data']['image/png']
                  for cell in nb.cells
                  for output in cell.get('outputs', [])
                  if 'image/png' in output.get('data', {})]
        if index >= len(images):
            raise SystemExit('%s has no figure %d' % (notebook, index))
        (GALLERY_DIR / filename).write_bytes(base64.b64decode(images[index]))
        written += 1
    print('  wrote %d gallery figures to %s' % (written, GALLERY_DIR))


def build():
    """Writes every notebook, executes it, and checks it kept its outputs."""
    from nbclient import NotebookClient
    from nbclient.exceptions import CellExecutionError

    OUT.mkdir(exist_ok=True)
    for name, nb in books.items():
        nbf.write(nb, OUT / name)

    failed = []
    for path in sorted(OUT.glob('*.ipynb')):
        nb = nbf.read(path, as_version=4)
        started = time.perf_counter()
        try:
            NotebookClient(
                nb, timeout=1800, kernel_name='python3',
                resources={'metadata': {'path': str(path.parent)}}).execute()
            nbf.write(nb, path)
            print('  executed %-40s %5.1f s'
                  % (path.name, time.perf_counter()-started))
        except CellExecutionError as error:
            failed.append(path.name)
            print('  FAILED   %s\n%s' % (path.name, error))

    if failed:
        raise SystemExit('notebooks failed to execute: %s'
                         % ', '.join(failed))

    # The same check CI applies, run here so a stripped notebook is caught
    # before it is committed rather than after it is pushed.
    bare = [path.name for path in sorted(OUT.glob('*.ipynb'))
            if not any(cell.get('outputs') for cell
                       in nbf.read(path, as_version=4).cells
                       if cell.cell_type == 'code')]
    if bare:
        raise SystemExit('notebooks carry no stored outputs: %s'
                         % ', '.join(bare))

    print('  all %d notebooks executed and carry stored outputs'
          % len(books))

    extract_gallery()


if __name__ == '__main__':
    build()
