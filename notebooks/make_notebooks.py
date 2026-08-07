"""Generate and execute the NuOscProbExact example notebooks.

Run from the repository root::

    python notebooks/make_notebooks.py

The notebooks are generated from this one script rather than edited by hand,
so that all twenty notebooks share a setup cell, a plotting style and a voice,
and so that regenerating them after an API change is a single command.

**Edit this file, not the notebooks.**  Anything written directly into a
``.ipynb`` here is lost the next time this runs.

Generation and execution are deliberately one step.  Generating rewrites every
notebook, discarding its stored outputs; executing puts them back.  Doing only
the first leaves twenty notebooks that still open and still run but show a
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
        md('## What the Hamiltonian is\n\n'
           'Every routine here takes a Hamiltonian — a Hermitian matrix — '
           'and a baseline. If you have not met the oscillation Hamiltonian '
           'before, this is the one piece of physics worth having up front, '
           'because everything else in these notebooks is a variation on '
           'it.\n\n'
           'Neutrinos are produced and detected as **flavor** states '
           '($\\nu_e, \\nu_\\mu, \\nu_\\tau$) but propagate as **mass** '
           'states ($\\nu_1, \\nu_2, \\nu_3$). The two bases are related by '
           'the PMNS matrix $U$ — three angles $\\theta_{12}, \\theta_{13}, '
           '\\theta_{23}$ and a CP-violating phase $\\delta_{CP}$ — and the '
           'Hamiltonian is diagonal in the mass basis, so in the flavor '
           'basis it reads\n\n'
           '$$ H_{\\rm vac} = \\frac{1}{2E}\\, U M^2 U^\\dagger , \\qquad '
           'M^2 = \\mathrm{diag}\\left(0,\\, \\Delta m^2_{21},\\, '
           '\\Delta m^2_{31}\\right) . $$\n\n'
           'Only mass-squared *differences* appear, which is why the first '
           'entry can be set to zero. Matter adds a potential that the '
           'electron flavor feels and the others do not,\n\n'
           '$$ H = H_{\\rm vac} + \\mathrm{diag}(V_{CC},\\, 0,\\, 0) , $$\n\n'
           'and that is the whole content. Non-standard interactions, '
           'Lorentz-invariance violation and a sterile state are each a '
           'different Hermitian matrix in that same slot, evaluated by the '
           'same routine.\n\n'
           '`H_VAC_3NU` in the setup cell is exactly $U M^2 U^\\dagger/2$ at '
           'the NuFit 4.0 best fit for the normal ordering. Building it by '
           'hand shows there is nothing hidden in the builder:'),
        code('# The mixing angles, as sin^2(theta) -- which is how they are\n'
             '# quoted -- and the phase in degrees.\n'
             'print("sin^2(th12) = %.3f" % gd.S12_NO_BF**2)\n'
             'print("sin^2(th23) = %.3f" % gd.S23_NO_BF**2)\n'
             'print("sin^2(th13) = %.5f" % gd.S13_NO_BF**2)\n'
             'print("delta_CP    = %.0f degrees" '
             '% np.degrees(gd.DCP_NO_BF))\n'
             'print("dm21^2      = %.3e eV^2" % gd.D21_NO_BF)\n'
             'print("dm31^2      = %.3e eV^2" % gd.D31_NO_BF)\n\n'
             '# Note: the builders take sin(theta), not theta itself.\n'
             'U = np.asarray(hamiltonians3nu.pmns_mixing_matrix(\n'
             '    gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, gd.DCP_NO_BF))\n'
             'M2 = np.diag([0.0, gd.D21_NO_BF, gd.D31_NO_BF]).astype(complex)\n\n'
             'by_hand = U @ M2 @ U.conj().T / 2.0\n'
             'reference = np.asarray(H_VAC_3NU)\n'
             'print("\\nU M^2 U^dag / 2  vs  H_VAC_3NU")\n'
             'print("  absolute : %.2e eV^2"\n'
             '      % np.max(np.abs(by_hand - reference)))\n'
             'print("  relative : %.2e"\n'
             '      % (np.max(np.abs(by_hand - reference))\n'
             '         / np.max(np.abs(reference))))'),
        md('## One three-flavor probability\n\n'
           'A 1 GeV neutrino travelling 1300 km in vacuum — roughly the DUNE '
           'baseline. The nine probabilities come back with the *initial* '
           'flavor varying slowest.'),
        code('energy = 1.0*GEV\n'
             'baseline = 1300.0*KM\n\n'
             '# The vacuum Hamiltonian at this energy\n'
             'H = np.asarray(H_VAC_3NU)/energy\n\n'
             'prob = oscprob3nu.probabilities_3nu(H, baseline)\n\n'
             'flavors = ["e", "mu", "tau"]\n'
             'for k, value in enumerate(prob):\n'
             '    a, b = divmod(k, 3)\n'
             '    print("P_%-7s = %.6f"\n'
             '          % (flavors[a]+flavors[b], value))'),
        md('### Reading the nine numbers\n\n'
           'The convention, used everywhere in this library and in every '
           'notebook after this one:\n\n'
           '$$ P_{\\alpha\\beta} \\equiv P(\\nu_\\alpha \\to \\nu_\\beta) '
           '= \\left|[U(L)]_{\\beta\\alpha}\\right|^2 $$\n\n'
           'so the **initial** flavor varies slowest. Later notebooks index '
           'the result directly — `p[:, 3]` for $P_{\\mu e}$, `p[:, 4]` for '
           '$P_{\\mu\\mu}$ — and this is the map they are using:'),
        code('print("index   from      to        symbol")\n'
             'for k in range(9):\n'
             '    a, b = divmod(k, 3)\n'
             '    print("  %d     nu_%-4s -> nu_%-5s  P_%s%s"\n'
             '          % (k, flavors[a], flavors[b],\n'
             '             flavors[a], flavors[b]))'),
        md('Each initial flavor conserves probability, which is the first '
           'thing to check of any oscillation code:'),
        code('prob = np.array(prob)\n'
             'for start, flavor in zip((0, 3, 6), ("e", "mu", "tau")):\n'
             '    print("sum over final flavors, from nu_%-4s = %.15f"\n'
             '          % (flavor, prob[start:start+3].sum()))'),
        md('## The trace does not matter\n\n'
           'Adding a multiple of the identity to $H$ shifts every eigenvalue '
           'by the same amount, so it contributes an overall phase that '
           'cancels in $|U_{\\beta\\alpha}|^2$. The library drops the trace '
           'internally and works with the traceless part throughout.\n\n'
           'This is worth knowing because it means an overall energy offset '
           'in your Hamiltonian — a flavor-universal potential, for instance '
           '— can simply be left out.'),
        code('shift = 3.0*np.max(np.abs(H))     # comparable to H itself\n\n'
             'p_plain = np.array(oscprob3nu.probabilities_3nu(H, baseline))\n'
             'p_shift = np.array(oscprob3nu.probabilities_3nu(\n'
             '    H + shift*np.eye(3), baseline))\n\n'
             'print("max |P(H) - P(H + c*1)| = %.2e"\n'
             '      % np.max(np.abs(p_plain - p_shift)))'),
        md('One caveat, since someone will try it: that holds *arithmetically* '
           'for any shift, but not numerically for an absurd one. The '
           'traceless part is recovered by subtracting numbers of order the '
           'shift, so adding $10^{10}\\,|H|$ leaves only a few significant '
           'digits of the physics and the agreement degrades to about '
           '$10^{-7}$. That is ordinary floating-point cancellation rather '
           'than anything about the method — but if your Hamiltonian carries '
           'a huge flavor-universal term, subtract it yourself before '
           'passing it in.'),
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
           'that matrix.\n\n'
           'The matrix below is not any physical scenario — it is a reminder '
           'that the routine neither knows nor needs to know where the '
           'numbers came from.'),
        code('H_arbitrary = np.array([[0.0, 1.0+2.0j, 0.5],\n'
             '                        [1.0-2.0j, 1.0, 0.0],\n'
             '                        [0.5, 0.0, -1.0]], dtype=complex)\n\n'
             'p_arb = oscprob3nu.probabilities_3nu(H_arbitrary, 1.0)\n\n'
             'for k, value in enumerate(p_arb):\n'
             '    a, b = divmod(k, 3)\n'
             '    print("P_%s%s = %.6f" % (flavors[a], flavors[b], value))'),
        md('## Pass arrays, do not loop\n\n'
           'The single most useful thing to know about performance here: the '
           'routines broadcast. A stack of Hamiltonians, an array of '
           'baselines, or both, are evaluated in one pass — tens of times '
           'faster than the equivalent Python loop, and the results agree to '
           'round-off. Notebook 09 measures both the speedup and the '
           'agreement.'),
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
    'two figures of [arXiv:1904.12391](https://arxiv.org/abs/1904.12391), '
    'produced inline and reproducible from the repository alone.',
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
           'A non-standard interaction (NSI) is a hypothetical extra '
           'neutrino–matter coupling. It enters as a potential with its own '
           'flavor structure, written relative to $V_{CC}$ as a matrix of '
           'dimensionless strengths $\\epsilon_{\\alpha\\beta}$:\n\n'
           '$$ H = H_{\\rm vac} + V_{CC}\\begin{pmatrix}'
           '1 + \\epsilon_{ee} & \\epsilon_{e\\mu} & \\epsilon_{e\\tau} \\\\'
           '\\epsilon_{e\\mu}^* & \\epsilon_{\\mu\\mu} & \\epsilon_{\\mu\\tau} \\\\'
           '\\epsilon_{e\\tau}^* & \\epsilon_{\\mu\\tau}^* & '
           '\\epsilon_{\\tau\\tau}\\end{pmatrix} $$\n\n'
           'The `1` recovers standard matter when every $\\epsilon$ '
           'vanishes. `globaldefs` carries an illustrative set — these are '
           'the values the curve below actually uses, not a fit to '
           'anything:'),
        code('names = ["eps_ee", "eps_emu", "eps_etau",\n'
             '         "eps_mumu", "eps_mutau", "eps_tautau"]\n'
             'for n, v in zip(names, gd.EPS_3):\n'
             '    print("  %-11s = %s" % (n, v))'),
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
           'A CPT-odd Lorentz-invariance-violating (LIV) background couples '
           'to neutrinos through a term that is *added* to the Hamiltonian '
           'with its own mixing structure. The physically important point is '
           'the energy dependence: the vacuum term falls as $1/E$ and this '
           'one does not, so it is invisible at low energy and dominant at '
           'high energy — which is why the curve below only separates at the '
           'right-hand end.\n\n'
           'The parameters below are illustrative defaults from '
           '`globaldefs`, shown because a plot of "LIV" means nothing '
           'without them:'),
        code('print("mixing, as sin(xi):")\n'
             'for n, v in (("sin(xi12)", gd.SXI12), ("sin(xi23)", gd.SXI23),\n'
             '             ("sin(xi13)", gd.SXI13)):\n'
             '    print("  %-10s = %.3f" % (n, v))\n'
             'print("phase  : %.3f rad" % gd.DXICP)\n'
             'print("eigenvalues b1, b2, b3 : %.1e, %.1e, %.1e"\n'
             '      % (gd.B1, gd.B2, gd.B3))\n'
             'print("scale  Lambda          : %.1e eV" % gd.LAMBDA)'),
        md('Two things in that list are worth reading rather than skipping '
           'past, and they are the reason it is printed at all.\n\n'
           'The LIV mixing angles are all **zero**, so the rotation is the '
           'identity and the added term is simply\n\n'
           '$$ \\frac{E}{\\Lambda}\\,\\mathrm{diag}(b_1, b_2, b_3) $$\n\n'
           'in the flavor basis — a background that shifts the flavor '
           'eigenvalues without introducing any mixing of its own. That is '
           'the simplest non-trivial case, and it is enough to move the '
           'curve. And since $b_1 = b_2$ here, only the difference '
           '$b_3 - b_1$ survives dropping the trace; a common $b$ would do '
           'nothing at all.\n\n'
           'The energy dependence is the signature. The vacuum term falls '
           'as $1/E$ while this one *rises* as $E$, so there is a crossover '
           'where they are equal, at\n\n'
           '$$ E_\\times = \\sqrt{\\frac{\\Delta m^2_{31}\\Lambda}{2b_3}} '
           '\\approx 0.8~\\mathrm{GeV} $$\n\n'
           'for these values. Below it the vacuum term wins and the curve is '
           'the ordinary one; above it the LIV term takes over, and by '
           '100 GeV it is four orders of magnitude larger. That growth with '
           'energy is what makes high-energy neutrinos the place to look for '
           'this kind of physics.'),
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
           'this is tens of times slower; **notebook 09** measures exactly '
           'how much, on whatever machine runs it.'),
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
        md('## What sets the size of the effect\n\n'
           'The ellipse is a picture; the number behind it is the '
           '**Jarlskog invariant**, the one parametrisation-independent '
           'measure of CP violation in the mixing matrix:\n\n'
           '$$ J_{CP} = c_{12}s_{12}\\,c_{23}s_{23}\\,c_{13}^2 s_{13}\\,'
           '\\sin\\delta_{CP} . $$\n\n'
           'It vanishes if any mixing angle is zero or if '
           '$\\sin\\delta_{CP} = 0$ — and when it does, the vacuum ellipse '
           'collapses onto the diagonal and there is nothing to measure. '
           'Every CP-odd term in the vacuum probabilities is proportional '
           'to it.\n\n'
           '`hamiltonians3nu.J` returns the quartic product '
           '$U^*_{\\alpha k}U_{\\beta k}U_{\\alpha j}U^*_{\\beta j}$, whose '
           'imaginary part is $J_{CP}$ up to a sign convention:'),
        code('def jarlskog(dcp):\n'
             '    # J_CP from the library quartic product\n'
             '    U = hamiltonians3nu.pmns_mixing_matrix(\n'
             '        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, dcp)\n'
             '    return -np.imag(hamiltonians3nu.J(U, 0, 1, 0, 1))\n\n\n'
             'def jarlskog_closed_form(dcp):\n'
             '    # The textbook expression, for comparison\n'
             '    s12, s23, s13 = gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF\n'
             '    c12, c23, c13 = [np.sqrt(1.0-s*s)\n'
             '                     for s in (s12, s23, s13)]\n'
             '    return c12*s12*c23*s23*c13*c13*s13*np.sin(dcp)\n\n\n'
             'for d in (0.0, np.pi/2, gd.DCP_NO_BF, np.pi):\n'
             '    print("dCP = %6.1f deg :  J = %+.5e   "\n'
             '          "(closed form %+.5e)"\n'
             '          % (np.degrees(d), jarlskog(d),\n'
             '             jarlskog_closed_form(d)))'),
        code('dcp_scan = np.linspace(0.0, 2.0*np.pi, 400)\n'
             'J_scan = np.array([jarlskog(d) for d in dcp_scan])\n\n'
             'fig, ax = plt.subplots(figsize=(6.4, 3.4))\n'
             'ax.plot(np.degrees(dcp_scan), J_scan)\n'
             'ax.axhline(0.0, color="0.7", lw=0.8)\n'
             'ax.axvline(np.degrees(gd.DCP_NO_BF), color="k", ls=":",\n'
             '           lw=1)\n'
             'ax.set_xlim(0.0, 360.0)\n'
             'ax.set_xlabel(r"$\\delta_{CP}$ [degrees]")\n'
             'ax.set_ylabel(r"$J_{CP}$")\n'
             'ax.set_title("The Jarlskog invariant "\n'
             '             "(dotted: the NuFit best fit)")\n'
             'plt.show()\n\n'
             'print("largest |J| over the scan : %.5e"\n'
             '      % np.max(np.abs(J_scan)))'),
        md('The zeros at $\\delta_{CP} = 0$ and $180^\\circ$ are exactly '
           'where the vacuum ellipse degenerates to a point on the '
           'diagonal. Note the scale: $J_{CP}$ peaks near $0.033$, small '
           'because it carries a factor of $s_{13}$ — the CP-violating '
           'effect is suppressed by the smallest mixing angle, which is a '
           'large part of why the measurement is hard.'),
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
        md('## From slabs to a probability\n\n'
           'The pieces are now in place, and this is where `slabs` earns its '
           'name. `earth` gives the widths and densities; '
           '`hamiltonians3nu.hamiltonian_3nu_matter` turns each density into '
           'a Hamiltonian; and `slabs.probabilities_3nu_slabs` solves each '
           'piece exactly and multiplies the evolution operators together, '
           'in order.\n\n'
           'Nothing here is approximated except the profile itself: the '
           'evolution *within* each slab is the same closed form as anywhere '
           'else in this library. `earth.probabilities_3nu_earth` is exactly '
           'this three-step recipe wrapped up, and doing it by hand once is '
           'the way to see what it does — and the way to build a profile of '
           'your own, which is what notebook 08 goes on to do.'),
        code('widths_km, densities = earth.earth_slabs(\n'
             '    -0.8, n_slabs_per_segment=6)\n\n'
             '# One Hamiltonian per slab, from that slab\'s density\n'
             'H_slabs = hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '    H_VAC_3NU, 10.0*GEV,\n'
             '    earth.matter_potential(densities))\n\n'
             '# Solve each exactly, multiply the operators in order\n'
             'p_by_hand = slabs.probabilities_3nu_slabs(\n'
             '    H_slabs, widths_km*KM)\n\n'
             '# The wrapper that does all three steps for you\n'
             'p_wrapped = earth.probabilities_3nu_earth(\n'
             '    H_VAC_3NU, 10.0*GEV, -0.8, n_slabs_per_segment=6)\n\n'
             'print("%d slabs, %.0f km total" % (len(widths_km), '
             'widths_km.sum()))\n'
             'print("P_mumu, by hand   = %.9f" % p_by_hand[4])\n'
             'print("P_mumu, wrapped   = %.9f" % p_wrapped[4])\n'
             'print("difference        = %.2e"\n'
             '      % abs(p_by_hand[4] - p_wrapped[4]))'),
        md('Identical, as they must be — the wrapper is the recipe, not a '
           'different calculation.\n\n'
           'If you need the evolution **operator** across the trajectory '
           'rather than the probabilities — to compose it with something '
           'else, or to propagate a density matrix — '
           '`slabs.evolution_operator_3nu_slabs` returns it from the same '
           'inputs.'),
        code('U = np.array(slabs.evolution_operator_3nu_slabs(\n'
             '    H_slabs, widths_km*KM))\n\n'
             'print("U has shape", U.shape)\n'
             'print("unitary to      %.2e"\n'
             '      % np.max(np.abs(U.conj().T @ U - np.eye(3))))\n'
             '# P_(alpha->beta) = |U_(beta,alpha)|^2\n'
             'print("|U[1,1]|^2      = %.9f  (= P_mumu above)"\n'
             '      % abs(U[1, 1])**2)'),
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
        md('## Asking for an accuracy instead\n\n'
           'The plot above answers *how does the error fall with `n`?* — '
           'but not the question a reader actually has, which is *what '
           '`n` do I need?* Those are different questions, because **the '
           'error at a fixed `n` is not a fixed number.** It depends '
           'strongly on the energy, and on the zenith angle through the '
           'chord length and which shells are crossed.\n\n'
           'So `rtol` and `atol` invert the relationship: give a '
           'tolerance instead of a subdivision and the chord is refined '
           'until the measured error meets it.'),
        code('atol = 1.0e-5\n'
             'energies_gev = np.array([1.0, 3.0, 10.0, 20.0, 40.0, 80.0])\n\n'
             'print("cos(theta_z) = -0.8, atol = %.0e\\n" % atol)\n'
             'print("%9s %14s %11s" % ("E [GeV]", "error at n=8", '
             '"n needed"))\n'
             'for e in energies_gev:\n'
             '    converged = np.array(earth.probabilities_3nu_earth(\n'
             '        H_VAC_3NU, e*GEV, -0.8, n_slabs_per_segment=512))\n'
             '    default = np.array(earth.probabilities_3nu_earth(\n'
             '        H_VAC_3NU, e*GEV, -0.8))\n'
             '    n = earth.slabs_for_tolerance(\n'
             '        H_VAC_3NU, e*GEV, -0.8, atol=atol)\n'
             '    print("%9.0f %14.2e %11d"\n'
             '          % (e, np.max(np.abs(default-converged)), n))'),
        md('The middle column is the point: at the default eight '
           'sub-slabs the error spans a factor of some four hundred '
           'across these energies alone — and this is one zenith angle. '
           'No single `n` gives a single accuracy. The right column is '
           'what a fixed accuracy actually costs, and it varies '
           'eightfold over the same range.\n\n'
           'Passing the tolerance directly does the same search and then '
           'the evaluation. `return_n_slabs` reports what it settled on, '
           'which is worth asking for: a tight tolerance can be '
           'expensive, and that should be visible rather than inferred '
           'from a slow notebook.'),
        code('prob, n = earth.probabilities_3nu_earth(\n'
             '    H_VAC_3NU, 10.0*GEV, -0.8, atol=1.0e-6, '
             'return_n_slabs=True)\n'
             'print("atol=1e-6 needed n = %d;  P_mumu = %.8f" % (n, '
             'prob[4]))\n\n'
             '# The tolerance binds on *every* returned probability, not '
             'on one\n'
             '# channel, so the subdivision is set by the '
             'worst-converging entry.\n'
             '# A tolerance that cannot be reached raises rather than '
             'quietly\n'
             '# returning less.\n'
             'try:\n'
             '    earth.slabs_for_tolerance(H_VAC_3NU, 10.0*GEV, -0.8,\n'
             '                              atol=1.0e-16, n_max=64)\n'
             'except ValueError as error:\n'
             '    print("\\nunreachable tolerance ->", '
             'str(error).split(";")[0])'),
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
           'experiment measures.\n\n'
           '**This used to have to loop, and no longer does.** Before '
           '1.12.0 `probabilities_3nu_earth` took one energy and one '
           'zenith angle, so a grid meant a Python loop over every point '
           'of it, and this notebook said so. Both arguments may now be '
           'arrays, and they broadcast: index the energies and the angles '
           'on different axes and the whole map comes back on the grid, '
           'in one call.\n\n'
           'What has *not* changed is `slabs.probabilities_3nu_slabs`, '
           'which still expects one Hamiltonian per slab rather than a '
           'stack of them per slab. The number and position of the slabs '
           'depend on the trajectory, so a set of energies through '
           '*different* zenith angles is genuinely not a rectangular '
           'array of Hamiltonians. `earth` deals with that by taking one '
           'angle at a time internally — every energy at a given angle '
           'shares the geometry and the matter potentials, which is where '
           'the saving comes from.'),
        code('n_e, n_c = 60, 60\n'
             'E_gev = np.logspace(0.0, 2.0, n_e)\n'
             'cz = np.linspace(-0.999, -0.05, n_c)\n\n'
             '# One call: energies down one axis, angles across the '
             'other.  The\n'
             '# trailing axis holds the nine probabilities; index 4 is '
             'P_mumu.\n'
             'grid = earth.probabilities_3nu_earth(\n'
             '    H_VAC_3NU, (E_gev*GEV)[:, None], cz[None, :],\n'
             '    n_slabs_per_segment=4)[..., 4]\n\n'
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
        md('### One subdivision for the whole map\n\n'
           'The map above fixes `n_slabs_per_segment=4`, chosen the way '
           'such numbers usually are — by eye. `slabs_for_tolerance` '
           'chooses it instead, from a stated accuracy. Handed the whole '
           'grid it returns one subdivision that covers every point of '
           'it, set by the worst-converging one; and because the search '
           'costs about twice a single evaluation, it is worth calling '
           'once for the map rather than per point.'),
        code('atol = 1.0e-4\n'
             'n = earth.slabs_for_tolerance(\n'
             '    H_VAC_3NU, (E_gev*GEV)[:, None], cz[None, :], '
             'atol=atol)\n'
             'print("one subdivision for the whole %dx%d map at '
             'atol=%.0e:  n = %d"\n'
             '      % (n_e, n_c, atol, n))\n\n'
             'refined = earth.probabilities_3nu_earth(\n'
             '    H_VAC_3NU, (E_gev*GEV)[:, None], cz[None, :], '
             'n)[..., 4]\n'
             'print("largest change against the n=4 map above:  %.2e"\n'
             '      % np.max(np.abs(refined-grid)))'),
        md('That difference is worth reading twice. The map above, drawn '
           'at `n_slabs_per_segment=4`, differs from the refined one by '
           'a few times $10^{-2}$ in $P_{\\mu\\mu}$ — visible on the '
           'colour scale, and far larger than anything an analysis would '
           'accept. It is fine for a picture of where the resonance '
           'sits, which is what this notebook wants; it is not fine for '
           'fitting. Choosing the subdivision from a tolerance rather '
           'than by eye is how that distinction stops being '
           'invisible.'),
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
    'and wins the most, and **the compiled backend**, which since 1.13.0 is '
    'installed with the package and needs nothing asked of it.\n\n'
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
        md('The two agree to round-off — the same arithmetic, organised '
           'differently. Not *bit for bit*: the batched path sums and '
           'multiplies in a different order, and floating-point addition is '
           'not associative, so the last digit or two can differ. The '
           'number below is the honest measure of that.'),
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
           '`numba` comes with the package as of 1.13.0, so the batched paths '
           'hand large stacks to a compiled kernel without anything being '
           'installed on purpose. It is switchable at runtime, which is how '
           'the regression suite checks the two paths against each other — '
           'and is what makes the comparison below possible in one '
           'process.'),
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
        md('## Four flavors cost more\n\n'
           '`oscprob4nu` broadcasts exactly like its three-flavor sibling, '
           'but two things make it dearer per element.\n\n'
           '1. **More algebra.** Fifteen generators rather than eight, a '
           'quartic characteristic equation rather than a cubic.\n'
           '2. **The root strategy.** `oscprob4nu.ROOT_STRATEGY` defaults '
           'to `\'double-double\'`, which is what keeps a stiff 3+1 '
           'spectrum accurate and costs roughly a fifth more than the '
           'eigensolver route behind it.\n\n'
           'Both flavor counts have a compiled kernel, so the ratio is '
           'measured twice below: once with both paths on NumPy and once '
           'with both compiled. Taking it across the two — four-flavor '
           'NumPy against three-flavor Numba — would blame the algebra for '
           'a difference in backend, which is the mistake this section had '
           'to warn about while four flavors had no kernel.'),
        code('H4_VAC = hamiltonians4nu.'
             'hamiltonian_4nu_vacuum_energy_independent(\n'
             '    gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF,\n'
             '    np.sqrt(0.10), np.sqrt(0.10), 0.0,\n'
             '    gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF, 1.0)\n\n'
             'n4 = 20000\n'
             'E4 = np.logspace(-1.0, 1.5, n4)*GEV\n'
             'H4 = hamiltonians4nu.hamiltonian_4nu_matter(\n'
             '    H4_VAC, E4, gd.VCC_EARTH_CRUST, gd.VNC_EARTH_CRUST)\n'
             'H3 = hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '    H_VAC_3NU, E4, gd.VCC_EARTH_CRUST)\n\n'
             '# Warm both paths before timing either\n'
             'oscprob3nu.probabilities_3nu(H3, L)\n'
             'oscprob4nu.probabilities_4nu(H4, L)\n\n'
             'fastkernels.USE_NUMBA = True\n'
             't3_numba = best_of(lambda: '
             'oscprob3nu.probabilities_3nu(H3, L))\n'
             't4_numba = best_of(lambda: '
             'oscprob4nu.probabilities_4nu(H4, L))\n'
             'fastkernels.USE_NUMBA = False\n'
             't3_numpy = best_of(lambda: '
             'oscprob3nu.probabilities_3nu(H3, L))\n'
             't4 = best_of(lambda: oscprob4nu.probabilities_4nu(H4, L))\n'
             'oscprob4nu.POLISH_ROOTS = False\n'
             't4_raw = best_of(lambda: '
             'oscprob4nu.probabilities_4nu(H4, L))\n'
             'oscprob4nu.POLISH_ROOTS = True\n'
             'fastkernels.USE_NUMBA = True\n\n'
             'print("%d energies" % n4)\n'
             'print("  3 flavors, NumPy path      : %7.1f ms" '
             '% (t3_numpy*1e3))\n'
             'print("  4 flavors, NumPy, polished : %7.1f ms" '
             '% (t4*1e3))\n'
             'print("  4 flavors, NumPy, no polish: %7.1f ms" '
             '% (t4_raw*1e3))\n'
             'print("  3 flavors, compiled kernel : %7.1f ms" '
             '% (t3_numba*1e3))\n'
             'print("  4 flavors, compiled kernel : %7.1f ms" '
             '% (t4_numba*1e3))\n'
             'print("\\n  4nu / 3nu, both on NumPy   : %.1fx" '
             '% (t4/t3_numpy))\n'
             'print("  4nu / 3nu, both compiled   : %.1fx" '
             '% (t4_numba/t3_numba))\n'
             'print("  what the kernel buys at 4nu: %.1fx" '
             '% (t4/t4_numba))\n'
             'print("  the refinement is           : %.0f%% of the "\n'
             '      "four-flavor NumPy call" % (100*(t4-t4_raw)/t4))'),
        md('So four flavors costs several times three on the NumPy path — '
           'the algebra and the refinement — but noticeably less than that '
           'once both are compiled. That is not the algebra getting '
           'cheaper: a good part of the four-flavor penalty was never '
           'algebra at all but the fixed cost of driving it as forty '
           'whole-array passes, a batched determinant and a Newton step '
           'among them, and the kernel does not pay that cost.\n\n'
           'Switching the refinement off buys back roughly a third of the '
           'NumPy call and gives up a factor of two on the roots. Changing '
           '`ROOT_STRATEGY` away from its default gives up rather more — an '
           'order of magnitude — for roughly a fifth of the time. Notebook '
           '16 measures both against a fifty-digit oracle, and is also '
           'where the catch lives: on a stiff spectrum neither choice is '
           'visible in the *probabilities*, because the accumulated phase '
           'dominates them.'),
        md('## What a tolerance costs\n\n'
           '`rtol` and `atol` are an accuracy feature rather than a '
           'speed one, but they have a cost worth knowing. The search '
           'evaluates the chord at a doubling sequence of subdivisions, '
           'which comes to roughly twice one evaluation at the '
           'subdivision it settles on — the series being dominated by '
           'its last term. Called once for a scan that is nothing. '
           'Called per point, it is paid at every point.'),
        code('import earth\n\n'
             'E_scan = np.logspace(0.0, 2.0, 200)*GEV\n'
             'atol = 1.0e-5\n\n'
             'n = earth.slabs_for_tolerance(H_VAC_3NU, E_scan, -0.8, '
             'atol=atol)\n'
             't_search = best_of(lambda: earth.slabs_for_tolerance(\n'
             '    H_VAC_3NU, E_scan, -0.8, atol=atol), repeat=3)\n'
             't_scan = best_of(lambda: earth.probabilities_3nu_earth(\n'
             '    H_VAC_3NU, E_scan, -0.8, n))\n\n'
             '# The trap, on a tenth of the points so the notebook '
             'still finishes.\n'
             'sample = E_scan[:20]\n'
             't_each = best_of(lambda: [earth.probabilities_3nu_earth(\n'
             '    H_VAC_3NU, e, -0.8, atol=atol) for e in sample], '
             'repeat=2)\n\n'
             'print("subdivision chosen           : %d" % n)\n'
             'print("search, once for the scan    : %7.1f ms" '
             '% (t_search*1e3))\n'
             'print("the scan itself, at that n   : %7.1f ms" '
             '% (t_scan*1e3))\n'
             'print("together                     : %7.1f ms"\n'
             '      % ((t_search+t_scan)*1e3))\n'
             'print("per point, scaled to 200     : %7.1f ms"\n'
             '      % (t_each*1e3*len(E_scan)/len(sample)))'),
        md('## What the palindrome is worth\n\n'
           'A neutrino crossing a spherically symmetric Earth meets '
           'every radius on the way in and again on the way out, so slab '
           '$j$ and slab $n-1-j$ carry the same Hamiltonian, the same '
           'width, and therefore the same operator. Composing from the '
           'centre outwards computes each of them once instead of '
           'twice.\n\n'
           'It is on by default and needs nothing from the caller. It is '
           'switchable because the two orderings round differently — by '
           'a few times $10^{-15}$ — so turning it off is how to ask for '
           'the plain left-to-right product when a comparison needs '
           'one.'),
        code('E_long = np.logspace(0.0, 2.0, 2000)*GEV\n\n'
             'fastkernels.USE_PALINDROME = True\n'
             't_on = best_of(lambda: earth.probabilities_3nu_earth(\n'
             '    H_VAC_3NU, E_long, -0.8))\n'
             'with_it = earth.probabilities_3nu_earth(H_VAC_3NU, '
             'E_long, -0.8)\n\n'
             'fastkernels.USE_PALINDROME = False\n'
             't_off = best_of(lambda: earth.probabilities_3nu_earth(\n'
             '    H_VAC_3NU, E_long, -0.8))\n'
             'without = earth.probabilities_3nu_earth(H_VAC_3NU, '
             'E_long, -0.8)\n\n'
             'fastkernels.USE_PALINDROME = True     # back to the '
             'default\n\n'
             'print("2000-energy Earth scan, three flavors")\n'
             'print("  composing the whole chord : %7.2f ms" '
             '% (t_off*1e3))\n'
             'print("  using the palindrome      : %7.2f ms" '
             '% (t_on*1e3))\n'
             'print("  speedup                   : %7.2fx" '
             '% (t_off/t_on))\n'
             'print("  largest difference        : %.2e"\n'
             '      % np.max(np.abs(with_it-without)))'),
        md('## What to take away\n\n'
           '1. Replace loops with array arguments. This is the large win, it '
           'needs no extra dependency, and the results agree to round-off.\n'
           '2. The compiled backend needs nothing installed — `numba` is a '
           'dependency as of 1.13.0. It is used at three and four flavors, '
           'and declined at two, where NumPy is already as quick.\n'
           '3. Do not bother for a handful of points: the library already '
           'takes the quicker path there on its own.\n'
           '4. At four flavors, expect several times the per-element cost '
           'of three — and leave `POLISH_ROOTS` alone unless you have '
           'measured that you can afford not to.'),
    ])

# ------------------------------------------------------- 10 the paper figures
books['10_paper_figures.ipynb'] = notebook(
    "The paper's figures",
    r'''Every figure in [arXiv:1904.12391](https://arxiv.org/abs/1904.12391) is produced here, and only here.

The plotting style is the group's standard `matplotlibrc`, inlined into the setup cell below rather than shipped as a separate file, with its sizes set for figures included at the paper's `\columnwidth`.

Running this notebook writes all eight PDFs. Set `NUOSC_PAPER_FIGDIR` to write them straight into the paper's directory.''',
    [
        code(r'''import earth
import slabs
from scipy.linalg import expm
import json
import shutil
from matplotlib.patches import (Rectangle, FancyArrowPatch,
                                Wedge, Polygon, Circle)
import matplotlib.patheffects as pe

# The paper's figures are drawn at \columnwidth -- 229.5 pt = 3.184 in for its
# two-column elsarticle layout -- so that each PDF is included at 1:1 and the
# labels come out at the size they are set at here.
COLW = 229.5/72.27

# The group's standard matplotlibrc, inlined so that the paper's figures do not
# depend on a resource file outside this repository.  The style is unchanged --
# Palatino set through LaTeX, ticks inward on all four sides, framed legends
# with a black border, Type 42 fonts -- but the sizes, which that file chooses
# for standalone 5-inch figures, are set here for figures included at
# \columnwidth beside 10 pt body text.
#
# Setting the labels through LaTeX needs a TeX installation, and Palatino comes
# from that installation rather than from matplotlib.  A machine that only wants
# to check the notebook still runs -- CI, most obviously -- has no reason to
# carry one, so this asks rather than assumes.  With LaTeX the figures are
# exactly the paper's; without it the labels fall back to matplotlib's own
# mathtext, which changes how they are typeset and nothing about what is
# plotted.  Every construct these figures use, \mathbb included, renders under
# both.
HAVE_LATEX = shutil.which('latex') is not None

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Palatino'],
    'text.usetex': HAVE_LATEX,
    'font.size': 8,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7.5,
    'axes.linewidth': 0.7,
    'lines.linewidth': 1.1,
    'xtick.top': True, 'xtick.bottom': True, 'xtick.direction': 'in',
    'ytick.left': True, 'ytick.right': True, 'ytick.direction': 'in',
    'xtick.major.size': 3.5, 'xtick.minor.size': 2.0,
    'ytick.major.size': 3.5, 'ytick.minor.size': 2.0,
    'xtick.major.width': 0.7, 'xtick.minor.width': 0.7,
    'ytick.major.width': 0.7, 'ytick.minor.width': 0.7,
    'xtick.major.pad': 3.0, 'ytick.major.pad': 3.0,
    'legend.frameon': True, 'legend.edgecolor': 'black',
    'legend.framealpha': 1.0, 'legend.fancybox': False,
    'legend.borderpad': 0.4, 'legend.handlelength': 1.9,
    'legend.handletextpad': 0.5, 'legend.columnspacing': 0.9,
    'legend.labelspacing': 0.35,
    # The shared setup cell above turns the grid on for the exploratory
    # notebooks; the paper's figures, like the group's matplotlibrc, leave
    # it off.
    'axes.grid': False,
    # Axes end exactly at the data in x.  y is left to autoscale: forcing 0-1
    # flattens the three-flavor panels, whose probabilities span a narrow range.
    'axes.xmargin': 0.0,
    'figure.dpi': 150, 'savefig.dpi': 300,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

if not HAVE_LATEX:
    # Palatino is supplied by the TeX installation, so asking for it without
    # one earns a missing-font warning per figure and DejaVu Serif anyway.
    plt.rcParams['font.serif'] = plt.rcParamsDefault['font.serif']
    print('latex not found: labels fall back to mathtext, and the serif face '
          'to matplotlib default.  The figures are the paper figures; the '
          'typesetting is not.')

# One-panel figures are square; the stacked ones keep the aspect the paper used.
FIGSIZE_SQUARE = (COLW, COLW)
FIGSIZE_STACK = (COLW, COLW*991.458/565.587)
STYLES = ['-', '--', ':', '-.']

def framed(ax, **kw):
    """A legend with a black frame at the weight of the axes spines."""
    leg = ax.legend(**kw)
    leg.get_frame().set_linewidth(0.7)
    return leg

# Set NUOSC_PAPER_FIGDIR to write the PDFs into the paper's directory.
FIGDIR = os.environ.get('NUOSC_PAPER_FIGDIR', '.')'''),
        md('''## Exactness, against an independent reference

The paper's opening figure, and its central claim. The probabilities computed
by the closed forms are compared with direct numerical exponentiation of the
*same* Hamiltonian, `scipy.linalg.expm`, at two, three and four flavors. The
lower panel is the difference. Agreement is at the level of double-precision
round-off, which is what "exact" means here.'''),
        code(r'''L_LBL = 1300.0*KM      # long baseline, for the two- and three-flavor cases
L_SBL = 0.6*KM         # short baseline, natural for an eV-scale sterile
E_gev = np.logspace(-0.5, 1.5, 300)
E = E_gev*GEV

h2 = hamiltonians2nu.hamiltonian_2nu_matter(
    np.asarray(hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(
        gd.S23_NO_BF, gd.D31_NO_BF)), E, gd.VCC_EARTH_CRUST)
h3 = hamiltonians3nu.hamiltonian_3nu_nsi(
    H_VAC_3NU, E, gd.VCC_EARTH_CRUST, gd.EPS_3)
h4 = hamiltonians4nu.hamiltonian_4nu_matter(
    hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(
        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF,
        np.sqrt(0.10), np.sqrt(0.10), 0.0,
        gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF, 1.0),
    E, gd.VCC_EARTH_CRUST, gd.VNC_EARTH_CRUST)

p2 = oscprob2nu.probabilities_2nu(h2, L_LBL)[:, 0]
p3 = oscprob3nu.probabilities_3nu(h3, L_LBL)[:, 0]
p4 = oscprob4nu.probabilities_4nu(h4, L_SBL)[:, 0]


def by_expm(h_stack, baseline):
    """P_ee from direct numerical exponentiation of the same Hamiltonian."""
    return np.array([abs(expm(-1j*np.asarray(h)*baseline)[0, 0])**2.0
                     for h in np.asarray(h_stack)])


r2 = by_expm(h2, L_LBL)
r3 = by_expm(h3, L_LBL)
r4 = by_expm(h4, L_SBL)
for n, p, r in ((2, p2, r2), (3, p3, r3), (4, p4, r4)):
    print("%d flavors: max |closed form - expm| = %.1e"
          % (n, np.abs(p-r).max()))

fig, (ax, axr) = plt.subplots(
    2, 1, figsize=(COLW, COLW*1.06), sharex=True,
    gridspec_kw={'height_ratios': [2.6, 1.0], 'hspace': 0.06})

for p, r, lab, ls in ((p2, r2, "Two flavors", "-"),
                      (p3, r3, "Three flavors", "--"),
                      (p4, r4, "Four flavors", ":")):
    # The reference curve lies under the closed-form one to round-off, so
    # plotting it would only overdraw; the lower panel carries the comparison.
    ax.semilogx(E_gev, p, ls, label=lab)
    axr.loglog(E_gev, np.abs(p-r) + 1.0e-18, ls, lw=0.9)

ax.set_ylabel(r"$P_{\nu_e \to \nu_e}$")
ax.set_ylim(0.0, 1.0)
framed(ax, loc="lower right")
axr.axhline(2.2e-16, color="0.45", ls=":", lw=0.7)
axr.set_ylabel(r"$|\Delta P|$")
axr.set_xlabel("Neutrino energy [GeV]")
axr.set_xlim(E_gev[0], E_gev[-1])
axr.set_ylim(1.0e-18, 1.0e-13)
axr.set_yticks([1.0e-17, 1.0e-15])
fig.savefig(os.path.join(FIGDIR, "validation.pdf"))
plt.show()'''),
        md('''## Three flavors, four scenarios

Vacuum, constant-density matter, matter with non-standard interactions, and a
CPT-odd Lorentz invariance-violating background, at $L = 1300$ km. All four use
the default parameters in `globaldefs`, so the figure is reproducible from the
repository alone.'''),
        code(r'''L = 1300.0*KM
E_gev = np.logspace(np.log10(0.5), np.log10(30.0), 400)
E = E_gev*GEV

p = {
    "Vacuum": oscprob3nu.probabilities_3nu(
        np.asarray(H_VAC_3NU)/E[:, None, None], L),
    "Matter": oscprob3nu.probabilities_3nu(
        hamiltonians3nu.hamiltonian_3nu_matter(
            H_VAC_3NU, E, gd.VCC_EARTH_CRUST), L),
    "NSI": oscprob3nu.probabilities_3nu(
        hamiltonians3nu.hamiltonian_3nu_nsi(
            H_VAC_3NU, E, gd.VCC_EARTH_CRUST, gd.EPS_3), L),
    "CPT-odd LIV": oscprob3nu.probabilities_3nu(
        hamiltonians3nu.hamiltonian_3nu_liv(
            H_VAC_3NU, E, gd.SXI12, gd.SXI23, gd.SXI13,
            gd.DXICP, gd.B1, gd.B2, gd.B3, gd.LAMBDA), L),
}

panels = [(0, r"$P_{\nu_e \to \nu_e}$"),
          (3, r"$P_{\nu_\mu \to \nu_e}$"),
          (4, r"$P_{\nu_\mu \to \nu_\mu}$")]

fig, axes = plt.subplots(3, 1, figsize=FIGSIZE_STACK, sharex=True)
for ax, (k, ylabel) in zip(axes, panels):
    for (name, prob), ls in zip(p.items(), STYLES):
        ax.semilogx(E_gev, prob[:, k], ls, label=name)
    ax.set_ylabel(ylabel)
    ax.set_xlim(E_gev[0], E_gev[-1])
framed(axes[0], loc="lower right")
axes[-1].set_xlabel("Neutrino energy [GeV]")
fig.tight_layout(pad=0.3)
fig.savefig(os.path.join(FIGDIR, "prob_3nu_vs_energy_compare.pdf"))
plt.show()'''),
        md('''## How the slab composition works

A diagram rather than a calculation: the neutrino crosses slabs of differing
width and density, each solved exactly, and the resulting operators are
multiplied. The ties cross, which is the point --- the slab crossed *first* is
applied first, and so stands *rightmost* in the product.'''),
        code(r'''W = [0.9, 1.7, 0.7, 1.1]          # slab widths, arbitrary units
SHADE = [0.16, 0.80, 0.52, 0.30]  # shading, standing in for density
edges = np.concatenate(([0.0], np.cumsum(W)))
total = edges[-1]
cmap = plt.get_cmap("Blues")
fill = [cmap(s) for s in SHADE]
ink = [cmap(np.clip(s + 0.46, 0.62, 0.95)) for s in SHADE]

fig, ax = plt.subplots(figsize=(COLW, COLW*0.56))
ax.set_axis_off()
ax.set_xlim(-0.60, total + 0.60)
ax.set_ylim(-1.16, 0.80)

for k in range(4):
    ax.add_patch(Rectangle((edges[k], 0.0), W[k], 0.40,
                           facecolor=fill[k], edgecolor="0.25", lw=0.6))
    xc = 0.5*(edges[k] + edges[k+1])
    ax.text(xc, 0.48, r"$\mathbb{H}_%d$" % (k+1), ha="center", va="bottom",
            fontsize=7.5)
    ax.text(xc, 0.20, r"$\rho_%d$" % (k+1), ha="center", va="center",
            fontsize=7, color="white" if SHADE[k] > 0.5 else "0.15")
    ax.plot([edges[k], edges[k+1]], [-0.15, -0.15], color="0.35", lw=0.6)
    for e in (edges[k], edges[k+1]):
        ax.plot([e, e], [-0.12, -0.18], color="0.35", lw=0.6)
    ax.text(xc, -0.24, r"$L_%d$" % (k+1), ha="center", va="top", fontsize=7)

for x0, dx in ((-0.52, 0.46), (total + 0.06, 0.46)):
    ax.add_patch(FancyArrowPatch((x0, 0.20), (x0 + dx, 0.20),
                                 arrowstyle="-|>", mutation_scale=6,
                                 lw=0.8, color="0.15"))
ax.text(-0.56, 0.20, r"$\nu_\alpha$", ha="right", va="center", fontsize=8)
ax.text(total + 0.56, 0.20, r"$\nu_\beta$", ha="left", va="center",
        fontsize=8)

y_eq = -0.98
ax.text(0.02, y_eq, r"$\mathbb{U} =$", ha="left", va="center", fontsize=8)
step = (total - 0.78)/4.0
fx = 0.78 + np.arange(4)*step
for k in range(4):
    j = 3 - k                     # factor 4 leftmost, factor 1 rightmost
    ax.text(fx[k], y_eq, r"$\mathbb{U}_%d(L_%d)$" % (j+1, j+1),
            ha="left", va="center", fontsize=8, color=ink[j])
    ax.annotate("", xy=(fx[k] + 0.5*step*0.62, y_eq + 0.14),
                xytext=(0.5*(edges[j] + edges[j+1]), -0.44),
                arrowprops=dict(arrowstyle="-|>", mutation_scale=5,
                                lw=0.6, color=ink[j], alpha=0.85,
                                shrinkA=0, shrinkB=0))

fig.savefig(os.path.join(FIGDIR, "slabs_composition.pdf"))
plt.show()'''),
        md('''## An oscillogram through the Earth

Muon-neutrino survival across energy and arrival direction, through the PREM
density profile. Each chord is decomposed into slabs, each solved exactly.'''),
        code(r'''n_e, n_c = 150, 150
E_gev = np.logspace(0.0, 2.0, n_e)
cz = np.linspace(-0.999, -0.05, n_c)

grid = np.empty((n_e, n_c))
for i, e in enumerate(E_gev):
    for j, c in enumerate(cz):
        grid[i, j] = earth.probabilities_3nu_earth(
            H_VAC_3NU, e*GEV, c, n_slabs_per_segment=4)[4]

fig, (axe, ax) = plt.subplots(
    2, 1, figsize=(COLW, COLW*1.78),
    gridspec_kw={"height_ratios": [1.30, 1.9], "hspace": 0.02})

# --- Earth in cutaway, with the PREM layers and a few chords.
# Coastlines are the public-domain GSHHS crude polygons, vendored beside
# this notebook so that the figure needs no mapping library.
with open('coastlines_crude.json') as handle:
    COASTS = json.load(handle)['polygons']

LAT0, LON0 = np.deg2rad(22.0), np.deg2rad(12.0)


def ortho(lon_deg, lat_deg):
    """Orthographic projection; the third return is the visibility mask."""
    lon, lat = np.deg2rad(lon_deg), np.deg2rad(lat_deg)
    cosc = (np.sin(LAT0)*np.sin(lat)
            + np.cos(LAT0)*np.cos(lat)*np.cos(lon - LON0))
    return (np.cos(lat)*np.sin(lon - LON0),
            np.cos(LAT0)*np.sin(lat) - np.sin(LAT0)*np.cos(lat)
            * np.cos(lon - LON0),
            cosc)


R = gd.EARTH_RADIUS
r_edge = np.concatenate(([0.0], earth.PREM_BOUNDARIES, [R]))
r_mid = 0.5*(r_edge[:-1] + r_edge[1:])
rho_shell = np.array([earth.density_prem(r) for r in r_mid])
shade = (rho_shell - rho_shell.min())/(rho_shell.max() - rho_shell.min())

axe.add_patch(Circle((0.0, 0.0), 1.0, facecolor="#bcd7e8", edgecolor="none",
                     zorder=1))
for seg in COASTS:
    seg = np.asarray(seg)
    x, y, cosc = ortho(seg[:, 0], seg[:, 1])
    if (cosc > 0.0).sum() < 3:
        continue
    axe.add_patch(Polygon(np.column_stack([x[cosc > 0.0], y[cosc > 0.0]]),
                          closed=True, facecolor="#a98055",
                          edgecolor="#8a6640", lw=0.15, zorder=2))
axe.add_patch(Circle((0.0, 0.0), 1.0, facecolor="none", edgecolor="0.35",
                     lw=0.7, zorder=6))

# The quarter cutaway, top right, showing the PREM layers.
cmap_e = plt.get_cmap("YlOrRd")
for rb, sh in zip(r_edge[1:][::-1], shade[::-1]):
    axe.add_patch(Wedge((0.0, 0.0), rb/R, 0.0, 90.0, zorder=3, lw=0.3,
                        facecolor=cmap_e(0.10 + 0.78*sh), edgecolor="0.40"))
axe.plot([0.0, 0.0], [0.0, 1.0], color="0.35", lw=0.6, zorder=4)
axe.plot([0.0, 1.0], [0.0, 0.0], color="0.35", lw=0.6, zorder=4)

for c, col in zip((-1.0, -0.8, -0.5, -0.2),
                  ["#d62728", "#1f77b4", "#ff7f0e", "#2ca02c"]):
    phi = 2.0*np.arcsin(min(1.0, -c))
    axe.annotate("", xy=(0.0, 1.0), xytext=(np.sin(phi), np.cos(phi)),
                 zorder=7,
                 arrowprops=dict(arrowstyle="-|>", mutation_scale=6, lw=1.2,
                                 color=col, shrinkA=0, shrinkB=2,
                                 path_effects=[
                                     pe.Stroke(linewidth=2.4,
                                               foreground="white"),
                                     pe.Normal()]))
    axe.plot([], [], "-", color=col, lw=1.2, label=r"$%.1f$" % c)
axe.plot(0.0, 1.0, "v", color="k", ms=4.2, zorder=8)
axe.text(0.0, 1.125, r"$\nu$", color="k", fontsize=8, ha="center",
         va="center", zorder=8)
axe.set_aspect("equal")
axe.set_axis_off()
axe.set_xlim(-1.10, 2.78)
axe.set_ylim(-1.10, 1.21)
leg_e = axe.legend(loc="center right", fontsize=7, borderpad=0.45,
                   handlelength=1.6, handletextpad=0.5, labelspacing=0.4,
                   title="Neutrino direction,\n" + r"$\cos\theta_z$")
leg_e.get_frame().set_linewidth(0.7)
leg_e.get_title().set_fontsize(7)
leg_e.get_title().set_multialignment("center")

mesh = ax.pcolormesh(cz, E_gev, grid, shading="auto", cmap="magma",
                     vmin=0.0, vmax=1.0, rasterized=True)
ax.set_yscale("log")
ax.set_xlabel(r"Cosine of zenith angle, $\cos\theta_z$")
ax.set_ylabel("Neutrino energy [GeV]")
cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
cbar.set_label(r"$P_{\nu_\mu \to \nu_\mu}$ through Earth "
               r"(standard $3\nu$)")
cbar.outline.set_linewidth(0.7)
fig.savefig(os.path.join(FIGDIR, "earth_oscillogram.pdf"))
plt.show()'''),
        md('''## The arrangement of matter, not only its mean

Four density profiles with the *same* mean density, and indeed the same
multiset of slabs, give different probabilities. The upper panel shows the
profiles; the lower one, what they do.'''),
        code(r'''def probabilities_profile(widths_km, densities, energy_ev, h_vac=None):
    """Nine probabilities through an arbitrary matter profile."""
    h_vac = H_VAC_3NU if h_vac is None else h_vac
    vcc = earth.matter_potential(np.asarray(densities, dtype=float))
    H = hamiltonians3nu.hamiltonian_3nu_matter(h_vac, energy_ev, vcc)
    return np.array(slabs.probabilities_3nu_slabs(
        H, np.asarray(widths_km, dtype=float)*KM))


total_km = 6000.0
n_slab = 24
rho_lo, rho_hi = 2.0, 8.0
mean_rho = 0.5*(rho_lo + rho_hi)

widths = np.full(n_slab, total_km/n_slab)
castle = np.where(np.arange(n_slab) % 2 == 0, rho_lo, rho_hi)
serrated = np.tile(np.linspace(rho_lo, rho_hi, 6), 4)
rng = np.random.default_rng(20260801)
# A *permutation* of the castle wall, not a fresh random draw: that guarantees
# the identical multiset of slabs, and so exactly the same mean density, which
# is the whole premise of the comparison.
random_wall = rng.permutation(castle)
uniform = np.full(n_slab, mean_rho)

profiles = [("Castle wall", castle), ("Serrated", serrated),
            ("Random wall", random_wall), ("Uniform", uniform)]
E_gev = np.logspace(-0.7, 1.7, 400)

fig = plt.figure(figsize=(COLW, COLW*1.62))
# Row 4 is left empty, as a gap between the profiles and the probabilities.
gs = fig.add_gridspec(6, 1, height_ratios=[1.0, 1.0, 1.0, 1.0, 0.55, 4.2],
                      hspace=0.18)
axes_p = [fig.add_subplot(gs[i]) for i in range(4)]
ax = fig.add_subplot(gs[5])

edges = np.concatenate(([0.0], np.cumsum(widths)))
for axp, (name, rho), ls, c in zip(axes_p, profiles, STYLES,
                                   ["C0", "C1", "C2", "C3"]):
    axp.step(edges, np.concatenate((rho[:1], rho)), ls, where="pre",
             lw=1.0, color=c, label=name)
    axp.set_xlim(0.0, edges[-1])
    axp.set_ylim(0.0, rho_hi*1.35)
    axp.set_yticks([0.0, 4.0, 8.0])
    leg = axp.legend(loc="upper left", handlelength=1.4, borderpad=0.25,
                     handletextpad=0.4, fontsize=5.5)
    leg.get_frame().set_linewidth(0.6)
    if axp is not axes_p[-1]:
        axp.set_xticklabels([])
axes_p[-1].set_xlabel("Distance travelled [km]")
# One y-label, centred on the four profile panels.
pos_top = axes_p[0].get_position(fig)
pos_bot = axes_p[-1].get_position(fig)
fig.text(0.005, 0.5*(pos_top.y1 + pos_bot.y0), r"Density [g cm$^{-3}$]",
         rotation="vertical", va="center", ha="left", fontsize=9)

for (name, rho), ls, c in zip(profiles, STYLES, ["C0", "C1", "C2", "C3"]):
    p = np.array([probabilities_profile(widths, rho, e*GEV)[3]
                  for e in E_gev])
    ax.semilogx(E_gev, p, ls, color=c, label=name)
ax.set_xlim(E_gev[0], E_gev[-1])
ax.set_xlabel("Neutrino energy [GeV]")
ax.set_ylabel(r"$P_{\nu_\mu \to \nu_e}$")
leg = ax.legend(loc="upper left", ncol=1, fontsize=7.5)
leg.get_frame().set_linewidth(0.7)
fig.savefig(os.path.join(FIGDIR, "density_arrangement.pdf"))
plt.show()'''),
        md('''## A 3+1 sterile resonance through the Earth

Four flavors, antineutrinos, through the Earth: the matter resonance that a
three-flavor treatment cannot produce at all.'''),
        code(r'''DM41 = 1.0                       # [eV^2]
S14 = S24 = np.sqrt(0.10)
S34 = 0.0

H_VAC_4NU = hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(
    gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, S14, S24, S34,
    gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF, DM41)

n_e, n_c = 300, 300
E_tev = np.logspace(-0.5, 1.5, n_e)          # 0.3 - 30 TeV
costhz = np.linspace(-1.0, -0.05, n_c)

# Antineutrinos: conjugate the vacuum term and flip both potentials.  One
# average mantle density here, so the whole grid is a single broadcast call.
rho_mantle = 4.5                              # [g cm^-3]
vcc = earth.matter_potential(rho_mantle)
H_bar = hamiltonians4nu.hamiltonian_4nu_matter(
    np.conj(H_VAC_4NU), E_tev[:, None]*1.0e3*GEV,
    -vcc, -gd.VNC_EARTH_CRUST*rho_mantle/3.0)

L_grid = np.array([earth.distance_traveled_inside_earth(c)
                   for c in costhz])*KM
grid = oscprob4nu.probabilities_4nu(H_bar, L_grid[None, :])

fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
mesh = ax.pcolormesh(costhz, E_tev, grid[:, :, 5], shading="auto",
                     cmap="viridis", vmin=0.0, vmax=1.0, rasterized=True)
ax.set_yscale("log")
ax.set_xlabel(r"Cosine of zenith angle, $\cos\theta_z$")
ax.set_ylabel("Antineutrino energy [TeV]")
cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
cbar.set_label(r"$P_{\bar\nu_\mu \to \bar\nu_\mu}$ through Earth "
               r"($3+1\nu$)")
cbar.outline.set_linewidth(0.7)
fig.tight_layout(pad=0.3)
fig.savefig(os.path.join(FIGDIR, "sterile_earth_oscillogram.pdf"))
plt.show()'''),
        md(r'''## How the pieces fit together

A map of the library: what each module is responsible for, and how a
Hamiltonian becomes a probability.'''),
    code(r'''from matplotlib.patches import FancyBboxPatch

# Two columns wide: TEXTW is the paper's \textwidth.
TEXTW = 522.0/72.27

fig, axd = plt.subplots(figsize=(TEXTW, TEXTW*0.34))
axd.set_axis_off()
axd.set_xlim(0.0, 30.0)
axd.set_ylim(0.0, 10.2)

C_IN, C_CORE, C_ACC, C_COMP = "#e8eef6", "#dce9dc", "#f6ecd9", "#efe2ee"
E_IN, E_CORE, E_ACC, E_COMP = "#4a6fa5", "#4a7a4a", "#b07d2a", "#8a5a86"


def box(x, y, w, h, face, edge, title, body):
    axd.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.12,rounding_size=0.25",
        facecolor=face, edgecolor=edge, lw=0.9, zorder=2))
    axd.text(x + w/2.0, y + h - 0.42, title, ha="center", va="top",
             fontsize=7.6, color=edge, zorder=3,
             fontfamily="monospace" if title[0].islower() else None)
    if body:
        axd.text(x + w/2.0, y + h - 1.30, body, ha="center", va="top",
                 fontsize=6.0, color="0.25", zorder=3, linespacing=1.5)


def arrow(x0, y0, x1, y1, color="0.35", ls="-"):
    axd.annotate("", xy=(x1, y1), xytext=(x0, y0), zorder=1,
                 arrowprops=dict(arrowstyle="-|>", mutation_scale=8, lw=0.9,
                                 color=color, linestyle=ls, shrinkA=2,
                                 shrinkB=2))


# Left: what goes in.  Middle: the library proper.  Right: layered matter.
box(0.3, 5.6, 7.4, 3.9, C_IN, E_IN, "Your Hamiltonian",
    "Any Hermitian $2{\\times}2$, $3{\\times}3$\nor $4{\\times}4$ array, or a whole\nstack of them")
box(0.3, 0.5, 7.4, 3.9, C_IN, E_IN, "hamiltonians$n$nu",
    "Vacuum, matter, NSI, LIV,\n3+1 --- the worked scenarios,\nbuilt from globaldefs")

# The green box was 5.0 tall for four lines of body, leaving its lower half
# empty, and its rounded corner met the orange box's: 2.6 against 2.5, minus
# 0.12 of padding each, overlapped them.  Shortened, and moved up to leave a
# clear 0.5 between the two.
box(9.6, 3.4, 9.4, 4.05, C_CORE, E_CORE, "oscprob$n$nu",
    "Invariants from traces\n"
    "$\\rightarrow$ eigenvalues in closed form\n"
    "$\\rightarrow$ coefficients $u_0,\\, u_k$\n"
    "$\\rightarrow$ the probabilities")
box(9.6, 0.5, 9.4, 2.2, C_ACC, E_ACC, "fastkernels",
    "Compiled kernels: the same\nnumbers, up to $12{\\times}$ faster")

box(21.0, 5.6, 8.7, 3.9, C_COMP, E_COMP, "slabs",
    "Composes the operators of\nadjacent constant-density\nlayers, in order")
box(21.0, 0.5, 8.7, 3.9, C_COMP, E_COMP, "earth",
    "PREM chords, named sites\nand zenith scans --- it\nbuilds the slabs")

arrow(7.9, 7.5, 9.4, 6.0)
arrow(7.9, 2.4, 9.4, 4.2)
arrow(14.3, 3.35, 14.3, 2.85, color=E_ACC, ls="--")
arrow(19.2, 6.0, 20.8, 7.0, color=E_COMP)
arrow(20.8, 4.6, 19.2, 4.4, color=E_COMP)
arrow(25.3, 4.5, 25.3, 5.5, color=E_COMP)
fig.savefig(os.path.join(FIGDIR, "architecture.pdf"),
            bbox_inches="tight", pad_inches=0.02)
plt.show()'''),
    md(r'''## Speed and accuracy, against five other codes

Every code here solves the *same* problem and is refereed by the *same*
arbitrary-precision reference, so that none of them is its own judge.

The external numbers are frozen in `tests/speed_accuracy.json`; the cell below
records exactly how each was obtained, so the comparison can be reproduced or
disputed without guessing what was run.'''),
    code(r'''# ---------------------------------------------------------------------------
# Provenance of the external numbers.  None of these codes is needed to run
# this notebook: the results are frozen in tests/speed_accuracy.json.  What
# follows is exactly what was done to produce them, on Linux, August 2026.
#
# TWO CONVENTIONS ARE MATCHED FIRST, or the comparison measures bookkeeping
# rather than physics:
#
#   Matter potential.  This library takes n_e = Y_e rho / m_bar, with m_bar the
#   mean free-nucleon mass; the others take rho N_A Y_e.  The ratio is the
#   nuclear mass defect,
#       scale = (CONV_G_TO_EV/((MASS_PROTON+MASS_NEUTRON)/2)) / N_A = 0.99209
#   and each external code is handed rho*scale = 2.97628 g/cm^3 so that all
#   propagate with the same V_CC.  (For GLoBES this factor was *found* by
#   scanning the density, and the minimum landed on 0.99209 to a part in 1e7 --
#   an independent check that the reasoning above is right.)
#
#   Length unit.  The codes disagree on hbar*c in the 7th digit:
#       this library  5.0677300000e9   eV^-1 per km
#       nuSQuIDS      5.0677307162e9
#       NuFast-LBL    5.0677302143e9   ( = 1e3 / 1.97327e-7 )
#   so each is given the baseline in km that reproduces *our* L in eV^-1.
#   NuFast and Prob3++ and GLoBES get 1299.999945 km instead of 1300.
#
# THE REFEREE.  Accuracy is measured against a 50-digit mpmath matrix
# exponential of the same Hamiltonian:
#     mp.mp.dps = 50;  U = mp.expm(-1j*M*mp.mpf(L));  P = abs(U[0,1])**2
# This library agrees with it to 4e-17, which is why it can be plotted at
# "double precision" rather than at zero.
#
# ---------------------------------------------------------------------------
# nuSQuIDS 1.13.3 -- manylinux wheel, no compiler needed:
#     pip install nusquids            # imports as `nuSQuIDS`
#   Multiple-energy mode, Set_rel_error = Set_abs_error swept over
#   1e-4 ... 1e-12.  Timed through the Python interface, including building
#   the solver and reading the results back.  Note its own default tolerance
#   is far tighter and ~40x slower than an explicit 1e-12.
#
# NuFast-LBL -- header-only C++, not on PyPI:
#     git clone https://github.com/PeterDenton/NuFast-LBL.git
#     g++ -Ofast -ffast-math driver.cpp -o driver
#   Probability_Matter_LBL(...), sweeping N_Newton = 0,1,2,3.  Timed inside
#   C++, best of five, with no Python in the loop.  rhoYe = 1.488737 (that is
#   3*0.5*0.99209).  N_Newton = 0 is its default and is a truncation, not an
#   error: two steps reach 8e-12.
#
# GLoBES 3.2.18 -- needs GSL (2.7.1 here):
#     curl -O https://www.mpi-hd.mpg.de/personalhomes/globes/download/globes-3.2.18.tar.gz
#     tar xzf globes-3.2.18.tar.gz && cd globes-3.2.18
#     ./configure --prefix=$PREFIX --with-gsl-prefix=/usr && make && make install
#     gcc -O2 -I$PREFIX/include drv.c -L$PREFIX/lib -lglobes -lgsl -lgslcblas -lm
#   glbConstantDensityProbability(2, 1, +1, E, L, rho), after glbDefineParams
#   with angles as arcsin(sqrt(sin^2 theta)) and glbSetOscillationParameters.
#
# Prob3++ -- C++; mosc.c and mosc3.c must be compiled as C, not C++:
#     git clone https://github.com/rogerwendell/Prob3plusplus.git
#     gcc -O2 -c mosc.c mosc3.c ; g++ -O2 -c EarthDensity.cc BargerPropagator.cc
#     g++ -O2 drv.cpp mosc.o mosc3.o EarthDensity.o BargerPropagator.o
#   BargerPropagator::SetMNS(s12sq, s13sq, s23sq, dm21, dmAtm, dcp, E, true, 1)
#   then propagateLinear(1, L, rho) and GetProb(2, 1).
#   CAUTION: the fifth argument is Dm2_32, NOT Dm2_31.  Passing Dm2_31 leaves a
#   1.7e-2 disagreement that no density rescaling will remove; with
#   Dm2_32 = 2.525e-3 - 7.39e-5 it drops to 4.5e-5.
#
# nuCraft r22 -- pure Python, installs by unpacking:
#     curl -O https://nucraft.hepforge.org/nuCraft-r22.tar.gz && tar xzf ...
#   NOT included here: its interface is (PDG type, energy, zenith angle) and it
#   propagates through its own Earth model, with no constant-density mode, so
#   it belongs in an Earth comparison rather than this one.  Two of its
#   conventions would need matching there: default electron fractions
#   y = (0.4957, 0.4656, 0.4656) rather than 0.5, and a default rICore of
#   1121.5 km against the 1221.5 km in its own PREM density table.
# ---------------------------------------------------------------------------

import json

with open(os.path.join('..', 'tests', 'speed_accuracy.json')) as handle:
    sa = json.load(handle)

# The same colours, markers and layout as the two Earth planes, so that a
# code keeps its identity across all three.  This library shows only the
# batched-plus-kernel point: the other two routes are what the performance
# figure is for, and here they would be three labels on one horizontal line.
STYLE = {"NuOscProbExact": ("-o", "C3", 4.0),
         "nuSQuIDS": ("-v", "C2", 3.6),
         "NuFast-LBL": ("-D", "C4", 3.2),
         "GLoBES": ("-*", "C6", 6.0),
         "Prob3++": ("-P", "C5", 4.4),
         "Second-order expansion": ("-s", "C1", 3.6)}
DRAWN = {"NuOscProbExact": ("Array + kernel",)}

fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
for series in sa["series"]:
    marker, colour, size = STYLE[series["name"]]
    keep = DRAWN.get(series["name"])
    points = [q for q in series["points"]
              if keep is None or q["label"] in keep]
    t = [q["us_per_probability"] for q in points]
    e = [q["max_abs_error"] for q in points]
    kw = dict(ms=size, color=colour, label=series["name"], zorder=4)
    if series["name"] == "NuOscProbExact":
        kw.update(mfc="white", mew=1.0, zorder=5)
    ax.loglog(t, e, marker, **kw)

# What is varied along each curve, at both ends where there are two.
for lab, x, y, dx, dy, c, ha in (
        (r"tol $10^{-4}$", 51.61, 1.83e-04, 0, 7, "C2", "center"),
        (r"$10^{-12}$", 189.37, 1.89e-08, 7, -1, "C2", "left"),
        (r"$N_{\rm Newton} = 0$", 0.044, 1.53e-05, 7, -2, "C4", "left"),
        (r"$3$", 0.066, 8.30e-12, 7, -2, "C4", "left")):
    ax.annotate(lab, xy=(x, y), xytext=(dx, dy), textcoords="offset points",
                fontsize=5.2, color=c, ha=ha)

ax.axhline(2.2e-16, color="0.5", ls=":", lw=0.7, zorder=1)
ax.text(2.6e-2, 3.2e-16, "Double precision", fontsize=5.4, color="0.4")
# This panel's data does not leave the same corners free as the Earth ones:
# the single NuOscProbExact point sits at the bottom left and nuSQuIDS runs
# along the top right, so the subtitle goes just above the former and the
# legend into the empty bottom right.
ax.text(0.03, 0.17, "Constant density:  $L = 1300$ km,\n"
        r"$E = 0.6$--$20$ GeV,  $\rho = 3$ g cm$^{-3}$",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=6.0,
        color="0.2", linespacing=1.4)
ax.annotate("Array + kernel", xy=(0.230, 9.71e-16), xytext=(8, -1),
            textcoords="offset points", fontsize=5.2, color="C3", ha="left")
ax.set_xlabel(r"Time per probability [$\mu$s]")
ax.set_ylabel(r"Error vs.\ a 50-digit reference,  "
              r"max $|\Delta P_{\nu_\mu \to \nu_e}|$")
ax.set_xlim(2.0e-2, 1.0e3)
ax.set_ylim(1.0e-16, 1.0e-2)
leg = ax.legend(loc="lower right", bbox_to_anchor=(0.995, 0.02),
                fontsize=6.0)
leg.get_frame().set_linewidth(0.7)
fig.tight_layout(pad=0.3)
fig.savefig(os.path.join(FIGDIR, "speed_accuracy.pdf"))
plt.show()'''),
    md(r'''## Speed and accuracy through the Earth

The same plane as above, on the PREM profile rather than at constant density,
at three flavors and at 3+1.

These two measure something different from the figure above, and the
difference is the whole point of them. At constant density there is an exact
answer, and the figure measures how closely each code evaluates it: the floor
is round-off. Through the Earth the density varies continuously, and *every*
code approximates the profile, this one included. So the reference here is a
**converged** PREM solution rather than an exact one, and each code sits where
its treatment of the profile puts it, not where the accuracy of its
probability formula would.

The consequence is that this library is a curve rather than a horizontal line,
with `n_slabs_per_segment` as its dial. That is more informative than a line:
the dial is what a user actually turns, and it keeps converging — second order
in the slab width — through the region where both external codes flatten out.

The external numbers are frozen in `tests/prem_speed_accuracy.json`.'''),
    code(r'''# ---------------------------------------------------------------------------
# Provenance.  Neither external code is needed to run this notebook: the
# results are frozen in tests/prem_speed_accuracy.json.  The generator,
# tests/prem_scan.py, carries the full account in its module docstring --
# every convention matched, how each was established, and the three
# Python 2 remnants in nuCraft r22.  In brief:
#
#   THE REFEREE.  The slab product converges at SECOND order in the slab
#   width: the density is sampled at slab midpoints and the PREM shell
#   boundaries are cut exactly, so the leading error is O(h^2).  Measured
#   exponents are 2.000 at 3, 10 and 40 GeV.  The reference is therefore
#       ref = (4*P(256) - P(128))/3
#   of a 30-digit mpmath slab product -- not the first-order 2*P(32)-P(16),
#   which for this problem amplifies the error instead of cancelling it.
#
#   REFEREEING THE REFEREE.  That extrapolation shares earth.earth_slabs
#   with the code it judges, so it is checked against a discretisation with
#   nothing in common: an adaptive DOP853 integration of
#   dpsi/dx = -i H(x) psi through the continuous profile, restarted at each
#   shell boundary so that no step straddles a discontinuity.  The two
#   agree to 2.3e-11 at three flavors and 1.2e-9 at 3+1.
#
#   nuSQuIDS 1.13.3.  Atmosphere height set to zero (it defaults to 22 km,
#   which lengthens the chord by 24.4 km); cos(theta_z) shrunk by 1.413e-7
#   so that its km reproduces ours in eV^-1 -- in vacuum that takes the
#   agreement from 3.8e-7 to 1.0e-15; this library's PREM handed over
#   through the documented (x, rho, ye) constructor at ye = 0.5 with the
#   usual 0.99209 mass defect, on a grid whose points land exactly on the
#   PREM shell boundaries.  Its ~1e-6 floor is not the integrator: it is
#   identical for RK4, RKF45, RKCK, RK8PD and MSADAMS, for any h_max, and
#   for every tolerance below 1e-9.
#
#   nuCraft r22.  CalcWeights takes the ZENITH ANGLE IN RADIANS while
#   ConstructMixingMatrix takes its angles in degrees; passing degrees to
#   both leaves a disagreement of 0.198.  atmMode 0 with atmHeight 0 and
#   detectorDepth 0 is surface to surface.  Its V_CC constant, 15.256e-5,
#   is the atomic-mass-unit value rather than the mean-nucleon one its own
#   source comment claims, so its density is scaled by 0.9926737; scanning
#   for the minimum residual returns 0.9926748, and the two codes then
#   agree to 3e-11 at constant density.
#
#   nuCraft is absent from the 3+1 panel, and not because it cannot do 3+1.
#   Its sterile and charged-current entries come from two independently
#   rounded constants whose ratio is 0.5016 where the isoscalar value is
#   exactly 0.5.  Scaling the density cannot separate them, so a 3.7e-4
#   floor survives every setting; forcing the ratio to 0.5 by hand drops the
#   same run to 2.8e-7.  Publishing the patched curve would misrepresent the
#   released code and publishing the unpatched one would misrepresent its
#   solver, so it appears only where the sterile entry never enters.
# ---------------------------------------------------------------------------

with open(os.path.join('..', 'tests', 'prem_speed_accuracy.json')) as handle:
    prem = json.load(handle)

# The colours and markers of the constant-density plane, so that a code
# keeps its identity between the two figures.
# The tolerance dial is the primary curve: it is the one where the user
# states the quantity the y-axis measures, and the one a new reader meets
# first.  The explicit slab count is the same code told the answer instead
# of asked to find it, and is drawn broken behind it.
PREM_STYLE = {"NuOscProbExact (tolerance)": ("-o", "C3", 4.0),
              "NuOscProbExact": ("--s", "C3", 3.2),
              # The style is keyed off the frozen name, not the legend
              # label, so the four-flavor slab curve needs the broken style
              # under its own name too.
              "NuOscProbExact (double-double)": ("--s", "C3", 3.2),
              "NuOscProbExact (eigensolver)": (":^", "C1", 3.8),
              "nuSQuIDS": ("-v", "C2", 3.6),
              "NuFast-Earth": ("-D", "C4", 3.2),
              "GLoBES": ("-*", "C6", 6.0),
              "Prob3++": ("-P", "C5", 4.4),
              "nuCraft": ("-s", "C0", 3.4)}


def prem_plane(panel, annotations, subtitle, xlim, ylim, outfile,
               dial_at=None, legend_loc="lower left",
               legend_anchor=None, only=None, relabel=None):
    """One speed-accuracy plane: time across, error against the referee up.

    `only` restricts which series are drawn, and `relabel` renames them in
    the legend.  Both exist for the four-flavor panel, where the two root
    strategies are measured and frozen but only the default is drawn: the
    curves coincide to the last bit, so plotting both put two labels on one
    line.
    """
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    wanted = {}
    relabel = relabel or {}
    order = list(PREM_STYLE)
    for series in sorted(panel["series"],
                         key=lambda x: order.index(x["name"])):
        if only is not None and series["name"] not in only:
            continue
        marker, colour, size = PREM_STYLE[series["name"]]
        t = [q["us_per_probability"] for q in series["points"]]
        e = [q["max_abs_error"] for q in series["points"]]
        kw = dict(ms=size, color=colour,
                  label=relabel.get(series["name"], series["name"]), zorder=4)
        if series["name"].startswith("NuOscProbExact"):
            kw.update(mfc="white", mew=1.0, zorder=5)
        ax.loglog(t, e, marker, **kw)
        for q in series["points"]:
            key = relabel.get(series["name"], series["name"])
            wanted[(key, q["label"])] = (
                q["us_per_probability"], q["max_abs_error"], colour)

    for (name, label), text, dx, dy, ha in annotations:
        x, y, colour = wanted[(name, label)]
        ax.annotate(text, xy=(x, y), xytext=(dx, dy),
                    textcoords="offset points", fontsize=5.2, color=colour,
                    ha=ha)

    # What the numbers along this library's curve are, said once -- but only
    # where they are unambiguous.  In the three-flavor panel three other
    # codes are also dialled by a shell count, so the caption says it there.
    if dial_at is not None:
        ax.text(dial_at[0], dial_at[1], "Slabs per\nsegment", fontsize=5.2,
                color="C0", ha="center", linespacing=1.3)
    # Bottom left is the one corner no curve reaches in either panel, and
    # the legend needs the top right: with four to seven entries it was
    # otherwise sitting on top of nuCraft.
    ax.text(0.03, 0.03, subtitle, transform=ax.transAxes, ha="left",
            va="bottom", fontsize=6.0, color="0.2", linespacing=1.4)
    ax.set_xlabel(r"Time per probability [$\mu$s]")
    ax.set_ylabel(r"Error vs.\ converged solution,  "
                  r"max $|\Delta P_{\nu_\mu \to \nu_\mu}|$")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    leg = ax.legend(loc=legend_loc, bbox_to_anchor=legend_anchor,
                    fontsize=6.0)
    leg.get_frame().set_linewidth(0.7)
    fig.tight_layout(pad=0.3)
    fig.savefig(os.path.join(FIGDIR, outfile))
    plt.show()


prem_plane(
    prem["three_flavor"],
    [(("NuOscProbExact, rtol", "3e+00"), r"rtol $=3$", 7, -2, "left"),
     (("NuOscProbExact, rtol", "1e-05"), r"$10^{-5}$", 0, 7, "center"),
     ((r"NuOscProbExact, $N_{\rm slabs}$", "1"), r"$N_{\rm slabs}=1$",
      -7, -2, "right"),
     ((r"NuOscProbExact, $N_{\rm slabs}$", "256"), "256", 6, -3, "left"),
     (("nuSQuIDS", "1e-03"), r"tol $10^{-3}$", -3, 7, "center"),
     (("nuSQuIDS", "1e-12"), r"$10^{-12}$", 7, -1, "left"),
     (("NuFast-Earth", "1"), r"$N_{\rm shells}=1$", 4, -8, "left"),
     (("NuFast-Earth", "256"), "256", 6, -2, "left")],
    "PREM, three flavors:  " + r"$\cos\theta_z = -0.9$," + "\n"
    r"$E = 3$--$40$ GeV,  $L = 11468$ km",
    (1.0e0, 1.0e5), (2.0e-8, 2.0e-1),
    "prem_speed_accuracy.pdf", None,
    legend_loc="upper right", legend_anchor=(0.995, 0.995),
    relabel={"NuOscProbExact": r"NuOscProbExact, $N_{\rm slabs}$",
             "NuOscProbExact (tolerance)": "NuOscProbExact, rtol"})

prem_plane(
    prem["sterile_3plus1"],
    [(("NuOscProbExact, rtol", "3e+00"), r"rtol $=3$", 7, -2, "left"),
     (("NuOscProbExact, rtol", "1e-05"), r"$10^{-5}$", 0, 7, "center"),
     ((r"NuOscProbExact, $N_{\rm slabs}$", "1"), r"$N_{\rm slabs}=1$",
      0, 7, "center"),
     ((r"NuOscProbExact, $N_{\rm slabs}$", "256"), "256", 6, -3, "left"),
     (("nuSQuIDS", "1e-03"), r"tol $10^{-3}$", -3, 7, "center"),
     (("nuSQuIDS", "1e-12"), r"$10^{-12}$", 7, -1, "left")],
    r"PREM, $3+1$:  $\cos\theta_z = -0.9$," + "\n"
    r"$E = 0.3$--$30$ TeV,  $\Delta m_{41}^2 = 1$ eV$^2$",
    (1.0e1, 1.0e5), (2.0e-8, 2.0e-1),
    "prem_speed_accuracy_3plus1.pdf", None,
    legend_loc="upper right", legend_anchor=(0.995, 0.995),
    # Both root strategies are frozen in the data; only the default is
    # drawn, because the two agree to the last bit and their curves lie on
    # top of one another.
    only={"NuOscProbExact (double-double)", "NuOscProbExact (tolerance)",
          "nuSQuIDS", "nuCraft"},
    relabel={"NuOscProbExact (double-double)":
             r"NuOscProbExact, $N_{\rm slabs}$",
             "NuOscProbExact (tolerance)": "NuOscProbExact, rtol"})'''),
    md(r'''## Performance

Three ways of evaluating the same scan: one point at a time, the whole stack
in one call, and the whole stack through the compiled kernel. The cost per
point is what a parameter scan actually pays.'''),
    code(r'''import time
import fastkernels


def best_of(func, repeat=5):
    # The minimum, not the mean: timing noise is one-sided, so the fastest
    # run is the one least polluted by whatever else the machine was doing.
    best = float("inf")
    for _ in range(repeat):
        t0 = time.perf_counter()
        func()
        best = min(best, time.perf_counter()-t0)
    return best


L = 1300.0*KM
sizes = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000]
rows = []
for n in sizes:
    Hs = hamiltonians3nu.hamiltonian_3nu_matter(
        H_VAC_3NU, np.logspace(-1.0, 1.5, n)*GEV, gd.VCC_EARTH_CRUST)
    fastkernels.USE_NUMBA = False
    t_loop = best_of(lambda: [oscprob3nu.probabilities_3nu(h, L)
                              for h in Hs], repeat=3)
    t_array = best_of(lambda: oscprob3nu.probabilities_3nu(Hs, L))
    if fastkernels.HAVE_NUMBA:
        fastkernels.USE_NUMBA = True
        oscprob3nu.probabilities_3nu(Hs, L)          # warm the compiler
        t_numba = best_of(lambda: oscprob3nu.probabilities_3nu(Hs, L))
        fastkernels.USE_NUMBA = False
    else:
        t_numba = np.nan
    rows.append((n, t_loop, t_array, t_numba))
rows = np.array(rows)
fastkernels.USE_NUMBA = fastkernels.HAVE_NUMBA

print("%8s %12s %12s %12s" % ("N", "loop [us/pt]", "array [us/pt]",
                              "numba [us/pt]"))
for n, tl, ta, tn in rows:
    print("%8d %12.3f %12.3f %12.3f"
          % (n, tl/n*1e6, ta/n*1e6, tn/n*1e6))
print("\nat the largest N: array is %.0fx the loop, numba a further %.1fx"
      % (rows[-1, 1]/rows[-1, 2], rows[-1, 2]/rows[-1, 3]))

with open(os.path.join('..', 'tests', 'timing_other_codes.json')) as handle:
    other = json.load(handle)

fig, (axt, ax) = plt.subplots(
    2, 1, figsize=(COLW, COLW*1.30), sharex=True,
    gridspec_kw={"height_ratios": [1.0, 1.0], "hspace": 0.07})

# Two of the three routes through this library: the compiled kernel, which
# is what an installation gets, and one point at a time, which is what the
# batching is worth against.  The plain array path lies between them and is
# left out to keep the legend readable.  Marker and colour match the two
# Earth planes, so a reader tracks one code across three figures.
for col, style, size, lab in (
        (3, "-o", 3.4, "NuOscProbExact, array + kernel"),
        (1, "-^", 3.6, "NuOscProbExact, one point at a time")):
    if col == 3 and not fastkernels.HAVE_NUMBA:
        continue
    kw = dict(ms=size, color="C3", mfc="white", mew=0.9, zorder=5)
    axt.loglog(rows[:, 0], rows[:, col]*1e3, style, label=lab, **kw)
    ax.loglog(rows[:, 0], rows[:, col]/rows[:, 0]*1e6, style, **kw)

# The external codes.  None of them batches: GLoBES, Prob3++ and NuFast-LBL
# take one energy per call, so their cost per probability is flat in N, and
# that flatness is the point of the lower panel.  nuSQuIDS has a
# multiple-energy mode and still does not fall, because what it amortises is
# the solver setup and not the integration.
for key, style, col, size, lab in (
        ("nusquids", "-v", "C2", 2.6, "nuSQuIDS"),
        ("nufast_lbl", "-D", "C4", 2.4, "NuFast-LBL"),
        ("globes", "-*", "C6", 4.4, "GLoBES"),
        ("prob3pp", "-P", "C5", 3.2, "Prob3++")):
    n = np.array(other[key]["sizes"], dtype=float)
    t = np.array(other[key]["seconds"])
    axt.loglog(n, t*1e3, style, ms=size, color=col, label=lab)
    ax.loglog(n, t/n*1e6, style, ms=size, color=col)

axt.set_ylabel("Total time [ms]")
axt.text(0.02, 0.97, "Constant density:  " + r"$L = 1300$ km," + "\n"
         r"$E = 0.1$--$32$ GeV,  $\rho = 3$ g cm$^{-3}$" + "\n"
         "All nine three-flavor probabilities" + "\n"
         "All codes at std.\\ settings",
         transform=axt.transAxes, ha="left", va="top", fontsize=5.6,
         color="0.2", linespacing=1.4)
# `mode="expand"` with the anchor spanning the axes makes the legend box
# exactly as wide as the panels it sits above.
leg = axt.legend(loc="lower left", bbox_to_anchor=(0.0, 1.015, 1.0, 0.1),
                 mode="expand", ncol=2, fontsize=5.6, borderaxespad=0.0,
                 columnspacing=1.2, handlelength=2.0)
leg.get_frame().set_linewidth(0.7)
ax.set_xlim(rows[0, 0], rows[-1, 0])
ax.set_xlabel("Number of energies in the scan")
ax.set_ylabel(r"Time per probability [$\mu$s]")
ax.set_ylim(2.0e-2, 1.0e3)
fig.savefig(os.path.join(FIGDIR, "performance.pdf"), bbox_inches="tight")
plt.show()'''),
    md(r'''## Comparison to the alternatives

Three ways of getting the same number at $L = 1300$ km, in constant matter:
the closed form here; **nuSQuIDS**, an independently written C++ code that
integrates the density matrix numerically; and the standard second-order
expansion in $\alpha = \Delta m^2_{21}/\Delta m^2_{31}$ and $s_{13}$ used
throughout long-baseline physics (Cervera *et al.*; Akhmedov *et al.*).

The nuSQuIDS curve is read from a frozen JSON file, generated separately by
`tests/nusquids_scan.py`, so this notebook does not need nuSQuIDS installed.'''),
    code(r'''# P(nu_mu -> nu_e) to second order in alpha and s13, constant density:
# Akhmedov, Johansson, Lindner, Ohlsson & Schwetz, hep-ph/0402175, Eq. (32).
# Evaluated directly as the analytic expression it is -- NuOscProbExact is
# not involved, and could not be: it returns the exact answer.
def p_mue_expansion(E_ev, L, s12, s13, s23, dcp, dm21, dm31, vcc):
    c12 = np.sqrt(1.0 - s12**2.0)
    c23 = np.sqrt(1.0 - s23**2.0)
    alpha = dm21/dm31
    Delta = dm31*L/(4.0*E_ev)
    A = 2.0*vcc*E_ev/dm31
    om = 1.0 - A
    return (4.0*s13**2.0*s23**2.0*np.sin(om*Delta)**2.0/om**2.0
            + alpha*8.0*s13*s12*c12*s23*c23*np.cos(Delta + dcp)
            * np.sin(A*Delta)/A*np.sin(om*Delta)/om
            + alpha**2.0*4.0*s12**2.0*c12**2.0*c23**2.0
            * np.sin(A*Delta)**2.0/A**2.0)


with open(os.path.join('..', 'tests', 'nusquids_scan.json')) as handle:
    ref = json.load(handle)
with open(os.path.join('..', 'tests', 'nufast_scan.json')) as handle:
    nuf = json.load(handle)

E_gev = np.array(ref['energy_gev'])
E = E_gev*GEV
L_km = ref['baseline_km']
rho = ref['density_g_cm3']
L = L_km*KM
vcc = earth.matter_potential(rho)

exact = oscprob3nu.probabilities_3nu(
    hamiltonians3nu.hamiltonian_3nu_matter(H_VAC_3NU, E, vcc), L)[:, 3]
nusq = np.array(ref['probability'])
nufast = np.array(nuf['probability_N_Newton_0'])
approx = p_mue_expansion(E, L, gd.S12_NO_BF, gd.S13_NO_BF, gd.S23_NO_BF,
                         gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF, vcc)

print("max |exact - nuSQuIDS|  = %.1e" % np.abs(exact-nusq).max())
print("max |exact - NuFast|    = %.1e" % np.abs(exact-nufast).max())
print("max |exact - expansion| = %.1e" % np.abs(exact-approx).max())

fig, (ax, axr) = plt.subplots(
    2, 1, figsize=(COLW, COLW*1.06), sharex=True,
    gridspec_kw={"height_ratios": [2.4, 1.2], "hspace": 0.06})

ax.semilogx(E_gev, exact, "-", color="C0", label="NuOscProbExact")
ax.semilogx(E_gev[::4], nusq[::4], "o", color="C2", ms=2.4, mfc="none",
            mew=0.7, label="nuSQuIDS")
ax.semilogx(E_gev[::4], nufast[::4], "s", color="C4", ms=2.2, mfc="none",
            mew=0.7, label="NuFast-LBL")
ax.semilogx(E_gev, approx, "--", color="C1", label=r"$\alpha$--$s_{13}$ expansion")
ax.set_ylabel(r"$P_{\nu_\mu \to \nu_e}$")
ax.set_ylim(0.0, 0.16)
framed(ax, loc="upper right", ncol=1)

axr.loglog(E_gev, np.abs(exact-nusq), "-", color="C2", lw=0.9,
           label="nuSQuIDS")
axr.loglog(E_gev, np.abs(exact-nufast), "-.", color="C4", lw=0.9,
           label="NuFast-LBL")
axr.loglog(E_gev, np.abs(exact-approx), "--", color="C1", lw=0.9,
           label="expansion")
axr.set_ylabel(r"$|\Delta P|$")
axr.set_xlabel("Neutrino energy [GeV]")
axr.set_xlim(E_gev[0], E_gev[-1])
axr.set_ylim(1.0e-10, 1.0e-1)
axr.set_yticks([1.0e-9, 1.0e-6, 1.0e-3])
fig.savefig(os.path.join(FIGDIR, "exact_vs_approximations.pdf"))
plt.show()'''),
        md(r'''## Two flavors, four scenarios

The two-neutrino counterpart of the three-flavor figure, driven by
$\Delta m^2_{31}$ and $\theta_{23}$, over a slightly wider energy range.

Note the first argument: the builders take $\sin\theta$, not $\theta$.'''),
        code(r'''E_gev = np.logspace(np.log10(0.5), np.log10(40.0), 400)
E = E_gev*GEV
L = 1300.0*KM

H2_VAC = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(
    gd.S23_NO_BF, gd.D31_NO_BF)

p2 = {
    "Vacuum": oscprob2nu.probabilities_2nu(
        np.asarray(H2_VAC)/E[:, None, None], L),
    "Matter": oscprob2nu.probabilities_2nu(
        hamiltonians2nu.hamiltonian_2nu_matter(
            H2_VAC, E, gd.VCC_EARTH_CRUST), L),
    "NSI": oscprob2nu.probabilities_2nu(
        hamiltonians2nu.hamiltonian_2nu_nsi(
            H2_VAC, E, gd.VCC_EARTH_CRUST, gd.EPS_2), L),
    "CPT-odd LIV": oscprob2nu.probabilities_2nu(
        hamiltonians2nu.hamiltonian_2nu_liv(
            H2_VAC, E, gd.SXI12, gd.B1, gd.B3, gd.LAMBDA), L),
}

panels2 = [(0, r"$P_{\nu_e \to \nu_e}$"),
           (2, r"$P_{\nu_\mu \to \nu_e}$"),
           (3, r"$P_{\nu_\mu \to \nu_\mu}$")]

fig, axes = plt.subplots(3, 1, figsize=FIGSIZE_STACK, sharex=True)
for ax, (k, ylabel) in zip(axes, panels2):
    for (name, prob), ls in zip(p2.items(), STYLES):
        ax.semilogx(E_gev, prob[:, k], ls, label=name)
    ax.set_ylabel(ylabel)
    ax.set_xlim(E_gev[0], E_gev[-1])
    ax.set_ylim(0.0, 1.0)
framed(axes[0], loc="lower right")
axes[-1].set_xlabel("Neutrino energy [GeV]")
fig.tight_layout(pad=0.3)
fig.savefig(os.path.join(FIGDIR, "prob_2nu_vs_energy_compare.pdf"))
plt.show()'''),
    ])


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
        md('The reason is in the two combinations of $\\theta_{23}$ that the '
           'channels carry, and it is worth doing the arithmetic rather '
           'than asserting it:'),
        code('for x in (0.45, 0.55):\n'
             '    print("sin^2(th23) = %.2f  ->  sin^2(2 th23) = %.4f"\n'
             '          % (x, 4.0*x*(1.0-x)))'),
        md('$\\sin^2 2\\theta_{23}$ is the *same* for both — that is the '
           'degeneracy — while $\\sin^2\\theta_{23}$ differs by about 20%. '
           'Disappearance is governed by the first at leading order and '
           'appearance by the second, so the leading term cancels in one '
           'channel and not the other.\n\n'
           'What that does *not* mean is that the disappearance curves lie '
           'on top of each other. Over the range plotted above the largest '
           'separation is about 0.017 in $P_{\\mu\\mu}$ and 0.013 in '
           '$P_{\\mu e}$ — comparable in absolute terms. The difference is '
           '*where* they sit: the disappearance gap is concentrated near '
           'the oscillation minimum, where $P_{\\mu\\mu}$ has fallen to a '
           'few per cent and there are correspondingly few events, whereas '
           'the appearance gap is a steady fraction of a signal that is '
           'small everywhere.\n\n'
           'So the honest statement is not that disappearance cannot see '
           'the octant, but that its sensitivity comes from subleading '
           'terms in the least favourable place, while appearance carries '
           'it at leading order. That is why the appearance channel drives '
           'the measurement.'),
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
        md('It converges at around ten thousand slabs, against an estimate of '
           'twenty thousand — the right order, which is all a '
           'back-of-the-envelope count of oscillation lengths can be '
           'expected to give. So far so good.\n\n'
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
           'two axes and let one call fill the grid.\n\n'
           'This uses a single average mantle density, which is what makes '
           'the grid one broadcast call. `earth.probabilities_4nu_earth` '
           'does the full PREM profile instead, slab by slab, one zenith '
           'angle at a time; the cell after this one checks the two against '
           'each other.'),
        code('import earth\n\n'
             'n_e, n_c = 150, 150\n'
             'E_tev = np.logspace(-0.5, 1.5, n_e)          # 0.3 - 30 TeV\n'
             'costhz = np.linspace(-1.0, -0.05, n_c)\n\n'
             '# Antineutrinos: conjugate the vacuum term and flip both\n'
             '# potentials.  One average mantle density here, so that the\n'
             '# whole grid is a single broadcast call.\n'
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
        md('### What the averaging costs\n\n'
           'The map above uses one average mantle density, which is what '
           'makes the whole grid a single broadcast call. '
           '`earth.probabilities_4nu_earth` instead cuts the chord at every '
           'PREM shell boundary, takes the density at the midpoint of each '
           'sub-slab, and multiplies the per-slab operators — the '
           'four-flavor case of the machinery notebooks 06 and 07 use at '
           'three.\n\n'
           'The comparison below is for **neutrinos**, not the antineutrinos '
           'of the map: `probabilities_4nu_earth` takes the potentials with '
           'the sign they have for neutrinos, exactly as its three-flavor '
           'sibling does. Each row uses the path-weighted mean PREM density '
           'along that very chord, so the only difference left is the '
           'averaging itself.'),
        code('print("costhz   rho_mean   averaged      PREM   difference")\n'
             'for c in (-0.2, -0.5, -0.8, -0.95, -1.0):\n'
             '    L = earth.distance_traveled_inside_earth(c)*KM\n'
             '    w, d = earth.earth_slabs(c, 8)\n'
             '    rho = float(np.sum(w*d)/np.sum(w))\n'
             '    H_avg = hamiltonians4nu.hamiltonian_4nu_matter(\n'
             '        H_VAC_4NU, 3.0e3*GEV,\n'
             '        earth.matter_potential(rho),\n'
             '        earth.matter_potential_nc(rho))\n'
             '    averaged = oscprob4nu.probabilities_4nu(H_avg, L)[5]\n'
             '    prem = earth.probabilities_4nu_earth(H_VAC_4NU, 3.0e3*GEV, c)[5]\n'
             '    print("%6.2f   %7.2f   %8.4f  %8.4f   %+9.4f"\n'
             '          % (c, rho, averaged, prem, prem-averaged))'),
        md('Through the mantle the two agree to a few parts in a hundred, '
           'which is the level at which replacing a varying density by its '
           'mean is defensible. The chord that grazes the core is several '
           'times worse than any of them: a single number cannot stand in '
           'for a profile that jumps by a factor of two at the core-mantle '
           'boundary, however well chosen the number is. Note that the '
           'diametric chord is *not* the worst case here — it is symmetric '
           'and spends most of its length in the core, so its mean is '
           'closer to representative than the grazing one\'s. That is the '
           'kind of thing `earth` exists to get right without being '
           'guessed at.'),
        md('The dark band is the sterile matter resonance: a region of '
           'energy and zenith angle where muon antineutrinos are strongly '
           'depleted into the sterile state. Its position depends on '
           '$\\Delta m^2_{41}$, which is what makes it a search channel.'),
        md('## Sterile states and new interactions together\n\n'
           '3+1 does not have to be the only new physics in the problem. '
           '`hamiltonians4nu.hamiltonian_4nu_nsi` adds non-standard '
           'interactions on top of the sterile state, and the flavor '
           'structure is the physically interesting part: NSI are a '
           'modification of a *standard-model* interaction, and a sterile '
           'state has none — so the $\\epsilon$ matrix acts on the active '
           'block only, and the sterile row and column stay '
           'untouched.\n\n'
           'The sterile entry keeps the $-V_{NC}$ it had before, since '
           'that came from removing the neutral-current potential of the '
           'active states, not from any interaction of the sterile one.'),
        code('H_nsi = hamiltonians4nu.hamiltonian_4nu_nsi(\n'
             '    H_VAC_4NU, energy, gd.VCC_EARTH_CRUST,\n'
             '    gd.VNC_EARTH_CRUST, gd.EPS_3)\n\n'
             '# What the NSI term added, on top of standard matter\n'
             'added = np.array(H_nsi) - np.array(H)\n'
             'print("NSI contribution [eV], real part:")\n'
             'print("%8s" % "" + "".join("%12s" % f '
             'for f in flavors))\n'
             'for i, a in enumerate(flavors):\n'
             '    print("%8s" % a\n'
             '          + "".join("%12.3e" % added[i, j].real\n'
             '                     for j in range(4)))\n\n'
             'print("\\nsterile row and column untouched : %s"\n'
             '      % bool(np.all(np.abs(added[3, :]) == 0.0)\n'
             '             and np.all(np.abs(added[:, 3]) == 0.0)))'),
        code('p_nsi = np.array(oscprob4nu.probabilities_4nu(\n'
             '    H_nsi, baseline)).reshape(4, 4)\n'
             'p_std = np.array(oscprob4nu.probabilities_4nu(\n'
             '    H, baseline)).reshape(4, 4)\n\n'
             'print("%-10s %12s %12s %12s"\n'
             '      % ("channel", "3+1", "3+1 with NSI", "change"))\n'
             'for i, a in enumerate(flavors):\n'
             '    for j, b in enumerate(flavors):\n'
             '        if i == j or (i, j) not in ((0, 1), (1, 0), (1, 3)):\n'
             '            continue\n'
             '        print("P(%s->%s)%s %12.6f %12.6f %+12.6f"\n'
             '              % (a, b, " "*(3-len(a)-len(b)),\n'
             '                 p_std[i, j], p_nsi[i, j],\n'
             '                 p_nsi[i, j]-p_std[i, j]))'),
        md('The sterile channel $P(\\nu_\\mu \\to \\nu_s)$ moves too, even '
           'though the sterile state feels no NSI directly — because NSI '
           'reshape the active block, and the active and sterile states '
           'are mixed. Nothing in this calculation is a special case: it '
           'is still one Hermitian matrix handed to the same routine.'),
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
           'There are two ways out of a bottleneck like that, and the '
           'library ships both. `oscprob4nu.ROOT_STRATEGY` chooses.\n\n'
           '**`\'eigensolver\'`** stops asking the invariants: it takes the '
           'roots from `numpy.linalg.eigvalsh`, which forms no invariants at '
           'all, and then refines them against the *matrix* with one Newton '
           'step on $\\chi(\\psi) = \\det(\\psi\\mathbb{1} - \\tilde H)$, '
           'reading the Hamiltonian entries at full precision. '
           '`POLISH_ROOTS` switches that step. This was the whole answer in '
           '1.12.0.\n\n'
           '**`\'double-double\'`**, the default since 1.13.0, widens the '
           'bottleneck instead. The compression only hurts because '
           '$10^{-16}$ is what a `float64` coefficient carries; forming '
           '$I_2, I_3, I_4$ as *pairs* of `float64` — some 32 digits — '
           'leaves the same $2.3\\times10^{9}$ amplification acting on '
           '$10^{-32}$, landing at $10^{-23}$. One Aberth sweep, also in '
           'double-double, then takes the quartic\'s roots to the last '
           '`float64` bit.\n\n'
           'Here is what each costs, and — the part worth reading twice — '
           'where the difference does and does not show up.'),
        md('One wrinkle first: **`float64` cannot referee this comparison.** '
           'The two routes differ at $4\\times10^{-17}$ against '
           '$4\\times10^{-16}$, and an `eigh` reference is itself only good '
           'to $\\sim10^{-15}$. So the oracle below is the same frozen one '
           'the test suite uses — `mpmath.eig` at fifty decimal digits on '
           'nine Hamiltonians, stored as hexadecimal floats in '
           '`tests/stiff_reference.json` so that a reader gets the exact '
           'bits it was computed from. It spans the physical 3+1 range out '
           'to $\\Delta m^2_{41} = 1000$ eV$^2$ and down to a pair of roots '
           'separated by $10^{-16}$.'),
        code(r'''import json
import time
import fastkernels

with open(os.path.join('..', 'tests', 'stiff_reference.json')) as handle:
    CASES = json.load(handle)['cases']


def frozen(case, key):
    """The stored fifty-digit values, off their exact bits."""
    return np.array([float.fromhex(x) for x in case[key]])


def traceless_of(case):
    entries = np.array([[float.fromhex(case['matrix_real'][i][j])
                         + 1j*float.fromhex(case['matrix_imag'][i][j])
                         for j in range(4)] for i in range(4)])
    return entries - np.trace(entries)/4.0*np.eye(4)


# A stack large enough to time stably, stiff at every point
H_timed = hamiltonians4nu.hamiltonian_4nu_matter(
    H_VAC_4NU, np.logspace(-1.0, 1.0, 50000)*GEV,
    gd.VCC_EARTH_CRUST, gd.VNC_EARTH_CRUST)

fastkernels.USE_NUMBA = True
oscprob4nu.probabilities_4nu(H_timed[:8], baseline)      # warm the kernel

print("%-15s  %-11s  %-11s  %s"
      % ("ROOT_STRATEGY", "worst root", "worst |dP|", "|dP| at 0.1 eV^2"))
for strategy in ("double-double", "eigensolver"):
    oscprob4nu.ROOT_STRATEGY = strategy
    root_error = 0.0
    prob_error = 0.0
    benign = 0.0
    for case in CASES:
        matrix = traceless_of(case)
        # `_latent_roots` is internal; it is used here because it is the
        # only way to see the roots alone, which is where the routes differ.
        psi = np.sort(oscprob4nu._latent_roots(matrix[None, ...])[0])
        exact = frozen(case, "roots")
        root_error = max(root_error, np.max(np.abs(psi - exact))
                         / np.max(np.abs(exact)))
        p = np.asarray(oscprob4nu.probabilities_4nu(
            matrix, float.fromhex(case["baseline"]))).ravel()
        gap = np.max(np.abs(p - frozen(case, "probabilities")))
        prob_error = max(prob_error, gap)
        if case["label"].endswith("0.1 eV^2"):
            benign = gap
    print("%-15s  %-11.2e  %-11.2e  %.2e"
          % (strategy, root_error, prob_error, benign))

# Timed in ALTERNATING PAIRS, reported as the median ratio.  Timing the two
# routes once each is dominated by whatever else the machine is doing: done
# that way on a shared runner, a real 25% difference reads as parity, and
# the faster route can even come out ahead.  Pairing cancels the drift.
ratios = []
for repetition in range(9):
    order = ("double-double", "eigensolver")
    if repetition % 2:
        order = order[::-1]
    elapsed = {}
    for strategy in order:
        oscprob4nu.ROOT_STRATEGY = strategy
        t0 = time.perf_counter()
        oscprob4nu.probabilities_4nu(H_timed, baseline)
        elapsed[strategy] = time.perf_counter() - t0
    ratios.append(elapsed["double-double"]/elapsed["eigensolver"])
ratios.sort()

print()
print("cost of 'double-double' over 'eigensolver', 50k points:")
print("  median of 9 alternated pairs : %.2fx" % ratios[len(ratios)//2])
print("  spread across those pairs    : %.2fx to %.2fx"
      % (ratios[0], ratios[-1]))

oscprob4nu.ROOT_STRATEGY = "double-double"'''),
        md('Three things in that output, and the second is the one that '
           'catches people.\n\n'
           '**On the roots the default wins by about an order of '
           'magnitude** — $3.6\\times10^{-17}$ against '
           '$3.9\\times10^{-16}$, which is under a fifth of a `float64` ulp '
           'against nearly two. That is the whole claim of the '
           'double-double route, and it is measurable only against an '
           'oracle carrying more digits than either.\n\n'
           '**On probabilities the two are indistinguishable, and both are '
           'poor, on the stiffest cases.** That is not a contradiction and '
           'not a defect in either route: rebuilding $U_4$ in Newton form '
           'takes *second* differences of $e^{-i\\psi L}$ over the roots, '
           'so a root error enters the coefficients twice and the '
           'probability error grows as $(\\psi_{\\max}L)^2$ — measured at '
           '$5\\times10^{-17}(\\psi_{\\max}L)^2$ across five decades of '
           'phase. At $\\Delta m^2_{41} = 1000$ eV$^2$ over 1300 km the '
           'accumulated phase is 2.5 million radians and `float64` simply '
           'cannot carry it; every route lands in the same place. The '
           '`|dP| at 0.1` column is the physical end of the range, where '
           'the phase is a few hundred radians and the answer is good to '
           '$10^{-12}$.\n\n'
           '**The default costs roughly a fifth more**, and notice how '
           'hard that is to measure: the spread across nine alternated '
           'pairs is wide enough that a single unpaired before-and-after '
           'would have been worthless, and could easily have shown the '
           'slower route winning. The cost is unsurprising once stated — '
           'the default runs the *same* eigensolver, because it needs a '
           'start that separates the cluster and that is exactly what a '
           'backward-stable Hermitian eigensolver is for, and then adds '
           'its invariants and one Aberth sweep on top of it.\n\n'
           'So: switch to `\'eigensolver\'` if four-flavor root-finding '
           'dominates your runtime and you are content with two ulp. Leave '
           'it alone otherwise. And if you are chasing accuracy in '
           '*probabilities* at large $\\Delta m^2_{41}L/E$, neither switch '
           'is your problem — the phase is.'),
        md('For the record, the alternatives measured and rejected along the '
           'way, all against the same fifty-digit oracle.\n\n'
           '| Strategy for the roots | Relative error | Closed form? |\n'
           '|---|---|---|\n'
           '| Closed form alone | 2.2e-07 | yes |\n'
           '| `numpy.linalg.eigvalsh` alone | 6.9e-16 | no |\n'
           '| `eigvalsh` + one Newton step | 3.9e-16 | no |\n'
           '| **Invariants in double-double** | **3.6e-17** | yes |\n'
           '| Closed form in `numpy.longdouble` | 4.5e-11 | yes |\n\n'
           'The Newton step is about **twice as accurate as LAPACK** alone, '
           'because `eigvalsh` reduces the matrix by similarity transforms '
           'that each carry a backward error $\\sim\\epsilon\\|H\\|$, while '
           'the step converges onto the root of $\\det(\\psi\\mathbb{1} - '
           '\\tilde H)$ for the matrix it was handed.\n\n'
           'Extended precision is the instructive failure. It buys under one '
           'digit rather than the three its extra mantissa suggests — the '
           'cluster amplifies whatever coefficient error survives — is '
           'slower for not being hardware-vectorised, and is silently '
           '`float64` on Apple Silicon and Windows, where it would appear '
           'to work while doing nothing. Double-double gets four orders '
           'more than that, portably, out of pairs of the `float64` every '
           'machine has.\n\n'
           'Two more were tried inside the double-double route itself and '
           'lost. Starting from the closed form rather than the eigensolver '
           'is twice as fast and exact on every stiff case here, yet fails '
           'on a near-degenerate pair, where Aberth converges *linearly* at '
           'ratio one half: from $10^{-7}$ five sweeps reach '
           '$3.8\\times10^{-9}$ where thirty are needed. And Durand-Kerner, '
           'at one division per root against Aberth\'s five, is non-monotone '
           'in the sweep count — 3.9e-16, then 9.7e-17, then 1.9e-16 — so it '
           'cannot be a default.'),
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

# --------------------------------------------------------- 17 cross-checks
books['17_cross_checks.ipynb'] = notebook(
    'Cross-checks with other codes',
    'The regression suite is thorough about algebra and, unavoidably, '
    'self-consistent about **conventions**. `scipy.linalg.expm` confirms '
    'that $U$ is the matrix exponential of the Hamiltonian — but of the '
    'Hamiltonian this library built. The standard oscillation formulas were '
    'transcribed from the same papers the builders were. A mixing-matrix '
    'ordering, a sign or a unit conversion that was wrong *consistently* '
    'across the package would sit comfortably inside all of it.\n\n'
    'Corroboration from outside is what closes that gap. This notebook '
    'collects two sources of it — an independently written code, and a '
    'published closed form reached by a different route. Both agree; the '
    'interesting part is what had to be lined up first.\n\n'
    '---\n\n'
    '### The answer, up front\n\n'
    'Once two conventions are matched — a length unit and an electron '
    'density, both differences in *definition* rather than errors in '
    'either code — this library and '
    '[nuSQuIDS](https://github.com/arguelles/nuSQuIDS) agree to between '
    '$4\\times10^{-16}$ and $3\\times10^{-10}$ across eleven '
    'configurations, and the Zaglauer–Schwarzer closed form reproduces '
    'the matter spectrum to $7\\times10^{-16}$ at any configuration '
    'asked of it.\n\n'
    '| Checked | By |\n'
    '|---|---|\n'
    '| $U$ is $e^{-iHL}$ | `scipy.linalg.expm`, in the suite |\n'
    '| The Hamiltonian we build is the right one | nuSQuIDS, and '
    'Zaglauer–Schwarzer for the matter spectrum |\n'
    '| Mixing conventions, $\\theta_{23}$, $\\delta_{CP}$ | nuSQuIDS, '
    'and the vacuum reference formulas |\n'
    '| The antineutrino rule — conjugate *and* flip | nuSQuIDS, three '
    'and four flavors, vacuum and matter |\n'
    '| The mass ordering sign | nuSQuIDS, inverted-ordering case |\n'
    '| Four-flavor probabilities | nuSQuIDS, at 3e-10 to 4e-16 |\n'
    '| Arbitrary matter configurations | Zaglauer–Schwarzer, live |\n\n'
    '**Not** covered, and worth knowing before leaning on any of it: '
    'the frozen set is eleven configurations rather than a parameter '
    'space, with $\\delta_{CP}$, $\\theta_{23}$ and the sterile angles '
    'each held at one value; the Zaglauer–Schwarzer check reaches '
    'eigenvalues but not mixing; and no external code here exercises '
    '`slabs` or `earth`.\n\n'
    'The rest of the notebook is how those numbers were arrived at. The '
    'closing section returns to this table with the caveats spelled '
    'out in full.',
    [
        md('## Two sources of corroboration\n\n'
           '**[nuSQuIDS](https://github.com/arguelles/nuSQuIDS)** evolves the '
           'neutrino density matrix numerically, in C++. It shares no lineage '
           'with the closed-form SU($n$) expansions here, and it is '
           'configured from mixing *angles* and a density in g/cm³ rather '
           'than from a Hamiltonian — so agreement exercises '
           '`hamiltonians3nu` and `hamiltonians4nu` as much as the '
           'expansions. It is a considerably more general tool than this '
           'one, which is what makes the agreement worth having.\n\n'
           'The two arrive at a probability by genuinely different '
           'routes, which is the point — a shared mistake is only '
           'possible where the routes share something, and these share '
           'almost nothing:\n\n'
           '| | NuOscProbExact | nuSQuIDS |\n'
           '|---|---|---|\n'
           '| You hand it | a Hermitian matrix $H$ and a baseline | '
           'mixing angles, $\\Delta m^2$, a body and a track |\n'
           '| Eigenvalues of $H$ | in closed form, from radicals | '
           'never formed explicitly |\n'
           '| Evolution | $U = e^{-iHL}$ built by interpolating the '
           'exponential over those eigenvalues | the density matrix '
           'integrated numerically along the track |\n'
           '| You get back | $\\lvert U_{\\beta\\alpha}\\rvert^2$ | the '
           'evolved '
           'flavor state, projected |\n'
           '| Error comes from | floating-point round-off | the ODE '
           'tolerance, plus round-off |\n\n'
           'So the comparison below is not one implementation of a '
           'formula against another. It is a closed form against a '
           'numerical integration, starting from different inputs.'
           '\n\n'
           'Its values are **frozen** into `tests/nusquids_reference.json` '
           'rather than computed here. nuSQuIDS ships manylinux wheels and '
           'would install in CI, but a comparison that runs in CI ties '
           'this suite to an external project\'s release cadence. '
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
           'For each frozen case, `ours()` below rebuilds the *same '
           'physical situation* from this library: it constructs the '
           'vacuum Hamiltonian from the mixing parameters nuSQuIDS was '
           'given, adds the matter potential, and evaluates the SU(3) '
           'or SU(4) expansion. nuSQuIDS\' answer for that situation '
           'is already in the JSON.\n\n'
           'The two are then compared entry by entry, as '
           '$\\max_{\\alpha\\beta}|P^{\\rm ours}_{\\alpha\\beta} - '
           'P^{\\rm nuSQuIDS}_{\\alpha\\beta}|$ over the full 3×3 or '
           '4×4 table — a worst case, not an average.\n\n'
           'Each case is computed **twice**: once with this library\'s '
           'own constants, and once with nuSQuIDS\' — so the effect of '
           'matching conventions is shown rather than asserted. Only '
           'the constants change between the two columns; the physics '
           'and the code path are identical.'),
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
        code('fig, ax = plt.subplots(figsize=(7.6, 5.4))\n'
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
             '# Headroom above the first bar so the legend sits clear\n'
             '# of it rather than on top of it.\n'
             'ax.set_ylim(len(names) - 0.5, -2.1)\n'
             'ax.legend(loc="upper right", frameon=True,\n'
             '          edgecolor="black", framealpha=1.0)\n'
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
           'stiffest case here. That last one is the stiff-spectrum limit '
           '`oscprob4nu.ROOT_STRATEGY` documents, arrived at independently; '
           'it is attributed at the end of this notebook rather than left '
           'hanging. Note what it is *not*: since 1.13.0 the roots '
           'themselves are good to a fifth of an ulp, and what remains here '
           'is the accumulated phase in the reconstruction, which grows as '
           'its square.'),
        md('## A second check: matter eigenvalues, live\n\n'
           'The nuSQuIDS comparison above is frozen at eleven '
           'configurations. This one runs at **any** configuration, with no '
           'external dependency — but it checks something narrower, and it '
           'is worth being exact about what.\n\n'
           '### What is compared with what\n\n'
           'Both sides produce the **three eigenvalues of the same object**: '
           'the effective mass-squared matrix in matter,\n\n'
           '$$M^2_{\\rm eff} = U M^2 U^\\dagger + '
           '\\mathrm{diag}(A, 0, 0), \\qquad A = 2EV_{CC},$$\n\n'
           'in units of eV$^2$. They reach them from opposite directions:\n\n'
           '| | Our side | Zaglauer–Schwarzer |\n'
           '|---|---|---|\n'
           '| Starts from | the $3\\times3$ matrix `hamiltonians3nu` builds '
           '| the numbers $\\theta_{12}, \\theta_{13}, \\Delta m^2_{21}, '
           '\\Delta m^2_{31}, A$ |\n'
           '| Is a matrix ever formed? | yes — that is the thing being '
           'tested | **no** |\n'
           '| Route to the eigenvalues | `numpy.linalg.eigvalsh`, i.e. '
           'LAPACK | an explicit cubic, solved in radicals |\n\n'
           'The Hamiltonian this library builds is $H = M^2_{\\rm '
           'eff}/(2E)$, so multiplying it by $2E$ recovers $M^2_{\\rm eff}$ '
           'and the two sides become directly comparable. That factor is the '
           'only manipulation applied.\n\n'
           '### What this does and does not exercise\n\n'
           'It tests **`hamiltonians3nu`, not `oscprob3nu`**. The SU(3) '
           'expansion never runs here: the left-hand side goes through '
           'LAPACK, not through the closed form. What is under test is '
           'whether the *matrix we assembled* is the right one — the PMNS '
           'parametrization, the placement of $V_{CC}$ on the $ee$ entry, '
           'and the relative normalization of the vacuum and matter terms. '
           'Those are conventions, and conventions are exactly what an '
           'internal test cannot catch.\n\n'
           'It also stops at eigenvalues. Eigenvalues fix the oscillation '
           '*frequencies* but not the amplitudes, so this says nothing about '
           'the mixing — that is what nuSQuIDS and the vacuum reference '
           'formulas cover.\n\n'
           'One detail that keeps it honest: the $|U_{ei}|^2$ in the formula '
           'below are written from the mixing angles directly, as '
           '$c_{12}^2c_{13}^2$, $s_{12}^2c_{13}^2$ and $s_{13}^2$, rather '
           'than read out of `pmns_mixing_matrix`. Had they been taken from '
           'our own mixing matrix, an error in it would appear on both sides '
           'and cancel.'),
        code('def zaglauer_schwarzer(s12, s13, dm21, dm31, A):\n'
             '    """Exact eigenvalues of M^2_eff, from the vacuum '
             'parameters."""\n'
             '    c13sq = 1.0 - s13*s13\n'
             '    ue1, ue2, ue3 = (1.0 - s12*s12)*c13sq, s12*s12*c13sq, '
             's13*s13\n'
             '    # lambda^3 - x lambda^2 + y lambda - z = 0, whose three\n'
             '    # roots are the eigenvalues of M^2_eff\n'
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
             '    E = e*GEV\n\n'
             '    # OUR SIDE: build the matter Hamiltonian, then\n'
             '    # diagonalise it numerically.  H = M^2_eff/(2E), so\n'
             '    # multiplying by 2E gives back M^2_eff in eV^2.\n'
             '    # Note this goes through LAPACK, not through the\n'
             '    # SU(3) expansion -- what is under test here is the\n'
             '    # matrix, not the solver.\n'
             '    H = np.asarray(hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '        H_VAC_3NU, E, VCC))\n'
             '    mine = np.sort(np.linalg.eigvalsh(2.0*E*H))\n\n'
             '    # THEIR SIDE: the same three eigenvalues, straight\n'
             '    # from the vacuum parameters and A = 2 E V_CC.  No\n'
             '    # matrix is built at any point.\n'
             '    theirs = zaglauer_schwarzer(gd.S12_NO_BF, gd.S13_NO_BF,\n'
             '                                gd.D21_NO_BF, gd.D31_NO_BF,\n'
             '                                2.0*E*VCC)\n\n'
             '    # Compare, scaled by the size of the spectrum\n'
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
        md('The upper panel is the physics: the three eigenvalues of '
           '$M^2_{\\rm eff}$ as the energy rises, computed from our '
           'matrix. Two of them approach near 10 GeV — the MSW '
           'resonance, marked — and then swap character, which is the '
           'level crossing that drives matter effects.\n\n'
           'The lower panel is the check: at each energy, the largest '
           'difference between those eigenvalues and the '
           'Zaglauer–Schwarzer ones, relative to the spectral scale. '
           'It sits at round-off across two and a half decades of '
           'energy, including through the resonance — which is where '
           'the two eigenvalues come closest and any error in how the '
           'matrix was assembled would show first.\n\n'           'The dotted line marks $A = \\Delta m^2_{31}\\cos2\\theta_{13}$, '
           'the two-flavor resonance condition; the true minimum-gap '
           'energy is about 1% lower, the difference being the solar '
           'term that condition omits.\n\n'
           'One transcription note, because it is the kind of mistake that '
           'usually produces a plausible wrong curve rather than an obvious '
           'one: the depressed cubic\'s constant term is '
           '$-2x^3/27 + xy/3 - z$. Written first with the signs of the first '
           'two terms flipped, it disagreed with the spectrum by 50% — large '
           'enough to be unmissable, which was luck. The corrected form is '
           'identical to the resolvent cubic already inside `oscprob4nu`, '
           'which is where it should have been checked against from the '
           'start.'),
        md('## What each code is for\n\n'
           'Agreement is one question; which tool suits a given calculation '
           'is another. The two are built for different things that happen '
           'to overlap.\n\n'
           '| If you need | Reach for |\n'
           '|---|---|\n'
           '| Constant or piecewise-constant density | either — this '
           'library gives it in closed form |\n'
           '| The evolution operator returned directly | this library |\n'
           '| Two, three or four flavors | either |\n'
           '| More than four flavors | nuSQuIDS |\n'
           '| A profile varying continuously along the trajectory | '
           'nuSQuIDS, or Magnus |\n'
           '| Neutrino–nucleon interactions, attenuation | nuSQuIDS |\n'
           '| Tau regeneration | nuSQuIDS |\n'
           '| Open-system, non-unitary evolution | nuSQuIDS |\n\n'
           'On the overlap the closed form is quicker, because it evaluates '
           'an expression where nuSQuIDS integrates a density matrix. That '
           'is a difference in what each is solving rather than a measure of '
           'either — the integration is what buys the rest of the column '
           'above, and a calculation that needs any of it is not one a '
           'closed form can serve.\n\n'
           'The same reasoning marks this library\'s own boundary. Where a '
           'profile varies appreciably over an oscillation length — the Sun, '
           'adiabatic MSW — slabbing stops being practical and a '
           'Magnus-type method is the right tool. Notebook 14 shows where '
           'that begins.'),
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
           'The antineutrino and mass-ordering rows are there deliberately. '
           'Those are the two sign conventions this library has a history '
           'of getting wrong — the 1.1.0 release fixed a two-flavor vacuum '
           'Hamiltonian built with the opposite sign, invisible in vacuum '
           'and wrong the moment matter was added — so a comparison meant '
           'to catch convention errors ought to exercise them. The '
           'three-flavor cases agree at $10^{-13}$ to $10^{-15}$; the '
           'four-flavor antineutrino case is the stiffest configuration '
           'here and is discussed just below.\n\n'
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

# ------------------------------------------- 18 the operator and the algebra
books['18_evolution_operator.ipynb'] = notebook(
    'The evolution operator, and the SU($n$) coefficients',
    'Every notebook so far has asked for probabilities. This one asks for '
    'the machinery underneath — the evolution operator itself, and the '
    'expansion coefficients the method is built on.\n\n'
    'Two reasons to care. If you want to **compose** evolution across '
    'segments, or propagate a **density matrix**, or feed an oscillation '
    'calculation into something larger, a probability is a dead end and the '
    'operator is not. And if you want to see *how* the closed form works '
    'rather than take it on trust, the coefficients are where it lives.',
    [
        md('## The operator, not the probability\n\n'
           'Every `probabilities_*` routine has an `evolution_operator_*` '
           'twin taking the same arguments. It returns\n\n'
           '$$ U(L) = e^{-i\\tilde{H}L} $$\n\n'
           'with $\\tilde H$ the traceless part of the Hamiltonian, since '
           'the trace is only an overall phase.\n\n'
           'The probabilities are read off it as '
           '$P_{\\alpha\\beta} = |U_{\\beta\\alpha}|^2$ — note the index '
           'order, which is the one place this is easy to get backwards: the '
           'operator is indexed (final, initial), the probabilities '
           '(initial, final).'),
        code('energy = 1.0*GEV\n'
             'baseline = 1300.0*KM\n'
             'H = np.asarray(H_VAC_3NU)/energy\n\n'
             'U = np.array(oscprob3nu.evolution_operator_3nu(H, baseline))\n'
             'p = np.array(oscprob3nu.probabilities_3nu(H, baseline))\n\n'
             'print("U is %s, complex" % (U.shape,))\n'
             'print("unitary to        %.2e"\n'
             '      % np.max(np.abs(U.conj().T @ U - np.eye(3))))\n\n'
             '# P[3*alpha + beta] = |U[beta, alpha]|^2\n'
             'from_U = np.array([abs(U[b, a])**2\n'
             '                   for a in range(3) for b in range(3)])\n'
             'print("max |P - |U|^2|   %.2e" % np.max(np.abs(p - from_U)))'),
        md('## The group property, and why slabbing is legitimate\n\n'
           'Because $\\tilde H$ does not depend on $L$, the operator obeys\n\n'
           '$$ U(L_1 + L_2) = U(L_2)\\, U(L_1) . $$\n\n'
           'This is not a detail — it is the licence for everything the '
           '`slabs` and `earth` modules do. A trajectory through varying '
           'matter is cut into pieces over which $H$ is taken constant, each '
           'piece is solved *exactly*, and the operators are multiplied. The '
           'only approximation is the cutting, never the solving.\n\n'
           'Note the order: the operator for the **first** segment travelled '
           'goes on the **right**.'),
        code('L1, L2 = 0.4*baseline, 0.6*baseline\n\n'
             'U_whole = np.array(oscprob3nu.evolution_operator_3nu(\n'
             '    H, L1 + L2))\n'
             'U_1 = np.array(oscprob3nu.evolution_operator_3nu(H, L1))\n'
             'U_2 = np.array(oscprob3nu.evolution_operator_3nu(H, L2))\n\n'
             'print("max |U(L1+L2) - U(L2) U(L1)| = %.2e"\n'
             '      % np.max(np.abs(U_2 @ U_1 - U_whole)))\n\n'
             '# Both segments share one Hamiltonian, so their operators\n'
             '# commute and the order happens not to matter *here*:\n'
             'print("with the factors swapped     = %.2e"\n'
             '      % np.max(np.abs(U_1 @ U_2 - U_whole)))'),
        md('That second line is not the general case, and it is worth being '
           'clear about why it came out the same. Both segments carry the '
           '*same* Hamiltonian, so their operators are functions of one '
           'matrix and commute. The moment the segments differ — which is '
           'what varying matter means — the order matters, and the next '
           'section shows it costing you a visibly different answer.'),
        md('`slabs.evolution_operator_3nu_slabs` does that multiplication '
           'for you, in order, given one Hamiltonian per slab:'),
        code('import earth\n'
             'import slabs\n\n'
             '# Deliberately asymmetric: the order is then observable\n'
             'widths = np.array([300.0, 600.0, 500.0])*KM\n'
             'densities = np.array([2.6, 4.5, 9.0])          # g/cm^3\n\n'
             'H_seg = hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '    H_VAC_3NU, energy, earth.matter_potential(densities))\n\n'
             'U_slabs = np.array(slabs.evolution_operator_3nu_slabs(\n'
             '    H_seg, widths))\n\n'
             '# The same product, written out\n'
             'U_manual = np.eye(3, dtype=complex)\n'
             'for k in range(len(widths)):\n'
             '    U_manual = np.array(oscprob3nu.evolution_operator_3nu(\n'
             '        H_seg[k], widths[k])) @ U_manual\n\n'
             'print("max |slabs - by hand| = %.2e"\n'
             '      % np.max(np.abs(U_slabs - U_manual)))\n\n'
             '# Now traverse the same three slabs in reverse\n'
             'U_rev = np.array(slabs.evolution_operator_3nu_slabs(\n'
             '    H_seg[::-1], widths[::-1]))\n'
             'print("max |forward - reversed| = %.2e" \n'
             '      % np.max(np.abs(U_slabs - U_rev)))\n'
             'print("P_mue forward  = %.6f" % abs(U_slabs[0, 1])**2)\n'
             'print("P_mue reversed = %.6f" % abs(U_rev[0, 1])**2)'),
        md('Same three slabs, same total length, same mean density — and a '
           'different answer, because matrix multiplication does not '
           'commute. A neutrino that meets dense matter last is not in the '
           'same state as one that met it first. Notebook 08 turns that '
           'observation into physics.'),
        md('## Inside the closed form\n\n'
           'The expansion writes the Hamiltonian in the basis of SU(3) '
           'generators (the Gell-Mann matrices $\\lambda^k$),\n\n'
           '$$ H = h_0\\mathbb{1} + h_k\\lambda^k , $$\n\n'
           'and the evolution operator in the same basis,\n\n'
           '$$ U_3(L) = u_0\\mathbb{1} + i\\,u_k\\lambda^k . $$\n\n'
           'Each of those steps is a public routine, so the method can be '
           'inspected rather than assumed.'),
        code('h = oscprob3nu.hamiltonian_3nu_coefficients(H)\n'
             'print("h_k : %d real coefficients" % len(h))\n'
             'print("     ", np.array2string(np.array(h), precision=3))\n\n'
             '# Two invariants control the whole result\n'
             'h2, h3 = oscprob3nu.su3_invariants(h)\n'
             'print("\\n|h|^2      = %.6e   (= Tr(H~^2)/2)" % h2)\n'
             'print("<h>        = %.6e   (= Tr(H~^3)/2)" % h3)'),
        md('The oscillation frequencies are the roots of the characteristic '
           'equation, which for SU(3) is a cubic and has a closed '
           'trigonometric solution. Up to a sign they are the eigenvalues '
           'of the traceless Hamiltonian — the quantities every oscillation '
           'phase is built from. The paper writes the characteristic '
           'equation for $-\\tilde H$, so $\\psi_m = -\\lambda_m$:'),
        code('psi = np.sort(np.array(oscprob3nu.psi_roots(h2, h3)))\n\n'
             '# The roots are those of the characteristic equation of -H~,\n'
             '# so they are minus the eigenvalues of H~ itself.\n'
             'H_traceless = H - np.trace(H)/3.0*np.eye(3)\n'
             'reference = np.sort(-np.linalg.eigvalsh(H_traceless))\n\n'
             'print("psi_m (closed form) :", np.array2string(psi, '
             'precision=4))\n'
             'print("-eigvalsh(H~)       :", np.array2string(reference, '
             'precision=4))\n'
             'print("max difference      : %.2e"\n'
             '      % np.max(np.abs(psi - reference)))'),
        md('And the coefficients of the operator, which rebuild it '
           'exactly. The Gell-Mann matrices are written out by hand '
           'below, and the reason is worth knowing if you go looking for '
           'them: `oscprob4nu` exposes `generators_su4()` because at four '
           'flavors every quantity it needs is a trace against the '
           'generators, whereas `oscprob3nu` contracts the $d$ tensor '
           'instead — see `oscprob3nu.tensor_d` — and so never builds the '
           'matrices. There is no `generators_su3()` to import.'),
        code('u = oscprob3nu.evolution_operator_3nu_u_coefficients(\n'
             '    H, baseline)\n'
             'print("u_0 and u_1..u_8 : %d complex coefficients" % len(u))\n\n'
             '# Rebuild U from them, using the Gell-Mann matrices\n'
             'lam = np.zeros((8, 3, 3), dtype=complex)\n'
             'lam[0, 0, 1] = lam[0, 1, 0] = 1.0\n'
             'lam[1, 0, 1], lam[1, 1, 0] = -1.0j, 1.0j\n'
             'lam[2, 0, 0], lam[2, 1, 1] = 1.0, -1.0\n'
             'lam[3, 0, 2] = lam[3, 2, 0] = 1.0\n'
             'lam[4, 0, 2], lam[4, 2, 0] = -1.0j, 1.0j\n'
             'lam[5, 1, 2] = lam[5, 2, 1] = 1.0\n'
             'lam[6, 1, 2], lam[6, 2, 1] = -1.0j, 1.0j\n'
             'lam[7] = np.diag([1.0, 1.0, -2.0])/np.sqrt(3.0)\n\n'
             'U_rebuilt = (u[0]*np.eye(3)\n'
             '             + 1.0j*np.einsum("k,kij->ij",\n'
             '                              np.array(u[1:]), lam))\n'
             'print("max |rebuilt - U| = %.2e"\n'
             '      % np.max(np.abs(U_rebuilt - U)))'),
        md('Seeing them vary is more instructive than one value of each. '
           'The $h_k$ are fixed by the Hamiltonian and do not move with '
           'distance; the $u_k$ carry the whole $L$ dependence, '
           'oscillating at the three frequencies $\\psi_m$ set by the '
           'latent roots. Every probability in every other notebook is '
           'built from these:'),
        code('L_km = np.linspace(0.0, 6000.0, 800)\n'
             'u_of_L = np.array([\n'
             '    oscprob3nu.evolution_operator_3nu_u_coefficients(\n'
             '        H, l*KM) for l in L_km])\n\n'
             'fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.6),\n'
             '                         sharex=True)\n'
             'axes[0].plot(L_km, u_of_L[:, 0].real, label=r"$u_0$")\n'
             'for k in (1, 2, 3):\n'
             '    axes[0].plot(L_km, u_of_L[:, k].real, lw=1.0,\n'
             '                 label=r"$u_%d$" % k)\n'
             'axes[0].set_ylabel("coefficient (real part)")\n'
             'axes[0].set_title("The SU(3) coefficients of $U_3(L)$")\n'
             '# Headroom so the legend clears the u_0 peak\n'
             'axes[0].set_ylim(-0.75, 1.55)\n'
             'axes[0].legend(ncol=4, fontsize=9, loc="upper center")\n\n'
             '# Rebuilt probabilities, for comparison\n'
             'p_scan = oscprob3nu.probabilities_3nu(H, L_km*KM)\n'
             'axes[1].plot(L_km, p_scan[:, 0], label=r"$P_{ee}$")\n'
             'axes[1].plot(L_km, p_scan[:, 1], label=r"$P_{e\\mu}$")\n'
             'axes[1].set_xlim(L_km[0], L_km[-1])\n'
             'axes[1].set_ylim(0.0, 1.0)\n'
             'axes[1].set_xlabel("Baseline [km]")\n'
             'axes[1].set_ylabel("probability")\n'
             'axes[1].legend(ncol=2)\n'
             'fig.tight_layout()\n'
             'plt.show()'),
        md('The lower panel is what the other notebooks plot; the upper '
           'panel is what it is made of. $u_0$ is the coefficient of the '
           'identity, so it stays near one where little has happened, and '
           'the $u_k$ grow as flavor content moves between states.'),
        md('## The same shape at two and four flavors\n\n'
           'Nothing above is specific to three. Each module carries the same '
           'four entry points — coefficients, invariants, roots, operator — '
           'over SU(2), SU(3) and SU(4).'),
        code('# Two flavors\n'
             'H2 = np.asarray(\n'
             '    hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(\n'
             '        gd.S23_NO_BF, gd.D31_NO_BF))/energy\n'
             'U2 = np.array(oscprob2nu.evolution_operator_2nu(H2, baseline))\n'
             'print("SU(2): U is %s, unitary to %.1e"\n'
             '      % (U2.shape,\n'
             '         np.max(np.abs(U2.conj().T @ U2 - np.eye(2)))))\n\n'
             '# Four flavors: fifteen generators, three invariants\n'
             'lam4 = oscprob4nu.generators_su4()\n'
             'H4 = np.diag([1.0, 2.0, 3.0, -6.0]).astype(complex)\n'
             'i2, i3, i4 = oscprob4nu.su4_invariants(H4)\n'
             'U4 = np.array(oscprob4nu.evolution_operator_4nu(H4, 1.0))\n'
             'print("SU(4): %d generators, invariants "\n'
             '      "%.3f %.3f %.3f" % (len(lam4), i2, i3, i4))\n'
             'print("       U is %s, unitary to %.1e"\n'
             '      % (U4.shape,\n'
             '         np.max(np.abs(U4.conj().T @ U4 - np.eye(4)))))'),
        md('## Which one to reach for\n\n'
           '| You want | Use |\n'
           '|---|---|\n'
           '| A probability | `probabilities_2nu` / `_3nu` / `_4nu` |\n'
           '| To compose across segments, or a density matrix | '
           '`evolution_operator_*` |\n'
           '| A trajectory through layered matter | `slabs.*`, `earth.*` |\n'
           '| The oscillation frequencies alone | `psi_roots`, '
           '`psi_roots_4nu` |\n'
           '| To inspect or extend the expansion | `*_coefficients`, '
           '`su3_invariants`, `su4_invariants` |\n\n'
           'The probability routines are the ones to use by default: they '
           'take the shortest path and avoid building the operator when it '
           'is not needed. Reach past them when you need what they discard.'),
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
    ('19_animations.ipynb', 0): 'gallery_anim_cp.png',
    ('19_animations.ipynb', 1): 'gallery_anim_sterile.png',
    ('19_animations.ipynb', 2): 'gallery_anim_earth.png',
    ('19_animations.ipynb', 3): 'gallery_anim_slabs.png',
    ('20_arbitrary_hamiltonian.ipynb', 5): 'gallery_long_range.png',
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


# ------------------------------------------------------------ 19 animations
books['19_animations.ipynb'] = notebook(
    'Animated scenes',
    'Four short scenes, each showing one thing the library does while a '
    'parameter sweeps. They began as the scripts behind a demonstration '
    'reel; this notebook is the same four, with the sweeps drawn as static '
    'filmstrips so that reading it costs nothing.\n\n'
    'Rendering them as animations is expensive — a few hundred frames each, '
    'most of the time spent in matplotlib rather than in the physics — so '
    'the stills are what this notebook draws, and the animation is left to '
    'an opt-in switch at the end.\n\n'
    'The four are chosen to differ, and the difference is the point: the '
    'first two recompute a whole map in one broadcast call, the third '
    're-slabs the Earth at every angle, and the fourth animates the one '
    'approximation the library makes.',
    [
        code('from matplotlib.patches import Circle\n\n'
             'import earth\n'
             'import slabs\n\n'
             '# Indices into the flavor-ordered probability tuples: the '
             'initial\n'
             '# flavor varies slowest, so with n flavors P[n*a + b] is '
             'P(nu_a -> nu_b).\n'
             'P3_MUE, P3_MUMU, P4_MUS = 3, 4, 7\n\n'
             "ACCENT, MARK, MUTED = '#1d4ed8', '#dc2626', '#94a3b8'\n\n\n"
             'def vacuum_3nu(dcp=gd.DCP_NO_BF):\n'
             '    return '
             'hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(\n'
             '        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, dcp,\n'
             '        gd.D21_NO_BF, gd.D31_NO_BF)\n\n\n'
             'def vacuum_4nu(d41):\n'
             '    return '
             'hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(\n'
             '        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF,\n'
             '        np.sqrt(0.10), np.sqrt(0.10), 0.0,\n'
             '        gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF, d41)'),

        md('## 1. The CP phase\n\n'
           'An oscillogram of $P_{\\mu e}$ in matter, recomputed at each '
           'phase, beside the bi-probability ellipse — which is the locus '
           'traced *by* the phase, so the ellipse is drawn once and the '
           'marker says where on it the map currently sits.\n\n'
           'Each map is 200 x 200 = 40 000 probabilities in **one** call: '
           'the Hamiltonians vary along one axis and the baselines along '
           'the other, and they broadcast.'),
        code('grid = 200\n'
             'energies_cp = np.logspace(-1.0, 1.0, grid)*GEV\n'
             'baselines_cp = np.linspace(50.0, 12000.0, grid)*KM\n\n\n'
             'def oscillogram(dcp):\n'
             '    stack = hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '        vacuum_3nu(dcp), energies_cp, gd.VCC_EARTH_CRUST)\n'
             '    return oscprob3nu.probabilities_3nu(\n'
             '        stack[:, None, :, :], '
             'baselines_cp[None, :])[:, :, P3_MUE]\n\n\n'
             'def ellipse_point(dcp, energy=0.8*GEV, baseline=1300.0):\n'
             '    length = baseline*KM\n'
             '    # Antineutrinos need *both* changes: conjugate the '
             'vacuum\n'
             '    # Hamiltonian and reverse the sign of the matter '
             'potential.\n'
             '    h_nu = hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '        vacuum_3nu(dcp), energy, gd.VCC_EARTH_CRUST)\n'
             '    h_bar = hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '        np.conj(vacuum_3nu(dcp)), energy, '
             '-gd.VCC_EARTH_CRUST)\n'
             '    return (oscprob3nu.probabilities_3nu(h_nu, '
             'length)[P3_MUE],\n'
             '            oscprob3nu.probabilities_3nu(h_bar, '
             'length)[P3_MUE])\n\n\n'
             'shown = [0.0, 2.0*np.pi/3.0, 4.0*np.pi/3.0]\n'
             'locus = np.array([ellipse_point(d)\n'
             '                  for d in np.linspace(0.0, 2.0*np.pi, '
             '240)])\n'
             '# The colour scale is fixed over the whole sweep: taking it '
             'from\n'
             '# one frame lets the others clip silently, and clipping '
             'reads as\n'
             '# structure rather than as saturation.\n'
             'ceiling_cp = max(float(oscillogram(d).max())\n'
             '                 for d in np.linspace(0.0, 2.0*np.pi, 12,\n'
             '                                      endpoint=False))\n\n'
             'fig, axes = plt.subplots(1, 4, figsize=(15.0, 3.4),\n'
             "                         gridspec_kw={'width_ratios': "
             '[1, 1, 1, 1.05]})\n'
             'for ax, dcp in zip(axes[:3], shown):\n'
             "    im = ax.imshow(oscillogram(dcp), origin='lower', "
             "aspect='auto',\n"
             "                   cmap='viridis', vmin=0.0, "
             "vmax=ceiling_cp,\n"
             '                   extent=[50.0, 12000.0, -1.0, 1.0])\n'
             "    ax.set_title(r'$\\delta_{CP} = %.2f\\pi$' % "
             '(dcp/np.pi), fontsize=11)\n'
             "    ax.set_xlabel('Baseline [km]')\n"
             '    ax.set_yticks([-1, 0, 1])\n'
             "    ax.set_yticklabels(['0.1', '1', '10'])\n"
             "    ax.grid(False)\n"
             "axes[0].set_ylabel('Energy [GeV]')\n"
             "fig.colorbar(im, ax=axes[2], "
             "pad=0.02).set_label(r'$P_{\\mu e}$')\n\n"
             'ax = axes[3]\n'
             'ax.plot(locus[:, 0], locus[:, 1], color=ACCENT, lw=1.6)\n'
             'lo, hi = float(min(locus.min(), 0.0)), '
             'float(locus.max())*1.10\n'
             "ax.plot([lo, hi], [lo, hi], ls=':', color='#64748b', "
             "lw=1.0,\n"
             "        label='CP symmetric')\n"
             'for dcp in shown:\n'
             "    ax.plot(*ellipse_point(dcp), 'o', ms=8, color=MARK, "
             'zorder=5)\n'
             'ax.set_xlim(lo, hi)\n'
             'ax.set_ylim(lo, hi)\n'
             "ax.set_aspect('equal', adjustable='box')\n"
             "ax.set_xlabel(r'$P(\\nu_\\mu \\to \\nu_e)$')\n"
             "ax.set_ylabel(r'$P(\\bar\\nu_\\mu \\to "
             "\\bar\\nu_e)$')\n"
             "ax.set_title('Bi-probability, 0.8 GeV, 1300 km', "
             'fontsize=11)\n'
             "ax.legend(loc='upper left', fontsize=9)\n\n"
             "fig.suptitle('The CP phase sweeping through $2\\pi$', "
             'fontsize=13)\n'
             'fig.tight_layout(rect=[0, 0, 1, 0.93])\n'
             'plt.show()'),

        md('## 2. A sterile state\n\n'
           'Four flavors, in matter of constant density so the whole map '
           'is again one call. The sterile state feels neither the '
           'charged- nor the neutral-current potential, so $V_{NC}$ stops '
           'cancelling between the flavors and sits on the sterile entry '
           '— which is what places the resonance that moves across the '
           'frame as $\\Delta m^2_{41}$ sweeps.'),
        code('grid = 180\n'
             'energies_s = np.logspace(0.0, 2.0, grid)*GEV\n'
             'baselines_s = np.linspace(100.0, 12742.0, grid)*KM\n\n\n'
             'def sterile_map(d41):\n'
             '    stack = hamiltonians4nu.hamiltonian_4nu_matter(\n'
             '        vacuum_4nu(d41), energies_s, gd.VCC_EARTH_CRUST,\n'
             '        gd.VNC_EARTH_CRUST)\n'
             '    return oscprob4nu.probabilities_4nu(\n'
             '        stack[:, None, :, :], '
             'baselines_s[None, :])[:, :, P4_MUS]\n\n\n'
             'splittings = [0.1, 1.0, 10.0]\n'
             'maps = [sterile_map(d) for d in splittings]\n'
             'ceiling_s = max(float(m.max()) for m in maps)\n\n'
             'fig, axes = plt.subplots(1, 4, figsize=(15.0, 3.4),\n'
             "                         gridspec_kw={'width_ratios': "
             '[1, 1, 1, 1.05]})\n'
             'for ax, d41, m in zip(axes[:3], splittings, maps):\n'
             "    im = ax.imshow(m, origin='lower', aspect='auto', "
             "cmap='magma',\n"
             '                   vmin=0.0, vmax=ceiling_s,\n'
             '                   extent=[100.0, 12742.0, 0.0, 2.0])\n'
             "    ax.set_title(r'$\\Delta m^2_{41} = %g$ eV$^2$' % d41, "
             'fontsize=11)\n'
             "    ax.set_xlabel('Baseline [km]')\n"
             '    ax.set_yticks([0, 1, 2])\n'
             "    ax.set_yticklabels(['1', '10', '100'])\n"
             "    ax.axvline(12742.0, color='white', lw=1.0, ls='--', "
             'alpha=0.7)\n'
             "    ax.grid(False)\n"
             "axes[0].set_ylabel('Energy [GeV]')\n"
             'fig.colorbar(im, ax=axes[2], pad=0.02).set_label(\n'
             "    r'$P(\\nu_\\mu\\to\\nu_s)$')\n\n"
             'ax = axes[3]\n'
             'for d41, m in zip(splittings, maps):\n'
             '    ax.plot(np.log10(energies_s/GEV), m[:, -1],\n'
             "            lw=1.8, label=r'$%g$ eV$^2$' % d41)\n"
             'ax.set_xlim(0.0, 2.0)\n'
             'ax.set_ylim(0.0, ceiling_s*1.05)\n'
             'ax.set_xticks([0, 1, 2])\n'
             "ax.set_xticklabels(['1', '10', '100'])\n"
             "ax.set_xlabel('Energy [GeV]')\n"
             "ax.set_ylabel(r'$P(\\nu_\\mu \\to \\nu_s)$')\n"
             "ax.set_title('Through the diameter, 12742 km', "
             'fontsize=11)\n'
             'ax.legend(fontsize=9)\n\n'
             "fig.suptitle('A sterile state appearing as "
             "$\\Delta m^2_{41}$ sweeps',\n"
             '             fontsize=13)\n'
             'fig.tight_layout(rect=[0, 0, 1, 0.93])\n'
             'plt.show()'),

        md('## 3. Through the Earth\n\n'
           'The PREM density in cross-section with the neutrino path '
           'drawn on it, beside the survival probability along that path, '
           'as the zenith angle swings from diametric to grazing.\n\n'
           '**This scene used to be the expensive one.** When these '
           'scripts were written `probabilities_3nu_earth` took a scalar '
           'energy — the chord, and the slabs cut from it, change with '
           'the angle — so a curve of a hundred energies was a hundred '
           'calls, and the scene said so in its own caption. Since 1.12.0 '
           'the energy may be an array, so each curve below is **one** '
           'call, with the geometry and the matter potentials built once '
           'for it rather than once per energy.'),
        code('angles = [-1.0, -0.6, -0.2]\n'
             'energies_e = np.logspace(0.0, 2.0, 100)*GEV\n\n'
             'fig, axes = plt.subplots(1, 4, figsize=(15.0, 3.6),\n'
             "                         gridspec_kw={'width_ratios': "
             '[1, 1, 1, 1.35]})\n\n'
             '# The density disk is the same in each panel: concentric '
             'circles\n'
             '# shaded by PREM, outermost first so the inner ones sit on '
             'top.\n'
             'radii = np.linspace(gd.EARTH_RADIUS, 0.0, 240)\n'
             'densities = earth.density_prem(radii)\n'
             'colours = plt.cm.YlOrBr(0.15 + '
             '0.75*densities/densities.max())\n\n'
             'for ax, costhz in zip(axes[:3], angles):\n'
             '    for radius, colour in zip(radii, colours):\n'
             '        ax.add_patch(Circle((0.0, 0.0), radius, '
             'facecolor=colour,\n'
             "                            edgecolor='none', zorder=1))\n"
             '    ax.add_patch(Circle((0.0, 0.0), gd.EARTH_RADIUS, '
             'fill=False,\n'
             "                        edgecolor='#334155', lw=1.2, "
             'zorder=3))\n'
             '    length = earth.distance_traveled_inside_earth(costhz)\n'
             '    sin_thz = np.sqrt(max(0.0, 1.0 - costhz*costhz))\n'
             '    ax.plot([length*sin_thz, 0.0],\n'
             '            [gd.EARTH_RADIUS + length*costhz, '
             'gd.EARTH_RADIUS],\n'
             '            color=MARK, lw=2.4, zorder=4)\n'
             "    ax.plot([0.0], [gd.EARTH_RADIUS], 'v', "
             "color='#0f172a', ms=9,\n"
             '            zorder=5)\n'
             '    span = gd.EARTH_RADIUS*1.08\n'
             '    ax.set_xlim(-span, span)\n'
             '    ax.set_ylim(-span, span)\n'
             "    ax.set_aspect('equal')\n"
             "    ax.axis('off')\n"
             "    ax.set_title(r'$\\cos\\theta_z = %+.2f$,  %d km'\n"
             '                 % (costhz, round(length)), fontsize=11)\n\n'
             'ax = axes[3]\n'
             'for costhz in angles:\n'
             '    # One call per curve: the energies are an array.\n'
             '    curve = earth.probabilities_3nu_earth(\n'
             '        H_VAC_3NU, energies_e, costhz)[:, P3_MUMU]\n'
             '    ax.plot(np.log10(energies_e/GEV), curve, lw=1.8,\n'
             "            label=r'$\\cos\\theta_z = %+.2f$' % "
             'costhz)\n'
             'ax.set_xlim(0.0, 2.0)\n'
             'ax.set_ylim(0.0, 1.02)\n'
             'ax.set_xticks([0, 1, 2])\n'
             "ax.set_xticklabels(['1', '10', '100'])\n"
             "ax.set_xlabel('Energy [GeV]')\n"
             "ax.set_ylabel(r'$P(\\nu_\\mu \\to \\nu_\\mu)$')\n"
             "ax.set_title('Survival along the chord', fontsize=11)\n"
             'ax.legend(fontsize=9)\n\n'
             "fig.suptitle('A chord through the Earth, as the zenith "
             "angle swings',\n"
             '             fontsize=13)\n'
             'fig.tight_layout(rect=[0, 0, 1, 0.92])\n'
             'plt.show()'),

        md('## 4. Cutting a profile into slabs\n\n'
           'The one approximation the library makes. Within a slab '
           'nothing is approximated — the expansion is exact for a '
           'constant Hamiltonian — so the only question is how finely a '
           'profile that really varies is sliced, and that is the '
           'caller\'s to answer. The reference here is the same '
           'calculation at 600 slabs.'),
        code('baseline_km, energy_s = 4000.0, 1.0*GEV\n'
             'total = baseline_km*KM\n\n\n'
             'def profile(x):\n'
             '    """A smooth, deliberately awkward density profile '
             '[g/cm^3].\n\n'
             '    The amplitudes sum to more than the offset, so the '
             'positivity\n'
             '    of this is a property of the phase offset and the '
             'frequency\n'
             '    ratio rather than of the coefficients: the two sines do '
             'not\n'
             '    reach their minima together, and the measured minimum '
             'is\n'
             '    0.28 g/cm^3.  A negative density would quietly reverse '
             'the\n'
             '    sign of the matter potential --- the antineutrino case, '
             'and\n'
             '    not what this is about --- so re-check the minimum '
             'before\n'
             '    changing any of the five numbers below.\n'
             '    """\n'
             '    return 3.0 + 2.2*np.sin(2.0*np.pi*x/baseline_km) \\\n'
             '        + 1.1*np.sin(6.0*np.pi*x/baseline_km + 0.7)\n\n\n'
             'def probability(n):\n'
             '    edges = np.linspace(0.0, baseline_km, n + 1)\n'
             '    mid = 0.5*(edges[:-1] + edges[1:])\n'
             '    h = hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '        vacuum_3nu(), energy_s, '
             'earth.matter_potential(profile(mid)))\n'
             '    return slabs.probabilities_3nu_slabs(\n'
             '        h, np.full(n, total/n))[P3_MUE]\n\n\n'
             'counts = np.unique(np.round(\n'
             '    np.logspace(0.0, np.log10(80.0), 40)).astype(int))\n'
             'values = [probability(int(n)) for n in counts]\n'
             'reference = probability(600)\n'
             '# Snapped to counts that were actually evaluated, so the '
             'markers\n'
             '# on the trace are the same three cuts the panels draw.\n'
             'shown_n = [int(counts[np.argmin(np.abs(counts - n))])\n'
             '           for n in (2, 8, 40)]\n\n'
             'fig, axes = plt.subplots(1, 4, figsize=(15.0, 3.4),\n'
             "                         gridspec_kw={'width_ratios': "
             '[1, 1, 1, 1.1]})\n'
             'fine = np.linspace(0.0, baseline_km, 600)\n'
             'for ax, n in zip(axes[:3], shown_n):\n'
             '    edges = np.linspace(0.0, baseline_km, n + 1)\n'
             '    mid = 0.5*(edges[:-1] + edges[1:])\n'
             '    ax.plot(fine, profile(fine), color=MUTED, lw=1.6,\n'
             "            label='the true profile')\n"
             '    ax.step(np.append(edges[:-1], baseline_km),\n'
             '            np.append(profile(mid), profile(mid[-1])),\n'
             "            where='post', color=ACCENT, lw=1.8, "
             "label='what is solved')\n"
             "    ax.set_title('%d slab%s' % (n, '' if n == 1 else 's'), "
             'fontsize=11)\n'
             "    ax.set_xlabel('Distance [km]')\n"
             "axes[0].set_ylabel(r'Density [g cm$^{-3}$]')\n"
             "axes[0].legend(loc='upper right', fontsize=8)\n\n"
             'ax = axes[3]\n'
             "ax.axhline(reference, color=MUTED, ls='--', lw=1.2, "
             "label='600 slabs')\n"
             'ax.plot(counts, values, color=ACCENT, lw=1.6)\n'
             'for n in shown_n:\n'
             '    index = int(np.argmin(np.abs(counts - n)))\n'
             "    ax.plot([counts[index]], [values[index]], 'o', ms=7, "
             'color=MARK,\n'
             '            zorder=5)\n'
             "ax.set_xscale('log')\n"
             'ax.set_xlim(1, 90)\n'
             "ax.set_xlabel('Number of slabs')\n"
             "ax.set_ylabel(r'$P_{\\mu e}$')\n"
             "ax.set_title('Converging on the answer', fontsize=11)\n"
             "ax.legend(loc='lower right', fontsize=9)\n\n"
             "fig.suptitle('Each slab is solved exactly; only the cutting "
             "is approximate',\n"
             '             fontsize=13)\n'
             'fig.tight_layout(rect=[0, 0, 1, 0.93])\n'
             'plt.show()'),

        md('## Rendering them as animations\n\n'
           'Each scene is a sweep over one parameter, so each animates '
           'directly: the panels below are redrawn frame by frame while '
           '$\\delta_{CP}$, $\\Delta m^2_{41}$, $\\cos\\theta_z$ '
           'or the slab count moves.\n\n'
           'The stills above are what this notebook draws by default, '
           'because rendering is expensive — a few hundred frames each, '
           'and most of the time goes into matplotlib rather than into '
           'the physics. Set `RENDER = True` to write all four as GIFs. '
           'Expect several minutes and a few megabytes each. '
           '`matplotlib`\'s Pillow writer is used rather than `ffmpeg`, '
           'so there is nothing to install — though the committed GIFs in '
           '`img/` went through an ffmpeg palette, which is markedly '
           'smaller.\n\n'
           'The scenes are drawn two panels wide here, as the reel has '
           'them, rather than as the four-panel filmstrips above.'),
        code('RENDER = False        # set True to write the four GIFs\n'
             'FPS = 24\n\n'
             'if RENDER:\n'
             '    from matplotlib.animation import FuncAnimation, '
             'PillowWriter\n\n'
             '    FIGSIZE, DPI = (11.0, 4.6), 120\n'
             "    OUT = os.path.join('..', 'img')\n\n"
             '    def caption(fig, note):\n'
             '        """The standing caption, and the mutable '
             'heading."""\n'
             "        fig.patch.set_facecolor('white')\n"
             '        fig.text(0.5, 0.028, note, ha=\'center\', '
             "fontsize=9.5,\n"
             "                 color='#334155')\n"
             "        heading = fig.suptitle('', fontsize=13)\n"
             '        fig.tight_layout(rect=[0, 0.06, 1, 0.94])\n'
             '        return heading\n\n'
             '    def write(fig, update, frames, name):\n'
             '        anim = FuncAnimation(fig, update, frames=frames, '
             'blit=False)\n'
             '        path = os.path.join(OUT, name)\n'
             "        anim.save(path, writer=PillowWriter(fps=FPS),\n"
             "                  savefig_kwargs={'facecolor': 'white'})\n"
             "        plt.close(fig)\n"
             "        print('wrote %s (%.1f MB)'\n"
             "              % (path, os.path.getsize(path)/1024.0/1024.0))\n\n"
             '    # ---- the CP phase '
             '-------------------------------------------\n'
             '    phases = np.linspace(0.0, 2.0*np.pi, 240, '
             'endpoint=False)\n'
             '    fig, (ax_map, ax_ell) = plt.subplots(\n'
             '        1, 2, figsize=FIGSIZE, dpi=DPI,\n'
             "        gridspec_kw={'width_ratios': [1.35, 1.0]})\n"
             '    image = ax_map.imshow(oscillogram(phases[0]), '
             "origin='lower',\n"
             "                          aspect='auto', cmap='viridis',\n"
             '                          vmin=0.0, vmax=ceiling_cp,\n'
             '                          extent=[50.0, 12000.0, -1.0, '
             '1.0])\n'
             "    ax_map.set_xlabel('Baseline [km]')\n"
             "    ax_map.set_ylabel('Energy [GeV]')\n"
             '    ax_map.set_yticks([-1, 0, 1])\n'
             "    ax_map.set_yticklabels(['0.1', '1', '10'])\n"
             '    ax_map.grid(False)\n'
             '    ax_ell.plot(locus[:, 0], locus[:, 1], color=ACCENT, '
             'lw=1.6)\n'
             "    marker, = ax_ell.plot([], [], 'o', ms=9, color=MARK, "
             'zorder=5)\n'
             '    ax_ell.set_xlim(lo, hi)\n'
             '    ax_ell.set_ylim(lo, hi)\n'
             "    ax_ell.set_aspect('equal', adjustable='box')\n"
             "    ax_ell.set_xlabel(r'$P(\\nu_\\mu \\to "
             "\\nu_e)$')\n"
             "    ax_ell.set_ylabel(r'$P(\\bar\\nu_\\mu \\to "
             "\\bar\\nu_e)$')\n"
             "    head = caption(fig, 'One broadcast call a frame  ---  "
             "40,000 probabilities')\n\n"
             '    def update_cp(i):\n'
             '        image.set_data(oscillogram(phases[i]))\n'
             '        marker.set_data(*[[v] for v in '
             'ellipse_point(phases[i])])\n'
             "        head.set_text(r'The CP phase:  $\\delta_{CP} = "
             "%.2f\\pi$'\n"
             '                      % (phases[i]/np.pi))\n'
             '        return image, marker\n\n'
             "    write(fig, update_cp, len(phases), 'anim_cp.gif')\n\n"
             '    # ---- the sterile state '
             '---------------------------------------\n'
             '    d41s = np.logspace(-1.0, 1.0, 180)\n'
             '    fig, (ax_map, ax_cut) = plt.subplots(\n'
             '        1, 2, figsize=FIGSIZE, dpi=DPI,\n'
             "        gridspec_kw={'width_ratios': [1.35, 1.0]})\n"
             '    image = ax_map.imshow(sterile_map(d41s[0]), '
             "origin='lower',\n"
             "                          aspect='auto', cmap='magma',\n"
             '                          vmin=0.0, vmax=ceiling_s,\n'
             '                          extent=[100.0, 12742.0, 0.0, '
             '2.0])\n'
             "    ax_map.set_xlabel('Baseline [km]')\n"
             "    ax_map.set_ylabel('Energy [GeV]')\n"
             '    ax_map.set_yticks([0, 1, 2])\n'
             "    ax_map.set_yticklabels(['1', '10', '100'])\n"
             '    ax_map.grid(False)\n'
             '    curve, = ax_cut.plot([], [], color=ACCENT, lw=1.8)\n'
             '    ax_cut.set_xlim(0.0, 2.0)\n'
             '    ax_cut.set_ylim(0.0, ceiling_s*1.05)\n'
             '    ax_cut.set_xticks([0, 1, 2])\n'
             "    ax_cut.set_xticklabels(['1', '10', '100'])\n"
             "    ax_cut.set_xlabel('Energy [GeV]')\n"
             "    ax_cut.set_ylabel(r'$P(\\nu_\\mu \\to "
             "\\nu_s)$')\n"
             "    ax_cut.set_title('Through the diameter, 12742 km', "
             'fontsize=11)\n'
             "    head = caption(fig, 'Four flavors, constant density  "
             "---  the sterile '\n"
             "                        'state feels neither potential')\n"
             '    axis_s = np.log10(energies_s/GEV)\n\n'
             '    def update_sterile(i):\n'
             '        now = sterile_map(d41s[i])\n'
             '        image.set_data(now)\n'
             '        curve.set_data(axis_s, now[:, -1])\n'
             "        head.set_text(r'A sterile state:  $\\Delta "
             "m^2_{41} = %.2f$ eV$^2$'\n"
             '                      % d41s[i])\n'
             '        return image, curve\n\n'
             "    write(fig, update_sterile, len(d41s), "
             "'anim_sterile.gif')\n\n"
             '    # ---- through the Earth '
             '---------------------------------------\n'
             '    cosines = np.linspace(-1.0, -0.05, 150)\n'
             '    fig, (ax_earth, ax_prob) = plt.subplots(\n'
             '        1, 2, figsize=FIGSIZE, dpi=DPI,\n'
             "        gridspec_kw={'width_ratios': [1.0, 1.35]})\n"
             '    for radius, colour in zip(radii, colours):\n'
             '        ax_earth.add_patch(Circle((0.0, 0.0), radius, '
             'facecolor=colour,\n'
             "                                  edgecolor='none', "
             'zorder=1))\n'
             '    ax_earth.add_patch(Circle((0.0, 0.0), gd.EARTH_RADIUS, '
             'fill=False,\n'
             "                              edgecolor='#334155', lw=1.2, "
             'zorder=3))\n'
             '    span = gd.EARTH_RADIUS*1.08\n'
             '    ax_earth.set_xlim(-span, span)\n'
             '    ax_earth.set_ylim(-span, span)\n'
             "    ax_earth.set_aspect('equal')\n"
             "    ax_earth.axis('off')\n"
             "    ax_earth.plot([0.0], [gd.EARTH_RADIUS], 'v', "
             "color='#0f172a',\n"
             '                  ms=9, zorder=5)\n'
             "    path_line, = ax_earth.plot([], [], color=MARK, lw=2.2, "
             'zorder=4)\n'
             '    prob_line, = ax_prob.plot([], [], color=ACCENT, '
             'lw=1.8)\n'
             '    ax_prob.set_xlim(0.0, 2.0)\n'
             '    ax_prob.set_ylim(0.0, 1.02)\n'
             '    ax_prob.set_xticks([0, 1, 2])\n'
             "    ax_prob.set_xticklabels(['1', '10', '100'])\n"
             "    ax_prob.set_xlabel('Energy [GeV]')\n"
             "    ax_prob.set_ylabel(r'$P(\\nu_\\mu \\to "
             "\\nu_\\mu)$')\n"
             "    head = caption(fig, 'The chord is re-cut into PREM slabs "
             "at every '\n"
             "                        'angle  ---  100 energies, one "
             "call')\n"
             '    axis_e = np.log10(energies_e/GEV)\n\n'
             '    def update_earth(i):\n'
             '        costhz = cosines[i]\n'
             '        length = '
             'earth.distance_traveled_inside_earth(costhz)\n'
             '        sin_thz = np.sqrt(max(0.0, 1.0 - costhz*costhz))\n'
             '        path_line.set_data(\n'
             '            [length*sin_thz, 0.0],\n'
             '            [gd.EARTH_RADIUS + length*costhz, '
             'gd.EARTH_RADIUS])\n'
             '        prob_line.set_data(axis_e, '
             'earth.probabilities_3nu_earth(\n'
             '            H_VAC_3NU, energies_e, costhz)[:, P3_MUMU])\n'
             "        head.set_text(r'Through the Earth:  "
             "$\\cos\\theta_z = %+.2f$, '\n"
             "                      r'%d km' % (costhz, round(length)))\n"
             '        return path_line, prob_line\n\n'
             "    write(fig, update_earth, len(cosines), "
             "'anim_earth.gif')\n\n"
             '    # ---- cutting the profile '
             '-------------------------------------\n'
             '    ramp = np.round(np.logspace(0.0, np.log10(80.0), '
             '120)).astype(int)\n'
             '    ramp_values = [probability(int(n)) for n in ramp]\n'
             '    fig, (ax_rho, ax_conv) = plt.subplots(\n'
             '        1, 2, figsize=FIGSIZE, dpi=DPI,\n'
             "        gridspec_kw={'width_ratios': [1.15, 1.0]})\n"
             '    ax_rho.plot(fine, profile(fine), color=MUTED, lw=1.6,\n'
             "                label='the true profile')\n"
             "    stair, = ax_rho.step([], [], where='post', "
             'color=ACCENT, lw=1.8,\n'
             "                         label='what is solved')\n"
             "    ax_rho.set_xlabel('Distance along the trajectory "
             "[km]')\n"
             "    ax_rho.set_ylabel(r'Density [g cm$^{-3}$]')\n"
             '    ax_rho.set_xlim(0.0, baseline_km)\n'
             '    ax_rho.set_ylim(profile(fine).min()-0.3, '
             'profile(fine).max()+0.3)\n'
             "    ax_rho.legend(loc='upper right', fontsize=9)\n"
             "    ax_conv.axhline(reference, color=MUTED, ls='--', "
             "lw=1.2,\n"
             "                    label='600 slabs')\n"
             '    trace, = ax_conv.plot([], [], color=ACCENT, lw=1.6)\n'
             "    dot, = ax_conv.plot([], [], 'o', ms=7, color=MARK, "
             'zorder=5)\n'
             "    ax_conv.set_xscale('log')\n"
             '    ax_conv.set_xlim(1, 90)\n'
             '    edge = 0.12*(max(max(ramp_values), reference)\n'
             '                 - min(min(ramp_values), reference)) or '
             '0.01\n'
             '    ax_conv.set_ylim(min(min(ramp_values), reference)-edge,\n'
             '                     max(max(ramp_values), '
             'reference)+edge)\n'
             "    ax_conv.set_xlabel('Number of slabs')\n"
             "    ax_conv.set_ylabel(r'$P_{\\mu e}$')\n"
             "    ax_conv.legend(loc='lower right', fontsize=9)\n"
             "    head = caption(fig, 'Each slab is solved exactly  ---  "
             "the only '\n"
             "                        'approximation is how finely the "
             "profile is cut')\n\n"
             '    def update_slabs(i):\n'
             '        n = int(ramp[i])\n'
             '        edges = np.linspace(0.0, baseline_km, n + 1)\n'
             '        mid = 0.5*(edges[:-1] + edges[1:])\n'
             '        stair.set_data(np.append(edges[:-1], '
             'baseline_km),\n'
             '                       np.append(profile(mid), '
             'profile(mid[-1])))\n'
             '        trace.set_data(ramp[:i+1], ramp_values[:i+1])\n'
             '        dot.set_data([n], [ramp_values[i]])\n'
             "        head.set_text(r'Layered matter:  %d slab%s,  "
             "$P_{\\mu e} = %.5f$'\n"
             "                      % (n, '' if n == 1 else 's', "
             'ramp_values[i]))\n'
             '        return stair, trace, dot\n\n'
             "    write(fig, update_slabs, len(ramp), 'anim_slabs.gif')\n\n"
             '    # Pillow writes a colour table per frame, so what it '
             'just wrote is\n'
             '    # some twenty times larger than it needs to be.  The '
             'same two-pass\n'
             '    # palette the command line uses brings it down; '
             'skipped without\n'
             '    # ffmpeg, since the GIFs are already usable.\n'
             '    import shutil\n\n'
             "    sys.path.insert(0, os.path.abspath(os.path.join('..', "
             "'tools')))\n"
             '    import make_demo_video\n\n'
             "    if shutil.which('ffmpeg') is None:\n"
             "        print('\\nffmpeg not found: the GIFs above are "
             "Pillow-sized.')\n"
             '    else:\n'
             "        for name in ('cp', 'sterile', 'earth', 'slabs'):\n"
             "            path = os.path.join(OUT, 'anim_%s.gif' % "
             'name)\n'
             '            before = os.path.getsize(path)/1024.0/1024.0\n'
             '            make_demo_video.shrink(path, path, fps=12,\n'
             '                                   width=860, colors=128)\n'
             '            after = os.path.getsize(path)/1024.0/1024.0\n'
             "            print('  %-16s %5.1f MB -> %4.1f MB' % (name, "
             'before, after))\n'
             'else:\n'
             "    print('RENDER is False: the stills above are the "
             "output,')\n"
             "    print('and img/anim_*.gif are the rendered scenes.')"),
        md('### After rendering\n\n'
           'The cell above writes the four GIFs and then shrinks them '
           'in place, through `tools/make_demo_video.py` — the same code '
           'the command line runs, so there is one implementation of the '
           'encoding rather than two.\n\n'
           'To do either step by hand:\n\n'
           '```shell\n'
           '# Join the four scenes into one reel\n'
           'python tools/make_demo_video.py --join img/anim_cp.gif \\\n'
           '    img/anim_sterile.gif img/anim_earth.gif '
           'img/anim_slabs.gif \\\n'
           '    --out ~/reel.mp4\n\n'
           '# Shrink any clip, which is what makes a GIF publishable\n'
           'python tools/make_demo_video.py --shrink ~/reel.mp4 '
           '--out ~/reel.gif\n'
           '```\n\n'
           'Both go through a **shared palette**: one pass computes a '
           'colour table for the whole clip, the second applies it. '
           'Writing a GIF directly gives every frame its own table, '
           'which is where the twenty-fold difference comes from. The '
           'knobs are `--fps`, `--width` and `--colors`, in that order '
           'of effect on size.\n\n'
           '**Two traps, both hit while writing this.** A GIF stores its '
           'frame delays in hundredths of a second, so ffmpeg reads a '
           '90 ms frame as roughly 100 fps and, without an explicit '
           'output rate, writes every frame nine times over — a '
           'four-hundred-frame reel became ninety megabytes and climbing '
           'before it was noticed. And `ffmpeg` installed as a **snap** '
           'has a private `/tmp` and cannot read or write under the real '
           'one; it fails naming a path that plainly exists. Work under '
           '`$HOME`.'),
    ])

# ------------------------------------------ 20 an arbitrary Hamiltonian
books['20_arbitrary_hamiltonian.ipynb'] = notebook(
    'An arbitrary Hamiltonian, through three profiles',
    'The other notebooks vary one thing at a time. `03` puts a '
    'non-standard Hamiltonian through matter of **constant** density; '
    '`08` and `14` put the **standard** Hamiltonian through profiles '
    'that vary; `01` evaluates an arbitrary Hermitian matrix once. '
    'Nothing so far carries a Hamiltonian of your own through a profile '
    'that changes under it.\n\n'
    'This notebook does, three times over — an invented body, an '
    'exponential one, and the Earth.\n\n'
    'The worked case is a **long-range force**: a gauged '
    '$L_e - L_\\mu$ symmetry with a very light mediator, following '
    '[arXiv:1808.02042](https://arxiv.org/abs/1808.02042). It is a good '
    'test of the machinery for one specific reason. $V_{CC}$ reads the '
    'density at the neutrino\'s own position; this potential integrates '
    'the electrons of the *whole body*, so it varies along a trajectory '
    'in a way that has nothing to do with the local density — and it '
    'cannot be expressed in the form the `earth` entry points build.',
    [
        code('import earth\nimport slabs'),
        md('## The Hamiltonian\n\n'
           'Gauge $L_e - L_\\mu$ and the mediator $Z\'$ couples to that '
           'charge. Electrons carry $L_e = 1$, so ordinary matter is a '
           'source; among the neutrinos, $\\nu_e$ carries $+1$, '
           '$\\nu_\\mu$ carries $-1$ and $\\nu_\\tau$ is neutral. So the '
           'new term enters the flavor basis as\n\n'
           '$$ H = \\underbrace{\\frac{H_{\\rm vac}}{E} + V_{CC}\\,'
           '\\mathrm{diag}(1,0,0)}_{\\text{the standard part}} '
           '+ V_{e\\mu}(r)\\,\\mathrm{diag}(1,-1,0) , $$\n\n'
           'with the potential a Yukawa integral over the electrons '
           'around the neutrino,\n\n'
           '$$ V_{e\\mu}(\\mathbf{r}) = \\frac{g\'^2}{4\\pi} \\int d^3r\' '
           '\\, n_e(\\mathbf{r}\')\\, '
           '\\frac{e^{-m_{Z\'}|\\mathbf{r}-\\mathbf{r}\'|}}'
           '{|\\mathbf{r}-\\mathbf{r}\'|} . $$\n\n'
           '**This is the whole point of the notebook.** That matrix is '
           'not $H_{\\rm vac}/E + V_{CC}P_{ee}$ for any $V_{CC}$, so '
           '`earth.probabilities_3nu_earth` — which builds exactly that '
           'form from a zenith angle — cannot produce it, however its '
           'arguments are chosen. The last section shows what to do '
           'instead.\n\n'
           '**Scope.** Only the electrons of the body being crossed are '
           'counted; outside it, the electron density is taken to be '
           'zero. The reference adds the rest of the Universe, which '
           'dominates for a mediator lighter than these, and is a '
           'different calculation.'),
        md('### The integral is one-dimensional\n\n'
           'A $3n$-dimensional integral per slab would be hopeless. But '
           'for a spherically symmetric $n_e$ the angular average of the '
           'Yukawa kernel over a shell of radius $r\'$ is elementary,\n\n'
           '$$ \\left\\langle \\frac{e^{-m d}}{d} \\right\\rangle_{\\rm '
           'shell} = \\frac{e^{-m|r-r\'|} - e^{-m(r+r\')}}{2 m r r\'} , '
           '$$\n\n'
           'and the result separates into an interior and an exterior '
           'piece:\n\n'
           '$$ V_{e\\mu}(r) = \\frac{g\'^2}{r}\\left[ e^{-mr} \\int_0^r '
           'dr\'\\, r\'^2 n_e(r\')\\,\\mathrm{shc}(m r\') '
           '+ r\\,\\mathrm{shc}(m r) \\int_r^R dr\'\\, r\' n_e(r\')\\, '
           'e^{-m r\'} \\right] $$\n\n'
           'with $\\mathrm{shc}(x) \\equiv \\sinh(x)/x$, which is $1$ at '
           'the origin. Two things follow, and both matter here. The '
           '$r$-dependence has left the integrands, so **one pass over '
           'the profile serves every point on the trajectory** — the '
           'potential is built once per body, not once per slab. And '
           'nothing diverges as $m \\to 0$: the bracket becomes the '
           'enclosed-charge form of an ordinary $1/r$ potential.'),
        code('def integral(f, x):\n'
             '    """Trapezoid rule, spelled the same in every NumPy."""\n'
             '    return np.sum(0.5*(f[1:]+f[:-1])*np.diff(x))\n\n\n'
             'def n_e_of_density(density):\n'
             '    """Electron number density [eV^3] from [g/cm^3]."""\n'
             '    # Routed through the library\'s own V_CC so that the two\n'
             '    # potentials cannot disagree about the electron fraction.\n'
             '    return earth.matter_potential(density)/(np.sqrt(2.0)*gd.GF)'
             '\n\n\n'
             'def shc(x):\n'
             '    """sinh(x)/x, continued to 1 at the origin."""\n'
             '    safe = np.where(x == 0.0, 1.0, x)\n'
             '    return np.where(x == 0.0, 1.0, np.sinh(safe)/safe)\n\n\n'
             'def yukawa_potential(r, radii, n_e, m_mediator, alpha):\n'
             '    """V_emu [eV] at radii `r`, from a ball of electrons.\n\n'
             '    `r`, `radii` and 1/`m_mediator` in eV^-1, `n_e` in eV^3,\n'
             '    and alpha = g\'^2/(4 pi).\n'
             '    """\n'
             '    f = radii*n_e\n\n'
             '    # The interior piece, accumulated outward from the centre.\n'
             '    a = f*radii*shc(m_mediator*radii)\n'
             '    inner = np.concatenate(([0.0], np.cumsum(\n'
             '        0.5*(a[1:]+a[:-1])*np.diff(radii))))\n\n'
             '    # The exterior piece, accumulated *inward* from the\n'
             '    # surface.  Taking it as (total - interior) instead is the\n'
             '    # obvious thing to write and is wrong: once exp(-m R) is\n'
             '    # small the tail is a difference of two nearly equal\n'
             '    # numbers, and it vanishes into round-off.\n'
             '    b = f*np.exp(-m_mediator*radii)\n'
             '    d = 0.5*(b[1:]+b[:-1])*np.diff(radii)\n'
             '    outer = np.concatenate((np.cumsum(d[::-1])[::-1], [0.0]))\n\n'
             '    r = np.asarray(r, dtype=float)\n'
             '    safe = np.where(r == 0.0, 1.0, r)   # the 1/r is removable\n'
             '    return 4.0*np.pi*alpha*np.where(\n'
             '        r == 0.0, np.interp(r, radii, outer),\n'
             '        np.exp(-m_mediator*safe)*np.interp(r, radii, inner)'
             '/safe\n'
             '        + shc(m_mediator*safe)*np.interp(r, radii, outer))\n\n\n'
             '# nu_e carries +1 and nu_mu carries -1 of the gauged charge.\n'
             'Q_LR = np.diag([1.0, -1.0, 0.0])\n\n'
             'R_E = gd.EARTH_RADIUS*KM        # the Earth\'s radius in eV^-1\n'
             'print("R_Earth = %.4e eV^-1" % R_E)\n'
             'print("a mediator of range R_Earth has m = %.4e eV" '
             '% (1.0/R_E))'),
        md('## Does the potential come out right?\n\n'
           'A uniform ball has a closed form at **any** mediator mass,\n\n'
           '$$ V_{e\\mu}(r) = \\frac{g\'^2 n_e}{m^2}\\left[ 1 - '
           '\\left(R + \\frac{1}{m}\\right) e^{-mR} '
           '\\frac{\\sinh(m r)}{r} \\right] , $$\n\n'
           'so this is a check against something external rather than a '
           'self-consistency test. It is worth doing before any physics: '
           'a quadrature that is quietly wrong would produce plots that '
           'look entirely reasonable.'),
        code('n_uniform = n_e_of_density(5.51)      # the Earth\'s mean\n'
             'radii = np.linspace(0.0, R_E, 20001)\n'
             'r = np.linspace(0.0, 1.0, 11)*R_E\n'
             'alpha_ref = 1.0e-52\n\n'
             'print("  m*R    max relative difference from the closed form")\n'
             'for m_r in (0.0, 1.0, 10.0, 100.0):\n'
             '    m = m_r/R_E\n'
             '    num = yukawa_potential(r, radii,\n'
             '                           n_uniform*np.ones_like(radii),\n'
             '                           m, alpha_ref)\n'
             '    if m_r == 0.0:\n'
             '        # The m -> 0 limit: an ordinary 1/r potential.\n'
             '        exact = (4.0*np.pi*alpha_ref*n_uniform\n'
             '                 * (R_E**2.0/2.0 - r**2.0/6.0))\n'
             '    else:\n'
             '        safe = np.where(r == 0.0, 1.0, r)\n'
             '        exact = (4.0*np.pi*alpha_ref*n_uniform/m**2.0)*np.where('
             '\n'
             '            r == 0.0, 1.0 - (R_E+1.0/m)*np.exp(-m*R_E)*m,\n'
             '            1.0 - (R_E+1.0/m)*np.exp(-m*R_E)'
             '*np.sinh(m*safe)/safe)\n'
             '    print("%6.1f   %.2e" % (m_r, np.max(np.abs(num/exact'
             '-1.0))))'),
        md('## The mediator\'s range decides what the potential looks '
           'like\n\n'
           'Now the Earth\'s own electrons, from PREM. The density jumps '
           'at a shell boundary, so the quadrature grid has to put a '
           'node on each one — the same rule that makes `earth_slabs` '
           'cut a chord at the boundaries, and the same rule behind '
           '`probabilities_3nu_profile`\'s warning about '
           'discontinuities. Ignoring it costs an order of accuracy: '
           'the error falls as $h$ rather than $h^2$.'),
        code('PREM_EDGES = np.concatenate(([0.0], earth.PREM_BOUNDARIES))\n\n\n'
             'def earth_electrons(n_per_shell=400):\n'
             '    """Radii [eV^-1] and n_e [eV^3] through the Earth."""\n'
             '    radii, n_e = [], []\n'
             '    for lo, hi in zip(PREM_EDGES[:-1], PREM_EDGES[1:]):\n'
             '        nodes = np.linspace(lo, hi, n_per_shell+1)\n'
             '        # Every boundary appears twice, once with the\n'
             '        # density on each side.  `density_prem` assigns a\n'
             '        # boundary radius to the shell *below* it, so the\n'
             '        # lower edge needs a nudge to land in this one.\n'
             '        inside = nodes.copy()\n'
             '        inside[0] = lo + 1.0e-9*(hi-lo)\n'
             '        radii.append(nodes*KM)\n'
             '        n_e.append(n_e_of_density(earth.density_prem(inside)))\n'
             '    return np.concatenate(radii), np.concatenate(n_e)\n\n\n'
             'R_EARTH_GRID, NE_EARTH = earth_electrons()\n'
             'N_E_EARTH = 4.0*np.pi*integral(R_EARTH_GRID**2.0*NE_EARTH,\n'
             '                               R_EARTH_GRID)\n'
             'print("quadrature nodes    : %d" % len(R_EARTH_GRID))\n'
             'print("electrons in the Earth : %.3e" % N_E_EARTH)'),
        md('The two limits are worth having in mind before looking at '
           'the figure.\n\n'
           '* **Short range**, $1/m \\ll R$: only the neighbourhood '
           'contributes and $V_{e\\mu} \\to g\'^2 n_e(r)/m^2$. The '
           'potential tracks the *local* density, which makes it '
           'degenerate with a non-standard interaction — `03`\'s '
           '$\\epsilon_{\\alpha\\beta}$ with a particular flavor '
           'structure.\n'
           '* **Long range**, $1/m \\gg R$: the whole body contributes '
           'and the potential is smooth, largest at the centre, and '
           'indifferent to where any individual electron sits.\n\n'
           'Only the second is genuinely new physics as far as this '
           'library is concerned, and only the second needs the profile '
           'machinery.'),
        code('ranges = [0.03, 0.1, 0.3, 1.0, 3.0, 10.0]\n'
             'r_plot = np.linspace(0.0, R_E, 400)\n\n'
             'fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.6))\n'
             'for k, f in enumerate(ranges):\n'
             '    v = yukawa_potential(r_plot, R_EARTH_GRID, NE_EARTH,\n'
             '                         1.0/(f*R_E), alpha_ref)\n'
             '    axes[0].plot(r_plot/R_E, v/v[0], color="C%d" % k,\n'
             '                 label=r"$1/m = %g\\,R$" % f)\n'
             'ne_plot = np.interp(r_plot, R_EARTH_GRID, NE_EARTH)\n'
             'axes[0].plot(r_plot/R_E, ne_plot/ne_plot[0], "k--", lw=1.4,\n'
             '             label=r"$n_e(r)/n_e(0)$")\n'
             'axes[0].set_xlim(0.0, 1.0)\n'
             'axes[0].set_xlabel("$r/R$")\n'
             'axes[0].set_ylabel(r"$V_{e\\mu}(r)/V_{e\\mu}(0)$")\n'
             'axes[0].set_title("Shape", fontsize=10)\n'
             'axes[0].legend(fontsize=7.5)\n\n'
             'spread = np.logspace(-1.5, 1.5, 60)\n'
             'centre = np.array([yukawa_potential(\n'
             '    np.array([0.0]), R_EARTH_GRID, NE_EARTH, 1.0/(f*R_E),\n'
             '    alpha_ref)[0] for f in spread])\n'
             'axes[1].loglog(spread, centre)\n'
             'axes[1].set_xlim(spread[0], spread[-1])\n'
             'axes[1].set_xlabel("Mediator range $1/m$ [$R$]")\n'
             'axes[1].set_ylabel(r"$V_{e\\mu}(0)$ [eV]")\n'
             'axes[1].set_title(r"Depth at the centre, '
             '$\\alpha\'=10^{-52}$", fontsize=10)\n'
             'fig.tight_layout()\n'
             'plt.show()'),
        md('The left panel is the whole argument in one picture. At '
           '$1/m = 0.03\\,R$ the potential has very nearly taken the '
           'shape of the dashed density curve, shell jumps and all; by '
           '$1/m = 10\\,R$ it has forgotten the profile entirely and is '
           'a smooth bowl. The right panel shows the crossover: '
           '$V_{e\\mu}(0)$ grows as $1/m^2$ while the range is short and '
           'saturates once the mediator reaches across the body, because '
           'there are no more electrons to enclose.\n\n'
           'Everything below uses $1/m = R$, one range per body, which '
           'is the awkward middle where neither limit applies.'),
        md('## Profile 1: an invented body\n\n'
           'A body of the Earth\'s radius with '
           '$\\rho(r) = \\rho_c\\left[1-(r/R)^2\\right]$, and $\\rho_c$ '
           'chosen so that its **mean density, and so its total electron '
           'count, match the Earth\'s**. Nothing in the library knows '
           'this profile; `slabs.probabilities_3nu_profile` takes it as '
           'a callable.\n\n'
           'Matching the electron count is what makes the comparison '
           'with PREM later worth making: the two bodies hold the same '
           'electrons, arranged differently.'),
        code('RHO_C = 5.51/0.4            # <rho> = 0.4 rho_c for this shape\n'
             'BODY_RADII = np.linspace(0.0, R_E, 4001)\n'
             'BODY_RHO = RHO_C*(1.0 - (BODY_RADII/R_E)**2.0)\n'
             'BODY_NE = n_e_of_density(BODY_RHO)\n\n'
             'print("central density : %.3f g/cm^3" % RHO_C)\n'
             'print("mean density    : %.3f g/cm^3"\n'
             '      % (integral(BODY_RADII**2.0*BODY_RHO, BODY_RADII)\n'
             '         / integral(BODY_RADII**2.0, BODY_RADII)))\n'
             'print("electrons, relative to the Earth : %.4f"\n'
             '      % (4.0*np.pi*integral(BODY_RADII**2.0*BODY_NE,\n'
             '                            BODY_RADII)/N_E_EARTH))'),
        md('### The callable\n\n'
           '`hamiltonian_of(positions)` receives **all** the midpoints '
           'of a refinement at once, as distances along the trajectory '
           'in eV$^{-1}$, and must return one $3\\times3$ Hamiltonian '
           'per position. Writing it vectorised is not an optimisation '
           'here: the routine calls it once per refinement, so a loop '
           'inside would be the only loop in the calculation.\n\n'
           'The trajectory is a **diameter**, so the radius is '
           '$|x - R|$ — which is also why the potential it sees is '
           'symmetric about the midpoint, the same palindrome the '
           'compiled Earth kernels exploit.'),
        code('def body_hamiltonian_of(energy_ev, alpha, m):\n'
             '    """Builds the callable that `..._profile` wants."""\n'
             '    def h_of(x):\n'
             '        r = np.abs(x - R_E)                  # a diameter\n'
             '        rho = RHO_C*(1.0 - (r/R_E)**2.0)\n'
             '        h = np.asarray(hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '            H_VAC_3NU, energy_ev, earth.matter_potential(rho)))\n'
             '        if alpha:\n'
             '            v = yukawa_potential(r, BODY_RADII, BODY_NE, m,\n'
             '                                 alpha)\n'
             '            h = h + v[:, None, None]*Q_LR\n'
             '        return h\n'
             '    return h_of\n\n\n'
             'E_BODY = 1.5*GEV\n'
             'M_BODY = 1.0/R_E\n\n'
             '# No folklore exists about how many slabs an invented profile\n'
             '# needs, which is exactly the argument for a tolerance.\n'
             'for a in (0.0, 1.0e-52):\n'
             '    p, n = slabs.probabilities_3nu_profile(\n'
             '        body_hamiltonian_of(E_BODY, a, M_BODY), 2.0*R_E,\n'
             '        rtol=1.0e-4, atol=1.0e-8, return_n_slabs=True)\n'
             '    print("alpha\' = %-8.0e n_slabs = %4d  P_mue = %.6f  '
             'sum = %.15f"\n'
             '          % (a, n, p[3], sum(p[0:3])))'),
        md('The probabilities still sum to one from each initial flavor, '
           'which is the first thing to check of any Hamiltonian you '
           'have built yourself: an algebra slip that breaks Hermiticity '
           'shows up here and nowhere else.\n\n'
           'And the sampling is second-order, as it is for any profile — '
           'the midpoint rule does not care that the Hamiltonian is '
           'unfamiliar:'),
        code('reference = np.array(slabs.probabilities_3nu_profile(\n'
             '    body_hamiltonian_of(E_BODY, 1.0e-52, M_BODY), 2.0*R_E,\n'
             '    n_slabs=8192))\n\n'
             'previous = None\n'
             'print("n_slabs   max error   order")\n'
             'for n in (16, 32, 64, 128, 256):\n'
             '    p = np.array(slabs.probabilities_3nu_profile(\n'
             '        body_hamiltonian_of(E_BODY, 1.0e-52, M_BODY),\n'
             '        2.0*R_E, n_slabs=n))\n'
             '    err = np.max(np.abs(p - reference))\n'
             '    print("%7d   %.3e   %s"\n'
             '          % (n, err, "-" if previous is None\n'
             '             else "%.2f" % np.log2(previous/err)))\n'
             '    previous = err'),
        code('fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.4),\n'
             '                         sharex=True)\n'
             'axes[0].plot(BODY_RADII/R_E, BODY_RHO, label="invented body")\n'
             'axes[0].plot(R_EARTH_GRID/R_E,\n'
             '             earth.density_prem(R_EARTH_GRID/KM), "--",\n'
             '             label="PREM")\n'
             'axes[0].set_ylabel(r"$\\rho$ [g cm$^{-3}$]")\n'
             'axes[0].set_ylim(0.0, 14.5)\n'
             'axes[0].legend()\n'
             'axes[0].set_title("Same electrons, arranged differently")\n\n'
             'axes[1].plot(BODY_RADII/R_E, yukawa_potential(\n'
             '    BODY_RADII, BODY_RADII, BODY_NE, M_BODY, alpha_ref),\n'
             '    label="invented body")\n'
             'axes[1].plot(R_EARTH_GRID/R_E, yukawa_potential(\n'
             '    R_EARTH_GRID, R_EARTH_GRID, NE_EARTH, M_BODY, alpha_ref),\n'
             '    "--", label="PREM")\n'
             'axes[1].set_xlim(0.0, 1.0)\n'
             'axes[1].set_xlabel("$r/R$")\n'
             'axes[1].set_ylabel(r"$V_{e\\mu}$ [eV]")\n'
             'axes[1].legend()\n'
             'fig.tight_layout()\n'
             'plt.show()'),
        code('v_body = yukawa_potential(BODY_RADII, BODY_RADII, BODY_NE,\n'
             '                          M_BODY, alpha_ref)\n'
             'v_prem = yukawa_potential(BODY_RADII, R_EARTH_GRID, NE_EARTH,\n'
             '                          M_BODY, alpha_ref)\n'
             'rho_prem = earth.density_prem(BODY_RADII[:-1]/KM)\n\n'
             'print("largest ratio between the two densities  : %.2f"\n'
             '      % np.max(BODY_RHO[:-1]/rho_prem))\n'
             'print("largest ratio between the two potentials : %.3f"\n'
             '      % np.max(v_body/v_prem))\n'
             'print("at the surface: %.4e vs %.4e eV" % (v_body[-1],\n'
             '                                          v_prem[-1]))'),
        md('The two densities differ by nearly a factor of two through '
           'the middle of the body — PREM\'s core-mantle jump against a '
           'smooth parabola — while the potentials they source differ by '
           'under a tenth, and at the surface they agree to one part in '
           'a hundred. '
           'The coincidence is not luck: at the surface there is no '
           'exterior piece left, so the potential is fixed by the '
           'enclosed electron count, which is what was matched. For a '
           'massless mediator it would be exact.\n\n'
           'That is the non-locality made quantitative. A long-range '
           'force is sensitive to how many electrons there are and only '
           'weakly to where they sit, which is the exact opposite of '
           '$V_{CC}$.'),
        code('alphas = np.logspace(-54.0, -51.0, 40)\n'
             'p_body = np.array([slabs.probabilities_3nu_profile(\n'
             '    body_hamiltonian_of(E_BODY, a, M_BODY), 2.0*R_E,\n'
             '    rtol=1.0e-4, atol=1.0e-8)[3] for a in alphas])\n'
             'p_standard = slabs.probabilities_3nu_profile(\n'
             '    body_hamiltonian_of(E_BODY, 0.0, M_BODY), 2.0*R_E,\n'
             '    rtol=1.0e-4, atol=1.0e-8)[3]\n\n'
             'fig, ax = plt.subplots()\n'
             'ax.semilogx(alphas, p_body, label="with the long-range term")\n'
             'ax.axhline(p_standard, color="0.5", ls="--",\n'
             '           label="standard matter only")\n'
             'ax.set_xlim(alphas[0], alphas[-1])\n'
             'ax.set_xlabel(r"$\\alpha\' = g\'^2/4\\pi$")\n'
             'ax.set_ylabel(r"$P_{\\mu e}$")\n'
             'ax.set_title("Through the invented body, E = 1.5 GeV, "\n'
             '             r"$1/m = R$")\n'
             'ax.legend()\n'
             'plt.show()'),
        md('The left-hand end of that curve lies exactly on the '
           'standard-matter line, which is why the line looks missing: '
           'it is the $\\alpha\' \\to 0$ check, and it passes.\n\n'
           'What happens after that is not simple suppression. '
           '$P_{\\mu e}$ falls by a factor of forty, bottoms out near '
           '$\\alpha\' = 6\\times10^{-52}$, and climbs again — the extra '
           'splitting is re-tuning the interference between the slabs, '
           'and past the minimum a different oscillation is coming into '
           'phase. A bound built by assuming the effect is monotonic in '
           'the coupling would go wrong here, which is the sort of thing '
           'an exact calculation is for.'),
        md('## Profile 2: an exponential body\n\n'
           'A solar-like body: $\\rho(r) = \\rho_0 e^{-r/r_0}$ with '
           '$\\rho_0 = 150$ g cm$^{-3}$ and $r_0 = R_\\odot/10.54$, the '
           'same fit `14` uses, and the trajectory is a radius from the '
           'centre outward.\n\n'
           '**How this differs from `14`.** That notebook asks what the '
           '*standard* Hamiltonian does in this profile, validates the '
           'slab calculation against the analytic adiabatic MSW result, '
           'and argues that a Magnus method is the better tool. All of '
           'that stands, and none of it is repeated here. What is new is '
           'that the Hamiltonian is no longer standard — and that the '
           'body is now interesting for a second reason. A long-range '
           'force integrates the whole object, so the Sun\'s enormous '
           'electron count matters in a way it does not for $V_{CC}$:'),
        code('R_SUN_KM = 6.957e5\n'
             'RHO_SUN_C = 150.0\n'
             'R_SCALE_KM = R_SUN_KM/10.54\n'
             'R_SUN = R_SUN_KM*KM\n\n'
             'SUN_RADII = np.linspace(0.0, R_SUN, 8001)\n'
             'SUN_NE = n_e_of_density(\n'
             '    RHO_SUN_C*np.exp(-SUN_RADII/(R_SCALE_KM*KM)))\n'
             'N_E_SUN = 4.0*np.pi*integral(SUN_RADII**2.0*SUN_NE, SUN_RADII)\n'
             'V_CC_SUN_CENTRE = earth.matter_potential(RHO_SUN_C)\n\n'
             'v_sun = yukawa_potential(np.array([0.0]), SUN_RADII, SUN_NE,\n'
             '                         1.0/R_SUN, alpha_ref)[0]\n'
             'v_earth = yukawa_potential(np.array([0.0]), R_EARTH_GRID,\n'
             '                           NE_EARTH, 1.0/R_E, alpha_ref)[0]\n\n'
             'print("central density, Sun / Earth\'s core : %.1f"\n'
             '      % (RHO_SUN_C/earth.density_prem(0.0)))\n'
             'print("electrons,      Sun / Earth         : %.3e"\n'
             '      % (N_E_SUN/N_E_EARTH))\n'
             'print("V_emu(0),       Sun / Earth         : %.0f"\n'
             '      % (v_sun/v_earth))\n'
             'print()\n'
             'print("V_emu(0) at alpha\' = 1e-52 : %.3e eV" % v_sun)\n'
             'print("V_CC at the centre         : %.3e eV"\n'
             '      % V_CC_SUN_CENTRE)'),
        md('Eleven times the central density, but seven thousand times '
           'the long-range potential. A local probe and a non-local one '
           'rank the two bodies completely differently, and that '
           'asymmetry is the reason long-range forces are hunted with '
           'astrophysical sources rather than in a laboratory.'),
        code('def sun_hamiltonian_of(energy_ev, alpha, m):\n'
             '    def h_of(x):\n'
             '        r = np.minimum(x, R_SUN)             # a radius\n'
             '        rho = RHO_SUN_C*np.exp(-r/(R_SCALE_KM*KM))\n'
             '        h = np.asarray(hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '            H_VAC_3NU, energy_ev,\n'
             '            earth.matter_potential(rho)))\n'
             '        if alpha:\n'
             '            v = yukawa_potential(r, SUN_RADII, SUN_NE, m,\n'
             '                                 alpha)\n'
             '            h = h + v[:, None, None]*Q_LR\n'
             '        return h\n'
             '    return h_of\n\n\n'
             'E_SUN = 1.0e7                        # 10 MeV\n'
             'p, n_tol = slabs.probabilities_3nu_profile(\n'
             '    sun_hamiltonian_of(E_SUN, 1.0e-54, 1.0/R_SUN), R_SUN,\n'
             '    atol=1.0e-4, n_max=1 << 17, return_n_slabs=True)\n'
             'print("atol = 1e-4 asks for %d slabs, and gives P_ee = %.6f"\n'
             '      % (n_tol, p[0]))'),
        md('Eight thousand slabs for one number — `14`\'s wall, and it '
           'is unmoved by the extra term. But that number is also '
           'meaningless on its own, for `14`\'s other reason: at 10 MeV '
           'the probability oscillates far faster than any detector, or '
           'any production region, resolves. The physical quantity is '
           'the average over a band, so that is what to plot against the '
           'coupling.\n\n'
           'One practical note on tolerances, since it looks like a '
           'contradiction. The scan below fixes `n_slabs=4096`, half of '
           'what the tolerance search just asked for, because an '
           '*average* over two hundred energies tolerates more error per '
           'point than a single reported probability does — the errors '
           'are not correlated across the band. Checking that 4096 and '
           '8192 give the same average is a two-line experiment and it '
           'is how that number was chosen.'),
        code('def averaged_p_ee(alpha, n_slabs=4096, n_energies=200):\n'
             '    """<P_ee> over a 2% band, on a fixed energy grid."""\n'
             '    band = np.linspace(0.99e7, 1.01e7, n_energies)\n'
             '    return np.mean([slabs.probabilities_3nu_profile(\n'
             '        sun_hamiltonian_of(e, alpha, 1.0/R_SUN), R_SUN,\n'
             '        n_slabs=n_slabs)[0] for e in band])\n\n\n'
             '# The analytic adiabatic MSW result, as in notebook 14, with\n'
             '# the three-flavor cos^4(theta13) correction.\n'
             'cos2th = 1.0 - 2.0*gd.S12_NO_BF**2.0\n'
             'x = 2.0*E_SUN*V_CC_SUN_CENTRE/gd.D21_NO_BF\n'
             'cos2th_m = (cos2th-x)/np.sqrt((cos2th-x)**2.0\n'
             '                              + (1.0-cos2th**2.0))\n'
             'c13sq = 1.0 - gd.S13_NO_BF**2.0\n'
             'p_adiabatic = (c13sq**2.0*0.5*(1.0 + cos2th_m*cos2th)\n'
             '               + gd.S13_NO_BF**4.0)\n\n'
             'p_sun_standard = averaged_p_ee(0.0)\n'
             'print("adiabatic MSW, with the c13^4 correction : %.4f"\n'
             '      % p_adiabatic)\n'
             'print("slabs, averaged, standard matter         : %.4f"\n'
             '      % p_sun_standard)\n\n'
             '# The claim that 4096 slabs is enough for an average,\n'
             '# checked rather than asserted -- and the second line\n'
             '# measures how much of the residual above is nothing but\n'
             '# the finite sampling of the band.\n'
             'print("the same average at 8192 slabs           : %.4f"\n'
             '      % averaged_p_ee(0.0, n_slabs=8192))\n'
             'print("at 4096 slabs but half the energies      : %.4f"\n'
             '      % averaged_p_ee(0.0, n_energies=100))'),
        md('The slab count is settled: 4096 and 8192 agree to the fourth '
           'decimal, so the discretisation is not what limits this. '
           'Halving the number of energies instead moves the average by '
           'a few parts in a thousand, comparable to the remaining gap '
           'to the analytic value — so the calculation agrees with the '
           'adiabatic MSW prediction as closely as a two-hundred-point '
           'sampling of the band can show, reproducing `14`\'s '
           'two-flavor validation with the $\\theta_{13}$ correction '
           'included. The machinery is right, so what it says about the '
           'new term can be believed.'),
        code('alphas_sun = np.logspace(-55.0, -52.0, 10)\n'
             'p_sun = np.array([averaged_p_ee(a) for a in alphas_sun])\n\n'
             'fig, ax = plt.subplots()\n'
             'ax.semilogx(alphas_sun, p_sun, "o-", ms=4,\n'
             '            label="with the long-range term")\n'
             'ax.axhline(p_sun_standard, color="0.5", ls="--",\n'
             '           label="standard matter only")\n'
             'ax.axhline(p_adiabatic, color="C1", ls=":",\n'
             '           label="analytic adiabatic MSW")\n'
             'ax.set_xlim(alphas_sun[0], alphas_sun[-1])\n'
             'ax.set_xlabel(r"$\\alpha\' = g\'^2/4\\pi$")\n'
             'ax.set_ylabel(r"$\\langle P_{ee} \\rangle$")\n'
             'ax.set_title("Solar-like body, 10 MeV averaged over a 2% '
             'band")\n'
             'ax.legend()\n'
             'plt.show()'),
        md('The average is far less sensitive than a single energy: it '
           'is flat until $\\alpha\' \\approx 5\\times10^{-54}$, where '
           '$V_{e\\mu}$ at the centre reaches about a third of $V_{CC}$ '
           'there, and only then does the resonance move enough to '
           'survive the averaging. That is worth stating plainly, '
           'because a single-energy calculation would suggest a reach '
           'roughly an order of magnitude better than the averaged one '
           'actually delivers.\n\n'
           'The direction is the physics. $\\langle P_{ee}\\rangle$ does '
           'not drift toward the vacuum-averaged value — it climbs past '
           'it, toward one. A large $V_{e\\mu}$ dominates the '
           '$e$–$\\mu$ block and makes $\\nu_e$ very nearly an '
           'eigenstate, so the MSW conversion is not merely detuned but '
           'switched off.\n\n'
           '**On the coupling values.** They are chosen to make the '
           'figures show something, exactly as `03`\'s NSI and LIV '
           'parameters are. Nothing here is a bound: a bound needs a '
           'real source spectrum, a detector response, and the rest of '
           'the Universe\'s electrons, which the reference has and this '
           'notebook does not.'),
        md('## Profile 3: the Earth\n\n'
           'And now the one that needs care. `earth.probabilities_3nu_earth` '
           'takes the energy-independent **vacuum** Hamiltonian and '
           'builds $H = H_{\\rm vac}/E + V_{CC} P_{ee}$ per slab '
           'itself. There is no argument to it that produces our extra '
           'term, so the convenience of a single call is not available. '
           'The path is four lines, and it is the honest answer to "how '
           'do I put *my* Hamiltonian through the Earth":\n\n'
           '```python\n'
           'widths_km, densities = earth.earth_slabs(costhz, n_per_segment)\n'
           'potentials = earth.matter_potential(densities)\n'
           'h = your_own_construction(potentials)      # (n_slabs, 3, 3)\n'
           'slabs.probabilities_3nu_slabs(h, widths_km*KM)\n'
           '```\n\n'
           '`earth_slabs` does the part worth reusing — it cuts the '
           'chord at every PREM boundary it crosses, so no slab '
           'straddles a discontinuity — and hands back widths and '
           'densities you are free to interpret. The radius at each slab '
           'midpoint, which the new potential needs and the standard one '
           'does not, comes from '
           '`earth_radial_distance_from_depth`.'),
        code('def earth_chord(costhz, energy_ev, alpha, m,\n'
             '                n_per_segment=8):\n'
             '    """The nine probabilities along one chord."""\n'
             '    widths_km, densities = earth.earth_slabs(costhz,\n'
             '                                             n_per_segment)\n'
             '    h = np.asarray(hamiltonians3nu.hamiltonian_3nu_matter(\n'
             '        H_VAC_3NU, energy_ev,\n'
             '        earth.matter_potential(densities)))\n\n'
             '    if alpha:\n'
             '        edges = np.concatenate(([0.0], np.cumsum(widths_km)))\n'
             '        mid_km = 0.5*(edges[:-1]+edges[1:])\n'
             '        r_mid = earth.earth_radial_distance_from_depth(\n'
             '            costhz, mid_km)*KM\n'
             '        v = yukawa_potential(r_mid, R_EARTH_GRID, NE_EARTH,\n'
             '                             m, alpha)\n'
             '        h = h + v[:, None, None]*Q_LR\n\n'
             '    return np.array(slabs.probabilities_3nu_slabs(\n'
             '        h, widths_km*KM))\n\n\n'
             '# Switch the new term off and this must reproduce the\n'
             '# library\'s own Earth routine.  Not "to within a tolerance":\n'
             '# it is the same construction, so it is the same arithmetic.\n'
             'worst = 0.0\n'
             'for cz in (-1.0, -0.8, -0.35, -0.05):\n'
             '    for e_gev in (1.0, 5.0, 10.0, 50.0):\n'
             '        mine = earth_chord(cz, e_gev*GEV, 0.0, 0.0)\n'
             '        theirs = np.array(earth.probabilities_3nu_earth(\n'
             '            H_VAC_3NU, e_gev*GEV, cz))\n'
             '        worst = max(worst, np.max(np.abs(mine-theirs)))\n'
             'print("max |difference| over 16 (angle, energy) pairs: %.1e"\n'
             '      % worst)'),
        md('Zero. The recipe is not an approximation to the library\'s '
           'Earth routine, it *is* the library\'s Earth routine with the '
           'Hamiltonian construction pulled out into the open — which is '
           'the only reason it is safe to add a term to.\n\n'
           '### How finely to cut the chord\n\n'
           '`earth.slabs_for_tolerance` sizes the subdivision for a '
           'given accuracy, but it sizes it for the *standard* '
           'Hamiltonian, since that is the only one it can build. It is '
           'still the right starting point — the new potential is '
           'smooth on the scale of a PREM shell, so it asks less of the '
           'discretisation than the density jumps do — and the way to '
           'confirm that is to double the subdivision and watch the '
           'answer not move.'),
        code('E_SCAN = np.logspace(0.0, 2.0, 300)*GEV\n'
             'ALPHA_EARTH = 1.0e-52\n'
             'M_EARTH = 1.0/R_E\n\n'
             'n_per = earth.slabs_for_tolerance(H_VAC_3NU, E_SCAN, -1.0,\n'
             '                                  rtol=1.0e-3)\n'
             'print("slabs_for_tolerance, rtol 1e-3 : %d per segment"\n'
             '      % n_per)\n\n'
             'coarse = earth_chord(-1.0, 10.0*GEV, ALPHA_EARTH, M_EARTH,\n'
             '                     n_per_segment=n_per)\n'
             'fine = earth_chord(-1.0, 10.0*GEV, ALPHA_EARTH, M_EARTH,\n'
             '                   n_per_segment=2*n_per)\n'
             'print("doubling it moves the nine probabilities by %.2e"\n'
             '      % np.max(np.abs(coarse-fine)))'),
        code('fig, ax = plt.subplots()\n'
             'for a in (0.0, 1.0e-53, 1.0e-52, 3.0e-52):\n'
             '    p = np.array([earth_chord(-1.0, e, a, M_EARTH,\n'
             '                              n_per_segment=n_per)[3]\n'
             '                  for e in E_SCAN])\n'
             '    ax.semilogx(E_SCAN/GEV, p,\n'
             '                lw=1.8 if a == 0.0 else 1.2,\n'
             '                ls="--" if a == 0.0 else "-",\n'
             '                color="0.35" if a == 0.0 else None,\n'
             '                label="standard matter" if a == 0.0\n'
             '                else r"$\\alpha\' = %.0e$" % a)\n'
             'ax.set_xlim(E_SCAN[0]/GEV, E_SCAN[-1]/GEV)\n'
             'ax.set_xlabel("Energy [GeV]")\n'
             'ax.set_ylabel(r"$P_{\\mu e}$")\n'
             'ax.set_title("Straight through the Earth, "\n'
             '             r"$\\cos\\theta_z = -1$, $1/m = R$")\n'
             'ax.legend()\n'
             'plt.show()'),
        md('The long-range term adds to the $e$–$\\mu$ splitting, so it '
           'bites where that splitting is comparable to the vacuum one: '
           'between roughly 2 and 10 GeV it shifts and reshapes the '
           'resonance structure, moving the main peak and lifting the '
           'shoulder above it. Above about 20 GeV nothing happens, and '
           'for a reason worth keeping: there the standard matter term '
           'already dominates and has already suppressed the '
           'transition, so one more diagonal term changes nothing.\n\n'
           'Note also what is *not* happening. Even at '
           '$\\alpha\' = 3\\times10^{-52}$ the peak is displaced rather '
           'than destroyed — the effect is a rearrangement of where the '
           'oscillation sits, not an overall damping.\n\n'
           'Every direction at once shows where the sensitivity actually '
           'lives. This is an oscillogram of the *difference*, and each '
           'point of it is a separate hand-built chord — one call per '
           'point, because the batched entry points build the standard '
           'Hamiltonian and cannot be borrowed here. It is the one place '
           'where declining the library\'s Earth routine has a real '
           'cost.'),
        code('cz_grid = np.linspace(-1.0, -0.05, 70)\n'
             'e_grid = np.logspace(0.0, 2.0, 70)*GEV\n\n'
             'shift = np.array([\n'
             '    [earth_chord(cz, e, ALPHA_EARTH, M_EARTH,\n'
             '                 n_per_segment=n_per)[3]\n'
             '     - earth_chord(cz, e, 0.0, 0.0, n_per_segment=n_per)[3]\n'
             '     for e in e_grid] for cz in cz_grid])\n\n'
             'fig, ax = plt.subplots(figsize=(7.2, 4.4))\n'
             'limit = np.max(np.abs(shift))\n'
             'mesh = ax.pcolormesh(e_grid/GEV, cz_grid, shift,\n'
             '                     cmap="RdBu_r", vmin=-limit, vmax=limit,\n'
             '                     shading="auto")\n'
             'ax.set_xscale("log")\n'
             'ax.set_xlabel("Energy [GeV]")\n'
             'ax.set_ylabel(r"$\\cos\\theta_z$")\n'
             'ax.set_title(r"Change in $P_{\\mu e}$ from a long-range "\n'
             '             r"force, $\\alpha\' = 10^{-52}$")\n'
             'ax.grid(False)\n'
             'fig.colorbar(mesh, ax=ax, label=r"$\\Delta P_{\\mu e}$")\n'
             'fig.tight_layout()\n'
             'plt.show()'),
        md('The signal is confined to the long chords and to a few GeV '
           '— core-crossing trajectories, where both the path and the '
           'enclosed electron count are largest. Short chords barely '
           'enter the body and see a small fraction of it.'),
        md('## Putting your own Hamiltonian through this\n\n'
           'Nothing above is specific to a long-range force. The pattern '
           'is the same whatever the extra term is:\n\n'
           '1. **Write the matrix.** Any $3\\times3$ Hermitian addition '
           'to $H_{\\rm vac}/E + V_{CC}P_{ee}$ will do. Only the '
           'traceless part matters (`01`), so a flavor-universal piece '
           'can be dropped — and if it is enormous, it *should* be, or '
           'it will eat the significant digits of the physics.\n'
           '2. **For a body, use '
           '`slabs.probabilities_{2,3,4}nu_profile`** with a vectorised '
           'callable, and let `rtol`/`atol` choose the subdivision. An '
           'invented profile has no folklore about how many slabs it '
           'needs.\n'
           '3. **For the Earth, use `earth.earth_slabs`** and hand the '
           'result to `slabs.probabilities_3nu_slabs` yourself. The '
           'boundary cuts are the part worth reusing; the Hamiltonian is '
           'yours.\n'
           '4. **Check it.** Probabilities sum to one from every initial '
           'flavor; switching the new term off must reproduce the '
           'standard routine exactly; and the discretisation must '
           'converge at second order. All three are cheap, and all three '
           'were run above.\n\n'
           'Where the profile varies smoothly over many oscillation '
           'lengths, `14`\'s advice still applies and applies here too: '
           'the slab count is set by the oscillation, not by the '
           'Hamiltonian, and a Magnus-type method is the better tool. '
           'The extra term changes none of that.'),
    ])


# ---------------------------------------------------------------------------
# Reading order, and the footer that makes it navigable
# ---------------------------------------------------------------------------
#
# Twenty notebooks described as "numbered in reading order" is only a
# reading order if a reader can follow it from inside one.  Before this, no
# notebook linked to any other, none linked to the documentation, and none
# said what to read next -- so the order existed in the README and nowhere a
# reader actually was.
#
# Declared once here rather than written into twenty footers by hand, so
# that inserting a notebook is one edit and cannot leave a dangling link.
# `test_file_tree.py` already guarantees the filenames exist; the assertion
# in `add_footers` guarantees this list matches them.

READING_ORDER = [
    ('01_basics.ipynb', 'Basics',
     'units, one probability, and why to pass arrays'),
    ('02_vacuum_oscillations.ipynb', 'Oscillations in vacuum',
     'the probabilities against baseline and against energy'),
    ('03_matter_nsi_liv.ipynb', 'Matter, NSI, and LIV',
     'constant-density matter and two kinds of new physics'),
    ('04_oscillogram.ipynb', 'Oscillograms',
     'a two-dimensional map in a single call'),
    ('05_biprobability.ipynb', 'Bi-probability plots',
     'CP violation, as an ellipse'),
    ('06_earth_and_prem.ipynb', 'The Earth: PREM, chords, and slabs',
     'how a varying profile becomes a sequence of exact pieces'),
    ('07_earth_probabilities.ipynb', 'Probabilities through the Earth',
     'zenith-angle scans, an Earth oscillogram, and real baselines'),
    ('08_unusual_density_profiles.ipynb', 'Unusual density profiles',
     'castle walls, and why the arrangement of matter matters'),
    ('09_performance.ipynb', 'Performance',
     'looping versus broadcasting, measured on your machine'),
    ('10_paper_figures.ipynb', "The paper's figures",
     'the two figures of arXiv:1904.12391, reproduced'),
    ('11_exact_vs_approximations.ipynb', 'Exact versus the approximations',
     'where the familiar formulas break, and by how much'),
    ('12_ordering_and_octant.ipynb', 'Mass ordering and the octant',
     'two open questions, and how they show up'),
    ('13_antineutrinos.ipynb', 'Antineutrinos, done properly',
     'conjugate *and* flip -- and two ways to get it half right'),
    ('14_solar_and_adiabatic_msw.ipynb', 'Solar neutrinos and the MSW '
     'resonance', 'the adiabatic resonance, and the limits of slabbing'),
    ('15_numerical_edge_cases.ipynb', 'Numerical edge cases',
     'degeneracies, and what returns a number instead of NaN'),
    ('16_four_neutrinos.ipynb', 'Four neutrinos, and a sterile state',
     'a 3+1 system through SU(4), and why the method stops there'),
    ('17_cross_checks.ipynb', 'Cross-checks with other codes',
     'corroboration from nuSQuIDS and Zaglauer-Schwarzer'),
    ('18_evolution_operator.ipynb',
     'The evolution operator and the SU(n) coefficients',
     'the machinery underneath, for composing and extending'),
    ('19_animations.ipynb', 'Animated scenes',
     'four sweeps, as stills, and the reel they came from'),
    ('20_arbitrary_hamiltonian.ipynb',
     'An arbitrary Hamiltonian, through three profiles',
     'a long-range force carried through an invented body, an '
     'exponential one, and the Earth'),
]

DOCS = 'https://mbustama.github.io/NuOscProbExact'


def add_footers():
    """Appends a navigation footer to every notebook.

    Each carries the previous and next notebook and a pointer to the API
    reference, so a reader who arrives at any one of them -- which is what
    happens when a search engine or a colleague sends them a link -- can
    find both the path through and the underlying documentation.
    """
    assert set(name for name, _, _ in READING_ORDER) == set(books), (
        'READING_ORDER and the notebooks built here have diverged')

    for index, (name, _, _) in enumerate(READING_ORDER):
        parts = []
        if index:
            previous, title, _ = READING_ORDER[index-1]
            parts.append('**Previous:** [%s](%s)' % (title, previous))
        if index + 1 < len(READING_ORDER):
            following, title, blurb = READING_ORDER[index+1]
            parts.append('**Next:** [%s](%s) --- %s'
                         % (title, following, blurb))

        parts.append('[API reference](%s/functions.html) &middot; '
                     '[Numerical recipes](%s/recipes.html) &middot; '
                     '[All notebooks](.)' % (DOCS, DOCS))

        books[name].cells.append(md('---\n\n' + '  \n'.join(parts)))


def build():
    """Writes every notebook, executes it, and checks it kept its outputs."""
    from nbclient import NotebookClient
    from nbclient.exceptions import CellExecutionError

    add_footers()

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
