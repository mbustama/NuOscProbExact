#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = "0.1"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@nbi.ku.dk"


"""
3nutestsuite.py:
    Routines to run a suite of tests

Created: 2019/04/17 17:14
Last modified: 2019/04/17 17:14
"""

import oscprob3nu
import oscprob3nuhamiltonians
from globaldefs import *


def plot_probability_3nu_vs_l(case, energy_nu=1.e6, output_format='pdf'):

    # Baselines, L
    log10_l_min = 0.0 # [km]
    log10_l_max = log10(5.e2) # [km]
    log10_l_npts = 6000
    log10_l_val = np.linspace(log10_l_min, log10_l_max, log10_l_npts)
    l_val =[10.**x for x in log10_l_val]

    if (case == 'vacuum'):
        hamiltonian = (1./energy_nu) \
            * hamiltonian_vacuum_energy_independent(s12, s23, s13, dCP, D21,
                                                    D31)
        label_case = r'Vacuum'
        filename_output = 'prob_3nu_vacuum_vs_l'

    # Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt
    lst_prob = [oscprob3nu.probabilities_3nu(   hamiltonian,
                                                l/CONV_KM_TO_INV_EV) \
                for l in l_val]
    lst_prob_ee = [x[0] for x in lst_prob]
    lst_prob_em = [x[1] for x in lst_prob]
    lst_prob_et = [x[2] for x in lst_prob]

    # Formatting
    mpl.rcParams['xtick.labelsize']=26
    mpl.rcParams['ytick.labelsize']=26
    mpl.rcParams['legend.fontsize']=26
    mpl.rcParams['legend.borderpad']=0.4
    mpl.rcParams['axes.labelpad']=10
    mpl.rcParams['ps.fonttype']=42
    mpl.rcParams['pdf.fonttype']=42

    fig = plt.figure(figsize=[9,9])
    ax = fig.add_subplot(1,1,1)

    ax.set_xlabel(r'Baseline $L$ [km]', fontsize=25)
    ax.set_ylabel(r'Probability', fontsize=25)

    yaxis_minor_locator = mpl.ticker.MultipleLocator(0.1)
    ax.yaxis.set_minor_locator(yaxis_minor_locator)

    ax.tick_params('both', length=10, width=2, which='major')
    ax.tick_params('both', length=5, width=1, which='minor')
    ax.tick_params(axis='both', which='major', pad=10, direction='in')
    ax.tick_params(axis='both', which='minor', pad=10, direction='in')
    ax.tick_params(axis='x', which='minor', bottom=True)
    ax.tick_params(axis='x', which='minor', top=True)
    ax.tick_params(axis='y', which='minor', left=True)
    ax.tick_params(axis='y', which='minor', right=True)
    ax.tick_params(bottom=True, top=True, left=True, right=True)

    ax.set_xlim([10.**log10_l_min, 10.**log10_l_max])
    ax.set_xscale('log')
    ax.set_ylim([0.0, 1.0])

    # Plot
    ax.plot(l_val, lst_prob_ee, label=r'$P_{\nu_e \to \nu_e}$', color='C0',
        zorder=1)
    ax.plot(l_val, lst_prob_em, label=r'$P_{\nu_e \to \nu_\mu}$', color='C1',
        zorder=1)
    ax.plot(l_val, lst_prob_et, label=r'$P_{\nu_e \to \nu_\tau}$', color='C2',
        zorder=1)

    # Legend
    ax.legend(loc='center left', frameon=False, ncol=1)

    pylab.savefig(filename_output+'.'+output_format, bbox_inches='tight',
        dpi=300)

    return




def plot_probability_3nu_vacuum_vs_l(output_format='pdf'):

    # Best-fit values of mixing parameters, normal ordering, from 1708.01186
    s12 = sqrt(3.21e-1)
    s23 = sqrt(4.30e-1)
    s13 = sqrt(2.155e-2)
    dCP = 1.4*np.pi # [rad]
    D21 = 7.56e-5 # [eV^2]
    D31 = 2.55e-3 # [eV^2]

    # Neutrino energy, E
    energy_nu = 1.e6 # [eV]

    # Baselines, L
    log10_l_min = 0.0 # [km]
    log10_l_max = log10(5.e2) # [km]
    log10_l_npts = 6000
    log10_l_val = np.linspace(log10_l_min, log10_l_max, log10_l_npts)
    l_val =[10.**x for x in log10_l_val]

    # Hamiltonian evaluated at a single energy; in the case of vacuum
    # oscillations we only need to evaluate it once
    hamiltonian_matrix = (1./energy_nu) \
        * hamiltonian_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31)

    # Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt
    lst_prob = [oscprob3nu.probabilities_3nu(   hamiltonian_matrix,
                                                l/CONV_KM_TO_INV_EV) \
                for l in l_val]
    lst_prob_ee = [x[0] for x in lst_prob]
    lst_prob_em = [x[1] for x in lst_prob]
    lst_prob_et = [x[2] for x in lst_prob]

    # Plot
    mpl.rcParams['xtick.labelsize']=26
    mpl.rcParams['ytick.labelsize']=26
    mpl.rcParams['legend.fontsize']=26
    mpl.rcParams['legend.borderpad']=0.4
    mpl.rcParams['axes.labelpad']=10
    mpl.rcParams['ps.fonttype']=42
    mpl.rcParams['pdf.fonttype']=42

    fig = plt.figure(figsize=[9,9])
    ax = fig.add_subplot(1,1,1)

    ax.set_xlabel(r'Baseline $L$ [km]', fontsize=25)
    ax.set_ylabel(r'Probability in vacuum ($E = 1$ MeV)', fontsize=25)

    yaxis_minor_locator = mpl.ticker.MultipleLocator(0.1)
    ax.yaxis.set_minor_locator(yaxis_minor_locator)

    ax.tick_params('both', length=10, width=2, which='major')
    ax.tick_params('both', length=5, width=1, which='minor')
    ax.tick_params(axis='both', which='major', pad=10, direction='in')
    ax.tick_params(axis='both', which='minor', pad=10, direction='in')
    ax.tick_params(axis='x', which='minor', bottom=True)
    ax.tick_params(axis='x', which='minor', top=True)
    ax.tick_params(axis='y', which='minor', left=True)
    ax.tick_params(axis='y', which='minor', right=True)
    ax.tick_params(bottom=True, top=True, left=True, right=True)

    ax.set_xlim([10.**log10_l_min, 10.**log10_l_max])
    ax.set_xscale('log')
    ax.set_ylim([0.0, 1.0])

    # Plot
    ax.plot(l_val, lst_prob_ee, label=r'$P_{\nu_e \to \nu_e}$', color='C0',
        zorder=1)
    ax.plot(l_val, lst_prob_em, label=r'$P_{\nu_e \to \nu_\mu}$', color='C1',
        zorder=1)
    ax.plot(l_val, lst_prob_et, label=r'$P_{\nu_e \to \nu_\tau}$', color='C2',
        zorder=1)

    # Legend
    ax.legend(loc='center left', frameon=False, ncol=1)

    pylab.savefig('prob_3nu_vacuum_vs_l.'+output_format, bbox_inches='tight',
        dpi=300)

    return


def plot_probability_3nu_vacuum_vs_l_std(output_format='pdf'):

    # Best-fit values of mixing parameters, normal ordering, from 1708.01186
    s12 = sqrt(3.21e-1)
    s23 = sqrt(4.30e-1)
    s13 = sqrt(2.155e-2)
    dCP = 1.4*np.pi # [rad]
    D21 = 7.56e-5 # [eV^2]
    D31 = 2.55e-3 # [eV^2]

    # Neutrino energy, E
    energy_nu = 1.e6 # [eV]

    # Baselines, L
    log10_l_min = 0.0 # [km]
    log10_l_max = log10(5.e2) # [km]
    log10_l_npts = 6000
    log10_l_val = np.linspace(log10_l_min, log10_l_max, log10_l_npts)
    l_val =[10.**x for x in log10_l_val]

    U = pmns_mixing_matrix(s12, s23, s13, dCP)

    # Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt
    lst_prob = [probabilities_3nu_std(  U, D21, D31, energy_nu,
                                        l/CONV_KM_TO_INV_EV) \
                for l in l_val]
    lst_prob_ee = [x[0] for x in lst_prob]
    lst_prob_em = [x[1] for x in lst_prob]
    lst_prob_et = [x[2] for x in lst_prob]

    # Plot
    mpl.rcParams['xtick.labelsize']=26
    mpl.rcParams['ytick.labelsize']=26
    mpl.rcParams['legend.fontsize']=26
    mpl.rcParams['legend.borderpad']=0.4
    mpl.rcParams['axes.labelpad']=10
    mpl.rcParams['ps.fonttype']=42
    mpl.rcParams['pdf.fonttype']=42

    fig = plt.figure(figsize=[9,9])
    ax = fig.add_subplot(1,1,1)

    ax.set_xlabel(r'Baseline $L$ [km]', fontsize=25)
    ax.set_ylabel(r'Probability in vacuum ($E = 1$ MeV)', fontsize=25)

    yaxis_minor_locator = mpl.ticker.MultipleLocator(0.1)
    ax.yaxis.set_minor_locator(yaxis_minor_locator)

    ax.tick_params('both', length=10, width=2, which='major')
    ax.tick_params('both', length=5, width=1, which='minor')
    ax.tick_params(axis='both', which='major', pad=10, direction='in')
    ax.tick_params(axis='both', which='minor', pad=10, direction='in')
    ax.tick_params(axis='x', which='minor', bottom=True)
    ax.tick_params(axis='x', which='minor', top=True)
    ax.tick_params(axis='y', which='minor', left=True)
    ax.tick_params(axis='y', which='minor', right=True)
    ax.tick_params(bottom=True, top=True, left=True, right=True)

    ax.set_xlim([10.**log10_l_min, 10.**log10_l_max])
    ax.set_xscale('log')
    ax.set_ylim([0.0, 1.0])

    # Plot
    ax.plot(l_val, lst_prob_ee, label=r'$P_{\nu_e \to \nu_e}$', color='C0',
        zorder=1)
    ax.plot(l_val, lst_prob_em, label=r'$P_{\nu_e \to \nu_\mu}$', color='C1',
        zorder=1)
    ax.plot(l_val, lst_prob_et, label=r'$P_{\nu_e \to \nu_\tau}$', color='C2',
        zorder=1)

    # Legend
    ax.legend(loc='center left', frameon=False, ncol=1)

    pylab.savefig('prob_3nu_vacuum_std_vs_l.'+output_format, bbox_inches='tight',
        dpi=300)

    return

def plot_probability_3nu_matter_vs_l(output_format='pdf'):

    # Best-fit values of mixing parameters, normal ordering, from 1708.01186
    s12 = sqrt(3.21e-1)
    s23 = sqrt(4.30e-1)
    s13 = sqrt(2.155e-2)
    dCP = 1.4*np.pi # [rad]
    D21 = 7.56e-5 # [eV^2]
    D31 = 2.55e-3 # [eV^2]

    # Neutrino energy, E
    energy_nu = 1.e6 # [eV]

    # Baselines, L
    log10_l_min = 0.0 # [km]
    log10_l_max = log10(5.e2) # [km]
    log10_l_npts = 6000
    log10_l_val = np.linspace(log10_l_min, log10_l_max, log10_l_npts)
    l_val =[10.**x for x in log10_l_val]

    # Hamiltonian evaluated at a single energy
    h_vacuum = (1./energy_nu) \
        * hamiltonian_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31)
    h_matter = hamiltonian_matter(h_vacuum, A_EARTH_CRUST)

    # Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt
    lst_prob = [oscprob3nu.probabilities_3nu(   hamiltonian_matrix,
                                                l/CONV_KM_TO_INV_EV) \
                for l in l_val]
    lst_prob_ee = [x[0] for x in lst_prob]
    lst_prob_em = [x[1] for x in lst_prob]
    lst_prob_et = [x[2] for x in lst_prob]

    # Plot
    mpl.rcParams['xtick.labelsize']=26
    mpl.rcParams['ytick.labelsize']=26
    mpl.rcParams['legend.fontsize']=26
    mpl.rcParams['legend.borderpad']=0.4
    mpl.rcParams['axes.labelpad']=10
    mpl.rcParams['ps.fonttype']=42
    mpl.rcParams['pdf.fonttype']=42

    fig = plt.figure(figsize=[9,9])
    ax = fig.add_subplot(1,1,1)

    ax.set_xlabel(r'Baseline $L$ [km]', fontsize=25)
    ax.set_ylabel(r'Probability in matter ($E = 1$ MeV)', fontsize=25)

    yaxis_minor_locator = mpl.ticker.MultipleLocator(0.1)
    ax.yaxis.set_minor_locator(yaxis_minor_locator)

    ax.tick_params('both', length=10, width=2, which='major')
    ax.tick_params('both', length=5, width=1, which='minor')
    ax.tick_params(axis='both', which='major', pad=10, direction='in')
    ax.tick_params(axis='both', which='minor', pad=10, direction='in')
    ax.tick_params(axis='x', which='minor', bottom=True)
    ax.tick_params(axis='x', which='minor', top=True)
    ax.tick_params(axis='y', which='minor', left=True)
    ax.tick_params(axis='y', which='minor', right=True)
    ax.tick_params(bottom=True, top=True, left=True, right=True)

    ax.set_xlim([10.**log10_l_min, 10.**log10_l_max])
    ax.set_xscale('log')
    ax.set_ylim([0.0, 1.0])

    # Plot
    ax.plot(l_val, lst_prob_ee, label=r'$P_{\nu_e \to \nu_e}$', color='C0',
        zorder=1)
    ax.plot(l_val, lst_prob_em, label=r'$P_{\nu_e \to \nu_\mu}$', color='C1',
        zorder=1)
    ax.plot(l_val, lst_prob_et, label=r'$P_{\nu_e \to \nu_\tau}$', color='C2',
        zorder=1)

    # Legend
    ax.legend(loc='center left', frameon=False, ncol=1)

    pylab.savefig('prob_3nu_matter_vs_l.'+output_format, bbox_inches='tight',
        dpi=300)

    return

