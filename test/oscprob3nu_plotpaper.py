#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = "0.1"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@nbi.ku.dk"


"""
3nuplotpaper.py:
    Routine to generate the plot of probability vs. neutrino energy
    included in the paper

Created: 2019/04/17 18:08
Last modified: 2019/04/17 18:08
"""

from numpy import *
import numpy as np
from pylab import *
from matplotlib import *
import matplotlib as mpl

import sys
sys.path.append('../src')

import oscprob3nu
import hamiltonians3nu
from globaldefs import *


def plot_probability_vs_energy_compare(output_format='pdf'):

    # Baseline (DUNE)
    l = 1.3e3*CONV_KM_TO_INV_EV # [eV^{-1}]

    # Neutrino energies
    log10_energy_nu_min = log10(1.e-1) # [GeV]
    log10_energy_nu_max = log10(1.e1) # [GeV]
    log10_energy_nu_npts = 200
    log10_energy_nu = np.linspace(  log10_energy_nu_min,
                                    log10_energy_nu_max,
                                    log10_energy_nu_npts)
    energy_nu = [10.**x for x in log10_energy_nu] # [GeV]

    # Plot formatting
    mpl.rcParams['xtick.labelsize']=28
    mpl.rcParams['ytick.labelsize']=28
    mpl.rcParams['legend.fontsize']=21
    mpl.rcParams['legend.borderpad']=0.4
    mpl.rcParams['axes.labelpad']=10
    mpl.rcParams['ps.fonttype']=42
    mpl.rcParams['pdf.fonttype']=42

    fig, axes = plt.subplots(3, 1, figsize=[8,12])
    fig.subplots_adjust(hspace=0.05, wspace=0.05)

    h_vacuum_energy_indep = \
        hamiltonians3nu.hamiltonian3nu_vacuum_energy_independent(   S12_BF,
                                                                    S23_BF,
                                                                    S13_BF,
                                                                    DCP_BF,
                                                                    D21_BF,
                                                                    D31_BF)

    # eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt
    eps = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

    prob_vacuum = [oscprob3nu.probabilities_3nu( \
                    np.multiply(1./x/1.e9, h_vacuum_energy_indep), l) \
                    for x in energy_nu]

    prob_matter = [oscprob3nu.probabilities_3nu( \
                    hamiltonians3nu.hamiltonian3nu_matter( \
                        np.multiply(1./x/1.e9, h_vacuum_energy_indep),
                                    VCC_EARTH_CRUST),
                        l)
                    for x in energy_nu]

    prob_nsi = [oscprob3nu.probabilities_3nu( \
                    hamiltonians3nu.hamiltonian3nu_nsi( \
                        np.multiply(1./x/1.e9, h_vacuum_energy_indep),
                                    A_EARTH_CRUST, eps),
                        l)
                for x in energy_nu]

    # Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt

    for i, ax in enumerate(np.array(axes).reshape((1,3))[0]):

        if (i == 0): # Pee

            p_vacuum = [x[0] for x in prob_vacuum]
            p_matter = [x[0] for x in prob_matter]
            p_nsi = [x[0] for x in prob_nsi]

            ylabel = r'$P_{\nu_e \to \nu_e}$'

        elif (i == 1): # Pme

            p_vacuum = [x[3] for x in prob_vacuum]
            p_matter = [x[3] for x in prob_matter]
            p_nsi = [x[3] for x in prob_nsi]

            ylabel = r'$P_{\nu_\mu \to \nu_e}$'

        elif (i == 2): # Pmm

            p_vacuum = [x[4] for x in prob_vacuum]
            p_matter = [x[4] for x in prob_matter]
            p_nsi = [x[4] for x in prob_nsi]

            ylabel = r'$P_{\nu_\mu \to \nu_\mu}$'
            ax.set_xlabel(r'Neutrino energy [GeV]', fontsize=25)

        ax.set_ylabel(ylabel, fontsize=25)

        ax.plot(energy_nu, p_vacuum, color='C0', ls='-', lw=2., zorder=1,
            label=r'Vacuum')
        ax.plot(energy_nu, p_matter, color='C1', ls='--', lw=2., zorder=1,
            label=r'Matter')
        ax.plot(energy_nu, p_nsi, color='C2', ls=':', lw=2.5, zorder=1,
            label=r'NSI')

        ax.tick_params('both', length=10, width=2, which='major')
        ax.tick_params('both', length=5, width=1, which='minor')
        ax.tick_params(axis='both', which='major', pad=10, direction='in')
        ax.tick_params(axis='both', which='minor', pad=10, direction='in')
        ax.tick_params(axis='x', which='minor', bottom=True)
        ax.tick_params(axis='x', which='minor', top=True)
        ax.tick_params(axis='y', which='minor', left=True)
        ax.tick_params(axis='y', which='minor', right=True)
        ax.tick_params(bottom=True, top=True, left=True, right=True)

        ax.set_xlim([10.**log10_energy_nu_min, 10.**log10_energy_nu_max])
        ax.set_xscale('log')

        if (i == 0):

            ax.set_xticklabels([])
            ax_yticks_major = np.array([0.85, 0.90, 0.95, 1.0])
            ax.set_yticks(ax_yticks_major, minor=False)
            ax_yticks_minor = np.array([0.86, 0.87, 0.88, 0.89, 0.91, 0.92,
                                        0.93, 0.94, 0.66, 0.97, 0.98, 0.99])
            ax.set_yticks(ax_yticks_minor, minor=True)
            ax.set_ylim([0.85, 1.0])
            ax.legend(loc='lower right', ncol=2, frameon=False,
                columnspacing=1.)

        elif (i == 1):

            ax.set_xticklabels([])
            ax_yticks_major = np.array([0.000, 0.025, 0.050, 0.075])
            ax.set_yticks(ax_yticks_major, minor=False)
            ax_yticks_minor = np.array([0.000, 0.005, 0.010, 0.015, 0.020,
                                        0.030, 0.035, 0.040, 0.045, 0.055,
                                        0.060, 0.065, 0.070, 0.080, 0.085,
                                        0.090, 0.095])
            ax.set_yticks(ax_yticks_minor, minor=True)
            ax.set_ylim([0.0, 0.1])

        elif (i == 2):

            ax_yticks_major = np.array([0.00, 0.25, 0.50, 0.75])
            ax.set_yticks(ax_yticks_major, minor=False)
            ax_yticks_minor = np.array([0.05, 0.10, 0.15, 0.20, 0.30, 0.35,
                                        0.40, 0.45, 0.55, 0.60, 0.65, 0.70,
                                        0.80, 0.85, 0.90, 0.95])
            ax.set_yticks(ax_yticks_minor, minor=True)
            ax.set_ylim([0.0, 1.0])

        pylab.savefig('../fig/prob_vs_energy_compare.'+output_format,
            bbox_inches='tight', dpi=300)

    return


# plot_probability_vs_energy_compare()

# help(oscprob3nu)
# print(oscprob3nu.__doc__)
