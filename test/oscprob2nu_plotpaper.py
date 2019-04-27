# -*- coding: utf-8 -*-
r"""Produce the plot of 2nu probabilities vs. energy shown in the paper.

Contains a routine to generate and save the plot of two-neutrino
probabilities vs. energy that is included in the paper.

Routine listings
----------------

    * plot_probability_2nu_vs_energy_compare - Generates, saves the plot

References
----------

.. [1] Mauricio Bustamante, "Exact neutrino oscillation probabilities
   with arbitrary time-independent Hamiltonians", arXiv:1904.XXXXX.

Created: 2019/04/26 21:20
Last modified: 2019/04/26 21:20
"""

__version__ = "0.1"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@nbi.ku.dk"


from numpy import *
import numpy as np
from pylab import *
from matplotlib import *
import matplotlib as mpl

import sys
sys.path.append('../src')

import oscprob2nu
import hamiltonians2nu
from globaldefs import *


def plot_probability_2nu_vs_energy_compare(output_format='pdf',
    output_path='./fig/'):
    r"""Generates and saves a plot of 2nu probabilities vs. energy.

    Generates and saves a plot of two-neutrino probabilities vs.
    energy for oscillations in vacuum, matter, with NSI, and with
    CPT-odd LIV.  This is the same plot that is included in the paper.

    Parameters
    ----------
    output_format : str, optional
        File extension of the plot to save (e.g., 'pdf', 'png', 'jpg').
    output_path : str, optional
        File path where to save the plot.

    Returns
    -------
    None
        The plot is generated and saved.
    """
    # Baseline (DUNE)
    l = 1.3e3*CONV_KM_TO_INV_EV # [eV^{-1}]

    # Neutrino energies
    log10_energy_nu_min = log10(5.e-1) # [GeV]
    log10_energy_nu_max = log10(4.e1) # [GeV]
    log10_energy_nu_npts = 400
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

    fig, axes = plt.subplots(3, 1, figsize=[8,15])
    fig.subplots_adjust(hspace=0.05, wspace=0.05)

    h_vacuum_energy_indep = \
        hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(  S23_BF,
                                                                    D31_BF)

    prob_vacuum = [oscprob2nu.probabilities_2nu( \
                    np.multiply(1./x/1.e9, h_vacuum_energy_indep), l) \
                    for x in energy_nu]

    # Uncomment to compare to the probability computed with the standard
    # ocillation formula in vacuum
    # prob_vacuum_std = [hamiltonians2nu.probabilities_2nu_vacuum_std( \
    #                     S23_BF, D31_BF, x, l/CONV_KM_TO_INV_EV) \
    #                 for x in energy_nu]

    prob_matter = [oscprob2nu.probabilities_2nu( \
                    hamiltonians2nu.hamiltonian_2nu_matter( \
                        h_vacuum_energy_indep, x*1.e9, VCC_EARTH_CRUST), l)
                    for x in energy_nu]

    # Uncomment to compare to the probability computed with the standard
    # ocillation formula in matter
    # prob_matter_std = [hamiltonians2nu.probabilities_2nu_matter_std( \
    #                     S23_BF, D31_BF, VCC_EARTH_CRUST, x, l/CONV_KM_TO_INV_EV) \
    #                 for x in energy_nu]

    prob_nsi = [oscprob2nu.probabilities_2nu( \
                    hamiltonians2nu.hamiltonian_2nu_nsi( \
                        h_vacuum_energy_indep, x*1.e9, VCC_EARTH_CRUST, EPS_2),
                        l)
                for x in energy_nu]

    prob_liv = [oscprob2nu.probabilities_2nu( \
                    hamiltonians2nu.hamiltonian_2nu_liv( \
                        h_vacuum_energy_indep, x*1.e9, SXI12, B1, B3, LAMBDA),
                        l)
                for x in energy_nu]

    # Pee, Pem, Pmm
    for i, ax in enumerate(np.array(axes).reshape((1,3))[0]):

        if (i == 0): # Pee

            p_vacuum = [x[0] for x in prob_vacuum]
            # p_vacuum_std = [x[0] for x in prob_vacuum_std]
            p_matter = [x[0] for x in prob_matter]
            # p_matter_std = [x[0] for x in prob_matter_std]
            p_nsi = [x[0] for x in prob_nsi]
            p_liv = [x[0] for x in prob_liv]

            ylabel = r'$P_{\nu_e \to \nu_e}$'

        elif (i == 1): # Pme

            p_vacuum = [x[2] for x in prob_vacuum]
            # p_vacuum_std = [x[2] for x in prob_vacuum_std]
            p_matter = [x[2] for x in prob_matter]
            # p_matter_std = [x[2] for x in prob_matter_std]
            p_nsi = [x[2] for x in prob_nsi]
            p_liv = [x[2] for x in prob_liv]

            ylabel = r'$P_{\nu_\mu \to \nu_e}$'

        elif (i == 2): # Pmm

            p_vacuum = [x[3] for x in prob_vacuum]
            # p_vacuum_std = [x[3] for x in prob_vacuum_std]
            p_matter = [x[3] for x in prob_matter]
            # p_matter_std = [x[3] for x in prob_matter_std]
            p_nsi = [x[3] for x in prob_nsi]
            p_liv = [x[3] for x in prob_liv]

            ylabel = r'$P_{\nu_\mu \to \nu_\mu}$'
            ax.set_xlabel(r'Neutrino energy [GeV]', fontsize=27)

        ax.set_ylabel(ylabel, fontsize=27)

        ax.plot(energy_nu, p_vacuum, color='C0', ls='-', lw=3.0, zorder=1,
            label=r'Vacuum')
        # ax.plot(energy_nu, p_vacuum_std, color='k', ls='-', lw=3.0, zorder=1,
            # label=r'Vacuum std.')
        ax.plot(energy_nu, p_matter, color='C1', ls='--', lw=3.0, zorder=1,
            label=r'Matter')
        # ax.plot(energy_nu, p_matter_std, color='b', ls='--', lw=3.0, zorder=1,
            # label=r'Matter std.')
        ax.plot(energy_nu, p_nsi, color='C2', ls=':', lw=3.0, zorder=1,
            label=r'NSI')
        ax.plot(energy_nu, p_liv, color='C3', ls='-.', lw=3.0, zorder=1,
            label=r'CPT-odd LIV')

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
            ax_yticks_major = np.array([0.00, 0.25, 0.50, 0.75, 1.00])
            ax.set_yticks(ax_yticks_major, minor=False)
            ax_yticks_minor = np.array([0.05, 0.10, 0.15, 0.20, 0.30, 0.35,
                                        0.40, 0.45, 0.55, 0.60, 0.65, 0.70,
                                        0.80, 0.85, 0.90, 0.95])
            ax.set_yticks(ax_yticks_minor, minor=True)
            ax.set_ylim([0.0, 1.0])
            ax.legend(loc='lower right', ncol=1, frameon=False,
                columnspacing=1.)

        elif (i == 1):

            ax.set_xticklabels([])
            ax_yticks_major = np.array([0.00, 0.25, 0.50, 0.75])
            ax.set_yticks(ax_yticks_major, minor=False)
            ax_yticks_minor = np.array([0.05, 0.10, 0.15, 0.20, 0.30, 0.35,
                                        0.40, 0.45, 0.55, 0.60, 0.65, 0.70,
                                        0.80, 0.85, 0.90, 0.95])
            ax.set_yticks(ax_yticks_minor, minor=True)
            ax.set_ylim([0.0, 1.])

        elif (i == 2):

            ax_yticks_major = np.array([0.00, 0.25, 0.50, 0.75])
            ax.set_yticks(ax_yticks_major, minor=False)
            ax_yticks_minor = np.array([0.05, 0.10, 0.15, 0.20, 0.30, 0.35,
                                        0.40, 0.45, 0.55, 0.60, 0.65, 0.70,
                                        0.80, 0.85, 0.90, 0.95])
            ax.set_yticks(ax_yticks_minor, minor=True)
            ax.set_ylim([0.0, 1.0])

        pylab.savefig(output_path+'prob_2nu_vs_energy_compare.'+output_format,
            bbox_inches='tight', dpi=300)

    return

plot_probability_2nu_vs_energy_compare( output_format='pdf',
                                        output_path='../fig/')

