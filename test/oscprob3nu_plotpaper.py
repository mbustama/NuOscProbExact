# -*- coding: utf-8 -*-
r"""Produce the plot of probabilities vs. energy shown in the paper.

Contains a routine to generate and save the plot of three-neutrino
probabilities vs. energy that is included in the paper.

Routine listings
----------------

    * plot_probability_vs_energy_compare - Generates and saves the plot

References
----------

.. [1] Mauricio Bustamante, "Exact neutrino oscillation probabilities
   with arbitrary time-independent Hamiltonians", arXiv:1904.XXXXX.

Created: 2019/04/17 18:08
Last modified: 2019/04/22 20:36
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

import oscprob3nu
import hamiltonians3nu
from globaldefs import *


def plot_probability_vs_energy_compare(output_format='pdf',
    output_path='./fig/'):
    r"""Generates and saves a plot of 3nu probabilities vs. energy.

    Generates and saves a plot of three-neutrino probabilities vs.
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

    fig, axes = plt.subplots(3, 1, figsize=[8,15])
    fig.subplots_adjust(hspace=0.05, wspace=0.05)

    h_vacuum_energy_indep = \
        hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(  S12_BF,
                                                                    S23_BF,
                                                                    S13_BF,
                                                                    DCP_BF,
                                                                    D21_BF,
                                                                    D31_BF)

    prob_vacuum = [oscprob3nu.probabilities_3nu( \
                    np.multiply(1./x/1.e9, h_vacuum_energy_indep), l) \
                    for x in energy_nu]

    prob_matter = [oscprob3nu.probabilities_3nu( \
                    hamiltonians3nu.hamiltonian_3nu_matter( \
                        h_vacuum_energy_indep, x*1.e9, VCC_EARTH_CRUST), l)
                    for x in energy_nu]

    # eps = eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt
    prob_nsi = [oscprob3nu.probabilities_3nu( \
                    hamiltonians3nu.hamiltonian_3nu_nsi( \
                        h_vacuum_energy_indep, x*1.e9, VCC_EARTH_CRUST, EPS_3),
                        l)
                for x in energy_nu]

    prob_liv = [oscprob3nu.probabilities_3nu( \
                    hamiltonians3nu.hamiltonian_3nu_liv( \
                        h_vacuum_energy_indep, x*1.e9, SXI12, SXI23, SXI13,
                        DXICP, B1, B2, B3, LAMBDA), l)
                for x in energy_nu]

    # Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt
    for i, ax in enumerate(np.array(axes).reshape((1,3))[0]):

        if (i == 0): # Pee

            p_vacuum = [x[0] for x in prob_vacuum]
            p_matter = [x[0] for x in prob_matter]
            p_nsi = [x[0] for x in prob_nsi]
            p_liv = [x[0] for x in prob_liv]

            ylabel = r'$P_{\nu_e \to \nu_e}$'

        elif (i == 1): # Pme

            p_vacuum = [x[3] for x in prob_vacuum]
            p_matter = [x[3] for x in prob_matter]
            p_nsi = [x[3] for x in prob_nsi]
            p_liv = [x[3] for x in prob_liv]

            ylabel = r'$P_{\nu_\mu \to \nu_e}$'

        elif (i == 2): # Pmm

            p_vacuum = [x[4] for x in prob_vacuum]
            p_matter = [x[4] for x in prob_matter]
            p_nsi = [x[4] for x in prob_nsi]
            p_liv = [x[4] for x in prob_liv]

            ylabel = r'$P_{\nu_\mu \to \nu_\mu}$'
            ax.set_xlabel(r'Neutrino energy [GeV]', fontsize=27)

        ax.set_ylabel(ylabel, fontsize=27)

        ax.plot(energy_nu, p_vacuum, color='C0', ls='-', lw=3.0, zorder=1,
            label=r'Vacuum')
        ax.plot(energy_nu, p_matter, color='C1', ls='--', lw=3.0, zorder=1,
            label=r'Matter')
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
            ax_yticks_major = np.array([0.90, 0.95, 1.0])
            ax.set_yticks(ax_yticks_major, minor=False)
            ax_yticks_minor = np.array([0.89, 0.91, 0.92,
                                        0.93, 0.94, 0.66, 0.97, 0.98, 0.99])
            ax.set_yticks(ax_yticks_minor, minor=True)
            ax.set_ylim([0.89, 1.0])
            ax.legend(loc='lower right', ncol=1, frameon=False,
                columnspacing=1.)

        elif (i == 1):

            ax.set_xticklabels([])
            ax_yticks_major = np.array([0.000, 0.025, 0.050, 0.075])
            ax.set_yticks(ax_yticks_major, minor=False)
            ax_yticks_minor = np.array([0.000, 0.005, 0.010, 0.015, 0.020,
                                        0.030, 0.035, 0.040, 0.045, 0.055,
                                        0.060, 0.065, 0.070, 0.080])
            ax.set_yticks(ax_yticks_minor, minor=True)
            ax.set_ylim([0.0, 0.085])

        elif (i == 2):

            ax_yticks_major = np.array([0.00, 0.25, 0.50, 0.75])
            ax.set_yticks(ax_yticks_major, minor=False)
            ax_yticks_minor = np.array([0.05, 0.10, 0.15, 0.20, 0.30, 0.35,
                                        0.40, 0.45, 0.55, 0.60, 0.65, 0.70,
                                        0.80, 0.85, 0.90, 0.95])
            ax.set_yticks(ax_yticks_minor, minor=True)
            ax.set_ylim([0.0, 1.0])

        pylab.savefig(output_path+'prob_vs_energy_compare.'+output_format,
            bbox_inches='tight', dpi=300)

    return


plot_probability_vs_energy_compare(output_format='pdf', output_path='../fig/')

# help(oscprob3nu)
# print(oscprob3nu.__doc__)
