# -*- coding: utf-8 -*-
r"""Routines to plot three-neutrino flavor-transition probabilities.

This module contains contains routines to plot three-neutrino
oscillation probabilities vs. the neutrino baseline and energy.  These
routines are used by run_testsuite.py to produce a suite of test plots.

Routine listings
----------------

    * plot_probability_3nu_vs_baseline - Plot probabilities vs. baseline
    * plot_probability_3nu_vs_energy - Plot probabilities vs. energy

Created: 2019/04/17 17:14
Last modified: 2019/04/22 18:01
"""

__version__ = "0.1"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


from numpy import *
import numpy as np
from pylab import *
from matplotlib import *
import matplotlib as mpl

import sys
sys.path.append('../src')
sys.path.append('../test')

import oscprob3nu
import hamiltonians3nu
from globaldefs import *


def plot_probability_3nu_vs_baseline(
                case, energy=1.e-1,
                log10_l_min=0.0, log10_l_max=3.0, log10_l_npts=6000,
                plot_prob_ee=True, plot_prob_em=True, plot_prob_et=True,
                plot_prob_me=False, plot_prob_mm=False, plot_prob_mt=False,
                plot_prob_te=False, plot_prob_tm=False, plot_prob_tt=False,
                output_filename='prob_vs_baseline', output_format='pdf',
                output_path='../fig/', legend_loc='center left',
                legend_ncol=1):
    r"""Generates and saves a plot of 3nu probabilities vs. baseline.

    Generates a plot of three-neutrino oscillation probabilities vs.
    baseline, for a fixed neutrino energy.  The probabilities to be
    plotted are turned on and off via the flags plot_prob_ee,
    plot_prob_em, etc.  (At least one of them must be True.)  The
    parameter 'case' selects between 'vacuum', 'matter', 'nsi', and
    'liv' (see below).  The plot is saved with the provided name and
    file format under the specified path.

    Parameters
    ----------
    case : str
        Not optional.  Must be one of the following: 'vacuum', 'matter',
        'nsi', or 'liv'.  In each case, the probabilities are computed
        using the default parameter values pulled from globaldefs.
    energy : float, optional
        Neutrino energy [GeV].
    log10_l_min : float, optional
        Log10 of the minimum baseline [km].
    log10_l_max : float, optional
        Log10 of the maximum baseline [km].
    log10_l_npts : int, optional
        Number of baseline values at which to compute the probabilities.
    plot_prob_ee : bool, optional
        True to plot Pee, False otherwise.
    plot_prob_em : bool, optional
        True to plot Pem, False otherwise.
    plot_prob_et : bool, optional
        True to plot Pet, False otherwise.
    plot_prob_me : bool, optional
        True to plot Pme, False otherwise.
    plot_prob_mm : bool, optional
        True to plot Pmm, False otherwise.
    plot_prob_mt : bool, optional
        True to plot Pmt, False otherwise.
    plot_prob_te : bool, optional
        True to plot Pte, False otherwise.
    plot_prob_tm : bool, optional
        True to plot Ptm, False otherwise.
    plot_prob_tt : bool, optional
        True to plot Ptt, False otherwise.
    output_filename : str, optional
        File name of plot to save (without the file extension).
    output_format : str, optional
        File extension of the plot to save (e.g., 'pdf', 'png', 'jpg').
    output_path : str, optional
        File path where to save the plot.
    legend_loc : str, optional
        Location of the legend in the plot.  Must be one of the allowed
        values of the plot routine of matplotlib.
    legend_ncol : int, optional
        Number of columns to include in the legend box.  Must be at
        least 1.

    Returns
    -------
    None
        The plot is generated and saved.
    """
    if (not plot_prob_ee) and (not plot_prob_em) and (not plot_prob_et) \
        and (not plot_prob_me) and (not plot_prob_mm) and (not plot_prob_mt) \
        and (not plot_prob_te) and (not plot_prob_tm) and (not plot_prob_tt):
        quit()

    # Baselines, L
    log10_l_val = np.linspace(log10_l_min, log10_l_max, log10_l_npts)
    l_val =[10.**x for x in log10_l_val]

    h_vacuum_energy_independent = \
        hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(  S12_BF,
                                                                    S23_BF,
                                                                    S13_BF,
                                                                    DCP_BF,
                                                                    D21_BF,
                                                                    D31_BF)

    if (case.lower() == 'vacuum'):

        hamiltonian = np.multiply(1./energy/1.e9, h_vacuum_energy_independent)
        label_case = r'Vacuum'

    elif (case.lower() == 'matter'):

        hamiltonian = hamiltonians3nu.hamiltonian_3nu_matter( \
                                                h_vacuum_energy_independent,
                                                energy*1.e9,
                                                VCC_EARTH_CRUST)
        label_case = r'Matter'

    elif (case.lower() == 'nsi'):

        hamiltonian = hamiltonians3nu.hamiltonian_3nu_nsi( \
                                                h_vacuum_energy_independent,
                                                energy*1.e9,
                                                VCC_EARTH_CRUST,
                                                EPS_3)
        label_case = r'NSI'

    elif (case.lower() == 'liv'):

        hamiltonian = hamiltonians3nu.hamiltonian_3nu_liv( \
                                                h_vacuum_energy_independent,
                                                energy*1.e9,
                                                SXI12, SXI23, SXI13, DXICP,
                                                B1, B2, B3, LAMBDA)
        label_case = r'CPT-odd LIV'


    # Each element of prob: [Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt]
    prob = [oscprob3nu.probabilities_3nu(   hamiltonian,
                                            l*CONV_KM_TO_INV_EV) \
            for l in l_val]
    prob_ee = [x[0] for x in prob]
    prob_em = [x[1] for x in prob]
    prob_et = [x[2] for x in prob]
    prob_me = [x[3] for x in prob]
    prob_mm = [x[4] for x in prob]
    prob_mt = [x[5] for x in prob]
    prob_te = [x[6] for x in prob]
    prob_tm = [x[7] for x in prob]
    prob_tt = [x[8] for x in prob]

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
    ax.set_ylabel(r'Three-neutrino probability', fontsize=25)

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
    if (plot_prob_ee):
        ax.plot(l_val, prob_ee, label=r'$P_{\nu_e \to \nu_e}$',
            color='C0', zorder=1)
    if (plot_prob_em):
        ax.plot(l_val, prob_em, label=r'$P_{\nu_e \to \nu_\mu}$',
            color='C1', zorder=1)
    if (plot_prob_et):
        ax.plot(l_val, prob_et, label=r'$P_{\nu_e \to \nu_\tau}$',
            color='C2', zorder=1)
    if (plot_prob_me):
        ax.plot(l_val, prob_me, label=r'$P_{\nu_\mu \to \nu_e}$',
            color='C3', zorder=1)
    if (plot_prob_mm):
        ax.plot(l_val, prob_mm, label=r'$P_{\nu_\mu \to \nu_\mu}$',
            color='C4', zorder=1)
    if (plot_prob_mt):
        ax.plot(l_val, prob_mt, label=r'$P_{\nu_\mu \to \nu_\tau}$',
            color='C5', zorder=1)
    if (plot_prob_te):
        ax.plot(l_val, prob_te, label=r'$P_{\nu_\tau \to \nu_e}$',
            color='C6', zorder=1)
    if (plot_prob_tm):
        ax.plot(l_val, prob_tm, label=r'$P_{\nu_\tau \to \nu_\mu}$',
            color='C7', zorder=1)
    if (plot_prob_tt):
        ax.plot(l_val, prob_tt, label=r'$P_{\nu_\tau \to \nu_\tau}$',
            color='C8', zorder=1)

    # Legend
    ax.legend(loc=legend_loc, frameon=False, ncol=legend_ncol)

    # Annotations
    ax.annotate( label_case, xy = (0.05, 0.86), \
        xycoords='axes fraction', color='k', fontsize=25,
        horizontalalignment='left', rotation=0, zorder=2 )
    ax.annotate( r'$\log_{10}(E/{\rm GeV}) = $' + \
        str(int(log10(energy)*100.)/100.),
        xy = (0.05, 0.80), xycoords='axes fraction', color='k', fontsize=25,
        horizontalalignment='left', rotation=0, zorder=2 )

    pylab.savefig(output_path+output_filename+'.'+output_format,
        bbox_inches='tight', dpi=300)

    plt.close()

    return


def plot_probability_3nu_vs_energy(
                case, baseline=1.3e3,
                log10_energy_min=-1.0, log10_energy_max=1.0,
                log10_energy_npts=200,
                plot_prob_ee=True, plot_prob_em=True, plot_prob_et=True,
                plot_prob_me=False, plot_prob_mm=False, plot_prob_mt=False,
                plot_prob_te=False, plot_prob_tm=False, plot_prob_tt=False,
                output_filename='prob_vs_energy', output_format='pdf',
                output_path='../fig/', legend_loc='center right',
                legend_ncol=1):
    r"""Generates and saves a plot of 3nu probabilities vs. energy.

    Generates a plot of three-neutrino oscillation probabilities vs.
    energy, for a fixed neutrino baseline.  The probabilities to be
    plotted are turned on and off via the flags plot_prob_ee,
    plot_prob_em, etc.  (At least one of them must be True.)  The
    parameter 'case' selects between 'vacuum', 'matter', 'nsi', and
    'liv' (see below).  The plot is saved with the provided name and
    file format under the specified path.

    Parameters
    ----------
    case : str
        Not optional.  Must be one of the following: 'vacuum', 'matter',
        'nsi', or 'liv'.  In each case, the probabilities are computed
        using the default parameter values pulled from globaldefs.
    baseline : float, optional
        Neutrino baseline [km].
    log10_energy_min : float, optional
        Log10 of the minimum energy [GeV].
    log10_energy_max : float, optional
        Log10 of the maximum energy [GeV].
    log10_energy_npts : int, optional
        Number of energy values at which to compute the probabilities.
    plot_prob_ee : bool, optional
        True to plot Pee, False otherwise.
    plot_prob_em : bool, optional
        True to plot Pem, False otherwise.
    plot_prob_et : bool, optional
        True to plot Pet, False otherwise.
    plot_prob_me : bool, optional
        True to plot Pme, False otherwise.
    plot_prob_mm : bool, optional
        True to plot Pmm, False otherwise.
    plot_prob_mt : bool, optional
        True to plot Pmt, False otherwise.
    plot_prob_te : bool, optional
        True to plot Pte, False otherwise.
    plot_prob_tm : bool, optional
        True to plot Ptm, False otherwise.
    plot_prob_tt : bool, optional
        True to plot Ptt, False otherwise.
    output_filename : str, optional
        File name of plot to save (without the file extension).
    output_format : str, optional
        File extension of the plot to save (e.g., 'pdf', 'png', 'jpg').
    output_path : str, optional
        File path where to save the plot.
    legend_loc : str, optional
        Location of the legend in the plot.  Must be one of the allowed
        values of the plot routine of matplotlib.
    legend_ncol : int, optional
        Number of columns to include in the legend box.  Must be at
        least 1.

    Returns
    -------
    None
        The plot is generated and saved.
    """
    if (not plot_prob_ee) and (not plot_prob_em) and (not plot_prob_et) \
        and (not plot_prob_me) and (not plot_prob_mm) and (not plot_prob_mt) \
        and (not plot_prob_te) and (not plot_prob_tm) and (not plot_prob_tt):
        quit()

    baseline = baseline*CONV_KM_TO_INV_EV # [eV^{-1}]

    # Neutrino energies
    log10_energy_val = np.linspace( log10_energy_min, log10_energy_max,
                                    log10_energy_npts)
    energy_val =[10.**x for x in log10_energy_val]

    h_vacuum_energy_independent = \
        hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(  S12_BF,
                                                                    S23_BF,
                                                                    S13_BF,
                                                                    DCP_BF,
                                                                    D21_BF,
                                                                    D31_BF)

    if (case.lower() == 'vacuum'):

        prob = [oscprob3nu.probabilities_3nu( \
                    np.multiply(1./energy/1.e9, h_vacuum_energy_independent),
                    baseline)  \
                for energy in energy_val]
        label_case = r'Vacuum'

    elif (case.lower() == 'matter'):

        prob = [oscprob3nu.probabilities_3nu( \
                    hamiltonians3nu.hamiltonian_3nu_matter( \
                                                h_vacuum_energy_independent,
                                                energy*1.e9,
                                                VCC_EARTH_CRUST),
                    baseline)  \
                for energy in energy_val]
        label_case = r'Matter'

    elif (case.lower() == 'nsi'):

        prob = [oscprob3nu.probabilities_3nu( \
                    hamiltonians3nu.hamiltonian_3nu_nsi( \
                                                h_vacuum_energy_independent,
                                                energy*1.e9,
                                                VCC_EARTH_CRUST,
                                                EPS_3),
                    baseline)  \
                for energy in energy_val]
        label_case = r'NSI'

    elif (case.lower() == 'liv'):

        prob = [oscprob3nu.probabilities_3nu( \
                    hamiltonians3nu.hamiltonian_3nu_liv( \
                                                h_vacuum_energy_independent,
                                                energy*1.e9,
                                                SXI12, SXI23, SXI13, DXICP,
                                                B1, B2, B3, LAMBDA),
                    baseline)  \
                for energy in energy_val]
        label_case = r'CPT-odd LIV'

    # Each element of prob: [Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt]
    prob_ee = [x[0] for x in prob]
    prob_em = [x[1] for x in prob]
    prob_et = [x[2] for x in prob]
    prob_me = [x[3] for x in prob]
    prob_mm = [x[4] for x in prob]
    prob_mt = [x[5] for x in prob]
    prob_te = [x[6] for x in prob]
    prob_tm = [x[7] for x in prob]
    prob_tt = [x[8] for x in prob]

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

    ax.set_xlabel(r'Neutrino energy $E$ [GeV]', fontsize=25)
    ax.set_ylabel(r'Three-neutrino probability', fontsize=25)

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

    ax.set_xlim([10.**log10_energy_min, 10.**log10_energy_max])
    ax.set_xscale('log')
    ax.set_ylim([0.0, 1.0])

    # Plot
    if (plot_prob_ee):
        ax.plot(energy_val, prob_ee, label=r'$P_{\nu_e \to \nu_e}$',
            color='C0', zorder=1)
    if (plot_prob_em):
        ax.plot(energy_val, prob_em, label=r'$P_{\nu_e \to \nu_\mu}$',
            color='C1', zorder=1)
    if (plot_prob_et):
        ax.plot(energy_val, prob_et, label=r'$P_{\nu_e \to \nu_\tau}$',
            color='C2', zorder=1)
    if (plot_prob_me):
        ax.plot(energy_val, prob_me, label=r'$P_{\nu_\mu \to \nu_e}$',
            color='C3', zorder=1)
    if (plot_prob_mm):
        ax.plot(energy_val, prob_mm, label=r'$P_{\nu_\mu \to \nu_\mu}$',
            color='C4', zorder=1)
    if (plot_prob_mt):
        ax.plot(energy_val, prob_mt, label=r'$P_{\nu_\mu \to \nu_\tau}$',
            color='C5', zorder=1)
    if (plot_prob_te):
        ax.plot(energy_val, prob_te, label=r'$P_{\nu_\tau \to \nu_e}$',
            color='C6', zorder=1)
    if (plot_prob_tm):
        ax.plot(energy_val, prob_tm, label=r'$P_{\nu_\tau \to \nu_\mu}$',
            color='C7', zorder=1)
    if (plot_prob_tt):
        ax.plot(energy_val, prob_tt, label=r'$P_{\nu_\tau \to \nu_\tau}$',
            color='C8', zorder=1)

    # Legend
    ax.legend(loc=legend_loc, frameon=False, ncol=legend_ncol)

    # Annotations
    ax.annotate( label_case, xy = (0.05, 0.86), \
        xycoords='axes fraction', color='k', fontsize=25,
        horizontalalignment='left', rotation=0, zorder=2 )
    ax.annotate( r'$\log_{10}(L/{\rm km}) = $' + \
        str(int(log10(baseline/CONV_KM_TO_INV_EV)*100.)/100.),
        xy = (0.05, 0.80), xycoords='axes fraction', color='k', fontsize=25,
        horizontalalignment='left', rotation=0, zorder=2 )

    pylab.savefig(output_path+output_filename+'.'+output_format,
        bbox_inches='tight', dpi=100)

    plt.close()

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

