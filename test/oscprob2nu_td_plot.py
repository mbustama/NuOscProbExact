# -*- coding: utf-8 -*-
r"""Routines to plot two-neutrino flavor-transition probabilities (t-dep. H).

This module contains contains routines to plot two-neutrino
oscillation probabilities vs. the neutrino baseline and energy for time-
dependent Hamiltonians.  These routines are used by run_testsuite.py to 
produce a suite of test plots.

Routine listings
----------------

    * plot_probability_2nu_td_vs_baseline - Plot probabilities vs. baseline
    * plot_probability_2nu_td_vs_energy - Plot probabilities vs. energy

Created: 2024/11/28 21:00
Last modified: 2024/11/29 17:26
"""


__version__ = "1.0"
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

import oscprob2nu
import hamiltonians2nu
import oscprob2nu_td
import hamiltonians2nu_td
import matter
from globaldefs import *


def plot_probability_2nu_td_vs_baseline(
                sector, energy=1.e-1,
                log10_l_min=0.0, log10_l_max=3.0, log10_l_npts=6000,
                plot_prob_ee=True, plot_prob_em=True, plot_prob_mm=False,
                output_filename='prob_td_vs_baseline', 
                output_format='pdf', output_path='../fig/', 
                legend_loc='center left', legend_ncol=1):
    r"""Generates and saves a plot of 2nu probabilities vs. baseline.

    Generates a plot of two-neutrino oscillation probabilities vs.
    baseline, for a fixed neutrino energy.  The probabilities to be
    plotted are turned on and off via the flags plot_prob_ee,
    plot_prob_em, etc.  (At least one of them must be True.)  The
    parameter 'case' selects between 'vacuum', 'matter', 'nsi', and
    'liv' (see below).  The plot is saved with the provided name and
    file format under the specified path.

    Parameters
    ----------
    sector : str
        Not optional.  Must be one of the following: '12' (for nu_e <-->
        nu_mu oscillations) of '23' (for nu_mu <--> nu_tau
        oscillations).
    energy : float, optional
        Neutrino energy [GeV].
    log10_l_min : float, optional
        Log10 of the minimum baseline [km].
    log10_l_max : float, optional
        Log10 of the maximum baseline [km].
    log10_l_npts : int, optional
        Number of baseline values at which to compute the probabilities.
    plot_prob_ee : bool, optional
        True to plot Pee (if sector == '12') or Pmm (if sector == '23),
        False otherwise.
    plot_prob_em : bool, optional
        True to plot Pem (if sector == '12') or Pmt (if sector == '23),
        False otherwise.
    plot_prob_mm : bool, optional
        True to plot Pmm (if sector == '12') or Ptt (if sector == '23),
        False otherwise.
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
    if (not plot_prob_ee) and (not plot_prob_em) and (not plot_prob_mm):
        quit()

    # Formatting
    mpl.rcParams['xtick.labelsize']=26
    mpl.rcParams['ytick.labelsize']=26
    mpl.rcParams['legend.fontsize']=26
    mpl.rcParams['legend.borderpad']=0.4
    mpl.rcParams['axes.labelpad']=10
    mpl.rcParams['ps.fonttype']=42
    mpl.rcParams['pdf.fonttype']=42

    # Open the plot
    fig = plt.figure(figsize=[18,9])
    ax = fig.add_subplot(1,1,1)

    # Baselines, L
    log10_l_val = np.linspace(log10_l_min, log10_l_max, log10_l_npts)
    l_val = [10.**x for x in log10_l_val]

    if sector == '12':
        sth = S12_NO_BF
        Dm2 = D21_NO_BF
        label_ee = r'$P_{\nu_e \to \nu_e}$'
        label_em = r'$P_{\nu_e \to \nu_\mu}$'
        label_mm = r'$P_{\nu_\mu \to \nu_\mu}$'
        color_ee = 'C0'
        color_em = 'C1'
        color_mm = 'C4'
    elif sector == '23':
        sth = S23_NO_BF
        Dm2 = D31_NO_BF
        label_ee = r'$P_{\nu_\mu \to \nu_\mu}$'
        label_em = r'$P_{\nu_\mu \to \nu_\tau}$'
        label_mm = r'$P_{\nu_\tau \to \nu_\tau}$'
        color_ee = 'C4'
        color_em = 'C5'
        color_mm = 'C8'

    for case in ['vacuum', 'matter_const', 'matter_exp_1', 'matter_exp_2']:

        if (case.lower() == 'vacuum'):

            # Define a potential that is only a function of position, l.
            # In this case, oscilllations are in vacuum, so this is only for
            # testing.
            h_vacuum_func = lambda l: np.multiply(1./energy/1.e9, 
                hamiltonians2nu_td.hamiltonian_2nu_vacuum_energy_independent_td(l,
                sth, Dm2))
            h_func = h_vacuum_func
            label_case = r'Vacuum ($t$-dep. calculation)'
            ls = '-'
            lc = 'k'

        elif (case.lower() == 'matter_const'):

            # Define a potential that is only a function of position, l.
            # In this case, the matter density is constant, so this is only
            # for testing.
            def VCC_func_const_density(l):
                if (sector == '12'):
                    density_matter_const = 10 # [g cm^{-3}]
                    label_case = r'Matter, constant $\rho_0 = 10$ g cm$^{-3}$'
                elif (sector == '23'):
                    density_matter_const = 100 # [g cm^{-3}]
                    label_case = r'Matter, constant $\rho_0 = 100$ g cm$^{-3}$'
                density_matter_func = lambda l: \
                    matter.density_matter_func_const(l, 
                        density_matter_const) # [g cm^{-3}]
                num_density_e_func = lambda l : \
                    matter.num_density_e_func(l,
                        density_matter_func, electron_fraction=0.5) # [eV^{-3}]
                return matter.VCC_func(l, num_density_e_func) # [eV]

            h_vacuum_energy_independent = \
                hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(sth,
                    Dm2)
            h_matter_func = lambda l: \
                hamiltonians2nu_td.hamiltonian_2nu_matter_td( \
                    h_vacuum_energy_independent, l, energy*1.e9, 
                    VCC_func_const_density)
            h_func = h_matter_func
            ls = '--'
            lc = 'C0'

        elif (case.lower() == 'matter_exp_1'):

            # Define a potential that is only a function of position, l
            def VCC_func_exp(l):
                if (sector == '12'):
                    density_matter_central = 10 # [g cm^{-3}]
                    l_scale = 500*CONV_KM_TO_INV_EV # [km]
                elif (sector == '23'):
                    density_matter_central = 100 # [g cm^{-3}]
                    l_scale = 100*CONV_KM_TO_INV_EV # [km]
                density_matter_func = lambda l: \
                    matter.density_matter_func_exp(l, 
                        density_matter_central, l_scale) # [g cm^{-3}]
                num_density_e_func = lambda l : \
                    matter.num_density_e_func(l,
                        density_matter_func, electron_fraction=0.5) # [eV^{-3}]
                return matter.VCC_func(l, num_density_e_func) # [eV]

            h_vacuum_energy_independent = \
                hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(sth,
                    Dm2)
            h_matter_func = lambda l: \
                hamiltonians2nu_td.hamiltonian_2nu_matter_td( \
                    h_vacuum_energy_independent, l, energy*1.e9, 
                    VCC_func_exp)
            h_func = h_matter_func
            label_case = r'Matter, $\rho = \rho_0~e^{-r/r_0}~(r_0 = 500~{\rm km})$'
            ls = '-.'
            lc = 'C1'

        elif (case.lower() == 'matter_exp_2'):

            # Define a potential that is only a function of position, l
            def VCC_func_exp(l):
                if (sector == '12'):
                    density_matter_central = 10 # [g cm^{-3}]
                    l_scale = 100*CONV_KM_TO_INV_EV # [km]
                elif (sector == '23'):
                    density_matter_central = 100 # [g cm^{-3}]
                    l_scale = 100*CONV_KM_TO_INV_EV # [km]
                density_matter_func = lambda l: \
                    matter.density_matter_func_exp(l, 
                        density_matter_central, l_scale) # [g cm^{-3}]
                num_density_e_func = lambda l : \
                    matter.num_density_e_func(l,
                        density_matter_func, electron_fraction=0.5) # [eV^{-3}]
                return matter.VCC_func(l, num_density_e_func) # [eV]

            h_vacuum_energy_independent = \
                hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(sth,
                    Dm2)
            h_matter_func = lambda l: \
                hamiltonians2nu_td.hamiltonian_2nu_matter_td( \
                    h_vacuum_energy_independent, l, energy*1.e9, 
                    VCC_func_exp)
            h_func = h_matter_func
            label_case = r'Matter, $\rho = \rho_0~e^{-r/r_0}~(r_0 = 100~{\rm km})$'
            ls = ':'
            lc = 'C2'

        # Each element of prob: [Pee, Pem, Pmm]
        prob = [oscprob2nu_td.probabilities_2nu_td(
            h_func, l_val[0]*CONV_KM_TO_INV_EV, l*CONV_KM_TO_INV_EV, 
            integration_method='quad', epsrel=1.e-8, epsabs=1.e-8) for l in l_val]
        prob_ee = [x[0] for x in prob]
        prob_em = [x[1] for x in prob]
        prob_mm = [x[3] for x in prob]

        # Plot
        if (plot_prob_ee):
            ax.plot(l_val, prob_ee, label=label_case,
                color=lc, ls=ls, zorder=1)
        if (plot_prob_em):
            ax.plot(l_val, prob_em, label=label_case,
                color=lc, ls=ls, zorder=1)
        if (plot_prob_mm):
            ax.plot(l_val, prob_mm, label=label_case,
                color=lc, ls=ls, zorder=1)

    ax.set_xlabel(r'Baseline $L$ [km]', fontsize=25)

    if (sector == '12'):
        if plot_prob_ee and (not plot_prob_em) and (not plot_prob_mm):
            y_label = r'Two-neutrino probability, $P_{\nu_e \to \nu_e}$'
        elif (not plot_prob_ee) and plot_prob_em and (not plot_prob_mm):
            y_label = r'Two-neutrino probability, $P_{\nu_e \to \nu_\mu}$'
        elif (not plot_prob_ee) and (not plot_prob_em) and plot_prob_mm:
            y_label = r'Two-neutrino probability, $P_{\nu_\mu \to \nu_\mu}$'
    elif (sector == '23'):
        if plot_prob_ee and (not plot_prob_em) and (not plot_prob_mm):
            y_label = r'Two-neutrino probability, $P_{\nu_\mu \to \nu_\mu}$'
        elif (not plot_prob_ee) and plot_prob_em and (not plot_prob_mm):
            y_label = r'Two-neutrino probability, $P_{\nu_\mu \to \nu_\tau}$'
        elif (not plot_prob_ee) and (not plot_prob_em) and plot_prob_mm:
            y_label = r'Two-neutrino probability, $P_{\nu_\tau \to \nu_\tau}$'
    else:
        y_label = r'Two-neutrino probability'
    ax.set_ylabel(y_label, fontsize=25)

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

    # Legend
    ax.legend(loc=legend_loc, frameon=True, ncol=legend_ncol,
        fontsize=20)

    # Annotations
    ax.annotate( r'$\log_{10}(E/{\rm GeV}) = $' + \
        str(int(log10(energy)*100.)/100.),
        xy = (0.98, 0.05), xycoords='axes fraction', color='k', fontsize=20,
        horizontalalignment='right', rotation=0, zorder=2 )

    pylab.savefig(output_path+output_filename+'.'+output_format,
        bbox_inches='tight', dpi=100)

    plt.close()

    return


def plot_oscillogram_earth_2nu_td(
                sector, prob_sel='ee',
                costhz_min=-1.0, costhz_max=0.0, costhz_npts=100,
                log10_Enu_min=0.0, log10_Enu_max=2.0, Enu_npts=100,
                output_filename='oscillogram_earth_2nu_td', 
                output_format='png', output_path='../fig/'):

    def VCC_func_prem(r):
        density_matter_func = lambda r: \
            matter.density_matter_func_prem(r) # [g cm^{-3}]
        num_density_e_func = lambda r : \
            matter.num_density_e_func(r, density_matter_func,
            electron_fraction=0.5) # [eV^{-3}]
        return matter.VCC_func(r, num_density_e_func) # [eV]

    if sector == '12':
        sth = S12_NO_BF
        Dm2 = D21_NO_BF
        if (prob_sel == 'ee'):
            prob_index = 0
            label = r'$P_{\nu_e \to \nu_e}$'
        elif (prob_sel == 'em'):
            prob_index = 1
            label = r'$P_{\nu_e \to \nu_\mu}$'
        elif (prob_sel == 'mm'):
            prob_index = 2
            label = r'$P_{\nu_\mu \to \nu_\mu}$'
    elif sector == '23':
        sth = S23_NO_BF
        Dm2 = D31_NO_BF
        if (prob_sel == 'mm'):
            prob_index = 0
            label = r'$P_{\nu_\mu \to \nu_\mu}$'
        elif (prob_sel == 'mt'):
            prob_index = 1
            label = r'$P_{\nu_\mu \to \nu_\tau}$'
        elif (prob_sel == 'tt'):
            prob_index = 2
            label = r'$P_{\nu_\tau \to \nu_\tau}$'

    # Cosine of zenith angle
    costhz_val = np.linspace(costhz_min, costhz_max, costhz_npts)

    # Baselines, L
    l_val = [matter.distance_traveled_inside_earth(costhz) \
        for costhz in costhz_val] # [km]

    # Neutrino energies, Enu
    Enu_val = np.logspace(log10_Enu_min, log10_Enu_max, Enu_npts) # [GeV]
    log10_Enu_val = np.log10(Enu_val)

    # Create the 2D array to store the probability
    prob_2d = np.zeros((costhz_npts, Enu_npts))

    # Compute vacuum, energy-independent Hamiltonian just once
    h_vacuum_energy_independent = \
        hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(sth, 
            Dm2)

    # Generate the probability for all combinations of costhz and Enu
    for costhz_index, costhz in enumerate(costhz_val):
        for Enu_index, Enu in enumerate(Enu_val):
            # print(costhz_index, Enu_index, costhz, 
            # l_val[costhz_index], Enu_val[Enu_index])
            print(costhz_index, Enu_index)
            def VCC_func_prem_wrapper(l):
                r = matter.earth_radial_distance_from_depth(costhz, 
                    l/CONV_KM_TO_INV_EV)
                return VCC_func_prem(r)
            # Hamiltonian function including matter effects
            h_func = lambda l: \
                hamiltonians2nu_td.hamiltonian_2nu_matter_td( \
                    h_vacuum_energy_independent, l, Enu*1.e9, 
                    VCC_func_prem_wrapper) # [eV]
            # Each element of prob, e.g., [Pee, Pem, Pmm]
            prob = oscprob2nu_td.probabilities_2nu_td(
                h_func, 0.0, l_val[costhz_index]*CONV_KM_TO_INV_EV, 
                integration_method='quad', epsrel=1.e-8, epsabs=1.e-8)
            prob_2d[Enu_index][costhz_index] = prob[prob_index]

    # Plot the contour plot
    mpl.rcParams['xtick.labelsize']=26
    mpl.rcParams['ytick.labelsize']=26
    mpl.rcParams['legend.fontsize']=26
    mpl.rcParams['legend.borderpad']=0.4
    mpl.rcParams['axes.labelpad']=10
    mpl.rcParams['ps.fonttype']=42
    mpl.rcParams['pdf.fonttype']=42

    fig = plt.figure(figsize=[9,9])
    ax = fig.add_subplot(1,1,1)

    ax.tick_params('both', length=10, width=2, which='major')
    ax.tick_params('both', length=5, width=1, which='minor')
    ax.tick_params(axis='both', which='major', pad=10, direction='in')
    ax.tick_params(axis='both', which='minor', pad=10, direction='in')
    ax.tick_params(axis='x', which='minor', bottom=True)
    ax.tick_params(axis='x', which='minor', top=True)
    ax.tick_params(axis='y', which='minor', left=True)
    ax.tick_params(axis='y', which='minor', right=True)
    ax.tick_params(bottom=True, top=True, left=True, right=True)
    
    ax.set_xlim([costhz_min, costhz_max])
    ax.set_ylim([log10_Enu_min, log10_Enu_max])

    xaxis_major_locator = mpl.ticker.MultipleLocator(0.2)
    ax.xaxis.set_major_locator(xaxis_major_locator)
    xaxis_minor_locator = mpl.ticker.MultipleLocator(0.02)
    ax.xaxis.set_minor_locator(xaxis_minor_locator)
    yaxis_major_locator = mpl.ticker.MultipleLocator(0.1)
    ax.yaxis.set_major_locator(yaxis_major_locator)
    yaxis_minor_locator = mpl.ticker.MultipleLocator(0.02)
    ax.yaxis.set_minor_locator(yaxis_minor_locator)

    ax.set_xlabel(r'Zenith angle, $\cos(\theta_z)$', fontsize=25)
    ax.set_ylabel(r'Neutrino energy, $\log_{10}(E_\nu/{\rm GeV})$', 
        fontsize=25)

    # print(prob_2d)

    # Plot
    cs = ax.contourf(costhz_val, log10_Enu_val, prob_2d,
        levels=120, cmap=mpl.cm.plasma)
    cbar = fig.colorbar(cs)
    cbar.ax.tick_params(labelsize=25) 
    cbar.set_label(label=r'Two-neutrino probability, '+label, fontsize=25)

    plt.savefig(output_path+output_filename+'.'+output_format,
        bbox_inches='tight', dpi=100)

    plt.close()

    return



# def plot_probability_2nu_td_vs_baseline(
#                 case, sector, energy=1.e-1,
#                 log10_l_min=0.0, log10_l_max=3.0, log10_l_npts=6000,
#                 plot_prob_ee=True, plot_prob_em=True, plot_prob_mm=False,
#                 output_filename='prob_td_vs_baseline', 
#                 output_format='pdf', output_path='../fig/', 
#                 legend_loc='center left', legend_ncol=1):
#     r"""Generates and saves a plot of 2nu probabilities vs. baseline.

#     Generates a plot of two-neutrino oscillation probabilities vs.
#     baseline, for a fixed neutrino energy.  The probabilities to be
#     plotted are turned on and off via the flags plot_prob_ee,
#     plot_prob_em, etc.  (At least one of them must be True.)  The
#     parameter 'case' selects between 'vacuum', 'matter', 'nsi', and
#     'liv' (see below).  The plot is saved with the provided name and
#     file format under the specified path.

#     Parameters
#     ----------
#     case : str
#         Not optional.  Must be one of the following: 'vacuum', 'matter',
#         'nsi', or 'liv'.  In each case, the probabilities are computed
#         using the default parameter values pulled from globaldefs.
#     sector : str
#         Not optional.  Must be one of the following: '12' (for nu_e <-->
#         nu_mu oscillations) of '23' (for nu_mu <--> nu_tau
#         oscillations).
#     energy : float, optional
#         Neutrino energy [GeV].
#     log10_l_min : float, optional
#         Log10 of the minimum baseline [km].
#     log10_l_max : float, optional
#         Log10 of the maximum baseline [km].
#     log10_l_npts : int, optional
#         Number of baseline values at which to compute the probabilities.
#     plot_prob_ee : bool, optional
#         True to plot Pee (if sector == '12') or Pmm (if sector == '23),
#         False otherwise.
#     plot_prob_em : bool, optional
#         True to plot Pem (if sector == '12') or Pmt (if sector == '23),
#         False otherwise.
#     plot_prob_mm : bool, optional
#         True to plot Pmm (if sector == '12') or Ptt (if sector == '23),
#         False otherwise.
#     output_filename : str, optional
#         File name of plot to save (without the file extension).
#     output_format : str, optional
#         File extension of the plot to save (e.g., 'pdf', 'png', 'jpg').
#     output_path : str, optional
#         File path where to save the plot.
#     legend_loc : str, optional
#         Location of the legend in the plot.  Must be one of the allowed
#         values of the plot routine of matplotlib.
#     legend_ncol : int, optional
#         Number of columns to include in the legend box.  Must be at
#         least 1.

#     Returns
#     -------
#     None
#         The plot is generated and saved.
#     """
#     if (not plot_prob_ee) and (not plot_prob_em) \
#         and (not plot_prob_me) and (not plot_prob_mm):
#         quit()

#     # Baselines, L
#     log10_l_val = np.linspace(log10_l_min, log10_l_max, log10_l_npts)
#     l_val =[10.**x for x in log10_l_val]

#     if sector == '12':
#         sth = S12_NO_BF
#         Dm2 = D21_NO_BF
#         label_ee = r'$P_{\nu_e \to \nu_e}$'
#         label_em = r'$P_{\nu_e \to \nu_\mu}$'
#         label_mm = r'$P_{\nu_\mu \to \nu_\mu}$'
#         color_ee = 'C0'
#         color_em = 'C1'
#         color_mm = 'C4'
#     elif sector == '23':
#         sth = S23_NO_BF
#         Dm2 = D31_NO_BF
#         label_ee = r'$P_{\nu_\mu \to \nu_\mu}$'
#         label_em = r'$P_{\nu_\mu \to \nu_\tau}$'
#         label_mm = r'$P_{\nu_\tau \to \nu_\tau}$'
#         color_ee = 'C4'
#         color_em = 'C5'
#         color_mm = 'C8'

#     h_vacuum_func = lambda l: np.multiply(1./energy/1.e9, 
#         hamiltonians2nu_td.hamiltonian_2nu_vacuum_energy_independent_td(l,
#         sth, Dm2))

#     if (case.lower() == 'vacuum'):

#         h_func = h_vacuum_func
#         label_case = r'Vacuum ($t$-dep. calc.)'

#     elif (case.lower() == 'matter'):

#         hamiltonian = hamiltonians2nu.hamiltonian_2nu_matter( \
#                                                 h_vacuum_energy_independent,
#                                                 energy*1.e9,
#                                                 VCC_EARTH_CRUST)
#         label_case = r'Matter (constant density, $t$-dep. calc.)'

#     elif (case.lower() == 'nsi'):

#         hamiltonian = hamiltonians2nu.hamiltonian_2nu_nsi( \
#                                                 h_vacuum_energy_independent,
#                                                 energy*1.e9,
#                                                 VCC_EARTH_CRUST,
#                                                 EPS_2)
#         label_case = r'NSI'

#     elif (case.lower() == 'liv'):

#         hamiltonian = hamiltonians2nu.hamiltonian_2nu_liv( \
#                                                 h_vacuum_energy_independent,
#                                                 energy*1.e9,
#                                                 SXI12,
#                                                 B1, B3, LAMBDA)
#         label_case = r'CPT-odd LIV'


#     # Each element of prob: [Pee, Pem, Pmm]
#     prob = [oscprob2nu_td.probabilities_2nu_td(
#         h_func, l_val[0]*CONV_KM_TO_INV_EV, l*CONV_KM_TO_INV_EV, 
#         integration_method='quad', epsrel=1.e-8, epsabs=1.e-8) for l in l_val]
#     prob_ee = [x[0] for x in prob]
#     prob_em = [x[1] for x in prob]
#     prob_mm = [x[3] for x in prob]

#     # Formatting
#     mpl.rcParams['xtick.labelsize']=26
#     mpl.rcParams['ytick.labelsize']=26
#     mpl.rcParams['legend.fontsize']=26
#     mpl.rcParams['legend.borderpad']=0.4
#     mpl.rcParams['axes.labelpad']=10
#     mpl.rcParams['ps.fonttype']=42
#     mpl.rcParams['pdf.fonttype']=42

#     fig = plt.figure(figsize=[9,9])
#     ax = fig.add_subplot(1,1,1)

#     ax.set_xlabel(r'Baseline $L$ [km]', fontsize=25)
#     ax.set_ylabel(r'Two-neutrino probability', fontsize=25)

#     yaxis_minor_locator = mpl.ticker.MultipleLocator(0.1)
#     ax.yaxis.set_minor_locator(yaxis_minor_locator)

#     ax.tick_params('both', length=10, width=2, which='major')
#     ax.tick_params('both', length=5, width=1, which='minor')
#     ax.tick_params(axis='both', which='major', pad=10, direction='in')
#     ax.tick_params(axis='both', which='minor', pad=10, direction='in')
#     ax.tick_params(axis='x', which='minor', bottom=True)
#     ax.tick_params(axis='x', which='minor', top=True)
#     ax.tick_params(axis='y', which='minor', left=True)
#     ax.tick_params(axis='y', which='minor', right=True)
#     ax.tick_params(bottom=True, top=True, left=True, right=True)

#     ax.set_xlim([10.**log10_l_min, 10.**log10_l_max])
#     ax.set_xscale('log')
#     ax.set_ylim([0.0, 1.0])

#     # Plot
#     if (plot_prob_ee):
#         ax.plot(l_val, prob_ee, label=label_ee,
#             color=color_ee, zorder=1)
#     if (plot_prob_em):
#         ax.plot(l_val, prob_em, label=label_em,
#             color=color_em, zorder=1)
#     if (plot_prob_mm):
#         ax.plot(l_val, prob_mm, label=label_mm,
#             color=color_mm, zorder=1)

#     # Legend
#     ax.legend(loc=legend_loc, frameon=False, ncol=legend_ncol)

#     # Annotations
#     ax.annotate( label_case, xy = (0.05, 0.86), \
#         xycoords='axes fraction', color='k', fontsize=25,
#         horizontalalignment='left', rotation=0, zorder=2 )
#     ax.annotate( r'$\log_{10}(E/{\rm GeV}) = $' + \
#         str(int(log10(energy)*100.)/100.),
#         xy = (0.05, 0.80), xycoords='axes fraction', color='k', fontsize=25,
#         horizontalalignment='left', rotation=0, zorder=2 )

#     pylab.savefig(output_path+output_filename+'.'+output_format,
#         bbox_inches='tight', dpi=100)

#     plt.close()

#     return


# def plot_probability_2nu_td_vs_energy(
#                 case, sector, baseline=1.3e3,
#                 log10_energy_min=-1.0, log10_energy_max=1.0,
#                 log10_energy_npts=200,
#                 plot_prob_ee=True, plot_prob_em=True, plot_prob_mm=False,
#                 output_filename='prob_vs_energy', output_format='pdf',
#                 output_path='../fig/', legend_loc='center right',
#                 legend_ncol=1):
#     r"""Generates and saves a plot of 2nu probabilities vs. energy.

#     Generates a plot of two-neutrino oscillation probabilities vs.
#     energy, for a fixed neutrino baseline.  The probabilities to be
#     plotted are turned on and off via the flags plot_prob_ee,
#     plot_prob_em, etc.  (At least one of them must be True.)  The
#     parameter 'case' selects between 'vacuum', 'matter', 'nsi', and
#     'liv' (see below).  The plot is saved with the provided name and
#     file format under the specified path.

#     Parameters
#     ----------
#     case : str
#         Not optional.  Must be one of the following: 'vacuum', 'matter',
#         'nsi', or 'liv'.  In each case, the probabilities are computed
#         using the default parameter values pulled from globaldefs.
#     sector : str
#         Not optional.  Must be one of the following: '12' (for nu_e <-->
#         nu_mu oscillations) of '23' (for nu_mu <--> nu_tau
#         oscillations).
#     baseline : float, optional
#         Neutrino baseline [km].
#     log10_energy_min : float, optional
#         Log10 of the minimum energy [GeV].
#     log10_energy_max : float, optional
#         Log10 of the maximum energy [GeV].
#     log10_energy_npts : int, optional
#         Number of energy values at which to compute the probabilities.
#     plot_prob_ee : bool, optional
#         True to plot Pee (if sector == '12') or Pmm (if sector == '23),
#         False otherwise.
#     plot_prob_em : bool, optional
#         True to plot Pem (if sector == '12') or Pmt (if sector == '23),
#         False otherwise.
#     plot_prob_mm : bool, optional
#         True to plot Pmm (if sector == '12') or Ptt (if sector == '23),
#         False otherwise.
#     output_filename : str, optional
#         File name of plot to save (without the file extension).
#     output_format : str, optional
#         File extension of the plot to save (e.g., 'pdf', 'png', 'jpg').
#     output_path : str, optional
#         File path where to save the plot.
#     legend_loc : str, optional
#         Location of the legend in the plot.  Must be one of the allowed
#         values of the plot routine of matplotlib.
#     legend_ncol : int, optional
#         Number of columns to include in the legend box.  Must be at
#         least 1.

#     Returns
#     -------
#     None
#         The plot is generated and saved.
#     """
#     if (not plot_prob_ee) and (not plot_prob_em) \
#         and (not plot_prob_me) and (not plot_prob_mm):
#         quit()

#     baseline = baseline*CONV_KM_TO_INV_EV # [eV^{-1}]

#     # Neutrino energies
#     log10_energy_val = np.linspace( log10_energy_min, log10_energy_max,
#                                     log10_energy_npts)
#     energy_val =[10.**x for x in log10_energy_val]

#     if sector == '12':
#         sth = S12_NO_BF
#         Dm2 = D21_NO_BF
#         label_ee = r'$P_{\nu_e \to \nu_e}$'
#         label_em = r'$P_{\nu_e \to \nu_\mu}$'
#         label_mm = r'$P_{\nu_\mu \to \nu_\mu}$'
#         color_ee = 'C0'
#         color_em = 'C1'
#         color_mm = 'C4'
#     elif sector == '23':
#         sth = S23_NO_BF
#         Dm2 = D31_NO_BF
#         label_ee = r'$P_{\nu_\mu \to \nu_\mu}$'
#         label_em = r'$P_{\nu_\mu \to \nu_\tau}$'
#         label_mm = r'$P_{\nu_\tau \to \nu_\tau}$'
#         color_ee = 'C4'
#         color_em = 'C5'
#         color_mm = 'C8'

#     h_vacuum_energy_independent_func = lambda l: \
#         hamiltonians2nu_td.hamiltonian_2nu_vacuum_energy_independent_td(l,
#             sth, Dm2)

#     if (case.lower() == 'vacuum'):

#         prob = [oscprob2nu_td.probabilities_2nu_td(lambda l: 
#             np.multiply(1./energy/1.e9, 
#                 hamiltonians2nu_td.hamiltonian_2nu_vacuum_energy_independent_td(l, sth, Dm2)), 
#             0, baseline, 
#             integration_method='quad', epsrel=1.e-8, epsabs=1.e-8) \
#         for energy in energy_val]
#         label_case = r'Vacuum ($t$-dep. calc.)'

#     elif (case.lower() == 'matter'):

#         prob = [oscprob2nu.probabilities_2nu( \
#                     hamiltonians2nu.hamiltonian_2nu_matter( \
#                                                 h_vacuum_energy_independent,
#                                                 energy*1.e9,
#                                                 VCC_EARTH_CRUST),
#                     baseline)  \
#                 for energy in energy_val]
#         label_case = r'Matter'

#     elif (case.lower() == 'nsi'):

#         prob = [oscprob2nu.probabilities_2nu( \
#                     hamiltonians2nu.hamiltonian_2nu_nsi( \
#                                                 h_vacuum_energy_independent,
#                                                 energy*1.e9,
#                                                 VCC_EARTH_CRUST,
#                                                 EPS_2),
#                     baseline)  \
#                 for energy in energy_val]
#         label_case = r'NSI'

#     elif (case.lower() == 'liv'):

#         prob = [oscprob2nu.probabilities_2nu( \
#                     hamiltonians2nu.hamiltonian_2nu_liv( \
#                                                 h_vacuum_energy_independent,
#                                                 energy*1.e9,
#                                                 SXI12,
#                                                 B1, B3, LAMBDA),
#                     baseline)  \
#                 for energy in energy_val]
#         label_case = r'CPT-odd LIV'

#     # Each element of prob: [Pee, Pem, Pmm]
#     prob_ee = [x[0] for x in prob]
#     prob_em = [x[1] for x in prob]
#     prob_mm = [x[3] for x in prob]

#     # Formatting
#     mpl.rcParams['xtick.labelsize']=26
#     mpl.rcParams['ytick.labelsize']=26
#     mpl.rcParams['legend.fontsize']=26
#     mpl.rcParams['legend.borderpad']=0.4
#     mpl.rcParams['axes.labelpad']=10
#     mpl.rcParams['ps.fonttype']=42
#     mpl.rcParams['pdf.fonttype']=42

#     fig = plt.figure(figsize=[9,9])
#     ax = fig.add_subplot(1,1,1)

#     ax.set_xlabel(r'Neutrino energy $E$ [GeV]', fontsize=25)
#     ax.set_ylabel(r'Two-neutrino probability', fontsize=25)

#     yaxis_minor_locator = mpl.ticker.MultipleLocator(0.1)
#     ax.yaxis.set_minor_locator(yaxis_minor_locator)

#     ax.tick_params('both', length=10, width=2, which='major')
#     ax.tick_params('both', length=5, width=1, which='minor')
#     ax.tick_params(axis='both', which='major', pad=10, direction='in')
#     ax.tick_params(axis='both', which='minor', pad=10, direction='in')
#     ax.tick_params(axis='x', which='minor', bottom=True)
#     ax.tick_params(axis='x', which='minor', top=True)
#     ax.tick_params(axis='y', which='minor', left=True)
#     ax.tick_params(axis='y', which='minor', right=True)
#     ax.tick_params(bottom=True, top=True, left=True, right=True)

#     ax.set_xlim([10.**log10_energy_min, 10.**log10_energy_max])
#     ax.set_xscale('log')
#     ax.set_ylim([0.0, 1.0])

#     # Plot
#     if (plot_prob_ee):
#         ax.plot(energy_val, prob_ee, label=label_ee,
#             color=color_ee, zorder=1)
#     if (plot_prob_em):
#         ax.plot(energy_val, prob_em, label=label_em,
#             color=color_em, zorder=1)
#     if (plot_prob_mm):
#         ax.plot(energy_val, prob_mm, label=label_mm,
#             color=color_mm, zorder=1)

#     # Legend
#     ax.legend(loc=legend_loc, frameon=False, ncol=legend_ncol)

#     # Annotations
#     ax.annotate( label_case, xy = (0.05, 0.86), \
#         xycoords='axes fraction', color='k', fontsize=25,
#         horizontalalignment='left', rotation=0, zorder=2 )
#     ax.annotate( r'$\log_{10}(L/{\rm km}) = $' + \
#         str(int(log10(baseline/CONV_KM_TO_INV_EV)*100.)/100.),
#         xy = (0.05, 0.80), xycoords='axes fraction', color='k', fontsize=25,
#         horizontalalignment='left', rotation=0, zorder=2 )

#     pylab.savefig(output_path+output_filename+'.'+output_format,
#         bbox_inches='tight', dpi=100)

#     plt.close()

#     return
