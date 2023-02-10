# -*- coding: utf-8 -*-
r"""Routines to plot two-neutrino flavor-transition probabilities.

This module contains contains routines to plot two-neutrino
oscillation probabilities vs. the neutrino baseline and energy.  These
routines are used by run_testsuite.py to produce a suite of test plots."""


from numpy import *
import numpy as np
from pylab import *
from matplotlib import *
import matplotlib as mpl
import oscprob2nu_slab

import sys
sys.path.append('../src')
sys.path.append('../test')

import hamiltonians2nu
from globaldefs import *


def plot_probability_2nu_vs_baseline(
                case, sector, energy=1.e-1,
                log10_l_min=0.0, log10_l_max=3.0, log10_l_npts=6000,
                plot_prob_ee=True, plot_prob_em=True, plot_prob_mm=False,
                output_filename='prob_vs_baseline', output_format='pdf',
                output_path='../fig/', legend_loc='center left',
                legend_ncol=1):
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
    case : str
        Not optional.  Must be one of the following: 'matterslabs', 'nsimatterslabs',
        or 'castlewall'.  In each case, the probabilities are computed using the default parameter 
        values pulled from globaldefs.'matterslabs' contains two distinct slabs of matter, each 
        occupying half the input length (see below); 'nsimatterslabs' contains two distinct slabs
        of matter with non-standard neutrino interactions, each occupying half the input length; 
        'castlewall' containts three slabs, one for the earth's crust (occupying 30% of input 
        length), one for earth's crust (next 40%), and one again for earth's crust (last 30%).
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
    if (not plot_prob_ee) and (not plot_prob_em) \
        and (not plot_prob_me) and (not plot_prob_mm):
        quit()

    # Baselines, L
    log10_l_val = np.linspace(log10_l_min, log10_l_max, log10_l_npts)
    l_val =[10.**x for x in log10_l_val]

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


    #Hamiltonians
        
    h_vacuum_energy_indep = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2)
    
    h_earthcrust = hamiltonians2nu.hamiltonian_2nu_matter(h_vacuum_energy_indep, energy*1.e9, VCC_EARTH_CRUST)

    h_earthcore = hamiltonians2nu.hamiltonian_2nu_matter(h_vacuum_energy_indep, energy*1.e9, 10*VCC_EARTH_CRUST)
        
    h_nsi_earthcrust = hamiltonians2nu.hamiltonian_2nu_nsi(h_vacuum_energy_indep, energy*1.e9, VCC_EARTH_CRUST,
                                                EPS_2)
    
    h_nsi_earthcore = hamiltonians2nu.hamiltonian_2nu_nsi(h_vacuum_energy_indep, energy*1.e9, 10*VCC_EARTH_CRUST,
                                                EPS_2)

    length = 10**log10_l_max - 10**log10_l_min

    if (case.lower() == 'matterslabs'):
        hamiltonian_matrices = [h_earthcrust, h_earthcore]
        slabs_initial = np.multiply(CONV_KM_TO_INV_EV,[10**log10_l_min, 10**log10_l_min + length/2])
        width_final_slab = CONV_KM_TO_INV_EV*length/2
        label_case = r'MatterSlabs'
        
    elif (case.lower() == 'nsimatterslabs'):
        hamiltonian_matrices = [h_nsi_earthcrust, h_nsi_earthcore]
        slabs_initial = np.multiply(CONV_KM_TO_INV_EV,[10**log10_l_min, 10**log10_l_min + length/2])
        width_final_slab = CONV_KM_TO_INV_EV*length/2
        label_case = r'NSIMatterSlabs'
        
    elif (case.lower() == 'castlewall'):
        hamiltonian_matrices = [h_earthcrust, h_earthcore, h_earthcrust]
        slabs_initial = np.multiply(CONV_KM_TO_INV_EV,[10**log10_l_min, 10**log10_l_min + length*3/10, 10**log10_l_min + length*7/10])
        width_final_slab = CONV_KM_TO_INV_EV*length*3/10
        label_case = r'Castlewall'
        

    # Each element of prob: [Pee, Pem, Pmm]
    prob = [oscprob2nu_slab.probabilities_2nu_slabs(   hamiltonian_matrices, slabs_initial, width_final_slab,
                                            l*CONV_KM_TO_INV_EV) \
            for l in l_val]
        
    for p in prob:
        for x in p:
            if not 0<=x<=1:
                print('error')
                
    prob_ee = [x[0] for x in prob]
    prob_em = [x[1] for x in prob]
    prob_mm = [x[3] for x in prob]

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
    ax.set_ylabel(r'Two-neutrino probability', fontsize=25)

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
        ax.plot(l_val, prob_ee, label=label_ee,
            color=color_ee, zorder=1)
    if (plot_prob_em):
        ax.plot(l_val, prob_em, label=label_em,
            color=color_em, zorder=1)
    if (plot_prob_mm):
        ax.plot(l_val, prob_mm, label=label_mm,
            color=color_mm, zorder=1)
        
    # Vertical lines
    for x in slabs_initial:
        x = 1/CONV_KM_TO_INV_EV*x
        ax.plot((x, x), (0, 1), color = '0.5', ls = '--')

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
        bbox_inches='tight', dpi=100)

    plt.close()

    return


def plot_probability_2nu_vs_energy(
                case, sector, baseline=1.e3,
                log10_energy_min=-1.0, log10_energy_max=1.0,
                log10_energy_npts=200,
                plot_prob_ee=True, plot_prob_em=True, plot_prob_mm=False,
                output_filename='prob_vs_energy', output_format='pdf',
                output_path='../fig/', legend_loc='center right',
                legend_ncol=1):
    r"""Generates and saves a plot of 2nu probabilities vs. energy.

    Generates a plot of two-neutrino oscillation probabilities vs.
    energy, for a fixed neutrino baseline.  The probabilities to be
    plotted are turned on and off via the flags plot_prob_ee,
    plot_prob_em, etc.  (At least one of them must be True.)  The
    parameter 'case' selects between 'vacuum', 'matter', 'nsi', and
    'liv' (see below).  The plot is saved with the provided name and
    file format under the specified path.

    Parameters
    ----------
    case : str
        Not optional.  Must be one of the following: 'matterslabs', 'nsimatterslabs',
        or 'castlewall'.  In each case, the probabilities are computed using the default parameter 
        values pulled from globaldefs.'matterslabs' contains two distinct slabs of matter, each 
        occupying half the input length (see below); 'nsimatterslabs' contains two distinct slabs
        of matter with non-standard neutrino interactions, each occupying half the input length; 
        'castlewall' containts three slabs, one for the earth's crust (occupying 30% of input 
        length), one for earth's crust (next 40%), and one again for earth's crust (last 30%).
    sector : str
        Not optional.  Must be one of the following: '12' (for nu_e <-->
        nu_mu oscillations) of '23' (for nu_mu <--> nu_tau
        oscillations).
    log10_l_min : float, optional
        Log10 of the minimum baseline [km].
    log10_l_max : float, optional
        Log10 of the maximum baseline [km].
    baseline : float, optional
        Neutrino baseline [km].
    log10_energy_min : float, optional
        Log10 of the minimum energy [GeV].
    log10_energy_max : float, optional
        Log10 of the maximum energy [GeV].
    log10_energy_npts : int, optional
        Number of energy values at which to compute the probabilities.
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
    if (not plot_prob_ee) and (not plot_prob_em) \
        and (not plot_prob_me) and (not plot_prob_mm):
        quit()
        
    baseline = baseline*CONV_KM_TO_INV_EV # [eV^{-1}]
    log10_l_min = 0
    log10_l_max = 3.1
    length = 10**log10_l_max - 10**log10_l_min

    # Neutrino energies
    log10_energy_val = np.linspace( log10_energy_min, log10_energy_max,
                                    log10_energy_npts)
    energy_val =[10.**x for x in log10_energy_val]

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

        
    h_vacuum_energy_indep = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(sth, Dm2)

    if (case.lower() == 'matterslabs'):
        
        slabs_initial = np.multiply(CONV_KM_TO_INV_EV,[10**log10_l_min, 10**log10_l_min + length/2])
        width_final_slab = CONV_KM_TO_INV_EV*length/2
        prob = [oscprob2nu_slab.probabilities_2nu_slabs([
                hamiltonians2nu.hamiltonian_2nu_matter(h_vacuum_energy_indep, e*1.e9, 
                                                       VCC_EARTH_CRUST),
                hamiltonians2nu.hamiltonian_2nu_matter(h_vacuum_energy_indep, e*1.e9, 
                                                       10*VCC_EARTH_CRUST)],
             slabs_initial, width_final_slab,
             baseline) for e in energy_val]
        
        label_case = r'MatterSlabs'
    
    if (case.lower() == 'nsimatterslabs'):
        
        slabs_initial = np.multiply(CONV_KM_TO_INV_EV,[10**log10_l_min, 10**log10_l_min + length/2])
        width_final_slab = CONV_KM_TO_INV_EV*length/2
        prob = [oscprob2nu_slab.probabilities_2nu_slabs([
                hamiltonians2nu.hamiltonian_2nu_nsi(h_vacuum_energy_indep, e*1.e9, 
                                                        VCC_EARTH_CRUST, EPS_2),
                hamiltonians2nu.hamiltonian_2nu_nsi(h_vacuum_energy_indep, e*1.e9,
                                                        10*VCC_EARTH_CRUST, EPS_2)],
             slabs_initial, width_final_slab,
             baseline) for e in energy_val]
        
        label_case = r'NSIMatterSlabs'
    
    if (case.lower() == 'castlewall'):
        
        slabs_initial = np.multiply(CONV_KM_TO_INV_EV,[10**log10_l_min, 10**log10_l_min + length*3/10,
                                                       10**log10_l_min + length*7/10])
        width_final_slab = CONV_KM_TO_INV_EV*length*3/10
        prob = [oscprob2nu_slab.probabilities_2nu_slabs([
                hamiltonians2nu.hamiltonian_2nu_matter(h_vacuum_energy_indep, e*1.e9,
                                                        VCC_EARTH_CRUST),
                hamiltonians2nu.hamiltonian_2nu_matter(h_vacuum_energy_indep, e*1.e9,
                                                        10*VCC_EARTH_CRUST),
                hamiltonians2nu.hamiltonian_2nu_matter(h_vacuum_energy_indep, e*1.e9,
                                                        10*VCC_EARTH_CRUST)],
            slabs_initial, width_final_slab,
            baseline) for e in energy_val]
        
        label_case = r'Castlewall'
        
    for p in prob:
        for x in p:
            if not 0<=x<=1:
                print('error')

    # Each element of prob: [Pee, Pem, Pmm]
    prob_ee = [x[0] for x in prob]
    prob_em = [x[1] for x in prob]
    prob_mm = [x[3] for x in prob]

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
    ax.set_ylabel(r'Two-neutrino probability', fontsize=25)

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
        ax.plot(energy_val, prob_ee, label=label_ee,
            color=color_ee, zorder=1)
    if (plot_prob_em):
        ax.plot(energy_val, prob_em, label=label_em,
            color=color_em, zorder=1)
    if (plot_prob_mm):
        ax.plot(energy_val, prob_mm, label=label_mm,
            color=color_mm, zorder=1)

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
