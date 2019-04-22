# -*- coding: utf-8 -*-
r"""Produce a suite of test plots of probabilities.

Running this module ('python run_testsuite.py') creates a number of
test plots of the two- and three-neutrino oscillation probabilities
vs. baseline and vs. energy.  Also generates the plot included in the
paper.

Created: 2019/04/22 16:23
Last modified: 2019/04/22 20:21
"""

from __future__ import print_function

__version__ = "0.1"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


from numpy import *

import sys
sys.path.append('./src')
sys.path.append('./test')

import oscprob2nu_plot
import oscprob3nu_plot
import oscprob3nu_plotpaper


print('NuOscProbExact: Running test suite (plots will be stored inside ./fig)')
print()

print('Generating plots of 2nu probability vs. baseline:')

for case in ['vacuum', 'matter', 'nsi', 'liv']:

    print('  Case: '+case)

    print('    Pee, Pem... ', end='')
    oscprob2nu_plot.plot_probability_2nu_vs_baseline(
                case, '12', energy=1.e-2,
                log10_l_min=0.0, log10_l_max=3.0, log10_l_npts=3000,
                plot_prob_ee=True, plot_prob_em=True, plot_prob_mm=False,
                output_filename='prob_2nu_'+case+'_vs_baseline_ee_em',
                output_format='png', output_path='./fig/',
                legend_loc='center left', legend_ncol=1)
    print('Done')

    print('    Pmm, Pmt... ', end='')
    oscprob2nu_plot.plot_probability_2nu_vs_baseline(
                case, '23', energy=1.e-2,
                log10_l_min=0.0, log10_l_max=3.0, log10_l_npts=3000,
                plot_prob_ee=True, plot_prob_em=True, plot_prob_mm=False,
                output_filename='prob_2nu_'+case+'_vs_baseline_mm_mt',
                output_format='png', output_path='./fig/',
                legend_loc='center left', legend_ncol=1)
    print('Done')

    print('  Done')

print()

print('Generating plots of 2nu probability vs. energy:')

for case in ['vacuum', 'matter', 'nsi', 'liv']:

    print('  Case: '+case)

    print('    Pee, Pem... ', end='')
    oscprob2nu_plot.plot_probability_2nu_vs_energy(
                case, '12', baseline=1.e3,
                log10_energy_min=-1.0, log10_energy_max=1.0,
                log10_energy_npts=600,
                plot_prob_ee=True, plot_prob_em=True, plot_prob_mm=False,
                output_filename='prob_2nu_'+case+'_vs_energy_ee_em',
                output_format='png', output_path='./fig/',
                legend_loc='center right', legend_ncol=1)
    print('Done')

    print('    Pee, Pem... ', end='')
    oscprob2nu_plot.plot_probability_2nu_vs_energy(
                case, '23', baseline=1.e3,
                log10_energy_min=-1.0, log10_energy_max=1.0,
                log10_energy_npts=600,
                plot_prob_ee=True, plot_prob_em=True, plot_prob_mm=False,
                output_filename='prob_2nu_'+case+'_vs_energy_mm_mt',
                output_format='png', output_path='./fig/',
                legend_loc='center right', legend_ncol=1)
    print('Done')


    print('  Done')

print()

print('Generating plots of 3nu probability vs. baseline:')

for case in ['vacuum', 'matter', 'nsi', 'liv']:

    print('  Case: '+case)

    print('    Pee, Pem, Pet... ', end='')
    # Pee, Pem, Pet
    oscprob3nu_plot.plot_probability_3nu_vs_baseline(
                case, energy=1.e-2,
                log10_l_min=0.0, log10_l_max=3.0, log10_l_npts=1000,
                plot_prob_ee=True, plot_prob_em=True, plot_prob_et=True,
                plot_prob_me=False, plot_prob_mm=False, plot_prob_mt=False,
                plot_prob_te=False, plot_prob_tm=False, plot_prob_tt=False,
                output_filename='prob_3nu_'+case+'_vs_baseline_ee_em_et',
                output_format='png', output_path='./fig/',
                legend_loc='center left', legend_ncol=1)
    print('Done')

    print('    Pme, Pmm, Pmt... ', end='')
    # Pme, Pmm, Pmt
    oscprob3nu_plot.plot_probability_3nu_vs_baseline(
                case, energy=1.e-2,
                log10_l_min=0.0, log10_l_max=3.0, log10_l_npts=1000,
                plot_prob_ee=False, plot_prob_em=False, plot_prob_et=False,
                plot_prob_me=True, plot_prob_mm=True, plot_prob_mt=True,
                plot_prob_te=False, plot_prob_tm=False, plot_prob_tt=False,
                output_filename='prob_3nu_'+case+'_vs_baseline_me_mm_mt',
                output_format='png', output_path='./fig/',
                legend_loc='center left', legend_ncol=1)
    print('Done')

    print('    Pte, Ptm, Ptt... ', end='')
    # Pte, Ptm, Ptt
    oscprob3nu_plot.plot_probability_3nu_vs_baseline(
                case, energy=1.e-2,
                log10_l_min=0.0, log10_l_max=3.0, log10_l_npts=1000,
                plot_prob_ee=False, plot_prob_em=False, plot_prob_et=False,
                plot_prob_me=False, plot_prob_mm=False, plot_prob_mt=False,
                plot_prob_te=True, plot_prob_tm=True, plot_prob_tt=True,
                output_filename='prob_3nu_'+case+'_vs_baseline_te_tm_tt',
                output_format='png', output_path='./fig/',
                legend_loc='center left', legend_ncol=1)
    print('Done')
    print('  Done')

print()

print('Generating plots of 3nu probability vs. energy:')

for case in ['vacuum', 'matter', 'nsi', 'liv']:

    print('  Case: '+case)

    print('    Pee, Pem, Pet... ', end='')
    # Pee, Pem, Pet
    oscprob3nu_plot.plot_probability_3nu_vs_energy(
                case, baseline=1.e3,
                log10_energy_min=-1.0, log10_energy_max=1.0,
                log10_energy_npts=200,
                plot_prob_ee=True, plot_prob_em=True, plot_prob_et=True,
                plot_prob_me=False, plot_prob_mm=False, plot_prob_mt=False,
                plot_prob_te=False, plot_prob_tm=False, plot_prob_tt=False,
                output_filename='prob_3nu_'+case+'_vs_energy_ee_em_et',
                output_format='png', output_path='./fig/',
                legend_loc='center right', legend_ncol=1)
    print('Done')

    print('    Pme, Pmm, Pmt... ', end='')
    # Pme, Pmm, Pmt
    oscprob3nu_plot.plot_probability_3nu_vs_energy(
                case, baseline=1.e3,
                log10_energy_min=-1.0, log10_energy_max=1.0,
                log10_energy_npts=200,
                plot_prob_ee=False, plot_prob_em=False, plot_prob_et=False,
                plot_prob_me=True, plot_prob_mm=True, plot_prob_mt=True,
                plot_prob_te=False, plot_prob_tm=False, plot_prob_tt=False,
                output_filename='prob_3nu_'+case+'_vs_energy_me_mm_mt',
                output_format='png', output_path='./fig/',
                legend_loc='center right', legend_ncol=1)
    print('Done')

    print('    Pte, Ptm, Ptt... ', end='')
    # Pte, Ptm, Ptt
    oscprob3nu_plot.plot_probability_3nu_vs_energy(
                case, baseline=1.e3,
                log10_energy_min=-1.0, log10_energy_max=1.0,
                log10_energy_npts=200,
                plot_prob_ee=False, plot_prob_em=False, plot_prob_et=False,
                plot_prob_me=False, plot_prob_mm=False, plot_prob_mt=False,
                plot_prob_te=True, plot_prob_tm=True, plot_prob_tt=True,
                output_filename='prob_3nu_'+case+'_vs_energy_te_tm_tt',
                output_format='png', output_path='./fig/',
                legend_loc='center right', legend_ncol=1)
    print('Done')
    print('  Done')

print()

print('Generating plot in paper... ', end='')
oscprob3nu_plotpaper.plot_probability_vs_energy_compare(output_format='png',
  output_path='./fig/')
print('Done')




