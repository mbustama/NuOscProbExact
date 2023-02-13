# -*- coding: utf-8 -*-
r"""Produce a suite of test plots of probabilities.

Running this module ('python run_testsuite.py') creates a number of
test plots of the two- and three-neutrino oscillation probabilities
vs. baseline and vs. energy.  Also generates the plot included in the
paper.

Created: 2019/04/22 16:23 (Mauricio Bustamante)
Last modified: 2023/02/11 15:17 (Mauricio Bustamante)
"""


from __future__ import print_function

__version__ = "1.1"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


from numpy import *

import sys
import os
import argparse
sys.path.append('./src')
sys.path.append('./test')

import oscprob2nu_plot
import oscprob3nu_plot
import oscprob2nu_plotpaper
import oscprob3nu_plotpaper
import oscprob2nu_slab_plot
import oscprob3nu_slab_plot


def plot_testsuite_probability_2nu(output_format='png'):

    os.makedirs('./fig/prob_2nu/vs_baseline', exist_ok=True)
    os.makedirs('./fig/prob_2nu/vs_energy', exist_ok=True)

    message_print = {
        '12': '    Pee, Pem... ',
        '23': '    Pmm, Pmt... '
    }

    filename_end = {
        '12': 'ee_em',
        '23': 'mm_mt'
    }

    print('Generating plots of 2nu probability vs. baseline:')
    
    for case in ['vacuum', 'matter', 'nsi', 'liv']:
        print('  Case: '+case)
        for prob in ['12', '23']:
            print(message_print[prob], end='')
            oscprob2nu_plot.plot_probability_2nu_vs_baseline(
                        case, prob, energy=1.e-2,
                        log10_l_min=0.0, log10_l_max=3.0, log10_l_npts=3000,
                        plot_prob_ee=True, plot_prob_em=True, plot_prob_mm=False,
                        output_filename='prob_2nu_'+case+'_vs_baseline_'
                            +filename_end[prob],
                        output_format=output_format, 
                        output_path='./fig/prob_2nu/vs_baseline/',
                        legend_loc='center left', legend_ncol=1)
            print('Done')
        print('  Done')

    print()

    print('Generating plots of 2nu probability vs. energy:')

    for case in ['vacuum', 'matter', 'nsi', 'liv']:
        print('  Case: '+case)
        for prob in ['12', '23']:
            print(message_print[prob], end='')
            oscprob2nu_plot.plot_probability_2nu_vs_energy(
                        case, prob, baseline=1.e3,
                        log10_energy_min=-1.0, log10_energy_max=1.0,
                        log10_energy_npts=600,
                        plot_prob_ee=True, plot_prob_em=True, plot_prob_mm=False,
                        output_filename='prob_2nu_'+case+'_vs_energy_'
                            +filename_end[prob],
                        output_format=output_format, 
                        output_path='./fig/prob_2nu/vs_energy/',
                        legend_loc='center right', legend_ncol=1)
            print('Done')
        print('  Done')

    print()

    return


def plot_testsuite_probability_3nu(output_format='png'):

    os.makedirs('./fig/prob_3nu/vs_baseline', exist_ok=True)
    os.makedirs('./fig/prob_3nu/vs_energy', exist_ok=True)

    message_print = {
        'e': '    Pee, Pem, Pet... ',
        'mu': '    Pme, Pmm, Pmt... ',
        'tau': '    Pte, Ptm, Ptt... '
    }

    filename_end = {
        'e': 'ee_em_et',
        'mu': 'me_mm_mt',
        'tau': 'te_tm_tt'
    }

    prob_flags = {
        'e': [True, True, True, False, False, False, False, False, False],
        'mu': [False, False, False, True, True, True, False, False, False],
        'tau': [False, False, False, False, False, False, True, True, True]
    }

    print('Generating plots of 3nu probability vs. baseline:')

    for case in ['vacuum', 'matter', 'nsi', 'liv']:
        print('  Case: '+case)
        for prob in ['e', 'mu', 'tau']:
            print(message_print[prob], end='')
            plot_prob_ee, plot_prob_em, plot_prob_et, \
                plot_prob_me, plot_prob_mm, plot_prob_mt, \
                plot_prob_te, plot_prob_tm, plot_prob_tt = prob_flags[prob]
            # (Pee, Pem, Pet), (Pme, Pmm, Pmt), (Pte, Ptm, Ptt)
            oscprob3nu_plot.plot_probability_3nu_vs_baseline(
                        case, energy=1.e-2,
                        log10_l_min=0.0, log10_l_max=3.0, log10_l_npts=1000,
                        plot_prob_ee=plot_prob_ee, 
                        plot_prob_em=plot_prob_em, 
                        plot_prob_et=plot_prob_et,
                        plot_prob_me=plot_prob_me, 
                        plot_prob_mm=plot_prob_mm, 
                        plot_prob_mt=plot_prob_mt,
                        plot_prob_te=plot_prob_te, 
                        plot_prob_tm=plot_prob_tm, 
                        plot_prob_tt=plot_prob_tt,
                        output_filename='prob_3nu_'+case+'_vs_baseline_'
                            +filename_end[prob],
                        output_format=output_format, 
                        output_path='./fig/prob_3nu/vs_baseline/',
                        legend_loc='center left', legend_ncol=1)
            print('Done')
        print('  Done')

    print()

    print('Generating plots of 3nu probability vs. energy:')

    for case in ['vacuum', 'matter', 'nsi', 'liv']:
        print('  Case: '+case)
        for prob in ['e', 'mu', 'tau']:
            print(message_print[prob], end='')
            plot_prob_ee, plot_prob_em, plot_prob_et, \
                plot_prob_me, plot_prob_mm, plot_prob_mt, \
                plot_prob_te, plot_prob_tm, plot_prob_tt = prob_flags[prob]
            # (Pee, Pem, Pet), (Pme, Pmm, Pmt), (Pte, Ptm, Ptt)
            oscprob3nu_plot.plot_probability_3nu_vs_energy(
                        case, baseline=1.e3,
                        log10_energy_min=-1.0, log10_energy_max=1.0,
                        log10_energy_npts=200,
                        plot_prob_ee=plot_prob_ee, 
                        plot_prob_em=plot_prob_em, 
                        plot_prob_et=plot_prob_et,
                        plot_prob_me=plot_prob_me, 
                        plot_prob_mm=plot_prob_mm, 
                        plot_prob_mt=plot_prob_mt,
                        plot_prob_te=plot_prob_te, 
                        plot_prob_tm=plot_prob_tm, 
                        plot_prob_tt=plot_prob_tt,
                        output_filename='prob_3nu_'+case+'_vs_energy_'
                            +filename_end[prob],
                        output_format=output_format, 
                        output_path='./fig/prob_3nu/vs_energy/',
                        legend_loc='center right', legend_ncol=1)
            print('Done')
        print('  Done')

    print()

    return


def plot_testsuite_plots_paper(output_format='png'):

    os.makedirs('./fig/paper', exist_ok=True)

    print('Generating 2nu plot in paper... ', end='')
    oscprob2nu_plotpaper.plot_probability_2nu_vs_energy_compare( \
        output_format=output_format, output_path='./fig/paper/')
    print('Done')

    print()

    print('Generating 3nu plot in paper... ', end='')
    oscprob3nu_plotpaper.plot_probability_3nu_vs_energy_compare( \
        output_format=output_format, output_path='./fig/paper/')
    print('Done')

    print()

    return


def plot_testsuite_probability_2nu_slabs(output_format='png'):

    # Need to split the directory creation like this to avoid problems with long
    # paths
    os.makedirs('./fig/prob_2nu_slabs', exist_ok=True)
    os.makedirs('./fig/prob_2nu_slabs/vs_baseline', exist_ok=True)
    os.makedirs('./fig/prob_2nu_slabs/vs_energy', exist_ok=True)

    message_print = {
        '12': '    Pee, Pem... ',
        '23': '    Pmm, Pmt... '
    }

    filename_end = {
        '12': 'ee_em',
        '23': 'mm_mt'
    }

    print('Generating plots of 2nu probability (slabs) vs. baseline:')

    for case in ['slabs_matter', 'slabs_nsi', 'slabs_castle_wall']:
        print('  Case: '+case)
        for prob in ['12', '23']:
            print(message_print[prob], end='')
            oscprob2nu_slab_plot.plot_probability_2nu_vs_baseline(
                        case, prob, energy=1.e-2,
                        log10_l_min=0.0, log10_l_max=3.0, log10_l_npts=3000,
                        plot_prob_ee=True, plot_prob_em=True, plot_prob_mm=False,
                        output_filename='prob_2nu_'+case+'_vs_baseline_'
                            +filename_end[prob],
                        output_format=output_format, 
                        output_path='./fig/prob_2nu_slabs/vs_baseline/',
                        legend_loc='center left', legend_ncol=1)
            print('Done')
        print('  Done')

    print()

    print('Generating plots of 2nu probability (slabs) vs. energy:')

    for case in ['slabs_matter', 'slabs_nsi', 'slabs_castle_wall']:
        print('  Case: '+case)
        for prob in ['12', '23']:
            print(message_print[prob], end='')
            oscprob2nu_slab_plot.plot_probability_2nu_vs_energy(
                        case, prob, baseline=1.e3,
                        log10_energy_min=-1.0, log10_energy_max=1.0,
                        log10_energy_npts=600,
                        plot_prob_ee=True, plot_prob_em=True, plot_prob_mm=False,
                        output_filename='prob_2nu_'+case+'_vs_energy_'
                            +filename_end[prob],
                        output_format=output_format, 
                        output_path='./fig/prob_2nu_slabs/vs_energy/',
                        legend_loc='center right', legend_ncol=1)
            print('Done')
        print('  Done')

    print()

    return


def plot_testsuite_probability_3nu_slabs(output_format='png'):

    # Need to split the directory creation like this to avoid problems with long
    # paths
    os.makedirs('./fig/prob_3nu_slabs/', exist_ok=True)
    os.makedirs('./fig/prob_3nu_slabs/vs_baseline', exist_ok=True)
    os.makedirs('./fig/prob_3nu_slabs/vs_energy', exist_ok=True)

    message_print = {
        'e': '    Pee, Pem, Pet... ',
        'mu': '    Pme, Pmm, Pmt... ',
        'tau': '    Pte, Ptm, Ptt... '
    }

    filename_end = {
        'e': 'ee_em_et',
        'mu': 'me_mm_mt',
        'tau': 'te_tm_tt'
    }

    prob_flags = {
        'e': [True, True, True, False, False, False, False, False, False],
        'mu': [False, False, False, True, True, True, False, False, False],
        'tau': [False, False, False, False, False, False, True, True, True]
    }

    print('Generating plots of 3nu probability (slabs) vs. baseline:')

    for case in ['slabs_matter', 'slabs_nsi', 'slabs_castle_wall']:
        print('  Case: '+case)
        for prob in ['e', 'mu', 'tau']:
            print(message_print[prob], end='')
            plot_prob_ee, plot_prob_em, plot_prob_et, \
                plot_prob_me, plot_prob_mm, plot_prob_mt, \
                plot_prob_te, plot_prob_tm, plot_prob_tt = prob_flags[prob]
            # (Pee, Pem, Pet), (Pme, Pmm, Pmt), (Pte, Ptm, Ptt)
            oscprob3nu_slab_plot.plot_probability_3nu_vs_baseline(
                        case, energy=1.e-2  ,
                        log10_l_min=0.0, log10_l_max=3.0, log10_l_npts=1000,
                        plot_prob_ee=plot_prob_ee, 
                        plot_prob_em=plot_prob_em, 
                        plot_prob_et=plot_prob_et,
                        plot_prob_me=plot_prob_me, 
                        plot_prob_mm=plot_prob_mm, 
                        plot_prob_mt=plot_prob_mt,
                        plot_prob_te=plot_prob_te, 
                        plot_prob_tm=plot_prob_tm, 
                        plot_prob_tt=plot_prob_tt,
                        output_filename='prob_3nu_'+case+'_vs_baseline_'
                             +filename_end[prob],
                        output_format=output_format, 
                        output_path='./fig/prob_3nu_slabs/vs_baseline/',
                        legend_loc='center left', legend_ncol=1)
            print('Done')
        print('  Done')

    print()

    print('Generating plots of 3nu probability (slabs) vs. energy:')

    for case in ['slabs_matter', 'slabs_nsi', 'slabs_castle_wall']:
        print('  Case: '+case)
        for prob in ['e', 'mu', 'tau']:
            print(message_print[prob], end='')
            plot_prob_ee, plot_prob_em, plot_prob_et, \
                plot_prob_me, plot_prob_mm, plot_prob_mt, \
                plot_prob_te, plot_prob_tm, plot_prob_tt = prob_flags[prob]
            # (Pee, Pem, Pet), (Pme, Pmm, Pmt), (Pte, Ptm, Ptt)
            oscprob3nu_slab_plot.plot_probability_3nu_vs_energy(
                        case, baseline=1.e3,
                        log10_energy_min=-1.0, log10_energy_max=1.0,
                        log10_energy_npts=200,
                        plot_prob_ee=plot_prob_ee, 
                        plot_prob_em=plot_prob_em, 
                        plot_prob_et=plot_prob_et,
                        plot_prob_me=plot_prob_me, 
                        plot_prob_mm=plot_prob_mm, 
                        plot_prob_mt=plot_prob_mt,
                        plot_prob_te=plot_prob_te, 
                        plot_prob_tm=plot_prob_tm, 
                        plot_prob_tt=plot_prob_tt,
                        output_filename='prob_3nu_'+case+'_vs_energy_'+
                            filename_end[prob],
                        output_format=output_format, 
                        output_path='./fig/prob_3nu_slabs/vs_energy/',
                        legend_loc='center right', legend_ncol=1)
            print('Done')
        print('  Done')

    print()

    return


def main():

    print('NuOscProbExact: Running test suite (plots will be stored in ./fig)')
    print()
    plot_testsuite_probability_2nu()
    plot_testsuite_probability_3nu()
    plot_testsuite_plots_paper()
    plot_testsuite_probability_2nu_slabs()
    plot_testsuite_probability_3nu_slabs()
    
    return


if __name__ == '__main__':
    sys.exit(main())