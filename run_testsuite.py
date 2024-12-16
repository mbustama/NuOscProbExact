# -*- coding: utf-8 -*-
r"""Produce a suite of test plots of probabilities.

Running this module ('python run_testsuite.py') creates a number of
test plots of the two- and three-neutrino oscillation probabilities
vs. baseline and vs. energy.  Also generates the plot included in the
paper.

Created: 2019/04/22 16:23
Last modified: 2024/12/03 18:39
"""


from __future__ import print_function

__version__ = "1.0"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


from numpy import *

import sys
sys.path.append('./src')
sys.path.append('./test')

import oscprob2nu_plot
import oscprob2nu_td_plot
import oscprob3nu_plot
import oscprob3nu_td_plot
import oscprob2nu_plotpaper
import oscprob3nu_plotpaper
from globaldefs import *


print('NuOscProbExact: Running test suite (plots will be stored inside ./fig)')
print()

# print('Generating plots of 2nu probability vs. baseline:')

# for case in ['vacuum', 'matter', 'nsi', 'liv']:

#     print('  Case: '+case)

#     print('    Pee, Pem... ', end='')
#     oscprob2nu_plot.plot_probability_2nu_vs_baseline(
#                 case, '12', energy=1.e-2,
#                 log10_l_min=0.0, log10_l_max=3.0, log10_l_npts=3000,
#                 plot_prob_ee=True, plot_prob_em=True, plot_prob_mm=False,
#                 output_filename='prob_2nu_'+case+'_vs_baseline_ee_em',
#                 output_format='png', output_path='./fig/',
#                 legend_loc='center left', legend_ncol=1)
#     print('Done')

#     print('    Pmm, Pmt... ', end='')
#     oscprob2nu_plot.plot_probability_2nu_vs_baseline(
#                 case, '23', energy=1.e-2,
#                 log10_l_min=0.0, log10_l_max=3.0, log10_l_npts=3000,
#                 plot_prob_ee=True, plot_prob_em=True, plot_prob_mm=False,
#                 output_filename='prob_2nu_'+case+'_vs_baseline_mm_mt',
#                 output_format='png', output_path='./fig/',
#                 legend_loc='center left', legend_ncol=1)
#     print('Done')

#     print('  Done')

# print()


# print('Generating plots of 2nu probability vs. baseline (time-dependent Hamiltonian):')

# # --- Uncomment this when validating vacuum calculations (see oscprob2nu_td_plot.y)
# # print('    Pee (time-dependent Hamiltonian calculation)... ', end='')
# # oscprob2nu_td_plot.plot_probability_2nu_td_vs_baseline(
# #             '12', energy=1.e-2,
# #             log10_l_min=log10(1.e6), log10_l_max=log10(1.e7), log10_l_npts=10000,
# #             plot_prob_ee=True, plot_prob_em=False, plot_prob_mm=False,
# #             output_filename='prob_2nu_td_vacuum_vs_baseline_ee_validation',
# #             output_format='png', output_path='./fig/',
# #             legend_loc='lower left', legend_ncol=1)
# # print('Done')

# print('    Pee (time-dependent Hamiltonian calculation)... ', end='')
# oscprob2nu_td_plot.plot_probability_2nu_td_vs_baseline(
#             '12', energy=1.e-2,
#             log10_l_min=log10(5.e1), log10_l_max=3.0, log10_l_npts=3000,
#             plot_prob_ee=True, plot_prob_em=False, plot_prob_mm=False,
#             output_filename='prob_2nu_td_vacuum_matter_vs_baseline_ee',
#             output_format='png', output_path='./fig/',
#             legend_loc='lower left', legend_ncol=1)
# print('Done')

# print('    Pem (time-dependent Hamiltonian calculation)... ', end='')
# oscprob2nu_td_plot.plot_probability_2nu_td_vs_baseline(
#             '12', energy=1.e-2,
#             log10_l_min=log10(5.e1), log10_l_max=3.0, log10_l_npts=3000,
#             plot_prob_ee=False, plot_prob_em=True, plot_prob_mm=False,
#             output_filename='prob_2nu_td_vacuum_matter_vs_baseline_em',
#             output_format='png', output_path='./fig/',
#             legend_loc='upper left', legend_ncol=1)
# print('Done')

# print('    Pmm (time-dependent Hamiltonian calculation)... ', end='')
# oscprob2nu_td_plot.plot_probability_2nu_td_vs_baseline(
#             '23', energy=1.e-1,
#             log10_l_min=1.0, log10_l_max=log10(3.e2), log10_l_npts=3000,
#             plot_prob_ee=True, plot_prob_em=False, plot_prob_mm=False,
#             output_filename='prob_2nu_td_vacuum_matter_vs_baseline_mm',
#             output_format='png', output_path='./fig/',
#             legend_loc='upper left', legend_ncol=1)
# print('Done')

# print('    Pmt (time-dependent Hamiltonian calculation)... ', end='')
# oscprob2nu_td_plot.plot_probability_2nu_td_vs_baseline(
#             '23', energy=1.e-1,
#             log10_l_min=1.0, log10_l_max=log10(3.e2), log10_l_npts=3000,
#             plot_prob_ee=False, plot_prob_em=True, plot_prob_mm=False,
#             output_filename='prob_2nu_td_vacuum_matter_vs_baseline_mt',
#             output_format='png', output_path='./fig/',
#             legend_loc='lower left', legend_ncol=1)
# print('Done')

# print('    Pee solar (time-dependent Hamiltonian calculation)... ', end='')
# oscprob2nu_td_plot.plot_probability_2nu_td_vs_baseline_sun(
#             '12', energy=5e-3,#1.e-2,
#             log10_l_min=-2, log10_l_max=log10(1e2), log10_l_npts=5000,
#             integration_method='simpson_log', 
#             epsrel=1.e-10, epsabs=1.e-10, num_pts_integration=1001,
#             plot_prob_ee=True, plot_prob_em=False, plot_prob_mm=False,
#             output_filename='prob_2nu_td_vacuum_matter_vs_baseline_ee_sun',
#             output_format='png', output_path='./fig/',
#             legend_loc='lower left', legend_ncol=1)
# print('Done')

# print('  Done')

# print()

# print('Generating oscillograms 2nu for Earth (time-dependent Hamiltonian):')

# print('    Pee (time-dependent Hamiltonian calculation)... ', end='')
# oscprob2nu_td_plot.plot_oscillogram_earth_2nu_td(
#             '12', prob_sel='ee',
#             costhz_min=-1.0, costhz_max=0.0, costhz_npts=200,
#             log10_Enu_min=0.0, log10_Enu_max=1.0, Enu_npts=200,
#             output_filename='oscillogram_earth_2nu_td_ee', 
#             output_format='png', output_path='./fig/')
# print('Done')

# print('    Pmm (time-dependent Hamiltonian calculation)... ', end='')
# oscprob2nu_td_plot.plot_oscillogram_earth_2nu_td(
#             '23', prob_sel='mm',
#             costhz_min=-1.0, costhz_max=0.0, costhz_npts=200,
#             log10_Enu_min=0.0, log10_Enu_max=1.0, Enu_npts=200,
#             output_filename='oscillogram_earth_2nu_td_mm', 
#             output_format='png', output_path='./fig/')
# print('Done')

# print('  Done')

# print()

# print('Generating plots of 2nu probability vs. energy:')

# for case in ['vacuum']:#, 'matter', 'nsi', 'liv']:

#     print('  Case: '+case)

#     print('    Pee, Pem... ', end='')
#     oscprob2nu_plot.plot_probability_2nu_vs_energy(
#                 case, '12', baseline=1.e3,
#                 log10_energy_min=-1.0, log10_energy_max=1.0,
#                 log10_energy_npts=600,
#                 plot_prob_ee=True, plot_prob_em=True, plot_prob_mm=False,
#                 output_filename='prob_2nu_'+case+'_vs_energy_ee_em',
#                 output_format='png', output_path='./fig/',
#                 legend_loc='center right', legend_ncol=1)
#     print('Done')

#     print('    Pmm, Pmt... ', end='')
#     oscprob2nu_plot.plot_probability_2nu_vs_energy(
#                 case, '23', baseline=1.e3,
#                 log10_energy_min=-1.0, log10_energy_max=1.0,
#                 log10_energy_npts=600,
#                 plot_prob_ee=True, plot_prob_em=True, plot_prob_mm=False,
#                 output_filename='prob_2nu_'+case+'_vs_energy_mm_mt',
#                 output_format='png', output_path='./fig/',
#                 legend_loc='center right', legend_ncol=1)
#     print('Done')

#     print('  Done')

# print()

# print('Generating plots of 2nu probability vs. density:')

# print('    Pee (time-dependent Hamiltonian calculation)... ', end='')
# oscprob2nu_td_plot.plot_probability_2nu_td_vs_density(
#             '12', 
#             l_min=0.0, l_max=SUN_RADIUS,
#             num_density_e_min=0, num_density_e_max=300, 
#             num_density_e_npts=1000,
#             integration_method='quad_linear', 
#             epsrel=1.e-10, epsabs=1.e-10, num_pts_integration=1001,
#             plot_prob_ee=True, plot_prob_em=False, plot_prob_mm=False,
#             output_filename='prob_2nu_td_vacuum_matter_vs_density_ee', 
#             output_format='png', output_path='./fig/', 
#             legend_loc='center left', legend_ncol=1)
# print('Done')

# print('  Done')

# print()

print('Generating plots of 2nu probability vs. energy (time-dependent Hamiltonian):')

print('    Pee solar (time-dependent Hamiltonian calculation)... ', end='')
oscprob2nu_td_plot.plot_probability_2nu_td_vs_energy_sun(
            '12', 
            # r_ini=1.e-3*SUN_RADIUS, r_fin=SUN_RADIUS,
            r_ini=1.e-3*SUN_RADIUS, r_fin=SUN_RADIUS,
            log10_energy_min=-1.0, log10_energy_max=2.0, 
            log10_energy_npts=10000,
            integration_method='quad_log', 
            epsrel=1.e-10, epsabs=1.e-10, num_pts_integration=1001,
            plot_prob_ee=True, plot_prob_em=False, plot_prob_mm=False,
            output_filename='prob_2nu_td_vacuum_matter_vs_energy_ee_sun', 
            output_format='png', output_path='./fig/', 
            legend_loc='center left', legend_ncol=1)
print('Done')

print('  Done')

print()


# print('Generating plots of 3nu probability vs. baseline:')

# for case in ['vacuum', 'matter', 'nsi', 'liv']:

#     print('  Case: '+case)

#     print('    Pee, Pem, Pet... ', end='')
#     # Pee, Pem, Pet
#     oscprob3nu_plot.plot_probability_3nu_vs_baseline(
#                 case, energy=1.e-2,
#                 log10_l_min=0.0, log10_l_max=3.0, log10_l_npts=1000,
#                 plot_prob_ee=True, plot_prob_em=True, plot_prob_et=True,
#                 plot_prob_me=False, plot_prob_mm=False, plot_prob_mt=False,
#                 plot_prob_te=False, plot_prob_tm=False, plot_prob_tt=False,
#                 output_filename='prob_3nu_'+case+'_vs_baseline_ee_em_et',
#                 output_format='png', output_path='./fig/',
#                 legend_loc='center left', legend_ncol=1)
#     print('Done')

#     print('    Pme, Pmm, Pmt... ', end='')
#     # Pme, Pmm, Pmt
#     oscprob3nu_plot.plot_probability_3nu_vs_baseline(
#                 case, energy=1.e-2,
#                 log10_l_min=0.0, log10_l_max=3.0, log10_l_npts=1000,
#                 plot_prob_ee=False, plot_prob_em=False, plot_prob_et=False,
#                 plot_prob_me=True, plot_prob_mm=True, plot_prob_mt=True,
#                 plot_prob_te=False, plot_prob_tm=False, plot_prob_tt=False,
#                 output_filename='prob_3nu_'+case+'_vs_baseline_me_mm_mt',
#                 output_format='png', output_path='./fig/',
#                 legend_loc='center left', legend_ncol=1)
#     print('Done')

#     print('    Pte, Ptm, Ptt... ', end='')
#     # Pte, Ptm, Ptt
#     oscprob3nu_plot.plot_probability_3nu_vs_baseline(
#                 case, energy=1.e-2,
#                 log10_l_min=0.0, log10_l_max=3.0, log10_l_npts=1000,
#                 plot_prob_ee=False, plot_prob_em=False, plot_prob_et=False,
#                 plot_prob_me=False, plot_prob_mm=False, plot_prob_mt=False,
#                 plot_prob_te=True, plot_prob_tm=True, plot_prob_tt=True,
#                 output_filename='prob_3nu_'+case+'_vs_baseline_te_tm_tt',
#                 output_format='png', output_path='./fig/',
#                 legend_loc='center left', legend_ncol=1)
#     print('Done')
#     print('  Done')

# print()

# print('Generating plots of 3nu probability vs. baseline (time-dependent Hamiltonian):')

# # --- Uncomment this when validating vacuum calculations (see oscprob3nu_td_plot.y)
# print('    Pem (time-dependent Hamiltonian calculation)... ', end='')
# oscprob3nu_td_plot.plot_probability_3nu_td_vs_baseline(
#             energy=1.e-2,
#             log10_l_min=log10(5.e1), log10_l_max=4.0, log10_l_npts=200,
#             plot_prob_ee=False, plot_prob_em=True, plot_prob_et=False,
#             plot_prob_me=False, plot_prob_mm=False, plot_prob_mt=False,
#             plot_prob_te=False, plot_prob_tm=False, plot_prob_tt=False,
#             output_filename='prob_3nu_td_vacuum_vs_baseline_em_validation',
#             output_format='png', output_path='./fig/',
#             legend_loc='lower left', legend_ncol=1)
# print('Done')


# print('    Pee (time-dependent Hamiltonian calculation)... ', end='')
# oscprob3nu_td_plot.plot_probability_3nu_td_vs_baseline(
#             energy=1.e-2,
#             log10_l_min=log10(5.e1), log10_l_max=3.0, log10_l_npts=3000,
#             plot_prob_ee=True, plot_prob_em=False, plot_prob_et=False,
#             plot_prob_me=False, plot_prob_mm=False, plot_prob_mt=False,
#             plot_prob_te=False, plot_prob_tm=False, plot_prob_tt=False,
#             output_filename='prob_3nu_td_vacuum_matter_vs_baseline_ee',
#             output_format='png', output_path='./fig/',
#             legend_loc='lower left', legend_ncol=1)
# print('Done')

# print('    Pem (time-dependent Hamiltonian calculation)... ', end='')
# oscprob3nu_td_plot.plot_probability_3nu_td_vs_baseline(
#             energy=1.e-2,
#             log10_l_min=log10(5.e1), log10_l_max=3.0, log10_l_npts=3000,
#             plot_prob_ee=False, plot_prob_em=True, plot_prob_et=False,
#             plot_prob_me=False, plot_prob_mm=False, plot_prob_mt=False,
#             plot_prob_te=False, plot_prob_tm=False, plot_prob_tt=False,
#             output_filename='prob_3nu_td_vacuum_matter_vs_baseline_em',
#             output_format='png', output_path='./fig/',
#             legend_loc='upper left', legend_ncol=1)
# print('Done')

# print('    Pet (time-dependent Hamiltonian calculation)... ', end='')
# oscprob3nu_td_plot.plot_probability_3nu_td_vs_baseline(
#             energy=1.e-2,
#             log10_l_min=log10(5.e1), log10_l_max=3.0, log10_l_npts=3000,
#             plot_prob_ee=False, plot_prob_em=False, plot_prob_et=True,
#             plot_prob_me=False, plot_prob_mm=False, plot_prob_mt=False,
#             plot_prob_te=False, plot_prob_tm=False, plot_prob_tt=False,
#             output_filename='prob_3nu_td_vacuum_matter_vs_baseline_et',
#             output_format='png', output_path='./fig/',
#             legend_loc='upper left', legend_ncol=1)
# print('Done')

# print('    Pmm (time-dependent Hamiltonian calculation)... ', end='')
# oscprob3nu_td_plot.plot_probability_3nu_td_vs_baseline(
#             energy=1.e-2,
#             log10_l_min=log10(5.e1), log10_l_max=3.0, log10_l_npts=3000,
#             plot_prob_ee=False, plot_prob_em=False, plot_prob_et=False,
#             plot_prob_me=False, plot_prob_mm=True, plot_prob_mt=False,
#             plot_prob_te=False, plot_prob_tm=False, plot_prob_tt=False,
#             output_filename='prob_3nu_td_vacuum_matter_vs_baseline_mm',
#             output_format='png', output_path='./fig/',
#             legend_loc='upper left', legend_ncol=1)
# print('Done')

# print('    Pmt (time-dependent Hamiltonian calculation)... ', end='')
# oscprob3nu_td_plot.plot_probability_3nu_td_vs_baseline(
#             energy=1.e-2,
#             log10_l_min=log10(5.e1), log10_l_max=3.0, log10_l_npts=3000,
#             plot_prob_ee=False, plot_prob_em=False, plot_prob_et=False,
#             plot_prob_me=False, plot_prob_mm=False, plot_prob_mt=True,
#             plot_prob_te=False, plot_prob_tm=False, plot_prob_tt=False,
#             output_filename='prob_3nu_td_vacuum_matter_vs_baseline_mt',
#             output_format='png', output_path='./fig/',
#             legend_loc='upper left', legend_ncol=1)
# print('Done')

# print('    Ptt (time-dependent Hamiltonian calculation)... ', end='')
# oscprob3nu_td_plot.plot_probability_3nu_td_vs_baseline(
#             energy=1.e-2,
#             log10_l_min=log10(5.e1), log10_l_max=3.0, log10_l_npts=3000,
#             plot_prob_ee=False, plot_prob_em=False, plot_prob_et=False,
#             plot_prob_me=False, plot_prob_mm=False, plot_prob_mt=False,
#             plot_prob_te=False, plot_prob_tm=False, plot_prob_tt=True,
#             output_filename='prob_3nu_td_vacuum_matter_vs_baseline_tt',
#             output_format='png', output_path='./fig/',
#             legend_loc='upper left', legend_ncol=1)
# print('Done')

# print()

# print('Generating oscillograms 3nu for Earth (time-dependent Hamiltonian):')

# print('    Pee (time-dependent Hamiltonian calculation)... ', end='')
# oscprob3nu_td_plot.plot_oscillogram_earth_3nu_td(
#             prob_sel='ee',
#             costhz_min=-1.0, costhz_max=0.0, costhz_npts=200,
#             log10_Enu_min=0.0, log10_Enu_max=1.0, Enu_npts=200,
#             output_filename='oscillogram_earth_3nu_td_ee', 
#             output_format='png', output_path='./fig/')
# print('Done')

# print('    Pem (time-dependent Hamiltonian calculation)... ', end='')
# oscprob3nu_td_plot.plot_oscillogram_earth_3nu_td(
#             prob_sel='em',
#             costhz_min=-1.0, costhz_max=0.0, costhz_npts=200,
#             log10_Enu_min=0.0, log10_Enu_max=1.0, Enu_npts=200,
#             output_filename='oscillogram_earth_3nu_td_em', 
#             output_format='png', output_path='./fig/')
# print('Done')

# print('    Pet (time-dependent Hamiltonian calculation)... ', end='')
# oscprob3nu_td_plot.plot_oscillogram_earth_3nu_td(
#             prob_sel='et',
#             costhz_min=-1.0, costhz_max=0.0, costhz_npts=200,
#             log10_Enu_min=0.0, log10_Enu_max=1.0, Enu_npts=200,
#             output_filename='oscillogram_earth_3nu_td_et', 
#             output_format='png', output_path='./fig/')
# print('Done')

# print('    Pmm (time-dependent Hamiltonian calculation)... ', end='')
# oscprob3nu_td_plot.plot_oscillogram_earth_3nu_td(
#             prob_sel='mm',
#             costhz_min=-1.0, costhz_max=0.0, costhz_npts=200,
#             log10_Enu_min=0.0, log10_Enu_max=1.0, Enu_npts=200,
#             output_filename='oscillogram_earth_3nu_td_mm', 
#             output_format='png', output_path='./fig/')
# print('Done')

# print('    Pmt (time-dependent Hamiltonian calculation)... ', end='')
# oscprob3nu_td_plot.plot_oscillogram_earth_3nu_td(
#             prob_sel='mt',
#             costhz_min=-1.0, costhz_max=0.0, costhz_npts=200,
#             log10_Enu_min=0.0, log10_Enu_max=1.0, Enu_npts=200,
#             output_filename='oscillogram_earth_3nu_td_mt', 
#             output_format='png', output_path='./fig/')
# print('Done')

# print('    Ptt (time-dependent Hamiltonian calculation)... ', end='')
# oscprob3nu_td_plot.plot_oscillogram_earth_3nu_td(
#             prob_sel='tt',
#             costhz_min=-1.0, costhz_max=0.0, costhz_npts=200,
#             log10_Enu_min=0.0, log10_Enu_max=1.0, Enu_npts=200,
#             output_filename='oscillogram_earth_3nu_td_tt', 
#             output_format='png', output_path='./fig/')
# print('Done')

# print('  Done')

# print()

# print('Generating plots of 3nu probability vs. energy:')

# for case in ['vacuum', 'matter', 'nsi', 'liv']:

#     print('  Case: '+case)

#     print('    Pee, Pem, Pet... ', end='')
#     # Pee, Pem, Pet
#     oscprob3nu_plot.plot_probability_3nu_vs_energy(
#                 case, baseline=1.e3,
#                 log10_energy_min=-1.0, log10_energy_max=1.0,
#                 log10_energy_npts=200,
#                 plot_prob_ee=True, plot_prob_em=True, plot_prob_et=True,
#                 plot_prob_me=False, plot_prob_mm=False, plot_prob_mt=False,
#                 plot_prob_te=False, plot_prob_tm=False, plot_prob_tt=False,
#                 output_filename='prob_3nu_'+case+'_vs_energy_ee_em_et',
#                 output_format='png', output_path='./fig/',
#                 legend_loc='center right', legend_ncol=1)
#     print('Done')

#     print('    Pme, Pmm, Pmt... ', end='')
#     # Pme, Pmm, Pmt
#     oscprob3nu_plot.plot_probability_3nu_vs_energy(
#                 case, baseline=1.e3,
#                 log10_energy_min=-1.0, log10_energy_max=1.0,
#                 log10_energy_npts=200,
#                 plot_prob_ee=False, plot_prob_em=False, plot_prob_et=False,
#                 plot_prob_me=True, plot_prob_mm=True, plot_prob_mt=True,
#                 plot_prob_te=False, plot_prob_tm=False, plot_prob_tt=False,
#                 output_filename='prob_3nu_'+case+'_vs_energy_me_mm_mt',
#                 output_format='png', output_path='./fig/',
#                 legend_loc='center right', legend_ncol=1)
#     print('Done')

#     print('    Pte, Ptm, Ptt... ', end='')
#     # Pte, Ptm, Ptt
#     oscprob3nu_plot.plot_probability_3nu_vs_energy(
#                 case, baseline=1.e3,
#                 log10_energy_min=-1.0, log10_energy_max=1.0,
#                 log10_energy_npts=200,
#                 plot_prob_ee=False, plot_prob_em=False, plot_prob_et=False,
#                 plot_prob_me=False, plot_prob_mm=False, plot_prob_mt=False,
#                 plot_prob_te=True, plot_prob_tm=True, plot_prob_tt=True,
#                 output_filename='prob_3nu_'+case+'_vs_energy_te_tm_tt',
#                 output_format='png', output_path='./fig/',
#                 legend_loc='center right', legend_ncol=1)
#     print('Done')
#     print('  Done')

# print()

# print('Generating 2nu plot in paper... ', end='')
# oscprob2nu_plotpaper.plot_probability_2nu_vs_energy_compare( \
#     output_format='png', output_path='./fig/')
# print('Done')

# print()

# print('Generating 3nu plot in paper... ', end='')
# oscprob3nu_plotpaper.plot_probability_3nu_vs_energy_compare( \
#     output_format='png', output_path='./fig/')
# print('Done')
