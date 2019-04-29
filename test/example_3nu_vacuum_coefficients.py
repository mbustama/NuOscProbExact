# -*- coding: utf-8 -*-
r"""Run the vacuum coefficients 3nu example shown in README.md.

Runs the three-neutrino example of coefficients for oscillations in
vacuum shown in README.md

References
----------

.. [1] Mauricio Bustamante, "Exact neutrino oscillation probabilities:
   a fast general-purpose computation method for two and three neutrino
   flavors", arXiv:1904.XXXXX.

Created: 2019/04/29 23:48
Last modified: 2019/04/29 23:48
"""


from __future__ import print_function

__version__ = "1.0"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import sys
sys.path.append('../src')

import numpy as np

import oscprob3nu
import hamiltonians3nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_energy_indep = \
    hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(  S12_NO_BF,
                                                                S23_NO_BF,
                                                                S13_NO_BF,
                                                                DCP_NO_BF,
                                                                D21_NO_BF,
                                                                D31_NO_BF)
h_vacuum = np.multiply(1./energy, h_vacuum_energy_indep)

h1, h2, h3, h4, h5, h6, h7, h8 = \
    oscprob3nu.hamiltonian_3nu_coefficients(h_vacuum)
print('h1: {:.4e}'.format(h1))
print('h2: {:.4e}'.format(h2))
print('h3: {:.4e}'.format(h3))
print('h4: {:.4e}'.format(h4))
print('h5: {:.4e}'.format(h5))
print('h6: {:.4e}'.format(h6))
print('h7: {:.4e}'.format(h7))
print('h8: {:.4e}'.format(h8))
print()

u0, u1, u2, u3, u4, u5, u6, u7, u8 = \
    oscprob3nu.evolution_operator_3nu_u_coefficients( \
                                                    h_vacuum,
                                                    baseline*CONV_KM_TO_INV_EV)
print('u0: {:.4f}'.format(u0))
print('u1: {:.4f}'.format(u1))
print('u2: {:.4f}'.format(u2))
print('u3: {:.4f}'.format(u3))
print('u4: {:.4f}'.format(u4))
print('u5: {:.4f}'.format(u5))
print('u6: {:.4f}'.format(u6))
print('u7: {:.4f}'.format(u7))
print('u8: {:.4f}'.format(u8))
print()

evol_operator = \
    oscprob3nu.evolution_operator_3nu(h_vacuum, baseline*CONV_KM_TO_INV_EV)
print('U3 = ')
with np.printoptions(precision=3, suppress=True):
    print(np.array(evol_operator))
