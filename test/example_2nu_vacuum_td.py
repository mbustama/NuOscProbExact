# -*- coding: utf-8 -*-
r"""Run the vacuum 2nu example shown in README.md.

Runs the two-neutrino example of oscillations in vacuum shown in
README.md

References
----------

.. [1] Mauricio Bustamante, "NuOscProbExact: a general-purpose code to
   compute exact two-flavor and three-flavor neutrino oscillation 
   probabilities", arXiv:1904.12391.

Created: 2024/11/28 21:00
Last modified: 2024/11/29 01:15
"""


from __future__ import print_function

__version__ = "1.0"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import sys
sys.path.append('../src')

import numpy as np

import oscprob2nu_td
import hamiltonians2nu_td
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_func = lambda l: np.multiply(1./energy, 
    hamiltonians2nu_td.hamiltonian_2nu_vacuum_energy_independent_td(l,
    S23_NO_BF, D31_NO_BF))

l_init, l_final = 0.0, baseline*CONV_KM_TO_INV_EV # [eV^{-1}]
Pee, Pem, Pme, Pmm = oscprob2nu_td.probabilities_2nu_td(
    h_vacuum_func, l_init, l_final, integration_method='quad',
    epsrel=1.e-8, epsabs=1.e-8)

print("Pee = %6.5f, Pem = %6.5f" % (Pee, Pem))
print("Pme = %6.5f, Pmm = %6.5f" % (Pme, Pmm))
