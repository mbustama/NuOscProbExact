# -*- coding: utf-8 -*-
r"""Compute the two-neutrino flavor-transition probability, for a time-
dependent (or position-dependent) Hamiltonian, using time slices.

This module contains all the necessary routines to compute the
two-neutrino flavor-transition probabilities using the SU(2)
exponential expansion, for a time-dependent Hamiltonian.

Routine listings
----------------

    * probabilities_2nu_td - Returns the oscillation probabilities

References
----------

.. [1] Mauricio Bustamante, "NuOscProbExact: a general-purpose code to
   compute exact two-flavor and three-flavor neutrino oscillation 
   probabilities", arXiv:1904.12391.

Created: 2024/12/16 10:36
Last modified: 2024/12/16 10:36
"""


__version__ = "1.0"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


from numpy import *
import numpy as np
import cmath
import cmath as cmath
import scipy as sp
import oscprob2nu
# from mpmath import mp


def evolution_operator_2nu_tslice(hamiltonian_matrix_func, l_init, 
    l_final, l_npts, use_td_calculation=False, **kwargs):

    if (l_init == l_final):
        U = np.identity(2)
    else:
        # l = np.linspace(l_init, l_final, l_npts)
        l = np.logspace(log10(l_init), log10(l_final), l_npts)
        if use_td_calculation:
            pass
        else:
            U = np.identity(2)
            for i in range(1, l_npts):
                Dl = l[i]-l[i-1]
                U_new = np.array(oscprob2nu.evolution_operator_2nu( \
                    hamiltonian_matrix_func(l[i-1]+Dl/2), L=Dl))
                U = np.matmul(U_new, U)

    return U


def probabilities_2nu_tslice(hamiltonian_matrix_func, l_init, l_final,
    l_npts, use_td_calculation=False, **kwargs):

    if (l_init == l_final):
        Pee, Pem, Pme, Pmm  = 1.0, 0.0, 0.0, 1.0
    else:
        # Compute the evolution operator in time (or distance) slices
        U = evolution_operator_2nu_tslice(hamiltonian_matrix_func,
            l_init, l_final, l_npts, 
            use_td_calculation=use_td_calculation, **kwargs)
        # Compute probability nu_a --> nu_b at the end of the interval
        nu_a = np.array([[1], [0]])
        nu_b = np.array([0, 1])
        Pem = abs(np.matmul(nu_b, np.matmul(U, nu_a))[0])**2
        Pme = Pem
        Pee = 1.0-Pem
        Pmm = 1.0-Pme

    return Pee, Pem, Pme, Pmm
