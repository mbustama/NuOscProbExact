# -*- coding: utf-8 -*-
r"""Compute the two-neutrino flavor-transition probability, for a time-
dependent (or position-dependent) Hamiltonian.

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

Created: 2024/11/28 21:00
Last modified: 2024/11/30 21:28
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


def probabilities_2nu_td(hamiltonian_matrix_func, l_init, l_final,
    integration_method='quad', epsrel=1.e-8, epsabs=1.e-8):
    r"""Returns the 2nu oscillation probability.

    Returns the two-neutrino oscillation probabilities
    Pee, Pem, Pme, Pmm.

    Parameters
    ----------
    hamiltonian_matrix_func : list
        Function that returns a 2x2 Hamiltonian, [[H11,H12],[H21,H22]],
        as a function of position, l.
    l_init : float
        Initial position of the neutrino.
    l_final : float
        Final position of the neutrino.
    integration_method : string
        Flag to choose the method to integrate the Hamiltonian in l 
        (currently, only the scipy.integrate.quad is implemented).
    epsrel : float
        Relative error for the quad integration routine.
    epsabs : float
        Absolute error for the quad integration routine.

    Returns
    -------
    list
        Two-neutrino probabilities Pee, Pem, Pme, Pmm

    Example [=== OUTDATED ===]
    -------
    >>> hamiltonian_matrix = [
    ...                   [1.0+0.0j, 0.0+2.0j],
    ...                   [0.0-2.0j, 3.0+0.0j]
    ... ]
    >>> L = 1.0
    >>> Pee, Pem, Pme, Pmm = \
    ...     probabilities_2nu(hamiltonian_matrix, 1.0)
    >>> print(Pee, Pem, Pme, Pmm)
    0.504820 0.495179 0.495179 0.504820
    """

    if (l_init == l_final):
        Pee, Pem, Pme, Pmm  = 1.0, 0.0, 0.0, 1.0
    else:
        # First integrate each of the components of the Hamiltonian
        h_integral = np.zeros((2,2))
        for i in range(4):
            if (integration_method == 'quad'):
                res, err = sp.integrate.quad(lambda l: \
                    hamiltonian_matrix_func(l)[i//2][i%2], 
                    l_init, l_final, epsrel=epsrel, epsabs=epsabs)
            else:
                print('Error: probabilities_2nu_td: ' + \
                    'value of integration_method not allowed')
                quit()
            h_integral[i//2][i%2] = res

        # [h1, h2, h3]
        h_coeffs = oscprob2nu.hamiltonian_2nu_coefficients(h_integral)

        # h_abs = |h|
        h_abs = oscprob2nu.modulus(h_coeffs)

        Pem = (abs(h_coeffs[0])**2.0 + abs(h_coeffs[1])**2.0) / h_abs**2.0 \
            * pow(sin(h_abs), 2.0)
        Pme = Pem
        Pee = 1.0-Pem
        Pmm = 1.0-Pme

    return Pee, Pem, Pme, Pmm