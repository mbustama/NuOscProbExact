# -*- coding: utf-8 -*-
r"""Compute the three-neutrino flavor-transition probability.

This module contains all the necessary routines to compute the
three-neutrino flavor-transition probabilities using the SU(3)
exponential expansion.

Routine listings
----------------

    * hamiltonian_3nu_coefficients - Returns coefficients of Hamiltonian
    * tensor_d - Returns the value of the tensor :math:`d_{i,jk}`
    * star - Returns the SU(3) product :math:`(h * h)_k`
    * su3_invariants - Returns the SU(3) invariants :math:`|h|^2, <L>`
    * psi_roots - Returns the roots of the characteristic equation
    * evolution_operator_3nu_u_coefficients - Returns the :math:`u_k`
    * evolution_operator_3nu - Returns evolution operator :math:`U_3`
    * probabilities_3nu - Returns the oscillation probabilities

References
----------

.. [1] A.J. MacFarlane, A. Sudbery, and P.H. Weisz, "On Gell-Mann's
   :math:`\lambda`-matrices, :math:`d`- and math`f`-tensors, octets, and
   parametrizations of SU(3)", Commun. Math. Phys. 11, 77 (1968).

.. [2] Mauricio Bustamante, "NuOscProbExact: a general-purpose code to
   compute exact two-flavor and three-flavor neutrino oscillation 
   probabilities", arXiv:1904.12391.

Created: 2024/12/03 21:16
Last modified: 2024/12/03 21:16
"""


__version__ = "1.0"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


from numpy import *
import numpy as np
import cmath
import cmath as cmath
import scipy as sp
import oscprob3nu


def probabilities_3nu_td(hamiltonian_matrix_func, l_init, l_final,
    integration_method='quad', epsrel=1.e-8, epsabs=1.e-8):
    r"""Returns the 3nu oscillation probability.

    Returns the three-neutrino oscillation probabilities
    Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt.

    Parameters
    ----------
    hamiltonian_matrix : list
        3x3 Hamiltonian, [[H11,H12,H13],[H21,H22,H23],[H31,H32,H33]].
    L : float
        Baseline.

    Returns
    -------
    list
        Three-neutrino probabilities
        Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt.

    Example
    -------
    >>> hamiltonian_matrix = [
    ...                   [1.0+0.0j, 0.0+2.0j, 0.0-1.0j],
    ...                   [0.0-2.0j, 3.0+0.0j, 3.0+0.0j],
    ...                   [0.0+1.0j, 3.0-0.0j, -5.0+0.0j]
    ... ]
    >>> L = 1.0
    >>> Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = \
    ...     probabilities_3nu(hamiltonian_matrix, 1.0)
    >>> print(Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt)
    0.342732 0.413691 0.243576 0.413691 0.0048504 0.58145 0.243576
    0.58145 0.174965
    """

    if (l_init == l_final):
        Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = \
            1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0
    else:
        # First integrate each of the components of the Hamiltonian
        h_integral = np.zeros((3,3))
        for i in range(9):
            if (integration_method == 'quad'):
                res, err = sp.integrate.quad(lambda l: \
                    hamiltonian_matrix_func(l)[i//3][i%3], 
                    l_init, l_final, epsrel=epsrel, epsabs=epsabs)
            else:
                print('Error: probabilities_3nu_td: ' + \
                    'value of integration_method not allowed')
                quit()
            h_integral[i//3][i%3] = res
        Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = \
             oscprob3nu.probabilities_3nu(h_integral, 1.0)

    return Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt
