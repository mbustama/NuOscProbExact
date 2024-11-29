# -*- coding: utf-8 -*-
r"""Compute two-neutrino Hamiltonians for selected scenarios.

This module contains the routines to compute the two-neutrino
Hamiltonians for the following scenarios: oscillations in vacuum, in
matter of constant density, in matter with non-standard interactions
(NSI), and in a CPT-odd Lorentz invariance-violating background (LIV).

Routine listings
----------------

    * hamiltonian_2nu_vacuum_energy_independent - Returns H_vac (no 1/E)
    * hamiltonian_2nu_matter - Returns H_matter
    * hamiltonian_2nu_nsi - Returns H_NSI
    * hamiltonian_2nu_liv - Returns H_LIV

Created: 2024/11/28 21:00
Last modified: 2024/11/29 01:24
"""


__version__ = "1.0"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


from numpy import *
import numpy as np
import cmath
import cmath as cmath
import copy as cp

import hamiltonians2nu
from globaldefs import *



def hamiltonian_2nu_vacuum_energy_independent_td(l, sth, Dm2,
    compute_matrix_multiplication=False):
#     r"""Returns the two-neutrino Hamiltonian for vacuum oscillations,
#     r"""as a function of distance, even if it does not depend on it.

#     Computes and returns the 2x2 real two-neutrino Hamiltonian for
#     oscillations in vacuum, as a function of distance, parametrized by
#     a single mixing angle theta and a single mass-squared difference 
#     Dm2.  The Hamiltonian is H = (1/2)*R.M2.R^dagger, with R the 2x2
#     rotation matrix and M2 the mass matrix.  The multiplicative factor
#     1/E is not applied.  The vacuum Hamiltonian does not depend on
#     distance in reality, but we include the dependence here as a way
#     to validate the routine to compute probabilities for time-
#     dependent Hamiltonians.

#     Parameters
#     ----------
#     l : float
#         Position at which the Hamiltonian is evaluated.
#     sth : float
#         Sin(theta).
#     Dm2 : float
#         Mass-squared difference Delta m^2.
#     compute_matrix_multiplication : bool, optional
#         If False (default), use the pre-computed expressions; otherwise,
#         multiply R.M2.R^dagger live.

#     Returns
#     -------
#     list
#         Hamiltonian 2x2 matrix.
#     """

    H = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(sth, 
        Dm2, compute_matrix_multiplication)

    return H