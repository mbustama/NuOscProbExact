# -*- coding: utf-8 -*-
r"""Compute three-neutrino Hamiltonians for selected scenarios.

This module contains the routines to compute the three-neutrino
Hamiltonians for the following scenarios: oscillations in vacuum, in
matter of constant density, in matter with non-standard interactions
(NSI), and in a CPT-odd Lorentz invariance-violating background (LIV).

Routine listings
----------------

    * hamiltonian_3nu_vacuum_energy_independent - Returns H_vac (no 1/E)
    * hamiltonian_3nu_matter - Returns H_matter

Created: 2024/11/03 21:29
Last modified: 2024/11/03 21:29
"""


__version__ = "1.0"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


from numpy import *
import numpy as np
import cmath
import cmath as cmath
import copy as cp

import hamiltonians3nu
from globaldefs import *



def hamiltonian_3nu_vacuum_energy_independent_td(s12, s23, s13, dCP, 
    D21, D31, compute_matrix_multiplication=False):
    r"""Returns the three-neutrino Hamiltonian for vacuum oscillations.

    Computes and returns the 3x3 complex three-neutrino Hamiltonian for
    oscillations in vacuum, parametrized by three mixing angles ---
    theta_12, theta_23, theta_13 --- one CP-violation phase --- delta_CP
    --- and two mass-squared difference --- Delta m^2_21, Delta m^2_31.
    The Hamiltonian is H = (1/2)*R.M2.R^dagger, with R the 3x3 PMNS
    matrix and M2 the mass matrix.  The multiplicative factor 1/E is not
    applied.

    Parameters
    ----------
    s12 : float
        Sin(theta_12).
    s23 : float
        Sin(theta_23).
    s13 : float
        Sin(theta_13).
    D21 : float
        Mass-squared difference Delta m^2_21.
    D31 : float
        Mass-squared difference Delta m^2_31.
    compute_matrix_multiplication : bool, optional
        If False (default), use the pre-computed expressions; otherwise,
        multiply R.M2.R^dagger live.

    Returns
    -------
    list
        Hamiltonian 3x3 matrix.
    """
    H = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent( \
        s12, s23, s13, dCP, D21, D31,
        compute_matrix_multiplication=compute_matrix_multiplication)

    return H


def hamiltonian_3nu_matter_td(h_vacuum_energy_independent, l, energy, 
    VCC_func):
    r"""Returns the three-neutrino Hamiltonian for matter oscillations.

    Computes and returns the 3x3 real three-neutrino Hamiltonian for
    oscillations in matter with constant density.

    Parameters
    ----------
    h_vacuum_energy_independent : list
        Energy-independent part of the three-neutrino Hamiltonian for
        oscillations in vacuum.  This is computed by the routine
        hamiltonian_3nu_vacuum_energy_independent.
    energy : float
        Neutrino energy.
    VCC : float
        Potential due to charged-current interactions of nu_e with
        electrons.

    Returns
    -------
    list
        Hamiltonian 3x3 matrix.
    """
    h_matter = hamiltonians3nu.hamiltonian_3nu_matter( \
        h_vacuum_energy_independent, energy, VCC_func(l))

    return h_matter


