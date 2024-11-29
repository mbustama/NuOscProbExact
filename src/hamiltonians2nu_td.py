# -*- coding: utf-8 -*-
r"""Compute two-neutrino time-dependent Hamiltonians for selected 
scenarios.

This module contains the routines to compute the two-neutrino
Hamiltonians for the following scenarios: oscillations in vacuum, in
matter of constant density, in matter with non-standard interactions
(NSI), and in a CPT-odd Lorentz invariance-violating background (LIV).

Routine listings
----------------

    * hamiltonian_2nu_vacuum_energy_independent - Returns H_vac (no 1/E)
    * hamiltonian_2nu_matter - Returns H_matter

Created: 2024/11/28 21:00
Last modified: 2024/11/29 17:46
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
    r"""Returns the two-neutrino Hamiltonian for vacuum oscillations,
    as a function of distance, even if it does not depend on it.

    Computes and returns the 2x2 real two-neutrino Hamiltonian for
    oscillations in vacuum, as a function of distance, parametrized by
    a single mixing angle theta and a single mass-squared difference 
    Dm2.  The Hamiltonian is H = (1/2)*R.M2.R^dagger, with R the 2x2
    rotation matrix and M2 the mass matrix.  The multiplicative factor
    1/E is not applied.  The vacuum Hamiltonian does not depend on
    distance in reality, but we include the dependence here as a way
    to validate the routine to compute probabilities for time-
    dependent Hamiltonians.

    Parameters
    ----------
    l : float
        Position at which the Hamiltonian is evaluated.
    sth : float
        Sin(theta).
    Dm2 : float
        Mass-squared difference Delta m^2.
    compute_matrix_multiplication : bool, optional
        If False (default), use the pre-computed expressions; otherwise,
        multiply R.M2.R^dagger live.

    Returns 
    -------
    list
        Hamiltonian 2x2 matrix.
    """

    H = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(sth, 
        Dm2, compute_matrix_multiplication)

    return H


def density_matter_func_const(l, 
    density_matter_const=DENSITY_MATTER_CRUST_G_PER_CM3):
    r"""Returns the matter density as a function of position, assuming a 
    constant density. Used for testing purposes.

    Returns the matter density as a function of position, assuming a 
    constant density. Used for testing purposes.

    Parameters
    ----------
    l : float
        Position at which the density profile is evaluated (in this
        case, the profile is uniform, so any value of l returns the same
        constant density).
    density_matter_const : float
        Matter density [g cm^{-3}]

    Returns
    -------
    float
        Matter density [g cm^{-3}]
    """

    return density_matter_const


def density_matter_func_exp(l, density_matter_central, l_scale):
    r"""Returns the matter density as a function of position, assuming  
    an exponentially decreasing density profile.

    Returns the matter density as a function of position, assuming  
    an exponentially decreasing density profile of the form
    rho(l) = density_matter_central*exp(-l/l_scale), for given values
    of density_matter_central and l_scale.

    Parameters
    ----------
    l : float
        Position at which the density profile is evaluated (in this
        case, the profile is uniform, so any value of l returns the same
        constant density).
    density_matter_central : float
        Matter density at the center of the profile (l = 0) [g cm^{-3}]
    l_scale : float
        Length scale of the exponential density decrease.

    Returns
    -------
    float
        Matter density [g cm^{-3}]
    """

    density = density_matter_central*exp(-l/l_scale)

    return density


def num_density_e_func(l, density_matter_func, electron_fraction=0.5):
    r"""Converts matter density [g cm^{-3}] to electron number density
    [eV^3], for a given matter density profile and position.

    Converts the matter density [g cm^{-3}] to electron number density
    [eV^3], for a given matter density profile, density_matter_func,
    and position, l. Matter is assumed to be isoscalar, with the
    fraction of electrons given by electron_fraction.

    Parameters
    ----------
    l : float
        Position at which the density profile is evaluated (in this
        case, the profile is uniform, so any value of l returns the same
        constant density).
    density_matter_funct : float(l)
        Matter density as a function of l [g cm^{-3}].
    electron_fraction : float
        Electron fraction.

    Returns
    -------
    float
        Number density of electrons [eV^3]
    """

    num_density_e = density_matter_func(l) * CONV_G_TO_EV \
                        / ((MASS_PROTON+MASS_NEUTRON)/2.0) \
                        * electron_fraction \
                        / pow(CONV_CM_TO_INV_EV, 3.0) # [eV^3]

    return num_density_e


def VCC_func(l, num_density_e_func):
    r"""Computes and returns the coherent forward electron potential, 
    V_CC, at position l, for a given electron number density, 
    num_density_e_func.

    Computes and returns the coherent forward electron potential, V_CC,
    at position l, for a given electron number density profile, 
    num_density_e_func.

    Parameters
    ----------
    l : float
        Position at which the density profile is evaluated (in this
        case, the profile is uniform, so any value of l returns the same
        constant density).
    num_density_e_func : float(l)
        Electron number density [eV^3].

    Returns
    -------
    float
        Coherent forward electron potntial, V_CC [eV]
    """

    VCC = sqrt(2.0)*GF*num_density_e_func(l) # [eV]

    return VCC


def hamiltonian_2nu_matter_td(h_vacuum_energy_independent, l, energy, 
    VCC_func):
    r"""Returns the two-neutrino Hamiltonian for matter oscillations,
    as a function of distance.

    Computes and returns the 2x2 real two-neutrino Hamiltonian for
    oscillations in matter with a given density as a function of
    position.

    Parameters
    ----------
    h_vacuum_energy_independent : list
        Energy-independent part of the two-neutrino Hamiltonian for
        oscillations in vacuum.  This is computed by the routine
        hamiltonian_2nu_vacuum_energy_independent.
    energy : float
        Neutrino energy.
    VCC : float
        Potential due to charged-current interactions of nu_e with
        electrons.

    Returns
    -------
    list
        Hamiltonian 2x2 matrix.
    """

    h_matter = hamiltonians2nu.hamiltonian_2nu_matter( \
        h_vacuum_energy_independent, energy, VCC_func(l))

    return h_matter
