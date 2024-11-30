# -*- coding: utf-8 -*-
r"""Contains helper functions to compute the oscillation probability in
matter.

This module contains routines to common matter density profiles (e.g.,
constant, exponentially decreasing), electron number density, and
coherent forward scattering potential.

Routine listings
----------------

    * density_matter_func_const - Returns XXXX
    * density_matter_func_exp - Returns H_matter
    * num_density_e_func - XXX
    * VCC_func - XXX

Created: 2024/11/30 15:42
Last modified: 2024/11/30 15:42
"""

__version__ = "1.0"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


from globaldefs import *


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

