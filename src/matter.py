# -*- coding: utf-8 -*-
r"""Contains helper functions to compute the oscillation probability in
matter.

This module contains routines to common matter density profiles (e.g.,
constant, exponentially decreasing), electron number density, and
coherent forward scattering potential.

Routine listings
----------------

    * density_matter_func_const - Returns the density for a constant 
           matter density profile
    * density_matter_func_exp - Returns the density for an exponentially
           decreasing matter density profile
    * density_matter_prem - Returns the density inside the Earth using
           the Preliminary Reference Earth Model
    * num_density_e_func - Converts a matter density to an electron
           number density
    * VCC_func - Returns the potential for coherent forward electron
           scattering

Created: 2024/11/30 15:42
Last modified: 2024/11/30 21:23
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


def density_matter_func_prem(r):
    r"""Returns the matter density inside the Earth according to the 
    Preliminary Reference Earth Model (PREM).
    
    Returns the matter density inside the Earth according to the PREM,
    for a given radial distance measured from the center of the Earth.

    Parameters
    ----------
    r : float
        Radial distance measured from the center of the Earth [km].

    Returns
    -------
    float
        Matter density [g cm^{-3}].

    References
    ----------

    .. [1] Adam M. Dziewonski & Don L. Anderson, "Preliminary Reference
        Earth Model", Physics of the Earth and Planetary Interiors, 25,
        297 (1981).
    """

    RADIUS_EARTH = 6371.0 # [km]

    if (r > RADIUS_EARTH):
        print('Error: density_matter_prem: value of ' + \
                'l cannot be > RADIUS_EARTH = 6371 km')
        quit()

    x = r/RADIUS_EARTH

    if (0 <= r <= 1221.5):
        density = 13.0885-8.8381*x*x
    elif (1221.5 < r <= 3480.0):
        density = 12.5815-1.2638*x-3.6426*x*x-5.5281*x*x*x
    elif (3480.0 < r <= 5701.0):
        density = 7.9565-6.4761*x+5.5283*x*x-3.0807*x*x*x
    elif (5701.0 < r <= 5771.0):
        density = 5.3197-1.4836*x
    elif (5771.0 < r <= 5971.0):
        density = 11.2494-8.0298*x
    elif (5971.0 < r <= 6151.0):
        density = 7.1089-3.8045*x
    elif (6151.0 < r <= 6346.6):
        density = 2.6910+0.6924*x
    elif (6346.6 < r <= 6356.0):
        density = 2.900
    elif (6356.0 < r <= 6368.0):
        density = 2.600
    elif (6368.0 < r <= RADIUS_EARTH):
        density = 1.020

    return density


def distance_traveled_inside_earth(costhz):
    r"""Returns the distance traveled by a neutrino inside the Earth,
    traveling with a cosine of zenith angle costhz.
    
    Returns the length of the path traveled by a neutrino from the 
    surface ot the Earth, through it, until it reaches a detector. The
    direction of the neutrino is parametrized by the zenith angle of the
    neutrino. Assumes that the neutrino detector is on the surface of 
    the Earth, not underground. As a result, the distance is zero for 
    all values of costhz > 0.

    Parameters
    ----------
    costhz : float
        Cosine of the zenith angle of the neutrino.

    Returns
    -------
    float
        Path length inside the Earth [km].
    """
    if (costhz > 0.0):
        d = 0.0
    else:
        d = -2.0 * EARTH_RADIUS * costhz

    return d


def earth_radial_distance_from_depth(costhz, l):
    r"""Returns the radial distance measured from the center of the
    Earth to a position inside the Earth, given by costhz and l.
    
    A neutrino with direction given by the cosine of the zenith angle, 
    costhz, travels from l=0 to l=distance_traveled_inside_earth,
    computed below. The routine returns the radial distance to the 
    neutrino when its distance from its point of entry into the Earth is
    l.  

    Parameters
    ----------
    costhz : float
        Cosine of the zenith angle of the neutrino.
    l : float
        Distance of the neutrino from its point of entry into the Earth.

    Returns
    -------
    float
        Radial distance from to the neutrino [km].
    """
    d = distance_traveled_inside_earth(costhz)

    if (l > d):
        print('Error: earth_radial_distance_from_depth: value of ' + \
                'l cannot be larger than the distance traveled ' + \
                'inside Earth for this value of costhz')
        quit()
    elif ((l == 0.0) and (costhz == 0.0)):
        r = 0.0
    else:
        r2 = EARTH_RADIUS*EARTH_RADIUS
        r2 += (d-l)**2
        r2 += 2*EARTH_RADIUS*(d-l)*costhz
        r = sqrt(r2)

    return abs(r)


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
    VCC = -sqrt(2.0)*GF*num_density_e_func(l) # [eV]

    return VCC