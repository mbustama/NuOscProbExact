# -*- coding: utf-8 -*-
r"""Compute two-neutrino Hamiltonians for selected scenarios.

This module contains the routines to compute the two-neutrino
Hamiltonians for the following scenarios: oscillations in vacuum, in
matter of constant density, in matter with non-standard interactions
(NSI), and in a CPT-odd Lorentz invariance-violating background (LIV).

Routine listings
----------------

    * mixing_matrix_2nu - Returns 2x2 rotation matrix
    * hamiltonian_2nu_vacuum_energy_independent - Returns H_vac (no 1/E)
    * hamiltonian_2nu_matter - Returns H_matter
    * hamiltonian_2nu_nsi - Returns H_NSI
    * hamiltonian_2nu_liv - Returns H_LIV

Created: 2019/04/21 15:00
Last modified: 2019/04/23 21:04
"""


__version__ = "1.0"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


from numpy import *
import numpy as np
import cmath
import cmath as cmath
import copy as cp

import oscprob3nu
from globaldefs import *


def mixing_matrix_2nu(sth):
    r"""Returns the 2x2 rotation matrix.

    Computes and returns a 2x2 real rotation matrix parametrized by a
    single rotation angle theta.

    Parameters
    ----------
    sth : float
        Sin(theta).

    Returns
    -------
    list
        Rotation matrix [[cth, sth], [-sth, cth]], with cth = cos(theta)
        and sth = sin(theta).
    """
    cth = sqrt(1.0-sth*sth)

    U00 = cth
    U01 = sth
    U10 = -sth
    U11 = cth

    return [[U00,U01],[U10,U11]]


def hamiltonian_2nu_vacuum_energy_independent(sth, Dm2,
    compute_matrix_multiplication=False):
    r"""Returns the two-neutrino Hamiltonian for vacuum oscillations.

    Computes and returns the 2x2 real two-neutrino Hamiltonian for
    oscillations in vacuum, parametrized by a single mixing angle theta
    and a single mass-squared difference Dm2.  The Hamiltonian is
    H = (1/2)*R.M2.R^dagger, with R the 2x2 rotation matrix and M2 the
    mass matrix.  The multiplicative factor 1/E is not applied.

    Parameters
    ----------
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
    th = np.arcsin(sth)
    c2th = cos(2.0*th)
    s2th = sin(2.0*th)

    f = 1./4.

    if not compute_matrix_multiplication:

        H00 = Dm2*c2th
        H01 = -Dm2*s2th
        H10 = H01
        H11 = -H00

        H = [[H00*f,H01*f], [H10*f,H11*f]]

    else:

        # PMNS matrix
        R = np.array(mixing_matrix_2nu(sth))
        # Mass matrix
        M2 = np.array([[Dm2, 0.0], [0.0, -Dm2]])
        # Hamiltonian
        H = list(f*np.matmul(R, np.matmul(M2, matrix.transpose(R))))

    return H


def probabilities_2nu_vacuum_std(sth, Dm2, energy, L):
    r"""Returns 2nu oscillation vacuum probabilities, std. computation.

    Returns the probabilities for two-neutrino oscillations in vacuum,
    computed using the standard analytical expression of the
    probabilities.

    Parameters
    ----------
    sth : float
        Sin(theta).
    Dm2 : float
        Mass-squared difference Delta m^2.
    energy : float
        Neutrino energy.
    L : float
        Baseline.

    Returns
    -------
    list
        List of probabilities [Pee, Pem, Pme, Pmm].
    """
    arg = 1.27*Dm2*L/energy#/4.0
    cth = sqrt(1.0-sth*sth)
    s2th = 2.0*sth*cth

    Pem = s2th*s2th * pow(sin(arg), 2.0)
    Pme = Pem
    Pee = 1.0-Pem
    Pmm = 1.0-Pme

    prob = [Pee, Pem, Pme, Pmm]

    return prob


def hamiltonian_2nu_matter(h_vacuum_energy_independent, energy, VCC):
    r"""Returns the two-neutrino Hamiltonian for matter oscillations.

    Computes and returns the 2x2 real two-neutrino Hamiltonian for
    oscillations in matter with constant density.

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
    h_matter = cp.deepcopy(h_vacuum_energy_independent)
    h_matter = np.multiply(1.0/energy, h_matter)

    # Add the matter potential to the ee term to find the matter
    # Hamiltonian
    h_matter[0][0] += VCC

    return h_matter


def probabilities_2nu_matter_std(sth, Dm2, VCC, energy, L):
    r"""Returns 2nu oscillation matter probabilities, std. computation.

    Returns the probabilities for two-neutrino oscillations in matter,
    computed using the standard analytical expression of the
    probabilities.

    Parameters
    ----------
    sth : float
        Sin(theta).
    Dm2 : float
        Mass-squared difference Delta m^2.
    VCC : float
        Potential due to charged-current interactions of nu_e with
        electrons.
    energy : float
        Neutrino energy.
    L : float
        Baseline.

    Returns
    -------
    list
        List of probabilities [Pee, Pem, Pme, Pmm].
    """
    x = 2.0*VCC*(energy*1.e9)/Dm2
    cth = sqrt(1.0-sth*sth)
    s2th = 2.0*sth*cth
    s2thsq = s2th*s2th
    c2th = sqrt(1.0-s2thsq)

    Dm2m = Dm2*sqrt(s2thsq+pow(c2th-x, 2.0))
    s2thmsq = s2thsq / (s2thsq+pow(c2th-x, 2.0))

    arg = 1.27*Dm2m*L/energy#/4.0

    Pem = s2thmsq * pow(sin(arg), 2.0)
    Pme = Pem
    Pee = 1.0-Pem
    Pmm = 1.0-Pme

    prob = [Pee, Pem, Pme, Pmm]

    return prob


def hamiltonian_2nu_nsi(h_vacuum_energy_independent, energy, VCC, eps):
    r"""Returns the two-neutrino Hamiltonian for oscillations with NSI.

    Computes and returns the 2x2 real two-neutrino Hamiltonian for
    oscillations with non-standard interactions (NSI) in matter with
    constant density.

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
    eps : list
        Vector of NSI strength parameters: eps = eps_ee, eps_em, eps_mm.

    Returns
    -------
    list
        Hamiltonian 2x2 matrix.
    """
    h_nsi = cp.deepcopy(h_vacuum_energy_independent)
    h_nsi = np.multiply(1.0/energy, h_nsi)

    eps_ee, eps_em, eps_mm = eps

    h_nsi[0][0] += VCC*(1.0+eps_ee)
    h_nsi[0][1] += VCC*eps_em
    h_nsi[1][0] += VCC*np.conj(eps_em)
    h_nsi[1][1] += VCC*eps_mm

    return h_nsi


def hamiltonian_2nu_liv(h_vacuum_energy_independent, energy, sxi,
    b1, b2, Lambda):
    r"""Returns the two-neutrino Hamiltonian for oscillations with LIV.

    Computes and returns the 2x2 real two-neutrino Hamiltonian for
    oscillations in a CPT-odd Lorentz invariance-violating background.

    Parameters
    ----------
    h_vacuum_energy_independent : list
        Energy-independent part of the two-neutrino Hamiltonian for
        oscillations in vacuum.  This is computed by the routine
        hamiltonian_2nu_vacuum_energy_independent.
    energy : float
        Neutrino energy.
    sxi : float
        Sin(xi), with xi the rotation angle between the space of the
        eigenvectors of B2 and the flavor states.
    b1 : float
        Eigenvalue b1 of the LIV operator B2.
    b2 : float
        Eigenvalue b2 of the LIV operator B2.
    Lambda : float
        Energy scale of the LIV operator B2.

    Returns
    -------
    list
        Hamiltonian 2x2 matrix.
    """
    h_liv = cp.deepcopy(h_vacuum_energy_independent)
    h_liv = np.multiply(1.0/energy, h_liv)

    f = energy/Lambda
    cxi = sqrt(1.0-sxi-sxi)

    h_liv[0][0] += f*(b1*cxi*cxi + b2*sxi*sxi)
    h_liv[0][1] += f*((-b1+b2)*cxi*sxi)
    h_liv[1][0] += f*((-b1+b2)*cxi*sxi)
    h_liv[1][1] += f*(b2*cxi*cxi + b1*sxi*sxi)

    return h_liv
