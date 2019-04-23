# -*- coding: utf-8 -*-
r"""Compute three-neutrino Hamiltonians for selected scenarios.

This module contains the routines to compute the three-neutrino
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

Created: 2019/04/17 17:14
Last modified: 2019/04/23 22:04
"""

__version__ = "0.1"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


from numpy import *
import numpy as np
import cmath
import cmath as cmath
import copy as cp

import oscprob3nu
from globaldefs import *


def pmns_mixing_matrix(s12, s23, s13, dCP):
    r"""Returns the 3x3 PMNS mixing matrix.

    Computes and returns the 3x3 complex PMNS mixing matrix
    parametrized by three rotation angles, theta_12, theta_23, theta_13,
    and one CP-violation phase, delta_CP.

    Parameters
    ----------
    s12 : float
        Sin(theta_12).
    s23 : float
        Sin(theta_23).
    s13 : float
        Sin(theta_13).
    dCP : float
        delta_CP [radian].

    Returns
    -------
    list
        3x3 PMNS mixing matrix.
    """
    c12 = sqrt(1.0-s12*s12)
    c23 = sqrt(1.0-s23*s23)
    c13 = sqrt(1.0-s13*s13)

    cdCP = cos(dCP)
    sdCP = sin(dCP)

    U00 = c12*c13
    U01 = s12*c13
    U02 = s13*complex(cdCP,-sdCP)
    U10 = -s12*c23 - c12*s23*s13*complex(cdCP,sdCP)
    U11 = c12*c23 - s12*s23*s13*complex(cdCP,sdCP)
    U12 = s23*c13
    U20 = s12*s23 - c12*c23*s13*complex(cdCP,sdCP)
    U21 = -c12*s23 - s12*c23*s13*complex(cdCP,sdCP)
    U22 = c23*c13

    return [[U00,U01,U02],[U10,U11,U12],[U20,U21,U22]]


def hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31,
    compute_matrix_multiplication=False):
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
    c12 = sqrt(1.0-s12*s12)
    c23 = sqrt(1.0-s23*s23)
    c13 = sqrt(1.0-s13*s13)

    f = 1./2.

    if not compute_matrix_multiplication:

        # All Hij have units of [eV^2]
        H00 = c13*c13*D21*s12*s12 + D31*s13*s13
        H01 = c12*c13*c23*D21*s12 + \
                c13*(D31-D21*s12*s12)*s13*s23*complex(cos(dCP),-sin(dCP))
        H02 = c13*c23*(D31-D21*s12*s12)*s13*complex(cos(dCP),-sin(dCP)) - \
                c12*c13*D21*s12*s23
        H10 = c12*c13*c23*D21*s12 + \
                c13*(D31-D21*s12*s12)*s13*s23*complex(cos(dCP),sin(dCP))
        H11 = c12*c12*c23*c23*D21 + (c13*c13*D31 + D21*s12*s12*s13*s13)*s23*s23 - \
                2.0*c12*c23*D21*s12*s13*s23*cos(dCP)
        H12 = c13*c13*c23*D31*s23 + \
                (c23*s12*s13*complex(cos(dCP),-sin(dCP)) + c12*s23) * \
                (-c12*c23*D21 + D21*s12*s13*s23*complex(cos(dCP),sin(dCP)))
        H20 = c13*c23*(D31-D21*s12*s12)*s13*complex(cos(dCP),sin(dCP)) - \
                c12*c13*D21*s12*s23
        H21 = c13*c13*c23*D31*s23 - \
                D21*(c23*s12*s13*complex(cos(dCP),sin(dCP)) + c12*s23) * \
                (c12*c23 - s12*s13*s23*complex(cos(dCP),-sin(dCP)))
        H22 = c23*c23*(c13*c13*D31 + D21*s12*s12*s13*s13) + c12*c12*D21*s23*s23 + \
                2.0*c12*c23*D21*s12*s13*s23*cos(dCP)

        H = [[H00*f,H01*f,H02*f], [H10*f,H11*f,H12*f], [H20*f,H21*f,H22*f]]

    else:

        # PMNS matrix
        R = np.array(pmns_mixing_matrix(s12, s23, s13, dCP))
        # Mass matrix
        M2 = np.array([[0.0, 0.0, 0.0], [0.0, D21, 0.0], [0.0, 0.0, D31]])
        # Hamiltonian
        H = list(f*np.matmul(R, np.matmul(M2, np.conj(matrix.transpose(R)))))

    return H


def delta(a, b):
    r"""Returns the Kronecker delta function.

    Returns the delta function delta(a, b) = 1 if a == b and 0 if
    a != b.

    Parameters
    ----------
    a : int
        First index.
    b : int
        Second index.

    Returns
    -------
    int
        delta(a, b).
    """
    if (a == b):
        return 1
    else:
        return 0


def J(U, alpha, beta, k, j):
    r"""Returns U*_ak * U_bk * U_aj * U*_bj, with U the PMNS matrix.

    Returns the product U*_ak * U_bk * U_aj * U*_bj, where U is the
    PMNS mixing matrix.  This product appears in the standard expression
    for the three-neutrino oscillation probability in vacuum.

    Parameters
    ----------
    U : list
        3x3 PMNS complex mixing matrix.
    alpha : int
        Index of the initial flavor (0: e, 1: mu, 2: tau).
    beta : int
        Index of the final flavor (0: e, 1: mu, 2: tau).
    k : int
        First index of the sum over mass eigenstates (k = 0, 1, 2).
    j : int
        First index of the sum over mass eigenstates (k = 0, 1, 2).

    Returns
    -------
    float
        J(U, alpha, beta, j, j)
    """
    return np.conj(U[alpha][k])*U[beta][k]*U[alpha][j]*np.conj(U[beta][j])


def probabilities_3nu_std(U, D21, D31, energy, L):
    r"""Returns 3nu oscillation vacuum probabilities, std. computation.

    Returns the probabilities for three-neutrino oscillations in vacuum,
    computed using the standard analytical expression of the
    probabilities.

    Parameters
    ----------
    U : list
        3x3 PMNS complex mixing matrix.
    D21 : float
        Mass-squared difference Delta m^2_21.
    D31 : float
        Mass-squared difference Delta m^2_31.
    energy : float
        Neutrino energy.
    L : float
        Baseline.

    Returns
    -------
    list
        List of probabilities [Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm,
        Ptt].
    """
    D32 = D31-D21
    arg21 = D21*L/energy/2.0
    arg31 = D31*L/energy/2.0
    arg32 = D32*L/energy/2.0
    s21 = sin(arg21)
    s31 = sin(arg31)
    s32 = sin(arg32)
    ss21 = pow(sin(arg21/2.0), 2.0)
    ss31 = pow(sin(arg31/2.0), 2.0)
    ss32 = pow(sin(arg32/2.0), 2.0)

    # Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt
    prob = [delta(alpha, beta) \
            - 4.0 * ( J(U, alpha, beta, 1, 0).real*ss21
                    + J(U, alpha, beta, 2, 0).real*ss31
                    + J(U, alpha, beta, 2, 1).real*ss32 ) \
            + 2.0 * ( J(U, alpha, beta, 1, 0).imag*s21
                    + J(U, alpha, beta, 2, 0).imag*s31
                    + J(U, alpha, beta, 2, 1).imag*s32 ) \
            for alpha in [0,1,2] for beta in [0,1,2]]

    return prob


def hamiltonian_3nu_matter(h_vacuum_energy_independent, energy, VCC):
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

    h_matter = cp.deepcopy(h_vacuum_energy_independent)
    h_matter = np.multiply(1.0/energy, h_matter)

    # Add the matter potential to the ee term to find the matter
    # Hamiltonian
    h_matter[0][0] += VCC

    return h_matter


def hamiltonian_3nu_nsi(h_vacuum_energy_independent, energy, VCC, eps):
    r"""Returns the three-neutrino Hamiltonian for oscillations w/ NSI.

    Computes and returns the 3x3 complex three-neutrino Hamiltonian for
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
        Vector of NSI strength parameters: eps = eps_ee, eps_em, eps_et,
        eps_mm, eps_mt, eps_tt.

    Returns
    -------
    list
        Hamiltonian 3x3 matrix.
    """
    h_nsi = cp.deepcopy(h_vacuum_energy_independent)
    h_nsi = np.multiply(1.0/energy, h_nsi)

    eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt = eps

    h_nsi[0][0] += VCC*(1.0+eps_ee)
    h_nsi[0][1] += VCC*eps_em
    h_nsi[0][2] += VCC*eps_et
    h_nsi[1][0] += VCC*np.conj(eps_em)
    h_nsi[1][1] += VCC*eps_mm
    h_nsi[1][2] += VCC*eps_mt
    h_nsi[2][0] += VCC*np.conj(eps_et)
    h_nsi[2][1] += VCC*np.conj(eps_mt)
    h_nsi[2][2] += VCC*eps_tt

    return h_nsi


def hamiltonian_3nu_liv(h_vacuum_energy_independent, energy, sxi12, sxi23,
                        sxi13, dxiCP, b1, b2, b3, Lambda):
    r"""Returns the three-neutrino Hamiltonian for oscillations w/ LIV.

    Computes and returns the 3x3 complex three-neutrino Hamiltonian for
    oscillations in a CPT-odd Lorentz invariance-violating background.

    Parameters
    ----------
    h_vacuum_energy_independent : list
        Energy-independent part of the two-neutrino Hamiltonian for
        oscillations in vacuum.  This is computed by the routine
        hamiltonian_2nu_vacuum_energy_independent.
    energy : float
        Neutrino energy.
    sxi12 : float
        Sin(xi_12), with xi_12 the one of the mixing angles between the
        space of the eigenvectors of B3 and the flavor states.
    sxi23 : float
        Sin(xi_23), with xi_23 the one of the mixing angles between the
        space of the eigenvectors of B3 and the flavor states.
    sxi13 : float
        Sin(xi_12), with xi_13 the one of the mixing angles between the
        space of the eigenvectors of B3 and the flavor states.
    dciCP : float
        CP-violation angle of the LIV operator B3 [radian].
    b1 : float
        Eigenvalue b1 of the LIV operator B3.
    b2 : float
        Eigenvalue b2 of the LIV operator B3.
    b3 : float
        Eigenvalue b3 of the LIV operator B3.
    Lambda : float
        Energy scale of the LIV operator B2.

    Returns
    -------
    list
        Hamiltonian 3x3 matrix.
    """

    h_liv = cp.deepcopy(h_vacuum_energy_independent)
    h_liv = np.multiply(1.0/energy, h_liv)

    f = energy/Lambda
    # PMNS-like mixing matrix
    R = np.array(pmns_mixing_matrix(sxi12, sxi23, sxi13, dxiCP))
    # B matrix
    B = np.array([[b1, 0.0, 0.0], [0.0, b2, 0.0], [0.0, 0.0, b3]])
    # LIV term
    H = list(f*np.matmul(R, np.matmul(B, np.conj(matrix.transpose(R)))))

    h_liv[0][0] += H[0][0]
    h_liv[0][1] += H[0][1]
    h_liv[0][2] += H[0][2]
    h_liv[1][0] += H[1][0]
    h_liv[1][1] += H[1][1]
    h_liv[1][2] += H[1][2]
    h_liv[2][0] += H[2][0]
    h_liv[2][1] += H[2][1]
    h_liv[2][2] += H[2][2]

    return h_liv
