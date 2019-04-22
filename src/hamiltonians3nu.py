#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = "0.1"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


"""
3nuhamiltonians.py:
    Routines to compute the three-neutrino flavor-transition probability
    in vacuum, matter, NSI, and LIV

Created: 2019/04/17 17:14
Last modified: 2019/04/17 17:14
"""


from numpy import *
import numpy as np
# from pylab import *
# from matplotlib import *
# import matplotlib as mpl
import cmath
import cmath as cmath
import copy as cp

import oscprob3nu
from globaldefs import *


########################################################################
# Oscillations in vacuum
########################################################################


def pmns_mixing_matrix(s12, s23, s13, dCP):

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

    if (a == b):
        return 1
    else:
        return 0


def J(U, alpha, beta, k, j):

    return np.conj(U[alpha][k])*U[beta][k]*U[alpha][j]*np.conj(U[beta][j])


def probabilities_3nu_std(U, D21, D31, energy_nu, l):

    D32 = D31-D21
    arg21 = D21*l/energy_nu/2.0
    arg31 = D31*l/energy_nu/2.0
    arg32 = D32*l/energy_nu/2.0
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



########################################################################
# Oscillations in matter
########################################################################


def hamiltonian_3nu_matter(h_vacuum_energy_independent, energy, VCC):

    h_matter = cp.deepcopy(h_vacuum_energy_independent)
    h_matter = np.multiply(1.0/energy, h_matter)

    # Add the matter potential to the ee term to find the matter
    # Hamiltonian
    h_matter[0][0] += VCC

    return h_matter



########################################################################
# Oscillations in matter with non-standard interactions
########################################################################


def hamiltonian_3nu_nsi(h_vacuum_energy_independent, energy, VCC, eps):

    h_nsi = cp.deepcopy(h_vacuum_energy_independent)
    h_nsi = np.multiply(1.0/energy, h_nsi)

    eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt = eps

    # Add the matter potential to the ee term to find the matter
    # Hamiltonian
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



########################################################################
# Oscillations in a Lorentz-violating background
########################################################################

def hamiltonian_3nu_liv(h_vacuum_energy_independent, energy, sxi12, sxi23,
                        sxi13, dxiCP, b1, b2, b3, Lambda):

    h_liv = cp.deepcopy(h_vacuum_energy_independent)
    h_liv = np.multiply(1.0/energy, h_liv)

    f = energy/Lambda
    # PMNS-like mixing matrix
    R = np.array(pmns_mixing_matrix(sxi12, sxi23, sxi13, dxiCP))
    # B matrix
    B = np.array([[b1, 0.0, 0.0], [0.0, b2, 0.0], [0.0, 0.0, b3]])
    # LIV term
    H = list(f*np.matmul(R, np.matmul(B, np.conj(matrix.transpose(R)))))

    # print(H)

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




