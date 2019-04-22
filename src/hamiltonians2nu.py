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


def mixing_matrix_2nu(sth):

    cth = sqrt(1.0-sth*sth)

    U00 = cth
    U01 = sth
    U10 = -sth
    U11 = cth

    return [[U00,U01],[U10,U11]]


def hamiltonian_2nu_vacuum_energy_independent(sth, Dm2,
    compute_matrix_multiplication=False):

    th = np.arcsin(sth)
    c2th = cos(2.0*th)
    s2th = sin(2.0*th)

    f = 1./2.

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



########################################################################
# Oscillations in matter
########################################################################


def hamiltonian_2nu_matter(h_vacuum_energy_independent, energy, VCC):

    h_matter = cp.deepcopy(h_vacuum_energy_independent)
    h_matter = np.multiply(1.0/energy, h_matter)

    # Add the matter potential to the ee term to find the matter Hamiltonian
    h_matter[0][0] += VCC

    return h_matter



########################################################################
# Oscillations in matter with non-standard interactions
########################################################################


def hamiltonian_2nu_nsi(h_vacuum_energy_independent, energy, VCC, eps):

    h_nsi = cp.deepcopy(h_vacuum_energy_independent)
    h_nsi = np.multiply(1.0/energy, h_nsi)

    eps_ee, eps_em, eps_mm = eps

    h_nsi[0][0] += VCC*(1.0+eps_ee)
    h_nsi[0][1] += VCC*eps_em
    h_nsi[1][0] += VCC*np.conj(eps_em)
    h_nsi[1][1] += VCC*eps_mm

    return h_nsi



########################################################################
# Oscillations in a Lorentz-violating background
########################################################################

def hamiltonian_2nu_liv(h_vacuum_energy_independent, energy, b1, b2, sxi,
    Lambda):

    h_liv = cp.deepcopy(h_vacuum_energy_independent)
    h_liv = np.multiply(1.0/energy, h_liv)

    f = energy/Lambda
    cxi = sqrt(1.0-sxi-sxi)

    # Add the matter potential to the ee term to find the matter
    # Hamiltonian
    h_liv[0][0] += f*(b1*cxi*cxi + b2*sxi*sxi)
    h_liv[0][1] += f*((-b1+b2)*cxi*sxi)
    h_liv[1][0] += f*((-b1+b2)*cxi*sxi)
    h_liv[1][1] += f*(b2*cxi*cxi + b1*sxi*sxi)

    return h_liv



