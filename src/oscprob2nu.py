# -*- coding: utf-8 -*-
r"""Compute the two-neutrino flavor-transition probability.

This module contains all the necessary routines to compute the
two-neutrino flavor-transition probabilities using the SU(2)
exponential expansion.

Routine listings
----------------

    * hamiltonian_2nu_coefficients - Returns coefficients of Hamiltonian
    * modulus - Returns the modulus of a vector
    * evolution_operator_2nu_u_coefficients - Returns the :math:`u_k`
    * evolution_operator_2nu - Returns evolution operator :math:`U_2`
    * probabilities_2nu - Returns the oscillation probabilities

References
----------

.. [1] Mauricio Bustamante, "Exact neutrino oscillation probabilities
   with arbitrary time-independent Hamiltonians", arXiv:1904.XXXXX.

Created: 2019/04/20 19:07
Last modified: 2019/04/22 20:33
"""

__version__ = "0.3"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


from numpy import *
import numpy as np
import cmath
import cmath as cmath

sigma_1 = [[0.0, 1.0], [1.0, 0.0]]
sigma_2 = [[0.0, -1.0j], [1.0j, 0]]
sigma_3 = [[1.0, 0.0], [0.0, -1.0]]
sigma = [sigma_1, sigma_2, sigma_3]
identity = [[1.0, 0.0], [0.0, 1.0]]
base = [identity, sigma_1, sigma_2, sigma_3]


def hamiltonian_2nu_coefficients(hamiltonian_matrix):
    r"""Returns the h_k of the SU(2)-expansion of the 2nu Hamiltonian.

    Computes the coefficients :math:`h_1, ..., h_3` in the SU(s)
    expansion of the provided three-flavor Hamiltonian
    `hamiltonian_matrix`, which is assumed to be given in the flavor
    basis.  The Hamiltonian is a :math:`2\times2` Hermitian matrix.

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Two-flavor Hamiltonian matrix, given as the list
        [[H11, H12], [H12*, H22], where the componentes Hij are complex
        numbers.

    Returns
    -------
    list
        List of coefficients [h1, h2, h3].  These are complex numbers,
        in general.

    Example
    -------
    >>> hamiltonian_matrix = [
    ...                   [1.0+0.0j, 0.0+2.0j],
    ...                   [0.0-2.0j, 3.0+0.0j]
    ... ]
    >>> h_coeffs = hamiltonian_2nu_coefficients(hamiltonian_matrix)
    >>> print(h_coeffs)
    [2j, 0.0, (-1+0j)]
    """
    H11 = hamiltonian_matrix[0][0]
    H12 = hamiltonian_matrix[0][1]
    H21 = hamiltonian_matrix[1][0]
    H22 = hamiltonian_matrix[1][1]

    # h0 = (H11+H22)/2.0  # Not used
    h1 = H12
    h2 = 0.0+0.0j
    h3 = (H11-H22)/2.0+0.0j

    return [h1, h2, h3]


def modulus(h_coeffs):
    r"""Returns the modulus of the vector of h_k coefficients, |h|.

    Returns the modulus of the vector of h_k coefficients of the SU(2)
    expansion of the two-neutrino Hamiltonian,
    |h| = sqrt(|h_1|^1+|h_2|^2+|h_3|^2)

    Parameters
    ----------
    h_coeffs : array_like
        Eight-component vector.

    Returns
    -------
    float
        Modulus |h|
    """
    # h_abs = |h|
    h_abs = sqrt(sum([abs(h)**2.0 for h in h_coeffs]))

    return h_abs


def evolution_operator_2nu_u_coefficients(hamiltonian_matrix, L):
    r"""Returns coefficients u0, ..., u3 of the 2nu evolution operator.

    Returns the four coefficients u0, ..., u3 of the two-neutrino
    time-evolution operator U2(L) in its SU(2) exponential expansion,
    i.e., U2 = u0*I + i*u_k*sigma^k.

    Parameters
    ----------
    hamiltonian_matrix : list
        2x2 Hamiltonian, [[H11,H12],[H21,H22]].
    L : float
        Baseline.

    Returns
    -------
    list
        The four coefficients [u0, u1, u2, u3]

    Example
    -------
    >>> hamiltonian_matrix = [
    ...                   [1.0+0.0j, 0.0+2.0j],
    ...                   [0.0-2.0j, 3.0+0.0j]
    ... ]
    >>> L = 1.0
    >>> u_coeffs = \
    ...     evolution_operator_2nu_u_coefficients(hamiltonian_matrix, L)
    >>> print(u_coeffs)
    [-0.6172728764571667, (-0-0.7036898157513979j), -0.0,
        (0.35184490787569894-0j)]
    """
    # [h1, h2, h3]
    h_coeffs = hamiltonian_2nu_coefficients(hamiltonian_matrix)

    # h_abs = |h|
    h_abs = modulus(h_coeffs)

    u0 = cos(h_abs*L)
    ss = -sin(h_abs*L)/h_abs
    uk = [h_coeffs[k]*ss for k in range(0,3)]

    # [u0, u1, u2, u3]
    u_coeffs = [u0]+uk

    return u_coeffs


def evolution_operator_2nu(hamiltonian_matrix, L):
    r"""Returns the 2nu time-evolution operator in its SU(2) expanstion.

    Returns the two-neutrino time-evolution operator U2(L) in its
    exponential SU(2) expansion U2(L) = u0*I + i*u_k*sigma^k.  This is
    a 2x2 unitary matrix.

    Parameters
    ----------
    hamiltonian_matrix : list
        2x2 Hamiltonian, [[H11,H12],[H21,H22]].
    L : float
        Baseline.

    Returns
    -------
    list
        The U2(L) time-evolution operator, a 2x2 unitary complex matrix.

    Example
    -------
    >>> hamiltonian_matrix = [
    ...                   [1.0+0.0j, 0.0+2.0j],
    ...                   [0.0-2.0j, 3.0+0.0j]
    ... ]
    >>> L = 1.0
    >>> U2 = evolution_operator_2nu(hamiltonian_matrix, L)
    >>> print(U2)
    [[ 0.312551-0.495018j -0.374011+0.523265j  0.170201-0.463257j]
     [ 0.374011-0.523265j -0.051331-0.047068j -0.384525+0.65848j ]
     [-0.170201+0.463257j -0.384525+0.65848j  -0.173265+0.380716j]]
    """
    u0, u1, u2, u3 = \
        evolution_operator_2nu_u_coefficients(hamiltonian_matrix, L)

    evolution_operator = [
                            [u0+1.j*u3, 1.j*u1+u2],
                            [1.j*u1-u2, u0-1.j*u3]
    ]

    return evolution_operator


def probabilities_2nu(hamiltonian_matrix, L):
    r"""Returns the 2nu oscillation probability.

    Returns the three-neutrino oscillation probabilities
    Pee, Pem, Pme, Pmm.

    Parameters
    ----------
    hamiltonian_matrix : list
        2x2 Hamiltonian, [[H11,H12],[H21,H22]].
    L : float
        Baseline.

    Returns
    -------
    list
        Two-neutrino probabilities Pee, Pem, Pme, Pmm

    Example
    -------
    >>> hamiltonian_matrix = [
    ...                   [1.0+0.0j, 0.0+2.0j],
    ...                   [0.0-2.0j, 3.0+0.0j]
    ... ]
    >>> L = 1.0
    >>> Pee, Pem, Pme, Pmm = \
    ...     probabilities_2nu(hamiltonian_matrix, 1.0)
    >>> print(Pee, Pem, Pme, Pmm)
    0.504820 0.495179 0.495179 0.504820
    """
    # [h1, h2, h3]
    h_coeffs = hamiltonian_2nu_coefficients(hamiltonian_matrix)

    # h_abs = |h|
    h_abs = modulus(h_coeffs)

    Pem = abs(h_coeffs[0])**2.0 / h_abs**2.0 * pow(sin(h_abs*L), 2.0)
    Pme = Pem
    Pee = 1.0-Pem
    Pmm = 1.0-Pme

    return Pee, Pem, Pme, Pmm




# hamiltonian_matrix = [
#                         [1.0+0.0j, 0.0+2.0j],
#                         [0.0-2.0j, 3.0+0.0j]
# ]


# h_coeffs = hamiltonian_2nu_coefficients(hamiltonian_matrix)
# print("h = ", h_coeffs)
# print()



# u_coeffs = evolution_operator_2nu_u_coefficients(hamiltonian_matrix, 1.0)
# print("u = ", u_coeffs)
# u0, u1, u2, u3 = u_coeffs
# print()


# u = np.array(evolution_operator_2nu(hamiltonian_matrix, 1.0))
# udag = np.conj(matrix.transpose(u))

# print(u)
# print(udag)
# print(np.matmul(u,udag))
# print(u @ udag)
# print()
# print(abs(u[0][0])**2.0)
# print(abs(u[0][1])**2.0)


# HH = np.array([np.multiply(h_coeffs[i], sigma[i]) for i in range(0,3)])
# print(HH.sum(axis=0))
# print()

# UU = np.array([ np.multiply(u_coeffs[0], identity),
#                 np.multiply(1.0j*u_coeffs[1], sigma_1),
#                 np.multiply(1.0j*u_coeffs[2], sigma_2),
#                 np.multiply(1.0j*u_coeffs[3], sigma_3)]).sum(axis=0)
# UUdag = np.conj(matrix.transpose(UU))
# print(UU)
# print(UU @ UUdag)
# print()


"""
Pee, Pem, Pme, Pmm = \
    probabilities_2nu(hamiltonian_matrix, 1.0)

print(Pee, Pem, Pme, Pmm)
print(Pee+Pem)
print(Pme+Pmm)
"""

