# -*- coding: utf-8 -*-
r"""Compute the three-neutrino flavor-transition probability.

This module contains all the necessary routines to compute the
three-neutrino flavor-transition probabilities using the SU(3)
exponential expansion.

Routine listings
----------------

    * hamiltonian_3nu_coefficients - Returns coefficients of Hamiltonian
    * tensor_d - Returns the value of the tensor :math:`d_{i,jk}`
    * star - Returns the SU(3) product :math:`(h * h)_k`
    * su3_invariants - Returns the SU(3) invariants :math:`|h|^2, <L>`
    * psi_roots - Returns the roots of the characteristic equation
    * evolution_operator_3nu__u_coefficients - Returns the :math:`u_k`
    * evolution_operator_3nu - Returns evolution operator :math:`U_3`
    * probabilities_3nu - Returns the oscillation probabilities

References
----------

.. [1] A.J. MacFarlane, A. Sudbery, and P.H. Weisz, "On Gell-Mann's
   :math:`\lambda`-matrices, :math:`d`- and math`f`-tensors, octets, and
   parametrizations of SU(3)", Commun. Math. Phys. 11, 77 (1968).

.. [2] Mauricio Bustamante, "Exact neutrino oscillation probabilities
   with arbitrary time-independent Hamiltonians", arXiv:1904.XXXXX.

Created: 2019/04/11 15:36
Last modified: 2019/04/15 22:21.
"""

__version__ = "0.3"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@nbi.ku.dk"


from numpy import *
import numpy as np
import cmath
import cmath as cmath


SQRT3_INV = 1./sqrt(3.0)
r"""float: Module-level constant

Constant equal to 1.0/sqrt(3.0).
"""

NEG_HALF_SQRT3_INV = -SQRT3_INV/2.0
r"""float: Module-level constant

Constant equal to -1.0/sqrt(3.0)/2.0.
"""


def hamiltonian_3nu_coefficients(hamiltonian_matrix):
    r"""Returns the h_k of the SU(3)-expansion of the 3nu Hamiltonian.

    Computes the coefficients :math:`h_1, ..., h_8` in the SU(3)
    expansion of the provided three-flavor Hamiltonian
    `hamiltonian_matrix`, which is assumed to be given in the flavor
    basis.  The Hamiltonian is a :math:`3\times3` Hermitian matrix.

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Three-flavor Hamiltonian matrix, given as the list
        [[H11, H12, H13], [H12*, H22, H23], [H13*, H23*, H33]], where
        the componentes Hij are complex numbers.

    Returns
    -------
    array_like
        List of coefficients [h1, h2, h3, h4, h5, h6, h7, h8].  These
        are complex numbers, in general.

    Example
    -------
    >>> hamiltonian_matrix = [
    ...                   [1.0+0.0j, 0.0+2.0j, 0.0-1.0j],
    ...                   [0.0-2.0j, 3.0+0.0j, 3.0+0.0j],
    ...                   [0.0+1.0j, 3.0-0.0j, -5.0+0.0j]
    ... ]
    >>> h_coeffs = hamiltonian_3nu_coefficients(hamiltonian_matrix)
    >>> print(h_coeffs)
    [0.0, -2.0, (-1+0j), 0.0, 1.0, 3.0, -0.0, (-1.7320508075688774+0j)]
    """
    H11 = hamiltonian_matrix[0][0]
    H12 = hamiltonian_matrix[0][1]
    H13 = hamiltonian_matrix[0][2]
    H21 = hamiltonian_matrix[1][0]
    H22 = hamiltonian_matrix[1][1]
    H23 = hamiltonian_matrix[1][2]
    H31 = hamiltonian_matrix[2][0]
    H32 = hamiltonian_matrix[2][1]
    H33 = hamiltonian_matrix[2][2]

    # h0 = (H11+H22+H33)/3.0  # Not used
    h1 = H12.real
    h2 = -H12.imag
    h3 = (H11-H22)/2.0
    h4 = H13.real
    h5 = -H13.imag
    h6 = H23.real
    h7 = -H23.imag
    h8 = (H11+H22-2.0*H33)*sqrt(3.0)/6.0

    return [h1, h2, h3, h4, h5, h6, h7, h8]


def tensor_d(i, j, k):

    ip1 = i+1
    jp1 = j+1
    kp1 = k+1
    jkp1 = (jp1, kp1)

    if (ip1 == 1):
        if jkp1 == (1,8): return SQRT3_INV
        if jkp1 == (4,6): return 0.5
        if jkp1 == (5,7): return 0.5
        if jkp1 == (6,4): return 0.5
        if jkp1 == (7,5): return 0.5
        if jkp1 == (8,1): return SQRT3_INV
        return 0.0
    elif (ip1 == 2):
        if jkp1 == (2,8): return SQRT3_INV
        if jkp1 == (4,7): return -0.5
        if jkp1 == (5,6): return 0.5
        if jkp1 == (6,5): return 0.5
        if jkp1 == (7,4): return -0.5
        if jkp1 == (8,2): return SQRT3_INV
        return 0.0
    elif (ip1 == 3):
        if jkp1 == (3,8): return SQRT3_INV
        if jkp1 == (4,4): return 0.5
        if jkp1 == (5,5): return 0.5
        if jkp1 == (6,6): return -0.5
        if jkp1 == (7,7): return -0.5
        if jkp1 == (8,3): return SQRT3_INV
        return 0.0
    elif (ip1 == 4):
        if jkp1 == (1,6): return 0.5
        if jkp1 == (2,7): return -0.5
        if jkp1 == (3,4): return 0.5
        if jkp1 == (4,3): return 0.5
        if jkp1 == (4,8): return NEG_HALF_SQRT3_INV
        if jkp1 == (6,1): return 0.5
        if jkp1 == (7,2): return -0.5
        if jkp1 == (8,4): return NEG_HALF_SQRT3_INV
        return 0.0
    elif (ip1 == 5):
        if jkp1 == (1,7): return 0.5
        if jkp1 == (2,6): return 0.5
        if jkp1 == (3,5): return 0.5
        if jkp1 == (5,3): return 0.5
        if jkp1 == (5,8): return NEG_HALF_SQRT3_INV
        if jkp1 == (6,2): return 0.5
        if jkp1 == (7,1): return 0.5
        if jkp1 == (8,5): return NEG_HALF_SQRT3_INV
        return 0.0
    elif (ip1 == 6):
        if jkp1 == (1,4): return 0.5
        if jkp1 == (2,5): return 0.5
        if jkp1 == (3,6): return -0.5
        if jkp1 == (4,1): return 0.5
        if jkp1 == (5,2): return 0.5
        if jkp1 == (6,3): return -0.5
        if jkp1 == (6,8): return NEG_HALF_SQRT3_INV
        if jkp1 == (8,6): return NEG_HALF_SQRT3_INV
        return 0.0
    elif (ip1 == 7):
        if jkp1 == (1,5): return 0.5
        if jkp1 == (2,4): return -0.5
        if jkp1 == (3,7): return -0.5
        if jkp1 == (4,2): return -0.5
        if jkp1 == (5,1): return 0.5
        if jkp1 == (7,3): return -0.5
        if jkp1 == (7,8): return NEG_HALF_SQRT3_INV
        if jkp1 == (8,7): return NEG_HALF_SQRT3_INV
        return 0.0
    elif (ip1 == 8):
        if jkp1 == (1,1): return SQRT3_INV
        if jkp1 == (2,2): return SQRT3_INV
        if jkp1 == (3,3): return SQRT3_INV
        if jkp1 == (4,4): return NEG_HALF_SQRT3_INV
        if jkp1 == (5,5): return NEG_HALF_SQRT3_INV
        if jkp1 == (6,6): return NEG_HALF_SQRT3_INV
        if jkp1 == (7,7): return NEG_HALF_SQRT3_INV
        if jkp1 == (8,8): return -SQRT3_INV
        return 0.0


def star(i, h_coeffs):

    res = sum([tensor_d(i,j,k)*h_coeffs[j]*h_coeffs[k]
            for j in range(0,8) for k in range(0,8)])

    return res


def su3_invariants(h_coeffs):

    # h2 = |h|^2
    h2 = sum([h*h for h in h_coeffs])

    # h3 = <h>
    h3 = sum([tensor_d(i,j,k)*h_coeffs[i]*h_coeffs[j]*h_coeffs[k]
            for i in range(0,8) for j in range(0,8) for k in range(0,8)])

    return h2, h3


def psi_roots(h2, h3):

    pre = 2.0*sqrt(h2)*SQRT3_INV
    chi = cmath.acos(-sqrt(3.0)*h3*pow(h2,-1.5))

    roots = [pre*cmath.cos((chi+2.*np.pi*m)/3.0) for m in [1,2,3]]

    return roots


def evolution_operator_3nu_u_coefficients(hamiltonian_matrix, L):

    # [h1, h2, h3, h4, h5, h6, h7, h8]
    h_coeffs = hamiltonian_3nu_coefficients(hamiltonian_matrix)

    # h2 = |h|^2, h3 = <h>
    h2, h3 = su3_invariants(h_coeffs)

    # [psi1, psi2, psi3]
    psi = psi_roots(h2, h3)

    # [e^{i*L*psi1}, e^{i*L*psi2}, e^{i*L*psi3}]
    exp_psi = [cmath.exp(1.j*L*x) for x in psi]

    u0 = sum([x for x in exp_psi])/3.
    uk = [ 1.j*sum([exp_psi[m]*(psi[m]*h_coeffs[k]-star(k,h_coeffs)) \
            /(3.*psi[m]*psi[m]-h2) for m in [0,1,2]]) for k in range(0,8)]

    # [u0, u1, u2, u3, u4, u5, u6, u7, u8]
    u_coeffs = [u0]+uk

    return u_coeffs


def evolution_operator_3nu(hamiltonian_matrix, L):

    u0, u1, u2, u3, u4, u5, u6, u7, u8 = \
        evolution_operator_3nu_u_coefficients(hamiltonian_matrix, L)

    evolution_operator = [
                            [u0+1.j*(u3+u8/sqrt(3.)), 1.j*u1+u2, 1.j*u4+u5],
                            [1.j*u1-u2, u0-1.j*(u3-u8/sqrt(3.)), 1.j*u6+u7],
                            [1.j*u4-u5, 1.j*u6-u7, u0-1.j*2.*u8/sqrt(3.)]
    ]

    return evolution_operator


def probabilities_3nu(hamiltonian_matrix, L):

    U = evolution_operator_3nu(hamiltonian_matrix, L)

    Pee = abs(U[0][0])**2.
    Pem = abs(U[1][0])**2.
    Pet = abs(U[2][0])**2.
    Pme = abs(U[0][1])**2.
    Pmm = abs(U[1][1])**2.
    Pmt = abs(U[2][1])**2.
    Pte = abs(U[0][2])**2.
    Ptm = abs(U[1][2])**2.
    Ptt = abs(U[2][2])**2.

    return Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt



# hamiltonian_matrix = [
#                         [1.0, -2.0j, 1.0j],
#                         [2.0j, 3.0, 3.0j],
#                         [-1.0j, -3.0j, 4.0]
# ]

# hamiltonian_matrix = [
#                         [1.0+0.0j, 0.0+2.0j, 0.0-1.0j],
#                         [0.0-2.0j, 3.0+0.0j, 3.0+0.0j],
#                         [0.0+1.0j, 3.0-0.0j, 5.0+0.0j]
# ]

# hamiltonian_matrix = [
#                         [1.0, 0.0, 0.0],
#                         [0.0, 3.0, 3.0],
#                         [0.0, 3.0, 5.0]
# ]

# h_coeffs = hamiltonian_3nu_coefficients(hamiltonian_matrix)
# print("h = ", h_coeffs)
# print()


# print(np.array(hamiltonian_matrix))
# print()

# HH = np.array([np.multiply(h_coeffs[i], lambda_matrix[i]) for i in range(0,8)])
# print(HH.sum(axis=0))
# print()

# quit()

# h2, h3 = su3_invariants(h_coeffs)
# print("h2 = ", h2)
# print("h3 = ", h3)
# print()

# # print([star(i, h_coeffs) for i in range(0,8)])

# phi = phi_roots(I2, I3)
# print(phi)


# u_coeffs = evolution_operator_3nu_u_coefficients(hamiltonian_matrix, 1.0)
# print("u = ", u_coeffs)
# print()


# u = np.array(evolution_operator_3nu(hamiltonian_matrix, 1.0))
# udag = np.conj(matrix.transpose(u))

# print(u)
# print(udag)
# print(np.matmul(u,udag))
# print(u @ udag)




"""
Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = \
    probabilities_3nu(hamiltonian_matrix, 1.0)

print(Pee+Pem+Pet)
print(Pme+Pmm+Pmt)
print(Pte+Ptm+Ptt)
"""


"""
hamiltonian = [
                [1.0+0.0j, 0.0+2.0j, 0.0-1.0j],
                [0.0-2.0j, 3.0+0.0j, 3.0+0.0j],
                [0.0+1.0j, 3.0-0.0j, 5.0+0.0j]
]

L = 1.0

Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = \
    probabilities_3nu(hamiltonian, L)

print(Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt)
"""
