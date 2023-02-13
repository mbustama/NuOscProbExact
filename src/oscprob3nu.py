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
    * evolution_operator_3nu_u_coefficients - Returns the :math:`u_k`
    * evolution_operator_3nu - Returns evolution operator :math:`U_3`
    * probabilities_3nu - Returns the oscillation probabilities

References
----------

.. [1] A.J. MacFarlane, A. Sudbery, and P.H. Weisz, "On Gell-Mann's
   :math:`\lambda`-matrices, :math:`d`- and math`f`-tensors, octets, and
   parametrizations of SU(3)", Commun. Math. Phys. 11, 77 (1968).

.. [2] Mauricio Bustamante, "NuOscProbExact: a general-purpose code to compute
   exact two-flavor and three-flavor neutrino oscillation probabilities", 
   arXiv:1904.12391.

Created: 2019/04/11 15:36 (Mauricio Bustamante)
Last modified: 2023/02/13 13:44 (Mauricio Bustamante)
"""


__version__ = "1.0"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


from numpy import *
import numpy as np
import cmath
import cmath as cmath


SQRT3 = sqrt(3.0)
r"""float: Module-level constant

Constant equal to sqrt(3.0).
"""

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
    list
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
    h8 = (H11+H22-2.0*H33)*SQRT3/6.0

    return [h1, h2, h3, h4, h5, h6, h7, h8]


def tensor_d(i, j, k):
    r"""Returns the tensor d_ijk of the SU(3) algebra.

    Returns the SU(3) tensor d_ijk.

    Parameters
    ----------
    i : int
        First index.
    j : int
        Second index.
    k : int
        Third index.

    Returns
    -------
    float
        Value of the tensor d_ijk.
    """
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
    r"""Returns the SU(3) star oroduct (h*h)_i.

    Returns the SU(3) star product (h*h)_i = d_ijk*h^j*h^k (summed over
    repeated indices).

    Parameters
    ----------
    i : int
        Index of the star product.
    h_coeffs : array_like
        Eight-component vector.

    Returns
    -------
    float
        Star product (h*h)_i.
    """
    res = sum([tensor_d(i,j,k)*h_coeffs[j]*h_coeffs[k]
            for j in range(0,8) for k in range(0,8)])

    return res


def su3_invariants(h_coeffs):
    r"""Returns the two SU(3) invariants, |h|^2 and <h>.

    Returns the two SU(3) invariants computed from the coefficients
    h_coeffs.

    Parameters
    ----------
    h_coeffs : array_like
        Eight-component vector.

    Returns
    -------
    h2 : float
        SU(3) invariant |h|^2.
    h3 : float
        SU(3) invariant <h>.
    """
    # h2 = |h|^2
    h2 = sum([h*h for h in h_coeffs])

    # h3 = <h>
    h3 = sum([tensor_d(i,j,k)*h_coeffs[i]*h_coeffs[j]*h_coeffs[k]
            for i in range(0,8) for j in range(0,8) for k in range(0,8)])

    return h2, h3


def psi_roots(h2, h3):
    r"""Returns the roots psi.

    Returns the three latent roots psi of the characteristic equation of
    -h_h*lambda^k.  These roots are L-independent.

    Parameters
    ----------
    h2 : float
        SU(3) invariant |h|^2.
    h3 : float
        SU(3) invariant <h>.

    Returns
    -------
    roots : list
        The three roots [psi1, psi2, psi3].
    """
    pre = 2.0*sqrt(h2)*SQRT3_INV
    chi = cmath.acos(-SQRT3*h3*pow(h2,-1.5))

    roots = [pre*cmath.cos((chi+2.*np.pi*m)/3.0) for m in [1,2,3]]

    return roots


def evolution_operator_3nu_u_coefficients(hamiltonian_matrix, L):
    r"""Returns coefficients u0, ..., u8 of the 3nu evolution operator.

    Returns the nine coefficients u0, ..., u8 of the three-neutrino
    time-evolution operator U3(L) in its SU(3) exponential expansion,
    i.e., U3 = u0*I + i*u_k*lambda^k.

    Parameters
    ----------
    hamiltonian_matrix : list
        3x3 Hamiltonian, [[H11,H12,H13],[H21,H22,H23],[H31,H32,H33]].
    L : float
        Baseline.

    Returns
    -------
    list
        The nine coefficients [u0, u1, u2, u3, u4, u5, u6, u7, u8].

    Example
    -------
    >>> hamiltonian_matrix = [
    ...                   [1.0+0.0j, 0.0+2.0j, 0.0-1.0j],
    ...                   [0.0-2.0j, 3.0+0.0j, 3.0+0.0j],
    ...                   [0.0+1.0j, 3.0-0.0j, -5.0+0.0j]
    ... ]
    >>> L = 1.0
    >>> u_coeffs = \
    ...     evolution_operator_3nu_u_coefficients(hamiltonian_matrix, L)
    >>> print(u_coeffs)
    [(0.02931819348583329-0.05379039810587031j), 0j,
     (-0.37401185730706943+0.5232654591567741j),
     (-0.22397496123471491-0.18194182213957913j), 0j,
     (0.17020182293481648-0.4632575258138263j),
     (0.6584815991291701+0.3845256294303862j), 0j,
     (-0.37629378652956447-0.17544272350921003j)]
    """
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
    r"""Returns the 3nu time-evolution operator in its SU(3) expanstion.

    Returns the three-neutrino time-evolution operator U3(L) in its
    exponential SU(3) expansion U3(L) = u0*I + i*u_k*lambda^k.  This is
    a 3x3 unitary matrix.

    Parameters
    ----------
    hamiltonian_matrix : list
        3x3 Hamiltonian, [[H11,H12,H13],[H21,H22,H23],[H31,H32,H33]].
    L : float
        Baseline.

    Returns
    -------
    list
        The U3(L) time-evolution operator, a 3x3 unitary complex matrix.

    Example
    -------
    >>> hamiltonian_matrix = [
    ...                   [1.0+0.0j, 0.0+2.0j, 0.0-1.0j],
    ...                   [0.0-2.0j, 3.0+0.0j, 3.0+0.0j],
    ...                   [0.0+1.0j, 3.0-0.0j, -5.0+0.0j]
    ... ]
    >>> L = 1.0
    >>> U3 = evolution_operator_3nu(hamiltonian_matrix, L)
    >>> print(U3)
    [[ 0.312551-0.495018j -0.374011+0.523265j  0.170201-0.463257j]
     [ 0.374011-0.523265j -0.051331-0.047068j -0.384525+0.65848j ]
     [-0.170201+0.463257j -0.384525+0.65848j  -0.173265+0.380716j]]
    """
    u0, u1, u2, u3, u4, u5, u6, u7, u8 = \
        evolution_operator_3nu_u_coefficients(hamiltonian_matrix, L)

    evolution_operator = [
                            [u0+1.j*(u3+u8/SQRT3), 1.j*u1+u2, 1.j*u4+u5],
                            [1.j*u1-u2, u0-1.j*(u3-u8/SQRT3), 1.j*u6+u7],
                            [1.j*u4-u5, 1.j*u6-u7, u0-1.j*2.*u8/SQRT3]
    ]

    return evolution_operator


def probabilities_3nu(hamiltonian_matrix, L):
    r"""Returns the 3nu oscillation probability.

    Returns the three-neutrino oscillation probabilities
    Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt.

    Parameters
    ----------
    hamiltonian_matrix : list
        3x3 Hamiltonian, [[H11,H12,H13],[H21,H22,H23],[H31,H32,H33]].
    L : float
        Baseline.

    Returns
    -------
    list
        Three-neutrino probabilities
        Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt.

    Example
    -------
    >>> hamiltonian_matrix = [
    ...                   [1.0+0.0j, 0.0+2.0j, 0.0-1.0j],
    ...                   [0.0-2.0j, 3.0+0.0j, 3.0+0.0j],
    ...                   [0.0+1.0j, 3.0-0.0j, -5.0+0.0j]
    ... ]
    >>> L = 1.0
    >>> Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = \
    ...     probabilities_3nu(hamiltonian_matrix, 1.0)
    >>> print(Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt)
    0.342732 0.413691 0.243576 0.413691 0.0048504 0.58145 0.243576
    0.58145 0.174965
    """
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
