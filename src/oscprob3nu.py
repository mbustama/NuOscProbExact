# -*- coding: utf-8 -*-
r"""Compute the three-neutrino flavor-transition probabilities.

This module contains the routines needed to compute three-neutrino
flavor-transition probabilities for an arbitrary time-independent
:math:`3\times3` Hermitian Hamiltonian, using the SU(3) exponential
expansion described in [2]_.

The Hamiltonian is expanded in the basis of Gell-Mann matrices, whose
structure constants and :math:`d` tensor follow the conventions of [1]_,

.. math:: H = h_0 \mathbb{1} + h_k \lambda^k ,

and the time-evolution operator in the same basis,

.. math:: U_3(L) = u_0 \mathbb{1} + i u_k \lambda^k .

The term :math:`h_0` contributes only an overall phase and is dropped;
all routines therefore work with the traceless part of the Hamiltonian,
which leaves the oscillation probabilities unchanged.

`evolution_operator_3nu` and `probabilities_3nu` accept either a single
Hamiltonian and baseline or a stack of them, in which case the whole
stack is evaluated at once.  See the notes on `probabilities_3nu`.

Units
-----

The routines are unit-agnostic: they require only that the Hamiltonian
and the baseline be given in reciprocal units, so that the product
:math:`H L` is dimensionless.  Elsewhere in **NuOscProbExact** the
Hamiltonian is in eV and the baseline in eV\ :sup:`-1`; the module
:mod:`globaldefs` provides ``CONV_KM_TO_INV_EV`` to convert a baseline
in km into eV\ :sup:`-1`.

Routine listings
----------------

    * hamiltonian_3nu_coefficients - Returns the :math:`h_k`
    * tensor_d - Returns the SU(3) tensor :math:`d_{ijk}`
    * star - Returns the SU(3) star product :math:`(h \star h)_i`
    * su3_invariants - Returns the SU(3) invariants :math:`|h|^2, \langle h \rangle`
    * psi_roots - Returns the roots of the characteristic equation
    * evolution_operator_3nu_u_coefficients - Returns the :math:`u_k`
    * evolution_operator_3nu - Returns the evolution operator :math:`U_3`
    * probabilities_3nu - Returns the oscillation probabilities

References
----------

.. [1] A.J. MacFarlane, A. Sudbery, and P.H. Weisz, "On Gell-Mann's
   :math:`\lambda`-matrices, :math:`d`- and :math:`f`-tensors, octets,
   and parametrizations of SU(3)", Commun. Math. Phys. 11, 77 (1968).

.. [2] Mauricio Bustamante, "Exact neutrino oscillation probabilities
   with arbitrary time-independent Hamiltonians", arXiv:1904.12391.

Created: 2019/04/11 15:36
Last modified: 2026/07/31
"""

from __future__ import print_function

__version__ = "1.1"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

__all__ = ['hamiltonian_3nu_coefficients', 'tensor_d', 'star',
           'su3_invariants', 'psi_roots',
           'evolution_operator_3nu_u_coefficients', 'evolution_operator_3nu',
           'probabilities_3nu']

import cmath

import numpy as np


SQRT3 = np.sqrt(3.0)
r"""float: Module-level constant equal to :math:`\sqrt{3}`."""

SQRT3_INV = 1./np.sqrt(3.0)
r"""float: Module-level constant equal to :math:`1/\sqrt{3}`."""

NEG_HALF_SQRT3_INV = -SQRT3_INV/2.0
r"""float: Module-level constant equal to :math:`-1/(2\sqrt{3})`."""

DEGENERACY_TOL = 1.e-12
r"""float: Module-level constant.

Relative tolerance below which two latent roots :math:`\psi` are treated
as degenerate.  The general expression for the :math:`u_k` divides by
:math:`3\psi_m^2 - |h|^2`, which vanishes at a repeated root, so a
degenerate spectrum is handled by a separate, exact expression.
"""


def hamiltonian_3nu_coefficients(hamiltonian_matrix):
    r"""Returns the :math:`h_k` of the SU(3) expansion of the Hamiltonian.

    Computes the coefficients :math:`h_1, \ldots, h_8` of the SU(3)
    expansion :math:`H = h_0 \mathbb{1} + h_k \lambda^k` of the
    three-flavor Hamiltonian `hamiltonian_matrix`, which is assumed to
    be given in the flavor basis.  The coefficient :math:`h_0`
    contributes only an overall phase to the evolution operator and is
    not returned.

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Three-flavor Hamiltonian, given as the nested list
        ``[[H11, H12, H13], [H21, H22, H23], [H31, H32, H33]]``.  It
        must be Hermitian.

    Returns
    -------
    list of float
        The eight coefficients ``[h1, h2, h3, h4, h5, h6, h7, h8]``.
        They are real, because the Hamiltonian is Hermitian.

    Examples
    --------
    >>> hamiltonian_matrix = [[1.0+0.0j, 0.0+2.0j, 0.0-1.0j],
    ...                       [0.0-2.0j, 3.0+0.0j, 3.0+0.0j],
    ...                       [0.0+1.0j, 3.0-0.0j, -5.0+0.0j]]
    >>> h_coeffs = hamiltonian_3nu_coefficients(hamiltonian_matrix)
    >>> print('  '.join(['%.6f' % (h+0.0) for h in h_coeffs]))
    0.000000  -2.000000  -1.000000  0.000000  1.000000  3.000000  0.000000  4.041452
    """
    H11 = hamiltonian_matrix[0][0]
    H12 = hamiltonian_matrix[0][1]
    H13 = hamiltonian_matrix[0][2]
    H22 = hamiltonian_matrix[1][1]
    H23 = hamiltonian_matrix[1][2]
    H33 = hamiltonian_matrix[2][2]

    # h0 = (H11+H22+H33)/3.0 is not returned: it multiplies the identity
    # and so contributes only an overall phase to U3, which cancels in
    # the oscillation probabilities.
    h1 = np.real(H12)
    h2 = -np.imag(H12)
    h3 = np.real(H11-H22)/2.0
    h4 = np.real(H13)
    h5 = -np.imag(H13)
    h6 = np.real(H23)
    h7 = -np.imag(H23)
    h8 = np.real(H11+H22-2.0*H33)*SQRT3/6.0

    return [float(h1), float(h2), float(h3), float(h4), float(h5),
            float(h6), float(h7), float(h8)]


def tensor_d(i, j, k):
    r"""Returns the tensor :math:`d_{ijk}` of the SU(3) algebra.

    Returns the totally symmetric SU(3) tensor
    :math:`d_{ijk} = \frac{1}{4}\mathrm{Tr}
    (\{\lambda_i, \lambda_j\} \lambda_k)`, defined in [1]_.

    Parameters
    ----------
    i : int
        First index, in the range 0--7 (i.e., :math:`d_{ijk}` is indexed
        from zero, so that ``i = 0`` corresponds to :math:`d_{1jk}`).
    j : int
        Second index, in the range 0--7.
    k : int
        Third index, in the range 0--7.

    Returns
    -------
    float
        The value of :math:`d_{ijk}`.

    Raises
    ------
    IndexError
        If any index lies outside the range 0--7.

    References
    ----------

    .. [1] A.J. MacFarlane, A. Sudbery, and P.H. Weisz, "On Gell-Mann's
       :math:`\lambda`-matrices, :math:`d`- and :math:`f`-tensors,
       octets, and parametrizations of SU(3)", Commun. Math. Phys. 11,
       77 (1968).

    Examples
    --------
    >>> print('%.6f' % tensor_d(0, 0, 7))
    0.577350
    >>> print('%.6f' % tensor_d(0, 1, 2))
    0.000000
    """
    for index in (i, j, k):
        if not 0 <= index <= 7:
            raise IndexError(
                'tensor_d: index %r is outside the range 0-7' % index)

    ip1 = i+1
    jkp1 = (j+1, k+1)

    if (ip1 == 1):
        if jkp1 == (1, 8): return SQRT3_INV
        if jkp1 == (4, 6): return 0.5
        if jkp1 == (5, 7): return 0.5
        if jkp1 == (6, 4): return 0.5
        if jkp1 == (7, 5): return 0.5
        if jkp1 == (8, 1): return SQRT3_INV
        return 0.0
    elif (ip1 == 2):
        if jkp1 == (2, 8): return SQRT3_INV
        if jkp1 == (4, 7): return -0.5
        if jkp1 == (5, 6): return 0.5
        if jkp1 == (6, 5): return 0.5
        if jkp1 == (7, 4): return -0.5
        if jkp1 == (8, 2): return SQRT3_INV
        return 0.0
    elif (ip1 == 3):
        if jkp1 == (3, 8): return SQRT3_INV
        if jkp1 == (4, 4): return 0.5
        if jkp1 == (5, 5): return 0.5
        if jkp1 == (6, 6): return -0.5
        if jkp1 == (7, 7): return -0.5
        if jkp1 == (8, 3): return SQRT3_INV
        return 0.0
    elif (ip1 == 4):
        if jkp1 == (1, 6): return 0.5
        if jkp1 == (2, 7): return -0.5
        if jkp1 == (3, 4): return 0.5
        if jkp1 == (4, 3): return 0.5
        if jkp1 == (4, 8): return NEG_HALF_SQRT3_INV
        if jkp1 == (6, 1): return 0.5
        if jkp1 == (7, 2): return -0.5
        if jkp1 == (8, 4): return NEG_HALF_SQRT3_INV
        return 0.0
    elif (ip1 == 5):
        if jkp1 == (1, 7): return 0.5
        if jkp1 == (2, 6): return 0.5
        if jkp1 == (3, 5): return 0.5
        if jkp1 == (5, 3): return 0.5
        if jkp1 == (5, 8): return NEG_HALF_SQRT3_INV
        if jkp1 == (6, 2): return 0.5
        if jkp1 == (7, 1): return 0.5
        if jkp1 == (8, 5): return NEG_HALF_SQRT3_INV
        return 0.0
    elif (ip1 == 6):
        if jkp1 == (1, 4): return 0.5
        if jkp1 == (2, 5): return 0.5
        if jkp1 == (3, 6): return -0.5
        if jkp1 == (4, 1): return 0.5
        if jkp1 == (5, 2): return 0.5
        if jkp1 == (6, 3): return -0.5
        if jkp1 == (6, 8): return NEG_HALF_SQRT3_INV
        if jkp1 == (8, 6): return NEG_HALF_SQRT3_INV
        return 0.0
    elif (ip1 == 7):
        if jkp1 == (1, 5): return 0.5
        if jkp1 == (2, 4): return -0.5
        if jkp1 == (3, 7): return -0.5
        if jkp1 == (4, 2): return -0.5
        if jkp1 == (5, 1): return 0.5
        if jkp1 == (7, 3): return -0.5
        if jkp1 == (7, 8): return NEG_HALF_SQRT3_INV
        if jkp1 == (8, 7): return NEG_HALF_SQRT3_INV
        return 0.0
    else:  # ip1 == 8
        if jkp1 == (1, 1): return SQRT3_INV
        if jkp1 == (2, 2): return SQRT3_INV
        if jkp1 == (3, 3): return SQRT3_INV
        if jkp1 == (4, 4): return NEG_HALF_SQRT3_INV
        if jkp1 == (5, 5): return NEG_HALF_SQRT3_INV
        if jkp1 == (6, 6): return NEG_HALF_SQRT3_INV
        if jkp1 == (7, 7): return NEG_HALF_SQRT3_INV
        if jkp1 == (8, 8): return -SQRT3_INV
        return 0.0


# The d tensor is constant, so it is tabulated once, at import time, as a
# dense 8x8x8 array.  Evaluating the star product and the SU(3)
# invariants from this table, rather than by calling tensor_d inside
# nested loops, is what makes repeated probability evaluations fast.
_TENSOR_D = np.array([[[tensor_d(i, j, k) for k in range(8)]
                       for j in range(8)] for i in range(8)])


def star(i, h_coeffs):
    r"""Returns the SU(3) star product :math:`(h \star h)_i`.

    Returns the SU(3) star product
    :math:`(h \star h)_i = d_{ijk} h^j h^k`, summed over repeated
    indices.

    Parameters
    ----------
    i : int
        Index of the star product, in the range 0--7.
    h_coeffs : array_like
        Eight-component vector of coefficients :math:`h_k`, as returned
        by `hamiltonian_3nu_coefficients`.

    Returns
    -------
    float
        The star product :math:`(h \star h)_i`.

    Examples
    --------
    >>> print('%.6f' % star(0, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
    0.000000
    >>> print('%.6f' % star(7, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
    0.577350
    """
    h = np.asarray(h_coeffs, dtype=float)

    return float(h @ _TENSOR_D[i] @ h)


def su3_invariants(h_coeffs):
    r"""Returns the two SU(3) invariants, :math:`|h|^2` and
    :math:`\langle h \rangle`.

    Returns the two invariants of the SU(3) expansion,
    :math:`|h|^2 = h_i h_i` and
    :math:`\langle h \rangle = d_{ijk} h_i h_j h_k`.  They equal,
    respectively, :math:`\mathrm{Tr}(H_0^2)/2` and
    :math:`\mathrm{Tr}(H_0^3)/2`, with :math:`H_0` the traceless part of
    the Hamiltonian.

    Parameters
    ----------
    h_coeffs : array_like
        Eight-component vector of coefficients :math:`h_k`, as returned
        by `hamiltonian_3nu_coefficients`.

    Returns
    -------
    h2 : float
        The SU(3) invariant :math:`|h|^2`.
    h3 : float
        The SU(3) invariant :math:`\langle h \rangle`.

    Examples
    --------
    >>> h2, h3 = su3_invariants([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    >>> print('%.6f  %.6f' % (h2, h3))
    1.000000  0.000000
    """
    h = np.asarray(h_coeffs, dtype=float)

    # h2 = |h|^2
    h2 = float(h @ h)

    # h3 = <h> = h_i (h*h)_i, using the star product computed once
    h3 = float(h @ (_TENSOR_D @ h @ h))

    return h2, h3


def psi_roots(h2, h3):
    r"""Returns the three latent roots :math:`\psi`.

    Returns the three latent roots :math:`\psi` of the characteristic
    equation :math:`\psi^3 - |h|^2 \psi - \frac{2}{3}\langle h \rangle
    = 0`, which are the eigenvalues of minus the traceless part of the
    Hamiltonian.  The roots are independent of the baseline.

    Parameters
    ----------
    h2 : float
        The SU(3) invariant :math:`|h|^2`.
    h3 : float
        The SU(3) invariant :math:`\langle h \rangle`.

    Returns
    -------
    list of float
        The three roots ``[psi1, psi2, psi3]``.  They are real, because
        the Hamiltonian is Hermitian.

    Notes
    -----
    For a Hermitian Hamiltonian the argument of the arc cosine lies in
    :math:`[-1, 1]`; it is clipped to that interval so that round-off
    cannot produce spurious complex roots, which would spoil the
    unitarity of the evolution operator.  When :math:`|h|^2 = 0` the
    Hamiltonian is proportional to the identity and all three roots
    vanish.

    Examples
    --------
    >>> psi = psi_roots(1.0, 0.0)
    >>> print('  '.join(['%.6f' % (round(p, 9)+0.0) for p in sorted(psi)]))
    -1.000000  0.000000  1.000000
    """
    if h2 <= 0.0:
        return [0.0, 0.0, 0.0]

    pre = 2.0*np.sqrt(h2)*SQRT3_INV
    chi = np.arccos(np.clip(-SQRT3*h3*pow(h2, -1.5), -1.0, 1.0))

    return [float(pre*np.cos((chi+2.*np.pi*m)/3.0)) for m in [1, 2, 3]]


def _u_coefficients_3nu_single(h, h2, h3, L):
    r"""Returns the nine :math:`u_k` for one Hamiltonian and baseline.

    The shared core of the scalar and the vectorised paths: given the
    SU(3) coefficients and invariants of a *single* Hamiltonian, it
    returns the coefficients of :math:`U_3(L)`.  The vectorised path
    calls it only for the rare elements whose spectrum is degenerate.

    Parameters
    ----------
    h : numpy.ndarray
        The eight coefficients :math:`h_k`, as a real array.
    h2 : float
        The SU(3) invariant :math:`|h|^2`.
    h3 : float
        The SU(3) invariant :math:`\langle h \rangle`.
    L : float
        Baseline.

    Returns
    -------
    list of complex
        The nine coefficients ``[u0, u1, ..., u8]``.
    """
    if h2 <= 0.0:
        # The Hamiltonian is proportional to the identity: U3 = 1
        return [1.0+0.j] + [0.j]*8

    # [psi1, psi2, psi3]
    psi = psi_roots(h2, h3)

    # (h*h)_k, computed once: it does not depend on the baseline
    star_coeffs = _TENSOR_D @ h @ h

    # Find the closest pair of latent roots, to detect degeneracy
    scale = np.sqrt(h2)
    pairs = [(0, 1, 2), (0, 2, 1), (1, 2, 0)]
    a, b, c = min(pairs, key=lambda p: abs(psi[p[0]]-psi[p[1]]))

    if abs(psi[a]-psi[b]) <= DEGENERACY_TOL*scale:
        # Doubly degenerate root: the general expression would divide by
        # zero, so use the two-projector form instead.
        psi_deg = (psi[a]+psi[b])/2.0
        exp_deg = cmath.exp(1.j*L*psi_deg)
        exp_odd = cmath.exp(1.j*L*psi[c])
        weight = (exp_odd-exp_deg)/(psi_deg-psi[c])
        u0 = exp_deg + weight*psi_deg
        uk = [-1.j*weight*h[k] for k in range(0, 8)]

        return [u0]+uk

    # [e^{i*L*psi1}, e^{i*L*psi2}, e^{i*L*psi3}]
    exp_psi = [cmath.exp(1.j*L*x) for x in psi]

    u0 = sum(exp_psi)/3.
    uk = [1.j*sum([exp_psi[m]*(psi[m]*h[k]-star_coeffs[k])
                   / (3.*psi[m]*psi[m]-h2) for m in [0, 1, 2]])
          for k in range(0, 8)]

    # [u0, u1, u2, u3, u4, u5, u6, u7, u8]
    return [u0]+uk


def _hamiltonian_3nu_coefficients_batch(h_matrix):
    r"""Returns the :math:`h_k` for a stack of Hamiltonians.

    The vectorised counterpart of `hamiltonian_3nu_coefficients`.

    Parameters
    ----------
    h_matrix : numpy.ndarray
        Hamiltonians, of shape ``(..., 3, 3)``.

    Returns
    -------
    numpy.ndarray
        The coefficients, of shape ``(..., 8)``, real.
    """
    return np.stack([
        h_matrix[..., 0, 1].real,
        -h_matrix[..., 0, 1].imag,
        (h_matrix[..., 0, 0]-h_matrix[..., 1, 1]).real/2.0,
        h_matrix[..., 0, 2].real,
        -h_matrix[..., 0, 2].imag,
        h_matrix[..., 1, 2].real,
        -h_matrix[..., 1, 2].imag,
        (h_matrix[..., 0, 0]+h_matrix[..., 1, 1]
         - 2.0*h_matrix[..., 2, 2]).real*SQRT3/6.0,
    ], axis=-1)


def _u_coefficients_3nu_batch(h, L):
    r"""Returns the nine :math:`u_k` for a stack of Hamiltonians.

    Every step of the SU(3) expansion --- the star product, the two
    invariants, the latent roots, and the Lagrange interpolation --- is
    evaluated for the whole stack at once.

    Parameters
    ----------
    h : numpy.ndarray
        The coefficients :math:`h_k`, of shape ``(..., 8)``, real.
    L : numpy.ndarray
        Baselines, of shape ``(...)``, already broadcast against `h`.

    Returns
    -------
    numpy.ndarray
        The coefficients, of shape ``(..., 9)``, complex.

    Notes
    -----
    Degenerate spectra are handled exactly, as in the scalar path, but
    the branch cannot be taken elementwise inside a vectorised
    expression.  The general formula is therefore evaluated everywhere
    with the vanishing denominators replaced by one --- which produces a
    finite but meaningless value for those elements --- and the affected
    elements are then recomputed with `_u_coefficients_3nu_single`.
    Degeneracy is measure-zero among floating-point Hamiltonians, so the
    fallback loop is empty in essentially every real use.
    """
    # Everything up to here depends on the Hamiltonian alone, so it is
    # evaluated at the shape of `h` and only then broadcast against the
    # baselines.  Scanning one Hamiltonian over N baselines therefore
    # solves the characteristic equation once, not N times.
    star = np.einsum('ijk,...j,...k->...i', _TENSOR_D, h, h)
    h2 = np.einsum('...i,...i->...', h, h)
    h3 = np.einsum('...i,...i->...', h, star)

    positive = h2 > 0.0
    safe_h2 = np.where(positive, h2, 1.0)

    pre = 2.0*np.sqrt(safe_h2)*SQRT3_INV
    chi = np.arccos(np.clip(-SQRT3*h3*safe_h2**-1.5, -1.0, 1.0))
    m = np.array([1.0, 2.0, 3.0])
    psi = pre[..., None]*np.cos((chi[..., None]+2.0*np.pi*m)/3.0)
    psi = np.where(positive[..., None], psi, 0.0)

    denom = 3.0*psi*psi - h2[..., None]
    safe_denom = np.where(denom != 0.0, denom, 1.0)
    numer = psi[..., :, None]*h[..., None, :] - star[..., None, :]

    # The baselines enter only here
    exp_psi = np.exp(1.j*L[..., None]*psi)

    u0 = exp_psi.sum(-1)/3.0
    uk = 1.j*np.einsum('...m,...mk->...k', exp_psi/safe_denom, numer)
    u = np.concatenate([u0[..., None], uk], axis=-1)

    # Recompute the elements whose spectrum is degenerate, using the same
    # criterion as the scalar path: the closest pair of latent roots.
    # Degeneracy is a property of the Hamiltonian, so the test is made at
    # the shape of `h` and then broadcast over the baselines.
    gaps = np.stack([np.abs(psi[..., 0]-psi[..., 1]),
                     np.abs(psi[..., 0]-psi[..., 2]),
                     np.abs(psi[..., 1]-psi[..., 2])], axis=-1)
    degenerate = ((gaps.min(-1) <= DEGENERACY_TOL*np.sqrt(safe_h2))
                  | ~positive)
    if degenerate.any():
        full = u.shape[:-1]
        h_full = np.broadcast_to(h, full+(8,))
        h2_full = np.broadcast_to(h2, full)
        h3_full = np.broadcast_to(h3, full)
        L_full = np.broadcast_to(L, full)
        for idx in zip(*np.nonzero(np.broadcast_to(degenerate, full))):
            u[idx] = _u_coefficients_3nu_single(
                h_full[idx], float(h2_full[idx]), float(h3_full[idx]),
                float(L_full[idx]))

    return u


def _evolution_operator_3nu_batch(h_matrix, L):
    r"""Returns :math:`U_3(L)` for a stack of Hamiltonians and baselines.

    Parameters
    ----------
    h_matrix : array_like
        Hamiltonians, of shape ``(..., 3, 3)``.
    L : array_like
        Baselines, broadcastable against the leading axes of `h_matrix`.

    Returns
    -------
    numpy.ndarray
        The evolution operators, of shape ``(..., 3, 3)``, complex.
    """
    h_matrix = np.asarray(h_matrix, dtype=complex)
    L = np.asarray(L, dtype=float)

    # Check that the two broadcast against each other, and fail here with
    # a clear message rather than deep inside the expansion
    np.broadcast_shapes(h_matrix.shape[:-2], L.shape)

    h = _hamiltonian_3nu_coefficients_batch(h_matrix)
    u = _u_coefficients_3nu_batch(h, L)

    u0 = u[..., 0]
    u1, u2, u3, u4, u5, u6, u7, u8 = [u[..., k] for k in range(1, 9)]

    return np.stack([
        np.stack([u0+1.j*(u3+u8/SQRT3), 1.j*u1+u2, 1.j*u4+u5], axis=-1),
        np.stack([1.j*u1-u2, u0-1.j*(u3-u8/SQRT3), 1.j*u6+u7], axis=-1),
        np.stack([1.j*u4-u5, 1.j*u6-u7, u0-1.j*2.*u8/SQRT3], axis=-1),
    ], axis=-2)


def _is_batched(hamiltonian_matrix, L):
    r"""Returns whether the arguments describe a stack of problems.

    A single Hamiltonian is an ``n``-by-``n`` matrix and a single
    baseline is a scalar, so anything with more axes than that is a
    stack, and the vectorised path applies.

    This runs on every scalar call, so it is written to be cheap: an
    exact type check short-circuits the common case, and
    ``numpy.ndim`` --- which would convert a nested list to an array
    every time --- is reached only for an argument that is neither a
    plain Python number nor a NumPy array.
    """
    if type(L) is not float and type(L) is not int:
        if np.ndim(L) > 0:
            return True

    if type(hamiltonian_matrix) is np.ndarray:
        return hamiltonian_matrix.ndim > 2

    return isinstance(hamiltonian_matrix[0][0], (list, tuple, np.ndarray))


def evolution_operator_3nu_u_coefficients(hamiltonian_matrix, L):
    r"""Returns the coefficients :math:`u_0, \ldots, u_8`.

    Returns the nine coefficients :math:`u_0, \ldots, u_8` of the
    three-neutrino time-evolution operator :math:`U_3(L)` in its SU(3)
    exponential expansion,
    :math:`U_3 = u_0 \mathbb{1} + i u_k \lambda^k`.

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Three-flavor Hermitian Hamiltonian, given as the nested list
        ``[[H11, H12, H13], [H21, H22, H23], [H31, H32, H33]]``.
    L : float
        Baseline, in units reciprocal to those of the Hamiltonian.

    Returns
    -------
    list of complex
        The nine coefficients ``[u0, u1, ..., u8]``.

    Notes
    -----
    The general expression divides by :math:`3\psi_m^2 - |h|^2`, which
    vanishes when two latent roots coincide.  Two degenerate cases are
    therefore handled separately, and exactly:

    * :math:`|h|^2 = 0`, when the Hamiltonian is proportional to the
      identity and :math:`U_3 = \mathbb{1}`;
    * a doubly degenerate root :math:`\psi_a = \psi_b \neq \psi_c`, when
      the spectral decomposition reduces to a single projector and
      :math:`U_3 = e^{i\psi_a L}\mathbb{1} +
      (e^{i\psi_c L} - e^{i\psi_a L}) P_c`, with
      :math:`P_c = (h_k\lambda^k + \psi_a \mathbb{1})/(\psi_a - \psi_c)`.

    Examples
    --------
    >>> hamiltonian_matrix = [[1.0+0.0j, 0.0+2.0j, 0.0-1.0j],
    ...                       [0.0-2.0j, 3.0+0.0j, 3.0+0.0j],
    ...                       [0.0+1.0j, 3.0-0.0j, -5.0+0.0j]]
    >>> u_coeffs = evolution_operator_3nu_u_coefficients(hamiltonian_matrix,
    ...                                                  1.0)
    >>> print('%+.6f%+.6fj' % (u_coeffs[0].real, u_coeffs[0].imag))
    +0.621522-0.047327j
    """
    # [h1, h2, h3, h4, h5, h6, h7, h8]
    h_coeffs = hamiltonian_3nu_coefficients(hamiltonian_matrix)
    h = np.asarray(h_coeffs, dtype=float)

    # h2 = |h|^2, h3 = <h>
    h2, h3 = su3_invariants(h_coeffs)

    return _u_coefficients_3nu_single(h, h2, h3, L)


def evolution_operator_3nu(hamiltonian_matrix, L):
    r"""Returns the three-neutrino time-evolution operator.

    Returns the three-neutrino time-evolution operator :math:`U_3(L)` in
    its SU(3) exponential expansion
    :math:`U_3(L) = u_0 \mathbb{1} + i u_k \lambda^k`.  This is a
    :math:`3\times3` unitary matrix.

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Three-flavor Hermitian Hamiltonian, given as the nested list
        ``[[H11, H12, H13], [H21, H22, H23], [H31, H32, H33]]``, or a
        stack of them, of shape ``(..., 3, 3)``.
    L : float or array_like
        Baseline, in units reciprocal to those of the Hamiltonian, or an
        array of baselines broadcastable against the leading axes of
        `hamiltonian_matrix`.

    Returns
    -------
    list of list of complex or numpy.ndarray
        For a single Hamiltonian and baseline, the time-evolution
        operator :math:`U_3(L)` --- a :math:`3\times3` unitary complex
        matrix --- as a nested list.  If either argument is a stack, an
        array of shape ``(..., 3, 3)``.

    See Also
    --------
    probabilities_3nu : Returns the probabilities built from this matrix.

    Examples
    --------
    >>> hamiltonian_matrix = [[1.0+0.0j, 0.0+2.0j, 0.0-1.0j],
    ...                       [0.0-2.0j, 3.0+0.0j, 3.0+0.0j],
    ...                       [0.0+1.0j, 3.0-0.0j, -5.0+0.0j]]
    >>> U3 = evolution_operator_3nu(hamiltonian_matrix, 1.0)
    >>> for row in U3:
    ...     print('  '.join(['%+.6f%+.6fj' % (z.real+0.0, z.imag+0.0)
    ...                      for z in row]))
    +0.546090-0.496423j  -0.600964-0.114920j  -0.278885+0.056655j
    +0.600964+0.114920j  +0.430462+0.614384j  -0.171381+0.183031j
    +0.278885-0.056655j  -0.171381+0.183031j  +0.888015-0.259943j
    """
    if _is_batched(hamiltonian_matrix, L):
        return _evolution_operator_3nu_batch(hamiltonian_matrix, L)

    u0, u1, u2, u3, u4, u5, u6, u7, u8 = \
        evolution_operator_3nu_u_coefficients(hamiltonian_matrix, L)

    return [
        [u0+1.j*(u3+u8/SQRT3), 1.j*u1+u2, 1.j*u4+u5],
        [1.j*u1-u2, u0-1.j*(u3-u8/SQRT3), 1.j*u6+u7],
        [1.j*u4-u5, 1.j*u6-u7, u0-1.j*2.*u8/SQRT3]
    ]


def probabilities_3nu(hamiltonian_matrix, L):
    r"""Returns the three-neutrino oscillation probabilities.

    Returns the three-neutrino flavor-transition probabilities
    :math:`P_{ee}, P_{e\mu}, P_{e\tau}, P_{\mu e}, P_{\mu\mu},
    P_{\mu\tau}, P_{\tau e}, P_{\tau\mu}, P_{\tau\tau}`, where
    :math:`P_{\alpha\beta} \equiv P(\nu_\alpha \to \nu_\beta)
    = |[U_3]_{\beta\alpha}|^2`.

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Three-flavor Hermitian Hamiltonian, given as the nested list
        ``[[H11, H12, H13], [H21, H22, H23], [H31, H32, H33]]``, or a
        stack of them, of shape ``(..., 3, 3)``.
    L : float or array_like
        Baseline, in units reciprocal to those of the Hamiltonian, or an
        array of baselines broadcastable against the leading axes of
        `hamiltonian_matrix`.

    Returns
    -------
    tuple of float or numpy.ndarray
        For a single Hamiltonian and baseline, the nine probabilities
        ``(Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt)`` as a tuple,
        ordered with the initial flavor varying slowest.  If either
        argument is a stack, an array of shape ``(..., 9)`` in the same
        order, with the leading axes given by broadcasting the two
        arguments together.

    Notes
    -----
    Passing arrays evaluates the whole stack at once, which is between
    one and two orders of magnitude faster than calling this routine in
    a Python loop.  The two common scans both broadcast naturally:

    * *versus baseline*, with one Hamiltonian and an array of
      baselines.  The characteristic equation depends only on the
      Hamiltonian, so it is solved once for the whole scan;
    * *versus energy*, with an array of Hamiltonians --- one per energy,
      since :math:`H \propto 1/E` --- and a single baseline.

    An oscillogram is the outer combination of the two, obtained by
    giving the Hamiltonians and the baselines separate axes, e.g.
    ``probabilities_3nu(H[:, None, :, :], L[None, :])``, which returns
    an array of shape ``(len(H), len(L), 9)``.

    Examples
    --------
    >>> hamiltonian_matrix = [[1.0+0.0j, 0.0+2.0j, 0.0-1.0j],
    ...                       [0.0-2.0j, 3.0+0.0j, 3.0+0.0j],
    ...                       [0.0+1.0j, 3.0-0.0j, -5.0+0.0j]]
    >>> prob = probabilities_3nu(hamiltonian_matrix, 1.0)
    >>> print('  '.join(['%.6f' % p for p in prob[0:3]]))
    0.544650  0.374364  0.080986
    >>> print('  '.join(['%.6f' % p for p in prob[3:6]]))
    0.374364  0.562764  0.062872
    >>> print('  '.join(['%.6f' % p for p in prob[6:9]]))
    0.080986  0.062872  0.856142
    """
    if _is_batched(hamiltonian_matrix, L):
        U = _evolution_operator_3nu_batch(hamiltonian_matrix, L)
        # P_ab = |U_ba|^2: the evolution operator is indexed
        # (final, initial), the probabilities (initial, final)
        prob = np.abs(U)**2.
        return np.swapaxes(prob, -1, -2).reshape(prob.shape[:-2]+(9,))

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
