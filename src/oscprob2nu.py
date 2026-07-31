# -*- coding: utf-8 -*-
r"""Compute the two-neutrino flavor-transition probabilities.

This module contains the routines needed to compute two-neutrino
flavor-transition probabilities for an arbitrary time-independent
:math:`2\times2` Hermitian Hamiltonian, using the SU(2) exponential
expansion described in [1]_.

The Hamiltonian is expanded in the basis of Pauli matrices,

.. math:: H = h_0 \mathbb{1} + h_k \sigma^k ,

and the time-evolution operator in the same basis,

.. math:: U_2(L) = u_0 \mathbb{1} + i u_k \sigma^k .

The term :math:`h_0` contributes only an overall phase and is dropped;
all routines therefore work with the traceless part of the Hamiltonian,
which leaves the oscillation probabilities unchanged.

`evolution_operator_2nu` and `probabilities_2nu` accept either a single
Hamiltonian and baseline or a stack of them, in which case the whole
stack is evaluated at once.

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

    * hamiltonian_2nu_coefficients - Returns the :math:`h_k`
    * modulus - Returns the modulus :math:`|h|` of a vector
    * evolution_operator_2nu_u_coefficients - Returns the :math:`u_k`
    * evolution_operator_2nu - Returns the evolution operator :math:`U_2`
    * probabilities_2nu - Returns the oscillation probabilities

References
----------

.. [1] Mauricio Bustamante, "Exact neutrino oscillation probabilities
   with arbitrary time-independent Hamiltonians", arXiv:1904.12391.

Created: 2019/04/20 19:07
Last modified: 2026/07/31
"""

from __future__ import print_function

__version__ = "1.1"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

__all__ = ['hamiltonian_2nu_coefficients', 'modulus',
           'evolution_operator_2nu_u_coefficients', 'evolution_operator_2nu',
           'probabilities_2nu']

import numpy as np


def hamiltonian_2nu_coefficients(hamiltonian_matrix):
    r"""Returns the :math:`h_k` of the SU(2) expansion of the Hamiltonian.

    Computes the coefficients :math:`h_1, h_2, h_3` of the SU(2)
    expansion :math:`H = h_0 \mathbb{1} + h_k \sigma^k` of the
    two-flavor Hamiltonian `hamiltonian_matrix`, which is assumed to be
    given in the flavor basis.  The coefficient :math:`h_0` contributes
    only an overall phase to the evolution operator and is not returned.

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Two-flavor Hamiltonian, given as the nested list
        ``[[H11, H12], [H21, H22]]``.  It must be Hermitian, i.e.
        ``H21 == conj(H12)`` and ``H11``, ``H22`` real.

    Returns
    -------
    list of float
        The three coefficients ``[h1, h2, h3]``.  They are real, because
        the Hamiltonian is Hermitian.

    See Also
    --------
    modulus : Returns the modulus :math:`|h|` of the returned vector.

    Examples
    --------
    >>> hamiltonian_matrix = [[1.0+0.0j, 0.0+2.0j],
    ...                       [0.0-2.0j, 3.0+0.0j]]
    >>> h_coeffs = hamiltonian_2nu_coefficients(hamiltonian_matrix)
    >>> print('%.6f  %.6f  %.6f' % tuple(h_coeffs))
    0.000000  -2.000000  -1.000000
    """
    H11 = hamiltonian_matrix[0][0]
    H12 = hamiltonian_matrix[0][1]
    H22 = hamiltonian_matrix[1][1]

    # h0 = (H11+H22)/2.0 is not returned: it multiplies the identity and
    # so contributes only an overall phase to U2, which cancels in the
    # oscillation probabilities.
    h1 = np.real(H12)
    h2 = -np.imag(H12)
    h3 = np.real(H11-H22)/2.0

    return [float(h1), float(h2), float(h3)]


def modulus(h_coeffs):
    r"""Returns the modulus :math:`|h|` of the vector of coefficients.

    Returns the modulus of the vector of coefficients :math:`h_k` of the
    SU(2) expansion of the two-neutrino Hamiltonian,
    :math:`|h| = \sqrt{|h_1|^2 + |h_2|^2 + |h_3|^2}`.

    Parameters
    ----------
    h_coeffs : array_like
        Three-component vector of coefficients :math:`h_k`, as returned
        by `hamiltonian_2nu_coefficients`.

    Returns
    -------
    float
        The modulus :math:`|h|`.

    Examples
    --------
    >>> print('%.6f' % modulus([0.0, -2.0, -1.0]))
    2.236068
    """
    return float(np.sqrt(sum([abs(h)**2.0 for h in h_coeffs])))


def _hamiltonian_2nu_coefficients_batch(h_matrix):
    r"""Returns the :math:`h_k` for a stack of Hamiltonians.

    The vectorised counterpart of `hamiltonian_2nu_coefficients`.

    Parameters
    ----------
    h_matrix : numpy.ndarray
        Hamiltonians, of shape ``(..., 2, 2)``.

    Returns
    -------
    numpy.ndarray
        The coefficients, of shape ``(..., 3)``, real.
    """
    return np.stack([
        h_matrix[..., 0, 1].real,
        -h_matrix[..., 0, 1].imag,
        (h_matrix[..., 0, 0]-h_matrix[..., 1, 1]).real/2.0,
    ], axis=-1)


def _u_coefficients_2nu_batch(h, L):
    r"""Returns the four :math:`u_k` for a stack of Hamiltonians.

    Parameters
    ----------
    h : numpy.ndarray
        The coefficients :math:`h_k`, of shape ``(..., 3)``, real.
    L : numpy.ndarray
        Baselines, of shape ``(...)``, already broadcast against `h`.

    Returns
    -------
    numpy.ndarray
        The coefficients, of shape ``(..., 4)``, real.
    """
    # |h| depends on the Hamiltonian alone; the baselines enter only in
    # the trigonometric factors below
    h_abs = np.sqrt(np.einsum('...i,...i->...', h, h))

    positive = h_abs > 0.0
    safe_h_abs = np.where(positive, h_abs, 1.0)

    phase = h_abs*L
    u0 = np.cos(phase)
    # The limit of -sin(|h|L)/|h| as |h| -> 0 is -L
    ss = np.where(positive, -np.sin(phase)/safe_h_abs,
                  -np.broadcast_to(L, np.shape(phase)))

    return np.concatenate([np.broadcast_to(u0, np.shape(ss))[..., None],
                           h*ss[..., None]], axis=-1)


def _evolution_operator_2nu_batch(h_matrix, L):
    r"""Returns :math:`U_2(L)` for a stack of Hamiltonians and baselines.

    Parameters
    ----------
    h_matrix : array_like
        Hamiltonians, of shape ``(..., 2, 2)``.
    L : array_like
        Baselines, broadcastable against the leading axes of `h_matrix`.

    Returns
    -------
    numpy.ndarray
        The evolution operators, of shape ``(..., 2, 2)``, complex.
    """
    h_matrix = np.asarray(h_matrix, dtype=complex)
    L = np.asarray(L, dtype=float)

    # Check that the two broadcast against each other, and fail here with
    # a clear message rather than deep inside the expansion
    np.broadcast_shapes(h_matrix.shape[:-2], L.shape)

    u = _u_coefficients_2nu_batch(
        _hamiltonian_2nu_coefficients_batch(h_matrix), L)
    u0, u1, u2, u3 = [u[..., k] for k in range(4)]

    return np.stack([
        np.stack([u0+1.j*u3, 1.j*u1+u2], axis=-1),
        np.stack([1.j*u1-u2, u0-1.j*u3], axis=-1),
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


def evolution_operator_2nu_u_coefficients(hamiltonian_matrix, L):
    r"""Returns the coefficients :math:`u_0, \ldots, u_3`.

    Returns the four coefficients :math:`u_0, \ldots, u_3` of the
    two-neutrino time-evolution operator :math:`U_2(L)` in its SU(2)
    exponential expansion,
    :math:`U_2 = u_0 \mathbb{1} + i u_k \sigma^k`.

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Two-flavor Hermitian Hamiltonian, given as the nested list
        ``[[H11, H12], [H21, H22]]``.
    L : float
        Baseline, in units reciprocal to those of the Hamiltonian.

    Returns
    -------
    list of float
        The four coefficients ``[u0, u1, u2, u3]``.  They are real,
        because the Hamiltonian is Hermitian; the factor :math:`i` that
        multiplies :math:`u_k` is part of the expansion, not of the
        coefficients.

    Notes
    -----
    When :math:`|h| = 0` the Hamiltonian is proportional to the
    identity, there is no flavor evolution, and the limit
    :math:`\sin(|h| L)/|h| \to L` is used.

    Examples
    --------
    >>> hamiltonian_matrix = [[1.0+0.0j, 0.0+2.0j],
    ...                       [0.0-2.0j, 3.0+0.0j]]
    >>> u_coeffs = evolution_operator_2nu_u_coefficients(hamiltonian_matrix,
    ...                                                  1.0)
    >>> print('%.6f  %.6f  %.6f  %.6f' % tuple(u+0.0 for u in u_coeffs))
    -0.617273  0.000000  0.703690  0.351845
    """
    # [h1, h2, h3]
    h_coeffs = hamiltonian_2nu_coefficients(hamiltonian_matrix)

    # h_abs = |h|
    h_abs = modulus(h_coeffs)

    u0 = np.cos(h_abs*L)
    # The limit of -sin(|h|L)/|h| as |h| -> 0 is -L
    ss = -L if h_abs == 0.0 else -np.sin(h_abs*L)/h_abs
    uk = [h_coeffs[k]*ss for k in range(0, 3)]

    # [u0, u1, u2, u3]
    return [u0]+uk


def evolution_operator_2nu(hamiltonian_matrix, L):
    r"""Returns the two-neutrino time-evolution operator.

    Returns the two-neutrino time-evolution operator :math:`U_2(L)` in
    its SU(2) exponential expansion
    :math:`U_2(L) = u_0 \mathbb{1} + i u_k \sigma^k`.  This is a
    :math:`2\times2` unitary matrix.

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Two-flavor Hermitian Hamiltonian, given as the nested list
        ``[[H11, H12], [H21, H22]]``, or a stack of them, of shape
        ``(..., 2, 2)``.
    L : float or array_like
        Baseline, in units reciprocal to those of the Hamiltonian, or an
        array of baselines broadcastable against the leading axes of
        `hamiltonian_matrix`.

    Returns
    -------
    list of list of complex or numpy.ndarray
        For a single Hamiltonian and baseline, the time-evolution
        operator :math:`U_2(L)` --- a :math:`2\times2` unitary complex
        matrix --- as a nested list.  If either argument is a stack, an
        array of shape ``(..., 2, 2)``.

    See Also
    --------
    probabilities_2nu : Returns the probabilities directly, more cheaply.

    Examples
    --------
    >>> hamiltonian_matrix = [[1.0+0.0j, 0.0+2.0j],
    ...                       [0.0-2.0j, 3.0+0.0j]]
    >>> U2 = evolution_operator_2nu(hamiltonian_matrix, 1.0)
    >>> for row in U2:
    ...     print('  '.join(['%+.6f%+.6fj' % (z.real+0.0, z.imag+0.0)
    ...                      for z in row]))
    -0.617273+0.351845j  +0.703690+0.000000j
    -0.703690+0.000000j  -0.617273-0.351845j
    """
    if _is_batched(hamiltonian_matrix, L):
        return _evolution_operator_2nu_batch(hamiltonian_matrix, L)

    u0, u1, u2, u3 = \
        evolution_operator_2nu_u_coefficients(hamiltonian_matrix, L)

    return [
        [u0+1.j*u3, 1.j*u1+u2],
        [1.j*u1-u2, u0-1.j*u3]
    ]


def probabilities_2nu(hamiltonian_matrix, L):
    r"""Returns the two-neutrino oscillation probabilities.

    Returns the two-neutrino flavor-transition probabilities
    :math:`P_{ee}, P_{e\mu}, P_{\mu e}, P_{\mu\mu}`, where
    :math:`P_{\alpha\beta} \equiv P(\nu_\alpha \to \nu_\beta)`.

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Two-flavor Hermitian Hamiltonian, given as the nested list
        ``[[H11, H12], [H21, H22]]``, or a stack of them, of shape
        ``(..., 2, 2)``.
    L : float or array_like
        Baseline, in units reciprocal to those of the Hamiltonian, or an
        array of baselines broadcastable against the leading axes of
        `hamiltonian_matrix`.

    Returns
    -------
    tuple of float or numpy.ndarray
        For a single Hamiltonian and baseline, the probabilities
        ``(Pee, Pem, Pme, Pmm)`` as a tuple.  If either argument is a
        stack, an array of shape ``(..., 4)`` in the same order.

    Notes
    -----
    Passing arrays evaluates the whole stack at once, which is roughly
    fifty times faster than calling this routine in a Python loop; see
    the notes on :func:`oscprob3nu.probabilities_3nu` for the two scans
    that broadcast naturally.

    The transition probability is

    .. math::
       P_{e\mu} = \frac{|h_1|^2 + |h_2|^2}{|h|^2} \sin^2(|h| L) ,

    i.e. :math:`|U_{\mu e}|^2 = u_1^2 + u_2^2`.  Both :math:`h_1` and
    :math:`h_2` contribute; :math:`h_2` vanishes only when the
    off-diagonal entry of the Hamiltonian is real.

    Examples
    --------
    >>> hamiltonian_matrix = [[1.0+0.0j, 0.0+2.0j],
    ...                       [0.0-2.0j, 3.0+0.0j]]
    >>> Pee, Pem, Pme, Pmm = probabilities_2nu(hamiltonian_matrix, 1.0)
    >>> print('%.6f  %.6f  %.6f  %.6f' % (Pee, Pem, Pme, Pmm))
    0.504821  0.495179  0.495179  0.504821
    """
    if _is_batched(hamiltonian_matrix, L):
        U = _evolution_operator_2nu_batch(hamiltonian_matrix, L)
        # P_ab = |U_ba|^2: the evolution operator is indexed
        # (final, initial), the probabilities (initial, final)
        prob = np.abs(U)**2.
        return np.swapaxes(prob, -1, -2).reshape(prob.shape[:-2]+(4,))

    # [h1, h2, h3]
    h_coeffs = hamiltonian_2nu_coefficients(hamiltonian_matrix)

    # h_abs = |h|
    h_abs = modulus(h_coeffs)

    if h_abs == 0.0:
        # The Hamiltonian is proportional to the identity: no flavor
        # transitions occur, whatever the baseline.
        Pem = 0.0
    else:
        Pem = (abs(h_coeffs[0])**2.0 + abs(h_coeffs[1])**2.0) / h_abs**2.0 \
            * pow(np.sin(h_abs*L), 2.0)

    Pme = Pem
    Pee = 1.0-Pem
    Pmm = 1.0-Pme

    return Pee, Pem, Pme, Pmm
