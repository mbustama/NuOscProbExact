# -*- coding: utf-8 -*-
r"""Oscillation probabilities across a sequence of adjacent slabs.

**NuOscProbExact** computes the evolution operator exactly for a
Hamiltonian that does not change along the trajectory.  A great many
interesting cases are *piecewise* constant instead: a neutrino crossing
the Earth, a beam through a layered detector, a castle-wall profile
built to enhance CP-violating effects.  This module handles those by
doing the only thing the exactness of the method allows --- solving each
slab exactly and multiplying the results.

For a trajectory divided into :math:`n` slabs, the slab :math:`k` having
Hamiltonian :math:`H_k` and width :math:`L_k`, the evolution operator is

.. math::

   U = U_n(L_n) \cdots U_2(L_2) U_1(L_1) ,

with the slab the neutrino meets first applied first --- rightmost,
since the operators act to the left on the initial state.  Each
:math:`U_k` is the exact SU(2), SU(3) or SU(4) expansion of
:mod:`oscprob2nu`,
:mod:`oscprob3nu` or :mod:`oscprob4nu`, so the only approximation in the
result is the one
the caller makes in choosing how finely to slice a continuously varying
profile.  Within each slab there is none.

The per-slab operators are evaluated in a single batched call, so the
cost of :math:`n` slabs is one vectorised evaluation plus :math:`n-1`
small matrix products rather than :math:`n` separate evaluations.

One thing carries over from the single-slab case and is easier to
overlook here.  The expansions return :math:`e^{-i H_0 L}`, with
:math:`H_0` the *traceless* part of the Hamiltonian, dropping the phase
:math:`e^{-i h_0 L}` that the trace contributes.  Each slab therefore
drops its own phase, and their product differs from
:math:`\prod_k e^{-i H_k L_k}` by the single scalar
:math:`\exp(i \sum_k h_0^{(k)} L_k)`.  That is still one overall phase,
so every probability is unaffected --- but a caller comparing the
returned operator against an independent matrix exponential must
compare against the traceless one, exactly as the single-slab tests do.

For the Earth specifically, :mod:`earth` builds the slabs for you from
the Preliminary Reference Earth Model.

Routine listings
----------------

    * evolution_operator_2nu_slabs - Two-flavor evolution operator
    * evolution_operator_3nu_slabs - Three-flavor evolution operator
    * evolution_operator_4nu_slabs - Four-flavor evolution operator
    * probabilities_2nu_slabs - Two-flavor probabilities
    * probabilities_3nu_slabs - Three-flavor probabilities
    * probabilities_4nu_slabs - Four-flavor probabilities
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

__all__ = ['evolution_operator_2nu_slabs', 'evolution_operator_3nu_slabs',
           'evolution_operator_4nu_slabs', 'probabilities_2nu_slabs',
           'probabilities_3nu_slabs', 'probabilities_4nu_slabs']

from typing import Tuple, Union

import numpy as np

import fastkernels
import oscprob2nu
import oscprob3nu
import oscprob4nu


def _check_slabs(
    hamiltonian_matrices: Union[list, np.ndarray],
    widths: Union[list, np.ndarray],
    n_flavors: int,
    caller: str
) -> Tuple[np.ndarray, np.ndarray]:
    r"""Validates and normalises a slab sequence.

    Parameters
    ----------
    hamiltonian_matrices : array_like
        Stack of Hamiltonians, of shape ``(n, n_flavors, n_flavors)``.
    widths : array_like
        Slab widths, of shape ``(n,)``.
    n_flavors : int
        Number of neutrino flavors, 2, 3, or 4.
    caller : str
        Name of the calling routine, used in error messages.

    Returns
    -------
    tuple of numpy.ndarray
        The Hamiltonians and widths as arrays.

    Raises
    ------
    ValueError
        If the two have different lengths, if either is empty, if the
        Hamiltonians are not square of the expected size, or if any
        width is negative.
    """
    h = np.asarray(hamiltonian_matrices, dtype=complex)
    w = np.asarray(widths, dtype=float)

    if h.ndim != 3 or h.shape[1:] != (n_flavors, n_flavors):
        raise ValueError(
            '%s: hamiltonian_matrices must have shape (n, %d, %d), got %s'
            % (caller, n_flavors, n_flavors, (h.shape,)))
    if w.ndim != 1:
        raise ValueError(
            '%s: widths must be one-dimensional, got shape %s'
            % (caller, (w.shape,)))
    if h.shape[0] != w.shape[0]:
        raise ValueError(
            '%s: got %d Hamiltonians but %d widths; there must be one '
            'width per slab' % (caller, h.shape[0], w.shape[0]))
    if h.shape[0] == 0:
        raise ValueError('%s: at least one slab is required' % caller)
    if np.any(w < 0.0):
        raise ValueError('%s: slab widths cannot be negative' % caller)

    return h, w


def _evolution_operator_slabs(
    hamiltonian_matrices: Union[list, np.ndarray],
    widths: Union[list, np.ndarray],
    n_flavors: int,
    caller: str
) -> np.ndarray:
    r"""Returns the evolution operator across a sequence of slabs.

    Parameters
    ----------
    hamiltonian_matrices : array_like
        Stack of Hamiltonians, of shape ``(n, n_flavors, n_flavors)``.
    widths : array_like
        Slab widths, of shape ``(n,)``.
    n_flavors : int
        Number of neutrino flavors, 2, 3, or 4.
    caller : str
        Name of the calling routine, used in error messages.

    Returns
    -------
    numpy.ndarray
        The evolution operator, of shape
        ``(n_flavors, n_flavors)``.
    """
    h, w = _check_slabs(hamiltonian_matrices, widths, n_flavors, caller)

    # The compiled path computes the operators *and* composes them in one
    # pass, so the stack is never materialised and the products never
    # leave registers.  Composing in Python was the largest single cost of
    # an Earth crossing once the operators themselves were compiled.
    if n_flavors == 2 and fastkernels.worthwhile_slabs(2, w.shape[0]):
        return fastkernels.slab_product_2nu_kernel(h, w)
    if n_flavors == 3 and fastkernels.worthwhile_slabs(3, w.shape[0]):
        return fastkernels.slab_product_3nu_kernel(h, w)
    if n_flavors == 4 and fastkernels.worthwhile_slabs(4, w.shape[0]):
        # The traceless part is what the expansion acts on, and what
        # `oscprob4nu` hands its own kernel; the dropped phase is per slab
        # and cancels in every probability, as the module docstring says.
        return fastkernels.slab_product_4nu_kernel(
            oscprob4nu._traceless_part(h), w, oscprob4nu.POLISH_ROOTS)

    # One batched call for all the slabs, rather than one call per slab:
    # the per-slab operators are independent of each other, and only their
    # product is not.
    if n_flavors == 2:
        u_slabs = np.asarray(oscprob2nu.evolution_operator_2nu(h, w))
    elif n_flavors == 3:
        u_slabs = np.asarray(oscprob3nu.evolution_operator_3nu(h, w))
    else:
        u_slabs = np.asarray(oscprob4nu.evolution_operator_4nu(h, w))

    # U = U_n ... U_1, the first slab crossed applied first.  The operator
    # is indexed (final, initial) and acts to the left on the initial
    # state, so the first slab ends up rightmost in the product.
    u_total = u_slabs[0]
    for k in range(1, u_slabs.shape[0]):
        u_total = u_slabs[k] @ u_total

    return u_total


def evolution_operator_2nu_slabs(
    hamiltonian_matrices: Union[list, np.ndarray],
    widths: Union[list, np.ndarray]
) -> np.ndarray:
    r"""Returns the two-flavor evolution operator across adjacent slabs.

    Returns the :math:`2\times2` evolution operator
    :math:`U_2 = U_2^{(n)}(L_n) \cdots U_2^{(1)}(L_1)` for a trajectory
    divided into slabs, each with its own constant Hamiltonian and
    width.  Each slab is solved exactly by the SU(2) expansion of
    :mod:`oscprob2nu`.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    hamiltonian_matrices : array_like
        Hamiltonians, of shape ``(n, 2, 2)``, one per slab and ordered
        along the trajectory, in units of eV.
    widths : array_like
        Slab widths, of shape ``(n,)``, in units of eV\ :sup:`-1`.  Use
        `globaldefs.CONV_KM_TO_INV_EV` to convert from km.

    Returns
    -------
    numpy.ndarray
        The evolution operator, of shape ``(2, 2)``, indexed
        ``(final, initial)``.

    Raises
    ------
    ValueError
        If the number of Hamiltonians and widths differ, if either is
        empty, if the Hamiltonians are not :math:`2\times2`, or if any
        width is negative.

    Examples
    --------
    .. jupyter-execute::

        import slabs

        import numpy as np
        H = np.array([[[0.0, 1.0], [1.0, 0.0]],
                      [[0.0, 0.5], [0.5, 0.0]]])
        U = slabs.evolution_operator_2nu_slabs(H, [0.3, 0.4])
        print('%.6f' % abs(U[0][0]))
    """
    return _evolution_operator_slabs(hamiltonian_matrices, widths, 2,
                                     'evolution_operator_2nu_slabs')


def evolution_operator_3nu_slabs(
    hamiltonian_matrices: Union[list, np.ndarray],
    widths: Union[list, np.ndarray]
) -> np.ndarray:
    r"""Returns the three-flavor evolution operator across adjacent slabs.

    Returns the :math:`3\times3` evolution operator
    :math:`U_3 = U_3^{(n)}(L_n) \cdots U_3^{(1)}(L_1)` for a trajectory
    divided into slabs, each with its own constant Hamiltonian and
    width.  Each slab is solved exactly by the SU(3) expansion of
    :mod:`oscprob3nu`.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    hamiltonian_matrices : array_like
        Hamiltonians, of shape ``(n, 3, 3)``, one per slab and ordered
        along the trajectory, in units of eV.
    widths : array_like
        Slab widths, of shape ``(n,)``, in units of eV\ :sup:`-1`.  Use
        `globaldefs.CONV_KM_TO_INV_EV` to convert from km.

    Returns
    -------
    numpy.ndarray
        The evolution operator, of shape ``(3, 3)``, indexed
        ``(final, initial)``.

    Raises
    ------
    ValueError
        If the number of Hamiltonians and widths differ, if either is
        empty, if the Hamiltonians are not :math:`3\times3`, or if any
        width is negative.

    Examples
    --------
    .. jupyter-execute::

        import slabs

        import numpy as np
        H = np.array([np.diag([1.0, 0.0, -1.0]),
                      np.diag([0.5, 0.0, -0.5])], dtype=complex)
        U = slabs.evolution_operator_3nu_slabs(H, [0.2, 0.3])
        print('%.6f' % abs(U[0][0]))
    """
    return _evolution_operator_slabs(hamiltonian_matrices, widths, 3,
                                     'evolution_operator_3nu_slabs')


def probabilities_2nu_slabs(
    hamiltonian_matrices: Union[list, np.ndarray],
    widths: Union[list, np.ndarray]
) -> Tuple[float, float, float, float]:
    r"""Returns the two-flavor probabilities across adjacent slabs.

    Returns :math:`P_{ee}, P_{e\mu}, P_{\mu e}, P_{\mu\mu}`, where
    :math:`P_{\alpha\beta} \equiv P(\nu_\alpha \to \nu_\beta) =
    |[U_2]_{\beta\alpha}|^2`, for a trajectory divided into slabs.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    hamiltonian_matrices : array_like
        Hamiltonians, of shape ``(n, 2, 2)``, one per slab and ordered
        along the trajectory, in units of eV.
    widths : array_like
        Slab widths, of shape ``(n,)``, in units of eV\ :sup:`-1`.

    Returns
    -------
    tuple of float
        The probabilities
        :math:`P_{ee}, P_{e\mu}, P_{\mu e}, P_{\mu\mu}`.

    Raises
    ------
    ValueError
        If the slab sequence is malformed; see
        `evolution_operator_2nu_slabs`.

    Examples
    --------
    .. jupyter-execute::

        import slabs

        import numpy as np
        H = np.array([[[0.0, 1.0], [1.0, 0.0]],
                      [[0.0, 0.5], [0.5, 0.0]]])
        Pee, Pem, Pme, Pmm = slabs.probabilities_2nu_slabs(H, [0.3, 0.4])
        print('%.6f  %.6f' % (Pee, Pem))
    """
    u = evolution_operator_2nu_slabs(hamiltonian_matrices, widths)

    # P_ab = |U_ba|^2: the evolution operator is indexed (final, initial)
    return (abs(u[0][0])**2.0, abs(u[1][0])**2.0,
            abs(u[0][1])**2.0, abs(u[1][1])**2.0)


def probabilities_3nu_slabs(
    hamiltonian_matrices: Union[list, np.ndarray],
    widths: Union[list, np.ndarray]
) -> Tuple[float, float, float, float, float, float, float, float, float]:
    r"""Returns the three-flavor probabilities across adjacent slabs.

    Returns :math:`P_{ee}, P_{e\mu}, P_{e\tau}, P_{\mu e}, P_{\mu\mu},
    P_{\mu\tau}, P_{\tau e}, P_{\tau\mu}, P_{\tau\tau}`, where
    :math:`P_{\alpha\beta} \equiv P(\nu_\alpha \to \nu_\beta) =
    |[U_3]_{\beta\alpha}|^2`, for a trajectory divided into slabs.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    hamiltonian_matrices : array_like
        Hamiltonians, of shape ``(n, 3, 3)``, one per slab and ordered
        along the trajectory, in units of eV.
    widths : array_like
        Slab widths, of shape ``(n,)``, in units of eV\ :sup:`-1`.

    Returns
    -------
    tuple of float
        The nine probabilities, with the initial flavor varying slowest.

    Raises
    ------
    ValueError
        If the slab sequence is malformed; see
        `evolution_operator_3nu_slabs`.

    Examples
    --------
    .. jupyter-execute::

        import slabs

        import numpy as np
        H = np.array([np.diag([1.0, 0.0, -1.0]),
                      np.diag([0.5, 0.0, -0.5])], dtype=complex)
        prob = slabs.probabilities_3nu_slabs(H, [0.2, 0.3])
        print('%.6f  %.6f' % (prob[0], prob[1]))
    """
    u = evolution_operator_3nu_slabs(hamiltonian_matrices, widths)

    # P_ab = |U_ba|^2: the evolution operator is indexed (final, initial)
    return (abs(u[0][0])**2.0, abs(u[1][0])**2.0, abs(u[2][0])**2.0,
            abs(u[0][1])**2.0, abs(u[1][1])**2.0, abs(u[2][1])**2.0,
            abs(u[0][2])**2.0, abs(u[1][2])**2.0, abs(u[2][2])**2.0)


def evolution_operator_4nu_slabs(
    hamiltonian_matrices: Union[list, np.ndarray],
    widths: Union[list, np.ndarray]
) -> np.ndarray:
    r"""Returns the four-flavor evolution operator across adjacent slabs.

    Returns the :math:`4\times4` evolution operator
    :math:`U_4 = U_4^{(n)}(L_n) \cdots U_4^{(1)}(L_1)` for a trajectory
    divided into slabs, each with its own constant Hamiltonian and
    width.  Each slab is solved exactly by the SU(4) expansion of
    :mod:`oscprob4nu`.

    This is what makes a 3+1 scenario propagable through layered matter:
    the sterile state's matter entry is constant within a slab like
    every other, so nothing about the composition changes at four
    flavors.  See :func:`hamiltonians4nu.hamiltonian_4nu_matter` for the
    entry itself, which is :math:`-V_{NC}` rather than zero.

    .. versionadded:: 1.11.0

    Parameters
    ----------
    hamiltonian_matrices : array_like
        Hamiltonians, of shape ``(n, 4, 4)``, one per slab and ordered
        along the trajectory, in units of eV.
    widths : array_like
        Slab widths, of shape ``(n,)``, in units of eV\ :sup:`-1`.  Use
        `globaldefs.CONV_KM_TO_INV_EV` to convert from km.

    Returns
    -------
    numpy.ndarray
        The evolution operator, of shape ``(4, 4)``, indexed
        ``(final, initial)``.

    Raises
    ------
    ValueError
        If the number of Hamiltonians and widths differ, if either is
        empty, if the Hamiltonians are not :math:`4\times4`, or if any
        width is negative.

    Examples
    --------
    .. jupyter-execute::

        import slabs

        import numpy as np
        H = np.array([np.diag([1.0, 0.0, -0.5, -0.5]),
                      np.diag([0.5, 0.0, -0.25, -0.25])], dtype=complex)
        U = slabs.evolution_operator_4nu_slabs(H, [0.2, 0.3])
        print('%.6f' % abs(U[0][0]))
    """
    return _evolution_operator_slabs(hamiltonian_matrices, widths, 4,
                                     'evolution_operator_4nu_slabs')


def probabilities_4nu_slabs(
    hamiltonian_matrices: Union[list, np.ndarray],
    widths: Union[list, np.ndarray]
) -> Tuple[float, ...]:
    r"""Returns the four-flavor probabilities across adjacent slabs.

    Returns the sixteen probabilities
    :math:`P_{\alpha\beta} \equiv P(\nu_\alpha \to \nu_\beta) =
    |[U_4]_{\beta\alpha}|^2`, ordered with the initial flavor varying
    slowest, for a trajectory divided into slabs.  With the fourth state
    read as sterile, the flavor order is
    :math:`(\nu_e, \nu_\mu, \nu_\tau, \nu_s)`.

    .. versionadded:: 1.11.0

    Parameters
    ----------
    hamiltonian_matrices : array_like
        Hamiltonians, of shape ``(n, 4, 4)``, one per slab and ordered
        along the trajectory, in units of eV.
    widths : array_like
        Slab widths, of shape ``(n,)``, in units of eV\ :sup:`-1`.

    Returns
    -------
    tuple of float
        The sixteen probabilities, with the initial flavor varying
        slowest.

    Raises
    ------
    ValueError
        If the slab sequence is malformed; see
        `evolution_operator_4nu_slabs`.

    Examples
    --------
    .. jupyter-execute::

        import slabs

        import numpy as np
        H = np.array([np.diag([1.0, 0.0, -0.5, -0.5]),
                      np.diag([0.5, 0.0, -0.25, -0.25])], dtype=complex)
        prob = slabs.probabilities_4nu_slabs(H, [0.2, 0.3])
        print('%.6f  %.6f' % (prob[0], prob[1]))
    """
    u = evolution_operator_4nu_slabs(hamiltonian_matrices, widths)

    # P_ab = |U_ba|^2: the evolution operator is indexed (final, initial)
    return tuple(abs(u[beta][alpha])**2.0
                 for alpha in range(4) for beta in range(4))
