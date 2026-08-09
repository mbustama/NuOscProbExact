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

With the compiled backend this goes further: the
whole trajectory is one pass, each slab's operator computed and
multiplied into the running product in registers, so the stack is never
materialised and the :math:`n-1` products are never dispatched.  That is
worth between seven and two hundred times the NumPy path depending on
the flavor count and the number of slabs --- see
:data:`fastkernels.MIN_SLAB_BATCH`, which is why the threshold here is
one rather than `fastkernels.MIN_BATCH`.  Until 1.12.0 there was no
compiled path at all for this module: the backend had probability
kernels only, and composing operators cannot use one.

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

Chords that share their geometry and differ only in their Hamiltonians
--- an energy scan across a fixed profile --- can be passed together, as
an array of shape ``(..., n_slabs, n_flavors, n_flavors)`` against the
one set of widths they share.  Every routine here takes that form and
returns one result per chord, composing the whole batch in a single
pass rather than one chord at a time.

The widths are the limit of what that buys.  A batch shares them, so it
is the right tool for varying the Hamiltonian at fixed geometry and the
wrong one for varying the geometry, where each chord has its own slab
widths and there is nothing for one call to share.

For the Earth, do not build the batch by hand at all.  :mod:`earth`
generates the slabs from the Preliminary Reference Earth Model and
applies this same batching internally, and its routines take energies
and zenith angles on separate axes --- ``probabilities_3nu_earth(h,
energies[None, :], costhz[:, None])`` returns a whole oscillogram,
handling per-angle geometry that no single batch here could express.
Reach for the routines in this module when the profile is one
:mod:`earth` does not know about: a castle wall, a solar model, a
hand-built layer sequence, or an Earth chord carrying an extra term.

Routine listings
----------------

    * evolution_operator_2nu_slabs - Two-flavor evolution operator
    * evolution_operator_3nu_slabs - Three-flavor evolution operator
    * evolution_operator_4nu_slabs - Four-flavor evolution operator
    * probabilities_2nu_slabs - Two-flavor probabilities
    * probabilities_3nu_slabs - Three-flavor probabilities
    * probabilities_4nu_slabs - Four-flavor probabilities
    * probabilities_2nu_profile - Two flavors, across a varying profile
    * probabilities_3nu_profile - Three flavors, across a varying profile
    * probabilities_4nu_profile - Four flavors, across a varying profile
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

__all__ = ['N_SLABS_MAX',
           'evolution_operator_2nu_slabs', 'evolution_operator_3nu_slabs',
           'evolution_operator_4nu_slabs', 'probabilities_2nu_slabs',
           'probabilities_3nu_slabs', 'probabilities_4nu_slabs',
           'probabilities_2nu_profile', 'probabilities_3nu_profile',
           'probabilities_4nu_profile']

from typing import Callable, Optional, Tuple, Union

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


def _check_slabs_batch(
    hamiltonian_matrices: Union[list, np.ndarray],
    widths: Union[list, np.ndarray],
    n_flavors: int,
    caller: str
) -> Tuple[np.ndarray, np.ndarray]:
    r"""Validates a batch of chords sharing one set of slab widths.

    The batched counterpart of `_check_slabs`.  Every chord in the batch
    crosses the same geometry --- that is what an energy scan at fixed
    zenith angle is --- so there is one width per slab rather than one
    per slab per chord.

    Parameters
    ----------
    hamiltonian_matrices : array_like
        Stack of Hamiltonians, of shape
        ``(..., n_slabs, n_flavors, n_flavors)``, with at least one
        leading batch axis.
    widths : array_like
        Slab widths, of shape ``(n_slabs,)``.
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
        If the shapes disagree, if either is empty, if the Hamiltonians
        are not square of the expected size, or if any width is
        negative.
    """
    h = np.asarray(hamiltonian_matrices, dtype=complex)
    w = np.asarray(widths, dtype=float)

    if h.ndim < 4 or h.shape[-2:] != (n_flavors, n_flavors):
        raise ValueError(
            '%s: hamiltonian_matrices must have shape (..., n, %d, %d) with '
            'at least one leading batch axis, got %s'
            % (caller, n_flavors, n_flavors, (h.shape,)))
    if w.ndim != 1:
        raise ValueError(
            '%s: widths must be one-dimensional, got shape %s'
            % (caller, (w.shape,)))
    if h.shape[-3] != w.shape[0]:
        raise ValueError(
            '%s: got %d Hamiltonians per chord but %d widths; there must be '
            'one width per slab' % (caller, h.shape[-3], w.shape[0]))
    if w.shape[0] == 0:
        raise ValueError('%s: at least one slab is required' % caller)
    if np.any(w < 0.0):
        raise ValueError('%s: slab widths cannot be negative' % caller)

    return h, w


def _evolution_operator_slabs_batch(
    hamiltonian_matrices: Union[list, np.ndarray],
    widths: Union[list, np.ndarray],
    n_flavors: int,
    caller: str
) -> np.ndarray:
    r"""Returns one evolution operator per chord in a batch.

    The batched counterpart of `_evolution_operator_slabs`.  The chords
    are independent of one another, so the compiled path spreads them
    over the available cores --- an axis the per-chord kernel does not
    have, since the product *along* a chord cannot be reordered.

    Parameters
    ----------
    hamiltonian_matrices : array_like
        Hamiltonians, of shape ``(..., n_slabs, n_flavors, n_flavors)``.
    widths : array_like
        Slab widths, of shape ``(n_slabs,)``, shared by every chord.
    n_flavors : int
        Number of neutrino flavors, 2, 3, or 4.
    caller : str
        Name of the calling routine, used in error messages.

    Returns
    -------
    numpy.ndarray
        The evolution operators, of shape
        ``(..., n_flavors, n_flavors)``.
    """
    h, w = _check_slabs_batch(hamiltonian_matrices, widths, n_flavors, caller)

    n_chords = int(np.prod(h.shape[:-3])) if h.shape[:-3] else 1
    if fastkernels.worthwhile_slabs(n_flavors, n_chords*w.shape[0]):
        if n_flavors == 2:
            return fastkernels.slab_product_2nu_batch_kernel(h, w)
        if n_flavors == 3:
            return fastkernels.slab_product_3nu_batch_kernel(h, w)
        # As in the per-chord path, the expansion acts on the traceless
        # part and the dropped per-slab phase cancels in every probability
        return fastkernels.slab_product_4nu_batch_kernel(
            oscprob4nu._traceless_part(h), w, oscprob4nu.POLISH_ROOTS)

    # One batched call for every slab of every chord at once: the widths
    # broadcast along the chord axes, since the geometry is shared.
    w_b = np.broadcast_to(w, h.shape[:-2])
    if n_flavors == 2:
        u_slabs = np.asarray(oscprob2nu.evolution_operator_2nu(h, w_b))
    elif n_flavors == 3:
        u_slabs = np.asarray(oscprob3nu.evolution_operator_3nu(h, w_b))
    else:
        u_slabs = np.asarray(oscprob4nu.evolution_operator_4nu(h, w_b))

    # U = U_n ... U_1 along the slab axis, for every chord at once.  The
    # slab axis is -3, so each step is one batched matrix product.
    u_total = u_slabs[..., 0, :, :]
    for k in range(1, u_slabs.shape[-3]):
        u_total = u_slabs[..., k, :, :] @ u_total

    return u_total


def _probabilities_slabs_batch(
    hamiltonian_matrices: Union[list, np.ndarray],
    widths: Union[list, np.ndarray],
    n_flavors: int,
    caller: str
) -> np.ndarray:
    r"""Returns the probabilities for a batch of chords.

    Parameters
    ----------
    hamiltonian_matrices : array_like
        Hamiltonians, of shape ``(..., n_slabs, n_flavors, n_flavors)``.
    widths : array_like
        Slab widths, of shape ``(n_slabs,)``, shared by every chord.
    n_flavors : int
        Number of neutrino flavors, 2, 3, or 4.
    caller : str
        Name of the calling routine, used in error messages.

    Returns
    -------
    numpy.ndarray
        The probabilities, of shape ``(..., n_flavors*n_flavors)``, with
        the initial flavor varying slowest --- the same ordering the
        per-chord routines return as a tuple.
    """
    u = _evolution_operator_slabs_batch(hamiltonian_matrices, widths,
                                        n_flavors, caller)

    # P_ab = |U_ba|^2, and the returned ordering runs over the initial
    # flavor slowest, so the transpose comes before the flattening
    p = np.abs(np.swapaxes(u, -1, -2))**2.0

    return p.reshape(p.shape[:-2] + (n_flavors*n_flavors,))


N_SLABS_MAX = 1024
r"""int: Module-level constant.

Default ceiling on the refinement a tolerance request may ask for.

A tolerance is a statement about the answer, not about the cost, so
without a ceiling a tolerance the discretisation cannot reach --- one
below the round-off of the arithmetic itself, say --- would refine until
it ran out of memory.  Refinement past a thousand sub-slabs per segment
is in practice a sign that the tolerance was mis-stated rather than that
the answer is nearly in reach, so that is where the routines stop and
raise.  Pass ``n_max`` to move it.

.. versionadded:: 1.12.0
"""


def _check_tolerances(
    rtol: Optional[float],
    atol: Optional[float],
    caller: str
) -> Tuple[float, float]:
    r"""Validates a tolerance pair and fills in the unset one with zero.

    Parameters
    ----------
    rtol : float or None
        Relative tolerance, or None if not given.
    atol : float or None
        Absolute tolerance, or None if not given.
    caller : str
        Name of the calling routine, used in error messages.

    Returns
    -------
    tuple of float
        The two tolerances, with an unset one returned as zero.

    Raises
    ------
    ValueError
        If either is negative, or if both are zero or unset, which asks
        for an exact answer from an approximation.
    """
    r = 0.0 if rtol is None else float(rtol)
    a = 0.0 if atol is None else float(atol)

    if r < 0.0 or a < 0.0:
        raise ValueError(
            '%s: tolerances cannot be negative; got rtol=%s, atol=%s'
            % (caller, rtol, atol))
    if r == 0.0 and a == 0.0:
        raise ValueError(
            '%s: a tolerance of zero cannot be met by a discretisation; '
            'give a positive rtol or atol, or pass neither to use '
            'n_slabs_per_segment as given' % caller)

    return r, a


def _n_for_tolerance(
    evaluate: Callable[[int], np.ndarray],
    rtol: Optional[float],
    atol: Optional[float],
    n_start: int,
    n_max: int,
    caller: str
) -> Tuple[int, np.ndarray]:
    r"""Returns the smallest tried subdivision meeting a tolerance.

    The discretisation is second-order accurate --- midpoint sampling
    within each segment --- so halving the sub-slab width quarters the
    error.  That is what makes the error *measurable* without knowing
    the exact answer: with :math:`e(n) = 4 e(2n)`, two evaluations
    differ by

    .. math::

       P(2n) - P(n) = e(n) - e(2n) = 3 e(2n) ,

    so a third of the gap between consecutive refinements estimates the
    error of the finer one, and four thirds of it estimates the error of
    the coarser.  Both are used: the coarser test is what lets a loose
    tolerance be met by ``n_start`` itself rather than by twice it.

    The subdivision doubles until the estimate passes, rather than
    solving the second-order law for the required :math:`n` in one jump.
    Doubling costs little --- the evaluations form a geometric series,
    so reaching :math:`n` costs about twice what evaluating at
    :math:`n` costs on its own --- and it means every returned value has
    had its error *measured* rather than extrapolated.  That matters
    here because the caller asked to be told when the tolerance cannot
    be met, and an extrapolated error cannot tell anyone that.
    Extrapolating downwards would be worse still: the law is asymptotic,
    and at one sub-slab per segment the observed errors depart from it
    by an order of magnitude.

    Parameters
    ----------
    evaluate : callable
        Takes a subdivision count and returns the probabilities at it,
        as an array.  Every entry must be comparable across calls.
    rtol : float or None
        Relative tolerance, against the finer evaluation.
    atol : float or None
        Absolute tolerance.
    n_start : int
        Coarsest subdivision to try, and the smallest that can be
        returned.
    n_max : int
        Largest subdivision to try.
    caller : str
        Name of the calling routine, used in error messages.

    Returns
    -------
    tuple
        The subdivision met, and the probabilities evaluated at it.

    Raises
    ------
    ValueError
        If the tolerances are invalid, if ``n_start`` is not positive,
        or if the tolerance is not met by ``n_max``.
    """
    r, a = _check_tolerances(rtol, atol, caller)

    n_start = int(n_start)
    n_max = int(n_max)
    if n_start < 1:
        raise ValueError('%s: n_start must be at least 1, got %d'
                         % (caller, n_start))
    # Two evaluations are what an error estimate costs, so a budget that
    # cannot afford the pair cannot answer the question at all
    if n_max < 2*n_start:
        raise ValueError(
            '%s: n_max (%d) must be at least twice n_start (%d); the error '
            'is estimated by comparing consecutive refinements, so there is '
            'nothing to compare below that' % (caller, n_max, n_start))

    n = n_start
    p_coarse = np.asarray(evaluate(n), dtype=float)
    coarse_untested = True
    worst = np.inf

    while 2*n <= n_max:
        p_fine = np.asarray(evaluate(2*n), dtype=float)
        gap = np.abs(p_fine - p_coarse)

        # Four thirds of the gap is the error of the coarser evaluation,
        # so a tolerance loose enough to be met by n_start is met by it
        # rather than by twice it.  Only worth asking on the first pass:
        # every later coarse value is one this loop has already refused.
        if coarse_untested and np.all(4.0*gap/3.0 <= a + r*np.abs(p_coarse)):
            return n, p_coarse
        coarse_untested = False

        if np.all(gap/3.0 <= a + r*np.abs(p_fine)):
            return 2*n, p_fine

        # The error of the finest evaluation made, kept for the message
        # below: once the loop ends, the coarse and fine values are the
        # same array and the gap can no longer be recovered from them.
        worst = float(np.max(gap/3.0))
        n, p_coarse = 2*n, p_fine

    raise ValueError(
        '%s: could not meet rtol=%s, atol=%s with at most %d slabs per '
        'segment; the largest error estimate at %d was %.3e.  Raise n_max '
        'if the tolerance is genuinely wanted, or loosen the tolerance'
        % (caller, rtol, atol, n_max, n, worst))


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
        Stack of Hamiltonians, of shape ``(n, n_flavors, n_flavors)``,
        or of shape ``(..., n, n_flavors, n_flavors)`` for a batch of
        chords sharing one geometry.
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
        ``(n_flavors, n_flavors)``, or ``(..., n_flavors, n_flavors)``
        for a batch.
    """
    # A leading axis beyond the slab axis means a batch of chords that
    # share their widths; the batched path composes all of them at once.
    if np.ndim(hamiltonian_matrices) > 3:
        return _evolution_operator_slabs_batch(hamiltonian_matrices, widths,
                                               n_flavors, caller)

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

    .. versionchanged:: 1.13.1
       Accepts a batch of chords that share one set of slab widths,
       of shape ``(..., n, 2, 2)``, and returns one result per
       chord.  A single chord returns exactly what it returned
       before.  The batch is composed in one pass rather than one
       chord at a time, which is worth roughly an order of
       magnitude on an energy scan across a fixed profile; it agrees
       with the per-chord
       result to round-off rather than bit for bit, because the two
       take different paths through the compiled backend.

    Parameters
    ----------
    hamiltonian_matrices : array_like
        Hamiltonians, of shape ``(n, 2, 2)``, one per slab and ordered
        along the trajectory, in units of eV.  May instead have
        shape ``(..., n, 2, 2)``, a batch of chords that share the
        widths below and differ only in their Hamiltonians --- an
        energy scan across a fixed profile, typically, asked for as
        ``hamiltonian_2nu_matter(h_vac, energies[:, None], vcc)``.
    widths : array_like
        Slab widths, of shape ``(n,)``, in units of eV\ :sup:`-1`.  Use
        `globaldefs.CONV_KM_TO_INV_EV` to convert from km.

    Returns
    -------
    numpy.ndarray
        The evolution operator, of shape ``(2, 2)``, indexed
        ``(final, initial)``, or of shape ``(..., 2, 2)`` carrying
        one such operator per chord when a batch is given.

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

    .. versionchanged:: 1.13.1
       Accepts a batch of chords that share one set of slab widths,
       of shape ``(..., n, 3, 3)``, and returns one result per
       chord.  A single chord returns exactly what it returned
       before.  The batch is composed in one pass rather than one
       chord at a time, which is worth roughly an order of
       magnitude on an energy scan across a fixed profile; it agrees
       with the per-chord
       result to round-off rather than bit for bit, because the two
       take different paths through the compiled backend.

    Parameters
    ----------
    hamiltonian_matrices : array_like
        Hamiltonians, of shape ``(n, 3, 3)``, one per slab and ordered
        along the trajectory, in units of eV.  May instead have
        shape ``(..., n, 3, 3)``, a batch of chords that share the
        widths below and differ only in their Hamiltonians --- an
        energy scan across a fixed profile, typically, asked for as
        ``hamiltonian_3nu_matter(h_vac, energies[:, None], vcc)``.
    widths : array_like
        Slab widths, of shape ``(n,)``, in units of eV\ :sup:`-1`.  Use
        `globaldefs.CONV_KM_TO_INV_EV` to convert from km.

    Returns
    -------
    numpy.ndarray
        The evolution operator, of shape ``(3, 3)``, indexed
        ``(final, initial)``, or of shape ``(..., 3, 3)`` carrying
        one such operator per chord when a batch is given.

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

    .. versionchanged:: 1.13.1
       Accepts a batch of chords that share one set of slab widths,
       of shape ``(..., n, 2, 2)``, and returns one result per
       chord.  A single chord returns exactly what it returned
       before.  The batch is composed in one pass rather than one
       chord at a time, which is worth roughly an order of
       magnitude on an energy scan across a fixed profile; it agrees
       with the per-chord
       result to round-off rather than bit for bit, because the two
       take different paths through the compiled backend.

    Parameters
    ----------
    hamiltonian_matrices : array_like
        Hamiltonians, of shape ``(n, 2, 2)``, one per slab and ordered
        along the trajectory, in units of eV.  May instead have
        shape ``(..., n, 2, 2)``, a batch of chords that share the
        widths below and differ only in their Hamiltonians --- an
        energy scan across a fixed profile, typically, asked for as
        ``hamiltonian_2nu_matter(h_vac, energies[:, None], vcc)``.
    widths : array_like
        Slab widths, of shape ``(n,)``, in units of eV\ :sup:`-1`.

    Returns
    -------
    tuple of float or numpy.ndarray
        The probabilities
        :math:`P_{ee}, P_{e\mu}, P_{\mu e}, P_{\mu\mu}`.  Given a batch,
        an array of shape ``(..., 4)`` instead, with the same
        ordering along the last axis.

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
    if np.ndim(hamiltonian_matrices) > 3:
        return _probabilities_slabs_batch(
            hamiltonian_matrices, widths, 2,
            'probabilities_2nu_slabs')

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

    .. versionchanged:: 1.13.1
       Accepts a batch of chords that share one set of slab widths,
       of shape ``(..., n, 3, 3)``, and returns one result per
       chord.  A single chord returns exactly what it returned
       before.  The batch is composed in one pass rather than one
       chord at a time, which is worth roughly an order of
       magnitude on an energy scan across a fixed profile; it agrees
       with the per-chord
       result to round-off rather than bit for bit, because the two
       take different paths through the compiled backend.

    Parameters
    ----------
    hamiltonian_matrices : array_like
        Hamiltonians, of shape ``(n, 3, 3)``, one per slab and ordered
        along the trajectory, in units of eV.  May instead have
        shape ``(..., n, 3, 3)``, a batch of chords that share the
        widths below and differ only in their Hamiltonians --- an
        energy scan across a fixed profile, typically, asked for as
        ``hamiltonian_3nu_matter(h_vac, energies[:, None], vcc)``.
    widths : array_like
        Slab widths, of shape ``(n,)``, in units of eV\ :sup:`-1`.

    Returns
    -------
    tuple of float or numpy.ndarray
        The nine probabilities, with the initial flavor varying
        slowest.  Given a batch, an array of shape ``(..., 9)``
        instead, with the same ordering along the last axis.

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

    Several chords across the same two slabs, in one call.  The widths
    are given once, because the geometry is shared; only the
    Hamiltonians differ, which is what an energy scan is:

    .. jupyter-execute::

        import slabs

        import numpy as np
        H = np.array([np.diag([1.0, 0.0, -1.0]),
                      np.diag([0.5, 0.0, -0.5])], dtype=complex)
        stack = np.array([1.0, 2.0, 3.0])[:, None, None, None]*H
        prob = slabs.probabilities_3nu_slabs(stack, [0.2, 0.3])
        print(prob.shape)
        print('%.6f  %.6f' % (prob[0, 0], prob[2, 0]))
    """
    if np.ndim(hamiltonian_matrices) > 3:
        return _probabilities_slabs_batch(
            hamiltonian_matrices, widths, 3,
            'probabilities_3nu_slabs')

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

    .. versionchanged:: 1.13.1
       Accepts a batch of chords that share one set of slab widths,
       of shape ``(..., n, 4, 4)``, and returns one result per
       chord.  A single chord returns exactly what it returned
       before.  The batch is composed in one pass rather than one
       chord at a time, which is worth roughly a few per cent on an
       energy scan across a fixed profile; it agrees with the
       per-chord
       result to round-off rather than bit for bit, because the two
       take different paths through the compiled backend.

    Parameters
    ----------
    hamiltonian_matrices : array_like
        Hamiltonians, of shape ``(n, 4, 4)``, one per slab and ordered
        along the trajectory, in units of eV.  May instead have
        shape ``(..., n, 4, 4)``, a batch of chords that share the
        widths below and differ only in their Hamiltonians --- an
        energy scan across a fixed profile, typically, asked for as
        ``hamiltonian_4nu_matter(h_vac, energies[:, None], vcc)``.
    widths : array_like
        Slab widths, of shape ``(n,)``, in units of eV\ :sup:`-1`.  Use
        `globaldefs.CONV_KM_TO_INV_EV` to convert from km.

    Returns
    -------
    numpy.ndarray
        The evolution operator, of shape ``(4, 4)``, indexed
        ``(final, initial)``, or of shape ``(..., 4, 4)`` carrying
        one such operator per chord when a batch is given.

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

    .. versionchanged:: 1.13.1
       Accepts a batch of chords that share one set of slab widths,
       of shape ``(..., n, 4, 4)``, and returns one result per
       chord.  A single chord returns exactly what it returned
       before.  The batch is composed in one pass rather than one
       chord at a time, which is worth roughly a few per cent on an
       energy scan across a fixed profile; it agrees with the
       per-chord
       result to round-off rather than bit for bit, because the two
       take different paths through the compiled backend.

    Parameters
    ----------
    hamiltonian_matrices : array_like
        Hamiltonians, of shape ``(n, 4, 4)``, one per slab and ordered
        along the trajectory, in units of eV.  May instead have
        shape ``(..., n, 4, 4)``, a batch of chords that share the
        widths below and differ only in their Hamiltonians --- an
        energy scan across a fixed profile, typically, asked for as
        ``hamiltonian_4nu_matter(h_vac, energies[:, None], vcc)``.
    widths : array_like
        Slab widths, of shape ``(n,)``, in units of eV\ :sup:`-1`.

    Returns
    -------
    tuple of float or numpy.ndarray
        The sixteen probabilities, with the initial flavor varying
        slowest.  Given a batch, an array of shape ``(..., 16)``
        instead, with the same ordering along the last axis.

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
    if np.ndim(hamiltonian_matrices) > 3:
        return _probabilities_slabs_batch(
            hamiltonian_matrices, widths, 4,
            'probabilities_4nu_slabs')

    u = evolution_operator_4nu_slabs(hamiltonian_matrices, widths)

    # P_ab = |U_ba|^2: the evolution operator is indexed (final, initial)
    return tuple(abs(u[beta][alpha])**2.0
                 for alpha in range(4) for beta in range(4))


def _probabilities_profile(
    hamiltonian_of: Callable,
    baseline: Union[int, float],
    n_flavors: int,
    n_slabs: int,
    rtol: Optional[float],
    atol: Optional[float],
    n_max: int,
    return_n_slabs: bool,
    caller: str
) -> Union[Tuple[float, ...], tuple]:
    r"""Returns the probabilities across a continuously varying profile.

    The common body of the three public profile routines.

    Parameters
    ----------
    hamiltonian_of : callable
        Takes an array of positions and returns one Hamiltonian per
        position.
    baseline : int or float
        Total length of the trajectory, in units of eV\ :sup:`-1`.
    n_flavors : int
        Number of neutrino flavors, 2, 3, or 4.
    n_slabs : int
        Number of equal slabs, or the coarsest to try when refining.
    rtol : float or None
        Relative tolerance, or None.
    atol : float or None
        Absolute tolerance, or None.
    n_max : int
        Largest number of slabs to try when refining.
    return_n_slabs : bool
        Whether to return the number of slabs used.
    caller : str
        Name of the calling routine, used in error messages.

    Returns
    -------
    tuple
        The probabilities, paired with the slab count when
        `return_n_slabs` is set.
    """
    if not callable(hamiltonian_of):
        raise ValueError('%s: hamiltonian_of must be callable, got %s'
                         % (caller, type(hamiltonian_of).__name__))
    baseline = float(baseline)
    if not baseline > 0.0:
        raise ValueError('%s: baseline must be positive, got %s'
                         % (caller, baseline))
    if int(n_slabs) < 1:
        raise ValueError('%s: n_slabs must be at least 1, got %s'
                         % (caller, n_slabs))

    routine = {2: probabilities_2nu_slabs,
               3: probabilities_3nu_slabs,
               4: probabilities_4nu_slabs}[n_flavors]

    def evaluate(n: int) -> np.ndarray:
        # Equal slabs sampled at their midpoints, which is second-order
        # accurate and so refines by the law `_n_for_tolerance` assumes.
        # Sampling at an end would be first-order and would make the
        # error estimate there wrong rather than merely pessimistic.
        edges = np.linspace(0.0, baseline, n+1)
        midpoints = (edges[:-1] + edges[1:])/2.0
        h = np.asarray(hamiltonian_of(midpoints), dtype=complex)
        if h.shape != (n, n_flavors, n_flavors):
            raise ValueError(
                '%s: hamiltonian_of returned shape %s for %d positions; it '
                'must return one %dx%d Hamiltonian per position, of shape '
                '(%d, %d, %d)' % (caller, (h.shape,), n, n_flavors,
                                  n_flavors, n, n_flavors, n_flavors))
        return np.asarray(routine(h, np.diff(edges)), dtype=float)

    if rtol is None and atol is None:
        p = evaluate(int(n_slabs))
        n_used = int(n_slabs)
    else:
        n_used, p = _n_for_tolerance(evaluate, rtol, atol, int(n_slabs),
                                     n_max, caller)

    probabilities = tuple(float(x) for x in p)

    return (probabilities, n_used) if return_n_slabs else probabilities


def probabilities_2nu_profile(
    hamiltonian_of: Callable,
    baseline: Union[int, float],
    n_slabs: int = 8,
    rtol: Optional[float] = None,
    atol: Optional[float] = None,
    n_max: int = N_SLABS_MAX,
    return_n_slabs: bool = False
) -> Union[Tuple[float, ...], tuple]:
    r"""Returns the two-flavor probabilities across a varying profile.

    .. versionadded:: 1.12.0

    See `probabilities_3nu_profile`, of which this is the two-flavor
    counterpart in every respect.

    Parameters
    ----------
    hamiltonian_of : callable
        Takes an array of positions along the trajectory, in units of
        eV\ :sup:`-1`, and returns the Hamiltonian at each, as an array
        of shape ``(len(positions), 2, 2)`` in units of eV.
    baseline : int or float
        Total length of the trajectory, in units of eV\ :sup:`-1`.
    n_slabs : int, optional
        Number of equal slabs.  Default: 8.
    rtol : float, optional
        Relative tolerance on every returned probability.  Default:
        None.
    atol : float, optional
        Absolute tolerance on every returned probability.  Default:
        None.
    n_max : int, optional
        Largest number of slabs the refinement may try.  Default:
        `N_SLABS_MAX`.
    return_n_slabs : bool, optional
        Whether to return the number of slabs used.  Default: False.

    Returns
    -------
    tuple of float
        The probabilities
        :math:`P_{ee}, P_{e\mu}, P_{\mu e}, P_{\mu\mu}`, paired with the
        number of slabs used when `return_n_slabs` is set.

    Raises
    ------
    ValueError
        As `probabilities_3nu_profile`.
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np
        import slabs

        baseline = 1.0e13
        H0 = np.diag([1.0e-13, -1.0e-13])
        H0[0, 1] = H0[1, 0] = 0.3e-13    # something for the profile to act on

        def H_of(x):
            # A potential that rises linearly along the trajectory.
            # Normalise by `baseline`, not by `x[-1]`: the midpoints move
            # as the refinement doubles, so dividing by the last one makes
            # the profile itself depend on `n_slabs`, which costs an order
            # of convergence and can leave the tolerance unreachable.
            h = np.broadcast_to(H0, (len(x), 2, 2)).copy()
            h[:, 0, 0] += 1.0e-13*x/baseline
            return h

        prob, n = slabs.probabilities_2nu_profile(
            H_of, baseline, atol=1.0e-8, return_n_slabs=True)
        print(n, '%.6f' % prob[0])
    """
    return _probabilities_profile(hamiltonian_of, baseline, 2, n_slabs,
                                  rtol, atol, n_max, return_n_slabs,
                                  'probabilities_2nu_profile')


def probabilities_3nu_profile(
    hamiltonian_of: Callable,
    baseline: Union[int, float],
    n_slabs: int = 8,
    rtol: Optional[float] = None,
    atol: Optional[float] = None,
    n_max: int = N_SLABS_MAX,
    return_n_slabs: bool = False
) -> Union[Tuple[float, ...], tuple]:
    r"""Returns the three-flavor probabilities across a varying profile.

    The general counterpart of `earth.probabilities_3nu_earth`: where
    that one knows about PREM, this one takes the profile as a callable
    and so serves any continuously varying Hamiltonian --- a hand-built
    density profile, a castle wall, a solar model, a matter potential
    that varies for a reason other than density.

    The trajectory is cut into equal slabs, the Hamiltonian is sampled
    at the midpoint of each, and the slabs are solved exactly and
    composed.  Midpoint sampling makes that second-order accurate, so
    the answer converges as the slabs are refined, and `rtol` and `atol`
    let the routine do the refining: it doubles the slab count until the
    measured error meets the tolerance.

    Where a profile has *discontinuities* --- a wall, a shell boundary
    --- equal slabs are the wrong tool, because no amount of refinement
    recovers a jump that straddles a slab.  Split the trajectory at the
    discontinuities yourself and call this once per piece, or hand the
    pieces to `probabilities_3nu_slabs` directly; that is exactly what
    :mod:`earth` does with the PREM shells.

    .. versionadded:: 1.12.0

    Parameters
    ----------
    hamiltonian_of : callable
        Takes an array of positions along the trajectory, in units of
        eV\ :sup:`-1`, measured from the start, and returns the
        Hamiltonian at each, as an array of shape
        ``(len(positions), 3, 3)`` in units of eV.  It is called once
        per refinement, with all the midpoints at once, so it should be
        vectorised rather than called in a loop.
    baseline : int or float
        Total length of the trajectory, in units of eV\ :sup:`-1`.  Use
        `globaldefs.CONV_KM_TO_INV_EV` to convert from km.
    n_slabs : int, optional
        Number of equal slabs, or the coarsest to try when a tolerance
        is given.  Default: 8.
    rtol : float, optional
        Relative tolerance on every returned probability.  Default:
        None, meaning `n_slabs` is used as given.
    atol : float, optional
        Absolute tolerance on every returned probability.  Default:
        None.  When both are given the threshold is
        ``atol + rtol*abs(P)``, the convention of `numpy.isclose`.
    n_max : int, optional
        Largest number of slabs the refinement may try before giving up
        and raising.  Default: `N_SLABS_MAX`.
    return_n_slabs : bool, optional
        Whether to return the number of slabs used alongside the
        probabilities.  Default: False.

    Returns
    -------
    tuple of float
        The nine probabilities, with the initial flavor varying slowest,
        paired with the number of slabs used when `return_n_slabs` is
        set.

    Raises
    ------
    ValueError
        If `hamiltonian_of` is not callable or does not return one
        Hamiltonian of the right size per position, if `baseline` is not
        positive, if `n_slabs` is not positive, if the tolerances are
        invalid, or if the tolerance is not met by `n_max`.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np
        import slabs

        baseline = 1.0e13
        H0 = np.diag([1.0e-13, 0.0, -1.0e-13])
        H0[0, 1] = H0[1, 0] = 0.3e-13    # something for the profile to act on

        def H_of(x):
            # A potential that rises linearly along the trajectory.
            # Normalise by `baseline`, not by `x[-1]`: the midpoints move
            # as the refinement doubles, so dividing by the last one makes
            # the profile itself depend on `n_slabs`, which costs an order
            # of convergence and can leave the tolerance unreachable.
            h = np.broadcast_to(H0, (len(x), 3, 3)).copy()
            h[:, 0, 0] += 1.0e-13*x/baseline
            return h

        prob, n = slabs.probabilities_3nu_profile(
            H_of, baseline, atol=1.0e-8, return_n_slabs=True)
        print(n, '%.6f' % prob[0])
    """
    return _probabilities_profile(hamiltonian_of, baseline, 3, n_slabs,
                                  rtol, atol, n_max, return_n_slabs,
                                  'probabilities_3nu_profile')


def probabilities_4nu_profile(
    hamiltonian_of: Callable,
    baseline: Union[int, float],
    n_slabs: int = 8,
    rtol: Optional[float] = None,
    atol: Optional[float] = None,
    n_max: int = N_SLABS_MAX,
    return_n_slabs: bool = False
) -> Union[Tuple[float, ...], tuple]:
    r"""Returns the four-flavor probabilities across a varying profile.

    .. versionadded:: 1.12.0

    See `probabilities_3nu_profile`, of which this is the four-flavor
    counterpart in every respect.

    Parameters
    ----------
    hamiltonian_of : callable
        Takes an array of positions along the trajectory, in units of
        eV\ :sup:`-1`, and returns the Hamiltonian at each, as an array
        of shape ``(len(positions), 4, 4)`` in units of eV.
    baseline : int or float
        Total length of the trajectory, in units of eV\ :sup:`-1`.
    n_slabs : int, optional
        Number of equal slabs.  Default: 8.
    rtol : float, optional
        Relative tolerance on every returned probability.  Default:
        None.
    atol : float, optional
        Absolute tolerance on every returned probability.  Default:
        None.
    n_max : int, optional
        Largest number of slabs the refinement may try.  Default:
        `N_SLABS_MAX`.
    return_n_slabs : bool, optional
        Whether to return the number of slabs used.  Default: False.

    Returns
    -------
    tuple of float
        The sixteen probabilities, with the initial flavor varying
        slowest, paired with the number of slabs used when
        `return_n_slabs` is set.

    Raises
    ------
    ValueError
        As `probabilities_3nu_profile`.
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np
        import slabs

        baseline = 1.0e13
        H0 = np.diag([2.0e-13, 1.0e-13, 0.0, -1.0e-13])
        H0[0, 1] = H0[1, 0] = 0.5e-13    # something for the profile to act on

        def H_of(x):
            # A potential that rises linearly along the trajectory.
            # Normalise by `baseline`, not by `x[-1]`: the midpoints move
            # as the refinement doubles, so dividing by the last one makes
            # the profile itself depend on `n_slabs`, which costs an order
            # of convergence and can leave the tolerance unreachable.
            h = np.broadcast_to(H0, (len(x), 4, 4)).copy()
            h[:, 0, 0] += 1.0e-13*x/baseline
            return h

        prob, n = slabs.probabilities_4nu_profile(
            H_of, baseline, atol=1.0e-8, return_n_slabs=True)
        print(n, '%.6f' % prob[0])
    """
    return _probabilities_profile(hamiltonian_of, baseline, 4, n_slabs,
                                  rtol, atol, n_max, return_n_slabs,
                                  'probabilities_4nu_profile')
