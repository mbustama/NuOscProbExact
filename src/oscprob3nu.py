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
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

__all__ = ['CHECK_HERMITICITY', 'SMALL_BATCH',
           'hamiltonian_3nu_coefficients', 'tensor_d', 'star',
           'su3_invariants', 'psi_roots',
           'evolution_operator_3nu_u_coefficients', 'evolution_operator_3nu',
           'probabilities_3nu']

from typing import List, Optional, Tuple, Union

import cmath
import math

import numpy as np

try:
    import fastkernels
except ImportError:                                       # pragma: no cover
    # Copying this file on its own into another project is a supported way
    # to use **NuOscProbExact**, and is documented as such --- but it stopped
    # working in 1.6.0, when the optional compiled backend was added and
    # imported unconditionally.  A lone copy raised ImportError on the first
    # line that mattered, which is a poor return for a promise the
    # documentation makes twice.
    #
    # The backend is optional by design, so its absence is answered the same
    # way its being switched off is: `worthwhile` says no, and the NumPy path
    # runs.  Nothing else in this module touches it.
    class _NoFastKernels:
        r"""Stands in for :mod:`fastkernels` when it is not importable."""

        HAVE_NUMBA = False
        USE_NUMBA = False

        @staticmethod
        def available():
            r"""Returns False: there is no compiled backend here."""
            return False

        @staticmethod
        def worthwhile(n_flavors, size):
            r"""Returns False: no stack is worth a backend that is absent."""
            return False

    fastkernels = _NoFastKernels()


SQRT3 = np.sqrt(3.0)
r"""float: Module-level constant equal to :math:`\sqrt{3}`."""

SQRT3_INV = 1./np.sqrt(3.0)
r"""float: Module-level constant equal to :math:`1/\sqrt{3}`."""

NEG_HALF_SQRT3_INV = -SQRT3_INV/2.0
r"""float: Module-level constant equal to :math:`-1/(2\sqrt{3})`."""

SMALL_BATCH = 10
r"""int: Module-level constant.

Stacks with at most this many elements are evaluated one at a time
through the scalar path.  A batched call carries a fixed cost of a
couple of hundred microseconds --- allocating and reducing a dozen small
arrays --- which for a short stack exceeds what the scalar path spends
on the whole job.  Measured crossover: eleven elements, so a single
Hamiltonian and baseline costs about a third of what the array path
would.  The two-flavor expansion does less work per element, so its
threshold is lower; see :data:`oscprob2nu.SMALL_BATCH`.
"""

DEGENERACY_TOL = 1.e-12
r"""float: Module-level constant.

Relative tolerance below which two latent roots :math:`\psi` are treated
as degenerate.  The general expression for the :math:`u_k` divides by
:math:`3\psi_m^2 - |h|^2`, which vanishes at a repeated root, so a
degenerate spectrum is handled by a separate, exact expression.
"""


CHECK_HERMITICITY = True
r"""bool: Module-level switch.

Whether to verify that the Hamiltonian handed in is Hermitian before
evaluating anything.  It is on by default, and the reason is that the
failure it catches is silent: a non-Hermitian matrix does not raise, and
does not produce obviously broken output either --- the probabilities it
returns still sum to one, so the usual sanity check a caller would apply
cannot tell that anything went wrong.  The expansion assumes
Hermiticity; without it the result is meaningless rather than merely
inaccurate.

Checking is not free, and the cost is stated here rather than buried,
because it is larger than one might expect.  Validating a stack is a
pass over it, which is the same order of work as evaluating it --- and
the compiled kernel has made evaluating it *fast*, so on a large stack
the check dominates.  Measured by interleaving the two settings and
taking the best of fifteen rounds each:

===============  ==========  ==========
Stack            2 flavors   4 flavors
===============  ==========  ==========
2 000 points     1.5x        1.3x
200 000 points   5.7x        3.2x
===============  ==========  ==========

Three flavors sits between them, at 1.8x and 3.9x.  So a scan that the
compiled backend was installed to speed up can spend most of its time
here instead.  For production scans whose Hamiltonians come from a
construction already known to be Hermitian --- everything
:mod:`hamiltonians3nu` builds is Hermitian to round-off, as the table
below records --- set this to ``False``.

The default is nevertheless ``True``, because the alternative is a
library that silently returns meaningless numbers to anyone who makes
this mistake once, and finding out costs far more than the check does.

The tolerance is relative to the largest entry of the Hamiltonian, so a
matrix assembled in floating point passes: the ones this library builds
are Hermitian to about :math:`2 \times 10^{-17}` relative, against a
tolerance of :math:`10^{-12}`.

.. versionadded:: 1.11.0
"""

_HERMITICITY_TOL = 1.e-12
r"""float: Module-level constant.

Relative tolerance for `CHECK_HERMITICITY`, measured against the largest
entry of the Hamiltonian.
"""


def _check_hermitian(h_matrix: np.ndarray, caller: str) -> None:
    r"""Raises unless `h_matrix` is Hermitian to `_HERMITICITY_TOL`.

    Compares the 3 independent pairs and the imaginary parts of the
    diagonal, rather than forming ``H - H^dagger``, which would allocate
    a temporary the size of the whole stack.  Real and imaginary parts
    are compared separately, on views rather than copies: the condition
    is :math:`\mathrm{Re}\,H_{ij} = \mathrm{Re}\,H_{ji}` together
    with :math:`\mathrm{Im}\,H_{ij} = -\mathrm{Im}\,H_{ji}`, which
    needs no complex arithmetic and, unlike :func:`numpy.abs` on a
    complex array, no square root per element.  That is worth about a
    factor of three on a large stack, where this check would otherwise
    cost several times the evaluation it guards.

    Parameters
    ----------
    h_matrix : numpy.ndarray
        Hamiltonian, or stack of them, of shape ``(..., 3, 3)``.
    caller : str
        Name of the calling routine, used in the error message.

    Returns
    -------
    None
        Nothing; the routine either returns or raises.

    Raises
    ------
    ValueError
        If any element of the stack is not Hermitian.
    """
    if h_matrix.size == 0:
        return

    real, imaginary = h_matrix.real, h_matrix.imag

    # `np.abs(...).max()` allocates a float array the size of the stack,
    # and is still the quickest way to the largest entry: replacing it with
    # four reductions over `real` and `imaginary`, which allocate nothing,
    # was measured 1.4x *slower* on two hundred thousand elements, because
    # those views are strided over the complex array while `np.abs` reads
    # it contiguously.  Tried, measured, reverted.
    scale = max(float(np.max(np.abs(real))), float(np.max(np.abs(imaginary))))

    # A non-finite entry has to be caught here rather than left to
    # propagate.  It would otherwise make `scale` infinite or nan, hence
    # `tolerance` infinite or nan, and every comparison below false ---
    # so a Hamiltonian that is both non-finite *and* non-Hermitian would
    # pass a check whose whole purpose is to refuse the second.
    if not np.isfinite(scale):
        raise ValueError(
            '%s: the Hamiltonian has a non-finite entry, so it is neither '
            'Hermitian nor usable.  Set %s.CHECK_HERMITICITY = False to '
            'skip this check.' % (caller, 'oscprob3nu'))

    tolerance = _HERMITICITY_TOL*scale if scale > 0.0 else _HERMITICITY_TOL

    complaint = (
        '%s: the Hamiltonian is not Hermitian%s.  The expansion assumes '
        'Hermiticity, and without it the probabilities returned are '
        'meaningless even though they still sum to one.  Set '
        'oscprob3nu.CHECK_HERMITICITY = False to skip this check.')

    for i in range(3):
        if np.any(np.abs(imaginary[..., i, i]) > tolerance):
            raise ValueError(complaint % (
                caller, ' --- the diagonal entry (%d, %d) has a non-zero '
                'imaginary part' % (i, i)))
        for j in range(i+1, 3):
            if (np.any(np.abs(real[..., i, j] - real[..., j, i]) > tolerance)
                    or np.any(np.abs(imaginary[..., i, j]
                                     + imaginary[..., j, i]) > tolerance)):
                raise ValueError(complaint % (
                    caller, ' --- entry (%d, %d) is not the complex '
                    'conjugate of entry (%d, %d)' % (i, j, j, i)))



def hamiltonian_3nu_coefficients(
    hamiltonian_matrix: Union[list, np.ndarray]
) -> List[float]:
    r"""Returns the :math:`h_k` of the SU(3) expansion of the Hamiltonian.

    Computes the coefficients :math:`h_1, \ldots, h_8` of the SU(3)
    expansion :math:`H = h_0 \mathbb{1} + h_k \lambda^k` of the
    three-flavor Hamiltonian `hamiltonian_matrix`, which is assumed to
    be given in the flavor basis.  The coefficient :math:`h_0`
    contributes only an overall phase to the evolution operator and is
    not returned.

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.1.0
       Returns real floats.  The coefficients of a Hermitian Hamiltonian
       are real by construction, but the routine previously returned a
       mixture of floats and complex numbers.

    .. versionchanged:: 1.4.0
       Faster, with identical results --- all 42 figures generated by
       ``run_testsuite.py`` are byte-for-byte those of 1.3.0. The scalar
       path stopped dispatching NumPy for single numbers:
       :func:`numpy.real`, :func:`numpy.imag`, :obj:`numpy.arccos`,
       :func:`numpy.clip` and :obj:`numpy.sqrt` on one number give way
       to attribute access and the :mod:`math` module.

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
    .. jupyter-execute::

        import oscprob3nu

        hamiltonian_matrix = [[1.0+0.0j, 0.0+2.0j, 0.0-1.0j],
                              [0.0-2.0j, 3.0+0.0j, 3.0+0.0j],
                              [0.0+1.0j, 3.0-0.0j, -5.0+0.0j]]
        h_coeffs = oscprob3nu.hamiltonian_3nu_coefficients(hamiltonian_matrix)
        print('  '.join(['%.6f' % (h+0.0) for h in h_coeffs]))
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
    h1 = H12.real
    h2 = -H12.imag
    h3 = (H11-H22).real/2.0
    h4 = H13.real
    h5 = -H13.imag
    h6 = H23.real
    h7 = -H23.imag
    h8 = (H11+H22-2.0*H33).real*SQRT3/6.0

    return [float(h1), float(h2), float(h3), float(h4), float(h5),
            float(h6), float(h7), float(h8)]


def tensor_d(i: int, j: int, k: int) -> float:
    r"""Returns the tensor :math:`d_{ijk}` of the SU(3) algebra.

    Returns the totally symmetric SU(3) tensor
    :math:`d_{ijk} = \frac{1}{4}\mathrm{Tr}
    (\{\lambda_i, \lambda_j\} \lambda_k)`, defined in [1]_.

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.1.0
       Validates its indices and raises :exc:`IndexError` on one outside
       0-7. The dispatch previously fell off the end of an ``elif``
       chain and returned ``None``, so the failure surfaced far from its
       cause.

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
    .. jupyter-execute::

        import oscprob3nu

        print('%.6f' % oscprob3nu.tensor_d(0, 0, 7))
        print('%.6f' % oscprob3nu.tensor_d(0, 1, 2))
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


def star(i: int, h_coeffs: Union[list, np.ndarray]) -> float:
    r"""Returns the SU(3) star product :math:`(h \star h)_i`.

    Returns the SU(3) star product
    :math:`(h \star h)_i = d_{ijk} h^j h^k`, summed over repeated
    indices.

    .. versionadded:: 1.0.0

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
    .. jupyter-execute::

        import oscprob3nu

        print('%.6f' % oscprob3nu.star(0, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
        print('%.6f' % oscprob3nu.star(7, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
    """
    h = np.asarray(h_coeffs, dtype=float)

    return float(h @ _TENSOR_D[i] @ h)


def su3_invariants(h_coeffs: Union[list, np.ndarray]) -> Tuple[float, float]:
    r"""Returns the two SU(3) invariants, :math:`|h|^2` and
    :math:`\langle h \rangle`.

    Returns the two invariants of the SU(3) expansion,
    :math:`|h|^2 = h_i h_i` and
    :math:`\langle h \rangle = d_{ijk} h_i h_j h_k`.  They equal,
    respectively, :math:`\mathrm{Tr}(H_0^2)/2` and
    :math:`\mathrm{Tr}(H_0^3)/2`, with :math:`H_0` the traceless part of
    the Hamiltonian.

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.4.0
       Faster, with identical results --- all 42 figures generated by
       ``run_testsuite.py`` are byte-for-byte those of 1.3.0. The star
       product is computed once and passed on, rather than once to form
       :math:`\langle h \rangle` and again inside the expansion, and the
       sparse expansion of :math:`d_{ijk} h_j h_k` replaces the dense
       table for the eight scalar components.

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
    .. jupyter-execute::

        import oscprob3nu

        h2, h3 = oscprob3nu.su3_invariants([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        print('%.6f  %.6f' % (h2, h3))
    """
    star = _star_all(h_coeffs)

    # h2 = |h|^2
    h2 = sum([x*x for x in h_coeffs])

    # h3 = <h> = h_i (h*h)_i, using the star product computed once
    h3 = sum([x*y for x, y in zip(h_coeffs, star)])

    return float(h2), float(h3)


def psi_roots(h2: Union[int, float], h3: Union[int, float]) -> List[float]:
    r"""Returns the three latent roots :math:`\psi`.

    Returns the three latent roots :math:`\psi` of the characteristic
    equation :math:`\psi^3 - |h|^2 \psi - \frac{2}{3}\langle h \rangle
    = 0`, which are the eigenvalues of minus the traceless part of the
    Hamiltonian.  The roots are independent of the baseline.

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.1.0
       No longer returns NaN when the traceless part of the Hamiltonian
       vanishes, and the arc-cosine argument is clipped to :math:`[-1,
       1]` so that round-off cannot make the roots complex and the
       evolution operator non-unitary.

    .. versionchanged:: 1.4.0
       Faster, with identical results --- all 42 figures generated by
       ``run_testsuite.py`` are byte-for-byte those of 1.3.0. The scalar
       path stopped dispatching NumPy for single numbers:
       :func:`numpy.real`, :func:`numpy.imag`, :obj:`numpy.arccos`,
       :func:`numpy.clip` and :obj:`numpy.sqrt` on one number give way
       to attribute access and the :mod:`math` module.

    .. versionchanged:: 1.5.0
       Faster, with identical results; the probabilities agree with
       1.4.0 to 1.6e-13 across every code path.  The prefactor takes
       :math:`\sqrt{|h|^2}` once rather than three times.  The
       arc-cosine argument still forms :math:`|h|^2` to the power of
       -1.5, which is what the code does and what an earlier version of
       this note wrongly described as a division.

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
    .. jupyter-execute::

        import oscprob3nu

        psi = oscprob3nu.psi_roots(1.0, 0.0)
        print('  '.join(['%.6f' % (round(p, 9)+0.0) for p in sorted(psi)]))
    """
    if h2 <= 0.0:
        return [0.0, 0.0, 0.0]

    pre = 2.0*math.sqrt(h2)*SQRT3_INV

    # For a Hermitian Hamiltonian the argument lies in [-1, 1]; clipping
    # keeps round-off from producing spurious complex roots.  The bounds
    # are written out rather than passed through numpy, which would
    # dispatch a whole array machinery for one number.
    arg = -SQRT3*h3*pow(h2, -1.5)
    if arg <= -1.0:
        chi = math.pi
    elif arg >= 1.0:
        chi = 0.0
    else:
        chi = math.acos(arg)

    return [pre*math.cos((chi+2.0*math.pi*m)/3.0) for m in (1, 2, 3)]


TWO_SQRT3_INV = 2.0*SQRT3_INV
r"""float: Module-level constant equal to :math:`2/\sqrt{3}`."""


def _star_all(h: Union[list, np.ndarray]) -> Tuple[float, ...]:
    r"""Returns all eight components of the star product at once.

    The sparse expansion of :math:`(h \star h)_i = d_{ijk} h_j h_k`,
    written out.  The :math:`d` tensor has at most eight non-zero
    entries per component, so contracting the dense
    :math:`8\times8\times8` table costs several microseconds of NumPy
    dispatch to do a few dozen multiplications --- about six times
    longer than the arithmetic itself.  On the batched path the array
    machinery pays for itself and the table is used instead.

    These expressions are generated from, and checked against,
    `tensor_d` by ``tests/test_su3_algebra.py``; edit them only in step
    with that table.

    Parameters
    ----------
    h : array_like
        The eight coefficients :math:`h_k`.

    Returns
    -------
    tuple of float
        The eight components :math:`(h \star h)_i`.
    """
    h0, h1, h2, h3, h4, h5, h6, h7 = h

    return (
        TWO_SQRT3_INV*h0*h7 + h3*h5 + h4*h6,
        TWO_SQRT3_INV*h1*h7 - h3*h6 + h4*h5,
        TWO_SQRT3_INV*h2*h7 + 0.5*(h3*h3 + h4*h4 - h5*h5 - h6*h6),
        h0*h5 - h1*h6 + h2*h3 - SQRT3_INV*h3*h7,
        h0*h6 + h1*h5 + h2*h4 - SQRT3_INV*h4*h7,
        h0*h3 + h1*h4 - h2*h5 - SQRT3_INV*h5*h7,
        h0*h4 - h1*h3 - h2*h6 - SQRT3_INV*h6*h7,
        SQRT3_INV*(h0*h0 + h1*h1 + h2*h2 - h7*h7)
        - SQRT3_INV/2.0*(h3*h3 + h4*h4 + h5*h5 + h6*h6),
    )


def _abs2(z: complex) -> float:
    r"""Returns :math:`|z|^2` without taking a square root.

    :func:`abs` on a complex number computes a hypotenuse, which is then
    squared again; forming the sum of squares directly skips the square
    root and its rounding.

    Parameters
    ----------
    z : complex
        The number whose squared modulus is wanted.

    Returns
    -------
    float
        The squared modulus :math:`|z|^2`.
    """
    return z.real*z.real + z.imag*z.imag


def _u_coefficients_3nu_single(
    h: Union[list, np.ndarray],
    h2: Union[int, float],
    h3: Union[int, float],
    L: Union[int, float],
    star_coeffs: Optional[Union[tuple, np.ndarray]] = None
) -> List[complex]:
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
    star_coeffs : numpy.ndarray, optional
        The star product :math:`(h \star h)_k`, if the caller has
        already computed it --- which it has whenever it obtained
        :math:`\langle h \rangle` from it.  Recomputed here if omitted.

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

    if star_coeffs is None:
        star_coeffs = _star_all(h)

    # Find the closest pair of latent roots, to detect degeneracy
    psi0, psi1, psi2 = psi
    gaps = (abs(psi0-psi1), abs(psi0-psi2), abs(psi1-psi2))
    smallest = min(gaps)
    a, b, c = ((0, 1, 2), (0, 2, 1), (1, 2, 0))[gaps.index(smallest)]

    if smallest <= DEGENERACY_TOL*math.sqrt(h2):
        # Doubly degenerate root: the general expression would divide by
        # zero, so use the two-projector form instead.
        psi_deg = (psi[a]+psi[b])/2.0
        exp_deg = cmath.rect(1.0, L*psi_deg)
        exp_odd = cmath.rect(1.0, L*psi[c])
        weight = (exp_odd-exp_deg)/(psi_deg-psi[c])
        u0 = exp_deg + weight*psi_deg
        uk = [-1.j*weight*h[k] for k in range(0, 8)]

        return [u0]+uk

    # e^{i*L*psi_m}, and the same divided by the denominators of the
    # Lagrange interpolation, which do not depend on k
    # cmath.rect(1, t) is e^{it} without forming the complex argument
    exp0 = cmath.rect(1.0, L*psi0)
    exp1 = cmath.rect(1.0, L*psi1)
    exp2 = cmath.rect(1.0, L*psi2)
    w0 = exp0/(3.*psi0*psi0-h2)
    w1 = exp1/(3.*psi1*psi1-h2)
    w2 = exp2/(3.*psi2*psi2-h2)

    # Only two combinations of the three terms survive the sum over m,
    #   u_k = i [ (sum_m w_m psi_m) h_k - (sum_m w_m) (h*h)_k ],
    # so they are formed once instead of inside the loop over k
    weighted = w0*psi0 + w1*psi1 + w2*psi2
    total = w0 + w1 + w2

    u0 = (exp0+exp1+exp2)/3.
    uk = [1.j*(weighted*h[k] - total*star_coeffs[k]) for k in range(0, 8)]

    # [u0, u1, u2, u3, u4, u5, u6, u7, u8]
    return [u0]+uk


def _hamiltonian_3nu_coefficients_batch(h_matrix: np.ndarray) -> np.ndarray:
    r"""Returns the :math:`h_k` for a stack of Hamiltonians.

    The vectorised counterpart of `hamiltonian_3nu_coefficients`.

    Parameters
    ----------
    h_matrix : numpy.ndarray
        Hamiltonians, of shape ``(..., 3, 3)``.

    Returns
    -------
    numpy.ndarray
        The coefficients, of shape ``(8, ...)``.

    Notes
    -----
    The component index is the *first* axis, not the last, so that each
    :math:`h_k` is a contiguous array.  Everything downstream works one
    component at a time --- the star product, the nine entries of
    :math:`U_3` --- and reading those out of a strided view of a
    ``(..., 8)`` array costs about a third more.  It also lets
    `_star_all`, which unpacks its argument, serve both paths unchanged.
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
    ], axis=0)


def _u_coefficients_3nu_batch(h: np.ndarray, L: np.ndarray) -> np.ndarray:
    r"""Returns the nine :math:`u_k` for a stack of Hamiltonians.

    Every step of the SU(3) expansion --- the star product, the two
    invariants, the latent roots, and the Lagrange interpolation --- is
    evaluated for the whole stack at once.

    Parameters
    ----------
    h : numpy.ndarray
        The coefficients :math:`h_k`, of shape ``(8, ...)``, real.
    L : numpy.ndarray
        Baselines, of shape ``(...)``, broadcastable against `h`.

    Returns
    -------
    u0 : numpy.ndarray
        The coefficient :math:`u_0`, of shape ``(...)``.
    uk : numpy.ndarray
        The coefficients :math:`u_1, \ldots, u_8`, of shape
        ``(8, ...)``.  They come back separately rather than
        concatenated, because every caller takes them apart again.

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
    # The component index leads, so the batch axes of `h` are the
    # trailing ones and NumPy right-aligns them against `L`.  Padding
    # `h` with as many length-one axes as `L` has extra keeps that
    # alignment, and costs nothing: reshape returns a view.
    full = np.broadcast_shapes(h.shape[1:], np.shape(L))
    extra = len(full) - (h.ndim - 1)
    if extra > 0:
        h = h.reshape(h.shape[:1] + (1,)*extra + h.shape[1:])

    # Everything up to the exponential depends on the Hamiltonian alone,
    # so it is evaluated at the shape of `h` and only then broadcast
    # against the baselines.  Scanning one Hamiltonian over N baselines
    # therefore solves the characteristic equation once, not N times.
    # The star product uses the same sparse expansion as the scalar path.  Contracting the
    # dense table with einsum instead costs an order of magnitude more:
    # without a path plan it walks the 8x8x8 tensor for every element,
    # where this is a few dozen array multiplications.
    star = np.stack(_star_all(h))
    h2 = (h*h).sum(0)
    h3 = (h*star).sum(0)

    positive = h2 > 0.0
    safe_h2 = np.where(positive, h2, 1.0)

    # sqrt(|h|^2) is wanted three times over -- for the prefactor, for
    # the |h|^-3 in the arc-cosine argument, and for the degeneracy
    # scale below -- so it is taken once.  Writing that argument as a
    # division by |h|^2 sqrt(|h|^2), rather than as a power of -1.5, is
    # seven times quicker for the same value.
    sqrt_h2 = np.sqrt(safe_h2)
    pre = 2.0*SQRT3_INV*sqrt_h2

    # Not clipped in place: for a single Hamiltonian this is a scalar,
    # which has nowhere to write to
    chi = np.arccos(np.clip((-SQRT3)*h3/(safe_h2*sqrt_h2), -1.0, 1.0))

    m = np.array([1.0, 2.0, 3.0]).reshape((3,)+(1,)*np.ndim(chi))
    psi = np.where(positive, pre*np.cos((chi+2.0*np.pi*m)/3.0), 0.0)

    denom = 3.0*psi*psi - h2
    safe_denom = np.where(denom != 0.0, denom, 1.0)

    # The baselines enter only here
    exp_psi = np.exp(1.j*L*psi)

    # Only two combinations of the three terms survive the sum over m,
    #   u_k = i [ (sum_m w_m psi_m) h_k - (sum_m w_m) (h*h)_k ],
    # so they are formed directly.  Contracting the (..., 3, 8) array of
    # numerators instead would build, and then immediately discard, an
    # intermediate eight times the size of the result.
    w = exp_psi/safe_denom

    u0 = exp_psi.sum(0)/3.0

    if h[0].size == 1 and w.ndim > 1:
        # A single Hamiltonian scanned over baselines.  The bracket is
        # then a fixed 8-by-3 matrix and the sum over m is one matrix
        # product, which BLAS does better than three broadcasts.
        bracket = psi.reshape(1, 3)*h.reshape(8, 1) - star.reshape(8, 1)
        uk = 1.j*(bracket @ w.reshape(3, -1)).reshape((8,)+full)
    else:
        # The factor of i rides on the two (...)-shaped sums rather than
        # on the (8, ...) result, which is eight times larger
        weighted = 1.j*(w*psi).sum(0)
        total = 1.j*w.sum(0)
        uk = weighted*h - total*star

    # Recompute the elements whose spectrum is degenerate, using the same
    # criterion as the scalar path: the closest pair of latent roots.
    # Degeneracy is a property of the Hamiltonian, so the test is made at
    # the shape of `h` and then broadcast over the baselines.  Taking the
    # smallest gap pairwise avoids stacking the three of them into an
    # array only to reduce it away again.
    smallest_gap = np.minimum(np.minimum(np.abs(psi[0]-psi[1]),
                                         np.abs(psi[0]-psi[2])),
                              np.abs(psi[1]-psi[2]))
    degenerate = (smallest_gap <= DEGENERACY_TOL*sqrt_h2) | ~positive
    if degenerate.any():
        u0 = np.broadcast_to(u0, full).copy()
        uk = np.broadcast_to(uk, (8,)+full).copy()
        h_full = np.broadcast_to(h, (8,)+full)
        star_full = np.broadcast_to(star, (8,)+full)
        h2_full = np.broadcast_to(h2, full)
        h3_full = np.broadcast_to(h3, full)
        L_full = np.broadcast_to(L, full)
        for idx in zip(*np.nonzero(np.broadcast_to(degenerate, full))):
            column = (slice(None),)+idx
            values = _u_coefficients_3nu_single(
                h_full[column], float(h2_full[idx]), float(h3_full[idx]),
                float(L_full[idx]), star_full[column])
            u0[idx] = values[0]
            uk[column] = values[1:]

    return u0, uk


def _u_to_entries_batch(
    u0: np.ndarray,
    uk: np.ndarray
) -> Tuple[np.ndarray, ...]:
    r"""Returns the nine entries of :math:`U_3` from the coefficients.

    Parameters
    ----------
    u0 : numpy.ndarray
        The coefficient :math:`u_0`, of shape ``(...)``.
    uk : numpy.ndarray
        The coefficients :math:`u_1, \ldots, u_8`, of shape
        ``(8, ...)``.

    Returns
    -------
    tuple of numpy.ndarray
        The entries ``(U00, U01, U02, U10, U11, U12, U20, U21, U22)``,
        each of shape ``(...)``.
    """
    u1, u2, u3, u4, u5, u6, u7, u8 = uk
    u8_over_sqrt3 = u8/SQRT3

    return (u0+1.j*(u3+u8_over_sqrt3), 1.j*u1+u2, 1.j*u4+u5,
            1.j*u1-u2, u0-1.j*(u3-u8_over_sqrt3), 1.j*u6+u7,
            1.j*u4-u5, 1.j*u6-u7, u0-2.j*u8_over_sqrt3)


def _probabilities_3nu_batch(
    h_matrix: Union[list, np.ndarray],
    L: Union[int, float, list, np.ndarray]
) -> np.ndarray:
    r"""Returns the nine probabilities for a stack, without forming U.

    :math:`U_3` is needed only through the modulus squared of its
    entries, so the entries are squared as they are produced.  Stacking
    them into a ``(..., 3, 3)`` array first, and then transposing and
    reshaping it into the order the probabilities are returned in,
    allocates two further arrays and copies both.

    Parameters
    ----------
    h_matrix : array_like
        Hamiltonians, of shape ``(..., 3, 3)``.
    L : array_like
        Baselines, broadcastable against the leading axes of `h_matrix`.

    Returns
    -------
    numpy.ndarray
        The probabilities, of shape ``(..., 9)``, ordered with the
        initial flavor varying slowest.
    """
    h_matrix = np.asarray(h_matrix, dtype=complex)
    L = np.asarray(L, dtype=float)
    batch = np.broadcast_shapes(h_matrix.shape[:-2], L.shape)
    size = int(np.prod(batch, dtype=np.int64))

    if CHECK_HERMITICITY:
        _check_hermitian(h_matrix, 'probabilities_3nu')

    if size > 0 and fastkernels.worthwhile(3, size):
        return fastkernels.probabilities_3nu_kernel(
            np.broadcast_to(h_matrix, batch+(3, 3)),
            np.broadcast_to(L, batch))

    if size <= SMALL_BATCH:
        # Below this size the fixed cost of the array machinery exceeds
        # what the scalar path spends on the whole stack
        flat_h = np.broadcast_to(h_matrix, batch+(3, 3)).reshape(-1, 3, 3)
        flat_l = np.broadcast_to(L, batch).reshape(-1)
        out = np.empty((flat_l.shape[0], 9))
        for n in range(flat_l.shape[0]):
            # .tolist() gives Python complex numbers, on which the scalar
            # path is quicker than on NumPy scalars, by more than the
            # conversion costs
            out[n] = probabilities_3nu(flat_h[n].tolist(), float(flat_l[n]))
        return out.reshape(batch+(9,))

    entries = _u_to_entries_batch(*_u_coefficients_3nu_batch(
        _hamiltonian_3nu_coefficients_batch(h_matrix), L))

    # P_ab = |U_ba|^2: the entries come out row by row, so taking them
    # in column order is what puts the initial flavor slowest
    return np.stack([np.abs(entries[k])**2.
                     for k in (0, 3, 6, 1, 4, 7, 2, 5, 8)], axis=-1)


def _evolution_operator_3nu_batch(
    h_matrix: Union[list, np.ndarray],
    L: Union[int, float, list, np.ndarray]
) -> np.ndarray:
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

    if CHECK_HERMITICITY:
        _check_hermitian(h_matrix, 'evolution_operator_3nu')

    h = _hamiltonian_3nu_coefficients_batch(h_matrix)
    entries = _u_to_entries_batch(*_u_coefficients_3nu_batch(h, L))

    return np.stack([np.stack(entries[0:3], axis=-1),
                     np.stack(entries[3:6], axis=-1),
                     np.stack(entries[6:9], axis=-1)], axis=-2)


def _is_batched(
    hamiltonian_matrix: Union[list, np.ndarray],
    L: Union[int, float, list, np.ndarray]
) -> bool:
    r"""Returns whether the arguments describe a stack of problems.

    A single Hamiltonian is an ``n``-by-``n`` matrix and a single
    baseline is a scalar, so anything with more axes than that is a
    stack, and the vectorised path applies.

    This runs on every scalar call, so it is written to be cheap: an
    exact type check short-circuits the common case, and
    ``numpy.ndim`` --- which would convert a nested list to an array
    every time --- is reached only for an argument that is neither a
    plain Python number nor a NumPy array.

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Hamiltonian, or stack of them.
    L : int or float or array_like
        Baseline, or array of baselines.

    Returns
    -------
    bool
        Whether the vectorised path applies.
    """
    if type(L) is not float and type(L) is not int:
        if np.ndim(L) > 0:
            return True

    if type(hamiltonian_matrix) is np.ndarray:
        return hamiltonian_matrix.ndim > 2

    return isinstance(hamiltonian_matrix[0][0], (list, tuple, np.ndarray))


def evolution_operator_3nu_u_coefficients(
    hamiltonian_matrix: Union[list, np.ndarray],
    L: Union[int, float]
) -> List[complex]:
    r"""Returns the coefficients :math:`u_0, \ldots, u_8`.

    Returns the nine coefficients :math:`u_0, \ldots, u_8` of the
    three-neutrino time-evolution operator :math:`U_3(L)` in its SU(3)
    exponential expansion,
    :math:`U_3 = u_0 \mathbb{1} + i u_k \lambda^k`.

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.1.0
       Degenerate Hamiltonians are handled exactly instead of returning
       NaN.  Lagrange interpolation over the latent roots divides by
       :math:`3\psi_m^2 - |h|^2`, which vanishes at a repeated root; the
       two degenerate cases are now taken in their confluent limit.

    .. versionchanged:: 1.4.0
       Faster, with identical results --- all 42 figures generated by
       ``run_testsuite.py`` are byte-for-byte those of 1.3.0. The sum
       over the latent roots forms its :math:`k`-independent factors
       once instead of inside the loop over the eight :math:`k`, and the
       star product is computed once and passed in rather than
       recomputed.

    .. versionchanged:: 1.5.0
       Faster, with identical results; the probabilities agree with
       1.4.0 to 1.6e-13 across every code path.  Around the latent
       roots, :math:`\sqrt{|h|^2}` is taken once rather than three
       times, the arc-cosine argument is a division rather than a power
       of -1.5, and the degeneracy test uses two minima instead of
       stacking three gap arrays.  The scalar exponentials use
       :func:`cmath.rect`.

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
    .. jupyter-execute::

        import oscprob3nu

        hamiltonian_matrix = [[1.0+0.0j, 0.0+2.0j, 0.0-1.0j],
                              [0.0-2.0j, 3.0+0.0j, 3.0+0.0j],
                              [0.0+1.0j, 3.0-0.0j, -5.0+0.0j]]
        u_coeffs = oscprob3nu.evolution_operator_3nu_u_coefficients(hamiltonian_matrix,
                                                         1.0)
        print('%+.6f%+.6fj' % (u_coeffs[0].real, u_coeffs[0].imag))
    """
    # [h1, h2, h3, h4, h5, h6, h7, h8]
    h = hamiltonian_3nu_coefficients(hamiltonian_matrix)

    # The star product gives <h> = h_i (h*h)_i, so computing it here and
    # passing it on spares the expansion a second identical contraction
    star_coeffs = _star_all(h)
    h2 = sum([x*x for x in h])
    h3 = sum([x*y for x, y in zip(h, star_coeffs)])

    return _u_coefficients_3nu_single(h, h2, h3, L, star_coeffs)


def evolution_operator_3nu(
    hamiltonian_matrix: Union[list, np.ndarray],
    L: Union[int, float, list, np.ndarray]
) -> Union[List[List[complex]], np.ndarray]:
    r"""Returns the three-neutrino time-evolution operator.

    Returns the three-neutrino time-evolution operator :math:`U_3(L)` in
    its SU(3) exponential expansion
    :math:`U_3(L) = u_0 \mathbb{1} + i u_k \lambda^k`.  This is a
    :math:`3\times3` unitary matrix.

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.1.0
       Degenerate Hamiltonians are handled exactly instead of returning
       NaN.  Lagrange interpolation over the latent roots divides by
       :math:`3\psi_m^2 - |h|^2`, which vanishes at a repeated root; the
       two degenerate cases are now taken in their confluent limit.

    .. versionchanged:: 1.2.0
       Accepts a stack of Hamiltonians of shape ``(..., n, n)``, an
       array of baselines, or both broadcast against each other,
       returning an array with the broadcast leading axes.  A single
       Hamiltonian with a scalar baseline returns exactly what it
       returned before.

    .. versionchanged:: 1.4.0
       Faster, with identical results --- all 42 figures generated by
       ``run_testsuite.py`` are byte-for-byte those of 1.3.0. The sum
       over the three latent roots is rewritten to form its
       :math:`k`-independent factors once, removing an ``(N, 3, 8)``
       intermediate from the batched path.

    .. versionchanged:: 1.5.0
       Faster, with identical results; the probabilities agree with
       1.4.0 to 1.6e-13 across every code path.  The batched star
       product uses the sparse expansion rather than an
       :func:`numpy.einsum` contraction of the dense :math:`d` tensor.
       The batched coefficients are laid out with the component index
       first, so that each :math:`h_k` is contiguous; every downstream
       step works one component at a time, and reading those from a
       strided view cost about a third more.

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
    .. jupyter-execute::

        import oscprob3nu

        hamiltonian_matrix = [[1.0+0.0j, 0.0+2.0j, 0.0-1.0j],
                              [0.0-2.0j, 3.0+0.0j, 3.0+0.0j],
                              [0.0+1.0j, 3.0-0.0j, -5.0+0.0j]]
        U3 = oscprob3nu.evolution_operator_3nu(hamiltonian_matrix, 1.0)
        for row in U3:
            print('  '.join(['%+.6f%+.6fj' % (z.real+0.0, z.imag+0.0)
                             for z in row]))
    """
    if _is_batched(hamiltonian_matrix, L):
        return _evolution_operator_3nu_batch(hamiltonian_matrix, L)

    if CHECK_HERMITICITY:
        _check_hermitian(np.asarray(hamiltonian_matrix, dtype=complex),
                         'evolution_operator_3nu')

    u0, u1, u2, u3, u4, u5, u6, u7, u8 = \
        evolution_operator_3nu_u_coefficients(hamiltonian_matrix, L)

    return [
        [u0+1.j*(u3+u8/SQRT3), 1.j*u1+u2, 1.j*u4+u5],
        [1.j*u1-u2, u0-1.j*(u3-u8/SQRT3), 1.j*u6+u7],
        [1.j*u4-u5, 1.j*u6-u7, u0-1.j*2.*u8/SQRT3]
    ]


def probabilities_3nu(
    hamiltonian_matrix: Union[list, np.ndarray],
    L: Union[int, float, list, np.ndarray]
) -> Union[Tuple[float, ...], np.ndarray]:
    r"""Returns the three-neutrino oscillation probabilities.

    Returns the three-neutrino flavor-transition probabilities
    :math:`P_{ee}, P_{e\mu}, P_{e\tau}, P_{\mu e}, P_{\mu\mu},
    P_{\mu\tau}, P_{\tau e}, P_{\tau\mu}, P_{\tau\tau}`, where
    :math:`P_{\alpha\beta} \equiv P(\nu_\alpha \to \nu_\beta)
    = |[U_3]_{\beta\alpha}|^2`.

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.1.0
       Degenerate Hamiltonians are handled exactly instead of returning
       NaN.  Lagrange interpolation over the latent roots divides by
       :math:`3\psi_m^2 - |h|^2`, which vanishes at a repeated root; the
       two degenerate cases are now taken in their confluent limit.

    .. versionchanged:: 1.2.0
       Accepts a stack of Hamiltonians of shape ``(..., n, n)``, an
       array of baselines, or both broadcast against each other,
       returning an array with the broadcast leading axes.  A single
       Hamiltonian with a scalar baseline returns exactly what it
       returned before.

    .. versionchanged:: 1.4.0
       Faster, with identical results --- all 42 figures generated by
       ``run_testsuite.py`` are byte-for-byte those of 1.3.0. The sum
       over the three latent roots is rewritten so that the two factors
       independent of :math:`k` are formed once rather than inside the
       loop over the eight :math:`k`, which also removes an ``(N, 3,
       8)`` intermediate from the batched path; the star product is
       computed once per call rather than twice; and :math:`|z|^2` skips
       the square root that squaring undoes.  A scalar three-flavor
       probability is 3.2x quicker.

    .. versionchanged:: 1.5.0
       Faster, with identical results; the probabilities agree with
       1.4.0 to 1.6e-13 across every code path.  The batched star
       product no longer contracts the dense 8x8x8 :math:`d` tensor
       through :func:`numpy.einsum`, which with no path plan walked the
       whole table for every element and was 70% of a 2000-point energy
       scan; the sparse expansion the scalar path already used
       vectorises unchanged.  The routine also forms the nine entries
       and squares them directly instead of building and transposing the
       evolution operator.  The batched coefficients are laid out with
       the component index first, so that each :math:`h_k` is
       contiguous; every downstream step works one component at a time,
       and reading those from a strided view cost about a third more.
       The scalar path is 3.2x quicker, a 2000-point energy scan 3.5x,
       and an oscillogram 1.8x.

    .. versionchanged:: 1.6.0
       Two dispatch decisions, neither changing the result.  A stack of
       at most :data:`SMALL_BATCH` elements is evaluated one at a time
       through the scalar path, because a batched call carries a fixed
       cost that a handful of points does not amortise.  Larger stacks
       are evaluated by a compiled kernel instead, if the optional
       `numba` extra is installed; for three flavors that kernel wins at
       every size.  See :mod:`fastkernels`.

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
    .. jupyter-execute::

        import oscprob3nu

        hamiltonian_matrix = [[1.0+0.0j, 0.0+2.0j, 0.0-1.0j],
                              [0.0-2.0j, 3.0+0.0j, 3.0+0.0j],
                              [0.0+1.0j, 3.0-0.0j, -5.0+0.0j]]
        prob = oscprob3nu.probabilities_3nu(hamiltonian_matrix, 1.0)
        print('  '.join(['%.6f' % p for p in prob[0:3]]))
        print('  '.join(['%.6f' % p for p in prob[3:6]]))
        print('  '.join(['%.6f' % p for p in prob[6:9]]))
    """
    if _is_batched(hamiltonian_matrix, L):
        return _probabilities_3nu_batch(hamiltonian_matrix, L)

    if CHECK_HERMITICITY:
        _check_hermitian(np.asarray(hamiltonian_matrix, dtype=complex),
                         'probabilities_3nu')

    U = evolution_operator_3nu(hamiltonian_matrix, L)

    row_e, row_m, row_t = U

    # P_ab = |U_ba|^2: the evolution operator is indexed (final, initial)
    return (_abs2(row_e[0]), _abs2(row_m[0]), _abs2(row_t[0]),
            _abs2(row_e[1]), _abs2(row_m[1]), _abs2(row_t[1]),
            _abs2(row_e[2]), _abs2(row_m[2]), _abs2(row_t[2]))
