# -*- coding: utf-8 -*-
r"""Compute the four-neutrino flavor-transition probabilities.

This module contains the routines needed to compute four-neutrino
flavor-transition probabilities for an arbitrary time-independent
:math:`4\times4` Hermitian Hamiltonian, using the SU(4) exponential
expansion.  It is the :math:`n = 4` member of the family whose
:math:`n = 2` and :math:`n = 3` members are :mod:`oscprob2nu` and
:mod:`oscprob3nu`, following the method of [1]_.  It is the last one: see
:ref:`why-four-is-the-end` in the methodology page.  The Cayley-Hamilton
route at :math:`n = 4` in constant matter was set out in [2]_.

The Hamiltonian is expanded in the basis of the fifteen generalized
Gell-Mann matrices,

.. math:: H = h_0 \mathbb{1} + h_a \lambda^a , \qquad a = 1, \ldots, 15 ,

and the time-evolution operator in the same basis,

.. math:: U_4(L) = u_0 \mathbb{1} + i u_a \lambda^a .

The term :math:`h_0` contributes only an overall phase and is dropped;
all routines therefore work with the traceless part of the Hamiltonian,
which leaves the oscillation probabilities unchanged.

What is new at four flavors
---------------------------

Three ingredients have no counterpart in :mod:`oscprob3nu`:

* **A third invariant.**  SU(4) has rank three, so the traceless part
  carries three independent invariants rather than two,

  .. math::
     I_2 = \tfrac12 \mathrm{Tr}\,\tilde{H}^2 , \qquad
     I_3 = \tfrac12 \mathrm{Tr}\,\tilde{H}^3 , \qquad
     I_4 = \tfrac12 \left(\mathrm{Tr}\,\tilde{H}^4 - I_2^2\right) ,

  and the cubic characteristic equation of the three-flavor case becomes
  the quartic

  .. math::
     \psi^4 - I_2 \psi^2 - \tfrac23 I_3 \psi
     + \tfrac14 \left(I_2^2 - 2 I_4\right) = 0 .

* **A quartic that still solves in closed form.**  Euler's method
  reduces it to a *resolvent cubic* whose three roots are real and
  non-negative because :math:`\tilde{H}` is Hermitian, so the same
  trigonometric formula that :func:`oscprob3nu.psi_roots` uses solves it.
  The SU(3) machinery is literally nested inside the SU(4) solution.

* **A longer star-product tower.**  The three-flavor identity
  :math:`(h \star h) \star h = \tfrac13 |h|^2 h` is a Cayley-Hamilton
  accident of :math:`n = 3` and is false for SU(4) --- it is off by some
  tens of per cent on a random Hamiltonian --- so
  :math:`((h \star h) \star h)_a` enters as independent data.

Accuracy
--------

For a generic Hamiltonian the expansion is exact to round-off, like its
two- and three-flavor counterparts.  A *stiff* spectrum is the case that
needs care, and the physically interesting 3+1 scenario is stiff: with
:math:`\Delta m^2_{41} \sim 1` eV\ :sup:`2` the eigenvalues span four
orders of magnitude, three of them clustering, and the information
needed to separate the cluster is destroyed when :math:`I_2, I_3, I_4`
are formed in double precision.  No amount of care in solving the
quartic recovers it.

The roots are therefore refined against the *matrix*, by one Newton step
on :math:`\chi(\psi) = \det(\psi \mathbb{1} - \tilde{H})`, which uses the
Hamiltonian entries rather than the three compressed invariants.  That
restores the roots to round-off and the probabilities from about
:math:`5 \times 10^{-7}` to :math:`10^{-9}`.

Neither figure is anywhere near a measurable effect --- probabilities
are confronted with data at the per-cent level --- so this matters for
the exactness claim, for error accumulating when :mod:`slabs` and
:mod:`earth` compose operators across layers, and for a regression suite
tight enough to catch a real mistake.  :data:`POLISH_ROOTS` has the
measured comparison against the alternatives, including why LAPACK and
extended precision both lose.

`evolution_operator_4nu` and `probabilities_4nu` accept either a single
Hamiltonian and baseline or a stack of them, in which case the whole
stack is evaluated at once, exactly as at two and three flavors.

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

    * generators_su4 - Returns the fifteen generalized Gell-Mann matrices
    * hamiltonian_4nu_coefficients - Returns the :math:`h_a`
    * su4_invariants - Returns the invariants :math:`I_2, I_3, I_4`
    * psi_roots_4nu - Returns the roots of the characteristic equation
    * evolution_operator_4nu_u_coefficients - Returns the :math:`u_a`
    * evolution_operator_4nu - Returns the evolution operator :math:`U_4`
    * probabilities_4nu - Returns the oscillation probabilities

References
----------

.. [1] Mauricio Bustamante, "Exact neutrino oscillation probabilities
   with arbitrary time-independent Hamiltonians", arXiv:1904.12391.

.. [2] S. Kamo et al., "Matter enhanced transitions of active and
   sterile neutrinos", Eur. Phys. J. C 28, 211 (2003), which applies the
   Cayley-Hamilton approach at :math:`n = 4` in constant matter.
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

__all__ = ['CHECK_HERMITICITY', 'LAMBDA_SU4', 'POLISH_ROOTS', 'SMALL_BATCH',
           'generators_su4', 'hamiltonian_4nu_coefficients',
           'su4_invariants', 'psi_roots_4nu',
           'evolution_operator_4nu_u_coefficients', 'evolution_operator_4nu',
           'probabilities_4nu']

from typing import List, Tuple, Union

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
    # The backend is a dependency, but this module is written to run without
    # it -- a lone copy like this one, or an environment where it was removed
    # on purpose.  Its absence is answered the same way its being switched
    # off is: `worthwhile` says no, and the NumPy path runs.  Nothing else in
    # this module touches it.
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


SMALL_BATCH = 10
r"""int: Module-level constant.

Nothing reads this.  At two and three flavors the constant of the same
name selects between a scalar path and a batched one; four flavors has
no separate scalar closed form to select, so every stack goes through
the array path whatever its length, and there is no short-stack case for
a threshold to name.  It is exported and documented as being what it is,
rather than quietly kept.

It also no longer matches :data:`oscprob3nu.SMALL_BATCH`, which it once
did and was documented as doing: that one was re-measured when the
scalar path was given a cheap Hermiticity check, and this one, governing
nothing, was not.
"""

POLISH_ROOTS = True
r"""bool: Module-level switch.

Whether to refine the latent roots against the Hamiltonian matrix, by
one Newton step on :math:`\chi(\psi) = \det(\psi \mathbb{1} - \tilde{H})`,
after solving the quartic in closed form.

This is on by default and should stay on.  The closed-form roots are
limited by the conditioning of :math:`I_2, I_3, I_4`, not by the quartic
solver: perturbing the three invariants at the :math:`10^{-16}` level
moves the roots of a stiff 3+1 spectrum by :math:`6 \times 10^{-11}`
relative, which is what the unrefined closed form achieves and what no
better root-finder can beat.  The Newton step reads the matrix entries
instead, and is not subject to that floor.

Measured against ``mpmath`` at fifty decimal digits, on stiff 3+1
Hamiltonians, with cost quoted for a 200 000-point scan:

======================================  ==========  =======
Strategy for the latent roots           Rel. error  Cost
======================================  ==========  =======
Closed form alone                       8.3e-11     0.17 s
Closed form + one Newton step           1.1e-16     0.41 s
``numpy.linalg.eigvalsh``               7.4e-16     0.17 s
Closed form in ``numpy.longdouble``     4.5e-11     0.43 s
======================================  ==========  =======

Note the second row beats the third: the Newton step is some seven times
more accurate than LAPACK, because ``eigvalsh`` reduces by similarity
transforms that each carry a backward error of order
:math:`\epsilon\|H\|`, while this converges onto the root of
:math:`\det(\psi\mathbb{1} - \tilde{H})` for the matrix it was given.
Extended precision was rejected for buying under a digit, being slower,
and silently being ``float64`` on Apple Silicon and Windows.

The refined figure has been confirmed from outside: against `nuSQuIDS
<https://github.com/arguelles/nuSQuIDS>`_, which integrates the density
matrix numerically, the four-flavor probabilities agree to
:math:`4 \times 10^{-16}` on a benign spectrum and
:math:`3 \times 10^{-10}` on the stiffest one tested.  See notebook 17.

In probabilities the difference is :math:`5 \times 10^{-7}` unrefined
against :math:`10^{-9}` refined.  Both are orders of magnitude below what
any experiment resolves; the reasons to want the smaller one are the
library's claim to exactness, error accumulating when :mod:`slabs` and
:mod:`earth` compose operators across many layers, and a regression
suite with no room for a bug to hide in.

A second Newton step changes nothing --- one already reaches the floor
--- so exactly one is taken.  It costs about 40% of the runtime, which
brings the four-flavor closed form to parity with a batched ``eigh``
rather than ahead of it.

The step is applied unconditionally rather than only where a spectrum
looks stiff, and that is a measured decision rather than a lazy one.
Two skip criteria were tried on 6300 Hamiltonians: the gap-based
amplification that perturbation theory suggests, which misjudges doubly
paired spectra by four orders of magnitude, and a matrix residual
comparing :math:`\prod_m \psi_m` with :math:`\det \tilde{H}`, which is
one constraint on four roots and misses errors that cancel in the
product.  Neither can safely skip any elements at all.  The reason is
structural: a criterion complete enough to certify four roots must
evaluate :math:`\chi` at four roots, which is this refinement --- the
check and the fix are the same computation.  See the methodology page.

Set to ``False`` to skip it, which is useful only for reproducing the
unrefined figures or for spectra known to be well separated.
"""

DEGENERACY_TOL = 1.e-12
r"""float: Module-level constant.

Relative tolerance below which two latent roots :math:`\psi` are treated
as degenerate.  Reconstructing :math:`U_4` from the roots divides by
their differences, which vanish at a repeated root, so a degenerate
spectrum is handled by a separate, exact expression.
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
here instead.  On a *single* matrix the check costs about 1.35x, and it
takes its own branch to manage that: the reductions above are all fixed
cost at one element, so without the branch one probability spent nine
tenths of its time here.  For production scans whose Hamiltonians come from a
construction already known to be Hermitian --- everything
:mod:`hamiltonians4nu` builds is Hermitian to round-off, as the table
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

    Compares the 6 independent pairs and the imaginary parts of the
    diagonal, rather than forming ``H - H^dagger``, which would allocate
    a temporary the size of the whole stack.  Real and imaginary parts
    are compared separately, on views rather than copies: the condition
    is :math:`\mathrm{Re}\,H_{ij} = \mathrm{Re}\,H_{ji}` together
    with :math:`\mathrm{Im}\,H_{ij} = -\mathrm{Im}\,H_{ji}`, which
    needs no complex arithmetic and, unlike :func:`numpy.abs` on a
    complex array, no square root per element.  That is worth about a
    factor of three on a large stack, where this check would otherwise
    cost several times the evaluation it guards.

    That describes the stack.  A single matrix takes its own branch and
    compares its entries as Python complex numbers, because the
    reductions above are all fixed cost at one element and would
    otherwise dominate a scalar probability.

    Parameters
    ----------
    h_matrix : numpy.ndarray
        Hamiltonian, or stack of them, of shape ``(..., 4, 4)``.
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

    non_finite = (
        '%s: the Hamiltonian has a non-finite entry, so it is neither '
        'Hermitian nor usable.  Set %s.CHECK_HERMITICITY = False to '
        'skip this check.')

    complaint = (
        '%s: the Hamiltonian is not Hermitian%s.  The expansion assumes '
        'Hermiticity, and without it the probabilities returned are '
        'meaningless even though they still sum to one.  Set '
        'oscprob4nu.CHECK_HERMITICITY = False to skip this check.')


    # A single matrix takes its own path, because everything below is
    # fixed cost at this size.  The reductions run about sixty
    # microseconds on one matrix against half a microsecond per element
    # on a stack of two thousand, and a scalar probability costs eight
    # microseconds in total --- so when this check arrived in 1.11.0 it
    # became nine tenths of the work, and made short stacks slower than
    # the batched path they exist to avoid.  Comparing the entries as
    # Python complex numbers is one conversion and a few scalar
    # operations, and reaches the same verdict.
    if h_matrix.ndim == 2:
        entries = h_matrix.tolist()

        # Tested per entry rather than by letting a non-finite value
        # reach `scale`: `max` keeps its running value when the
        # comparison is with a nan, since every comparison against one is
        # false, so a nan would never arrive and `isfinite` would pass.
        # The array path has no such hole, `np.max` propagating nan --- so
        # the first draft of this branch returned probabilities for a
        # Hamiltonian the batched path refuses, which is the divergence
        # this check exists to prevent.
        scale = 0.0
        for row in entries:
            for entry in row:
                real, imaginary = abs(entry.real), abs(entry.imag)
                if not (math.isfinite(real) and math.isfinite(imaginary)):
                    raise ValueError(non_finite % (caller, 'oscprob4nu'))
                scale = max(scale, real, imaginary)

        tolerance = (_HERMITICITY_TOL*scale if scale > 0.0
                     else _HERMITICITY_TOL)

        for i in range(4):
            if abs(entries[i][i].imag) > tolerance:
                raise ValueError(complaint % (
                    caller, ' --- the diagonal entry (%d, %d) has a non-zero '
                    'imaginary part' % (i, i)))
            for j in range(i+1, 4):
                upper, lower = entries[i][j], entries[j][i]
                if (abs(upper.real - lower.real) > tolerance
                        or abs(upper.imag + lower.imag) > tolerance):
                    raise ValueError(complaint % (
                        caller, ' --- entry (%d, %d) is not the complex '
                        'conjugate of entry (%d, %d)' % (i, j, j, i)))
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
        raise ValueError(non_finite % (caller, 'oscprob4nu'))

    tolerance = _HERMITICITY_TOL*scale if scale > 0.0 else _HERMITICITY_TOL

    for i in range(4):
        if np.any(np.abs(imaginary[..., i, i]) > tolerance):
            raise ValueError(complaint % (
                caller, ' --- the diagonal entry (%d, %d) has a non-zero '
                'imaginary part' % (i, i)))
        for j in range(i+1, 4):
            if (np.any(np.abs(real[..., i, j] - real[..., j, i]) > tolerance)
                    or np.any(np.abs(imaginary[..., i, j]
                                     + imaginary[..., j, i]) > tolerance)):
                raise ValueError(complaint % (
                    caller, ' --- entry (%d, %d) is not the complex '
                    'conjugate of entry (%d, %d)' % (i, j, j, i)))



def generators_su4() -> np.ndarray:
    r"""Returns the fifteen generalized Gell-Mann matrices of SU(4).

    The generators are traceless, Hermitian, and normalized as
    :math:`\mathrm{Tr}(\lambda^a \lambda^b) = 2 \delta^{ab}`, which is
    the convention :mod:`oscprob3nu` uses at :math:`n = 3`.

    They are ordered as six symmetric and six antisymmetric off-diagonal
    matrices, one pair for each of the six index pairs
    :math:`(j, k)` with :math:`j < k` in lexicographic order, followed by
    the three diagonal matrices that span the Cartan subalgebra.

    .. versionadded:: 1.9.0

    Returns
    -------
    numpy.ndarray
        Complex array of shape ``(15, 4, 4)``, the generators
        :math:`\lambda^1, \ldots, \lambda^{15}`.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import oscprob4nu

        lam = oscprob4nu.generators_su4()
        gram = np.einsum('aij,bji->ab', lam, lam).real/2.0

        print('shape:', lam.shape)
        print('orthonormal to %.1e' % np.max(np.abs(gram - np.eye(15))))
    """
    generators = []

    for j in range(4):
        for k in range(j+1, 4):
            symmetric = np.zeros((4, 4), dtype=complex)
            symmetric[j, k] = symmetric[k, j] = 1.0
            antisymmetric = np.zeros((4, 4), dtype=complex)
            antisymmetric[j, k] = -1.j
            antisymmetric[k, j] = 1.j
            generators += [symmetric, antisymmetric]

    for level in range(1, 4):
        diagonal = np.zeros((4, 4), dtype=complex)
        diagonal[np.arange(level), np.arange(level)] = 1.0
        diagonal[level, level] = -float(level)
        generators.append(np.sqrt(2.0/(level*(level+1)))*diagonal)

    return np.asarray(generators)


LAMBDA_SU4 = generators_su4()
r"""numpy.ndarray: Module-level constant.

The fifteen generalized Gell-Mann matrices, of shape ``(15, 4, 4)``,
tabulated once at import time.  At three flavors the analogous table is
the :math:`d` tensor of :func:`oscprob3nu.tensor_d`; here the generators
themselves are tabulated, because every quantity this module needs is a
trace against them and none needs :math:`d_{abc}` explicitly.
"""

_TRACELESS_SUBTRACTION = np.eye(4)


def hamiltonian_4nu_coefficients(
    hamiltonian_matrix: Union[list, np.ndarray]
) -> List[float]:
    r"""Returns the :math:`h_a` of the SU(4) expansion of the Hamiltonian.

    Computes the coefficients :math:`h_1, \ldots, h_{15}` of the SU(4)
    expansion :math:`H = h_0 \mathbb{1} + h_a \lambda^a` of the
    four-flavor Hamiltonian `hamiltonian_matrix`, which is assumed to be
    given in the flavor basis.  The coefficient :math:`h_0` contributes
    only an overall phase to the evolution operator and is not returned.

    They follow from the normalization of the generators as
    :math:`h_a = \tfrac12 \mathrm{Tr}(H \lambda^a)`, and are real for a
    Hermitian Hamiltonian.

    .. versionadded:: 1.9.0

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Four-flavor Hermitian Hamiltonian, given as a nested list or an
        array of shape ``(4, 4)``.

    Returns
    -------
    list of float
        The fifteen coefficients :math:`h_1, \ldots, h_{15}`.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import oscprob4nu

        hamiltonian = np.diag([1.0, 2.0, 3.0, -6.0]).astype(complex)
        h = oscprob4nu.hamiltonian_4nu_coefficients(hamiltonian)

        print('%d coefficients' % len(h))
        print('the three Cartan ones: %.4f, %.4f, %.4f' % tuple(h[12:]))
    """
    matrix = np.asarray(hamiltonian_matrix, dtype=complex)
    traceless = matrix - (np.trace(matrix).real/4.0)*_TRACELESS_SUBTRACTION

    return list(0.5*np.einsum('aij,ji->a', LAMBDA_SU4, traceless).real)


def su4_invariants(
    hamiltonian_matrix: Union[list, np.ndarray]
) -> Tuple[float, float, float]:
    r"""Returns the three SU(4) invariants of the Hamiltonian.

    Returns :math:`I_2 = |h|^2`, :math:`I_3 = \langle h \rangle` and
    :math:`I_4 = |h \star h|^2` of the traceless part
    :math:`\tilde{H}` of `hamiltonian_matrix`, computed from traces of
    its powers,

    .. math::
       I_2 = \tfrac12 \mathrm{Tr}\,\tilde{H}^2 , \qquad
       I_3 = \tfrac12 \mathrm{Tr}\,\tilde{H}^3 , \qquad
       I_4 = \tfrac12 \left(\mathrm{Tr}\,\tilde{H}^4 - I_2^2\right) .

    SU(4) has rank three, so there are three of them; SU(3) has rank two,
    which is why :func:`oscprob3nu.su3_invariants` returns two.  Taking
    them from traces avoids ever building the :math:`d` tensor of SU(4),
    a :math:`15 \times 15 \times 15` table that nothing else here needs.

    .. versionadded:: 1.9.0

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Four-flavor Hermitian Hamiltonian, given as a nested list or an
        array of shape ``(4, 4)``.

    Returns
    -------
    tuple of float
        The invariants ``(I2, I3, I4)``.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import oscprob4nu

        hamiltonian = np.diag([1.0, 2.0, 3.0, -6.0]).astype(complex)
        i2, i3, i4 = oscprob4nu.su4_invariants(hamiltonian)

        print('I2 = %.4f' % i2)
        print('I3 = %.4f' % i3)
        print('I4 = %.4f' % i4)
    """
    matrix = np.asarray(hamiltonian_matrix, dtype=complex)
    traceless = matrix - (np.trace(matrix).real/4.0)*_TRACELESS_SUBTRACTION

    squared = traceless @ traceless
    invariant_2 = 0.5*np.trace(squared).real
    invariant_3 = 0.5*np.trace(squared @ traceless).real
    invariant_4 = 0.5*(np.einsum('ij,ji->', squared, squared).real
                       - invariant_2*invariant_2)

    return invariant_2, invariant_3, invariant_4


def _resolvent_cubic_roots(
    coeff_2: Union[int, float],
    coeff_1: Union[int, float],
    coeff_0: Union[int, float]
) -> np.ndarray:
    r"""Returns the three real roots of a depressed-friendly cubic.

    Solves :math:`z^3 + c_2 z^2 + c_1 z + c_0 = 0` by the trigonometric
    (Viete) construction, which is the same one
    :func:`oscprob3nu.psi_roots` uses at three flavors.  All three roots
    are real for the resolvent cubic of a Hermitian Hamiltonian.

    Parameters
    ----------
    coeff_2 : int or float
        Coefficient of :math:`z^2`.
    coeff_1 : int or float
        Coefficient of :math:`z`.
    coeff_0 : int or float
        Constant coefficient.

    Returns
    -------
    numpy.ndarray
        The three roots, of shape ``(3,)``.
    """
    depressed_p = coeff_1 - coeff_2*coeff_2/3.0
    depressed_q = (2.0*coeff_2**3/27.0 - coeff_2*coeff_1/3.0 + coeff_0)
    shift = -coeff_2/3.0

    if depressed_p >= 0.0:
        # A triple root, up to round-off: the three roots coincide and
        # the trigonometric form degenerates.
        return np.full(3, np.cbrt(-depressed_q)) + shift

    scale = 2.0*np.sqrt(-depressed_p/3.0)
    argument = np.clip(3.0*depressed_q/(depressed_p*scale), -1.0, 1.0)
    angle = np.arccos(argument)

    return scale*np.cos((angle + 2.0*np.pi*np.arange(3))/3.0) + shift


def psi_roots_4nu(
    invariant_2: Union[int, float],
    invariant_3: Union[int, float],
    invariant_4: Union[int, float]
) -> List[float]:
    r"""Returns the roots of the four-flavor characteristic equation.

    Returns the four real roots :math:`\psi_m` of

    .. math::
       \psi^4 - I_2 \psi^2 - \tfrac23 I_3 \psi
       + \tfrac14 \left(I_2^2 - 2 I_4\right) = 0 ,

    the characteristic equation of the traceless part of the
    Hamiltonian, obtained by Euler's solution of the quartic: the
    resolvent cubic

    .. math:: z^3 - 2 I_2 z^2 + 2 I_4 z - \tfrac49 I_3^2 = 0

    has roots :math:`z_i = (\psi_i + \psi_j)^2`, real and non-negative
    because the Hamiltonian is Hermitian, and then

    .. math::
       \psi_m = \tfrac12 \left(s_1 \sqrt{z_1} + s_2 \sqrt{z_2}
                               + s_3 \sqrt{z_3}\right) ,

    with signs :math:`s_i = \pm1` fixed by
    :math:`s_1 s_2 s_3 \sqrt{z_1 z_2 z_3} = \tfrac23 I_3`.

    These roots carry only the accuracy that :math:`I_2, I_3, I_4`
    carry.  For a stiff spectrum that is not enough, which is why the
    routines that use them refine them against the matrix; see
    :data:`POLISH_ROOTS`.

    .. versionadded:: 1.9.0

    Parameters
    ----------
    invariant_2 : int or float
        The invariant :math:`I_2 = |h|^2`.
    invariant_3 : int or float
        The invariant :math:`I_3 = \langle h \rangle`.
    invariant_4 : int or float
        The invariant :math:`I_4 = |h \star h|^2`.

    Returns
    -------
    list of float
        The four roots, in ascending order.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import oscprob4nu

        hamiltonian = np.diag([1.0, 2.0, 3.0, -6.0]).astype(complex)
        i2, i3, i4 = oscprob4nu.su4_invariants(hamiltonian)
        psi = oscprob4nu.psi_roots_4nu(i2, i3, i4)

        print('roots    :', np.round(psi, 6))
        print('reference:', np.round(np.linalg.eigvalsh(hamiltonian) - 0.0, 6))
    """
    quadratic = -invariant_2
    linear = -(2.0/3.0)*invariant_3
    constant = 0.25*(invariant_2*invariant_2 - 2.0*invariant_4)

    resolvent = _resolvent_cubic_roots(2.0*quadratic,
                                       quadratic*quadratic - 4.0*constant,
                                       -linear*linear)
    root = np.sqrt(np.clip(resolvent, 0.0, None))

    if linear > 0.0:
        root = root*np.array([1.0, 1.0, -1.0])

    return list(np.sort(0.5*np.array([
        root[0] + root[1] + root[2],
        root[0] - root[1] - root[2],
        -root[0] + root[1] - root[2],
        -root[0] - root[1] + root[2],
    ])))


def _polish_roots(
    traceless: np.ndarray,
    psi: np.ndarray
) -> np.ndarray:
    r"""Refines the latent roots against the Hamiltonian matrix.

    One Newton step on :math:`\chi(\psi) = \det(\psi\mathbb{1} -
    \tilde{H})`, with the derivative taken as
    :math:`\chi'(\psi_m) = \prod_{l \neq m}(\psi_m - \psi_l)` from the
    current estimates.  The determinant is evaluated from the matrix, so
    unlike the closed-form roots it is not limited by the conditioning of
    the three invariants.

    Works on a stack: `traceless` of shape ``(..., 4, 4)`` and `psi` of
    shape ``(..., 4)``.

    A Newton step is only taken where it is meaningful.  The derivative
    is a product of gaps, so a root that nearly coincides with another
    divides by a very small number and is thrown across the spectrum ---
    a pair separated by one unit in the last place gives a derivative of
    order :math:`10^{-16}` and a step of order one.  The guard is the
    standard one for polishing polynomial roots: a step for a simple
    root should never move it more than halfway to its nearest
    neighbour, and a step that wants to is evidence that the root
    belongs to a cluster, where the closed form is already the best
    estimate available and refining against a nearly singular
    :math:`\chi'` only destroys it.

    Testing that ``derivative != 0.0``, which is what this did before,
    is the same guard with its threshold at exactly zero, and that is a
    knife edge: whether a degenerate pair lands on identical bits or on
    adjacent ones is decided by the last bit of a square root taken near
    zero, and it changes between a stack and a scalar call, between
    machines, and between this and the compiled kernel.  On synthetic
    spectra with pair separations drawn between :math:`10^{-16}` and
    :math:`10^{-6}` relative, that version returned roots wrong by more
    than :math:`10^{-6}` relative for about 10% of them, the worst by
    4.8 --- a root thrown clean across a spectrum of order one.

    What the guard does *not* do is make such a spectrum accurate, and
    the distinction matters.  A nearly degenerate pair is poorly
    resolved by the closed form to begin with: Euler's reduction
    recovers the pair's separation as :math:`\sqrt{z}` for a resolvent
    root :math:`z` that vanishes as the pair closes, so the
    :math:`\epsilon` carried by :math:`z` becomes :math:`\sqrt{\epsilon}`
    on the separation, of order :math:`10^{-8}` relative.  No Newton step
    against :math:`\chi` recovers that, and the guarded step does not
    try.  What it guarantees is the weaker and correct thing: refining
    never leaves the roots worse than the closed form left them, and a
    test asserts that ordering case by case.  Over the same sweep the
    guard refuses the step outright for about forty per cent of the
    spectra, and on one machine the worst refined error came out about
    half the worst unrefined one --- but which spectrum carries the
    worst error is decided by the last bit and is not a property to
    depend on.

    It costs one comparison per root and does not fire otherwise: over
    three thousand ordinary random Hermitian spectra it changed nothing,
    bit for bit, and on the stiff 3+1 spectrum the refinement exists for
    it leaves the refined error at :math:`5.5 \times 10^{-16}`.

    Parameters
    ----------
    traceless : numpy.ndarray
        Traceless part of the Hamiltonian, of shape ``(..., 4, 4)``.
    psi : numpy.ndarray
        Current estimates of the roots, of shape ``(..., 4)``.

    Returns
    -------
    numpy.ndarray
        The refined roots, sorted ascending along the last axis.
    """
    shifted = (psi[..., :, None, None]*np.eye(4)
               - traceless[..., None, :, :])
    chi = np.linalg.det(shifted).real

    gaps = psi[..., :, None] - psi[..., None, :]
    diagonal = np.arange(4)
    gaps[..., diagonal, diagonal] = 1.0
    derivative = np.prod(gaps, axis=-1)

    step = np.where(derivative != 0.0,
                    chi/np.where(derivative == 0.0, 1.0, derivative), 0.0)

    # The distance from each root to its nearest neighbour, which is what
    # the step is not allowed to cross.  The diagonal of `gaps` was set to
    # one for the product above and has to be excluded here, or a spectrum
    # whose roots are all further apart than one would be capped by it.
    nearest = np.min(np.where(np.eye(4, dtype=bool), np.inf, np.abs(gaps)),
                     axis=-1)
    step = np.where(np.abs(step) > 0.5*nearest, 0.0, step)

    return np.sort(psi - step, axis=-1)


def _latent_roots(traceless: np.ndarray) -> np.ndarray:
    r"""Returns the latent roots of a stack of traceless Hamiltonians.

    Parameters
    ----------
    traceless : numpy.ndarray
        Traceless Hamiltonians, of shape ``(..., 4, 4)``.

    Returns
    -------
    numpy.ndarray
        The four roots per element, of shape ``(..., 4)``.
    """
    squared = traceless @ traceless
    invariant_2 = 0.5*np.einsum('...ii->...', squared).real
    invariant_3 = 0.5*np.einsum('...ii->...', squared @ traceless).real
    invariant_4 = 0.5*(np.einsum('...ij,...ji->...', squared, squared).real
                       - invariant_2*invariant_2)

    quadratic = -invariant_2
    linear = -(2.0/3.0)*invariant_3
    constant = 0.25*(invariant_2*invariant_2 - 2.0*invariant_4)

    coeff_2 = 2.0*quadratic
    coeff_1 = quadratic*quadratic - 4.0*constant
    coeff_0 = -linear*linear

    depressed_p = coeff_1 - coeff_2*coeff_2/3.0
    depressed_q = 2.0*coeff_2**3/27.0 - coeff_2*coeff_1/3.0 + coeff_0
    shift = -coeff_2/3.0

    scale = 2.0*np.sqrt(np.maximum(-depressed_p, 0.0)/3.0)
    denominator = depressed_p*scale
    argument = np.clip(3.0*depressed_q/np.where(denominator != 0.0,
                                                denominator, 1.0), -1.0, 1.0)
    angle = np.arccos(argument)

    index = np.arange(3)
    resolvent = (scale[..., None]
                 * np.cos((angle[..., None] + 2.0*np.pi*index)/3.0)
                 + shift[..., None])
    root = np.sqrt(np.clip(resolvent, 0.0, None))
    root = np.where((linear > 0.0)[..., None] & (index == 2), -root, root)

    psi = 0.5*np.stack([
        root[..., 0] + root[..., 1] + root[..., 2],
        root[..., 0] - root[..., 1] - root[..., 2],
        -root[..., 0] + root[..., 1] - root[..., 2],
        -root[..., 0] - root[..., 1] + root[..., 2],
    ], axis=-1)
    psi = np.sort(psi, axis=-1)

    if POLISH_ROOTS:
        psi = _polish_roots(traceless, psi)

    return psi


def _traceless_part(
    hamiltonian_matrix: Union[list, np.ndarray]
) -> np.ndarray:
    r"""Returns the traceless part of a Hamiltonian or a stack of them.

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Hamiltonian, or stack of them, of shape ``(..., 4, 4)``.

    Returns
    -------
    numpy.ndarray
        The traceless part, of the same shape.
    """
    matrix = np.asarray(hamiltonian_matrix, dtype=complex)
    trace = np.einsum('...ii->...', matrix).real/4.0

    return matrix - trace[..., None, None]*np.eye(4)


def _divided_differences(
    psi: np.ndarray,
    baseline: np.ndarray
) -> np.ndarray:
    r"""Returns the divided differences of :math:`e^{-i\psi L}`.

    Builds the Newton coefficients :math:`c_j = f[\psi_1, \ldots,
    \psi_{j+1}]` of the exponential over the four latent roots, taking
    the confluent limit wherever two nodes are closer than
    :data:`DEGENERACY_TOL` relative to the spectral scale.  For the
    exponential the confluent value is exact rather than approximated:
    :math:`f[\psi, \ldots, \psi]` over :math:`k+1` equal nodes is
    :math:`(-iL)^k e^{-i\psi L}/k!`.

    Parameters
    ----------
    psi : numpy.ndarray
        Latent roots, sorted ascending, of shape ``(..., 4)``.
    baseline : numpy.ndarray
        Baselines, of shape ``(...)``.

    Returns
    -------
    numpy.ndarray
        The four coefficients, of shape ``(..., 4)``.
    """
    scale = np.max(np.abs(psi), axis=-1)
    tolerance = DEGENERACY_TOL*np.where(scale > 0.0, scale, 1.0)

    phase = np.exp(-1.j*psi*baseline[..., None])
    minus_i_l = -1.j*baseline

    table = phase
    coefficients = [table[..., 0]]

    for order in range(1, 4):
        separation = psi[..., order:] - psi[..., :-order]
        # The confluent value, used wherever the nodes have merged.
        confluent = (minus_i_l[..., None]**order/math.factorial(order)
                     * phase[..., :-order])
        regular = ((table[..., 1:] - table[..., :-1])
                   / np.where(np.abs(separation) > tolerance[..., None],
                              separation, 1.0))
        table = np.where(np.abs(separation) > tolerance[..., None],
                         regular, confluent)
        coefficients.append(table[..., 0])

    return np.stack(coefficients, axis=-1)


def _evolution_operator_4nu_array(
    hamiltonian_matrix: Union[list, np.ndarray],
    L: Union[int, float, list, np.ndarray]
) -> np.ndarray:
    r"""Returns :math:`U_4` for a stack, as an array.

    Interpolates :math:`e^{-i\psi L}` through the four latent roots in
    Newton form,

    .. math::
       U_4 = c_0 \mathbb{1} + c_1 (\tilde{H} - \psi_1)
             + c_2 (\tilde{H} - \psi_1)(\tilde{H} - \psi_2)
             + c_3 (\tilde{H} - \psi_1)(\tilde{H} - \psi_2)
                   (\tilde{H} - \psi_3) ,

    with the :math:`c_j` the divided differences of the exponential over
    the roots.  This is the Cayley-Hamilton form written in the basis
    that tolerates repeated nodes: a *confluent* divided difference is
    just a derivative, and for the exponential that derivative is known
    in closed form, :math:`f^{(k)}(\psi)/k! = (-iL)^k e^{-i\psi L}/k!`.

    The obvious alternative --- solving the Vandermonde system
    :math:`\sum_j c_j \psi_m^j = e^{-i\psi_m L}` --- was tried first and
    rejected: it is singular the moment two roots coincide, which is not
    an exotic case but includes a Hamiltonian proportional to the
    identity, a zero Hamiltonian, and any triply degenerate spectrum.
    The Newton form handles all of them without a special branch, and
    needs no eigenvectors either.

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Hamiltonian, or stack of them, of shape ``(..., 4, 4)``.
    L : int or float or array_like
        Baseline, or array of baselines.

    Returns
    -------
    numpy.ndarray
        The evolution operator, of shape ``(..., 4, 4)``.
    """
    if CHECK_HERMITICITY:
        _check_hermitian(np.asarray(hamiltonian_matrix, dtype=complex),
                         'oscprob4nu')

    traceless = _traceless_part(hamiltonian_matrix)
    baseline = np.asarray(L, dtype=float)

    leading = np.broadcast_shapes(traceless.shape[:-2], baseline.shape)
    traceless = np.broadcast_to(traceless, leading + (4, 4))
    baseline = np.broadcast_to(baseline, leading)

    # The compiled path, where there is one.  This is what `slabs` and
    # `earth` reach at four flavors: they compose operators, so they never
    # call `probabilities_4nu` and saw no kernel until this existed.  The
    # traceless part is passed on, since that is what the kernel expands
    # and what the NumPy path below uses.
    size = int(np.prod(leading, dtype=np.int64))
    if size > 0 and fastkernels.worthwhile(4, size):
        return fastkernels.evolution_operator_4nu_kernel(
            traceless, baseline, POLISH_ROOTS)

    psi = _latent_roots(traceless)
    coefficients = _divided_differences(psi, baseline)

    identity = np.eye(4)
    factor = np.broadcast_to(identity, leading + (4, 4)).copy()
    operator = coefficients[..., 0, None, None]*identity

    for order in range(1, 4):
        factor = factor @ (traceless
                           - psi[..., order-1, None, None]*identity)
        operator = operator + coefficients[..., order, None, None]*factor

    return operator


def evolution_operator_4nu_u_coefficients(
    hamiltonian_matrix: Union[list, np.ndarray],
    L: Union[int, float]
) -> List[complex]:
    r"""Returns the :math:`u_a` of the SU(4) expansion of :math:`U_4`.

    Returns the sixteen coefficients :math:`u_0, u_1, \ldots, u_{15}` of
    the expansion :math:`U_4(L) = u_0 \mathbb{1} + i u_a \lambda^a` of the
    time-evolution operator, in the same convention as
    :func:`oscprob3nu.evolution_operator_3nu_u_coefficients`.

    The four-flavor analogue of Eqs. (10)-(11) of arXiv:1904.12391 is

    .. math::
       u_0 = \frac14 \sum_m e^{-i \psi_m L} , \qquad
       i u_a = \sum_m e^{-i \psi_m L}\,
       \frac{\left(\psi_m^2 - \tfrac12 I_2\right) h_a
             + \psi_m (h \star h)_a + ((h \star h) \star h)_a}
            {\chi'(\psi_m)} ,

    with :math:`\chi'(\psi_m) = 4\psi_m^3 - 2 I_2 \psi_m - \tfrac23 I_3`.
    Note the third term in the numerator: at three flavors the tower
    closes on itself and no such term appears.

    .. versionadded:: 1.9.0

    .. versionchanged:: 1.10.1
       Results changed for a spectrum with two nearly coincident latent
       roots.  The Newton refinement of :data:`POLISH_ROOTS` divides by
       a product of gaps, and was guarded only against a gap of exactly
       zero; a pair separated by one unit in the last place passed that
       guard and was thrown across the spectrum.  The step is now
       refused whenever it would carry a root more than halfway to its
       nearest neighbour.  Nothing else moved: ordinary and stiff
       spectra are bit-for-bit what 1.10.0 gave.

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Four-flavor Hermitian Hamiltonian, given as a nested list or an
        array of shape ``(4, 4)``.
    L : int or float
        Baseline, in units reciprocal to those of the Hamiltonian.

    Returns
    -------
    list of complex
        The sixteen coefficients ``[u0, u1, ..., u15]``.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import oscprob4nu

        hamiltonian = np.diag([1.0, 2.0, 3.0, -6.0]).astype(complex)
        u = oscprob4nu.evolution_operator_4nu_u_coefficients(hamiltonian, 1.0)

        print('%d coefficients' % len(u))
        print('u0 = %.6f%+.6fj' % (u[0].real, u[0].imag))
    """
    operator = _evolution_operator_4nu_array(hamiltonian_matrix, L)

    u_0 = np.trace(operator)/4.0
    u_a = 0.5*np.einsum('aij,ji->a', LAMBDA_SU4, operator)/1.j

    return [complex(u_0)] + [complex(value) for value in u_a]


def evolution_operator_4nu(
    hamiltonian_matrix: Union[list, np.ndarray],
    L: Union[int, float, list, np.ndarray]
) -> Union[List[List[complex]], np.ndarray]:
    r"""Returns the four-neutrino time-evolution operator :math:`U_4`.

    Returns :math:`U_4(L) = e^{-i \tilde{H} L}`, with :math:`\tilde{H}`
    the traceless part of `hamiltonian_matrix`.  The discarded trace
    contributes only the overall phase :math:`e^{-i h_0 L}`, which
    cancels in the probabilities.

    .. versionadded:: 1.9.0

    .. versionchanged:: 1.10.1
       Results changed for a spectrum with two nearly coincident latent
       roots.  The Newton refinement of :data:`POLISH_ROOTS` divides by
       a product of gaps, and was guarded only against a gap of exactly
       zero; a pair separated by one unit in the last place passed that
       guard and was thrown across the spectrum.  The step is now
       refused whenever it would carry a root more than halfway to its
       nearest neighbour.  Nothing else moved: ordinary and stiff
       spectra are bit-for-bit what 1.10.0 gave.

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Four-flavor Hermitian Hamiltonian, given as a nested list of
        shape ``(4, 4)``, or a stack of them, of shape ``(..., 4, 4)``.
    L : int or float or array_like
        Baseline, in units reciprocal to those of the Hamiltonian, or an
        array of baselines broadcastable against the leading axes of
        `hamiltonian_matrix`.

    Returns
    -------
    list of list of complex or numpy.ndarray
        For a single Hamiltonian and baseline, the :math:`4\times4`
        evolution operator as a nested list.  If either argument is a
        stack, an array of shape ``(..., 4, 4)``.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import oscprob4nu

        hamiltonian = np.diag([1.0, 2.0, 3.0, -6.0]).astype(complex)
        operator = oscprob4nu.evolution_operator_4nu(hamiltonian, 1.0)

        unitarity = np.asarray(operator).conj().T @ np.asarray(operator)
        print('|U| unitary to %.1e' % np.max(np.abs(unitarity - np.eye(4))))
    """
    operator = _evolution_operator_4nu_array(hamiltonian_matrix, L)

    if operator.ndim == 2:
        return [[complex(entry) for entry in row] for row in operator]

    return operator


def probabilities_4nu(
    hamiltonian_matrix: Union[list, np.ndarray],
    L: Union[int, float, list, np.ndarray]
) -> Union[Tuple[float, ...], np.ndarray]:
    r"""Returns the four-neutrino oscillation probabilities.

    Returns the sixteen flavor-transition probabilities
    :math:`P_{\alpha\beta} \equiv P(\nu_\alpha \to \nu_\beta)
    = |[U_4]_{\beta\alpha}|^2`, ordered with the initial flavor varying
    slowest, exactly as at two and three flavors.  With the fourth state
    read as sterile, the flavor order is
    :math:`(\nu_e, \nu_\mu, \nu_\tau, \nu_s)`.

    .. versionadded:: 1.9.0

    .. versionchanged:: 1.10.1
       Results changed for a spectrum with two nearly coincident latent
       roots.  The Newton refinement of :data:`POLISH_ROOTS` divides by
       a product of gaps, and was guarded only against a gap of exactly
       zero; a pair separated by one unit in the last place passed that
       guard and was thrown across the spectrum.  The step is now
       refused whenever it would carry a root more than halfway to its
       nearest neighbour.  Nothing else moved: ordinary and stiff
       spectra are bit-for-bit what 1.10.0 gave.

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Four-flavor Hermitian Hamiltonian, given as a nested list of
        shape ``(4, 4)``, or a stack of them, of shape ``(..., 4, 4)``.
    L : int or float or array_like
        Baseline, in units reciprocal to those of the Hamiltonian, or an
        array of baselines broadcastable against the leading axes of
        `hamiltonian_matrix`.

    Returns
    -------
    tuple of float or numpy.ndarray
        For a single Hamiltonian and baseline, the sixteen
        probabilities as a tuple, ordered
        :math:`P_{ee}, P_{e\mu}, P_{e\tau}, P_{es}, P_{\mu e}, \ldots,
        P_{ss}`.  If either argument is a stack, an array of shape
        ``(..., 16)`` in the same order.

    Notes
    -----
    Passing arrays evaluates the whole stack at once, which is far
    faster than calling this routine in a Python loop, exactly as at two
    and three flavors; see :ref:`scanning`.

    The accuracy deserves a word, because four flavors is where it stops
    being automatic.  Against :func:`numpy.linalg.eigh`, a generic
    Hamiltonian agrees to about :math:`2 \times 10^{-13}`.  A stiff 3+1
    spectrum, with :math:`\Delta m^2_{41}` four orders of magnitude above
    :math:`\Delta m^2_{21}`, agrees to about :math:`10^{-9}` --- limited
    not by this expansion but by what double precision retains when the
    invariants are formed.  See :data:`POLISH_ROOTS`, which is what keeps
    that figure from being :math:`5 \times 10^{-7}`.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import oscprob4nu

        hamiltonian = np.diag([1.0, 2.0, 3.0, -6.0]).astype(complex)
        prob = oscprob4nu.probabilities_4nu(hamiltonian, 1.0)

        print('%d probabilities' % len(prob))
        print('P_ee = %.6f' % prob[0])
        print('they sum to %.6f' % sum(prob[0:4]))
    """
    matrix = np.asarray(hamiltonian_matrix, dtype=complex)
    baseline = np.asarray(L, dtype=float)

    # A single Hamiltonian and baseline must keep returning a tuple, and
    # is not worth the shape arithmetic either, so it leaves before it
    if matrix.ndim > 2 or baseline.ndim > 0:
        batch = np.broadcast_shapes(matrix.shape[:-2], baseline.shape)
        size = int(np.prod(batch, dtype=np.int64))

        if size > 0 and fastkernels.worthwhile(4, size):
            # The compiled kernel returns without going through
            # `_evolution_operator_4nu_array`, so it is checked here
            # instead --- and only here, so that no path checks twice
            if CHECK_HERMITICITY:
                _check_hermitian(matrix, 'probabilities_4nu')
            return fastkernels.probabilities_4nu_kernel(
                np.broadcast_to(matrix, batch+(4, 4)),
                np.broadcast_to(baseline, batch),
                POLISH_ROOTS)

    operator = _evolution_operator_4nu_array(matrix, baseline)

    # P[alpha][beta] = |U[beta][alpha]|^2, so the initial flavor varies
    # slowest once the last two axes are swapped and flattened.
    probabilities = np.abs(np.swapaxes(operator, -1, -2))**2

    if probabilities.ndim == 2:
        return tuple(float(value) for value in probabilities.reshape(16))

    return probabilities.reshape(probabilities.shape[:-2] + (16,))
