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

  **Sign convention, and a trap.**  Here the :math:`\psi_m` are the
  eigenvalues of :math:`\tilde{H}` itself.  :func:`oscprob3nu.psi_roots`
  returns those of :math:`-\tilde{H}`, the opposite sign.  Both modules
  are internally consistent and the choice is invisible from outside
  either, but the quartic above depends on it: under
  :math:`\psi \to -\psi` the linear term is the only one that flips.
  The expression is therefore right as written and wrong if carried into
  the other convention -- which has twice produced an error in the
  accompanying paper, precisely because it looks correct on its own.

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

So the invariants are formed in *double-double* arithmetic instead ---
each number a pair of ``float64`` limbs, together carrying some 32
digits --- which leaves the cluster's separation intact where double
precision destroyed it.  One Aberth sweep then takes the quartic's roots
to the last ``float64`` bit: :math:`3.6 \times 10^{-17}` worst over the
nine reference Hamiltonians, under a fifth of an ulp.  :data:`ROOT_STRATEGY`
has the measured comparison, and switches to the older route, which
escaped the same bottleneck by refining against the *matrix* with one
Newton step on :math:`\chi(\psi) = \det(\psi \mathbb{1} - \tilde{H})` and
reached :math:`3.9 \times 10^{-16}`.

Neither figure is anywhere near a measurable effect --- probabilities
are confronted with data at the per-cent level --- so this matters for
the exactness claim, for error accumulating when :mod:`slabs` and
:mod:`earth` compose operators across layers, and for a regression suite
tight enough to catch a real mistake.  :data:`POLISH_ROOTS` records why
LAPACK and extended precision both lose, the latter being the reason the
extra precision here is built out of pairs of ``float64``.

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

__all__ = ['CHECK_HERMITICITY', 'LAMBDA_SU4', 'POLISH_ROOTS',
           'ROOT_STRATEGY', 'SMALL_BATCH',
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

ROOT_STRATEGY = 'double-double'
r"""str: Module-level switch.

How the four latent roots :math:`\psi_m` are obtained.  Either
``'double-double'``, the default, or ``'eigensolver'``.

``'double-double'`` forms :math:`I_2, I_3, I_4` in double-double
arithmetic --- each number a pair of ``float64`` limbs, carrying some 32
digits --- and refines the quartic's roots with one Aberth sweep, also in
double-double.  The invariants are the problem this solves: they compress
the matrix into three numbers, and the amplification from those
coefficients to the roots was measured at :math:`2.3 \times 10^9`, so a
:math:`10^{-16}` coefficient becomes a :math:`10^{-7}` root.  At
:math:`10^{-32}` the same amplification lands at :math:`10^{-23}`, well
under ``float64``, and the answer is limited by the final rounding rather
than by the algebra.

``'eigensolver'`` calls :func:`numpy.linalg.eigvalsh` on the traceless
matrix, which forms no invariants at all, and then takes the Newton step
:data:`POLISH_ROOTS` asks for.  It is the strategy the library shipped in
1.12.0, kept because it depends on nothing but LAPACK, and because a
second route that agrees to one ulp is worth having when a result is in
doubt.

Root errors are against ``mpmath`` at fifty decimal digits over
``tests/stiff_reference.json``, worst of its nine cases, which run from
the physical 3+1 range to :math:`\Delta m^2_{41} = 1000\ \mathrm{eV}^2`
and down to a pair separated by :math:`10^{-16}`.  Cost is a full
:func:`probabilities_4nu` over a 100 000-point stiff 3+1 stack through the
compiled kernel, relative to the route 1.12.0 shipped:

========================================  ==========  =========
Strategy for the latent roots             Rel. error  Cost
========================================  ==========  =========
``'double-double'``                       3.6e-17     1.15-1.3x
``'eigensolver'``, with the Newton step   3.9e-16     1.00x
``'eigensolver'``, without it             6.9e-16     0.9-1.0x
Closed form alone (`psi_roots_4nu`)       2.2e-07     ---
========================================  ==========  =========

`psi_roots_4nu` still solves the quartic in closed form and its contract is
unchanged, but nothing in the probability path is built on it any more, so
there is no comparable figure to time.  Its error is quoted on the same
nine Hamiltonians as the rest, which is what it is there to show.

Note also what the third row costs against the second: the Newton step is
*within the noise*, under a tenth, not the fifth an earlier single-shot
measurement here suggested.

The errors are exact figures, reproducible to the digits shown.  The costs
are not, and are given as ranges on purpose: timed in alternated pairs to
cancel machine drift, the double-double row lands at 1.18x, 1.20x and
1.25x by three different methods, while the *per-pair* spread runs from
0.76x to 1.72x.  Anything quoted more precisely than this would be
describing one afternoon's load average.  The last row's own cost is
absent because it is no longer a route a caller can select --- see below.

So the default buys an order of magnitude on the roots for roughly a fifth
more time.  Two things it is worth being clear about.

The cost is not double-double's alone.  Both eigensolver rows call
:func:`numpy.linalg.eigvalsh` once per element *inside* the compiled
kernel, across every core, which is why the middle row costs about what a
single batched :func:`numpy.linalg.eigvalsh` over the same stack costs
while also building sixteen probabilities from it.  Double-double adds its
invariants and one Aberth sweep on top of that same eigensolver call,
which it needs for the start --- see `_latent_roots_dd` for why the start
cannot instead be the closed form, however much faster that is.

On the NumPy path the ratio is worse, of order 1.5-2x rather than 1.2x,
the double-double primitives being elementwise NumPy operations with
nothing to amortise them over.  That path is the fallback for an installation
without Numba, which since 1.13.0 is not the usual one; accuracy is
deliberately the same on both, so that a result never depends on whether a
compiler was present.

.. versionadded:: 1.13.0
"""

POLISH_ROOTS = True
r"""bool: Module-level switch.

Whether to refine the latent roots against the Hamiltonian matrix, by
one Newton step on :math:`\chi(\psi) = \det(\psi \mathbb{1} - \tilde{H})`.

Read only when :data:`ROOT_STRATEGY` is ``'eigensolver'``, and ignored
entirely by the default double-double route, which reaches the
``float64`` floor on its own and leaves a Newton step nothing to find.
On the eigensolver route it is on by default and should stay on: it is
worth a factor of about two on the roots, 3.9e-16 against 6.9e-16 over
``tests/stiff_reference.json``, for a cost too small to measure reliably
here --- under a tenth, and inside the run-to-run spread.

The reason it helps at all is that ``eigvalsh`` reduces by similarity
transforms which each carry a backward error of order
:math:`\epsilon\|H\|`, while this converges onto the root of
:math:`\det(\psi\mathbb{1} - \tilde{H})` for the matrix it was actually
given.

Its history is worth keeping, because it is what the default route was
built to finish.  This switch was introduced when the roots came from
Euler's reduction of the quartic, where the invariants
:math:`I_2, I_3, I_4` --- not the quartic solver --- set the floor:
perturbing those three numbers at the :math:`10^{-16}` level moves a
stiff 3+1 spectrum's roots by :math:`10^{-7}` relative, an amplification
of :math:`2.3 \times 10^9` that no better root-finder can beat.  One
Newton step against the matrix escaped that floor, from 2.2e-07 to
3.9e-16, and for 1.12.0 that was the answer.  It is not the only one:
:data:`ROOT_STRATEGY`'s default keeps the invariants and removes the
floor instead, by forming them in double-double, where the same
amplification lands at :math:`10^{-23}`.

Extended precision was tried on the way and rejected --- ``longdouble``
bought under a digit, was slower, and is silently ``float64`` on Apple
Silicon and Windows, which is precisely why the default route implements
its own extra precision out of pairs of ``float64`` instead.

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

A second Newton step changes nothing --- one already reaches this route's
floor --- so exactly one is taken.

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

That conclusion survived being retested for the default route, and turns
out to cut the other way there: because Aberth's step *is* an evaluation
of :math:`\chi` at four roots, its own magnitude certifies the answer for
nothing --- about :math:`10^{-32}` relative once converged against
:math:`10^{-9}` while still crawling up a cluster.  A cheaper start
gated on that certificate was built and measured, and dropped for being
slower in the end rather than for being unsound; see `_latent_roots_dd`.

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

    non_finite_message = (
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
                    raise ValueError(
                        non_finite_message % (caller, 'oscprob4nu'))
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

    # The compiled scan, where there is one.  This check was cheap against
    # the expansion it guards and is no longer: once the kernels made the
    # expansion some ten times quicker, the NumPy comparisons below came to
    # most of a whole probability call on a large stack.  The compiled scan
    # reaches the same verdict, tolerance and all, in one pass per stage
    # instead of a dozen reductions that each allocate a temporary the size
    # of the stack.  See `fastkernels.MIN_BATCH_HERMITICITY`.
    #
    # The flavor count is read off the array rather than written as a
    # literal, so that this block is textually identical in all three
    # modules --- which `tests/test_edge_cases.py` requires of these copies.
    n_flavors = h_matrix.shape[-1]
    if (fastkernels.available()
            and h_matrix.size//(n_flavors*n_flavors)
            >= fastkernels.MIN_BATCH_HERMITICITY):
        non_finite, element, i, j = fastkernels.hermiticity_offender(
            h_matrix, n_flavors, _HERMITICITY_TOL)
        if non_finite:
            raise ValueError(
                non_finite_message % (caller, 'oscprob4nu'))
        if element >= 0:
            if i == j:
                raise ValueError(complaint % (
                    caller, ' --- the diagonal entry (%d, %d) has a non-zero '
                    'imaginary part' % (i, i)))
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
        raise ValueError(non_finite_message % (caller, 'oscprob4nu'))

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

    What feeds this changed in 1.13.0 and the argument below did not.  It
    is reached only when :data:`ROOT_STRATEGY` is ``'eigensolver'``, so
    `psi` now arrives from :func:`numpy.linalg.eigvalsh` rather than from
    Euler's reduction of the quartic.  Everything about the guard carries
    over unaltered --- it is a statement about :math:`\chi'` being a
    product of gaps, and about what a step means when two of them nearly
    coincide, neither of which cares where the estimates came from.  The
    :math:`\sqrt{\epsilon}` discussion further down is the exception: that
    one is specifically about Euler's reduction, and it is kept because it
    is the clearest statement of the limit a Newton step cannot cross, and
    therefore of why :data:`ROOT_STRATEGY`'s default gives the invariants
    more precision instead of trying to repair them afterwards.

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


DD_SPLIT = 134217729.0
r"""float: Module-level constant.

:math:`2^{27} + 1`, Dekker's splitting constant.  Multiplying by it and
subtracting isolates the leading 26 bits of a ``float64``, narrow enough
that the product of two halves is exact and `_two_prod` can hand back the
rounding error it would otherwise have discarded.
"""

DD_SWEEPS = 1
r"""int: Module-level constant.

Aberth sweeps taken over the four latent roots by the double-double route.

One, measured: every case in ``tests/stiff_reference.json`` gives an
identical answer at one, two and three sweeps, to the last bit.  A single
sweep converges because the start is already good --- the eigensolver's
quartet with the residual-trace shift removed *in double-double*, so the
low limb the exact traceless-ing computed is not thrown away before the
iteration begins.

This constant was two for exactly as long as that subtraction was done in
``float64``.  Discarding the low limb put the start a whole ulp out and
cost a second sweep to recover, and the second sweep was then mistaken for
something the iteration needed --- the error even looked like evidence for
it, since one sweep did measure :math:`3.9 \times 10^{-16}`.  It is worth
knowing which way that dependency runs: this is not a tolerance to raise
when an answer looks inaccurate, it is a count that a correct start makes
sufficient.

`fastkernels` carries the same constant for its compiled kernel, as it
does :data:`DEGENERACY_TOL`, because a module imported by this one cannot
import back.
"""


def _dd_reciprocal_of_three() -> tuple:
    r"""Returns 1/3 as a double-double pair, exact to 3e-33.

    Formed once at import, by the same compensated arithmetic the route uses:
    :math:`3 h` is split exactly into ``p + e`` by `_two_prod`, so what 1/3
    is missing after its first limb is ``((1 - p) - e)/3``.

    It exists because :math:`2/3` is the one constant the quartic needs that
    binary cannot hold, and lifting a *rounded* 2/3 into double-double would
    zero the low limb and drag :math:`c_1` back to ``float64`` accuracy ---
    the :math:`1.1 \times 10^{-7}` stall this route exists to remove, wearing
    the disguise of a solver failure.  Multiplying by this pair costs no
    division, where `_dd_div` costs three.

    Returns
    -------
    tuple
        The high and low limbs, as Python floats.
    """
    hi = 1.0/3.0
    product, error = _two_prod(3.0, hi)
    return hi, ((1.0 - product) - error)/3.0


def _strategy_code() -> int:
    r"""Returns the two root switches as the one integer a kernel can take.

    ``0`` for double-double, ``1`` for the eigensolver with the Newton step,
    ``2`` for the eigensolver alone.  A compiled kernel cannot read a Python
    global at call time without recompiling, so the choice has to travel as
    an argument; collapsing two switches into one keeps that argument, and
    the number of kernel specialisations, to one.

    Returns
    -------
    int
        The code `fastkernels` expects.

    Raises
    ------
    ValueError
        If :data:`ROOT_STRATEGY` names neither route.
    """
    if ROOT_STRATEGY == 'double-double':
        return 0
    if ROOT_STRATEGY == 'eigensolver':
        return 1 if POLISH_ROOTS else 2
    raise ValueError(
        "ROOT_STRATEGY must be 'double-double' or 'eigensolver', not %r"
        % (ROOT_STRATEGY,))


def _two_sum(a: np.ndarray, b: np.ndarray) -> tuple:
    r"""Returns ``a + b`` rounded, and the error the rounding dropped.

    Knuth's compensated addition.  The second element is exactly
    :math:`a + b - s`, so the pair together represent the sum with nothing
    lost.  Elementwise, so every primitive here vectorises over a stack.
    """
    s = a + b
    bb = s - a
    return s, (a - (s - bb)) + (b - bb)


def _quick_two_sum(a: np.ndarray, b: np.ndarray) -> tuple:
    r"""Returns ``a + b`` and its error, given :math:`|a| \geq |b|`."""
    s = a + b
    return s, b - (s - a)


def _two_prod(a: np.ndarray, b: np.ndarray) -> tuple:
    r"""Returns ``a*b`` rounded, and the error the rounding dropped.

    Dekker's algorithm: split both operands at :data:`DD_SPLIT` into halves
    narrow enough to multiply exactly, then assemble the four partial
    products and subtract the rounded one.
    """
    p = a*b
    t = DD_SPLIT*a
    a_hi = t - (t - a)
    a_lo = a - a_hi
    t = DD_SPLIT*b
    b_hi = t - (t - b)
    b_lo = b - b_hi
    return p, ((a_hi*b_hi - p) + a_hi*b_lo + a_lo*b_hi) + a_lo*b_lo


def _dd_add(x: tuple, y: tuple) -> tuple:
    r"""Returns the double-double sum of two ``(hi, lo)`` pairs."""
    s, e = _two_sum(x[0], y[0])
    return _quick_two_sum(s, e + x[1] + y[1])


def _dd_sub(x: tuple, y: tuple) -> tuple:
    r"""Returns the double-double difference, ``x - y``."""
    return _dd_add(x, (-y[0], -y[1]))


def _dd_mul(x: tuple, y: tuple) -> tuple:
    r"""Returns the double-double product, ``x*y``."""
    p, e = _two_prod(x[0], y[0])
    return _quick_two_sum(p, e + x[0]*y[1] + x[1]*y[0])


def _dd_div(x: tuple, y: tuple) -> tuple:
    r"""Returns the double-double quotient, ``x/y``.

    Three ``float64`` divisions, each correcting the remainder the one
    before left.  Much the most expensive primitive here, which is why the
    iteration below is Aberth's with its divisions counted rather than
    anything more elaborate.
    """
    q1 = x[0]/y[0]
    r = _dd_sub(x, _dd_mul(y, (q1, np.zeros_like(q1))))
    q2 = r[0]/y[0]
    r = _dd_sub(r, _dd_mul(y, (q2, np.zeros_like(q2))))
    q3 = r[0]/y[0]
    s, e = _quick_two_sum(q1, q2)
    return _dd_add((s, e), (q3, np.zeros_like(q3)))


def _dd_scale(x: tuple, factor: float) -> tuple:
    r"""Returns ``factor*x`` for a `factor` that is a power of two.

    Scaling both limbs by a power of two only shifts exponents, so this is
    exact and costs two flops where `_dd_mul` costs some twenty.  Every
    factor of 2, 4, 1/2 and 1/4 in the invariants and in Aberth's step goes
    through here; do not pass it anything else.
    """
    return (factor*x[0], factor*x[1])


THIRD_HI, THIRD_LO = _dd_reciprocal_of_three()
r"""tuple: Module-level constant.

The two ``float64`` limbs of 1/3, from `_dd_reciprocal_of_three`.  Bound
here rather than beside :data:`DD_SPLIT` because forming it needs
`_two_prod`, which is defined below it.
"""


def _latent_roots_dd(traceless: np.ndarray) -> np.ndarray:
    r"""Returns the latent roots to the last ``float64`` bit.

    Three invariants compress a :math:`4 \times 4` matrix into three
    numbers, and in ``float64`` that compression is where a stiff spectrum
    loses what separates its clustered roots: the amplification from
    coefficients to roots was measured at :math:`2.3 \times 10^9`, so a
    :math:`10^{-16}` coefficient becomes a :math:`10^{-7}` root.  Forming
    the invariants in double-double arithmetic --- each number a pair of
    ``float64`` limbs, carrying some 32 digits --- leaves a
    :math:`10^{-32}` coefficient error to amplify to :math:`10^{-23}`, and
    the roots are then limited by the final rounding to ``float64`` instead
    of by the algebra: :math:`3.6 \times 10^{-17}` worst over
    ``tests/stiff_reference.json``, against :math:`3.9 \times 10^{-16}` for
    the eigensolver with a Newton step and :math:`2.2 \times 10^{-7}` for
    the closed form the invariants used to be handed to.

    Parameters
    ----------
    traceless : numpy.ndarray
        Traceless Hamiltonians, of shape ``(..., 4, 4)``.

    Returns
    -------
    numpy.ndarray
        The four roots per element, ascending, of shape ``(..., 4)``.

    Notes
    -----
    Four choices here were each made against measurement, and each would be
    easy to undo by accident.

    `traceless` arrives with its trace removed in ``float64``, which leaves
    a residue :math:`\tau` of order :math:`10^{-23}`.  Removing that
    residue exactly in double-double is not a refinement but the difference
    between a quartic that describes this matrix and one that does not:
    otherwise the invariants belong to a matrix whose trace is not quite
    zero while the quartic pins its cubic coefficient to exactly zero, and
    the two disagree by up to :math:`3.8 \times 10^{-7}`.

    :math:`\tilde{H}` is Hermitized exactly, as
    :math:`(\tilde{H} + \tilde{H}^\dagger)/2`, which is representable
    because `_two_sum` of two ``float64`` loses nothing and halving only
    shifts exponents.  That makes :math:`\tilde{H}^2` exactly Hermitian, so
    only its upper triangle is formed --- ten entries rather than sixteen,
    a quarter off the work.  Mirroring *without* the Hermitization looks
    equivalent and costs :math:`I_3` and :math:`I_4` about
    :math:`10^{-17}`, because a Hamiltonian built as
    :math:`U M^2 U^\dagger` is Hermitian only to rounding and the mirror
    discards the asymmetry silently.  Discarding it deliberately is free:
    an anti-Hermitian perturbation moves a *real* eigenvalue only at second
    order, since :math:`\langle v|\delta|v\rangle` is imaginary.

    :math:`\mathrm{Tr}(\tilde{H}^2)` is taken off :math:`\tilde{H}` itself
    as :math:`\sum_{ij}|\tilde{H}_{ij}|^2`, and
    :math:`\mathrm{Tr}(\tilde{H}^4)` off :math:`\tilde{H}^2` as
    :math:`\sum_{ij}|S_{ij}|^2`.  Both are sums of squares, with no
    cancellation to lose digits to.

    The start is the eigensolver's, because what a start must supply is not
    proximity but four distinct basins.  Euler's closed form is twice as
    fast and exact on every stiff case, and still wrong here: on a cluster
    Aberth converges *linearly*, at ratio one half, so from :math:`10^{-7}`
    five sweeps reach :math:`3.8 \times 10^{-9}` and it would need some
    thirty.  A backward-stable Hermitian eigensolver separates a cluster as
    well as ``float64`` allows, which together with removing the
    residual-trace shift in double-double is what makes :data:`DD_SWEEPS`
    one rather than many.  Durand-Kerner was measured as well, at one division
    per root against Aberth's five, and rejected for being non-monotone in
    the sweep count: :math:`3.9 \times 10^{-16}`, then
    :math:`9.7 \times 10^{-17}`, then :math:`1.9 \times 10^{-16}`.
    """
    traceless = np.asarray(traceless, dtype=complex)
    zero = (np.zeros(traceless.shape[:-2]), np.zeros(traceless.shape[:-2]))

    def lift(a):
        # np.copy, not np.ascontiguousarray: the latter promotes the 0-d
        # array a single Hamiltonian's entry slices down to into shape (1,),
        # which would put a phantom leading axis on the returned roots and
        # on everything reconstructed from them.  Both give a contiguous
        # copy of the strided slice a stack produces.
        return (np.copy(a), np.zeros_like(a))

    def const(value):
        return (np.full_like(zero[0], value), zero[1])

    # Exact Hermitization.  A Hermitian diagonal is real, so its imaginary
    # limbs are exactly zero and never enter a product below.
    real = [[None]*4 for _ in range(4)]
    imag = [[None]*4 for _ in range(4)]
    for i in range(4):
        real[i][i] = lift(traceless[..., i, i].real)
        imag[i][i] = zero
        for j in range(i+1, 4):
            s, e = _two_sum(np.copy(traceless[..., i, j].real),
                            np.copy(traceless[..., j, i].real))
            real[i][j] = (0.5*s, 0.5*e)
            real[j][i] = real[i][j]
            s, e = _two_sum(np.copy(traceless[..., i, j].imag),
                            -np.copy(traceless[..., j, i].imag))
            imag[i][j] = (0.5*s, 0.5*e)
            imag[j][i] = (-imag[i][j][0], -imag[i][j][1])

    # The residual trace, removed exactly.  Dividing by four is exact.
    tau = real[0][0]
    for i in range(1, 4):
        tau = _dd_add(tau, real[i][i])
    shift = _dd_scale(tau, 0.25)
    for i in range(4):
        real[i][i] = _dd_sub(real[i][i], shift)

    trace_2 = zero
    for i in range(4):
        trace_2 = _dd_add(trace_2, _dd_mul(real[i][i], real[i][i]))
        for j in range(i+1, 4):
            term = _dd_add(_dd_mul(real[i][j], real[i][j]),
                           _dd_mul(imag[i][j], imag[i][j]))
            trace_2 = _dd_add(trace_2, _dd_scale(term, 2.0))

    square_re = [[zero]*4 for _ in range(4)]
    square_im = [[zero]*4 for _ in range(4)]
    for i in range(4):
        for j in range(i, 4):
            accumulated_re, accumulated_im = zero, zero
            for k in range(4):
                accumulated_re = _dd_add(
                    accumulated_re,
                    _dd_sub(_dd_mul(real[i][k], real[k][j]),
                            _dd_mul(imag[i][k], imag[k][j])))
                accumulated_im = _dd_add(
                    accumulated_im,
                    _dd_add(_dd_mul(real[i][k], imag[k][j]),
                            _dd_mul(imag[i][k], real[k][j])))
            square_re[i][j] = accumulated_re
            square_im[i][j] = accumulated_im

    # Tr(H~^3) = sum_i S_ii H~_ii + 2 sum_{i<j} Re(S_ij conj(H~_ij))
    # Tr(H~^4) = sum_i S_ii^2   + 2 sum_{i<j} |S_ij|^2
    trace_3 = zero
    trace_4 = zero
    for i in range(4):
        trace_3 = _dd_add(trace_3, _dd_mul(square_re[i][i], real[i][i]))
        trace_4 = _dd_add(trace_4, _dd_mul(square_re[i][i], square_re[i][i]))
        for j in range(i+1, 4):
            term = _dd_add(_dd_mul(square_re[i][j], real[i][j]),
                           _dd_mul(square_im[i][j], imag[i][j]))
            trace_3 = _dd_add(trace_3, _dd_scale(term, 2.0))
            term = _dd_add(_dd_mul(square_re[i][j], square_re[i][j]),
                           _dd_mul(square_im[i][j], square_im[i][j]))
            trace_4 = _dd_add(trace_4, _dd_scale(term, 2.0))

    invariant_2 = _dd_scale(trace_2, 0.5)
    invariant_3 = _dd_scale(trace_3, 0.5)
    invariant_4 = _dd_scale(_dd_sub(trace_4, _dd_mul(invariant_2,
                                                     invariant_2)), 0.5)

    # chi(psi) = psi^4 + c2 psi^2 + c1 psi + c0.  The two thirds is formed
    # by division: 2/3 is not representable in float64, so lifting a rounded
    # 2/3 into double-double zeroes the low limb and poisons c1 back to
    # float64 accuracy --- which is precisely the 1.1e-07 stall this route
    # exists to remove, and it looks like a solver failure rather than a
    # constant.
    c2 = (-invariant_2[0], -invariant_2[1])
    c1 = _dd_mul(_dd_scale(invariant_3, 2.0),
                 (const(THIRD_HI)[0], const(THIRD_LO)[0]))
    c1 = (-c1[0], -c1[1])
    c0 = _dd_scale(_dd_sub(_dd_mul(invariant_2, invariant_2),
                           _dd_scale(invariant_4, 2.0)), 0.25)

    start = np.linalg.eigvalsh(traceless)
    roots = [_dd_sub(lift(start[..., k]), shift) for k in range(4)]

    for _ in range(DD_SWEEPS):
        for i in range(4):
            squared = _dd_mul(roots[i], roots[i])
            chi = _dd_add(_dd_add(_dd_mul(squared, squared),
                                  _dd_mul(c2, squared)),
                          _dd_add(_dd_mul(c1, roots[i]), c0))

            # chi'(psi) = 4 psi^3 + 2 c2 psi + c1
            derivative = _dd_add(
                _dd_scale(_dd_mul(squared, roots[i]), 4.0),
                _dd_add(_dd_scale(_dd_mul(c2, roots[i]), 2.0), c1))

            # Aberth's correction is Newton's step divided by one minus the
            # pull of the other three roots, and that pull is what keeps a
            # cluster's members from converging onto each other --- it is
            # why no half-gap guard is needed on this route.  Written out,
            #
            #     step = (chi/chi') / (1 - (chi/chi') sum_j 1/(psi_i-psi_j))
            #
            # takes five double-double divisions per root: one for the
            # ratio, three for the reciprocal gaps, one for the step.
            # Clearing the denominators collapses all five into one,
            #
            #     step = chi P / (chi' P - chi S)
            #
            # for P and S the product and the second elementary symmetric
            # function of the three gaps, which cost three multiplications.
            # The two forms are the same expression; a dd division is three
            # float64 divisions where a dd multiplication is none, and this
            # loop was over half the route's cost.
            gaps = [_dd_sub(roots[i], roots[j])
                    for j in range(4) if j != i]
            pair = _dd_mul(gaps[1], gaps[2])
            product = _dd_mul(gaps[0], pair)
            symmetric = _dd_add(pair, _dd_mul(gaps[0],
                                              _dd_add(gaps[1], gaps[2])))

            denominator = _dd_sub(_dd_mul(derivative, product),
                                  _dd_mul(chi, symmetric))
            standing = denominator[0] != 0.0
            denominator = (np.where(standing, denominator[0], 1.0),
                           np.where(standing, denominator[1], 0.0))
            step = _dd_div(_dd_mul(chi, product), denominator)

            # A vanishing gap now leaves the step at exactly zero through
            # P rather than needing to be skipped: two roots that have
            # landed on identical bits stay where the eigensolver put them,
            # which is the same refusal the guarded Newton step makes.
            roots[i] = (np.where(standing, roots[i][0] - step[0],
                                 roots[i][0]),
                        np.where(standing, roots[i][1] - step[1],
                                 roots[i][1]))

    # Only here do the two limbs collapse into one float64
    psi = np.stack([roots[k][0] + roots[k][1] for k in range(4)], axis=-1)
    return np.sort(psi, axis=-1)


def _latent_roots(traceless: np.ndarray) -> np.ndarray:
    r"""Returns the latent roots of a stack of traceless Hamiltonians.

    Which route is taken is :data:`ROOT_STRATEGY`'s to decide, and
    :data:`POLISH_ROOTS`' only on the eigensolver route.

    Parameters
    ----------
    traceless : numpy.ndarray
        Traceless Hamiltonians, of shape ``(..., 4, 4)``.

    Returns
    -------
    numpy.ndarray
        The four roots per element, of shape ``(..., 4)``.

    Raises
    ------
    ValueError
        If :data:`ROOT_STRATEGY` names neither route.
    """
    if ROOT_STRATEGY == 'double-double':
        return _latent_roots_dd(traceless)
    if ROOT_STRATEGY != 'eigensolver':
        raise ValueError(
            "ROOT_STRATEGY must be 'double-double' or 'eigensolver', not %r"
            % (ROOT_STRATEGY,))

    # Eigenvalues of the matrix, not roots of its invariants.  Both are the
    # latent roots in exact arithmetic; they are not the same computation in
    # floating point.  `su4_invariants` compresses the matrix into three
    # numbers, and what separates a cluster of eigenvalues enters them at
    # relative order (gap/scale)^2 --- at a gap of 1e-6 that is 1e-12, four
    # digits above the noise.  A Hermitian eigensolver never forms them,
    # which is why this route beats the closed form without needing dd; it
    # stops one ulp short of `_latent_roots_dd`, which removes the
    # compression instead of avoiding it.
    #
    # `psi_roots_4nu` keeps the closed form, and its contract, for callers
    # who want the quartic solved.
    psi = np.linalg.eigvalsh(traceless)

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

    # cos and sin into the two halves, rather than `np.exp` of an
    # imaginary argument: identical to the bit, and about a third cheaper,
    # the complex exponential paying for an exp(0) it does not need.  See
    # `oscprob3nu._u_coefficients_3nu_batch`, where the same applies.
    angle = -psi*baseline[..., None]
    phase = np.empty(angle.shape, dtype=complex)
    np.cos(angle, out=phase.real)
    np.sin(angle, out=phase.imag)
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
    L: Union[int, float, list, np.ndarray],
    caller: str = 'evolution_operator_4nu'
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
    caller : str, optional
        Name of the public routine to blame in a Hermiticity complaint,
        since three of them reach this helper and none of them is called
        what this helper is called.

    Returns
    -------
    numpy.ndarray
        The evolution operator, of shape ``(..., 4, 4)``.
    """
    # Named by the caller rather than by this helper.  Three public
    # routines reach it, and one of them --- `probabilities_4nu` --- runs
    # the check itself before dispatching to the kernel, so labelling the
    # message from here made one user error report `probabilities_4nu` with
    # Numba installed and `oscprob4nu` without it.  The two- and
    # three-flavor modules name the public routine at their equivalent
    # sites; this is how a shared helper does the same.
    if CHECK_HERMITICITY:
        _check_hermitian(np.asarray(hamiltonian_matrix, dtype=complex),
                         caller)

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
            traceless, baseline, _strategy_code())

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
    operator = _evolution_operator_4nu_array(
        hamiltonian_matrix, L, 'evolution_operator_4nu_u_coefficients')

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
                _strategy_code())

    operator = _evolution_operator_4nu_array(matrix, baseline,
                                             'probabilities_4nu')

    # P[alpha][beta] = |U[beta][alpha]|^2, so the initial flavor varies
    # slowest once the last two axes are swapped and flattened.
    probabilities = np.abs(np.swapaxes(operator, -1, -2))**2

    if probabilities.ndim == 2:
        return tuple(float(value) for value in probabilities.reshape(16))

    return probabilities.reshape(probabilities.shape[:-2] + (16,))
