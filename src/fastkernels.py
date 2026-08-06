# -*- coding: utf-8 -*-
r"""Numba-compiled kernels for the batched evaluation paths.

`Numba <https://numba.pydata.org>`_ is a dependency of
**NuOscProbExact** as of 1.13.0, so this module compiles the two-, three-
and four-neutrino expansions into fused machine-code loops and
:mod:`oscprob2nu`, :mod:`oscprob3nu` and :mod:`oscprob4nu` use them for
large stacks.

It is still written to do without.  If the import fails --- an
environment where the dependency was removed on purpose, a platform with
no wheel --- ``HAVE_NUMBA`` is ``False``, nothing here is defined, and
the NumPy path is used instead.  Nothing else in the library changes
either way, and the results agree to round-off --- see
``tests/test_fastkernels.py``, which runs both paths against each other
whichever is available.  That fallback is not vestigial: it is the
independent implementation these kernels are checked against, and
``.github/workflows/tests.yml`` uninstalls Numba in one job to keep it
honest.

Why it is worth compiling
-------------------------

The NumPy path evaluates the expansion as a sequence of whole-array
operations, so a stack of N Hamiltonians makes roughly fifteen passes
over N-element arrays, each writing a temporary that the next pass
reads back.  The compiled kernel does the same arithmetic one element at
a time, keeping every intermediate in registers, and spreads the
elements over the available cores.  Measured against the NumPy path on
this library's own benchmarks, best of seven runs with the two paths
interleaved:

===============================  ==========
Stack                            Speedup
===============================  ==========
200 000 energies, four flavors   ~19x
20 000 energies, four flavors    ~18x
200 000 energies, three flavors  ~15x
20 000 energies, three flavors   ~9x
100 x 100 oscillogram            ~3.5x
200 000 baselines, two flavors   ~1.5x
===============================  ==========

Those are the *probability* kernels.  Composing operators across slabs
is a separate route with its own economics, and until 1.12.0 it had no
compiled path at all --- the backend offered probability kernels only,
which :mod:`slabs` and :mod:`earth` cannot use, so an Earth crossing ran
the NumPy path however this module was configured.  Against that path,
best of five rounds with the two interleaved:

===============================  ==========  ==========
Slab sequence                    1 slab      256 slabs
===============================  ==========  ==========
Two flavors                      ~137x       ~59x
Three flavors                    ~225x       ~13x
Four flavors                     ~187x       ~7x
===============================  ==========  ==========

The margin narrows with length, which is the reverse of the table above,
because what is being avoided here is the *caller's* fixed cost --- a
dispatched matrix product per slab --- rather than the kernel's.  A whole
Earth crossing at 120 slabs comes out ~12x, ~12x and ~9x quicker at two,
three and four flavors.

Since 1.12.0 that is not the whole of it, and for a scan it is not even
the larger part.  An Earth scan takes its energies as an array and its
zenith angles as another, so the geometry and the matter potentials are
built once rather than once per point; the chord kernels build each
slab's Hamiltonian as they go, so the stack that made the batched path
memory-bound is never allocated; and a chord is a palindrome, so half its
operators were being computed twice.  Against a Python loop over the same
2000 energies, an Earth scan is ~38x, ~10x and ~11x quicker at two, three
and four flavors, and ~460x, ~120x and ~87x against that loop on the
NumPy path.  The palindrome is ~1.4x, ~1.5x and ~1.8x of that; see
`USE_PALINDROME`.

Four flavors gains the most, and not because the kernel is cleverer
there: the NumPy path has the furthest to fall.  Its expansion needs a
quartic, a Newton refinement of the four roots against the matrix, and
a Newton-form reconstruction, which as whole-array operations is some
forty passes over the stack; done one element at a time none of it
leaves the registers.

These are one machine on one day, and they move by tens of per cent
between runs; read them as the shape of the gain rather than as
constants.  The figures quoted for 1.6.0 in ``CHANGELOG.md`` came from a
different session and differ by up to a factor of two --- which is why
notebook 09 measures the comparison when it runs, on whatever machine
is running it, rather than repeating a number from here.

Costs, so that the trade is visible
-----------------------------------

* importing Numba takes about 140 ms, against 65 ms for NumPy alone;
* the first call compiles, which takes a few seconds.  The kernels are
  declared with ``cache=True``, so that cost is paid once per machine
  and later runs load the compiled code from disk in milliseconds.

Both are why this was an optional extra until 1.13.0.  What changed the
argument is that neither cost falls where it would be felt: they are
paid once, and never by a scalar call at all.  The scalar path is
deliberately left alone --- a single probability takes about
8 microseconds and never enters a kernel --- so the compilation pause
belongs only to the batched paths, which are the ones that go on to save
whole seconds.

Turning it off
--------------

Set ``fastkernels.USE_NUMBA = False`` to force the NumPy path even when
Numba is installed --- useful for checking that the two agree, which is
what the test suite does.

Routine listings
----------------

    * available - Whether the compiled kernels can be used at all
    * worthwhile - Whether a stack is large enough to be worth compiling
    * worthwhile_slabs - The same question for a slab sequence
    * probabilities_2nu_kernel - Two-flavor probabilities for a stack
    * probabilities_3nu_kernel - Three-flavor probabilities for a stack
    * probabilities_4nu_kernel - Four-flavor probabilities for a stack
    * evolution_operator_3nu_kernel - Three-flavor U_3 for a stack
    * slab_product_3nu_kernel - U_3 composed across slabs
    * evolution_operator_4nu_kernel - Four-flavor U_4 for a stack
    * slab_product_4nu_kernel - U_4 composed across slabs
    * evolution_operator_2nu_kernel - Two-flavor U_2 for a stack
    * slab_product_2nu_kernel - U_2 composed across slabs
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

__all__ = ['HAVE_NUMBA', 'USE_NUMBA', 'USE_PALINDROME',
           'MIN_BATCH', 'MIN_SLAB_BATCH', 'MIN_MIRROR_SLABS',
           'PARALLEL_THRESHOLD',
           'available', 'worthwhile', 'worthwhile_slabs',
           'worthwhile_mirror', 'palindromic',
           'probabilities_2nu_kernel', 'probabilities_3nu_kernel',
           'probabilities_4nu_kernel', 'evolution_operator_3nu_kernel',
           'slab_product_3nu_kernel',
           'evolution_operator_4nu_kernel',
           'slab_product_4nu_kernel',
           'evolution_operator_2nu_kernel',
           'slab_product_2nu_kernel',
           'slab_product_2nu_batch_kernel',
           'slab_product_3nu_batch_kernel',
           'slab_product_4nu_batch_kernel',
           'earth_chords_2nu_kernel', 'earth_chords_3nu_kernel',
           'earth_chords_4nu_kernel']

from typing import Callable

import cmath
import math

import numpy as np

try:
    from numba import njit, prange
    HAVE_NUMBA = True
except ImportError:                                       # pragma: no cover
    HAVE_NUMBA = False

USE_NUMBA = True
r"""bool: Module-level switch.

Set to ``False`` to force the NumPy path even when Numba is installed.
`available` reports the two together.
"""

USE_PALINDROME = True
r"""bool: Module-level switch.

Whether a slab sequence that reads the same from either end may be
composed at roughly two thirds of the cost, by computing each distinct
operator once instead of twice.  True by default: every Earth chord
qualifies, and the saving is between 1.3x and 1.8x.

Set it to ``False`` to compose every sequence in full.  That is not a
correctness switch --- the two agree to a few times 1e-15, and the
mirrored composer's departure from unitarity is if anything slightly
smaller --- but the two orderings round differently, so it is the way to
ask for the plain left-to-right product when a comparison needs one, and
the way to establish that a discrepancy is or is not the palindrome's
doing.

Whether a given sequence *is* a palindrome is a separate question, which
`palindromic` answers; this switch only decides whether it is worth
asking.

.. versionadded:: 1.12.0
"""

MIN_BATCH_HERMITICITY = 256
r"""int: Module-level constant.

The smallest stack for which `hermiticity_offender` beats the NumPy
comparisons it replaces.  Below it the compiled call's own dispatch
dominates, so :func:`oscprob3nu._check_hermitian` and its siblings stay on
NumPy; above it the saving is most of the check.

Why the check needed a kernel at all is worth recording.  It was written
against the cost of the *expansion*, and was carefully made cheap by that
measure --- comparing the independent pairs on real and imaginary views
rather than forming :math:`H - H^\dagger`, worth about three times.  Then
the compiled kernels made the expansion roughly ten times cheaper and left
the check untouched, at which point validating a 100 000-element stack cost
70% of a three-flavor probability call and 82% of a two-flavor one: more
than the physics it guards.  The check's own docstring had predicted
exactly that failure, about the implementation it replaced.
"""

MIN_BATCH = {2: 50000, 3: 1, 4: 1}
r"""dict: Module-level constant.

The smallest stack for which the compiled kernel is worth using, by
number of flavors.  A backend that is sometimes slower than the path it
replaces is worse than no backend, so these are measured rather than
assumed.

For three flavors the kernel wins at every size, by between two and
sixteen times, so the threshold is one.  Four flavors is the same story
only more so, and for a reason worth stating: the NumPy path there has
no short-stack shortcut to fall back on --- :mod:`oscprob4nu` has no
separate scalar closed form, so even a stack of one pays for the whole
array machinery, a batched determinant and all.  Measured by alternating
the two paths through :func:`oscprob4nu.probabilities_4nu` and taking the
best of nine rounds each, the kernel leads by 15x at a single element,
falls to 5x just below `PARALLEL_THRESHOLD` where it is still
single-threaded, and settles at 18x once the threads are in use.  It is
never behind, so the threshold is one.

For two flavors it does not:
that expansion reduces to a square root and a sine per element, which
NumPy already does about as well as compiled code can, and the kernel
additionally has to materialise the Hamiltonian stack --- which for a
scan over baselines is the same matrix repeated, costing 2.5 ms to copy
at two hundred thousand points.  Measured by alternating the two
paths and taking the best of nine rounds each, the crossover sits
between twenty and fifty thousand elements: at twenty thousand NumPy is
still ahead by a few per cent, at fifty thousand the kernel leads by
1.3x and it grows slowly from there.  The threshold is set at the first
size where the kernel is unambiguously ahead, since the region around
the crossover is broad and varies between machines.
"""

MIN_SLAB_BATCH = {2: 1, 3: 1, 4: 1}
r"""dict: Module-level constant.

The smallest slab sequence for which the compiled *product* kernel is
worth using --- one, at every flavor count, which is to say always.

This is deliberately not `MIN_BATCH`, and the difference is the point.
`MIN_BATCH` weighs compiled arithmetic against NumPy arithmetic, which
is why two flavors sits at fifty thousand: that expansion reduces to a
square root and a sine per element, and NumPy does those about as well as
compiled code can.  The slab product is not that comparison.  It replaces
compiled arithmetic *and a Python loop of dispatched matrix products*
over a stack that had to be materialised first, so the alternative
carries a cost per slab that has nothing to do with the flavor count.

Measured by alternating the two paths, best of five rounds each, at
sequences of 1 to 256 slabs: the kernel leads by 137x, 225x and 187x at
a single slab, and by 59x, 13x and 7x at 256, at two, three and four
flavors.  It is never behind anywhere in that range, and the margin
*narrows* with length --- the opposite of `MIN_BATCH`'s kernels, because
here the fixed cost being avoided is the caller's rather than the
kernel's.

Reusing `MIN_BATCH` here would have left two flavors on the NumPy path
for every Earth crossing, since a chord is a hundred-odd slabs and the
threshold is fifty thousand.

.. versionadded:: 1.12.0
"""

MIN_MIRROR_SLABS = {2: 16, 3: 16, 4: 16}
r"""dict: Module-level constant.

The shortest chord for which composing a palindrome at half cost is
worth it, by number of flavors.

The mirrored composer computes each distinct operator once and uses it
twice, which halves the expansions --- about two thirds of a slab's
work --- but it accumulates two running products rather than one, so the
matrix multiplications stay at one per slab rather than falling with
them.  That is a good trade only once there are enough slabs for the
expansions to dominate the fixed cost of carrying the second product.

Measured over chords of 7 to 960 slabs at three flavors, interleaved
against the ordinary composer: behind or level at 7, ahead by 1.10x at
15, and between 1.43x and 1.79x from 56 slabs upwards.  The threshold is
set at sixteen, the first length where it is clearly ahead rather than
within noise.  An Earth chord has at least 7 slabs at the coarsest
subdivision and 120 at the default, so this is about protecting the
shallow, coarse corner rather than the common case.

.. versionadded:: 1.12.0
"""

PARALLEL_THRESHOLD = 256
r"""int: Module-level constant.

Stacks with at least this many elements are spread over the available
cores; smaller ones run in a single thread, because below roughly this
size the cost of waking the thread pool exceeds what it saves.
"""


def available() -> bool:
    r"""Returns whether the compiled kernels can be used at all.

    True when Numba was imported successfully *and* `USE_NUMBA` has not
    been turned off.  Whether they are *worth* using for a given stack
    is a separate question; see `worthwhile`.

    .. versionadded:: 1.6.0

    Returns
    -------
    bool
        Whether `probabilities_2nu_kernel`, `probabilities_3nu_kernel`
        and `probabilities_4nu_kernel` may be called.
    """
    return HAVE_NUMBA and USE_NUMBA


def worthwhile_slabs(n_flavors: int, size: int) -> bool:
    r"""Returns whether the compiled product kernel should be used.

    .. versionadded:: 1.12.0

    The slab counterpart of `worthwhile`, against `MIN_SLAB_BATCH`
    rather than `MIN_BATCH`.  The two thresholds answer different
    questions and one is not a good default for the other; see
    `MIN_SLAB_BATCH`.

    Parameters
    ----------
    n_flavors : int
        Number of neutrino flavors, 2, 3, or 4.
    size : int
        Number of slabs in the sequence.

    Returns
    -------
    bool
        Whether `slab_product_2nu_kernel` and its siblings may be used.
    """
    return available() and size >= MIN_SLAB_BATCH.get(n_flavors, 1)


def palindromic(*arrays: np.ndarray) -> bool:
    r"""Returns whether every array given reads the same both ways.

    .. versionadded:: 1.12.0

    A sequence of slabs whose Hamiltonians and widths are both
    palindromic has ``U_j = U_{n-1-j}``, so half its operators are
    recomputations of the other half and the product can be composed at
    roughly two thirds of the cost --- see
    `earth_chords_3nu_kernel`.  A chord through a spherically symmetric
    Earth is always such a sequence, because it meets every radius
    twice, but nothing here is specific to the Earth: a castle wall
    built symmetrically, or any hand-built profile that reads the same
    from either end, qualifies too.

    The comparison is exact, deliberately.  The saving relies on the two
    mirrored operators being *identical*, which follows from identical
    inputs and from nothing weaker; a tolerance here would silently
    return a different answer for a nearly-symmetric profile, which is
    the one thing an optimisation must never do.  It is
    `earth._earth_slabs_cached`'s business to make the Earth's chords
    exactly symmetric rather than nearly so, and it does.

    Parameters
    ----------
    arrays : numpy.ndarray
        Arrays to test, given as separate arguments and each reversed
        along its first axis.  An empty call, or any array shorter than
        two entries, is trivially palindromic.

    Returns
    -------
    bool
        Whether every array equals its own reverse exactly.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np
        import fastkernels

        print(fastkernels.palindromic(np.array([1.0, 2.0, 1.0])))
        print(fastkernels.palindromic(np.array([1.0, 2.0, 3.0])))
    """
    for array in arrays:
        a = np.asarray(array)
        if a.shape[0] > 1 and not np.array_equal(a, a[::-1]):
            return False

    return True


def worthwhile_mirror(n_flavors: int, size: int) -> bool:
    r"""Returns whether the palindromic composer should be used.

    .. versionadded:: 1.12.0

    Whether the chord *is* a palindrome is a separate question, and
    `palindromic` answers it; this one asks only whether the saving is
    worth having at this length.  It is not below a handful of slabs:
    the mirrored composer accumulates two running products instead of
    one, so it trades a matrix multiplication per slab against half an
    expansion, and that trade only pays once there are enough slabs for
    the expansions to dominate.  Measured across 7 to 960 slabs, it
    leads from about fifteen slabs upwards, by 1.4x to 1.8x, and is
    within noise below that.

    Parameters
    ----------
    n_flavors : int
        Number of neutrino flavors, 2, 3, or 4.
    size : int
        Number of slabs in the chord.

    Returns
    -------
    bool
        Whether to use the mirrored composer.
    """
    return (USE_PALINDROME and available()
            and size >= MIN_MIRROR_SLABS.get(n_flavors, 1))


def worthwhile(n_flavors: int, size: int) -> bool:
    r"""Returns whether the compiled kernel should be used for a stack.

    The kernels are only used where they have been measured to win.
    Below the per-flavor threshold in `MIN_BATCH` the NumPy path is
    quicker, and using the kernel anyway would make the compiled backend
    a pessimisation for those calls --- which matters more now that it
    arrives with the package rather than being asked for.

    .. versionadded:: 1.6.0

    Parameters
    ----------
    n_flavors : int
        Number of neutrino flavors, 2, 3, or 4.
    size : int
        Number of elements in the stack.

    Returns
    -------
    bool
        Whether to call the corresponding kernel.
    """
    return available() and size >= MIN_BATCH.get(n_flavors, 1)


if HAVE_NUMBA:                                          # pragma: no branch

    SQRT3 = math.sqrt(3.0)
    SQRT3_INV = 1.0/SQRT3
    TWO_SQRT3_INV = 2.0*SQRT3_INV
    DEGENERACY_TOL = 1.0e-12

    @njit(cache=True, parallel=True)
    def _hermiticity_scale(h_stack, n_flavors):
        r"""Returns ``(largest magnitude, count of non-finite entries)``.

        The tolerance is relative to one global scale rather than to each
        element --- which is what the NumPy version this replaces does, and
        the reason the comparisons cannot be fused into this pass.

        Non-finite entries are counted rather than left to poison the scale.
        Comparisons against a nan are all false, so a nan would never
        survive ``if real > local`` to reach a finiteness test afterwards,
        and a Hamiltonian both non-finite *and* non-Hermitian would pass a
        check whose whole purpose is to refuse the second.  The NumPy path
        gets this free from ``np.max`` propagating nan; here it is explicit.
        """
        biggest = 0.0
        non_finite = 0
        for element in prange(h_stack.shape[0]):
            local = 0.0
            local_bad = 0
            for i in range(n_flavors):
                for j in range(n_flavors):
                    entry = h_stack[element, i, j]
                    real = abs(entry.real)
                    imaginary = abs(entry.imag)
                    if not (math.isfinite(real) and math.isfinite(imaginary)):
                        local_bad = 1
                    if real > local:
                        local = real
                    if imaginary > local:
                        local = imaginary
            biggest = max(biggest, local)
            non_finite += local_bad
        return biggest, non_finite

    @njit(cache=True)
    def _hermiticity_first_bad(h_stack, n_flavors, tolerance):
        r"""Returns ``element*16 + i*4 + j`` for the first offender, or -1.

        Serial on purpose, and it costs nothing to be.  It returns on the
        first offending element, so a stack that fails is settled almost
        immediately, and a stack that passes is one sweep of reads with no
        allocation --- against the dozen NumPy reductions this replaces,
        each of which allocated a temporary the size of the stack and then
        reduced it away.  Running it across cores would need a minimum
        reduction guarded by a conditional to keep the reported element
        independent of the core count, which Numba's parfor pass rejects,
        and would buy a pass that is already not the expensive one.

        ``i == j`` means a diagonal entry whose imaginary part does not
        vanish.
        """
        for element in range(h_stack.shape[0]):
            for i in range(n_flavors):
                if abs(h_stack[element, i, i].imag) > tolerance:
                    return element*16 + i*4 + i
                for j in range(i+1, n_flavors):
                    upper = h_stack[element, i, j]
                    lower = h_stack[element, j, i]
                    if (abs(upper.real - lower.real) > tolerance
                            or abs(upper.imag + lower.imag) > tolerance):
                        return element*16 + i*4 + j
        return -1

    def hermiticity_offender(
        h_stack: np.ndarray,
        n_flavors: int,
        relative_tol: float
    ) -> tuple:
        r"""Returns ``(non_finite, element, i, j)`` for a stack.

        The compiled counterpart of the stack branch of
        :func:`oscprob3nu._check_hermitian` and its siblings, which own the
        error messages; this only finds the offender.  The verdict is
        identical, including the global relative tolerance and the refusal
        of non-finite entries --- what changes is that the scan happens in
        one compiled pass per stage rather than a dozen NumPy reductions,
        each allocating a temporary the size of the stack.

        Parameters
        ----------
        h_stack : numpy.ndarray
            Hamiltonians, of shape ``(..., n, n)``, with at least one
            leading axis.
        n_flavors : int
            The number of flavors, 2, 3 or 4.
        relative_tol : float
            The caller's ``_HERMITICITY_TOL``, applied relative to the
            largest magnitude anywhere in the stack.

        Returns
        -------
        tuple
            ``(non_finite, element, i, j)``.  `non_finite` is True if any
            entry is not finite, in which case the indices are ``-1``.
            Otherwise `element` is ``-1`` when the stack is Hermitian, and
            the flat index of the first offending element when it is not,
            with `i` and `j` its entry --- equal for a diagonal entry whose
            imaginary part does not vanish.
        """
        flat = h_stack.reshape(-1, n_flavors, n_flavors)
        biggest, non_finite = _hermiticity_scale(flat, n_flavors)
        if non_finite > 0:
            return True, -1, -1, -1

        tolerance = relative_tol*biggest if biggest > 0.0 else relative_tol
        code = int(_hermiticity_first_bad(flat, n_flavors, tolerance))
        if code < 0:
            return False, -1, -1, -1
        return False, code//16, (code % 16)//4, code % 4

    DD_SPLIT = 134217729.0
    r"""float: :math:`2^{27} + 1`, Dekker's splitting constant.

    Multiplying by it and subtracting isolates the leading 26 bits of a
    float64, so that the product of two halves is representable exactly and
    `_two_prod` can return the rounding error it discarded.
    """

    def _dd_reciprocal_of_three() -> tuple:
        r"""Returns 1/3 as a double-double pair, exact to 3e-33.

        Formed once at import, by the same compensated arithmetic the kernel
        uses: ``3*hi`` is split exactly into ``p + e`` by `_two_prod`, so
        what 1/3 is missing after the first limb is ``((1 - p) - e)/3``.

        It exists because ``2/3`` is the one irrational-in-binary constant
        the quartic needs, and lifting a *rounded* 2/3 into double-double
        would zero the low limb and drag c1 back to float64 accuracy --- the
        1.1e-07 stall this whole route exists to remove, wearing the
        disguise of a solver failure.  Multiplying by this pair costs no
        division at all, where `_dd_div` costs three.
        """
        hi = 1.0/3.0
        p, e = _two_prod(3.0, hi)
        return hi, ((1.0 - p) - e)/3.0

    DD_SWEEPS = 1
    r"""int: Aberth sweeps taken over the four latent roots.

    One, measured: every case in `tests/stiff_reference.json` gives an
    identical answer at one, two and three sweeps, to the last bit.  A
    single sweep suffices because the start is already good --- the
    eigensolver's quartet with the residual-trace shift removed in
    double-double rather than in float64, so the low limb the exact
    traceless-ing computed survives into the iteration.

    This was two for exactly as long as that subtraction was done in
    float64, which put the start an ulp out and cost a sweep to recover.
    See :data:`oscprob4nu.DD_SWEEPS`, which records the mistake, because
    one sweep did genuinely measure 3.9e-16 back then and the number
    looked like evidence the iteration needed the second.
    """

    @njit(cache=True, inline='always')
    def _two_sum(a, b):
        r"""Returns ``a + b`` as a rounded sum and the error it dropped.

        Knuth's compensated addition: `s` is the float64 sum and `e` is
        exactly :math:`a + b - s`, so the pair represents the sum with
        nothing lost.  Six flops, and no assumption about which operand is
        larger --- `_quick_two_sum` is the cheaper version that does assume
        it.
        """
        s = a + b
        bb = s - a
        return s, (a - (s - bb)) + (b - bb)

    @njit(cache=True, inline='always')
    def _quick_two_sum(a, b):
        r"""Returns ``a + b`` and its error, given :math:`|a| \geq |b|`."""
        s = a + b
        return s, b - (s - a)

    @njit(cache=True, inline='always')
    def _two_prod(a, b):
        r"""Returns ``a*b`` as a rounded product and the error it dropped.

        Dekker's algorithm: split both operands at `DD_SPLIT` into halves
        narrow enough to multiply exactly, then assemble the four partial
        products and subtract the rounded one.
        """
        p = a*b
        t = DD_SPLIT*a
        ah = t - (t - a)
        al = a - ah
        t = DD_SPLIT*b
        bh = t - (t - b)
        bl = b - bh
        return p, ((ah*bh - p) + ah*bl + al*bh) + al*bl

    @njit(cache=True, inline='always')
    def _dd_add(xh, xl, yh, yl):
        r"""Returns the double-double sum of ``(xh, xl)`` and ``(yh, yl)``."""
        s, e = _two_sum(xh, yh)
        return _quick_two_sum(s, e + xl + yl)

    @njit(cache=True, inline='always')
    def _dd_sub(xh, xl, yh, yl):
        r"""Returns the double-double difference, ``x - y``."""
        return _dd_add(xh, xl, -yh, -yl)

    @njit(cache=True, inline='always')
    def _dd_mul(xh, xl, yh, yl):
        r"""Returns the double-double product, ``x*y``."""
        p, e = _two_prod(xh, yh)
        return _quick_two_sum(p, e + xh*yl + xl*yh)

    @njit(cache=True, inline='always')
    def _dd_div(xh, xl, yh, yl):
        r"""Returns the double-double quotient, ``x/y``.

        Three float64 divisions, each correcting the remainder left by the
        one before.  The most expensive dd primitive by some margin, which
        is why the root iteration below is Aberth's with the divisions
        counted rather than anything more elaborate.
        """
        q1 = xh/yh
        th, tl = _dd_mul(yh, yl, q1, 0.0)
        rh, rl = _dd_sub(xh, xl, th, tl)
        q2 = rh/yh
        th, tl = _dd_mul(yh, yl, q2, 0.0)
        rh, rl = _dd_sub(rh, rl, th, tl)
        q3 = rh/yh
        sh, sl = _quick_two_sum(q1, q2)
        return _dd_add(sh, sl, q3, 0.0)

    # Bound here rather than beside `DD_SWEEPS` because forming it needs
    # `_two_prod`, which is defined above it.  Numba bakes a module global
    # into the compiled code as a constant, which is exactly what is wanted.
    THIRD_HI, THIRD_LO = _dd_reciprocal_of_three()

    @njit(cache=True, inline='always')
    def _entries_3nu(h_matrix, L):
        r"""Returns the nine entries of :math:`U_3(L)`, row by row.

        A transcription of the scalar path in :mod:`oscprob3nu`: the
        SU(3) coefficients, the sparse star product, the two invariants,
        the latent roots with the same degeneracy handling, and the nine
        entries of the evolution operator.

        Factored out of `_one_3nu` so that the probability kernel and the
        evolution-operator kernel share it.  The probability kernel used
        to compute these entries and square them on the next line, which
        meant `slabs` and `earth` --- which need the *operators*, to
        multiply them --- had no compiled path at all and silently ran
        the NumPy one.  The arithmetic here is unchanged; only its last
        step now has two callers.

        Returned as a tuple rather than written into an array because
        ``inline='always'`` folds it into both call sites, so nothing is
        materialised and the probability path costs exactly what it did.
        """
        h0 = h_matrix[0, 1].real
        h1 = -h_matrix[0, 1].imag
        h2 = (h_matrix[0, 0] - h_matrix[1, 1]).real/2.0
        h3 = h_matrix[0, 2].real
        h4 = -h_matrix[0, 2].imag
        h5 = h_matrix[1, 2].real
        h6 = -h_matrix[1, 2].imag
        h7 = (h_matrix[0, 0] + h_matrix[1, 1]
              - 2.0*h_matrix[2, 2]).real*SQRT3/6.0

        # (h*h)_k, the sparse expansion of d_ijk h_j h_k
        s0 = TWO_SQRT3_INV*h0*h7 + h3*h5 + h4*h6
        s1 = TWO_SQRT3_INV*h1*h7 - h3*h6 + h4*h5
        s2 = TWO_SQRT3_INV*h2*h7 + 0.5*(h3*h3 + h4*h4 - h5*h5 - h6*h6)
        s3 = h0*h5 - h1*h6 + h2*h3 - SQRT3_INV*h3*h7
        s4 = h0*h6 + h1*h5 + h2*h4 - SQRT3_INV*h4*h7
        s5 = h0*h3 + h1*h4 - h2*h5 - SQRT3_INV*h5*h7
        s6 = h0*h4 - h1*h3 - h2*h6 - SQRT3_INV*h6*h7
        s7 = (SQRT3_INV*(h0*h0 + h1*h1 + h2*h2 - h7*h7)
              - SQRT3_INV/2.0*(h3*h3 + h4*h4 + h5*h5 + h6*h6))

        # |h|^2 and <h>
        hsq = (h0*h0 + h1*h1 + h2*h2 + h3*h3
               + h4*h4 + h5*h5 + h6*h6 + h7*h7)
        hcu = (h0*s0 + h1*s1 + h2*s2 + h3*s3
               + h4*s4 + h5*s5 + h6*s6 + h7*s7)

        if hsq <= 0.0:
            # Proportional to the identity: U3 = 1
            u0 = 1.0 + 0.0j
            c1 = c2 = c3 = c4 = c5 = c6 = c7 = c8 = 0.0j
        else:
            root = math.sqrt(hsq)
            pre = 2.0*SQRT3_INV*root

            arg = -SQRT3*hcu/(hsq*root)
            if arg <= -1.0:
                chi = math.pi
            elif arg >= 1.0:
                chi = 0.0
            else:
                chi = math.acos(arg)

            psi0 = pre*math.cos((chi + 2.0*math.pi)/3.0)
            psi1 = pre*math.cos((chi + 4.0*math.pi)/3.0)
            psi2 = pre*math.cos((chi + 6.0*math.pi)/3.0)

            gap01 = abs(psi0-psi1)
            gap02 = abs(psi0-psi2)
            gap12 = abs(psi1-psi2)
            smallest = gap01
            if gap02 < smallest:
                smallest = gap02
            if gap12 < smallest:
                smallest = gap12

            if smallest <= DEGENERACY_TOL*root:
                # Doubly degenerate root: the two-projector form
                if gap01 <= gap02 and gap01 <= gap12:
                    psi_deg = 0.5*(psi0+psi1)
                    psi_odd = psi2
                elif gap02 <= gap12:
                    psi_deg = 0.5*(psi0+psi2)
                    psi_odd = psi1
                else:
                    psi_deg = 0.5*(psi1+psi2)
                    psi_odd = psi0
                exp_deg = cmath.rect(1.0, L*psi_deg)
                exp_odd = cmath.rect(1.0, L*psi_odd)
                weight = (exp_odd-exp_deg)/(psi_deg-psi_odd)
                u0 = exp_deg + weight*psi_deg
                factor = -1.0j*weight
                c1 = factor*h0
                c2 = factor*h1
                c3 = factor*h2
                c4 = factor*h3
                c5 = factor*h4
                c6 = factor*h5
                c7 = factor*h6
                c8 = factor*h7
            else:
                exp0 = cmath.rect(1.0, L*psi0)
                exp1 = cmath.rect(1.0, L*psi1)
                exp2 = cmath.rect(1.0, L*psi2)
                w0 = exp0/(3.0*psi0*psi0 - hsq)
                w1 = exp1/(3.0*psi1*psi1 - hsq)
                w2 = exp2/(3.0*psi2*psi2 - hsq)
                weighted = w0*psi0 + w1*psi1 + w2*psi2
                total = w0 + w1 + w2
                u0 = (exp0+exp1+exp2)/3.0
                c1 = 1.0j*(weighted*h0 - total*s0)
                c2 = 1.0j*(weighted*h1 - total*s1)
                c3 = 1.0j*(weighted*h2 - total*s2)
                c4 = 1.0j*(weighted*h3 - total*s3)
                c5 = 1.0j*(weighted*h4 - total*s4)
                c6 = 1.0j*(weighted*h5 - total*s5)
                c7 = 1.0j*(weighted*h6 - total*s6)
                c8 = 1.0j*(weighted*h7 - total*s7)

        eighth = c8/SQRT3
        u_ee = u0 + 1.0j*(c3 + eighth)
        u_em = 1.0j*c1 + c2
        u_et = 1.0j*c4 + c5
        u_me = 1.0j*c1 - c2
        u_mm = u0 - 1.0j*(c3 - eighth)
        u_mt = 1.0j*c6 + c7
        u_te = 1.0j*c4 - c5
        u_tm = 1.0j*c6 - c7
        u_tt = u0 - 2.0j*eighth

        return (u_ee, u_em, u_et, u_me, u_mm, u_mt, u_te, u_tm, u_tt)

    @njit(cache=True, inline='always')
    def _one_3nu(h_matrix, L, out, n):
        r"""Writes the nine probabilities for one Hamiltonian into
        ``out[n]``.
        """
        (u_ee, u_em, u_et,
         u_me, u_mm, u_mt,
         u_te, u_tm, u_tt) = _entries_3nu(h_matrix, L)

        # P_ab = |U_ba|^2, initial flavor slowest
        out[n, 0] = u_ee.real*u_ee.real + u_ee.imag*u_ee.imag
        out[n, 1] = u_me.real*u_me.real + u_me.imag*u_me.imag
        out[n, 2] = u_te.real*u_te.real + u_te.imag*u_te.imag
        out[n, 3] = u_em.real*u_em.real + u_em.imag*u_em.imag
        out[n, 4] = u_mm.real*u_mm.real + u_mm.imag*u_mm.imag
        out[n, 5] = u_tm.real*u_tm.real + u_tm.imag*u_tm.imag
        out[n, 6] = u_et.real*u_et.real + u_et.imag*u_et.imag
        out[n, 7] = u_mt.real*u_mt.real + u_mt.imag*u_mt.imag
        out[n, 8] = u_tt.real*u_tt.real + u_tt.imag*u_tt.imag

    @njit(cache=True, inline='always')
    def _one_3nu_u(h_matrix, L, out, n):
        r"""Writes the nine entries of :math:`U_3(L)` into ``out[n]``.

        Row-major and indexed ``(final, initial)``, so that reshaping to
        ``(3, 3)`` gives the same matrix `oscprob3nu` returns --- *not*
        the flavor order the probabilities use, which runs the initial
        index slowest.  The two orderings differ by a transpose, and
        conflating them is the obvious way to get this wrong.
        """
        (u_ee, u_em, u_et,
         u_me, u_mm, u_mt,
         u_te, u_tm, u_tt) = _entries_3nu(h_matrix, L)

        out[n, 0] = u_ee
        out[n, 1] = u_em
        out[n, 2] = u_et
        out[n, 3] = u_me
        out[n, 4] = u_mm
        out[n, 5] = u_mt
        out[n, 6] = u_te
        out[n, 7] = u_tm
        out[n, 8] = u_tt

    @njit(cache=True, inline='always')
    def _one_2nu(h_matrix, L, out, n):
        r"""Writes the four probabilities for one Hamiltonian into
        ``out[n]``.

        The coefficients of a Hermitian 2x2 Hamiltonian are real, so the
        transition probability follows from the Hamiltonian directly and
        the survival probability is its complement.
        """
        h0 = h_matrix[0, 1].real
        h1 = -h_matrix[0, 1].imag
        h2 = (h_matrix[0, 0] - h_matrix[1, 1]).real/2.0

        hsq = h0*h0 + h1*h1 + h2*h2
        if hsq <= 0.0:
            # Proportional to the identity: no flavor transitions
            p_em = 0.0
        else:
            sin_phase = math.sin(math.sqrt(hsq)*L)
            p_em = (h0*h0 + h1*h1)/hsq*sin_phase*sin_phase

        out[n, 0] = 1.0 - p_em
        out[n, 1] = p_em
        out[n, 2] = p_em
        out[n, 3] = 1.0 - p_em

    @njit(cache=True, inline='always')
    def _chi_4nu(traceless, psi, scratch):
        r"""Returns the characteristic polynomial
        :math:`\chi(\psi) = \det(\psi\mathbb{1} - \tilde{H})`.

        Gaussian elimination with partial pivoting, written out for a
        four-by-four in the caller's scratch buffer: the same
        factorisation LAPACK's ``zgetrf`` performs for
        :func:`numpy.linalg.det` on the NumPy path, with no allocation
        and no call.

        The obvious alternative --- a Laplace expansion in the six
        two-by-two minors of the first two rows, thirty products, no
        pivoting and no branches --- was written first and measured
        against ``mpmath`` at sixty digits.  It is **5.9x cheaper**, and
        it was still rejected.

        The reason is that this determinant is evaluated *at a root*,
        where it is meant to vanish.  On the stiff 3+1 spectrum the true
        value sits some seventeen orders of magnitude below the products
        being summed, so an expansion that cancels them only at the end
        has no significant digits left, while elimination cancels while
        the entries are still full precision.  On the clustered roots,
        where :math:`\chi'` is :math:`6 \times 10^{-35}`, the expansion
        was a thousand times the less accurate, and it refined those
        roots to :math:`4 \times 10^{-15}` relative against
        :math:`6 \times 10^{-16}` here; on a spectrum whose cluster is
        :math:`10^{-3}` wide the gap widens to 54x.  The refined figure
        :data:`oscprob4nu.POLISH_ROOTS` tabulates is
        :math:`1.1 \times 10^{-16}`, and a backend that quietly delivers
        forty times that whenever an optional dependency happens to be
        installed makes that table false.

        What the measurement did **not** show is any of this reaching
        the probabilities, and the honest record is that no test here
        distinguishes the two.  Below :math:`\psi L \sim 1` both sit on
        the one-ulp floor; above it the Newton-form reconstruction
        cancels by :math:`\sim 10^6` and swamps them both, and which
        scores better is then noise --- on the stiff spectrum at 1300 km
        the *rejected* expansion won, 3.1e-11 against 1.5e-10.  The case
        for elimination is fidelity to the roots the NumPy path
        computes, not a demonstrated gain in the numbers handed back.
        It costs about 40% of the kernel's serial runtime, which is
        cheap against the 18x the kernel wins overall.

        The result is real for a Hermitian argument, so only the real
        part is returned, exactly as the NumPy path takes ``.real`` of
        :func:`numpy.linalg.det`.
        """
        for i in range(4):
            for j in range(4):
                scratch[i, j] = -traceless[i, j]
            scratch[i, i] = psi - traceless[i, i]

        sign = 1.0
        for k in range(3):
            # The pivot LAPACK would choose: ``izamax`` ranks by
            # |Re| + |Im| rather than by the modulus
            pivot_row = k
            largest = abs(scratch[k, k].real) + abs(scratch[k, k].imag)
            for i in range(k+1, 4):
                candidate = abs(scratch[i, k].real) + abs(scratch[i, k].imag)
                if candidate > largest:
                    largest = candidate
                    pivot_row = i

            if pivot_row != k:
                sign = -sign
                for j in range(k, 4):
                    swap = scratch[k, j]
                    scratch[k, j] = scratch[pivot_row, j]
                    scratch[pivot_row, j] = swap

            pivot = scratch[k, k]
            if pivot == 0.0:
                return 0.0

            for i in range(k+1, 4):
                multiplier = scratch[i, k]/pivot
                for j in range(k+1, 4):
                    scratch[i, j] -= multiplier*scratch[k, j]

        return sign*(scratch[0, 0]*scratch[1, 1]
                     * scratch[2, 2]*scratch[3, 3]).real

    @njit(cache=True)
    def _latent_roots_dd(traceless):
        r"""Returns the four latent roots, to the last float64 bit.

        The transcription of :func:`oscprob4nu._latent_roots_dd`.  Three
        invariants compress a matrix into three numbers, and in float64 that
        compression is where a stiff spectrum loses the separation between
        its clustered roots: the amplification from coefficients to roots
        was measured at 2.3e9, so a 1e-16 coefficient becomes a 1e-7 root.
        The invariants are therefore formed in double-double arithmetic,
        where a 1e-32 coefficient error amplifies to 1e-23 and the roots are
        limited by float64 output rounding instead --- 3.6e-17 worst over
        `tests/stiff_reference.json`, against 3.9e-16 for the eigensolver
        with a Newton step and 2.2e-07 for the closed form the invariants
        used to be handed to.

        Four things here were each chosen against measurement, and each
        would be easy to undo by accident:

        `traceless` arrives with its trace removed in float64, which leaves
        a residue :math:`\tau` of order 1e-23.  Removing that residue
        exactly in dd is not a refinement, it is the difference between a
        quartic that describes this matrix and one that does not: the
        invariants would otherwise belong to a matrix whose trace is not
        quite zero while the quartic pins its cubic coefficient to exactly
        zero, and the two disagree by up to 3.8e-7.

        :math:`\tilde{H}` is Hermitized exactly --- ``(h[i,j] +
        conj(h[j,i]))/2``, representable because `_two_sum` of two float64s
        loses nothing and halving only shifts exponents.  That makes
        :math:`\tilde{H}^2` exactly Hermitian, so only its upper triangle is
        formed, ten entries instead of sixteen.  Mirroring *without* the
        Hermitization looks equivalent and costs :math:`I_3` and
        :math:`I_4` about 1e-17, because a Hamiltonian built as
        :math:`U M^2 U^\dagger` is Hermitian only to rounding and the mirror
        discards the asymmetry silently.  Discarding it deliberately is
        free: an anti-Hermitian perturbation moves a *real* eigenvalue only
        at second order, since :math:`\langle v|\delta|v \rangle` is
        imaginary.

        :math:`\mathrm{Tr}(\tilde{H}^2)` is taken off :math:`\tilde{H}`
        itself as :math:`\sum_{ij} |\tilde{H}_{ij}|^2` and
        :math:`\mathrm{Tr}(\tilde{H}^4)` off :math:`\tilde{H}^2` as
        :math:`\sum_{ij} |S_{ij}|^2`.  Both are sums of squares, with no
        cancellation to lose digits to.

        The start is the eigensolver's, because what a start must supply is
        neither accuracy nor proximity but four distinct basins.  Euler's
        closed form is twice as fast and exact on every stiff case, and
        still wrong here:
        on a cluster Aberth converges *linearly*, at ratio one half, so from
        1e-7 five sweeps reach 3.8e-9 and it would need some thirty.  A
        backward-stable Hermitian eigensolver separates a cluster as well as
        float64 allows, which together with removing the residual-trace shift
        in double-double is what makes `DD_SWEEPS` one rather than many.

        Durand-Kerner was measured too, at one dd division per root against
        Aberth's five, and rejected for being non-monotone in the sweep
        count: 3.9e-16, then 9.7e-17, then 1.9e-16.
        """
        re_hi = np.empty((4, 4))
        re_lo = np.empty((4, 4))
        im_hi = np.empty((4, 4))
        im_lo = np.empty((4, 4))

        # Exact Hermitization.  A Hermitian diagonal is real, so its
        # imaginary limbs are exactly zero and never enter a product.
        for i in range(4):
            re_hi[i, i] = traceless[i, i].real
            re_lo[i, i] = 0.0
            im_hi[i, i] = 0.0
            im_lo[i, i] = 0.0
            for j in range(i+1, 4):
                s, e = _two_sum(traceless[i, j].real, traceless[j, i].real)
                re_hi[i, j] = 0.5*s
                re_lo[i, j] = 0.5*e
                re_hi[j, i] = re_hi[i, j]
                re_lo[j, i] = re_lo[i, j]
                s, e = _two_sum(traceless[i, j].imag, -traceless[j, i].imag)
                im_hi[i, j] = 0.5*s
                im_lo[i, j] = 0.5*e
                im_hi[j, i] = -im_hi[i, j]
                im_lo[j, i] = -im_lo[i, j]

        # The residual trace, removed exactly.  Dividing by four is exact.
        tau_hi, tau_lo = re_hi[0, 0], re_lo[0, 0]
        for i in range(1, 4):
            tau_hi, tau_lo = _dd_add(tau_hi, tau_lo, re_hi[i, i], re_lo[i, i])
        shift_hi, shift_lo = 0.25*tau_hi, 0.25*tau_lo
        for i in range(4):
            re_hi[i, i], re_lo[i, i] = _dd_sub(re_hi[i, i], re_lo[i, i],
                                               shift_hi, shift_lo)

        trace_2_hi, trace_2_lo = 0.0, 0.0
        for i in range(4):
            p_hi, p_lo = _dd_mul(re_hi[i, i], re_lo[i, i],
                                 re_hi[i, i], re_lo[i, i])
            trace_2_hi, trace_2_lo = _dd_add(trace_2_hi, trace_2_lo,
                                             p_hi, p_lo)
            for j in range(i+1, 4):
                p_hi, p_lo = _dd_mul(re_hi[i, j], re_lo[i, j],
                                     re_hi[i, j], re_lo[i, j])
                q_hi, q_lo = _dd_mul(im_hi[i, j], im_lo[i, j],
                                     im_hi[i, j], im_lo[i, j])
                p_hi, p_lo = _dd_add(p_hi, p_lo, q_hi, q_lo)
                trace_2_hi, trace_2_lo = _dd_add(trace_2_hi, trace_2_lo,
                                                 2.0*p_hi, 2.0*p_lo)

        sq_re_hi = np.zeros((4, 4))
        sq_re_lo = np.zeros((4, 4))
        sq_im_hi = np.zeros((4, 4))
        sq_im_lo = np.zeros((4, 4))
        for i in range(4):
            for j in range(i, 4):
                x_hi, x_lo = 0.0, 0.0
                y_hi, y_lo = 0.0, 0.0
                for k in range(4):
                    p_hi, p_lo = _dd_mul(re_hi[i, k], re_lo[i, k],
                                         re_hi[k, j], re_lo[k, j])
                    q_hi, q_lo = _dd_mul(im_hi[i, k], im_lo[i, k],
                                         im_hi[k, j], im_lo[k, j])
                    d_hi, d_lo = _dd_sub(p_hi, p_lo, q_hi, q_lo)
                    x_hi, x_lo = _dd_add(x_hi, x_lo, d_hi, d_lo)
                    p_hi, p_lo = _dd_mul(re_hi[i, k], re_lo[i, k],
                                         im_hi[k, j], im_lo[k, j])
                    q_hi, q_lo = _dd_mul(im_hi[i, k], im_lo[i, k],
                                         re_hi[k, j], re_lo[k, j])
                    d_hi, d_lo = _dd_add(p_hi, p_lo, q_hi, q_lo)
                    y_hi, y_lo = _dd_add(y_hi, y_lo, d_hi, d_lo)
                sq_re_hi[i, j], sq_re_lo[i, j] = x_hi, x_lo
                sq_im_hi[i, j], sq_im_lo[i, j] = y_hi, y_lo

        # Tr(H~^3) = sum_i S_ii H~_ii + 2 sum_{i<j} Re(S_ij conj(H~_ij))
        # Tr(H~^4) = sum_i S_ii^2   + 2 sum_{i<j} |S_ij|^2
        trace_3_hi, trace_3_lo = 0.0, 0.0
        trace_4_hi, trace_4_lo = 0.0, 0.0
        for i in range(4):
            p_hi, p_lo = _dd_mul(sq_re_hi[i, i], sq_re_lo[i, i],
                                 re_hi[i, i], re_lo[i, i])
            trace_3_hi, trace_3_lo = _dd_add(trace_3_hi, trace_3_lo,
                                             p_hi, p_lo)
            p_hi, p_lo = _dd_mul(sq_re_hi[i, i], sq_re_lo[i, i],
                                 sq_re_hi[i, i], sq_re_lo[i, i])
            trace_4_hi, trace_4_lo = _dd_add(trace_4_hi, trace_4_lo,
                                             p_hi, p_lo)
            for j in range(i+1, 4):
                p_hi, p_lo = _dd_mul(sq_re_hi[i, j], sq_re_lo[i, j],
                                     re_hi[i, j], re_lo[i, j])
                q_hi, q_lo = _dd_mul(sq_im_hi[i, j], sq_im_lo[i, j],
                                     im_hi[i, j], im_lo[i, j])
                p_hi, p_lo = _dd_add(p_hi, p_lo, q_hi, q_lo)
                trace_3_hi, trace_3_lo = _dd_add(trace_3_hi, trace_3_lo,
                                                 2.0*p_hi, 2.0*p_lo)
                p_hi, p_lo = _dd_mul(sq_re_hi[i, j], sq_re_lo[i, j],
                                     sq_re_hi[i, j], sq_re_lo[i, j])
                q_hi, q_lo = _dd_mul(sq_im_hi[i, j], sq_im_lo[i, j],
                                     sq_im_hi[i, j], sq_im_lo[i, j])
                p_hi, p_lo = _dd_add(p_hi, p_lo, q_hi, q_lo)
                trace_4_hi, trace_4_lo = _dd_add(trace_4_hi, trace_4_lo,
                                                 2.0*p_hi, 2.0*p_lo)

        invariant_2_hi, invariant_2_lo = 0.5*trace_2_hi, 0.5*trace_2_lo
        invariant_3_hi, invariant_3_lo = 0.5*trace_3_hi, 0.5*trace_3_lo
        p_hi, p_lo = _dd_mul(invariant_2_hi, invariant_2_lo,
                             invariant_2_hi, invariant_2_lo)
        d_hi, d_lo = _dd_sub(trace_4_hi, trace_4_lo, p_hi, p_lo)
        invariant_4_hi, invariant_4_lo = 0.5*d_hi, 0.5*d_lo

        # chi(psi) = psi^4 + c2 psi^2 + c1 psi + c0.  The two thirds is
        # formed by dd division: 2/3 is not representable in float64, and
        # lifting a rounded 2/3 into dd zeroes the low limb and poisons c1
        # back to float64 accuracy --- which is exactly the 1.1e-07 stall
        # this whole route exists to remove.
        c2_hi, c2_lo = -invariant_2_hi, -invariant_2_lo
        p_hi, p_lo = _dd_mul(2.0*invariant_3_hi, 2.0*invariant_3_lo,
                             THIRD_HI, THIRD_LO)
        c1_hi, c1_lo = -p_hi, -p_lo
        p_hi, p_lo = _dd_mul(invariant_2_hi, invariant_2_lo,
                             invariant_2_hi, invariant_2_lo)
        d_hi, d_lo = _dd_sub(p_hi, p_lo,
                             2.0*invariant_4_hi, 2.0*invariant_4_lo)
        c0_hi, c0_lo = 0.25*d_hi, 0.25*d_lo

        latent = np.linalg.eigvalsh(traceless)
        z_hi = np.empty(4)
        z_lo = np.empty(4)
        for k in range(4):
            # In double-double, not float64.  Subtracting the residual-trace
            # shift with one float64 operation discards the low limb the
            # exact traceless-ing above went to the trouble of computing,
            # and the start then arrives a whole ulp out --- which cost a
            # second Aberth sweep to recover, and was mistaken for the
            # iteration needing it.
            z_hi[k], z_lo[k] = _dd_sub(latent[k], 0.0, shift_hi, shift_lo)

        for _ in range(DD_SWEEPS):
            for i in range(4):
                sq_hi, sq_lo = _dd_mul(z_hi[i], z_lo[i], z_hi[i], z_lo[i])
                a_hi, a_lo = _dd_mul(sq_hi, sq_lo, sq_hi, sq_lo)
                b_hi, b_lo = _dd_mul(c2_hi, c2_lo, sq_hi, sq_lo)
                a_hi, a_lo = _dd_add(a_hi, a_lo, b_hi, b_lo)
                b_hi, b_lo = _dd_mul(c1_hi, c1_lo, z_hi[i], z_lo[i])
                b_hi, b_lo = _dd_add(b_hi, b_lo, c0_hi, c0_lo)
                chi_hi, chi_lo = _dd_add(a_hi, a_lo, b_hi, b_lo)

                # chi'(psi) = 4 psi^3 + 2 c2 psi + c1
                a_hi, a_lo = _dd_mul(sq_hi, sq_lo, z_hi[i], z_lo[i])
                b_hi, b_lo = _dd_mul(c2_hi, c2_lo, z_hi[i], z_lo[i])
                b_hi, b_lo = _dd_add(2.0*b_hi, 2.0*b_lo, c1_hi, c1_lo)
                der_hi, der_lo = _dd_add(4.0*a_hi, 4.0*a_lo, b_hi, b_lo)

                # Aberth's correction is Newton's step divided by one minus
                # the pull of the other three roots, which is what keeps a
                # cluster's members from converging onto each other and why
                # no half-gap guard is needed here.  Written as
                #
                #   step = (chi/chi')/(1 - (chi/chi') sum_j 1/(psi_i-psi_j))
                #
                # it takes five dd divisions per root.  Clearing the
                # denominators leaves the same expression as
                #
                #   step = chi P/(chi' P - chi S)
                #
                # for P and S the product and second elementary symmetric
                # function of the three gaps, and one dd division.  A dd
                # division is three float64 divisions where a dd multiply is
                # none, and this loop was over half the route's cost.
                gap_1_hi, gap_1_lo = 0.0, 0.0
                gap_2_hi, gap_2_lo = 0.0, 0.0
                gap_3_hi, gap_3_lo = 0.0, 0.0
                for j in range(4):
                    if j != i:
                        g_hi, g_lo = _dd_sub(z_hi[i], z_lo[i],
                                             z_hi[j], z_lo[j])
                        if j == (i + 1) % 4:
                            gap_1_hi, gap_1_lo = g_hi, g_lo
                        elif j == (i + 2) % 4:
                            gap_2_hi, gap_2_lo = g_hi, g_lo
                        else:
                            gap_3_hi, gap_3_lo = g_hi, g_lo

                pair_hi, pair_lo = _dd_mul(gap_2_hi, gap_2_lo,
                                           gap_3_hi, gap_3_lo)
                prod_hi, prod_lo = _dd_mul(gap_1_hi, gap_1_lo,
                                           pair_hi, pair_lo)
                s_hi, s_lo = _dd_add(gap_2_hi, gap_2_lo, gap_3_hi, gap_3_lo)
                s_hi, s_lo = _dd_mul(gap_1_hi, gap_1_lo, s_hi, s_lo)
                s_hi, s_lo = _dd_add(pair_hi, pair_lo, s_hi, s_lo)

                a_hi, a_lo = _dd_mul(der_hi, der_lo, prod_hi, prod_lo)
                b_hi, b_lo = _dd_mul(chi_hi, chi_lo, s_hi, s_lo)
                den_hi, den_lo = _dd_sub(a_hi, a_lo, b_hi, b_lo)
                if den_hi == 0.0:
                    continue

                # A vanishing gap leaves the step at exactly zero through P
                # rather than needing to be skipped: two roots that landed
                # on identical bits stay where the eigensolver put them.
                a_hi, a_lo = _dd_mul(chi_hi, chi_lo, prod_hi, prod_lo)
                step_hi, step_lo = _dd_div(a_hi, a_lo, den_hi, den_lo)
                z_hi[i], z_lo[i] = _dd_sub(z_hi[i], z_lo[i],
                                           step_hi, step_lo)

        # Only here do the two limbs collapse into one float64
        psi_0 = z_hi[0] + z_lo[0]
        psi_1 = z_hi[1] + z_lo[1]
        psi_2 = z_hi[2] + z_lo[2]
        psi_3 = z_hi[3] + z_lo[3]

        # Ascending, by the five-comparator network for four elements.  The
        # eigensolver's order is ascending and Aberth moves each root by far
        # less than a gap, but `_divided_differences` below is written for
        # sorted roots and the sort is four comparisons.
        if psi_0 > psi_1:
            psi_0, psi_1 = psi_1, psi_0
        if psi_2 > psi_3:
            psi_2, psi_3 = psi_3, psi_2
        if psi_0 > psi_2:
            psi_0, psi_2 = psi_2, psi_0
        if psi_1 > psi_3:
            psi_1, psi_3 = psi_3, psi_1
        if psi_1 > psi_2:
            psi_1, psi_2 = psi_2, psi_1
        return psi_0, psi_1, psi_2, psi_3

    @njit(cache=True, inline='always')
    def _operator_4nu(h_matrix, L, strategy, work):
        r"""Builds :math:`U_4(L)` for one Hamiltonian in ``work[1]``.

        Factored out of `_one_4nu`, which computed the operator and
        then squared it, so that the probability kernel and the
        evolution-operator kernel share the expensive part.  Unlike
        three flavors the operator is already materialised here, in
        the caller's scratch, so the split costs nothing at all.

        A transcription of :func:`oscprob4nu._evolution_operator_4nu_array`
        for a single element: the traceless part, its four latent roots by
        whichever route `strategy` names, the divided differences of the
        exponential over them, and the Newton-form reconstruction of
        :math:`U_4`.

        ``work`` is scratch space of shape ``(5, 4, 4)``, supplied by the
        caller so that the loop over a stack allocates nothing.
        `_latent_roots_dd` breaks that rule: it allocates eight small real
        arrays for its limbs rather than taking them from `work`, which
        would mean widening the scratch that all eleven call sites pass for
        the sake of a route one of them may not take.  Whether Numba hoists
        those allocations out of the parallel loop has *not* been measured,
        and if the double-double route's cost ever wants attacking this is
        the first place to look.
        """
        traceless = work[0]
        operator = work[1]
        first = work[2]
        second = work[3]
        shifted = work[4]

        trace = (h_matrix[0, 0] + h_matrix[1, 1]
                 + h_matrix[2, 2] + h_matrix[3, 3]).real/4.0
        for i in range(4):
            for j in range(4):
                traceless[i, j] = h_matrix[i, j]
            traceless[i, i] = h_matrix[i, i] - trace

        # `second` no longer holds H~^2 on the way in; it is scratch for the
        # reconstruction below, which writes every entry before reading any.
        if strategy == 0:
            psi_0, psi_1, psi_2, psi_3 = _latent_roots_dd(traceless)
        else:
            # Eigenvalues of the matrix, not roots of its invariants: in
            # float64 those three numbers no longer hold a stiff spectrum's
            # cluster apart, while a Hermitian eigensolver never forms them.
            # LAPACK returns them ascending, so the five-comparator sorting
            # network the closed form needed goes with it.
            latent = np.linalg.eigvalsh(traceless)
            psi_0 = latent[0]
            psi_1 = latent[1]
            psi_2 = latent[2]
            psi_3 = latent[3]

        if strategy == 1:
            # One Newton step on chi, with chi'(psi_m) taken as the
            # product of the gaps to the other three roots, and refused
            # wherever it would carry a root more than halfway to its
            # nearest neighbour --- see oscprob4nu._polish_roots, whose
            # guard this is
            gap_01 = psi_0 - psi_1
            gap_02 = psi_0 - psi_2
            gap_03 = psi_0 - psi_3
            gap_12 = psi_1 - psi_2
            gap_13 = psi_1 - psi_3
            gap_23 = psi_2 - psi_3

            near_0 = min(abs(gap_01), abs(gap_02), abs(gap_03))
            near_1 = min(abs(gap_01), abs(gap_12), abs(gap_13))
            near_2 = min(abs(gap_02), abs(gap_12), abs(gap_23))
            near_3 = min(abs(gap_03), abs(gap_13), abs(gap_23))

            derivative = gap_01*gap_02*gap_03
            if derivative != 0.0:
                step = _chi_4nu(traceless, psi_0, shifted)/derivative
                if abs(step) <= 0.5*near_0:
                    psi_0 -= step
            derivative = -gap_01*gap_12*gap_13
            if derivative != 0.0:
                step = _chi_4nu(traceless, psi_1, shifted)/derivative
                if abs(step) <= 0.5*near_1:
                    psi_1 -= step
            derivative = gap_02*gap_12*gap_23
            if derivative != 0.0:
                step = _chi_4nu(traceless, psi_2, shifted)/derivative
                if abs(step) <= 0.5*near_2:
                    psi_2 -= step
            derivative = -gap_03*gap_13*gap_23
            if derivative != 0.0:
                step = _chi_4nu(traceless, psi_3, shifted)/derivative
                if abs(step) <= 0.5*near_3:
                    psi_3 -= step

            if psi_0 > psi_1:
                psi_0, psi_1 = psi_1, psi_0
            if psi_2 > psi_3:
                psi_2, psi_3 = psi_3, psi_2
            if psi_0 > psi_2:
                psi_0, psi_2 = psi_2, psi_0
            if psi_1 > psi_3:
                psi_1, psi_3 = psi_3, psi_1
            if psi_1 > psi_2:
                psi_1, psi_2 = psi_2, psi_1

        # Divided differences of exp(-i psi L) over the four roots,
        # taking the confluent value wherever two nodes have merged
        spectral = abs(psi_0)
        if abs(psi_1) > spectral:
            spectral = abs(psi_1)
        if abs(psi_2) > spectral:
            spectral = abs(psi_2)
        if abs(psi_3) > spectral:
            spectral = abs(psi_3)
        tolerance = DEGENERACY_TOL*(spectral if spectral > 0.0 else 1.0)

        phase_0 = cmath.rect(1.0, -psi_0*L)
        phase_1 = cmath.rect(1.0, -psi_1*L)
        phase_2 = cmath.rect(1.0, -psi_2*L)
        phase_3 = cmath.rect(1.0, -psi_3*L)
        minus_i_l = complex(0.0, -L)

        table_0 = phase_0
        table_1 = phase_1
        table_2 = phase_2
        table_3 = phase_3
        coeff_0th = table_0

        weight = minus_i_l
        separation = psi_1 - psi_0
        if abs(separation) > tolerance:
            new_0 = (table_1 - table_0)/separation
        else:
            new_0 = weight*phase_0
        separation = psi_2 - psi_1
        if abs(separation) > tolerance:
            new_1 = (table_2 - table_1)/separation
        else:
            new_1 = weight*phase_1
        separation = psi_3 - psi_2
        if abs(separation) > tolerance:
            new_2 = (table_3 - table_2)/separation
        else:
            new_2 = weight*phase_2
        table_0, table_1, table_2 = new_0, new_1, new_2
        coeff_1st = table_0

        weight = minus_i_l*minus_i_l/2.0
        separation = psi_2 - psi_0
        if abs(separation) > tolerance:
            new_0 = (table_1 - table_0)/separation
        else:
            new_0 = weight*phase_0
        separation = psi_3 - psi_1
        if abs(separation) > tolerance:
            new_1 = (table_2 - table_1)/separation
        else:
            new_1 = weight*phase_1
        table_0, table_1 = new_0, new_1
        coeff_2nd = table_0

        weight = minus_i_l*minus_i_l*minus_i_l/6.0
        separation = psi_3 - psi_0
        if abs(separation) > tolerance:
            coeff_3rd = (table_1 - table_0)/separation
        else:
            coeff_3rd = weight*phase_0

        # U_4 = c_0 + c_1 (H~ - psi_0) + c_2 (H~ - psi_0)(H~ - psi_1)
        #           + c_3 (H~ - psi_0)(H~ - psi_1)(H~ - psi_2)
        for i in range(4):
            for j in range(4):
                first[i, j] = traceless[i, j]
            first[i, i] = traceless[i, i] - psi_0

        for i in range(4):
            for j in range(4):
                operator[i, j] = coeff_1st*first[i, j]
            operator[i, i] += coeff_0th

        for i in range(4):
            for j in range(4):
                shifted[i, j] = traceless[i, j]
            shifted[i, i] = traceless[i, i] - psi_1

        for i in range(4):
            for j in range(4):
                entry = 0.0j
                for k in range(4):
                    entry += first[i, k]*shifted[k, j]
                second[i, j] = entry
                operator[i, j] += coeff_2nd*entry

        for i in range(4):
            for j in range(4):
                shifted[i, j] = traceless[i, j]
            shifted[i, i] = traceless[i, i] - psi_2

        for i in range(4):
            for j in range(4):
                entry = 0.0j
                for k in range(4):
                    entry += second[i, k]*shifted[k, j]
                operator[i, j] += coeff_3rd*entry


    @njit(cache=True, inline='always')
    def _one_4nu(h_matrix, L, out, n, strategy, work):
        r"""Writes the sixteen probabilities into ``out[n]``."""
        _operator_4nu(h_matrix, L, strategy, work)
        operator = work[1]

        # P_ab = |U_ba|^2, initial flavor slowest
        for alpha in range(4):
            for beta in range(4):
                entry = operator[beta, alpha]
                out[n, 4*alpha + beta] = (entry.real*entry.real
                                          + entry.imag*entry.imag)

    @njit(cache=True, inline='always')
    def _one_4nu_u(h_matrix, L, out, n, strategy, work):
        r"""Writes the sixteen entries of :math:`U_4(L)` into ``out[n]``.

        Row-major and indexed ``(final, initial)``, so reshaping to
        ``(4, 4)`` gives the matrix `oscprob4nu` returns --- not the
        flavor order the probabilities use, which runs the initial index
        slowest.  The two differ by a transpose.
        """
        _operator_4nu(h_matrix, L, strategy, work)
        operator = work[1]

        for i in range(4):
            for j in range(4):
                out[n, 4*i + j] = operator[i, j]

    @njit(cache=True)
    def _run_3nu_serial(h_stack, l_stack, out):
        for n in range(h_stack.shape[0]):
            _one_3nu(h_stack[n], l_stack[n], out, n)

    @njit(cache=True, parallel=True)
    def _run_3nu_parallel(h_stack, l_stack, out):
        for n in prange(h_stack.shape[0]):
            _one_3nu(h_stack[n], l_stack[n], out, n)

    @njit(cache=True, inline='always')
    def _entries_2nu(h_matrix, L):
        r"""Returns the four entries of :math:`U_2(L)`, row by row.

        Unlike three and four flavors this is *not* a refactor of the
        probability kernel.  `_one_2nu` never forms the operator: at two
        flavors :math:`P_{e\mu}` has a closed form in the coefficients
        alone, so the kernel goes straight to it and the survival
        probability is its complement.  That shortcut is worth keeping,
        so the operator is written here separately rather than factored
        out of it.

        The expansion is :math:`U_2 = u_0 \mathbb{1} + i u_k \sigma^k`
        with :math:`u_0 = \cos(|h|L)` and
        :math:`u_k = -h_k \sin(|h|L)/|h|`, transcribed from the scalar
        path in :mod:`oscprob2nu`.  When :math:`|h| = 0` the Hamiltonian
        is proportional to the identity and the limit
        :math:`\sin(|h|L)/|h| \to L` is taken, which is the difference
        between an exact identity and a NaN.
        """
        h0 = h_matrix[0, 1].real
        h1 = -h_matrix[0, 1].imag
        h2 = (h_matrix[0, 0] - h_matrix[1, 1]).real/2.0

        hsq = h0*h0 + h1*h1 + h2*h2
        h_abs = math.sqrt(hsq)
        u0 = math.cos(h_abs*L)
        if h_abs == 0.0:
            ss = -L
        else:
            ss = -math.sin(h_abs*L)/h_abs

        u1 = h0*ss
        u2 = h1*ss
        u3 = h2*ss

        return (u0 + 1.0j*u3, 1.0j*u1 + u2,
                1.0j*u1 - u2, u0 - 1.0j*u3)

    @njit(cache=True, inline='always')
    def _one_2nu_u(h_matrix, L, out, n):
        r"""Writes the four entries of :math:`U_2(L)` into ``out[n]``."""
        u_ee, u_em, u_me, u_mm = _entries_2nu(h_matrix, L)
        out[n, 0] = u_ee
        out[n, 1] = u_em
        out[n, 2] = u_me
        out[n, 3] = u_mm

    @njit(cache=True)
    def _run_2nu_u_serial(h_stack, l_stack, out):
        for n in range(l_stack.shape[0]):
            _one_2nu_u(h_stack[n], l_stack[n], out, n)

    @njit(cache=True, parallel=True)
    def _run_2nu_u_parallel(h_stack, l_stack, out):
        for n in prange(l_stack.shape[0]):
            _one_2nu_u(h_stack[n], l_stack[n], out, n)

    @njit(cache=True)
    def _slab_product_2nu(h_stack, widths, out):
        r"""Multiplies the per-slab two-flavor operators into ``out``.

        ``U = U_n ... U_1``, first slab crossed rightmost, as at three
        and four flavors.
        """
        a00, a01, a10, a11 = _entries_2nu(h_stack[0], widths[0])

        for k in range(1, widths.shape[0]):
            u00, u01, u10, u11 = _entries_2nu(h_stack[k], widths[k])
            b00 = u00*a00 + u01*a10
            b01 = u00*a01 + u01*a11
            b10 = u10*a00 + u11*a10
            b11 = u10*a01 + u11*a11
            a00 = b00
            a01 = b01
            a10 = b10
            a11 = b11

        out[0, 0] = a00
        out[0, 1] = a01
        out[1, 0] = a10
        out[1, 1] = a11

    @njit(cache=True)
    def _run_3nu_u_serial(h_stack, l_stack, out):
        for n in range(l_stack.shape[0]):
            _one_3nu_u(h_stack[n], l_stack[n], out, n)

    @njit(cache=True, parallel=True)
    def _run_3nu_u_parallel(h_stack, l_stack, out):
        for n in prange(l_stack.shape[0]):
            _one_3nu_u(h_stack[n], l_stack[n], out, n)

    @njit(cache=True)
    def _slab_product_3nu(h_stack, widths, out):
        r"""Multiplies the per-slab operators into ``out``, in order.

        Sequential by nature --- the product does not commute --- so
        there is no parallel counterpart.  The win is not threads: it is
        that the 120-odd operators of an Earth crossing are never
        materialised, and the 119 matrix products never leave registers,
        where the NumPy path allocated a stack of them and then made one
        dispatched call per multiplication.

        ``U = U_n ... U_2 U_1``: the slab crossed first is applied first
        and so stands rightmost, because the operator is indexed
        ``(final, initial)`` and acts to the left on the initial state.
        Accumulating ``acc <- U_k @ acc`` in increasing k is exactly that
        order, and getting it backwards is the classic way to be wrong
        here by something that still looks like a probability.
        """
        acc = np.empty((3, 3), dtype=np.complex128)
        tmp = np.empty((3, 3), dtype=np.complex128)

        (u_ee, u_em, u_et,
         u_me, u_mm, u_mt,
         u_te, u_tm, u_tt) = _entries_3nu(h_stack[0], widths[0])
        acc[0, 0] = u_ee
        acc[0, 1] = u_em
        acc[0, 2] = u_et
        acc[1, 0] = u_me
        acc[1, 1] = u_mm
        acc[1, 2] = u_mt
        acc[2, 0] = u_te
        acc[2, 1] = u_tm
        acc[2, 2] = u_tt

        for k in range(1, widths.shape[0]):
            (u_ee, u_em, u_et,
             u_me, u_mm, u_mt,
             u_te, u_tm, u_tt) = _entries_3nu(h_stack[k], widths[k])
            for j in range(3):
                a0 = acc[0, j]
                a1 = acc[1, j]
                a2 = acc[2, j]
                tmp[0, j] = u_ee*a0 + u_em*a1 + u_et*a2
                tmp[1, j] = u_me*a0 + u_mm*a1 + u_mt*a2
                tmp[2, j] = u_te*a0 + u_tm*a1 + u_tt*a2
            for i in range(3):
                for j in range(3):
                    acc[i, j] = tmp[i, j]

        for i in range(3):
            for j in range(3):
                out[i, j] = acc[i, j]

    @njit(cache=True)
    def _slab_product_2nu_batch_serial(h_stack, widths, out):
        for c in range(h_stack.shape[0]):
            _slab_product_2nu(h_stack[c], widths, out[c])

    @njit(cache=True, parallel=True)
    def _slab_product_2nu_batch_parallel(h_stack, widths, out):
        for c in prange(h_stack.shape[0]):
            _slab_product_2nu(h_stack[c], widths, out[c])

    @njit(cache=True)
    def _slab_product_3nu_batch_serial(h_stack, widths, out):
        r"""Composes one chord per leading index, sequentially.

        The product *along* a chord does not commute and so stays
        serial, but the chords themselves are independent: an energy
        scan at fixed zenith angle is the same geometry evaluated at
        many energies, and nothing couples one energy to another.  That
        is the axis the parallel counterpart spreads over, and it is why
        batching buys threads that the per-chord kernel could not use.
        """
        for c in range(h_stack.shape[0]):
            _slab_product_3nu(h_stack[c], widths, out[c])

    @njit(cache=True, parallel=True)
    def _slab_product_3nu_batch_parallel(h_stack, widths, out):
        for c in prange(h_stack.shape[0]):
            _slab_product_3nu(h_stack[c], widths, out[c])

    @njit(cache=True)
    def _slab_product_4nu_batch_serial(h_stack, widths, strategy, out):
        for c in range(h_stack.shape[0]):
            _slab_product_4nu(h_stack[c], widths, strategy, out[c])

    @njit(cache=True, parallel=True)
    def _slab_product_4nu_batch_parallel(h_stack, widths, strategy, out):
        for c in prange(h_stack.shape[0]):
            _slab_product_4nu(h_stack[c], widths, strategy, out[c])

    @njit(cache=True)
    def _build_h_2nu(h_vac, inv_e, potentials, k, h_work):
        r"""Builds slab ``k``'s two-flavor Hamiltonian into ``h_work``."""
        for i in range(2):
            for j in range(2):
                h_work[i, j] = h_vac[i, j]*inv_e
        h_work[0, 0] += potentials[k]

    @njit(cache=True)
    def _build_h_3nu(h_vac, inv_e, potentials, k, h_work):
        r"""Builds slab ``k``'s three-flavor Hamiltonian into ``h_work``.

        The whole point of the fused path, and it computes exactly what
        `hamiltonians3nu.hamiltonian_3nu_matter` computes --- this is a
        compiled mirror of that, not a different scheme.
        ``H = H_vac/E + V P_ee``, so the only thing that varies from
        slab to slab is one real number, and materialising a stack of
        3x3 matrices to carry it was what made the batched scan
        memory-bound.
        """
        for i in range(3):
            for j in range(3):
                h_work[i, j] = h_vac[i, j]*inv_e
        h_work[0, 0] += potentials[k]

    @njit(cache=True)
    def _build_h_4nu(h_vac, inv_e, potentials, potentials_nc, k, h_work):
        r"""Builds slab ``k``'s four-flavor Hamiltonian into ``h_work``.

        Two potentials rather than one: a sterile state does not feel
        the neutral current, so ``V_NC`` no longer cancels between the
        flavors.  The mirror of
        `hamiltonians4nu.hamiltonian_4nu_matter`.

        The result is made traceless here, because that is what
        `_operator_4nu` expects; the dropped phase is per slab and
        cancels in every probability.
        """
        for i in range(4):
            for j in range(4):
                h_work[i, j] = h_vac[i, j]*inv_e
        h_work[0, 0] += potentials[k]
        h_work[3, 3] -= potentials_nc[k]

        trace = (h_work[0, 0] + h_work[1, 1]
                 + h_work[2, 2] + h_work[3, 3])/4.0
        for i in range(4):
            h_work[i, i] -= trace

    @njit(cache=True)
    def _earth_chord_2nu(h_vac, inv_e, potentials, widths, out):
        r"""Composes one two-flavor chord, Hamiltonians built inline."""
        h_work = np.empty((2, 2), dtype=np.complex128)

        _build_h_2nu(h_vac, inv_e, potentials, 0, h_work)
        a00, a01, a10, a11 = _entries_2nu(h_work, widths[0])

        for k in range(1, widths.shape[0]):
            _build_h_2nu(h_vac, inv_e, potentials, k, h_work)
            u00, u01, u10, u11 = _entries_2nu(h_work, widths[k])
            b00 = u00*a00 + u01*a10
            b01 = u00*a01 + u01*a11
            b10 = u10*a00 + u11*a10
            b11 = u10*a01 + u11*a11
            a00 = b00
            a01 = b01
            a10 = b10
            a11 = b11

        out[0, 0] = a00
        out[0, 1] = a01
        out[1, 0] = a10
        out[1, 1] = a11

    @njit(cache=True)
    def _earth_chord_3nu(h_vac, inv_e, potentials, widths, out):
        r"""Composes one three-flavor chord, Hamiltonians built inline.

        ``U = U_n ... U_1``, the slab crossed first applied first and so
        standing rightmost, exactly as `_slab_product_3nu` orders it.
        """
        acc = np.empty((3, 3), dtype=np.complex128)
        tmp = np.empty((3, 3), dtype=np.complex128)
        h_work = np.empty((3, 3), dtype=np.complex128)

        _build_h_3nu(h_vac, inv_e, potentials, 0, h_work)
        (u_ee, u_em, u_et,
         u_me, u_mm, u_mt,
         u_te, u_tm, u_tt) = _entries_3nu(h_work, widths[0])
        acc[0, 0] = u_ee
        acc[0, 1] = u_em
        acc[0, 2] = u_et
        acc[1, 0] = u_me
        acc[1, 1] = u_mm
        acc[1, 2] = u_mt
        acc[2, 0] = u_te
        acc[2, 1] = u_tm
        acc[2, 2] = u_tt

        for k in range(1, widths.shape[0]):
            _build_h_3nu(h_vac, inv_e, potentials, k, h_work)
            (u_ee, u_em, u_et,
             u_me, u_mm, u_mt,
             u_te, u_tm, u_tt) = _entries_3nu(h_work, widths[k])
            for j in range(3):
                a0 = acc[0, j]
                a1 = acc[1, j]
                a2 = acc[2, j]
                tmp[0, j] = u_ee*a0 + u_em*a1 + u_et*a2
                tmp[1, j] = u_me*a0 + u_mm*a1 + u_mt*a2
                tmp[2, j] = u_te*a0 + u_tm*a1 + u_tt*a2
            for i in range(3):
                for j in range(3):
                    acc[i, j] = tmp[i, j]

        for i in range(3):
            for j in range(3):
                out[i, j] = acc[i, j]

    @njit(cache=True)
    def _earth_chord_4nu(h_vac, inv_e, potentials, potentials_nc, widths,
                         strategy, out):
        r"""Composes one four-flavor chord, Hamiltonians built inline."""
        work = np.empty((5, 4, 4), dtype=np.complex128)
        acc = np.empty((4, 4), dtype=np.complex128)
        tmp = np.empty((4, 4), dtype=np.complex128)
        h_work = np.empty((4, 4), dtype=np.complex128)

        _build_h_4nu(h_vac, inv_e, potentials, potentials_nc, 0, h_work)
        _operator_4nu(h_work, widths[0], strategy, work)
        for i in range(4):
            for j in range(4):
                acc[i, j] = work[1][i, j]

        for k in range(1, widths.shape[0]):
            _build_h_4nu(h_vac, inv_e, potentials, potentials_nc, k, h_work)
            _operator_4nu(h_work, widths[k], strategy, work)
            for i in range(4):
                for j in range(4):
                    total = 0.0 + 0.0j
                    for m in range(4):
                        total += work[1][i, m]*acc[m, j]
                    tmp[i, j] = total
            for i in range(4):
                for j in range(4):
                    acc[i, j] = tmp[i, j]

        for i in range(4):
            for j in range(4):
                out[i, j] = acc[i, j]

    @njit(cache=True)
    def _palindromic_stack(h_stack, widths):
        r"""Returns whether a supplied slab sequence is a palindrome.

        The compiled counterpart of `palindromic`, and it exists because
        the NumPy one is too slow to be worth calling here.  Comparing a
        materialised ``(n, 3, 3)`` complex stack against its reverse
        costs about 6 microseconds through NumPy --- reversed views,
        temporaries, and a full pass whatever the answer --- against the
        6 microseconds the halved composition saves on a 120-slab chord.
        Measured, that check turned a 1.55x win into a 0.98x loss.

        This walks only the first half, compares each slab against its
        mirror, and returns at the first disagreement, which is the
        common case for a sequence that is not symmetric.
        """
        n = widths.shape[0]
        for k in range(n//2):
            j = n - 1 - k
            if widths[k] != widths[j]:
                return False
            for a in range(h_stack.shape[1]):
                for b in range(h_stack.shape[2]):
                    if h_stack[k, a, b] != h_stack[j, a, b]:
                        return False

        return True

    @njit(cache=True)
    def _slab_product_2nu_mirrored(h_stack, widths, out):
        r"""Composes a palindromic two-flavor slab sequence at half cost."""
        n = widths.shape[0]
        m = n//2
        a00 = 1.0 + 0.0j
        a01 = 0.0 + 0.0j
        a10 = 0.0 + 0.0j
        a11 = 1.0 + 0.0j
        b00 = 1.0 + 0.0j
        b01 = 0.0 + 0.0j
        b10 = 0.0 + 0.0j
        b11 = 1.0 + 0.0j

        for k in range(m):
            u00, u01, u10, u11 = _entries_2nu(h_stack[k], widths[k])
            t00 = u00*b00 + u01*b10
            t01 = u00*b01 + u01*b11
            t10 = u10*b00 + u11*b10
            t11 = u10*b01 + u11*b11
            b00 = t00
            b01 = t01
            b10 = t10
            b11 = t11
            t00 = a00*u00 + a01*u10
            t01 = a00*u01 + a01*u11
            t10 = a10*u00 + a11*u10
            t11 = a10*u01 + a11*u11
            a00 = t00
            a01 = t01
            a10 = t10
            a11 = t11

        if n % 2 == 1:
            u00, u01, u10, u11 = _entries_2nu(h_stack[m], widths[m])
            t00 = u00*b00 + u01*b10
            t01 = u00*b01 + u01*b11
            t10 = u10*b00 + u11*b10
            t11 = u10*b01 + u11*b11
            b00 = t00
            b01 = t01
            b10 = t10
            b11 = t11

        out[0, 0] = a00*b00 + a01*b10
        out[0, 1] = a00*b01 + a01*b11
        out[1, 0] = a10*b00 + a11*b10
        out[1, 1] = a10*b01 + a11*b11

    @njit(cache=True)
    def _slab_product_3nu_mirrored(h_stack, widths, out):
        r"""Composes a palindromic three-flavor slab sequence at half cost.

        The counterpart of `_earth_chord_3nu_mirrored` for a sequence
        whose Hamiltonians are supplied rather than built, which is what
        `slabs` composes: a single Earth crossing, a symmetric castle
        wall, any profile that reads the same from either end.
        """
        acc_a = np.empty((3, 3), dtype=np.complex128)
        acc_b = np.empty((3, 3), dtype=np.complex128)
        tmp_a = np.empty((3, 3), dtype=np.complex128)
        tmp_b = np.empty((3, 3), dtype=np.complex128)
        u = np.empty((3, 3), dtype=np.complex128)

        for i in range(3):
            for j in range(3):
                acc_a[i, j] = 1.0 + 0.0j if i == j else 0.0 + 0.0j
                acc_b[i, j] = 1.0 + 0.0j if i == j else 0.0 + 0.0j

        n = widths.shape[0]
        m = n//2
        for k in range(m):
            (u[0, 0], u[0, 1], u[0, 2],
             u[1, 0], u[1, 1], u[1, 2],
             u[2, 0], u[2, 1], u[2, 2]) = _entries_3nu(h_stack[k], widths[k])
            for i in range(3):
                for j in range(3):
                    sb = 0.0 + 0.0j
                    sa = 0.0 + 0.0j
                    for t in range(3):
                        sb += u[i, t]*acc_b[t, j]
                        sa += acc_a[i, t]*u[t, j]
                    tmp_b[i, j] = sb
                    tmp_a[i, j] = sa
            for i in range(3):
                for j in range(3):
                    acc_b[i, j] = tmp_b[i, j]
                    acc_a[i, j] = tmp_a[i, j]

        if n % 2 == 1:
            (u[0, 0], u[0, 1], u[0, 2],
             u[1, 0], u[1, 1], u[1, 2],
             u[2, 0], u[2, 1], u[2, 2]) = _entries_3nu(h_stack[m], widths[m])
            for i in range(3):
                for j in range(3):
                    s = 0.0 + 0.0j
                    for t in range(3):
                        s += u[i, t]*acc_b[t, j]
                    tmp_b[i, j] = s
            for i in range(3):
                for j in range(3):
                    acc_b[i, j] = tmp_b[i, j]

        for i in range(3):
            for j in range(3):
                s = 0.0 + 0.0j
                for t in range(3):
                    s += acc_a[i, t]*acc_b[t, j]
                out[i, j] = s

    @njit(cache=True)
    def _slab_product_4nu_mirrored(h_stack, widths, strategy, out):
        r"""Composes a palindromic four-flavor slab sequence at half cost."""
        work = np.empty((5, 4, 4), dtype=np.complex128)
        acc_a = np.empty((4, 4), dtype=np.complex128)
        acc_b = np.empty((4, 4), dtype=np.complex128)
        tmp_a = np.empty((4, 4), dtype=np.complex128)
        tmp_b = np.empty((4, 4), dtype=np.complex128)

        for i in range(4):
            for j in range(4):
                acc_a[i, j] = 1.0 + 0.0j if i == j else 0.0 + 0.0j
                acc_b[i, j] = 1.0 + 0.0j if i == j else 0.0 + 0.0j

        n = widths.shape[0]
        m = n//2
        for k in range(m):
            _operator_4nu(h_stack[k], widths[k], strategy, work)
            for i in range(4):
                for j in range(4):
                    sb = 0.0 + 0.0j
                    sa = 0.0 + 0.0j
                    for t in range(4):
                        sb += work[1][i, t]*acc_b[t, j]
                        sa += acc_a[i, t]*work[1][t, j]
                    tmp_b[i, j] = sb
                    tmp_a[i, j] = sa
            for i in range(4):
                for j in range(4):
                    acc_b[i, j] = tmp_b[i, j]
                    acc_a[i, j] = tmp_a[i, j]

        if n % 2 == 1:
            _operator_4nu(h_stack[m], widths[m], strategy, work)
            for i in range(4):
                for j in range(4):
                    s = 0.0 + 0.0j
                    for t in range(4):
                        s += work[1][i, t]*acc_b[t, j]
                    tmp_b[i, j] = s
            for i in range(4):
                for j in range(4):
                    acc_b[i, j] = tmp_b[i, j]

        for i in range(4):
            for j in range(4):
                s = 0.0 + 0.0j
                for t in range(4):
                    s += acc_a[i, t]*acc_b[t, j]
                out[i, j] = s

    @njit(cache=True)
    def _earth_chord_2nu_mirrored(h_vac, inv_e, potentials, widths, out):
        r"""Composes one two-flavor chord, using its palindrome."""
        n = widths.shape[0]
        m = n//2
        h_work = np.empty((2, 2), dtype=np.complex128)
        a00 = 1.0 + 0.0j
        a01 = 0.0 + 0.0j
        a10 = 0.0 + 0.0j
        a11 = 1.0 + 0.0j
        b00 = 1.0 + 0.0j
        b01 = 0.0 + 0.0j
        b10 = 0.0 + 0.0j
        b11 = 1.0 + 0.0j

        for k in range(m):
            _build_h_2nu(h_vac, inv_e, potentials, k, h_work)
            u00, u01, u10, u11 = _entries_2nu(h_work, widths[k])
            # B <- U B
            t00 = u00*b00 + u01*b10
            t01 = u00*b01 + u01*b11
            t10 = u10*b00 + u11*b10
            t11 = u10*b01 + u11*b11
            b00 = t00
            b01 = t01
            b10 = t10
            b11 = t11
            # A <- A U
            t00 = a00*u00 + a01*u10
            t01 = a00*u01 + a01*u11
            t10 = a10*u00 + a11*u10
            t11 = a10*u01 + a11*u11
            a00 = t00
            a01 = t01
            a10 = t10
            a11 = t11

        if n % 2 == 1:
            _build_h_2nu(h_vac, inv_e, potentials, m, h_work)
            u00, u01, u10, u11 = _entries_2nu(h_work, widths[m])
            t00 = u00*b00 + u01*b10
            t01 = u00*b01 + u01*b11
            t10 = u10*b00 + u11*b10
            t11 = u10*b01 + u11*b11
            b00 = t00
            b01 = t01
            b10 = t10
            b11 = t11

        out[0, 0] = a00*b00 + a01*b10
        out[0, 1] = a00*b01 + a01*b11
        out[1, 0] = a10*b00 + a11*b10
        out[1, 1] = a10*b01 + a11*b11

    @njit(cache=True)
    def _earth_chord_3nu_mirrored(h_vac, inv_e, potentials, widths, out):
        r"""Composes one three-flavor chord, using its palindrome.

        A chord through a spherically symmetric Earth meets every radius
        twice, so slab ``j`` and slab ``n-1-j`` carry the same
        Hamiltonian and the same width and therefore the *same*
        operator.  Half the SU(3) expansions in a chord are recomputing
        the other half.

        Writing ``U = U_{n-1} ... U_0`` and splitting at the centre,

        ``U = (U_0 U_1 ... U_{m-1}) (U_{m-1} ... U_0) = A B``

        for even ``n = 2m``, and ``U = A U_m B`` for odd ``n = 2m+1``.
        Both products accumulate in one pass over the first half, so each
        expansion is computed once and used twice; only the matrix
        products still number ``n``.  The expansion is about two thirds
        of a slab's cost, which is why this is worth roughly 1.5x rather
        than 2x.

        The caller decides whether the chord is a palindrome, by exact
        equality --- see `worthwhile_mirror`.  Nothing here checks, so
        nothing here may be called on a chord that is not one.
        """
        acc_a = np.empty((3, 3), dtype=np.complex128)
        acc_b = np.empty((3, 3), dtype=np.complex128)
        tmp_a = np.empty((3, 3), dtype=np.complex128)
        tmp_b = np.empty((3, 3), dtype=np.complex128)
        h_work = np.empty((3, 3), dtype=np.complex128)
        u = np.empty((3, 3), dtype=np.complex128)

        for i in range(3):
            for j in range(3):
                acc_a[i, j] = 1.0 + 0.0j if i == j else 0.0 + 0.0j
                acc_b[i, j] = 1.0 + 0.0j if i == j else 0.0 + 0.0j

        n = widths.shape[0]
        m = n//2
        for k in range(m):
            _build_h_3nu(h_vac, inv_e, potentials, k, h_work)
            (u[0, 0], u[0, 1], u[0, 2],
             u[1, 0], u[1, 1], u[1, 2],
             u[2, 0], u[2, 1], u[2, 2]) = _entries_3nu(h_work, widths[k])
            for i in range(3):
                for j in range(3):
                    sb = 0.0 + 0.0j
                    sa = 0.0 + 0.0j
                    for t in range(3):
                        sb += u[i, t]*acc_b[t, j]
                        sa += acc_a[i, t]*u[t, j]
                    tmp_b[i, j] = sb
                    tmp_a[i, j] = sa
            for i in range(3):
                for j in range(3):
                    acc_b[i, j] = tmp_b[i, j]
                    acc_a[i, j] = tmp_a[i, j]

        if n % 2 == 1:
            # The middle slab is its own mirror and so is applied once
            _build_h_3nu(h_vac, inv_e, potentials, m, h_work)
            (u[0, 0], u[0, 1], u[0, 2],
             u[1, 0], u[1, 1], u[1, 2],
             u[2, 0], u[2, 1], u[2, 2]) = _entries_3nu(h_work, widths[m])
            for i in range(3):
                for j in range(3):
                    s = 0.0 + 0.0j
                    for t in range(3):
                        s += u[i, t]*acc_b[t, j]
                    tmp_b[i, j] = s
            for i in range(3):
                for j in range(3):
                    acc_b[i, j] = tmp_b[i, j]

        for i in range(3):
            for j in range(3):
                s = 0.0 + 0.0j
                for t in range(3):
                    s += acc_a[i, t]*acc_b[t, j]
                out[i, j] = s

    @njit(cache=True)
    def _earth_chord_4nu_mirrored(h_vac, inv_e, potentials, potentials_nc,
                                  widths, strategy, out):
        r"""Composes one four-flavor chord, using its palindrome."""
        work = np.empty((5, 4, 4), dtype=np.complex128)
        acc_a = np.empty((4, 4), dtype=np.complex128)
        acc_b = np.empty((4, 4), dtype=np.complex128)
        tmp_a = np.empty((4, 4), dtype=np.complex128)
        tmp_b = np.empty((4, 4), dtype=np.complex128)
        h_work = np.empty((4, 4), dtype=np.complex128)

        for i in range(4):
            for j in range(4):
                acc_a[i, j] = 1.0 + 0.0j if i == j else 0.0 + 0.0j
                acc_b[i, j] = 1.0 + 0.0j if i == j else 0.0 + 0.0j

        n = widths.shape[0]
        m = n//2
        for k in range(m):
            _build_h_4nu(h_vac, inv_e, potentials, potentials_nc, k, h_work)
            _operator_4nu(h_work, widths[k], strategy, work)
            for i in range(4):
                for j in range(4):
                    sb = 0.0 + 0.0j
                    sa = 0.0 + 0.0j
                    for t in range(4):
                        sb += work[1][i, t]*acc_b[t, j]
                        sa += acc_a[i, t]*work[1][t, j]
                    tmp_b[i, j] = sb
                    tmp_a[i, j] = sa
            for i in range(4):
                for j in range(4):
                    acc_b[i, j] = tmp_b[i, j]
                    acc_a[i, j] = tmp_a[i, j]

        if n % 2 == 1:
            _build_h_4nu(h_vac, inv_e, potentials, potentials_nc, m, h_work)
            _operator_4nu(h_work, widths[m], strategy, work)
            for i in range(4):
                for j in range(4):
                    s = 0.0 + 0.0j
                    for t in range(4):
                        s += work[1][i, t]*acc_b[t, j]
                    tmp_b[i, j] = s
            for i in range(4):
                for j in range(4):
                    acc_b[i, j] = tmp_b[i, j]

        for i in range(4):
            for j in range(4):
                s = 0.0 + 0.0j
                for t in range(4):
                    s += acc_a[i, t]*acc_b[t, j]
                out[i, j] = s

    @njit(cache=True)
    def _earth_chords_2nu_mirrored_serial(h_vac, inv_energies, potentials,
                                          widths, out):
        for c in range(inv_energies.shape[0]):
            _earth_chord_2nu_mirrored(h_vac, inv_energies[c], potentials,
                                      widths, out[c])

    @njit(cache=True, parallel=True)
    def _earth_chords_2nu_mirrored_parallel(h_vac, inv_energies, potentials,
                                            widths, out):
        for c in prange(inv_energies.shape[0]):
            _earth_chord_2nu_mirrored(h_vac, inv_energies[c], potentials,
                                      widths, out[c])

    @njit(cache=True)
    def _earth_chords_3nu_mirrored_serial(h_vac, inv_energies, potentials,
                                          widths, out):
        for c in range(inv_energies.shape[0]):
            _earth_chord_3nu_mirrored(h_vac, inv_energies[c], potentials,
                                      widths, out[c])

    @njit(cache=True, parallel=True)
    def _earth_chords_3nu_mirrored_parallel(h_vac, inv_energies, potentials,
                                            widths, out):
        for c in prange(inv_energies.shape[0]):
            _earth_chord_3nu_mirrored(h_vac, inv_energies[c], potentials,
                                      widths, out[c])

    @njit(cache=True)
    def _earth_chords_4nu_mirrored_serial(h_vac, inv_energies, potentials,
                                          widths, potentials_nc, strategy,
                                          out):
        for c in range(inv_energies.shape[0]):
            _earth_chord_4nu_mirrored(h_vac, inv_energies[c], potentials,
                                      potentials_nc, widths, strategy, out[c])

    @njit(cache=True, parallel=True)
    def _earth_chords_4nu_mirrored_parallel(h_vac, inv_energies, potentials,
                                            widths, potentials_nc, strategy,
                                            out):
        for c in prange(inv_energies.shape[0]):
            _earth_chord_4nu_mirrored(h_vac, inv_energies[c], potentials,
                                      potentials_nc, widths, strategy, out[c])

    @njit(cache=True)
    def _earth_chords_2nu_serial(h_vac, inv_energies, potentials, widths, out):
        for c in range(inv_energies.shape[0]):
            _earth_chord_2nu(h_vac, inv_energies[c], potentials, widths,
                             out[c])

    @njit(cache=True, parallel=True)
    def _earth_chords_2nu_parallel(h_vac, inv_energies, potentials, widths,
                                   out):
        for c in prange(inv_energies.shape[0]):
            _earth_chord_2nu(h_vac, inv_energies[c], potentials, widths,
                             out[c])

    @njit(cache=True)
    def _earth_chords_3nu_serial(h_vac, inv_energies, potentials, widths, out):
        for c in range(inv_energies.shape[0]):
            _earth_chord_3nu(h_vac, inv_energies[c], potentials, widths,
                             out[c])

    @njit(cache=True, parallel=True)
    def _earth_chords_3nu_parallel(h_vac, inv_energies, potentials, widths,
                                   out):
        for c in prange(inv_energies.shape[0]):
            _earth_chord_3nu(h_vac, inv_energies[c], potentials, widths,
                             out[c])

    # The widths come before the four-flavor extras so that every chord
    # kernel takes the same first four arguments, which is what lets
    # `_run_earth_chords` dispatch all three through one call site.
    @njit(cache=True)
    def _earth_chords_4nu_serial(h_vac, inv_energies, potentials, widths,
                                 potentials_nc, strategy, out):
        for c in range(inv_energies.shape[0]):
            _earth_chord_4nu(h_vac, inv_energies[c], potentials,
                             potentials_nc, widths, strategy, out[c])

    @njit(cache=True, parallel=True)
    def _earth_chords_4nu_parallel(h_vac, inv_energies, potentials, widths,
                                   potentials_nc, strategy, out):
        for c in prange(inv_energies.shape[0]):
            _earth_chord_4nu(h_vac, inv_energies[c], potentials,
                             potentials_nc, widths, strategy, out[c])

    @njit(cache=True)
    def _run_2nu_serial(h_stack, l_stack, out):
        for n in range(h_stack.shape[0]):
            _one_2nu(h_stack[n], l_stack[n], out, n)

    @njit(cache=True, parallel=True)
    def _run_2nu_parallel(h_stack, l_stack, out):
        for n in prange(h_stack.shape[0]):
            _one_2nu(h_stack[n], l_stack[n], out, n)

    @njit(cache=True)
    def _run_4nu_serial(h_stack, l_stack, out, strategy):
        work = np.empty((5, 4, 4), dtype=np.complex128)
        for n in range(h_stack.shape[0]):
            _one_4nu(h_stack[n], l_stack[n], out, n, strategy, work)

    @njit(cache=True, parallel=True)
    def _run_4nu_parallel(h_stack, l_stack, out, strategy):
        for n in prange(h_stack.shape[0]):
            work = np.empty((5, 4, 4), dtype=np.complex128)
            _one_4nu(h_stack[n], l_stack[n], out, n, strategy, work)

    @njit(cache=True)
    def _run_4nu_u_serial(h_stack, l_stack, out, strategy):
        work = np.empty((5, 4, 4), dtype=np.complex128)
        for n in range(h_stack.shape[0]):
            _one_4nu_u(h_stack[n], l_stack[n], out, n, strategy, work)

    @njit(cache=True, parallel=True)
    def _run_4nu_u_parallel(h_stack, l_stack, out, strategy):
        for n in prange(h_stack.shape[0]):
            work = np.empty((5, 4, 4), dtype=np.complex128)
            _one_4nu_u(h_stack[n], l_stack[n], out, n, strategy, work)

    @njit(cache=True)
    def _slab_product_4nu(h_stack, widths, strategy, out):
        r"""Multiplies the per-slab four-flavor operators into ``out``.

        The four-flavor counterpart of `_slab_product_3nu`, and the same
        ordering: ``U = U_n ... U_1``, first slab crossed rightmost.
        """
        work = np.empty((5, 4, 4), dtype=np.complex128)
        acc = np.empty((4, 4), dtype=np.complex128)
        tmp = np.empty((4, 4), dtype=np.complex128)

        _operator_4nu(h_stack[0], widths[0], strategy, work)
        for i in range(4):
            for j in range(4):
                acc[i, j] = work[1][i, j]

        for k in range(1, widths.shape[0]):
            _operator_4nu(h_stack[k], widths[k], strategy, work)
            for i in range(4):
                for j in range(4):
                    total = 0.0 + 0.0j
                    for m in range(4):
                        total += work[1][i, m]*acc[m, j]
                    tmp[i, j] = total
            for i in range(4):
                for j in range(4):
                    acc[i, j] = tmp[i, j]

        for i in range(4):
            for j in range(4):
                out[i, j] = acc[i, j]

    def _run(
        h_stack: np.ndarray,
        l_stack: np.ndarray,
        width: int,
        serial: Callable,
        parallel: Callable,
        extra: tuple = (),
        dtype: type = float
    ) -> np.ndarray:
        r"""Flattens, dispatches, and restores the batch shape.

        Parameters
        ----------
        h_stack : numpy.ndarray
            Hamiltonians, of shape ``(..., width, width)``.
        l_stack : numpy.ndarray
            Baselines, of shape ``(...)``.
        width : int
            Number of flavors, 2, 3, or 4.
        serial : Callable
            Kernel to use below `PARALLEL_THRESHOLD`.
        parallel : Callable
            Kernel to use at or above it.
        extra : tuple, optional
            Further arguments passed on to the kernel after the output
            array.  Empty at two and three flavors; at four it carries
            the root-polishing switch, which an ``@njit`` function cannot
            read from module state at call time.
        dtype : type, optional
            Element type of the output.  ``float`` for probabilities,
            ``complex`` for evolution operators.

        Returns
        -------
        numpy.ndarray
            The probabilities, of shape ``(..., width*width)``.
        """
        batch = l_stack.shape
        flat_h = np.ascontiguousarray(h_stack).reshape(-1, width, width)
        flat_l = np.ascontiguousarray(l_stack).reshape(-1)
        out = np.empty((flat_l.shape[0], width*width), dtype=dtype)
        if flat_l.shape[0] >= PARALLEL_THRESHOLD:
            parallel(flat_h, flat_l, out, *extra)
        else:
            serial(flat_h, flat_l, out, *extra)
        return out.reshape(batch + (width*width,))

    def probabilities_3nu_kernel(
        h_stack: np.ndarray,
        l_stack: np.ndarray
    ) -> np.ndarray:
        r"""Returns the nine probabilities for a stack of Hamiltonians.

        .. versionadded:: 1.6.0

        Parameters
        ----------
        h_stack : numpy.ndarray
            Hamiltonians, of shape ``(..., 3, 3)``, already broadcast
            against `l_stack`.
        l_stack : numpy.ndarray
            Baselines, of shape ``(...)``.

        Returns
        -------
        numpy.ndarray
            The probabilities, of shape ``(..., 9)``, ordered with the
            initial flavor varying slowest --- the same ordering, and
            the same values to round-off, as the NumPy path.
        """
        return _run(h_stack, l_stack, 3, _run_3nu_serial, _run_3nu_parallel)

    def evolution_operator_3nu_kernel(
        h_stack: np.ndarray,
        l_stack: np.ndarray
    ) -> np.ndarray:
        r"""Returns :math:`U_3(L)` for a stack of Hamiltonians.

        .. versionadded:: 1.12.0

        The companion to `probabilities_3nu_kernel`, and the reason it
        exists: :mod:`slabs` and :mod:`earth` compose *operators* across
        adjacent slabs, so they cannot use a kernel that returns
        probabilities.  Without this they ran the NumPy path however the
        backend was configured --- having Numba installed bought an Earth
        crossing nothing at all.

        Parameters
        ----------
        h_stack : numpy.ndarray
            Hamiltonians, of shape ``(..., 3, 3)``, already broadcast
            against `l_stack`.
        l_stack : numpy.ndarray
            Baselines, of shape ``(...)``.

        Returns
        -------
        numpy.ndarray
            The evolution operators, of shape ``(..., 3, 3)``, indexed
            ``(final, initial)`` --- the same convention, and the same
            values to round-off, as the NumPy path.
        """
        flat = _run(h_stack, l_stack, 3,
                    _run_3nu_u_serial, _run_3nu_u_parallel,
                    dtype=complex)
        return flat.reshape(flat.shape[:-1] + (3, 3))

    def evolution_operator_2nu_kernel(
        h_stack: np.ndarray,
        l_stack: np.ndarray
    ) -> np.ndarray:
        r"""Returns :math:`U_2(L)` for a stack of Hamiltonians.

        .. versionadded:: 1.12.0

        Parameters
        ----------
        h_stack : numpy.ndarray
            Hamiltonians, of shape ``(..., 2, 2)``, already broadcast
            against `l_stack`.
        l_stack : numpy.ndarray
            Baselines, of shape ``(...)``.

        Returns
        -------
        numpy.ndarray
            The evolution operators, of shape ``(..., 2, 2)``, indexed
            ``(final, initial)``.
        """
        flat = _run(h_stack, l_stack, 2,
                    _run_2nu_u_serial, _run_2nu_u_parallel,
                    dtype=complex)
        return flat.reshape(flat.shape[:-1] + (2, 2))

    def slab_product_2nu_kernel(
        h_stack: np.ndarray,
        widths: np.ndarray
    ) -> np.ndarray:
        r"""Returns the composed :math:`U_2` across a sequence of slabs.

        .. versionadded:: 1.12.0

        Parameters
        ----------
        h_stack : numpy.ndarray
            One Hamiltonian per slab, of shape ``(n, 2, 2)``, ordered
            along the trajectory.
        widths : numpy.ndarray
            Slab widths, of shape ``(n,)``.

        Returns
        -------
        numpy.ndarray
            The product ``U_n ... U_1``, of shape ``(2, 2)``.
        """
        out = np.empty((2, 2), dtype=complex)
        h = np.ascontiguousarray(h_stack, dtype=complex)
        w = np.ascontiguousarray(widths, dtype=float)
        if worthwhile_mirror(2, w.shape[0]) and _palindromic_stack(h, w):
            _slab_product_2nu_mirrored(h, w, out)
        else:
            _slab_product_2nu(h, w, out)
        return out

    def evolution_operator_4nu_kernel(
        h_stack: np.ndarray,
        l_stack: np.ndarray,
        strategy: int = 0
    ) -> np.ndarray:
        r"""Returns :math:`U_4(L)` for a stack of Hamiltonians.

        .. versionadded:: 1.12.0

        Parameters
        ----------
        h_stack : numpy.ndarray
            Hamiltonians, of shape ``(..., 4, 4)``, already broadcast
            against `l_stack`.
        l_stack : numpy.ndarray
            Baselines, of shape ``(...)``.
        strategy : int, optional
            How to obtain the latent roots, as
            :data:`oscprob4nu.ROOT_STRATEGY` and
            :data:`oscprob4nu.POLISH_ROOTS` select between on the NumPy
            path: ``0`` for double-double invariants with Aberth
            refinement, ``1`` for the eigensolver with a Newton step,
            ``2`` for the eigensolver alone.

        Returns
        -------
        numpy.ndarray
            The evolution operators, of shape ``(..., 4, 4)``, indexed
            ``(final, initial)``.
        """
        flat = _run(h_stack, l_stack, 4,
                    _run_4nu_u_serial, _run_4nu_u_parallel,
                    (int(strategy),), dtype=complex)
        return flat.reshape(flat.shape[:-1] + (4, 4))

    def slab_product_4nu_kernel(
        h_stack: np.ndarray,
        widths: np.ndarray,
        strategy: int = 0
    ) -> np.ndarray:
        r"""Returns the composed :math:`U_4` across a sequence of slabs.

        .. versionadded:: 1.12.0

        Parameters
        ----------
        h_stack : numpy.ndarray
            One Hamiltonian per slab, of shape ``(n, 4, 4)``, ordered
            along the trajectory.
        widths : numpy.ndarray
            Slab widths, of shape ``(n,)``.
        strategy : int, optional
            How to obtain the latent roots: ``0`` for double-double
            invariants with Aberth refinement, ``1`` for the eigensolver
            with a Newton step, ``2`` for the eigensolver alone.

        Returns
        -------
        numpy.ndarray
            The product ``U_n ... U_1``, of shape ``(4, 4)``.
        """
        out = np.empty((4, 4), dtype=complex)
        h = np.ascontiguousarray(h_stack, dtype=complex)
        w = np.ascontiguousarray(widths, dtype=float)
        if worthwhile_mirror(4, w.shape[0]) and _palindromic_stack(h, w):
            _slab_product_4nu_mirrored(h, w, int(strategy), out)
        else:
            _slab_product_4nu(h, w, int(strategy), out)
        return out

    def slab_product_3nu_kernel(
        h_stack: np.ndarray,
        widths: np.ndarray
    ) -> np.ndarray:
        r"""Returns the composed :math:`U_3` across a sequence of slabs.

        .. versionadded:: 1.12.0

        What :mod:`slabs` and :mod:`earth` actually want.  Computing the
        per-slab operators in a kernel and then multiplying them in a
        Python loop left the loop as the dominant cost of an Earth
        crossing; this does both in one pass.

        Parameters
        ----------
        h_stack : numpy.ndarray
            One Hamiltonian per slab, of shape ``(n, 3, 3)``, ordered
            along the trajectory.
        widths : numpy.ndarray
            Slab widths, of shape ``(n,)``, in units reciprocal to the
            Hamiltonian.

        Returns
        -------
        numpy.ndarray
            The product ``U_n ... U_1``, of shape ``(3, 3)``.
        """
        out = np.empty((3, 3), dtype=complex)
        h = np.ascontiguousarray(h_stack, dtype=complex)
        w = np.ascontiguousarray(widths, dtype=float)
        # A sequence that reads the same from either end has every
        # operator twice over; see `palindromic`.
        if worthwhile_mirror(3, w.shape[0]) and _palindromic_stack(h, w):
            _slab_product_3nu_mirrored(h, w, out)
        else:
            _slab_product_3nu(h, w, out)
        return out

    def _run_slab_batch(
        h_stack: np.ndarray,
        widths: np.ndarray,
        width: int,
        serial: Callable,
        parallel: Callable,
        extra: tuple = ()
    ) -> np.ndarray:
        r"""Flattens, dispatches, and restores the chord batch shape.

        The parallel decision is made on ``n_chords*n_slabs`` rather than
        on the number of chords, because the work a chord costs is
        proportional to how many slabs it has: sixteen chords of
        nineteen slabs and four of a hundred and twenty are the same
        amount of arithmetic, and measurement puts the crossover at the
        same place for both.  Using the chord count alone would have put
        it in three different places for three geometries.

        Parameters
        ----------
        h_stack : numpy.ndarray
            Hamiltonians, of shape ``(..., n_slabs, width, width)``.
        widths : numpy.ndarray
            Slab widths, of shape ``(n_slabs,)``, shared by every chord.
        width : int
            Number of flavors, 2, 3, or 4.
        serial : Callable
            Kernel to use below `PARALLEL_THRESHOLD`.
        parallel : Callable
            Kernel to use at or above it.
        extra : tuple, optional
            Further arguments passed on before the output array; at four
            flavors it carries the root-polishing switch.

        Returns
        -------
        numpy.ndarray
            The composed operators, of shape ``(..., width, width)``.
        """
        batch = h_stack.shape[:-3]
        n_slabs = h_stack.shape[-3]
        flat_h = np.ascontiguousarray(h_stack, dtype=complex).reshape(
            -1, n_slabs, width, width)
        flat_w = np.ascontiguousarray(widths, dtype=float)
        out = np.empty((flat_h.shape[0], width, width), dtype=complex)
        if flat_h.shape[0]*n_slabs >= PARALLEL_THRESHOLD:
            parallel(flat_h, flat_w, *extra, out)
        else:
            serial(flat_h, flat_w, *extra, out)
        return out.reshape(batch + (width, width))

    def _run_earth_chords(
        h_vacuum: np.ndarray,
        energies: np.ndarray,
        potentials: np.ndarray,
        widths: np.ndarray,
        width: int,
        serial: Callable,
        parallel: Callable,
        extra: tuple = (),
        mirrored: tuple = ()
    ) -> np.ndarray:
        r"""Dispatches a fused chord kernel and restores the batch shape.

        Parameters
        ----------
        h_vacuum : numpy.ndarray
            Energy-independent vacuum Hamiltonian, of shape
            ``(width, width)``.
        energies : numpy.ndarray
            Neutrino energies, of any shape.
        potentials : numpy.ndarray
            Matter potentials, of shape ``(n_slabs,)``, one per slab.
        widths : numpy.ndarray
            Slab widths, of shape ``(n_slabs,)``.
        width : int
            Number of flavors, 2, 3, or 4.
        serial : Callable
            Kernel to use below `PARALLEL_THRESHOLD`.
        parallel : Callable
            Kernel to use at or above it.
        extra : tuple, optional
            Further arguments passed on before the output array.
        mirrored : tuple, optional
            The serial and parallel composers to use when the chord is a
            palindrome and long enough to be worth it.  Empty to always
            compose the whole chord.

        Returns
        -------
        numpy.ndarray
            The evolution operators, of shape
            ``(..., width, width)``.
        """
        batch = np.shape(energies)
        flat_e = np.ascontiguousarray(energies, dtype=float).reshape(-1)
        # The kernel wants a reciprocal per chord, and computing it here
        # keeps a division out of the innermost loop.
        inv_e = 1.0/flat_e
        h_vac = np.ascontiguousarray(h_vacuum, dtype=complex)
        pot = np.ascontiguousarray(potentials, dtype=float)
        wid = np.ascontiguousarray(widths, dtype=float)
        out = np.empty((flat_e.shape[0], width, width), dtype=complex)

        # A palindromic chord costs about two thirds as much, and every
        # Earth chord is one; the check is a pass over two small arrays
        # against a per-slab expansion, so it is free at this scale.
        if (mirrored and worthwhile_mirror(width, wid.shape[0])
                and palindromic(pot, wid, *extra[:1])):
            serial, parallel = mirrored

        if flat_e.shape[0]*wid.shape[0] >= PARALLEL_THRESHOLD:
            parallel(h_vac, inv_e, pot, wid, *extra, out)
        else:
            serial(h_vac, inv_e, pot, wid, *extra, out)

        return out.reshape(batch + (width, width))

    def earth_chords_2nu_kernel(
        h_vacuum: np.ndarray,
        energies: np.ndarray,
        potentials: np.ndarray,
        widths: np.ndarray
    ) -> np.ndarray:
        r"""Returns one composed :math:`U_2` per energy, Hamiltonians
        built inline.

        .. versionadded:: 1.12.0

        Parameters
        ----------
        h_vacuum : numpy.ndarray
            Energy-independent vacuum Hamiltonian, of shape ``(2, 2)``.
        energies : numpy.ndarray
            Neutrino energies, in units of eV.
        potentials : numpy.ndarray
            Charged-current potentials, of shape ``(n_slabs,)``.
        widths : numpy.ndarray
            Slab widths, of shape ``(n_slabs,)``.

        Returns
        -------
        numpy.ndarray
            The products, of shape ``(..., 2, 2)``.
        """
        return _run_earth_chords(h_vacuum, energies, potentials, widths, 2,
                                 _earth_chords_2nu_serial,
                                 _earth_chords_2nu_parallel,
                                 mirrored=(_earth_chords_2nu_mirrored_serial,
                                           _earth_chords_2nu_mirrored_parallel))

    def earth_chords_3nu_kernel(
        h_vacuum: np.ndarray,
        energies: np.ndarray,
        potentials: np.ndarray,
        widths: np.ndarray
    ) -> np.ndarray:
        r"""Returns one composed :math:`U_3` per energy, Hamiltonians
        built inline.

        .. versionadded:: 1.12.0

        The fused counterpart of `slab_product_3nu_batch_kernel`, and
        what an Earth energy scan actually calls.  The batch kernel
        takes a stack of Hamiltonians, which for a scan means
        materialising one 3x3 matrix per slab per energy and streaming
        it back --- 17 KB per chord, against the two kilobyte-scale
        arrays this reads, shared by every chord.  That is the
        difference between a memory-bound kernel and a compute-bound
        one.

        Parameters
        ----------
        h_vacuum : numpy.ndarray
            Energy-independent vacuum Hamiltonian, of shape ``(3, 3)``,
            in units of eV\ :sup:`2`.
        energies : numpy.ndarray
            Neutrino energies, in units of eV, of any shape.
        potentials : numpy.ndarray
            Charged-current potentials at the slab midpoints, of shape
            ``(n_slabs,)``.  Shared by every energy, the potential
            depending on the geometry alone.
        widths : numpy.ndarray
            Slab widths, of shape ``(n_slabs,)``, in units of
            eV\ :sup:`-1`.

        Returns
        -------
        numpy.ndarray
            The products ``U_n ... U_1``, of shape ``(..., 3, 3)``.
        """
        return _run_earth_chords(h_vacuum, energies, potentials, widths, 3,
                                 _earth_chords_3nu_serial,
                                 _earth_chords_3nu_parallel,
                                 mirrored=(_earth_chords_3nu_mirrored_serial,
                                           _earth_chords_3nu_mirrored_parallel))

    def earth_chords_4nu_kernel(
        h_vacuum: np.ndarray,
        energies: np.ndarray,
        potentials: np.ndarray,
        potentials_nc: np.ndarray,
        widths: np.ndarray,
        strategy: int = 0
    ) -> np.ndarray:
        r"""Returns one composed :math:`U_4` per energy, Hamiltonians
        built inline.

        .. versionadded:: 1.12.0

        Parameters
        ----------
        h_vacuum : numpy.ndarray
            Energy-independent vacuum Hamiltonian, of shape ``(4, 4)``.
        energies : numpy.ndarray
            Neutrino energies, in units of eV.
        potentials : numpy.ndarray
            Charged-current potentials, of shape ``(n_slabs,)``.
        potentials_nc : numpy.ndarray
            Neutral-current potentials, of the same shape.
        widths : numpy.ndarray
            Slab widths, of shape ``(n_slabs,)``.
        strategy : int, optional
            How to obtain the latent roots: ``0`` for double-double
            invariants with Aberth refinement, ``1`` for the eigensolver
            with a Newton step, ``2`` for the eigensolver alone.

        Returns
        -------
        numpy.ndarray
            The products, of shape ``(..., 4, 4)``.
        """
        return _run_earth_chords(
            h_vacuum, energies, potentials, widths, 4,
            _earth_chords_4nu_serial, _earth_chords_4nu_parallel,
            (np.ascontiguousarray(potentials_nc, dtype=float),
             int(strategy)),
            mirrored=(_earth_chords_4nu_mirrored_serial,
                      _earth_chords_4nu_mirrored_parallel))

    def slab_product_2nu_batch_kernel(
        h_stack: np.ndarray,
        widths: np.ndarray
    ) -> np.ndarray:
        r"""Returns one composed :math:`U_2` per chord.

        .. versionadded:: 1.12.0

        Parameters
        ----------
        h_stack : numpy.ndarray
            Hamiltonians, of shape ``(..., n_slabs, 2, 2)``.
        widths : numpy.ndarray
            Slab widths, of shape ``(n_slabs,)``, shared by every chord.

        Returns
        -------
        numpy.ndarray
            The products, of shape ``(..., 2, 2)``.
        """
        return _run_slab_batch(h_stack, widths, 2,
                               _slab_product_2nu_batch_serial,
                               _slab_product_2nu_batch_parallel)

    def slab_product_3nu_batch_kernel(
        h_stack: np.ndarray,
        widths: np.ndarray
    ) -> np.ndarray:
        r"""Returns one composed :math:`U_3` per chord.

        .. versionadded:: 1.12.0

        What an energy scan across the Earth wants.  `earth` calls this
        once for a whole array of energies, where before it made one
        `slab_product_3nu_kernel` call per energy from Python and
        rebuilt the energy-independent matter potentials every time.

        Parameters
        ----------
        h_stack : numpy.ndarray
            Hamiltonians, of shape ``(..., n_slabs, 3, 3)``, one chord
            per leading index and one Hamiltonian per slab along it.
        widths : numpy.ndarray
            Slab widths, of shape ``(n_slabs,)``, shared by every chord.
            A scan over energy at fixed zenith angle crosses the same
            geometry every time, which is what makes one array of widths
            enough.

        Returns
        -------
        numpy.ndarray
            The products ``U_n ... U_1``, of shape ``(..., 3, 3)``.
        """
        return _run_slab_batch(h_stack, widths, 3,
                               _slab_product_3nu_batch_serial,
                               _slab_product_3nu_batch_parallel)

    def slab_product_4nu_batch_kernel(
        h_stack: np.ndarray,
        widths: np.ndarray,
        strategy: int = 0
    ) -> np.ndarray:
        r"""Returns one composed :math:`U_4` per chord.

        .. versionadded:: 1.12.0

        Parameters
        ----------
        h_stack : numpy.ndarray
            Traceless Hamiltonians, of shape ``(..., n_slabs, 4, 4)``,
            as `slab_product_4nu_kernel` also expects.
        widths : numpy.ndarray
            Slab widths, of shape ``(n_slabs,)``, shared by every chord.
        strategy : int, optional
            How to obtain the latent roots: ``0`` for double-double
            invariants with Aberth refinement, ``1`` for the eigensolver
            with a Newton step, ``2`` for the eigensolver alone.

        Returns
        -------
        numpy.ndarray
            The products, of shape ``(..., 4, 4)``.
        """
        return _run_slab_batch(h_stack, widths, 4,
                               _slab_product_4nu_batch_serial,
                               _slab_product_4nu_batch_parallel,
                               (int(strategy),))

    def probabilities_2nu_kernel(
        h_stack: np.ndarray,
        l_stack: np.ndarray
    ) -> np.ndarray:
        r"""Returns the four probabilities for a stack of Hamiltonians.

        .. versionadded:: 1.6.0

        Parameters
        ----------
        h_stack : numpy.ndarray
            Hamiltonians, of shape ``(..., 2, 2)``, already broadcast
            against `l_stack`.
        l_stack : numpy.ndarray
            Baselines, of shape ``(...)``.

        Returns
        -------
        numpy.ndarray
            The probabilities, of shape ``(..., 4)``, ordered
            ``(Pee, Pem, Pme, Pmm)``.
        """
        return _run(h_stack, l_stack, 2, _run_2nu_serial, _run_2nu_parallel)

    def probabilities_4nu_kernel(
        h_stack: np.ndarray,
        l_stack: np.ndarray,
        strategy: int = 0
    ) -> np.ndarray:
        r"""Returns the sixteen probabilities for a stack of Hamiltonians.

        .. versionadded:: 1.10.0

        Parameters
        ----------
        h_stack : numpy.ndarray
            Hamiltonians, of shape ``(..., 4, 4)``, already broadcast
            against `l_stack`.
        l_stack : numpy.ndarray
            Baselines, of shape ``(...)``.
        strategy : int, optional
            How to obtain the latent roots, as
            :data:`oscprob4nu.ROOT_STRATEGY` and
            :data:`oscprob4nu.POLISH_ROOTS` ask the NumPy path to:
            ``0`` for double-double invariants with Aberth refinement,
            ``1`` for the eigensolver with a Newton step, ``2`` for the
            eigensolver alone.  It is an argument rather than a module
            constant because a compiled kernel cannot read a Python
            global at call time without recompiling.

        Returns
        -------
        numpy.ndarray
            The probabilities, of shape ``(..., 16)``, ordered with the
            initial flavor varying slowest --- the same ordering, and
            the same values to round-off, as the NumPy path.
        """
        return _run(h_stack, l_stack, 4, _run_4nu_serial, _run_4nu_parallel,
                    (int(strategy),))
