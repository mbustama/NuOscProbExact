# -*- coding: utf-8 -*-
r"""Optional Numba-compiled kernels for the batched evaluation paths.

**NuOscProbExact** needs only NumPy.  If `Numba <https://numba.pydata.org>`_
happens to be installed, this module compiles the two-, three- and
four-neutrino expansions into fused machine-code loops and
:mod:`oscprob2nu`, :mod:`oscprob3nu` and :mod:`oscprob4nu` use them for
large stacks; if it is not, ``HAVE_NUMBA`` is ``False``, nothing here is
defined, and the NumPy path is used instead.  Nothing else in the
library changes either way, and the results agree to round-off --- see
``tests/test_fastkernels.py``, which runs both paths against each other
whichever is available.

Install the optional dependency with::

    pip install "nuoscprobexact[fast]"

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
Earth crossing at 120 slabs comes out 13.9x, 9.6x and 6.6x quicker at
two, three and four flavors.

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

Both are why this is an optional extra rather than a dependency, and why
the scalar path is deliberately left alone: a single probability takes
about 8 microseconds, which is not worth a compilation pause.

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

__all__ = ['HAVE_NUMBA', 'USE_NUMBA', 'MIN_BATCH', 'MIN_SLAB_BATCH',
           'PARALLEL_THRESHOLD',
           'available', 'worthwhile', 'worthwhile_slabs',
           'probabilities_2nu_kernel', 'probabilities_3nu_kernel',
           'probabilities_4nu_kernel', 'evolution_operator_3nu_kernel',
           'slab_product_3nu_kernel',
           'evolution_operator_4nu_kernel',
           'slab_product_4nu_kernel',
           'evolution_operator_2nu_kernel',
           'slab_product_2nu_kernel']

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


def worthwhile(n_flavors: int, size: int) -> bool:
    r"""Returns whether the compiled kernel should be used for a stack.

    The kernels are only used where they have been measured to win.
    Below the per-flavor threshold in `MIN_BATCH` the NumPy path is
    quicker, and using the kernel anyway would make installing the
    optional extra a pessimisation for those calls.

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

    @njit(cache=True, inline='always')
    def _operator_4nu(h_matrix, L, polish, work):
        r"""Builds :math:`U_4(L)` for one Hamiltonian in ``work[1]``.

        Factored out of `_one_4nu`, which computed the operator and
        then squared it, so that the probability kernel and the
        evolution-operator kernel share the expensive part.  Unlike
        three flavors the operator is already materialised here, in
        the caller's scratch, so the split costs nothing at all.

        A transcription of :func:`oscprob4nu._evolution_operator_4nu_array`
        for a single element: the traceless part, the three invariants
        from traces of powers, the quartic by Euler's reduction, the
        Newton refinement of the roots against the matrix, the divided
        differences of the exponential over them, and the Newton-form
        reconstruction of :math:`U_4`.

        ``work`` is scratch space of shape ``(5, 4, 4)``, supplied by the
        caller so that the loop over a stack allocates nothing.
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

        # The invariants, from traces of powers of the traceless part.
        # `second` holds H~^2 until the reconstruction needs it back.
        for i in range(4):
            for j in range(4):
                entry = 0.0j
                for k in range(4):
                    entry += traceless[i, k]*traceless[k, j]
                second[i, j] = entry

        trace_2 = 0.0j
        trace_3 = 0.0j
        trace_4 = 0.0j
        for i in range(4):
            trace_2 += second[i, i]
            row_3 = 0.0j
            row_4 = 0.0j
            for j in range(4):
                row_3 += second[i, j]*traceless[j, i]
                row_4 += second[i, j]*second[j, i]
            trace_3 += row_3
            trace_4 += row_4

        invariant_2 = 0.5*trace_2.real
        invariant_3 = 0.5*trace_3.real
        invariant_4 = 0.5*(trace_4.real - invariant_2*invariant_2)

        # Euler's reduction: the resolvent cubic, solved trigonometrically
        quadratic = -invariant_2
        linear = -(2.0/3.0)*invariant_3
        constant = 0.25*(invariant_2*invariant_2 - 2.0*invariant_4)

        coeff_2 = 2.0*quadratic
        coeff_1 = quadratic*quadratic - 4.0*constant
        coeff_0 = -linear*linear

        depressed_p = coeff_1 - coeff_2*coeff_2/3.0
        depressed_q = (2.0*coeff_2*coeff_2*coeff_2/27.0
                       - coeff_2*coeff_1/3.0 + coeff_0)
        shift = -coeff_2/3.0

        scale = 2.0*math.sqrt(max(-depressed_p, 0.0)/3.0)
        denominator = depressed_p*scale
        if denominator != 0.0:
            argument = 3.0*depressed_q/denominator
        else:
            argument = 3.0*depressed_q
        if argument < -1.0:
            argument = -1.0
        elif argument > 1.0:
            argument = 1.0
        angle = math.acos(argument)

        root_0 = math.sqrt(max(scale*math.cos(angle/3.0) + shift, 0.0))
        root_1 = math.sqrt(max(scale*math.cos((angle + 2.0*math.pi)/3.0)
                               + shift, 0.0))
        root_2 = math.sqrt(max(scale*math.cos((angle + 4.0*math.pi)/3.0)
                               + shift, 0.0))
        if linear > 0.0:
            root_2 = -root_2

        psi_0 = 0.5*(root_0 + root_1 + root_2)
        psi_1 = 0.5*(root_0 - root_1 - root_2)
        psi_2 = 0.5*(-root_0 + root_1 - root_2)
        psi_3 = 0.5*(-root_0 - root_1 + root_2)

        # Ascending, by the five-comparator network for four elements
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

        if polish:
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
    def _one_4nu(h_matrix, L, out, n, polish, work):
        r"""Writes the sixteen probabilities into ``out[n]``."""
        _operator_4nu(h_matrix, L, polish, work)
        operator = work[1]

        # P_ab = |U_ba|^2, initial flavor slowest
        for alpha in range(4):
            for beta in range(4):
                entry = operator[beta, alpha]
                out[n, 4*alpha + beta] = (entry.real*entry.real
                                          + entry.imag*entry.imag)

    @njit(cache=True, inline='always')
    def _one_4nu_u(h_matrix, L, out, n, polish, work):
        r"""Writes the sixteen entries of :math:`U_4(L)` into ``out[n]``.

        Row-major and indexed ``(final, initial)``, so reshaping to
        ``(4, 4)`` gives the matrix `oscprob4nu` returns --- not the
        flavor order the probabilities use, which runs the initial index
        slowest.  The two differ by a transpose.
        """
        _operator_4nu(h_matrix, L, polish, work)
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
    def _slab_product_4nu_batch_serial(h_stack, widths, polish, out):
        for c in range(h_stack.shape[0]):
            _slab_product_4nu(h_stack[c], widths, polish, out[c])

    @njit(cache=True, parallel=True)
    def _slab_product_4nu_batch_parallel(h_stack, widths, polish, out):
        for c in prange(h_stack.shape[0]):
            _slab_product_4nu(h_stack[c], widths, polish, out[c])

    @njit(cache=True)
    def _run_2nu_serial(h_stack, l_stack, out):
        for n in range(h_stack.shape[0]):
            _one_2nu(h_stack[n], l_stack[n], out, n)

    @njit(cache=True, parallel=True)
    def _run_2nu_parallel(h_stack, l_stack, out):
        for n in prange(h_stack.shape[0]):
            _one_2nu(h_stack[n], l_stack[n], out, n)

    @njit(cache=True)
    def _run_4nu_serial(h_stack, l_stack, out, polish):
        work = np.empty((5, 4, 4), dtype=np.complex128)
        for n in range(h_stack.shape[0]):
            _one_4nu(h_stack[n], l_stack[n], out, n, polish, work)

    @njit(cache=True, parallel=True)
    def _run_4nu_parallel(h_stack, l_stack, out, polish):
        for n in prange(h_stack.shape[0]):
            work = np.empty((5, 4, 4), dtype=np.complex128)
            _one_4nu(h_stack[n], l_stack[n], out, n, polish, work)

    @njit(cache=True)
    def _run_4nu_u_serial(h_stack, l_stack, out, polish):
        work = np.empty((5, 4, 4), dtype=np.complex128)
        for n in range(h_stack.shape[0]):
            _one_4nu_u(h_stack[n], l_stack[n], out, n, polish, work)

    @njit(cache=True, parallel=True)
    def _run_4nu_u_parallel(h_stack, l_stack, out, polish):
        for n in prange(h_stack.shape[0]):
            work = np.empty((5, 4, 4), dtype=np.complex128)
            _one_4nu_u(h_stack[n], l_stack[n], out, n, polish, work)

    @njit(cache=True)
    def _slab_product_4nu(h_stack, widths, polish, out):
        r"""Multiplies the per-slab four-flavor operators into ``out``.

        The four-flavor counterpart of `_slab_product_3nu`, and the same
        ordering: ``U = U_n ... U_1``, first slab crossed rightmost.
        """
        work = np.empty((5, 4, 4), dtype=np.complex128)
        acc = np.empty((4, 4), dtype=np.complex128)
        tmp = np.empty((4, 4), dtype=np.complex128)

        _operator_4nu(h_stack[0], widths[0], polish, work)
        for i in range(4):
            for j in range(4):
                acc[i, j] = work[1][i, j]

        for k in range(1, widths.shape[0]):
            _operator_4nu(h_stack[k], widths[k], polish, work)
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
        backend was configured --- installing the optional extra bought
        an Earth crossing nothing at all.

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
        _slab_product_2nu(np.ascontiguousarray(h_stack, dtype=complex),
                          np.ascontiguousarray(widths, dtype=float), out)
        return out

    def evolution_operator_4nu_kernel(
        h_stack: np.ndarray,
        l_stack: np.ndarray,
        polish: bool = True
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
        polish : bool, optional
            Whether to refine the latent roots against the Hamiltonian,
            as :data:`oscprob4nu.POLISH_ROOTS` does on the NumPy path.

        Returns
        -------
        numpy.ndarray
            The evolution operators, of shape ``(..., 4, 4)``, indexed
            ``(final, initial)``.
        """
        flat = _run(h_stack, l_stack, 4,
                    _run_4nu_u_serial, _run_4nu_u_parallel,
                    (bool(polish),), dtype=complex)
        return flat.reshape(flat.shape[:-1] + (4, 4))

    def slab_product_4nu_kernel(
        h_stack: np.ndarray,
        widths: np.ndarray,
        polish: bool = True
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
        polish : bool, optional
            Whether to refine the latent roots against the Hamiltonian.

        Returns
        -------
        numpy.ndarray
            The product ``U_n ... U_1``, of shape ``(4, 4)``.
        """
        out = np.empty((4, 4), dtype=complex)
        _slab_product_4nu(np.ascontiguousarray(h_stack, dtype=complex),
                          np.ascontiguousarray(widths, dtype=float),
                          bool(polish), out)
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
        _slab_product_3nu(np.ascontiguousarray(h_stack, dtype=complex),
                          np.ascontiguousarray(widths, dtype=float), out)
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
        polish: bool = True
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
        polish : bool, optional
            Whether to refine the latent roots against the Hamiltonian.

        Returns
        -------
        numpy.ndarray
            The products, of shape ``(..., 4, 4)``.
        """
        return _run_slab_batch(h_stack, widths, 4,
                               _slab_product_4nu_batch_serial,
                               _slab_product_4nu_batch_parallel,
                               (bool(polish),))

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
        polish: bool = True
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
        polish : bool, optional
            Whether to refine the latent roots against the Hamiltonian
            matrix, as :data:`oscprob4nu.POLISH_ROOTS` asks the NumPy
            path to.  It is an argument rather than a module constant
            because a compiled kernel cannot read a Python global at
            call time without recompiling.

        Returns
        -------
        numpy.ndarray
            The probabilities, of shape ``(..., 16)``, ordered with the
            initial flavor varying slowest --- the same ordering, and
            the same values to round-off, as the NumPy path.
        """
        return _run(h_stack, l_stack, 4, _run_4nu_serial, _run_4nu_parallel,
                    (bool(polish),))
