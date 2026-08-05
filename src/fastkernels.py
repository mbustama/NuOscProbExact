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
    * probabilities_2nu_kernel - Two-flavor probabilities for a stack
    * probabilities_3nu_kernel - Three-flavor probabilities for a stack
    * probabilities_4nu_kernel - Four-flavor probabilities for a stack
    * evolution_operator_3nu_kernel - Three-flavor U_3 for a stack
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

__all__ = ['HAVE_NUMBA', 'USE_NUMBA', 'MIN_BATCH', 'PARALLEL_THRESHOLD',
           'available', 'worthwhile',
           'probabilities_2nu_kernel', 'probabilities_3nu_kernel',
           'probabilities_4nu_kernel', 'evolution_operator_3nu_kernel']

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
    def _one_4nu(h_matrix, L, out, n, polish, work):
        r"""Writes the sixteen probabilities for one Hamiltonian into
        ``out[n]``.

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

        # P_ab = |U_ba|^2, initial flavor slowest
        for alpha in range(4):
            for beta in range(4):
                entry = operator[beta, alpha]
                out[n, 4*alpha + beta] = (entry.real*entry.real
                                          + entry.imag*entry.imag)

    @njit(cache=True)
    def _run_3nu_serial(h_stack, l_stack, out):
        for n in range(h_stack.shape[0]):
            _one_3nu(h_stack[n], l_stack[n], out, n)

    @njit(cache=True, parallel=True)
    def _run_3nu_parallel(h_stack, l_stack, out):
        for n in prange(h_stack.shape[0]):
            _one_3nu(h_stack[n], l_stack[n], out, n)

    @njit(cache=True)
    def _run_3nu_u_serial(h_stack, l_stack, out):
        for n in range(l_stack.shape[0]):
            _one_3nu_u(h_stack[n], l_stack[n], out, n)

    @njit(cache=True, parallel=True)
    def _run_3nu_u_parallel(h_stack, l_stack, out):
        for n in prange(l_stack.shape[0]):
            _one_3nu_u(h_stack[n], l_stack[n], out, n)

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
