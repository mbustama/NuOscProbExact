# -*- coding: utf-8 -*-
r"""Optional Numba-compiled kernels for the batched evaluation paths.

**NuOscProbExact** needs only NumPy.  If `Numba <https://numba.pydata.org>`_
happens to be installed, this module compiles the two- and three-neutrino
expansions into fused machine-code loops and :mod:`oscprob2nu` and
:mod:`oscprob3nu` use them for large stacks; if it is not, ``HAVE_NUMBA``
is ``False``, nothing here is defined, and the NumPy path is used
instead.  Nothing else in the library changes either way, and the
results agree to round-off --- see ``tests/test_fastkernels.py``, which
runs both paths against each other whichever is available.

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
200 000 energies, three flavors  ~15x
20 000 energies, three flavors   ~9x
100 x 100 oscillogram            ~3.5x
200 000 baselines, two flavors   ~1.5x
===============================  ==========

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
about 16 microseconds, which is not worth a compilation pause.

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
"""

__version__ = "1.6"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

__all__ = ['HAVE_NUMBA', 'USE_NUMBA', 'MIN_BATCH', 'PARALLEL_THRESHOLD',
           'available', 'worthwhile',
           'probabilities_2nu_kernel', 'probabilities_3nu_kernel']

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

MIN_BATCH = {2: 50000, 3: 1}
r"""dict: Module-level constant.

The smallest stack for which the compiled kernel is worth using, by
number of flavors.  A backend that is sometimes slower than the path it
replaces is worse than no backend, so these are measured rather than
assumed.

For three flavors the kernel wins at every size, by between two and
sixteen times, so the threshold is one.  For two flavors it does not:
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
        Whether `probabilities_2nu_kernel` and
        `probabilities_3nu_kernel` may be called.
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
        Number of neutrino flavors, 2 or 3.
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
    def _one_3nu(h_matrix, L, out, n):
        r"""Writes the nine probabilities for one Hamiltonian into
        ``out[n]``.

        A transcription of the scalar path in :mod:`oscprob3nu`: the
        SU(3) coefficients, the sparse star product, the two invariants,
        the latent roots with the same degeneracy handling, and the nine
        moduli squared.
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

    @njit(cache=True)
    def _run_3nu_serial(h_stack, l_stack, out):
        for n in range(h_stack.shape[0]):
            _one_3nu(h_stack[n], l_stack[n], out, n)

    @njit(cache=True, parallel=True)
    def _run_3nu_parallel(h_stack, l_stack, out):
        for n in prange(h_stack.shape[0]):
            _one_3nu(h_stack[n], l_stack[n], out, n)

    @njit(cache=True)
    def _run_2nu_serial(h_stack, l_stack, out):
        for n in range(h_stack.shape[0]):
            _one_2nu(h_stack[n], l_stack[n], out, n)

    @njit(cache=True, parallel=True)
    def _run_2nu_parallel(h_stack, l_stack, out):
        for n in prange(h_stack.shape[0]):
            _one_2nu(h_stack[n], l_stack[n], out, n)

    def _run(
        h_stack: np.ndarray,
        l_stack: np.ndarray,
        width: int,
        serial: Callable,
        parallel: Callable
    ) -> np.ndarray:
        r"""Flattens, dispatches, and restores the batch shape.

        Parameters
        ----------
        h_stack : numpy.ndarray
            Hamiltonians, of shape ``(..., width, width)``.
        l_stack : numpy.ndarray
            Baselines, of shape ``(...)``.
        width : int
            Number of flavors, 2 or 3.
        serial : Callable
            Kernel to use below `PARALLEL_THRESHOLD`.
        parallel : Callable
            Kernel to use at or above it.

        Returns
        -------
        numpy.ndarray
            The probabilities, of shape ``(..., width*width)``.
        """
        batch = l_stack.shape
        flat_h = np.ascontiguousarray(h_stack).reshape(-1, width, width)
        flat_l = np.ascontiguousarray(l_stack).reshape(-1)
        out = np.empty((flat_l.shape[0], width*width))
        if flat_l.shape[0] >= PARALLEL_THRESHOLD:
            parallel(flat_h, flat_l, out)
        else:
            serial(flat_h, flat_l, out)
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
