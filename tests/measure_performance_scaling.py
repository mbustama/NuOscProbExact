# -*- coding: utf-8 -*-
r"""Measures the three routes through this library, once, and freezes them.

The scaling curves of the performance figure used to be measured inside the
notebook, so every rebuild re-timed them and the figure moved.  Across a
single evening's rebuilds the array gain came out 29x, 33x and 35x, and the
kernel gain between 8.0x and 9.6x, while the paper quoted fixed bands in two
places.  A figure whose numbers change when nothing changed cannot be
shipped, and prose cannot be kept true against it.

So this runs the measurement deliberately, on a machine chosen to be quiet,
and writes the result beside the other frozen datasets the figures read:
``tests/speed_accuracy.json``, ``tests/timing_other_codes.json``, and
``tests/prem_speed_accuracy.json``.  The notebook then draws it without
timing anything.

    python tests/measure_performance_scaling.py > tests/performance_scaling.json

The three routes are the ones a user can actually choose between: a Python
loop over scalar calls, one batched call through {\tt NumPy}, and the same
batched call through the compiled kernel.
"""

import json
import os
import platform
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'src'))

import globaldefs as gd                                       # noqa: E402
import hamiltonians3nu                                        # noqa: E402
import oscprob3nu                                             # noqa: E402
import fastkernels                                            # noqa: E402

SIZES = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000]
L_KM = 1300.0


def best_of(func, repeat=5):
    r"""The minimum over `repeat` runs.

    The minimum, not the mean: timing noise is one-sided, so the fastest run
    is the one least polluted by whatever else the machine was doing.
    """
    best = float('inf')
    for _ in range(repeat):
        t0 = time.perf_counter()
        func()
        best = min(best, time.perf_counter() - t0)
    return best


def main():
    h_vac = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, gd.DCP_NO_BF,
        gd.D21_NO_BF, gd.D31_NO_BF)
    baseline = L_KM*gd.CONV_KM_TO_INV_EV

    loop, array, numba = [], [], []
    for n in SIZES:
        stack = hamiltonians3nu.hamiltonian_3nu_matter(
            h_vac, np.logspace(-1.0, 1.5, n)*1.0e9, gd.VCC_EARTH_CRUST)

        fastkernels.USE_NUMBA = False
        # The batched routes are quick enough to afford many repeats, and
        # they need them: at one repeat the cost per probability came out
        # NON-monotonic in the stack size, which is contention rather than
        # any property of the code.  The loop is slow, so it gets fewer.
        loop.append(best_of(lambda: [oscprob3nu.probabilities_3nu(h, baseline)
                                     for h in stack], repeat=3))
        array.append(best_of(lambda: oscprob3nu.probabilities_3nu(stack,
                                                                 baseline),
                             repeat=15))
        if fastkernels.HAVE_NUMBA:
            fastkernels.USE_NUMBA = True
            oscprob3nu.probabilities_3nu(stack, baseline)   # warm the compiler
            numba.append(best_of(lambda: oscprob3nu.probabilities_3nu(
                stack, baseline), repeat=15))
            fastkernels.USE_NUMBA = False
        else:
            numba.append(None)
    fastkernels.USE_NUMBA = fastkernels.HAVE_NUMBA

    record = {
        'generated_by': 'tests/measure_performance_scaling.py',
        'note': ('Frozen so that the performance figure does not move when it '
                 'is redrawn.  Re-run this deliberately, on a quiet machine, '
                 'and update the numbers quoted in the text with it.'),
        'baseline_km': L_KM,
        'energies_gev': '10^-1 to 10^1.5, log-spaced, n of them',
        'flavors': 3,
        'environment': {
            'python': platform.python_version(),
            'numpy': np.__version__,
            'platform': platform.platform(),
            'have_numba': bool(fastkernels.HAVE_NUMBA),
        },
        'sizes': SIZES,
        'seconds': {'loop': loop, 'array': array, 'numba': numba},
    }
    json.dump(record, sys.stdout, indent=1)
    sys.stdout.write('\n')


if __name__ == '__main__':
    main()
