# -*- coding: utf-8 -*-
r"""What the machine was doing, and whether a run may be believed.

Absolute times do not transfer between machines, so every artifact records the
one it came from.  That much is bookkeeping.  The harder problem is a run that
*looks* fine and is not: timings on this hardware swing by a quarter between
sessions, and the audit that prompted this pipeline found a published pair of
numbers whose ratio was noise dressed as a result.

So a run has to earn belief.  Three layers, all recorded:

* the environment, captured rather than assumed;
* a rejection rule --- a block whose coefficient of variation is too large, or
  whose median sits too far above its minimum, means the machine was busy and
  no artifact is written;
* a canary, this library's own kernel timed at the start, middle and end of a
  session, so that drift across an hours-long sweep is visible afterwards
  rather than invisible forever.
"""

import json
import os
import platform
import statistics
import subprocess
import time

#: A block whose spread exceeds this was measured on a busy machine.
MAX_BLOCK_CV = 0.10

#: The median of a well-behaved set of blocks sits close to its minimum.
MAX_MEDIAN_OVER_MIN = 1.25

#: Canary drift beyond this across a session invalidates its artifacts.
MAX_CANARY_DRIFT = 0.05

#: The canary's fixed workload.  Never change it: its whole value is being
#: comparable to the same number taken in another session.
CANARY_POINTS = 10000


def _first_line(command):
    try:
        out = subprocess.run(command, shell=True, capture_output=True,
                             text=True, timeout=10).stdout.strip()
        return out.splitlines()[0].strip() if out else 'unknown'
    except Exception:                                    # noqa: BLE001
        return 'unknown'


def environment():
    r"""Returns what this run can be attributed to.

    Thread counts are recorded per code by the caller; what is captured here is
    the machine and the policy, because a cross-code speed claim is meaningless
    without both.
    """
    load1, load5, load15 = os.getloadavg()
    return {
        'cpu': _first_line("grep -m1 'model name' /proc/cpuinfo | cut -d: -f2"),
        'cores': os.cpu_count(),
        'governor': _first_line(
            'cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor'),
        'load_average': {'1min': load1, '5min': load5, '15min': load15},
        'python': platform.python_version(),
        'cxx': _first_line('g++ --version'),
        'thread_policy': 'natural per code; counts recorded per series',
        'captured_at': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
    }


def canary():
    r"""Returns microseconds per probability for a fixed reference workload.

    This library's three-flavor kernel on `CANARY_POINTS` energies, timed the
    same way every time.  It is not a result and belongs in no figure; it
    exists so two artifacts can be checked for comparability before their
    numbers are allowed into one sentence.
    """
    import numpy as np
    import globaldefs as gd
    import hamiltonians3nu
    import oscprob3nu

    vacuum = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF,
        gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF)
    energies = np.logspace(-1.0, 1.0, CANARY_POINTS)*1.0e9
    stack = hamiltonians3nu.hamiltonian_3nu_matter(
        vacuum, energies, gd.VCC_EARTH_CRUST)
    baseline = 1300.0*gd.CONV_KM_TO_INV_EV

    oscprob3nu.probabilities_3nu(stack, baseline)          # warm
    best = float('inf')
    for _ in range(5):
        start = time.perf_counter()
        oscprob3nu.probabilities_3nu(stack, baseline)
        best = min(best, time.perf_counter() - start)
    return best/CANARY_POINTS*1.0e6


def admissible(blocks):
    r"""Returns ``(ok, why)`` for a set of per-block timings.

    A run that fails this writes no artifact.  Refusing to record is the point:
    an unbelievable number that reaches a file will eventually reach a figure.
    """
    if len(blocks) < 2:
        return False, 'fewer than two blocks'
    mean = statistics.fmean(blocks)
    if mean <= 0.0:
        return False, 'non-positive mean'
    cv = statistics.stdev(blocks)/mean
    if cv > MAX_BLOCK_CV:
        return False, 'block CV %.3f exceeds %.2f' % (cv, MAX_BLOCK_CV)
    ratio = statistics.median(blocks)/min(blocks)
    if ratio > MAX_MEDIAN_OVER_MIN:
        return False, ('median/min %.3f exceeds %.2f'
                       % (ratio, MAX_MEDIAN_OVER_MIN))
    return True, 'ok'


def admissible_stats(mean, sd, minimum, n):
    r"""Returns ``(ok, why)`` for one cell, from the summary it stored.

    The same judgement as :func:`admissible`, made from ``us_per_point``
    rather than from the raw blocks, so it can be applied to an artifact
    after the fact.  ``mean/min`` stands in for ``median/min``; both ask
    whether the typical block sits close to the best one, and the mean is
    the more conservative of the two because an outlier moves it further.

    This is what should decide whether a timing is believable.  The canary
    cannot: it runs this library on every thread the machine has, so it
    reports the state of a twelve-thread workload and then passes verdict
    on cells that were pinned to one thread, or that were compiled code
    linking no threading library at all.  A cell knows its own spread.
    """
    if not n or n < 2:
        return False, 'fewer than two blocks'
    if not mean or mean <= 0.0:
        return False, 'non-positive mean'
    cv = sd/mean
    if cv > MAX_BLOCK_CV:
        return False, 'block CV %.3f exceeds %.2f' % (cv, MAX_BLOCK_CV)
    if minimum and minimum > 0.0:
        ratio = mean/minimum
        if ratio > MAX_MEDIAN_OVER_MIN:
            return False, ('mean/min %.3f exceeds %.2f'
                           % (ratio, MAX_MEDIAN_OVER_MIN))
    return True, 'ok'


def canary_drift(start, mid, end):
    r"""Returns ``(ok, why)`` for a session's three canary readings.

    NO LONGER ENFORCED.  `run_all` records the canary readings but does not
    gate on this verdict: `MAX_CANARY_DRIFT` is 5 per cent while the canary's
    own repeatability on the development machine is 10.7 per cent across four
    consecutive best-of-five readings on an idle machine, so a before/after
    pair clears the threshold by chance most of the time.  It rejected eleven
    of twelve runs. Kept because it is still the right shape for the test;
    it needs a canary whose spread is narrower than the drift being tested
    for, which means a median of several readings rather than one.
    """
    lo, hi = min(start, mid, end), max(start, mid, end)
    drift = (hi - lo)/lo if lo > 0 else float('inf')
    if drift > MAX_CANARY_DRIFT:
        return False, 'canary drifted %.1f%% across the session' % (drift*100)
    return True, 'ok'


if __name__ == '__main__':
    import sys
    sys.path.insert(0, 'src')
    env = environment()
    print(json.dumps(env, indent=2))
    value = canary()
    print('canary: %.4f us per probability over %d points'
          % (value, CANARY_POINTS))
    ok, why = admissible([value, value*1.01, value*0.99])
    print('admissible(synthetic well-behaved blocks): %s (%s)' % (ok, why))
    ok, why = admissible([value, value*2.0, value*0.9])
    print('admissible(synthetic busy machine):        %s (%s)' % (ok, why))
