# -*- coding: utf-8 -*-
r"""The benchmark harness for the Python codes -- the mirror of bench.hpp.

bench.hpp owns main(), every clock and the statistics on the C++ side; this
file owns them on the Python side, and it is the ONLY place Python timing
happens.  An adapter under ``tests/bench/adapters/`` supplies physics and
nothing else: its source may not name a clock, which is checked at load the
way the C++ side's duplicate main() is a link error rather than a
convention.

Same protocols, same statistics, same JSON shape as bench.hpp:

* ``AMORTIZED`` -- the timed region is a delta_CP scan of ``--steps`` steps,
  each step = configure(dcp) + evaluate(whole grid), preceded by one untimed
  warm-up pass.  Everything invariant under the scan ran in ``setup()``,
  before any clock.
* ``THROUGHPUT`` -- one request for the whole grid, autoranged to a block of
  at least 0.05 s, repeated ``--samples`` times after one untimed warm-up.
  ``reset()`` is called before every ``evaluate()``, because the repetition
  is measurement scaffold rather than a workload -- nobody asks the same
  question five thousand times against unchanged state -- and the mean is
  only the cost of a request if every repetition costs what the first did.
  AMORTIZED never calls it: there the caching is the thing being measured.

An adapter is a class exposed as the module attribute ``ADAPTER``, with::

    name                     # str
    capabilities() -> dict   # batches_energy, batches_zenith, batch_symbol,
                             # knob_name, knob_domain
    setup(problem)           # hoisted; never timed
    configure(dcp)           # timed, once per scan step
    evaluate() -> float      # timed; returns a checksum
    reset()                  # timed; makes the next evaluate cold
    probabilities() -> list  # untimed; the nu_mu row per grid point, in
                             # grid order: numu->nue, numu->numu, numu->nutau

Oscillation parameters come from :mod:`conversions` (which reads the
manifest) and the grids come from the manifest too, so the Python side
cannot drift from what bench.hpp compiled in.

Usage::

    python tests/bench/runner.py <code> [--protocol amortized|throughput]
        [--grid CHORD/12x1] [--knob N] [--samples 30] [--steps 25]
        [--n-energies N] [--n-zenith N] [--json out.json]
"""

import argparse
import importlib
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, 'adapters'))

import conversions  # noqa: E402

#: code name on the command line -> adapter module under adapters/.
CODES = {
    'nuoscprobexact': 'nuoscprobexact',
    'nusquids':       'nusquids',
    'nucraft':        'nucraft',
}

#: Substrings no adapter source may contain: the harness times, the adapter
#: does not.  The C++ side enforces its half by owning main(); this is the
#: Python half of the same mechanism.
FORBIDDEN = ('perf_counter', 'time.time', 'monotonic', 'clock_gettime',
             'process_time', 'import timeit', 'time.clock')


class Problem(object):
    r"""The mirror of ``bench::Problem`` -- one shared parameter set.

    The oscillation parameters are filled from :mod:`conversions`, never
    typed here, for the same reason bench.hpp takes them from the generated
    ``conversions.h``: a second copy is how "every code gets the same
    numbers" stops being true.
    """

    def __init__(self):
        params = conversions.oscillation_parameters()
        self.energies_gev = []
        self.costhz = []                  # empty for constant-density grids
        self.L_km = 1300.0                # as bench.hpp's Problem defaults
        self.density = 3.0                # them; overridden from the
        self.ye = 0.5                     # manifest grid entry when present
        self.knob = 0
        self.s12sq = params['s12sq']
        self.s13sq = params['s13sq']
        self.s23sq = params['s23sq']
        self.dcp = conversions.for_code('generic')['dcp_rad']
        self.dm21 = params['dmsq21_ev2']
        self.dm31 = params['dmsq31_ev2']

    def points(self):
        return len(self.energies_gev)*max(len(self.costhz), 1)


def logspace(lo, hi, n):
    r"""bench.hpp's logspace, formula for formula."""
    return [lo if n == 1 else lo*math.pow(hi/lo, float(i)/(n - 1))
            for i in range(n)]


def linspace(lo, hi, n):
    r"""bench.hpp's linspace, formula for formula."""
    return [lo if n == 1 else lo + (hi - lo)*float(i)/(n - 1)
            for i in range(n)]


def build_problem(grid, knob, n_e, n_z):
    r"""Returns the Problem for `grid`, read from the manifest.

    bench.hpp compiles the same four grids in; here they are read from
    ``manifest.json`` so the two harnesses share one definition.  The
    N-sweep's point count is a runtime dial in both (``--n-energies``),
    defaulting to 1000 as bench.hpp does.
    """
    manifest = json.load(open(os.path.join(HERE, 'manifest.json')))
    grids = manifest['grids']
    if grid not in grids:
        grid = 'CONST/N-sweep'
    entry = grids[grid]

    p = Problem()
    p.knob = knob
    p.L_km = float(entry.get('L_km', p.L_km))
    p.density = float(entry.get('density_g_cm3', p.density))

    e = entry['energies_gev']
    lo, hi = e['log']
    n = e['n'] if isinstance(e['n'], int) else 1000
    p.energies_gev = logspace(lo, hi, n_e if n_e else n)

    z = entry.get('costhz')
    if isinstance(z, list):
        p.costhz = list(z)
    elif isinstance(z, dict):
        lo, hi = z['lin']
        p.costhz = linspace(lo, hi, n_z if n_z else z['n'])
    return p


def reduce_stats(values):
    r"""Mean, sample sd, min and n -- bench.hpp's reduce, guard included."""
    n = len(values)
    if n == 0:
        return {'mean': 0.0, 'sd': 0.0, 'min': 0.0, 'n': 0}
    mean = sum(values)/n
    var = sum(x*x for x in values)/n - mean*mean
    sd = math.sqrt(var*n/(n - 1 if n - 1 > 0 else 1)) if var > 0.0 else 0.0
    return {'mean': mean, 'sd': sd, 'min': min(values), 'n': n}


def amortized(driver, problem, samples, steps):
    r"""The scan is the timed region; the mirror of bench.hpp's amortized."""
    d0, dd = problem.dcp, 0.2/steps
    sink = 0.0
    for k in range(steps):            # untimed warm-up pass
        driver.configure(d0 + k*dd)
        sink += driver.evaluate()
    per_point = []
    for _ in range(samples):
        t0 = time.perf_counter()
        for k in range(steps):
            driver.configure(d0 + k*dd)
            sink += driver.evaluate()
        t1 = time.perf_counter()
        per_point.append((t1 - t0)/(steps*problem.points())*1.0e6)
    return reduce_stats(per_point), sink


def throughput(driver, problem, samples, min_block):
    r"""One request for the whole grid, repeated; bench.hpp's throughput."""
    driver.configure(problem.dcp)
    driver.reset()
    sink = driver.evaluate()          # untimed warm-up
    reps = 1                          # autorange to a stable block
    while True:
        t0 = time.perf_counter()
        for _ in range(reps):
            driver.reset()
            sink += driver.evaluate()
        dt = time.perf_counter() - t0
        if dt >= min_block or reps > (1 << 24):
            break
        reps = int(reps*(min_block/dt)*1.25) + 1 if dt > 0.0 else reps*8
    per_point = []
    for _ in range(samples):
        t0 = time.perf_counter()
        for _ in range(reps):
            driver.reset()
            sink += driver.evaluate()
        t1 = time.perf_counter()
        per_point.append((t1 - t0)/(float(reps)*problem.points())*1.0e6)
    return reduce_stats(per_point), sink


def load_adapter(code):
    r"""Imports the adapter and checks that its source names no clock."""
    module = importlib.import_module(CODES[code])
    source = open(module.__file__).read()
    hits = [t for t in FORBIDDEN if t in source]
    if hits:
        raise SystemExit(
            'adapter %s names a clock (%s); the harness times, the adapter '
            'does not' % (module.__file__, ', '.join(hits)))
    adapter = module.ADAPTER()
    # The Python mirror of the C++ side's link error: an adapter with no
    # reset() cannot be driven under THROUGHPUT without silently reusing
    # whatever the previous repetition left behind.
    if not hasattr(adapter, 'reset'):
        raise SystemExit('adapter %s defines no reset(); see runner.__doc__'
                         % module.__file__)
    # The Python mirror of bench.hpp's link error: an adapter that cannot be
    # checked for accuracy should not be usable for speed either.
    if not hasattr(adapter, 'probabilities'):
        raise SystemExit('adapter %s defines no probabilities(); see '
                         'runner.__doc__' % module.__file__)
    return adapter


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('code', choices=sorted(CODES))
    ap.add_argument('--protocol', default='amortized',
                    choices=['amortized', 'throughput', 'accuracy'])
    ap.add_argument('--grid', default='CHORD/12x1')
    ap.add_argument('--knob', type=int, default=0)
    ap.add_argument('--samples', type=int, default=30)
    ap.add_argument('--steps', type=int, default=25)
    ap.add_argument('--n-energies', type=int, default=0)
    ap.add_argument('--n-zenith', type=int, default=0)
    ap.add_argument('--json', dest='out', default='')
    args = ap.parse_args(argv)

    problem = build_problem(args.grid, args.knob, args.n_energies,
                            args.n_zenith)
    driver = load_adapter(args.code)
    driver.setup(problem)
    cap = driver.capabilities()

    # ACCURACY is untimed by construction, exactly as in bench.hpp: configure
    # once, ask for the probabilities, print them.  No clock is read on this
    # path, so a number produced here can never be mistaken for a speed, and
    # it can be produced on a machine that is not idle.
    if args.protocol == 'accuracy':
        driver.configure(problem.dcp)
        probs = driver.probabilities()
        record = {
            'code': driver.name,
            'protocol': {'name': 'accuracy', 'grid': args.grid},
            'knob': {cap.get('knob_name') or 'none': args.knob},
            'conventions': 'own-reference',
            'profile_basis': 'continuous',
            'environment': (driver.environment()
                            if hasattr(driver, 'environment') else {}),
            'channels': ['numu->nue', 'numu->numu', 'numu->nutau'],
            'n_points': problem.points(),
            'probabilities': list(probs),
        }
        text = json.dumps(record, indent=2)
        print(text)
        if args.out:
            with open(args.out, 'w') as handle:
                handle.write(text + '\n')
        return

    if args.protocol == 'throughput':
        stats, sink = throughput(driver, problem, args.samples, 0.05)
    else:
        stats, sink = amortized(driver, problem, args.samples, args.steps)

    record = {
        'code': driver.name,
        'protocol': {'name': args.protocol, 'grid': args.grid},
        'knob': {cap.get('knob_name') or 'none': args.knob},
        'conventions': 'own-reference',
        'profile_basis': 'continuous',
        'environment': (driver.environment()
                        if hasattr(driver, 'environment') else {}),
        'n_points': problem.points(),
        'us_per_point': stats,
        'batched': {'energy': bool(cap.get('batches_energy')),
                    'zenith': bool(cap.get('batches_zenith')),
                    'symbol': cap.get('batch_symbol', '')},
        'checksum': sink,
    }
    text = json.dumps(record, indent=2)
    print(text)
    if args.out:
        with open(args.out, 'w') as handle:
            handle.write(text + '\n')


if __name__ == '__main__':
    main()
