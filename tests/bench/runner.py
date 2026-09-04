# -*- coding: utf-8 -*-
r"""The benchmark harness for the Python codes -- the mirror of bench.hpp.

bench.hpp owns main(), every clock and the statistics on the C++ side; this
file owns them on the Python side, and it is the ONLY place Python timing
happens.  An adapter under ``tests/bench/adapters/`` supplies physics and
nothing else: its source may not name a clock, which is checked at load the
way the C++ side's duplicate main() is a link error rather than a
convention.

Same protocols, same statistics, same JSON shape as bench.hpp:

* ``AMORTIZED`` -- the timed region is a delta_CP scan, each step =
  configure(dcp) + evaluate(whole grid).  The step count autoranges up from
  ``--steps`` until the block lasts ``--min-block``, so the scan is finer on
  a cheap code rather than shorter; the first pass is the untimed warm-up.
  Everything invariant under the scan ran in ``setup()``, before any clock.
* ``THROUGHPUT`` -- one request for the whole grid, autoranged to a block of
  at least ``--min-block``, repeated ``--samples`` times after one untimed
  warm-up.
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
import hashlib
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
    'second_order':   'second_order',
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
        # Which parameter the amortized scan turns; bench.hpp's Problem::scan.
        self.scan = 'dcp'

    def scan_base(self):
        return self.dm31 if self.scan == 'dmsq31' else self.dcp

    def scan_width(self):
        r"""How far the scan travels.

        Two per cent of Dmsq31 is far inside its measured uncertainty and far
        outside anything a cache could call unchanged; 0.2 rad is what the
        delta_CP scan has always used.
        """
        return 0.02*self.dm31 if self.scan == 'dmsq31' else 0.2

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


def block_samples(block, min_block, lo, hi):
    r"""Same wall-clock budget for every cell; the mirror of bench.hpp's.

    A cell cheap enough to have been autoranged gets `hi` blocks; one whose
    single pass already outlasts the budget gets `lo`.  The rule reads the
    clock and nothing else -- no code is named in it, so none can be
    favoured by it.
    """
    if block <= 0.0:
        return hi
    return max(lo, min(hi, int(hi*min_block/block)))


def amortized(driver, problem, samples, steps, min_block, max_samples):
    r"""The scan is the timed region; the mirror of bench.hpp's amortized.

    The step count autoranges until the timed region reaches `min_block`.  A
    fixed count gave the cheapest cells a block of well under a millisecond,
    and what came back was the clock's spread rather than the code's.

    The scan covers the same delta_CP interval whatever the step count, so a
    longer block is a finer scan and never a repeated one: no code is handed
    the same delta twice, so none can look fast for having cached it.
    """
    d0, span = problem.scan_base(), problem.scan_width()
    sink = 0.0
    total = steps if steps > 0 else 1
    # Warm-up on its own and discarded; see bench.hpp.  Folding it into the
    # autorange made a code with an expensive first call calibrate on that
    # call and then time a block far below min_block.
    dd0 = 0.2/total
    for k in range(total):
        driver.configure(d0 + k*dd0)
        sink += driver.evaluate()
    while True:
        dd = span/total
        t0 = time.perf_counter()
        for k in range(total):
            driver.configure(d0 + k*dd)
            sink += driver.evaluate()
        dt = time.perf_counter() - t0
        if dt >= min_block or total > (1 << 24):
            break
        total = int(total*(min_block/dt)*1.25) + 1 if dt > 0.0 else total*8
    n = block_samples(dt, min_block, samples, max_samples)
    # The scan continues across blocks rather than restarting in each one, so
    # the cell sweeps the interval exactly once and no delta_CP is evaluated
    # twice.  Where one pass already outlasts a block -- nuCraft on the
    # oscillogram -- a block holds one step, and restarting would time the
    # same delta over and over.  See bench.hpp for the same reasoning.
    dd = span/(total*n)
    step = 0
    per_point = []
    for _ in range(n):
        t0 = time.perf_counter()
        for _k in range(total):
            driver.configure(d0 + step*dd)
            sink += driver.evaluate()
            step += 1
        t1 = time.perf_counter()
        per_point.append((t1 - t0)/(total*problem.points())*1.0e6)
    stats = reduce_stats(per_point)
    stats['reps'] = total
    return stats, sink


def throughput(driver, problem, samples, min_block, max_samples):
    r"""One request for the whole grid, repeated; bench.hpp's throughput."""
    driver.configure(problem.scan_base())
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
    n = block_samples(dt, min_block, samples, max_samples)
    per_point = []
    for _ in range(n):
        t0 = time.perf_counter()
        for _ in range(reps):
            driver.reset()
            sink += driver.evaluate()
        t1 = time.perf_counter()
        per_point.append((t1 - t0)/(float(reps)*problem.points())*1.0e6)
    stats = reduce_stats(per_point)
    stats['reps'] = reps
    return stats, sink


def thread_environment():
    r"""The threading environment of this process, for every code alike.

    ``environment()`` is optional on an adapter, so for five of the six
    external codes the artifact recorded ``{}`` and a cross-code speed claim
    had nothing to stand on but prose.  The bench documentation says the
    effective count is recorded per run; this is what makes that true.  It
    is gathered here rather than in the adapter for the same reason the
    clock is: a code cannot be trusted to report its own advantage.

    ``blas_threads`` matters for the codes driven from Python, whose real
    parallelism is whatever NumPy's BLAS decides and not anything the
    adapter sets.
    """
    env = {'cpu_count': os.cpu_count()}
    for var in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS'):
        env[var.lower()] = os.environ.get(var)
    try:
        import threadpoolctl
        env['blas_threads'] = [
            {'library': p.get('user_api') or p.get('internal_api'),
             'threads': p.get('num_threads')}
            for p in threadpoolctl.threadpool_info()]
    except Exception:                                        # noqa: BLE001
        env['blas_threads'] = 'threadpoolctl not installed'
    return env


def manifest_sha():
    r"""The hash of the manifest these numbers were measured under.

    The compiled harness stamps this into every artifact it writes; the
    Python adapters did not, so a hundred and fifty-seven artifacts --- all
    three Python codes, both protocols --- could not be tied to the
    manifest that pins their versions.  That manifest is the answer to the
    objection about taking whatever is on GitHub today, and half the
    artifacts were not pointing at it.
    """
    path = os.path.join(HERE, 'manifest.json')
    try:
        with open(path, 'rb') as handle:
            return hashlib.sha256(handle.read()).hexdigest()
    except OSError:
        return ''


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
    ap.add_argument('--steps', type=int, default=1)
    ap.add_argument('--min-block', type=float, default=0.25,
                    dest='min_block')
    ap.add_argument('--max-samples', type=int, default=100,
                    dest='max_samples')
    ap.add_argument('--scan', default='dcp', choices=('dcp', 'dmsq31'),
                    help='which parameter the amortized scan turns')
    # The knob is an integer and the tolerance dial was encoded as
    # rtol = 10**knob, which can express 1e-3 and 1e-4 and nothing between
    # them.  The paper's curve runs 3, 1, 0.3, 0.03, 0.01, 1e-3, 1e-4,
    # 3e-5, 1e-5 -- five of those nine are not powers of ten, so six of its
    # points were unreachable and the curve lost its whole loose end.
    ap.add_argument('--rtol', type=float, default=0.0,
                    help='explicit tolerance, overriding the knob encoding')
    ap.add_argument('--n-energies', type=int, default=0)
    ap.add_argument('--n-zenith', type=int, default=0)
    ap.add_argument('--json', dest='out', default='')
    args = ap.parse_args(argv)

    problem = build_problem(args.grid, args.knob, args.n_energies,
                            args.n_zenith)
    problem.scan = args.scan
    if args.rtol:
        problem.rtol_override = args.rtol
    driver = load_adapter(args.code)
    driver.setup(problem)
    cap = driver.capabilities()

    # ACCURACY is untimed by construction, exactly as in bench.hpp: configure
    # once, ask for the probabilities, print them.  No clock is read on this
    # path, so a number produced here can never be mistaken for a speed, and
    # it can be produced on a machine that is not idle.
    if args.protocol == 'accuracy':
        driver.configure(problem.scan_base())
        probs = driver.probabilities()
        record = {
            'code': driver.name,
            'protocol': {'name': 'accuracy', 'grid': args.grid},
            'knob': {cap.get('knob_name') or 'none': args.knob},
            'conventions': 'own-reference',
            'manifest_sha256': manifest_sha(),
            'scan': problem.scan,
            'profile_basis': 'continuous',
            'environment': dict(thread_environment(),
                                **(driver.environment()
                                   if hasattr(driver, 'environment') else {})),
            # From the adapter, not assumed: a code may answer for fewer
            # channels than the row holds, and claiming three when it
            # computed one would be the figure lying about its own data.
            'channels': cap.get('channels',
                                ['numu->nue', 'numu->numu', 'numu->nutau']),
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
        stats, sink = throughput(driver, problem, args.samples,
                                 args.min_block, args.max_samples)
    else:
        stats, sink = amortized(driver, problem, args.samples, args.steps,
                                args.min_block, args.max_samples)
    reps = stats.pop('reps', 0)

    record = {
        'code': driver.name,
        'protocol': {'name': args.protocol, 'grid': args.grid},
        'knob': {cap.get('knob_name') or 'none': args.knob},
        'conventions': 'own-reference',
        'manifest_sha256': manifest_sha(),
        'scan': problem.scan,
        'profile_basis': 'continuous',
        'environment': dict(thread_environment(),
                            **(driver.environment()
                               if hasattr(driver, 'environment') else {})),
        'n_points': problem.points(),
        'us_per_point': stats,
        'timing': {'min_block_s': args.min_block, 'block_reps': reps,
                   'block_seconds': stats['mean']*1e-6*reps*problem.points(),
                   'samples': stats['n']},
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
