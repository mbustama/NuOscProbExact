# -*- coding: utf-8 -*-
r"""The run matrix: which code runs which grid, under which protocol, at which
knob, and what the artifact is called.

This file exists because the matrix was undefined.  A cold reader given this
directory could reach the harness and the adapters and then had to *invent*
which cells to run, and an invented matrix is not reproducible: the next person
invents a different one and the two sets of numbers cannot be compared.  So the
matrix is written down, and it is the same object the run consumes.

What this does NOT do is decide fairness.  Every fairness decision lives
elsewhere and is merely obeyed here: the conventions in ``manifest.json``, the
references in ``reference.py``, the protocol definitions in ``bench.hpp`` and
``runner.py``.  This is the schedule.

Two properties are enforced rather than trusted, because both were violated by
the comparison this pipeline replaces:

* **A cell states its own configuration.**  Nothing is defaulted silently.  If
  NuFast-Earth is run with midpoint shell densities, the cell says so, because
  the same code with mean densities is a different measurement and the two may
  not share an axes.
* **A timed cell refuses an untrustworthy machine.**  The canary is read before
  and after, and if it drifted the artifacts are written with ``believable``
  false rather than quietly kept.

Usage::

    python tests/bench/run_all.py --dry-run        # print the matrix, run nothing
    python tests/bench/run_all.py --accuracy-only  # the untimed half
    python tests/bench/run_all.py                  # everything; needs an idle machine
"""

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
BIN = os.path.join(ROOT, '.bench-build', 'bin')
sys.path.insert(0, HERE)

#: Codes that can be posed each problem kind.  Getting this wrong is silent:
#: asked for a chord, a constant-density-only code answers a different question
#: and returns perfectly good numbers for it.
EARTH_CODES = ['NuFast-Earth', 'Prob3++', 'GLoBES', 'nuSQuIDS', 'nuCraft',
               'NuOscProbExact']
CONST_CODES = ['NuFast-LBL', 'NuFast-Earth', 'Prob3++', 'GLoBES', 'nuSQuIDS',
               'NuOscProbExact']

BINARY = {
    'NuFast-LBL':   'bench_nufast_lbl',
    'NuFast-Earth': 'bench_nufast_earth',
    'Prob3++':      'bench_prob3',
    'GLoBES':       'bench_globes',
}
ADAPTER = {
    'nuSQuIDS':       'nusquids',
    'nuCraft':        'nucraft',
    'NuOscProbExact': 'nuoscprobexact',
}

#: Codes whose spatial discretisation can be swept, and what dial reaches it.
#:
#: NuFast-LBL has no entry because it has no spatial cut at all -- it
#: propagates a constant density over one baseline, so there is nothing to
#: refine.  nuSQuIDS and nuCraft are deliberately absent for a different
#: reason: their discretisation is the node count of the table THIS harness
#: builds and hands them (200 per PREM shell), not a dial their own code
#: exposes.  Sweeping it would measure our input choice rather than their
#: algorithm, which is a legitimate question but a different one, and it would
#: be dishonest to present it beside dials the codes actually own.
DISCRETISATION_DIAL = {
    'NuFast-Earth':   'n_layers (separate; its knob is eigenvalue precision)',
    'Prob3++':        'n_shells (its knob)',
    'GLoBES':         'n_shells (its knob)',
    'NuOscProbExact': 'n_slabs_per_segment (its knob, positive values)',
}

#: The precision knob swept per code, for the accuracy axis.
KNOBS = {
    'NuFast-LBL':     [-1, 0, 1, 2, 3],
    'NuFast-Earth':   [-1, 0, 1, 2, 3],
    'Prob3++':        [1, 2, 4, 8, 16, 32, 64, 128, 256],
    'GLoBES':         [1, 2, 4, 8, 16, 32, 64, 128, 256],
    'nuSQuIDS':       [3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
    'nuCraft':        [2, 3, 4, 5, 6, 7, 8, 9, 10],
    'NuOscProbExact': [1, 2, 4, 8, 16, 32, 64, 128, 256, -3, -4, -5],
}

#: The knob a THROUGHPUT or oscillogram cell holds fixed.  Explicit, with the
#: measurement that justifies each entry, because deriving it went wrong three
#: times: pooling both grids picked settings that are meaningless on one of
#: them, and breaking ties on speed fitted the 35 per cent block noise of the
#: cheapest cells.
#:
#: The rule is: each code's best WORKING setup -- the setting that reaches its
#: best measured accuracy on that grid.  Accuracy is untimed and reliable, so
#: nothing here depends on a timing.  Where several settings tie, the most
#: exact is kept: a difference at the round-off floor is not an accuracy
#: difference, and pretending otherwise is how a code ends up timed at a
#: setting its author did not ask for.
#:
#: NuFast-LBL and NuFast-Earth keep -1 for exactly that reason.  On constant
#: density -1, 2 and 3 all sit at 1e-16; on the chord -1, 1, 2 and 3 all tie
#: at the discretisation floor.  -1 is never measurably worse, and it is the
#: exact-eigenvalue mode whose absence was the strongest objection raised
#: against the earlier comparison.  It stays.
THROUGHPUT_KNOB = {
    'NuFast-LBL':     -1,   # exact; ties with 3 at 7.6e-16, i.e. round-off
    'NuFast-Earth':   -1,   # exact; ties with 2 and 3 at 1e-16
    'Prob3++':        256,  # dial inert on constant density: one accuracy
    'GLoBES':         256,  # dial inert on constant density
    'nuSQuIDS':       12,   # dial inert on constant density
    'nuCraft':        10,   # no constant-density cells; it refuses that grid
    'NuOscProbExact': 256,  # dial inert on constant density
}

#: The oscillogram runs on the chord, where the dial is NOT inert, so two of
#: these differ from the constant-density choice -- and both differ by margins
#: far above round-off.
OSCILLOGRAM_KNOB = {
    'NuFast-Earth':   -1,   # ties with 1, 2, 3 at the 3.25e-6 floor; exact kept
    'Prob3++':        256,  # 3.25e-6, its best and reached only here
    'GLoBES':         256,  # 3.25e-6, likewise
    'nuSQuIDS':       12,   # 6.45e-9, its best
    'nuCraft':        6,    # 1.39e-6 against knob 10's 1.73e-6, AND 10x faster
    # 256 fixed slabs, and the tolerance dial deliberately NOT used here.
    #
    # On the twelve-point chord the tolerance dial wins: rtol=1e-5 converges
    # at 1024 slabs per segment and reaches 5.47e-8, against 8.75e-7 for the
    # fixed 256.  On this grid it does not run at all.  Measured, not
    # assumed: rtol=1e-5 and rtol=1e-4 both raise ValueError, unable to meet
    # a RELATIVE tolerance within the library's 1024-slab ceiling at the
    # grid points where the probability is near zero; and rtol=1e-3, which
    # does run, costs 429 us/point against the fixed 256's 106 -- four times
    # slower for sixteen times worse accuracy, because the adaptive search
    # pays for every trial subdivision below the one the worst point of ten
    # thousand demands.
    #
    # So the fixed dial IS this library's best working setup on the
    # oscillogram, and its 8.75e-7 is better than four of the five codes it
    # is timed against here.  Raising n_max to force rtol=1e-5 through would
    # need some 7000 slabs per segment for a precision no other code on this
    # grid is delivering, which would measure a self-imposed handicap rather
    # than the library.
    'NuOscProbExact': 256,
}

#: AMORTIZED sweeps the SAME knob domain the accuracy axis sweeps, so that
#: every timed point can be paired with the accuracy the code reached at that
#: exact setting.  That pairing is the speed-accuracy plane, and it is the
#: only fair way to compare these codes on speed.
#:
#: Pinning one knob per code, which this file did first, compares codes doing
#: different amounts of work for answers of different quality.  Measured at
#: the settings that version would have used: Prob3++, GLoBES and NuFast-Earth
#: all at 3.25e-6 on the chord, this library at 8.75e-7, and nuSQuIDS at
#: 2.86e-5 -- an order worse than the others, and four orders worse than the
#: 6.45e-9 it reaches at its tightest.  On constant density the spread was
#: eight orders.  A speed number is not comparable without the accuracy it
#: bought, and choosing that accuracy per code by hand is a thumb on the
#: scale.  The reader compares at equal accuracy instead.
#:
#: This is also what the previous generation of this comparison did, in
#: tests/speed_accuracy.json, and its defects are instructive: GLoBES and
#: Prob3++ were both recorded as exposing no dial when both expose n_shells,
#: and NuFast-LBL was swept 0..3 with its exact mode absent -- which is
#: objection LBL-3, visible in the shape of the artifact.


def cells(accuracy_only=False, speed_only=False):
    r"""Yields every cell of the matrix, as explicit dicts.

    A cell carries everything needed to reproduce it and everything needed to
    name it.  Nothing is implied.
    """
    out = []

    if not speed_only:
        for grid, codes in (('CONST/60E', CONST_CODES),
                            ('CHORD/12x1', EARTH_CODES)):
            for code in codes:
                for knob in KNOBS[code]:
                    out.append({'kind': 'accuracy', 'code': code, 'grid': grid,
                                'knob': knob, 'timed': False,
                                'shell_density': 'midpoint'})

    if not accuracy_only:
        # AMORTIZED: the cost inside a scan, with each code's caching working.
        # OSC/100x100 is here because a 12x1 grid hides NuFast-Earth's zenith
        # caching entirely -- that is objection Earth-3, and the grid exists
        # to expose it.
        # OSC/100x100 is tiered, and the reason is not only cost.
        #
        # That grid exists for one objection: a 12x1 chord hides NuFast-Earth's
        # caching across zenith angles, worth roughly eightfold per
        # probability, so the oscillogram is where that advantage is visible.
        # It is ONE comparison at one precision setting, not a sweep -- and it
        # could not be a sweep even if the time were free, because the
        # accuracy axis was only measured on CONST/60E and CHORD/12x1, so
        # every knob variant here would join to a null error and produce a
        # speed-accuracy point with no accuracy in it.
        #
        # Sweeping it anyway costs 114 hours of the 115-hour total, 99 per
        # cent of the run, for points that cannot enter the plane.  One knob
        # per code, at each code's best precision, stated.
        #
        # The sample count is NOT cut.  It was -- ten blocks here and three
        # for nuCraft -- on the reasoning that ten thousand grid points
        # already average away the per-point noise.  That reasoning confused
        # two different averages: many points per block makes the block's
        # MEAN precise, and says nothing about the spread BETWEEN blocks,
        # which is what the standard deviation reports and what three blocks
        # cannot estimate.  The harness now sets the count itself from the
        # measured block length, so a cell that can afford many blocks takes
        # many and one that cannot takes the floor of thirty.
        for code in EARTH_CODES:
            cell = {'kind': 'amortized', 'code': code,
                    'grid': 'OSC/100x100', 'knob': OSCILLOGRAM_KNOB[code],
                    'timed': True,
                    'tier': 'oscillogram: one knob, Earth-3',
                    'shell_density': 'midpoint'}
            # nuCraft needs its own tier on this grid, and the reason is a
            # measurement rather than a guess: the first attempt at this
            # matrix recorded NuFast-Earth at 0.2 s, GLoBES at 144, Prob3++
            # at 164, nuSQuIDS at 392 -- and nuCraft exceeded two hours for
            # the same cell and had to be killed.  Its cost per point is not
            # grid-independent, which the calibration on a twelve-point chord
            # had assumed; on ten thousand points it is some four times worse
            # than that extrapolation.
            #
            # It is not dropped: a code being slow IS the measurement, and
            # excluding it because it is slow would be the plainest
            # unfairness in this whole comparison.  It gets the same thirty
            # blocks as everything else -- roughly half an hour of wall clock
            # at this knob -- so its standard deviation means what every
            # other code's means, and only its timeout is raised.
            if code == 'nuCraft':
                cell.update({'timeout': 10800,
                             'tier': 'oscillogram: full samples, long timeout'})
            out.append(cell)

        for grid, codes in (('CHORD/12x1', EARTH_CODES),
                            ('CONST/60E', CONST_CODES)):
            for code in codes:
                for knob in KNOBS[code]:
                    out.append({'kind': 'amortized', 'code': code,
                                'grid': grid, 'knob': knob, 'timed': True,
                                'plane': True, 'shell_density': 'midpoint'})
                    # This library, and only this library, is also run on ONE
                    # thread.  The four compiled codes link no threading
                    # library at all -- verified with ldd, no libgomp, no
                    # threaded BLAS -- so they are single-threaded by
                    # construction, while this one spreads seventeen prange
                    # kernels over every core it is given.
                    #
                    # Forcing everyone to one thread would be the wrong fix:
                    # it strips this library's parallelism while leaving the
                    # others untouched, and it publishes a number no user of
                    # this library would ever see.  Reporting both thread
                    # counts and stopping there is also not enough, because a
                    # reader given a twelvefold gap cannot tell how much of it
                    # is the algorithm and how much is the core count.
                    #
                    # So both are measured.  The single-thread series is the
                    # algorithmic comparison, like for like against codes that
                    # have no choice; the parallel series is what a user
                    # actually gets.  Neither is the honest number on its own.
                    if code == 'NuOscProbExact':
                        out.append({'kind': 'amortized', 'code': code,
                                    'grid': grid, 'knob': knob, 'timed': True,
                                    'plane': True, 'threads': 1,
                                    'shell_density': 'midpoint'})

        # THROUGHPUT: what one request for N points costs, with every
        # repetition started afresh.  The N-sweep is the figure's x axis.
        for n in (1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000):
            for code in CONST_CODES:
                out.append({'kind': 'throughput', 'code': code,
                            'grid': 'CONST/N-sweep', 'knob': THROUGHPUT_KNOB[code],
                            'n_energies': n, 'timed': True,
                            'shell_density': 'midpoint'})
        # Objection LBL-1's control: the same code, the same entry point, one
        # energy per call, so that what batching buys is measured rather than
        # asserted.
        for n in (1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000):
            out.append({'kind': 'throughput', 'code': 'NuFast-LBL',
                        'grid': 'CONST/N-sweep', 'knob': THROUGHPUT_KNOB['NuFast-LBL'],
                        'n_energies': n, 'timed': True, 'loop': True,
                        'shell_density': 'midpoint'})
        # The discretisation sweep, for every code that has one.
        #
        # Objection Earth-1 asked for this of NuFast-Earth, but restricting it
        # there would have been unfair twice over: Prob3++ and GLoBES expose
        # the same dial and would have shown a single point beside a
        # competitor's curve, and THIS LIBRARY exposes it too and would have
        # been the single point beside three curves.  A discretisation dial is
        # the same kind of dial whoever owns it, so all four are swept.
        #
        # The dial is reached differently per code, which is why this is a
        # table rather than a loop: NuFast-Earth spends its knob on eigenvalue
        # precision and takes the shell count separately, while for the other
        # three the shell or slab count IS the knob.
        for n in (1, 4, 16, 64, 256):
            out.append({'kind': 'amortized', 'code': 'NuFast-Earth',
                        'grid': 'CHORD/12x1', 'knob': THROUGHPUT_KNOB['NuFast-Earth'],
                        'n_layers': n, 'timed': True, 'sweep': 'discretisation',
                        'shell_density': 'midpoint'})
            for code in ('Prob3++', 'GLoBES', 'NuOscProbExact'):
                out.append({'kind': 'amortized', 'code': code,
                            'grid': 'CHORD/12x1', 'knob': n, 'timed': True,
                            'sweep': 'discretisation',
                            'shell_density': 'midpoint'})
    return out


def artifact_name(cell):
    r"""A cell's artifact name, derived from the cell and nothing else."""
    bits = [cell['kind'], cell['code'].replace('+', 'p').replace('-', '_'),
            cell['grid'].replace('/', '_'), 'knob%s' % cell['knob']]
    if cell.get('sweep'):
        bits.append(cell['sweep'])
    if cell.get('n_energies'):
        bits.append('nE%d' % cell['n_energies'])
    if cell.get('n_layers'):
        bits.append('L%d' % cell['n_layers'])
    if cell.get('loop'):
        bits.append('looped')
    if cell.get('threads') == 1:
        bits.append('1thread')
    return '_'.join(bits) + '.json'


def command(cell):
    r"""The exact command line a cell runs.  Printed by --dry-run."""
    code = cell['code']
    if code in BINARY:
        exe = os.path.join(BIN, BINARY[code]
                           + ('_accuracy' if cell['kind'] == 'accuracy'
                              and code != 'GLoBES' else ''))
        cmd = [exe]
    else:
        cmd = [sys.executable, os.path.join(HERE, 'runner.py'), ADAPTER[code]]
    cmd += ['--protocol', cell['kind'], '--grid', cell['grid'],
            '--knob', str(cell['knob'])]
    if cell.get('samples'):
        cmd += ['--samples', str(cell['samples'])]
    if cell.get('steps'):
        cmd += ['--steps', str(cell['steps'])]
    if cell.get('n_energies'):
        cmd += ['--n-energies', str(cell['n_energies'])]
    if cell.get('n_layers') and code in BINARY:
        cmd += ['--n-layers', str(cell['n_layers'])]
    if cell.get('loop') and code in BINARY:
        cmd += ['--loop']
    if cell.get('shell_density') == 'mean' and code == 'NuFast-Earth':
        cmd += ['--mean-density']
    return cmd


def join_plane(outdir):
    r"""Pairs each timed point with the accuracy reached at the same knob.

    This is the speed-accuracy plane, and it is emitted in the shape the
    previous generation of this comparison used --- one series per code,
    carrying the name of its dial and a list of points, each point holding
    both the cost and the error.  That shape was right; what was wrong was
    its contents, and the artifact showed it: GLoBES and Prob3++ were both
    recorded as exposing no dial when each exposes a shell count, and
    NuFast-LBL was swept from zero to three Newton steps with its exact mode
    -- the setting objection LBL-3 is about -- simply absent.

    Speed and accuracy are measured by different runs on different machines'
    terms: accuracy is untimed and indifferent to load, speed needs a quiet
    machine.  Joining them on (code, grid, knob) is what lets each be taken
    when it can be taken, instead of forcing both into one window.
    """
    accuracy = {}
    for grid, name in (('CONST/60E', 'accuracy_const.json'),
                       ('CHORD/12x1', 'accuracy_chord.json')):
        path = os.path.join(HERE, name)
        if os.path.exists(path):
            accuracy[grid] = json.load(open(path))

    plane = {'schema': 'bench-speed-accuracy/1',
             'generated_by': 'tests/bench/run_all.py --join',
             'manifest_sha256': manifest_sha(),
             'note': ('Speed against accuracy, per code, over that code\'s own '
                      'precision dial.  Every point pairs a timed measurement '
                      'with the accuracy the SAME setting reached against that '
                      'code\'s own 50-digit reference on the continuous PREM. '
                      'Comparing two codes means comparing at equal accuracy, '
                      'not at a knob somebody chose for them.'),
             'series': {}}

    for cell in cells():
        if not cell.get('plane'):
            continue
        path = os.path.join(outdir, artifact_name(cell))
        if not os.path.exists(path):
            continue
        timed = json.load(open(path))
        grid, code, knob = cell['grid'], cell['code'], str(cell['knob'])
        acc = accuracy.get(grid, {}).get('series', {}).get(code, {})
        entry = acc.get('by_knob', {}).get(knob, {})
        key = '%s|%s%s' % (code, grid,
                           '|1thread' if cell.get('threads') == 1 else '')
        series = plane['series'].setdefault(
            key, {'code': code, 'grid': grid,
                  'dial': DISCRETISATION_DIAL.get(code, 'precision knob'),
                  'threads': cell.get('threads',
                                       timed.get('environment', {})
                                       .get('numba_threads', 1)),
                  'points': []})
        series['points'].append({
            'knob': cell['knob'],
            'us_per_probability': timed['us_per_point']['mean'],
            'us_sd': timed['us_per_point']['sd'],
            'max_abs_error': entry.get('max_abs_deviation'),
            'accuracy_measured': 'max_abs_deviation' in entry,
        })
    for series in plane['series'].values():
        series['points'].sort(key=lambda q: q['knob'])
    return plane


def manifest_sha():
    with open(os.path.join(HERE, 'manifest.json'), 'rb') as handle:
        return hashlib.sha256(handle.read()).hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--accuracy-only', action='store_true')
    ap.add_argument('--speed-only', action='store_true')
    ap.add_argument('--outdir', default=os.path.join(HERE, 'artifacts'))
    # A full run is hours.  The first execution of it should never BE the
    # full run: --limit takes the first N cells so that artifact writing, the
    # canary bracketing and every command line can be proven end to end for
    # the price of a minute.
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--only', default='',
                    help='comma-separated substrings; keep only cells whose '
                         'artifact name contains one of them')
    ap.add_argument('--force', action='store_true',
                    help='recompute cells whose artifact already exists')
    ap.add_argument('--join', action='store_true',
                    help='join existing artifacts into the speed-accuracy '
                         'plane and exit, measuring nothing')
    args = ap.parse_args(argv)

    if args.join:
        plane = join_plane(args.outdir)
        path = os.path.join(HERE, 'speed_accuracy_plane.json')
        with open(path, 'w') as handle:
            json.dump(plane, handle, indent=2, sort_keys=True)
            handle.write('\n')
        n = sum(len(v['points']) for v in plane['series'].values())
        print('%d series, %d points -> %s'
              % (len(plane['series']), n, os.path.relpath(path, ROOT)))
        # And the three files the paper's figures actually draw from.  The
        # plane above is the whole measurement; these are its projections in
        # the shapes the notebook plots, regenerated here so that no figure
        # is fed by a file with no run behind it.
        import emit_figures
        emit_figures.main(args.outdir)
        return

    matrix = cells(args.accuracy_only, args.speed_only)
    if args.only:
        # Substring match on the artifact name, so a single cell can be
        # re-measured without re-running its neighbours.  --limit takes the
        # first N of the matrix, which is matrix order and not cost order:
        # the first six happen to be the six most expensive cells there are.
        wanted = [w.strip() for w in args.only.split(',') if w.strip()]
        matrix = [c for c in matrix
                  if any(w in artifact_name(c) for w in wanted)]
    if args.limit:
        matrix = matrix[:args.limit]
    timed = [c for c in matrix if c['timed']]
    print('%d cells: %d untimed, %d timed' % (len(matrix),
                                              len(matrix) - len(timed),
                                              len(timed)))
    if args.dry_run:
        for cell in matrix:
            print('  %-40s %s' % (artifact_name(cell),
                                  ' '.join(os.path.basename(x)
                                           for x in command(cell))))
        return

    if timed:
        import machine
        env = machine.environment()
        if env['governor'] != 'performance':
            print('WARNING: governor is %r, not performance.  Timed cells will '
                  'be recorded with the canary and may be rejected.'
                  % env['governor'])
        # Discard one reading first.  The canary's own warm-up -- numba
        # compiling, and the CPU climbing out of idle -- makes the FIRST call
        # in a fresh process run about twelve per cent slow, while every
        # reading after it agrees to under two.  Comparing that cold reading
        # against a warm one at the end is a guaranteed drift far above the
        # five per cent threshold, so every run would have been marked
        # unbelievable for a reason that has nothing to do with the machine.
        # Measured on an idle box at governor=performance: 0.0446, then
        # 0.0392, 0.0399, 0.0398, 0.0398, 0.0400.
        machine.canary()
        before = machine.canary()
        print('canary before: %.4f us/probability (after one discarded '
              'warm-up reading)' % before)

    os.makedirs(args.outdir, exist_ok=True)
    written = 0
    skipped = 0
    failures = []
    measured_here = []
    for cell in matrix:
        path = os.path.join(args.outdir, artifact_name(cell))

        # Resume.  A run that died at cell 109 of 328 should not repeat the
        # 108 that succeeded -- the first attempt lost two hours to a single
        # cell and had nothing to show for the rest.
        #
        # But existence is not enough.  The artifacts are written with a plain
        # fopen, so a run killed mid-write leaves a truncated file, and a
        # resume that trusted os.path.exists would skip it forever and carry a
        # corrupt cell into the figures.  Parse it before believing it.
        if os.path.exists(path) and not args.force:
            try:
                with open(path) as handle:
                    got = json.load(handle)
                intact = 'code' in got and ('probabilities' in got
                                            or 'us_per_point' in got)
            except (ValueError, OSError):
                intact = False
            if intact:
                skipped += 1
                written += 1
                continue
            print('%-46s %-7s  truncated or unparseable; recomputing'
                  % (artifact_name(cell), 'REDO'), flush=True)

        cmd = command(cell) + ['--json', path]
        started = time.time()
        env = dict(os.environ)
        if cell.get('threads') == 1:
            # numba reads this at import, so it must be in the child's
            # environment rather than set from inside the adapter.
            env['NUMBA_NUM_THREADS'] = '1'
            env['OMP_NUM_THREADS'] = '1'

        # One cell may not kill the run.  It did: nuCraft on the oscillogram
        # exceeded the two-hour cell timeout, TimeoutExpired propagated out of
        # main(), and every remaining cell was lost along with the run record
        # for the hundred and eight that had already succeeded.  A cell that
        # cannot finish is a fact about that cell, recorded and stepped over.
        limit = cell.get('timeout', 7200)
        try:
            result = subprocess.run(cmd, capture_output=True, text=True,
                                    timeout=limit, env=env)
            ok = result.returncode == 0
            note = '' if ok else (result.stderr or '').strip().splitlines()[-1:] and \
                   (result.stderr or '').strip().splitlines()[-1][:60] or ''
        except subprocess.TimeoutExpired:
            ok, note = False, 'exceeded %ds' % limit
            failures.append({'cell': artifact_name(cell),
                             'reason': 'timeout after %ds' % limit})
        except Exception as exc:                               # noqa: BLE001
            ok, note = False, type(exc).__name__
            failures.append({'cell': artifact_name(cell),
                             'reason': '%s: %s' % (type(exc).__name__, exc)})
        else:
            if not ok:
                failures.append({'cell': artifact_name(cell),
                                 'reason': note or 'exit %d' % result.returncode})

        print('%-46s %-7s %7.1fs %s'
              % (artifact_name(cell), 'ok' if ok else 'FAILED',
                 time.time() - started, note), flush=True)
        if ok:
            written += 1
            measured_here.append(artifact_name(cell))
            # Judge the cell by its own spread and record the verdict beside
            # the number.  The canary brackets the session; it cannot speak
            # for an individual cell, least of all one pinned to a single
            # thread while the canary itself ran on twelve.
            if cell.get('timed'):
                try:
                    import machine as _m
                    with open(path) as handle:
                        art = json.load(handle)
                    u = art.get('us_per_point') or {}
                    good, why = _m.admissible_stats(u.get('mean'), u.get('sd'),
                                                    u.get('min'), u.get('n'))
                    art['admissible'] = {
                        'ok': bool(good), 'why': why,
                        'block_cv': (u['sd']/u['mean']
                                     if u.get('mean') else None)}
                    with open(path, 'w') as handle:
                        json.dump(art, handle, indent=2, sort_keys=True)
                        handle.write('\n')
                    if not good:
                        print('%-46s %-7s %s'
                              % (artifact_name(cell), 'SPREAD', why),
                              flush=True)
                except (ValueError, OSError, KeyError):
                    pass

    # The run record.  The per-cell artifacts are written by the binaries,
    # which know nothing about the canary, so the verdict on whether this
    # session may be believed has to live beside them rather than inside
    # them.  Without this the drift was printed to a terminal and lost, and
    # an artifact directory carried no evidence of the machine that made it.
    record = {'schema': 'bench-run/1',
              'generated_by': 'tests/bench/run_all.py',
              'generated_at': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
              'manifest_sha256': manifest_sha(),
              'cells': len(matrix), 'artifacts_written': written,
              'artifacts_reused': skipped, 'failures': failures,
              'timed_cells': len(timed)}
    if timed:
        import machine
        after = machine.canary()
        believable, why = machine.canary_drift(before, before, after)
        print('canary after:  %.4f us/probability -> %s (%s)'
              % (after, 'believable' if believable else 'REJECTED', why))
        record['environment'] = machine.environment()
        record['canary'] = {'before': before, 'after': after,
                            'baseline_us_per_probability': 0.0508}
        record['believable'] = bool(believable)
        record['why'] = why
        if not believable:
            print('  the timed artifacts in this directory record a machine '
                  'that moved under them; they are kept, and marked.')

    # Records accumulate rather than overwrite.  A paused-and-resumed matrix
    # is several runs, each with its own canary bracket and its own machine,
    # and keeping only the last one would throw away the evidence for every
    # cell measured before the pause.  run_record.json is the latest; the
    # stamped copies are the history.
    record['cells_measured_this_run'] = measured_here
    stamp = time.strftime('%Y%m%dT%H%M%S')
    for name in ('run_record.json', 'run_record_%s.json' % stamp):
        with open(os.path.join(args.outdir, name), 'w') as handle:
            json.dump(record, handle, indent=2, sort_keys=True)
            handle.write('\n')
    if skipped:
        print('NOTE: %d cells were reused from earlier runs.  Timed cells from '
              'different runs were measured on different machine states; each '
              'run_record_*.json carries its own canary bracket, and a '
              'cross-code speed claim should not mix them without checking.'
              % skipped)

    print('%d artifacts in %s (%d reused), %d failed, manifest %s'
          % (written, os.path.relpath(args.outdir, ROOT), skipped,
             len(failures), manifest_sha()[:12]))
    for f in failures:
        print('  FAILED %-46s %s' % (f['cell'], f['reason']))


if __name__ == '__main__':
    main()
