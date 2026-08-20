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

#: The knob a timed cell holds fixed.  Speed is reported at each code's best
#: available precision, so that no code is timed in a configuration its authors
#: would not use, and the accuracy axis says what that precision buys.
SPEED_KNOB = {
    'NuFast-LBL': -1, 'NuFast-Earth': -1, 'Prob3++': 256, 'GLoBES': 256,
    'nuSQuIDS': 7, 'nuCraft': 5, 'NuOscProbExact': 256,
}


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
        for grid in ('CHORD/12x1', 'OSC/100x100'):
            for code in EARTH_CODES:
                out.append({'kind': 'amortized', 'code': code, 'grid': grid,
                            'knob': SPEED_KNOB[code], 'timed': True,
                            'shell_density': 'midpoint'})
        for code in CONST_CODES:
            out.append({'kind': 'amortized', 'code': code, 'grid': 'CONST/60E',
                        'knob': SPEED_KNOB[code], 'timed': True,
                        'shell_density': 'midpoint'})

        # THROUGHPUT: what one request for N points costs, with every
        # repetition started afresh.  The N-sweep is the figure's x axis.
        for n in (1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000):
            for code in CONST_CODES:
                out.append({'kind': 'throughput', 'code': code,
                            'grid': 'CONST/N-sweep', 'knob': SPEED_KNOB[code],
                            'n_energies': n, 'timed': True,
                            'shell_density': 'midpoint'})
        # Objection LBL-1's control: the same code, the same entry point, one
        # energy per call, so that what batching buys is measured rather than
        # asserted.
        for n in (1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000):
            out.append({'kind': 'throughput', 'code': 'NuFast-LBL',
                        'grid': 'CONST/N-sweep', 'knob': SPEED_KNOB['NuFast-LBL'],
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
                        'grid': 'CHORD/12x1', 'knob': SPEED_KNOB['NuFast-Earth'],
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
    if cell.get('n_energies'):
        cmd += ['--n-energies', str(cell['n_energies'])]
    if cell.get('n_layers') and code in BINARY:
        cmd += ['--n-layers', str(cell['n_layers'])]
    if cell.get('loop') and code in BINARY:
        cmd += ['--loop']
    if cell.get('shell_density') == 'mean' and code == 'NuFast-Earth':
        cmd += ['--mean-density']
    return cmd


def manifest_sha():
    with open(os.path.join(HERE, 'manifest.json'), 'rb') as handle:
        return hashlib.sha256(handle.read()).hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--accuracy-only', action='store_true')
    ap.add_argument('--speed-only', action='store_true')
    ap.add_argument('--outdir', default=os.path.join(HERE, 'artifacts'))
    args = ap.parse_args(argv)

    matrix = cells(args.accuracy_only, args.speed_only)
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
        before = machine.canary()
        print('canary before: %.4f us/probability' % before)

    os.makedirs(args.outdir, exist_ok=True)
    written = 0
    for cell in matrix:
        path = os.path.join(args.outdir, artifact_name(cell))
        cmd = command(cell) + ['--json', path]
        started = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True,
                                timeout=7200)
        ok = result.returncode == 0
        print('%-44s %s  %.1fs' % (artifact_name(cell),
                                   'ok' if ok else 'FAILED', time.time() - started))
        if ok:
            written += 1

    if timed:
        import machine
        after = machine.canary()
        believable, why = machine.canary_drift(before, before, after)
        print('canary after:  %.4f us/probability -> %s (%s)'
              % (after, 'believable' if believable else 'REJECTED', why))

    print('%d artifacts in %s, manifest %s'
          % (written, os.path.relpath(args.outdir, ROOT), manifest_sha()[:12]))


if __name__ == '__main__':
    main()
