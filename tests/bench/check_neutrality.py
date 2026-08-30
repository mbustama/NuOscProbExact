# -*- coding: utf-8 -*-
r"""Measures whether the timing rebuild moved the numbers it should not have.

The block-length fix was meant to cure a spread, not to change a speed.  For
the cells that were ALREADY clean under the old harness --- a long block, a
per-cent spread, nothing for the fix to repair --- the new harness must
reproduce the old number.  If it does not, the change is not the neutral
repair it was argued to be, and the full re-run should wait until that is
understood.

This has to run BEFORE ``run_all.py --speed-only --force``, which overwrites
the old artifacts it compares against.  It also has to run on an idle
machine: attempted under a load of four it reported discrepancies of up to
1.45x that were the other session's test suite, not the harness.

Run from the repository root::

    python tests/bench/check_neutrality.py

A ratio near one for every cell means the rebuild changed the spread and left
the speed alone, which is the whole claim.
"""

import glob
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

#: A cell qualifies as a witness only if its OLD reading was trustworthy:
#: a spread this tight over a block this long leaves nothing for the fix to
#: repair, so any movement is the harness and not the repair.
MAX_OLD_CV = 0.05
MIN_OLD_BLOCK_S = 0.05
OLD_STEPS = 25

BIN = {'NuFast-LBL': 'bench_nufast_lbl', 'NuFast-Earth': 'bench_nufast_earth',
       'Prob3++': 'bench_prob3', 'GLoBES': 'bench_globes'}
ADAPTER = {'nuSQuIDS': 'nusquids', 'nuCraft': 'nucraft',
           'NuOscProbExact': 'nuoscprobexact'}


def witnesses(limit):
    r"""Returns the cleanest old cells, cheapest first so the check is quick."""
    out = []
    for path in glob.glob(os.path.join(HERE, 'artifacts', 'amortized_*.json')):
        try:
            d = json.load(open(path))
        except (ValueError, OSError):
            continue
        u = d.get('us_per_point') or {}
        if not u.get('mean') or '1thread' in path or 'looped' in path:
            continue
        block = u['mean']*1e-6*d['n_points']*OLD_STEPS
        if u['sd']/u['mean'] > MAX_OLD_CV or block < MIN_OLD_BLOCK_S:
            continue
        knob = list(d['knob'].values())[0]
        out.append((u['mean']*d['n_points'], d['code'], d['protocol']['grid'],
                    knob, u['mean'], u['sd']/u['mean']))
    out.sort()
    return out[:limit]


def measure(code, grid, knob):
    if code in BIN:
        cmd = [os.path.join(ROOT, '.bench-build', 'bin', BIN[code])]
    else:
        cmd = [sys.executable, os.path.join(HERE, 'runner.py'), ADAPTER[code]]
    cmd += ['--protocol', 'amortized', '--grid', grid, '--knob', str(knob)]
    got = subprocess.run(cmd, capture_output=True, text=True, timeout=5400)
    if got.returncode != 0 or '"us_per_point"' not in got.stdout:
        return None
    return json.loads(got.stdout[got.stdout.index('{'):])


def main():
    load = os.getloadavg()[0]
    print('load average %.2f -- this check is meaningless above about 0.5\n'
          % load)
    rows = witnesses(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
    print('%-16s %-12s %-6s %11s %11s %8s %7s'
          % ('CODE', 'GRID', 'KNOB', 'OLD us/pt', 'NEW us/pt', 'RATIO', 'CV NEW'))
    worst = 0.0
    for _, code, grid, knob, old, oldcv in rows:
        got = measure(code, grid, knob)
        if not got:
            print('%-16s %-12s %-6s      FAILED' % (code, grid, knob))
            continue
        u = got['us_per_point']
        ratio = u['mean']/old
        worst = max(worst, abs(ratio - 1.0))
        print('%-16s %-12s %-6s %11.4f %11.4f %7.3fx %6.1f%%'
              % (code, grid, knob, old, u['mean'], ratio,
                 100*u['sd']/u['mean']))
    print()
    if worst <= 0.10:
        print('NEUTRAL: every witness within %.0f%% of its old value.  The '
              'rebuild changed the spread and left the speed alone; the full '
              're-run can proceed.' % (worst*100))
    else:
        print('NOT NEUTRAL: worst witness moved %.0f%%.  Either the machine '
              'is not idle (check the load above) or the rebuild changed a '
              'number it should not have.  Do not start the full re-run until '
              'this is understood.' % (worst*100))


if __name__ == '__main__':
    main()
