# -*- coding: utf-8 -*-
r"""Measures every code's accuracy against its own reference, over every knob.

The accuracy axis, end to end.  Nothing here reads a clock --- the protocol it
drives is untimed by construction --- so this can run on a machine that is
busy, and its numbers are unaffected by what else is happening.  That matters:
the speed axis needs an idle machine and this one does not, so the two need
not compete for the same window.

The expensive part is the references, and they are cached for a reason worth
stating.  A reference depends on the code, the chord and the energy; it does
not depend on the precision knob, because the knob changes what the *code*
does and not what the right answer is.  So a nine-point knob sweep costs one
reference, not nine.  Without that the sweep would be an order of magnitude
more work and would tempt whoever ran it into shortening the grid.

The cache is on disk and keyed by everything that can change an answer, so a
run that dies half way resumes instead of restarting, and a run against
changed conventions does not silently reuse stale numbers.

Usage::

    python tests/bench/sweep_accuracy.py [--grid CHORD/12x1] [--n-energies 12]
        [--codes NuFast-LBL,Prob3++] [--out artifact.json]
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
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, 'adapters'))
sys.path.insert(0, os.path.join(ROOT, 'src'))

import conversions                                             # noqa: E402
import reference as ref                                        # noqa: E402

CACHE = os.path.join(HERE, 'reference_cache.json')

#: Which binary or adapter answers for each code, and the knob domain swept.
COMPILED = {
    'NuFast-LBL':   'bench_nufast_lbl_accuracy',
    'NuFast-Earth': 'bench_nufast_earth_accuracy',
    'Prob3++':      'bench_prob3_accuracy',
    'GLoBES':       'bench_globes',
}
PYTHON = {
    'Second-order expansion': 'second_order',
    'nuSQuIDS':       'nusquids',
    'nuCraft':        'nucraft',
    'NuOscProbExact': 'nuoscprobexact',
}
KNOBS = {
    # No dial: the series is truncated at second order and that is that.
    'Second-order expansion': [0],
    'NuFast-LBL':     [-1, 0, 1, 2, 3],
    'NuFast-Earth':   [-1, 0, 1, 2, 3],
    # To 1024: the curves were stopping at 256 because this list did, not
    # because any code did.  All three accept more.
    'Prob3++':        [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096,
                      # Not a floor, though it was reported as one: 8192
                      # and 16384 read 8.0e-08 and 7.8e-08, but 12288
                      # between them is 4.3e-08 and 65536 reaches 2.4e-09.
                      # The error oscillates with the shell holding the
                      # turning point, so these rungs are deliberately not
                      # powers of four -- a sparse ladder misreads it.
                      8192, 12288, 24576, 32768, 49152, 65536],
    'GLoBES':         [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096,
                      # Not a floor, though it was reported as one: 8192
                      # and 16384 read 8.0e-08 and 7.8e-08, but 12288
                      # between them is 4.3e-08 and 65536 reaches 2.4e-09.
                      # The error oscillates with the shell holding the
                      # turning point, so these rungs are deliberately not
                      # powers of four -- a sparse ladder misreads it.
                      8192, 12288, 24576, 32768, 49152, 65536],
    'nuSQuIDS':       [3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
    'nuCraft':        [2, 3, 4, 5, 6, 7, 8, 9, 10],
    'NuOscProbExact': [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048,
                       4096, 8192, 12288, 24576, 32768, 49152, 65536,
                       # Reaching past -5 needs `n_max` raised above the
                       # shipped 1024; the adapter passes it.
                       -3, -4, -5, -6, -7, -8],
}

#: Codes that cannot be posed one of the two problem kinds, and why.
#:
#: Getting this wrong is silent rather than loud, which is why it is a table
#: rather than a comment.  Asked for a chord, NuFast-LBL's adapter ignores the
#: zenith angle and propagates the constant-density baseline instead -- it
#: returns numbers, they are simply answers to a different question, and
#: differencing them against a chord reference produced a knob-independent
#: 0.72 that looked like a catastrophically inaccurate code rather than a
#: mis-posed problem.
CONST_REFUSED = {'nuCraft': 'propagates through its Earth model only'}
CHORD_REFUSED = {'NuFast-LBL': 'constant-density baselines only; it has no '
                               'Earth model and no zenith axis'}


def logspace(lo, hi, n):
    return [lo if n == 1 else lo*(hi/lo)**(float(i)/(n - 1)) for i in range(n)]


def grid_points(grid, n_e, n_z):
    r"""Returns ``(energies_gev, costhz_list)`` for `grid`, as the harness does."""
    manifest = json.load(open(os.path.join(HERE, 'manifest.json')))
    entry = manifest['grids'][grid]
    e = entry['energies_gev']
    lo, hi = e['log']
    n = n_e or (e['n'] if isinstance(e['n'], int) else 12)
    energies = logspace(lo, hi, n)
    z = entry.get('costhz')
    if isinstance(z, list):
        costhz = list(z)
    elif isinstance(z, dict):
        a, b = z['lin']
        m = n_z or z['n']
        costhz = [a if m == 1 else a + (b - a)*float(i)/(m - 1) for i in range(m)]
    else:
        costhz = []
    return energies, costhz


#: Bumped when the SHAPE of a cached reference changes, not just its value.
#: v2 stores the nu_mu row (three channels) where v1 stored one; a v1 entry
#: reused under v2 would silently supply a scalar where a row is expected.
CACHE_SCHEMA = 'v2-numu-row'


def _conventions_key(code):
    r"""Everything about a code's conventions that can move its reference."""
    return '%s|%s|%.17g|%.17g|%s' % (
        CACHE_SCHEMA,
        code, conversions.matter_constant(code),
        conversions.km_to_inv_ev(code),
        json.dumps(conversions.oscillation_parameters(), sort_keys=True))


def load_cache():
    if os.path.exists(CACHE):
        with open(CACHE) as handle:
            return json.load(handle)
    return {}


def reference_for(code, energy_gev, costhz, ye, cache, l_km=None,
                  density=None):
    r"""Returns the reference probability, computing it only if not cached.

    The key carries the conventions, so a reference built before a constant
    changed is never silently reused after it.
    """
    key = hashlib.sha256(('%s|%.17g|%s|%.17g|%s|%s'
                          % (_conventions_key(code), energy_gev,
                             'const' if costhz is None else '%.17g' % costhz,
                             ye, l_km, density)).encode()).hexdigest()[:32]
    if key in cache:
        return cache[key]['value'], cache[key]['error'], True

    params = conversions.for_code(code)
    if costhz is None:
        value = ref.constant_density(code, [energy_gev], l_km, density, ye,
                                     params)[0]
        error = float(ref.mp.mpf(10)**(-ref.DPS + 5))
    else:
        value, err = ref.earth_chord_reference(code, energy_gev, costhz, ye,
                                               params)
        error = float(err)
    row = [str(v) for v in value]
    cache[key] = {'value': row, 'error': error, 'code': code,
                  'energy_gev': energy_gev, 'costhz': costhz}
    with open(CACHE, 'w') as handle:
        json.dump(cache, handle, indent=1, sort_keys=True)
    return row, error, False


#: The tolerance dial's loose end, which the integer knob cannot reach.
#: rtol = 10**knob expresses 1e-3 and 1e-4 and nothing between them, so six
#: of the paper's nine points -- 3, 1, 0.3, 0.03, 0.01, 3e-5 -- had no knob.
# 3e-6 is NOT reachable.  It succeeds at the central delta_CP, which is
# all an accuracy probe evaluates, and fails inside the amortized scan:
# somewhere in the 0.2 rad sweep the tolerance needs more than the
# 1024-slab ceiling and the call raises.  It would buy nothing anyway --
# measured, it reaches 5.468e-08, the same as 1e-5, because both are
# already pinned at that ceiling.  The dial stops at 1e-5.
RTOL_SWEEP = [3.0, 1.0, 0.3, 0.03, 0.01, 3.0e-5]


def probabilities(code, grid, knob, n_e, n_z, n_layers=0, rtol=0.0):
    r"""Runs the untimed accuracy protocol and returns the scored channel."""
    if code in COMPILED:
        cmd = [os.path.join(ROOT, '.bench-build', 'bin', COMPILED[code]),
               '--protocol', 'accuracy', '--grid', grid, '--knob', str(knob)]
    else:
        cmd = [sys.executable, os.path.join(HERE, 'runner.py'), PYTHON[code],
               '--protocol', 'accuracy', '--grid', grid, '--knob', str(knob)]
    if n_layers:
        cmd += ['--n-layers', str(n_layers)]
    if rtol:
        cmd += ['--rtol', repr(rtol)]
    if n_e:
        cmd += ['--n-energies', str(n_e)]
    if n_z:
        cmd += ['--n-zenith', str(n_z)]
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if out.returncode != 0 or '"probabilities"' not in out.stdout:
        return None
    return json.loads(out.stdout)['probabilities']


#: Codes whose PROFILE dial is separate from their precision dial.
#:
#: The Earth figure plots each code against how finely it cuts the PREM
#: profile, because on a chord that is what dominates the error.  For
#: Prob3++, GLoBES and this library the shell or slab count IS the
#: precision knob, so :data:`KNOBS` already measured it.  NuFast-Earth is
#: the exception: its knob is the eigenvalue precision, and its layer
#: count is a separate argument, so without this its discretisation cells
#: are timed against an accuracy that was never measured -- five points on
#: a speed axis with nothing to plot them against.
LAYER_SWEEP = {
    'NuFast-Earth': {'knob': -1,
                     'layers': [1, 4, 16, 64, 256, 1024, 4096, 8192,
                                12288, 24576, 32768, 49152, 65536]},
}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--grid', default='CHORD/12x1')
    ap.add_argument('--n-energies', type=int, default=0)
    ap.add_argument('--n-zenith', type=int, default=0)
    ap.add_argument('--codes', default='')
    ap.add_argument('--out', default=os.path.join(HERE, 'accuracy_sweep.json'))
    ap.add_argument('--rtol-sweep', action='store_true', dest='rtol_sweep',
                    help='sweep RTOL_SWEEP instead of the integer knob')
    ap.add_argument('--layers', action='store_true',
                    help='sweep the PROFILE dial (n_layers) instead of the '
                         'precision knob, for the codes in LAYER_SWEEP')
    args = ap.parse_args(argv)

    codes = ([c.strip() for c in args.codes.split(',') if c.strip()]
             or list(COMPILED) + list(PYTHON))
    energies, costhz = grid_points(args.grid, args.n_energies, args.n_zenith)
    constant = not costhz
    ye = 0.5
    l_km, density = (1300.0, 3.0) if constant else (None, None)

    cache = load_cache()
    record = {
        'schema': 'bench-accuracy-sweep/1',
        'generated_by': 'tests/bench/sweep_accuracy.py',
        'generated_at': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
        'grid': args.grid,
        'protocol': 'accuracy (untimed)',
        'dial': ('rtol' if args.rtol_sweep
                 else 'n_layers' if args.layers else 'precision knob'),
        'profile_basis': 'continuous',
        'ye': ye,
        'energies_gev': energies,
        'costhz': costhz,
        'series': {},
    }

    for code in codes:
        refused = (CONST_REFUSED if constant else CHORD_REFUSED).get(code)
        if refused:
            record['series'][code] = {'skipped': refused}
            print('%-16s skipped: %s' % (code, refused), flush=True)
            continue

        refs, ref_errs = [], []
        for cz in (costhz or [None]):
            for e in energies:
                value, error, hit = reference_for(code, e, cz, ye, cache,
                                                  l_km, density)
                refs.append(value)
                ref_errs.append(error)
                if not hit:
                    print('  %-16s reference E=%-8.3f cz=%-6s err %.0e'
                          % (code, e, cz, error), flush=True)

        entry = {'reference_error_max': max(ref_errs),
                 'channels': ['numu->nue', 'numu->numu', 'numu->nutau'],
                 'reference': [[str(v) for v in row] for row in refs],
                 'by_knob': {}}
        if args.rtol_sweep:
            if code != 'NuOscProbExact':
                record['series'][code] = {'skipped': 'no tolerance dial'}
                continue
            dials = [(0, 0, r) for r in RTOL_SWEEP]
        elif args.layers:
            spec = LAYER_SWEEP.get(code)
            if not spec:
                record['series'][code] = {
                    'skipped': 'its profile dial IS its precision knob; the '
                               'main sweep already measured it'}
                print('%-16s skipped: profile dial is its precision knob'
                      % code, flush=True)
                continue
            dials = [(spec['knob'], n, 0.0) for n in spec['layers']]
        else:
            dials = [(k, 0, 0.0) for k in KNOBS[code]]

        for knob, n_layers, rtol in dials:
            probs = probabilities(code, args.grid, knob, args.n_energies,
                                  args.n_zenith, n_layers, rtol)
            dial = ('%g' % rtol if args.rtol_sweep
                    else str(n_layers) if args.layers else str(knob))
            # A code may answer for fewer channels than the row holds: the
            # second-order expansion is an APPEARANCE formula and gives
            # P(numu->nue) only.  Scoring it against three would mean
            # inventing two it never computed.
            n_ch = len(probs)//len(refs) if probs else 0
            if probs is None or n_ch not in (1, 3):
                entry['by_knob'][dial] = {'failed': True}
                print('%-16s dial %-5s FAILED' % (code, dial), flush=True)
                continue
            from decimal import Decimal, getcontext
            getcontext().prec = 60
            # Three channels per grid point, flattened the same way in both,
            # so the comparison is elementwise and the per-channel split is
            # recoverable by striding.  The probabilities are stored as well
            # as the deviations, because a figure about the appearance
            # channel cannot be drawn from a disappearance summary and would
            # otherwise need its own data path -- which is how the previous
            # figures drifted from the pipeline that fed them.
            flat_refs = ([r[0] for r in refs] if n_ch == 1
                         else [r for row in refs for r in row])
            dev = [abs(Decimal(repr(p)) - Decimal(r)) for p, r in
                   zip(probs, flat_refs)]
            names = (('numu->nue',) if n_ch == 1
                     else ('numu->nue', 'numu->numu', 'numu->nutau'))
            per_channel = {name: [float(d) for d in dev[c::n_ch]]
                           for c, name in enumerate(names)}
            entry['by_knob'][dial] = {
                'max_abs_deviation': float(max(dev)),
                'median_abs_deviation': float(sorted(dev)[len(dev)//2]),
                'max_by_channel': {k: max(v) for k, v in per_channel.items()},
                'probabilities': probs,
            }
            print('%-16s dial %-5s max %.3e  median %.3e'
                  % (code, dial, float(max(dev)),
                     float(sorted(dev)[len(dev)//2])), flush=True)
        record['series'][code] = entry

    with open(args.out, 'w') as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
        handle.write('\n')
    print('wrote %s' % os.path.relpath(args.out, ROOT))


if __name__ == '__main__':
    main()
