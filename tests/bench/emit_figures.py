# -*- coding: utf-8 -*-
r"""Turns the measured artifacts into the files the paper figures draw from.

The pipeline measures; the notebook plots; this is the join between them, and
it exists because that join was previously a set of hand-maintained files
with no generator.  Two of them --- ``tests/speed_accuracy.json`` and
``tests/timing_other_codes.json`` --- carried timings of other people's codes
with no recorded run behind them, and both embodied the errors NuFast's
author wrote in about: no exact-eigenvalue point anywhere, GLoBES and Prob3++
recorded as exposing no precision dial when each exposes nine settings, and a
NuFast-LBL N_Newton=2 error of 8.3e-12 that was our unit conversion rather
than his algorithm.  Regenerating them from artifacts is what stops that
happening again: every number now has a run, a manifest hash and a date.

The shapes here are the ones the notebook already plots.  That is deliberate:
the plotting code is tuned --- axis limits, annotation offsets, a colour per
code held across three figures --- and rewriting it to a new shape would risk
the figures for no gain.  So the data moves and the drawing does not.

Written by ``run_all.py --join``; also runnable alone for a dry rehearsal::

    python tests/bench/emit_figures.py
"""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)

STAMP = ('Regenerated from measured artifacts by tests/bench/emit_figures.py. '
         'Every point pairs a timed measurement with the accuracy the SAME '
         'setting reached against that code\'s own 50-digit reference, in '
         'that code\'s own conventions, on the continuous PREM.')


def _tol(knob):
    r"""A tolerance dial's label, in the form the annotations use."""
    return '%.0e' % (10.0**-abs(knob))


#: How each code's knob reads on the figure, per grid.  A count is printed as
#: a count and a tolerance as a tolerance, because the annotation lookups key
#: on these strings and a reader has to recognise them.
def label_for(code, knob):
    if code == 'nuSQuIDS' or code == 'nuCraft':
        return _tol(knob)
    if code == 'NuOscProbExact' and knob < 0:
        return _tol(knob)
    return str(knob)


def load_accuracy():
    out = {}
    for grid, name in (('CONST/60E', 'accuracy_const.json'),
                       ('CHORD/12x1', 'accuracy_chord.json'),
                       ('CHORD/12x1|layers', 'accuracy_discretisation.json')):
        path = os.path.join(HERE, name)
        if os.path.exists(path):
            out[grid] = json.load(open(path))
    return out


def timed(outdir, name):
    path = os.path.join(outdir, name)
    if not os.path.exists(path):
        return None
    try:
        d = json.load(open(path))
    except (ValueError, OSError):
        return None
    return d if (d.get('us_per_point') or {}).get('mean') else None


#: Settings that are drawn whether or not anything ties them on speed.
#:
#: The "fastest at this accuracy" rule is right for an inert dial and wrong
#: for a sentinel.  NuFast-LBL reaches 7.607e-16 at both N_Newton = 3 and at
#: the exact-eigenvalue mode, and 3 is the quicker of the two, so the rule
#: dropped -1 from the figure -- the one setting objection LBL-3 is entirely
#: about, removed from the figure meant to answer it.  A named mode is a
#: claim about what the code can do, not a candidate for being outrun.
ALWAYS_DRAWN = {'NuFast-LBL': (-1,), 'NuFast-Earth': (-1,)}


def mark_best(points, code=None):
    r"""Flags, per distinct accuracy, the fastest point that reached it.

    On constant density four of the codes have an inert dial: every setting
    lands on the same accuracy because there is no profile to subdivide, so
    the sweep is a stack of points at one height.  Drawing all of them says
    nothing and drawing the slowest would misreport the code; the fastest at
    a given accuracy is the one a user would have.
    """
    best = {}
    for q in points:
        key = None if q['max_abs_error'] is None else '%.6e' % q['max_abs_error']
        if key not in best or q['us_per_probability'] < best[key]['us_per_probability']:
            best[key] = q
    keep = ALWAYS_DRAWN.get(code, ())
    for q in points:
        q['best_at_this_accuracy'] = (
            q is best.get(None if q['max_abs_error'] is None
                          else '%.6e' % q['max_abs_error'])
            or q['knob'] in keep)
    return points


def const_plane(cells_fn, artifact_name, outdir, acc):
    r"""Fig. 10: speed against accuracy at constant density."""
    series = {}
    for cell in cells_fn():
        if not (cell.get('plane') and cell['grid'] == 'CONST/60E'):
            continue
        d = timed(outdir, artifact_name(cell))
        if d is None:
            continue
        code, knob = cell['code'], cell['knob']
        name = code + (' (1 thread)' if cell.get('threads') == 1 else '')
        entry = (acc.get('CONST/60E', {}).get('series', {}).get(code, {})
                 .get('by_knob', {}).get(str(knob), {}))
        series.setdefault(name, {
            'name': name, 'code': code,
            'dial': 'N_Newton' if code.startswith('NuFast') else
                    ('solver tolerance' if code == 'nuSQuIDS' else
                     'n_shells' if code in ('Prob3++', 'GLoBES') else
                     'n_slabs_per_segment | rtol'),
            'threads': cell.get('threads'), 'points': []})['points'].append({
                'label': label_for(code, knob), 'knob': knob,
                'us_per_probability': d['us_per_point']['mean'],
                'us_sd': d['us_per_point']['sd'],
                'block_cv': (d['us_per_point']['sd']/d['us_per_point']['mean']
                             if d['us_per_point']['mean'] else None),
                'max_abs_error': entry.get('max_abs_deviation')})
    for s in series.values():
        s['points'].sort(key=lambda q: q['knob'])
        mark_best(s['points'], s.get('code'))
    return {'generated_by': 'tests/bench/emit_figures.py',
            'note': 'Constant density, L = 1300 km, rho = 3 g/cm^3, three '
                    'flavors, nu_mu row. ' + STAMP,
            'series': [series[k] for k in sorted(series)]}


#: Fig. 11 draws each code against how finely it cuts the PREM profile.  For
#: three of them the shell or slab count IS the precision knob; NuFast-Earth
#: dials eigenvalue precision instead and takes its layer count separately,
#: so its curve comes from the discretisation cells and the layer sweep.
EARTH_FROM_LAYERS = {'NuFast-Earth': -1}


def earth_plane(cells_fn, artifact_name, outdir, acc):
    r"""Fig. 11, three flavors: speed against accuracy through the Earth."""
    series = {}

    def add(name, code, dial, label, knob, d, err):
        series.setdefault(name, {'name': name, 'code': code, 'dial': dial,
                                 'points': []})['points'].append({
            'label': label, 'knob': knob,
            'us_per_probability': d['us_per_point']['mean'],
            'us_sd': d['us_per_point']['sd'],
            'max_abs_error': err})

    for cell in cells_fn():
        if cell['grid'] != 'CHORD/12x1' or not cell.get('timed'):
            continue
        if cell.get('threads') == 1:
            continue
        code = cell['code']
        d = timed(outdir, artifact_name(cell))
        if d is None:
            continue

        if cell.get('sweep') == 'discretisation' and code in EARTH_FROM_LAYERS:
            layers = cell.get('n_layers')
            # Accuracy is shared between the two scans: the untimed protocol
            # evaluates at the central value of whichever parameter moves, so
            # both lines are the same physics measured at the same settings.
            err = (acc.get('CHORD/12x1|layers', {}).get('series', {})
                   .get(code, {}).get('by_knob', {})
                   .get(str(layers), {}).get('max_abs_deviation'))
            # Two lines, because one number cannot describe this code here.
            # Scanning delta_CP invalidates none of its layer work, so that
            # curve is flat at 0.057 us whatever the layer count; scanning
            # Dmsq31 costs it 107.  The unqualified name goes to the second,
            # since that is the cost when a fit moves a parameter the cache
            # does not cover.
            name = (code if cell.get('scan') == 'dmsq31'
                    else code + ' (dCP only)')
            add(name, code, 'n_layers', str(layers), layers, d, err)
            continue
        if not cell.get('plane') or code in EARTH_FROM_LAYERS:
            continue

        knob = cell['knob']
        # An explicit tolerance is its own dial, keyed separately: the knob
        # is an integer and cannot carry 3, 0.3 or 3e-5, which is why six of
        # the paper's nine tolerance points had gone missing.
        if cell.get('rtol'):
            err = (acc.get('CHORD/12x1', {}).get('series', {}).get(code, {})
                   .get('by_rtol', {}).get('%g' % cell['rtol'], {})
                   .get('max_abs_deviation'))
            add('NuOscProbExact (tolerance)', code, 'rtol',
                '%.0e' % cell['rtol'], -cell['rtol'], d, err)
            continue
        err = (acc.get('CHORD/12x1', {}).get('series', {}).get(code, {})
               .get('by_knob', {}).get(str(knob), {}).get('max_abs_deviation'))
        # This library turns two different dials, and the figure has always
        # drawn them as two curves: the slab count is the code told what to
        # do, the tolerance is the same code asked to find out.
        if code == 'NuOscProbExact' and knob < 0:
            add('NuOscProbExact (tolerance)', code, 'rtol',
                label_for(code, knob), knob, d, err)
        else:
            add(code, code,
                'solver tolerance' if code == 'nuSQuIDS' else
                'numPrec' if code == 'nuCraft' else
                'eigenvalue precision' if code.startswith('NuFast') else
                'n_shells' if code in ('Prob3++', 'GLoBES') else
                'n_slabs_per_segment',
                label_for(code, knob), knob, d, err)

    for s in series.values():
        # By cost.  Sorting by error made a non-monotone code -- Prob3++ is
        # worse at 512 shells than at 256, which is its radial shell stack
        # and not noise -- draw as a line doubling back on itself.  Ordered
        # by cost the curve runs left to right and the rise is visible as a
        # rise, which is what it is.
        s['points'].sort(key=lambda q: q['us_per_probability'])
        mark_best(s['points'], s.get('code'))
    return {'generated_by': 'tests/bench/emit_figures.py',
            'note': 'PREM, cos(theta_z) = -0.9, E = 3-40 GeV, three flavors, '
                    'nu_mu row. ' + STAMP,
            'series': [series[k] for k in sorted(series)]}


def throughput(cells_fn, artifact_name, outdir):
    r"""Fig. 8: cost per probability against how many are asked for at once.

    The old file behind this figure recorded NuFast-LBL unbatched, and the
    figure's own text said that neither it, GLoBES nor Prob3++ batches at
    all.  That was true of GLoBES and Prob3++ and false of NuFast-LBL, which
    has taken a vector of energies since v2.0.0 --- objection LBL-1.  Here it
    is driven through that entry point, and the looped control is measured
    beside it rather than assumed, so what batching is worth is a number in
    the figure instead of a claim about it.
    """
    series = {}
    for cell in cells_fn():
        if cell.get('kind') != 'throughput':
            continue
        d = timed(outdir, artifact_name(cell))
        if d is None:
            continue
        name = cell['code'] + (' (looped)' if cell.get('loop') else '')
        n = cell.get('n_energies') or d['n_points']
        entry = series.setdefault(name, {
            'code': cell['code'], 'knob': cell['knob'],
            # From the code's OWN reported capability, not from how this
            # cell was driven.  Deriving it from the --loop flag marked
            # GLoBES and Prob3++ as batching, which they do not: the
            # harness loops over energies for them because neither exposes
            # an entry point that takes more than one.
            'batched': bool((d.get('batched') or {}).get('energy'))
                       and not cell.get('loop'),
            'batch_symbol': (d.get('batched') or {}).get('symbol', ''),
            'sizes': [], 'seconds': [], 'us_per_probability': []})
        entry['sizes'].append(n)
        entry['seconds'].append(d['us_per_point']['mean']*1e-6*n)
        entry['us_per_probability'].append(d['us_per_point']['mean'])
    for s in series.values():
        order = sorted(range(len(s['sizes'])), key=lambda i: s['sizes'][i])
        for k in ('sizes', 'seconds', 'us_per_probability'):
            s[k] = [s[k][i] for i in order]
    return {'generated_by': 'tests/bench/emit_figures.py',
            'note': 'One request for N energies at constant density, every '
                    'repetition cold. ' + STAMP,
            'series': series}


def scan_sensitivity(cells_fn, artifact_name, outdir):
    r"""How much each code's cost depends on WHICH parameter the fit moves.

    Objection Earth-2 named Dmsq31 and said different parameters cost
    different amounts of recomputation.  Each Dmsq31 cell here is paired with
    the delta_CP cell at the same code, grid and knob -- the two differ in
    nothing else -- so the ratio is the answer to that, per code, measured.
    """
    out = {}
    for cell in cells_fn():
        if cell.get('scan') != 'dmsq31':
            continue
        slow = timed(outdir, artifact_name(cell))
        control = dict(cell)
        control.pop('scan', None)
        fast = timed(outdir, artifact_name(control))
        if not slow or not fast:
            continue
        grid = cell['grid']
        out.setdefault(grid, []).append({
            'code': cell['code'], 'knob': cell['knob'],
            'dcp_us': fast['us_per_point']['mean'],
            'dcp_sd': fast['us_per_point']['sd'],
            'dmsq31_us': slow['us_per_point']['mean'],
            'dmsq31_sd': slow['us_per_point']['sd'],
            'ratio': slow['us_per_point']['mean']/fast['us_per_point']['mean']})
    for grid in out:
        out[grid].sort(key=lambda q: -q['ratio'])
    return {'generated_by': 'tests/bench/emit_figures.py',
            'note': 'Cost per probability when a fit moves delta_CP against '
                    'when it moves Dmsq31, at the same knob, differing in '
                    'nothing else. ' + STAMP,
            'grids': out}


def input_manifest_hashes(outdir):
    r"""Every distinct ``manifest_sha256`` carried by the artifacts in `outdir`.

    Returned sorted, so the emitted files are byte-stable across runs.
    """
    seen = set()
    for name in os.listdir(outdir):
        if not name.endswith('.json'):
            continue
        try:
            with open(os.path.join(outdir, name)) as handle:
                sha = json.load(handle).get('manifest_sha256')
        except (ValueError, OSError):
            continue
        if sha:
            seen.add(sha)
    return sorted(seen)


def main(outdir=None):
    import run_all
    outdir = outdir or os.path.join(HERE, 'artifacts')
    acc = load_accuracy()
    written = []
    for path, payload in (
            (os.path.join(ROOT, 'tests', 'speed_accuracy.json'),
             const_plane(run_all.cells, run_all.artifact_name, outdir, acc)),
            (os.path.join(HERE, 'earth_plane.json'),
             earth_plane(run_all.cells, run_all.artifact_name, outdir, acc)),
            (os.path.join(ROOT, 'tests', 'timing_other_codes.json'),
             throughput(run_all.cells, run_all.artifact_name, outdir)),
            (os.path.join(HERE, 'scan_sensitivity.json'),
             scan_sensitivity(run_all.cells, run_all.artifact_name, outdir))):
        payload['manifest_sha256'] = run_all.manifest_sha()
        # That stamp is the manifest as it stands now, which is not the same
        # thing as the manifest each input was measured under: a sweep resumed
        # after the manifest changed leaves artifacts carrying an older hash.
        # Recording the hashes actually present stops a file from claiming a
        # provenance its inputs do not share.
        seen = input_manifest_hashes(outdir)
        payload['manifest_sha256_inputs'] = seen
        if len(seen) > 1:
            payload['manifest_note'] = (
                'Inputs were measured under more than one manifest; the hashes '
                'are listed in manifest_sha256_inputs.')
        with open(path, 'w') as handle:
            json.dump(payload, handle, indent=1, sort_keys=True)
            handle.write('\n')
        body = payload.get('series', payload.get('grids', []))
        n = len(body) if isinstance(body, (list, dict)) else 0
        written.append((os.path.relpath(path, ROOT), n))
    for rel, n in written:
        print('  %-44s %2d series' % (rel, n))
    return written


if __name__ == '__main__':
    main()
