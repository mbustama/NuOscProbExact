# -*- coding: utf-8 -*-
r"""Checks that the measured output answers each objection, one at a time.

Written after the check it performs was skipped once.  The pipeline had been
verified for shape, coverage and crashes --- every plane cell joins, every
notebook executes --- and none of that noticed that the constant-density
figure dropped NuFast-LBL's exact-eigenvalue point, because N_Newton = 3 ties
it on accuracy and beats it on speed.  That is the one setting objection
LBL-3 is about, removed by a filter from the figure meant to answer it.

Shape and coverage are not the same question as whether the evidence is
there.  This asks the second question, per objection, against the files the
figures actually read, and it is meant to be re-run after every measurement
rather than reasoned about once.

Run from the repository root::

    python tests/bench/audit_objections.py
"""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)

ART = os.path.join(HERE, 'artifacts')


def _load(path):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (ValueError, OSError):
        return None


def _artifacts(predicate):
    out = []
    for name in sorted(os.listdir(ART)):
        if not name.endswith('.json') or name.startswith('run_record'):
            continue
        d = _load(os.path.join(ART, name))
        if d and predicate(name, d):
            out.append((name, d))
    return out


def _series(payload, name):
    for s in payload.get('series', []):
        if s.get('name') == name:
            return s
    return None


def check_lbl1():
    r"""LBL-1: NuFast-LBL batches over energy, and is driven that way."""
    fig8 = _load(os.path.join(ROOT, 'tests', 'timing_other_codes.json'))
    if not fig8:
        return False, 'tests/timing_other_codes.json missing'
    s = fig8.get('series', {})
    lbl, loop = s.get('NuFast-LBL'), s.get('NuFast-LBL (looped)')
    if not lbl or not lbl.get('batched'):
        return False, 'NuFast-LBL is not recorded as batched'
    if not lbl.get('batch_symbol'):
        return False, 'no batch entry point recorded'
    if not loop:
        return False, 'no looped control, so what batching buys is unmeasured'
    return True, ('batched through %s, with a looped control over %d sizes'
                  % (lbl['batch_symbol'].split('(')[0], len(loop['sizes'])))


def check_lbl2():
    r"""LBL-2: the Newton sweep, in his conventions, reaching double precision."""
    fig10 = _load(os.path.join(ROOT, 'tests', 'speed_accuracy.json'))
    s = _series(fig10 or {}, 'NuFast-LBL')
    if not s:
        return False, 'NuFast-LBL absent from the constant-density plane'
    by = {q['knob']: q['max_abs_error'] for q in s['points']}
    for k in (0, 1, 2, 3):
        if k not in by or by[k] is None:
            return False, 'N_Newton = %d has no accuracy' % k
    if not by[3] < 1e-15:
        return False, 'N_Newton = 3 does not reach double precision: %.2e' % by[3]
    if not by[0] > 1e-6:
        return False, 'N_Newton = 0 looks too good: %.2e' % by[0]
    return True, ('N_Newton 0,1,2,3 -> %.2e, %.2e, %.2e, %.2e; he stated '
                  '~1e-15 at 2 and double precision at 3'
                  % (by[0], by[1], by[2], by[3]))


def check_lbl3():
    r"""LBL-3: the exact mode is measured AND survives into the figure."""
    fig10 = _load(os.path.join(ROOT, 'tests', 'speed_accuracy.json'))
    notes = []
    for name in ('NuFast-LBL', 'NuFast-Earth'):
        s = _series(fig10 or {}, name)
        if not s:
            return False, '%s absent from the constant-density plane' % name
        exact = [q for q in s['points'] if q['knob'] == -1]
        if not exact:
            return False, '%s has no exact-eigenvalue point' % name
        if not exact[0].get('best_at_this_accuracy'):
            return False, ('%s exact point is measured but filtered out of '
                           'the figure' % name)
        notes.append('%s -1 drawn at %.2e' % (name, exact[0]['max_abs_error']))
    return True, '; '.join(notes)


def check_earth1():
    r"""Earth-1: whether 256 was per region or in total is now recorded."""
    bad = _artifacts(lambda n, d: d.get('n_layers') is not None
                     and d.get('shells_total') is None)
    if bad:
        return False, '%d artifacts record n_layers but not the total' % len(bad)
    got = _artifacts(lambda n, d: d.get('shells_total'))
    if not got:
        return False, 'no artifact records shells_total'
    ex = got[0][1]
    return True, ('recorded in %d artifacts, e.g. n_layers=%s -> %s shells '
                  'in total' % (len(got), ex['n_layers'], ex['shells_total']))


def check_earth3():
    r"""Earth-3: the realistic grid, not twelve energies at one zenith."""
    osc = _artifacts(lambda n, d: 'OSC_100x100' in n and d.get('us_per_point'))
    codes = sorted({d['code'] for _, d in osc})
    if len(codes) < 6:
        return False, 'oscillogram measured for only %d codes: %s' % (
            len(codes), ', '.join(codes))
    pts = osc[0][1]['n_points']
    return True, '%d codes on a %d-point grid' % (len(codes), pts)


def check_earth4():
    r"""Earth-4: NuFast-Earth does have a constant-density mode, and it is used."""
    got = _artifacts(lambda n, d: d['code'] == 'NuFast-Earth'
                     and 'CONST' in n)
    if not got:
        return False, 'NuFast-Earth never posed the constant-density problem'
    return True, 'measured on the constant-density grid in %d cells' % len(got)


def check_other1():
    r"""Other-1: every number ties to the manifest that pins the versions."""
    missing = _artifacts(lambda n, d: 'manifest_sha256' not in d)
    if missing:
        return False, ('%d artifacts carry no manifest hash (e.g. %s)'
                       % (len(missing), missing[0][0]))
    shas = {d['manifest_sha256'] for _, d in _artifacts(lambda n, d: True)}
    if len(shas) > 1:
        return False, ('artifacts span %d different manifests; a run must '
                       'leave exactly one' % len(shas))
    return True, 'every artifact names the one manifest %s' % list(shas)[0][:12]


def check_other2():
    r"""Other-2: enough repetitions, and a spread that means something."""
    timed = _artifacts(lambda n, d: d.get('us_per_point', {}).get('mean'))
    if not timed:
        return False, 'no timed artifacts'
    no_stamp = [n for n, d in timed if 'timing' not in d]
    if no_stamp:
        return False, ('%d timed artifacts do not record how long a block '
                       'was (e.g. %s)' % (len(no_stamp), no_stamp[0]))
    few = [n for n, d in timed if d['us_per_point']['n'] < 30]
    wide = [n for n, d in timed
            if d['us_per_point']['sd']/d['us_per_point']['mean'] > 0.10]
    reps = min(d['timing']['block_reps'] for _, d in timed)
    if few:
        return False, '%d cells have fewer than 30 blocks (e.g. %s)' % (
            len(few), few[0])
    if wide:
        return False, ('%d of %d cells still spread beyond 10 per cent '
                       '(e.g. %s)' % (len(wide), len(timed), wide[0]))
    return True, ('%d cells, every one at 30+ blocks, none above 10 per cent, '
                  'fewest repetitions in any block %d' % (len(timed), reps))


def check_earth2():
    r"""Earth-2: the loop is realistic, and BOTH named parameters are scanned.

    The objection asked for a loop a fit would actually run and named Dmsq31.
    Scanning only delta_CP answers half of it: for a code that caches by
    parameter, delta_CP is the cheapest thing a fit can move, so reporting
    that alone as "the" amortized cost flatters exactly the code the
    objection came from.  Both are measured now, at the same knob, differing
    in nothing else.
    """
    pairs, notes = 0, []
    for name, d in _artifacts(lambda n, d: n.endswith('_scan-dmsq31.json')
                              and d.get('us_per_point', {}).get('mean')):
        base = os.path.join(ART, name.replace('_scan-dmsq31', ''))
        b = _load(base)
        if not b or not (b.get('us_per_point') or {}).get('mean'):
            return False, 'no delta_CP control for %s' % name
        pairs += 1
        r = d['us_per_point']['mean']/b['us_per_point']['mean']
        if r > 1.5 or r < 0.67:
            notes.append('%s %s %.0fx' % (d['code'], d['protocol']['grid'], r))
    if not pairs:
        return False, 'no Dmsq31 scan measured; only delta_CP, his fast path'
    return True, ('%d paired scans; parameter-sensitive: %s'
                  % (pairs, ', '.join(notes) if notes
                     else 'none beyond 50 per cent'))


CHECKS = [
    ('LBL-1  batches over energy', check_lbl1),
    ('LBL-2  Newton sweep, his conventions', check_lbl2),
    ('LBL-3  exact mode measured AND drawn', check_lbl3),
    ('Earth-1 shells: per region or total', check_earth1),
    ('Earth-2 both parameters scanned', check_earth2),
    ('Earth-3 realistic 100x100 grid', check_earth3),
    ('Earth-4 constant-density mode used', check_earth4),
    ('Other-1 pinned versions', check_other1),
    ('Other-2 repetitions and spread', check_other2),
]


def main():
    print('%-40s %-4s %s' % ('OBJECTION', 'OK', 'EVIDENCE'))
    print('-'*118)
    bad = 0
    for label, fn in CHECKS:
        try:
            ok, why = fn()
        except Exception as exc:                                # noqa: BLE001
            ok, why = False, '%s: %s' % (type(exc).__name__, exc)
        bad += not ok
        print('%-40s %-4s %s' % (label, 'yes' if ok else 'NO', why))
    print()
    print('One half of Earth-2 stays unmeasurable from an artifact: that the\n'
          '        scan hoists every invariant into setup() before any clock\n'
          '        starts.  That is a property of the protocol, not a number.\n'
          '        Which parameter moves inside the loop IS measured, above.')
    print()
    print('%d of %d checkable objections answered' % (len(CHECKS)-bad,
                                                      len(CHECKS)))
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
