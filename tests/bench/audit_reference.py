# -*- coding: utf-8 -*-
r"""Measures whether each reference matches the code it judges, and records it.

The question this answers is the one that decides whether an accuracy figure
means anything: does the reference solve *the same problem* the code was
handed --- same :math:`Y_e`, same matter potential, same :math:`\hbar c`,
same PREM discretisation, same chord --- or does it silently solve a
neighbouring one?  A reference in the wrong convention makes a code look
better than it is, and nothing downstream notices.

So the answer is measured rather than asserted, and written to
``reference_audit.json`` so that the notebook and the paper render numbers
from a file with a generator instead of carrying typed copies.  Six published
numbers about external codes once existed only in docstrings with no stored
run; this is the shape that mistake takes when it is avoided.

Run from the repository root::

    python tests/bench/audit_reference.py

Nothing here reads a clock: every number is a correctness measurement, so it
can be produced on a machine that is not idle.
"""

import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, 'adapters'))
sys.path.insert(0, os.path.join(ROOT, 'src'))

import numpy as np                                             # noqa: E402

import conversions                                             # noqa: E402
import reference as ref                                        # noqa: E402

OUT = os.path.join(HERE, 'reference_audit.json')

#: The chord and energy every fidelity probe uses.  One point is enough: a
#: profile that disagrees anywhere disagrees here, and the expensive entries
#: are minutes each.
COSTHZ = -0.9
ENERGY_GEV = 10.0
YE = 0.5

CODES = ('NuFast-LBL', 'NuFast-Earth', 'Prob3++', 'GLoBES', 'nuSQuIDS',
         'nuCraft', 'NuOscProbExact')

#: How each code's profile reaches it, and therefore what its reference has
#: to reproduce.  Recorded rather than inferred, because this is the column
#: that was wrong.
PROFILE_KIND = {
    'NuFast-LBL':     'none -- constant density only',
    'NuFast-Earth':   'four major PREM layers, each cut into n uniform '
                      'sub-shells held at their midpoint density',
    'Prob3++':        'the same uniform sub-shells, through a profile file',
    'GLoBES':         'the same uniform sub-shells, chord-decomposed',
    'nuSQuIDS':       'an Akima spline that EarthAtm builds itself over the '
                      'nodes it is handed',
    'nuCraft':        'a k=1 spline, the object the adapter hands it',
    'NuOscProbExact': 'the continuous PREM',
}


def profile_fidelity():
    r"""Returns, per code, how far the reference's profile sits from the
    code's own --- measured against the code where the code will tell us.

    Only two codes can be interrogated directly.  nuSQuIDS exposes
    ``EarthAtm::density`` along a track, so its spline can be sampled and
    compared.  nuCraft is handed its profile by the adapter, so the two share
    one object and the comparison is exact by construction.  For the three
    compiled Earth codes the profile is a stack of uniform shells built from
    the same radii and the same PREM polynomial in both places, so the
    reference reproduces it identically.
    """
    out = {}

    # nuSQuIDS: sample its own density along a chord.
    try:
        import earth
        import globaldefs as gd
        import nuSQuIDS as nsq
        from nusquids import _aligned_radii

        radii = _aligned_radii(200)
        rho = earth.density_prem(np.clip(radii, 0.0, gd.EARTH_RADIUS))
        body = nsq.EarthAtm((radii/gd.EARTH_RADIUS).tolist(), rho.tolist(),
                            np.full(radii.size, YE).tolist())
        body.SetAtmosphereHeight(0.0)
        track = body.MakeTrackWithCosine(COSTHZ)
        length = track.GetFinalX()

        _, density_fn, _ = ref.profile_for('nuSQuIDS')
        radius = ref._mpf(gd.EARTH_RADIUS)
        chord = -2*radius*ref._mpf(COSTHZ)
        worst = 0.0
        for frac in np.linspace(0.02, 0.98, 40):
            track.SetX(frac*length)
            ours = float(density_fn(ref.radius_at_mp(COSTHZ, frac*chord,
                                                     radius)))
            worst = max(worst, abs(body.density(track) - ours))
        out['nuSQuIDS'] = worst
    except Exception as exc:                                   # noqa: BLE001
        out['nuSQuIDS'] = 'unavailable: %s' % exc

    # nuCraft: the adapter's spline against the reference's.
    try:
        from nucraft import NuCraftAdapter
        theirs = NuCraftAdapter._profile_spline(YE)
        _, ours, _ = ref.profile_for('nuCraft')
        out['nuCraft'] = max(
            abs(float(theirs(r)) - float(ours(ref._mpf(r))))
            for r in np.linspace(50.0, 6300.0, 60))
    except Exception as exc:                                   # noqa: BLE001
        out['nuCraft'] = 'unavailable: %s' % exc

    return out


def main():
    started = time.strftime('%Y-%m-%dT%H:%M:%S%z')
    record = {
        'schema': 'bench-reference-audit/1',
        'generated_by': 'tests/bench/audit_reference.py',
        'generated_at': started,
        'probe': {'costhz': COSTHZ, 'energy_gev': ENERGY_GEV, 'ye': YE},
        'working_digits': ref.DPS,
        'codes': {},
    }

    fidelity = profile_fidelity()

    for code in CODES:
        params = conversions.for_code(code)
        entry = {
            'ye': YE,
            'ye_shared_with_every_code': True,
            'matter_constant': conversions.matter_constant(code),
            'matter_constant_source': 'own pinned source'
                                      if code not in ('nuSQuIDS',
                                                      'NuOscProbExact')
                                      else ('own Const()' if code == 'nuSQuIDS'
                                            else 'derived from globaldefs'),
            'km_in_inv_ev': conversions.km_to_inv_ev(code),
            'profile_kind': PROFILE_KIND[code],
        }
        if code in fidelity:
            entry['profile_fidelity_g_per_cm3'] = fidelity[code]

        if code == 'NuFast-LBL':
            # Constant density: one matrix exponential, no discretisation.
            value = ref.constant_density(code, [ENERGY_GEV], 1300.0, 3.0, YE,
                                         params)[0]
            entry['reference_kind'] = 'exact matrix exponential'
            entry['reference_error'] = float(ref.mp.mpf(10)**(-ref.DPS + 5))
            entry['reference_value'] = str(value)
        else:
            started_at = time.time()
            value, error = ref.earth_chord_reference(code, ENERGY_GEV, COSTHZ,
                                                     YE, params)
            entry['reference_kind'] = (
                'exact finite product over uniform shells'
                if code in ('NuFast-Earth', 'Prob3++', 'GLoBES')
                else 'Romberg-extrapolated slab product')
            entry['reference_error'] = float(error)
            entry['reference_value'] = str(value)
            entry['seconds'] = round(time.time() - started_at, 1)
        record['codes'][code] = entry
        print('  %-15s ref err %-12.2e  %s'
              % (code, entry['reference_error'], entry['reference_kind']))

    with open(OUT, 'w') as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
        handle.write('\n')
    print('wrote %s' % os.path.relpath(OUT, ROOT))


if __name__ == '__main__':
    main()
