# -*- coding: utf-8 -*-
r"""The one place any physical conversion factor is computed.

Every code carries its own rounding of the matter constant and its own
:math:`\hbar c`, so comparing them requires converting one shared physical
input into each code's convention -- exactly, and once.  Doing it per driver
produced three different values of a single NuFast-Earth factor: 0.99209238
applied in the driver, 0.9920928368 implied by its own comment, and 0.9920938
printed in the README.

Each constant below is read from the *pinned* source in ``.bench-build/src``
rather than transcribed, so a factor cannot drift from the code it converts
for.  ``gen_conversions.py`` turns this into ``conversions.h`` for the C++
adapters, and ``test_bench_conversions.py`` checks each factor twice: once by
derivation, and once against a vacuum residual.
"""

import os
import re

BUILD = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), '.bench-build', 'src')

def _ours_matter_constant():
    r"""Returns this library's matter constant in the codes' normalisation.

    Derived from :mod:`globaldefs`, never transcribed.  The four external
    literals express the potential as ``V = k * rhoYe / 2e9`` with ``rhoYe`` in
    g cm^-3, so the same quantity for this library is
    ``sqrt(2) G_F N_e * 2e9 / rhoYe`` at any density.

    This function exists because the first version of this module opened with
    ``'matter': 1.514423e-4`` --- a seven-figure transcription of the value from
    an old driver's comment --- which is the very sin the module was written to
    prevent.  It was 2.4e-8 off, and every factor derived from it inherited that
    error.  Caught by checking a derived factor against the potential computed
    from the library's own constants, which is now a test.
    """
    import math
    import os
    import sys

    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), 'src'))
    import globaldefs as gd

    rho_ye = gd.DENSITY_MATTER_CRUST_G_PER_CM3*gd.ELECTRON_FRACTION_EARTH_CRUST
    potential = math.sqrt(2.0)*gd.GF*gd.NUM_DENSITY_E_EARTH_CRUST
    return potential*2.0e9/rho_ye


#: This library's own values.  The matter constant is derived; the length
#: conversion is `globaldefs.CONV_KM_TO_INV_EV` itself rather than a copy.
def _ours_km_to_inv_ev():
    import os
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), 'src'))
    import globaldefs as gd
    return gd.CONV_KM_TO_INV_EV


OURS = {
    'matter': _ours_matter_constant(),
    'km_to_inv_ev': _ours_km_to_inv_ev(),
}

#: Where each external constant lives in its pinned source.
#:
#: GLoBES defines GLB_V_FACTOR twice, old value behind ``#ifdef
#: GLB_OLD_CONSTANTS`` and current value in the ``#else`` branch; the
#: pattern anchors on the ``#else`` so it reads the one that compiles.
#: GLB_Ne_MANTLE has a commented-out 0.497 on the line above the live 0.5,
#: so its pattern requires the ``#define`` to start the line.
#: nuCraft's constant is the CC entry of EarthModel's potential vector,
#: ``A = array([15.256e-5*self.y[0], ...])`` -- same YerhoE2a convention as
#: the NuFast codes, with Y_e carried separately by ``y[0]``.
SITES = {
    'NuFast-LBL':   ('nufast-lbl/NuFast_LBL.cpp',      r'YerhoE2a\s*=\s*([0-9.eE+-]+)'),
    'NuFast-Earth': ('nufast-earth/src/Oscillation.cpp', r'YerhoE2a\s*=\s*([0-9.eE+-]+)'),
    'Prob3++':      ('prob3/mosc.c',                   r'tworttwoGf\s*=\s*([0-9.eE+-]+)'),
    'GLoBES':       ('globes-3.2.18/source/glb_probability.h',
                     r'#else\s*#define GLB_V_FACTOR\s+([0-9.eE+-]+)'),
    'GLoBES-Ne':    ('globes-3.2.18/source/glb_probability.h',
                     r'\n#define GLB_Ne_MANTLE\s+([0-9.eE+-]+)'),
    'nuCraft':      ('nucraft/NuCraft.py',
                     r'self\.A\s*=\s*array\(\[([0-9.eE+-]+)\*self\.y\[0\]'),
}


def extract(code):
    r"""Returns the matter constant a code carries, read from its own source.

    Raises
    ------
    RuntimeError
        If the constant is not where the manifest says it is -- which is the
        signal that the pin moved and the factor must be re-derived.
    """
    rel, pattern = SITES[code]
    path = os.path.join(BUILD, rel)
    with open(path) as handle:
        text = handle.read()
    found = re.search(pattern, text)
    if not found:
        raise RuntimeError('%s: %s not found in %s' % (code, pattern, path))
    return float(found.group(1)), path


def mass_defect(code):
    r"""Returns ours/theirs for `code`, the factor its density must carry.

    GLoBES's constant is in a different convention from the others'
    ``YerhoE2a``: GLB_V_FACTOR is the potential V in eV per (g/cm^3) of
    electron density, where the others carry :math:`A = 2 E V` per GeV.
    The bridge is :math:`2 \times 10^9` -- the 2 from the definition of A
    and the :math:`10^9` from GeV to eV, both exact unit factors rather
    than measured constants.
    """
    theirs, _ = extract(code)
    if code == 'GLoBES':
        theirs = 2.0e9*theirs
    return OURS['matter']/theirs


def globes_ne_mantle():
    r"""Returns the electron fraction GLoBES multiplies in by itself.

    ``glb_probability.c`` computes ``density * GLB_V_FACTOR *
    GLB_Ne_MANTLE``, so a density handed to GLoBES must divide our
    :math:`Y_e` by this rather than assume the two are equal.
    """
    value, _ = extract('GLoBES-Ne')
    return value


def hbar_c_cosine_scale(their_hbar_c_ev_m=1.97327e-7):
    r"""Returns the factor a chord's cosine carries to match a code's length.

    All three compiled Earth codes hard-code the same
    :math:`\hbar c = 1.97327 \times 10^{-7}` eV m, which makes a kilometre
    5.0677302143e9 eV^-1 against this library's 5.06773e9.  The chord is
    :math:`L = -2 R \cos\theta_z`, so the mismatch is absorbed by handing them
    a cosine shrunk by this factor rather than by changing the geometry.
    """
    theirs_km = 1.0e3/their_hbar_c_ev_m
    return OURS['km_to_inv_ev']/theirs_km


def _manifest():
    import json
    here = os.path.dirname(os.path.abspath(__file__))
    return json.load(open(os.path.join(here, 'manifest.json')))


def oscillation_parameters():
    r"""Returns the one parameter set every code is handed.

    Read from the manifest, so there is no second copy to drift.  Values are
    returned in this library's parameterisation; `for_code` converts.
    """
    return dict(_manifest()['oscillation_parameters'])


def for_code(code):
    r"""Returns the same physics in `code`'s own parameterisation.

    Prob3++ takes :math:`\Delta m^2_{32}` where the others take
    :math:`\Delta m^2_{31}`, and codes that want angles get
    :math:`\theta = \arcsin\sqrt{\sin^2\theta}`.  Both are *derived* here
    rather than typed a second time, because a hand-typed 2.4511e-3 beside a
    2.525e-3 is indistinguishable from a physics difference.
    """
    import math

    p = oscillation_parameters()
    out = {
        's12sq': p['s12sq'], 's13sq': p['s13sq'], 's23sq': p['s23sq'],
        'dcp_rad': p['dcp_deg']*math.pi/180.0,
        'dmsq21_ev2': p['dmsq21_ev2'], 'dmsq31_ev2': p['dmsq31_ev2'],
    }
    if code == 'Prob3++':
        out['dmsq32_ev2'] = p['dmsq31_ev2'] - p['dmsq21_ev2']
    if code in ('nuCraft',):
        for k in ('12', '13', '23'):
            out['theta%s_rad' % k] = math.asin(math.sqrt(p['s%ssq' % k]))
    return out
