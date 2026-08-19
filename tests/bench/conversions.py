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

#: This library's own values, from :mod:`globaldefs`.
OURS = {
    'matter': 1.514423e-4,      # 2 sqrt(2) G_F N_A, in the codes' units
    'km_to_inv_ev': 5.06773e9,
}

#: Where each external constant lives in its pinned source.
SITES = {
    'NuFast-LBL':   ('nufast-lbl/NuFast_LBL.cpp',      r'YerhoE2a\s*=\s*([0-9.eE+-]+)'),
    'NuFast-Earth': ('nufast-earth/src/Oscillation.cpp', r'YerhoE2a\s*=\s*([0-9.eE+-]+)'),
    'Prob3++':      ('prob3/mosc.c',                   r'tworttwoGf\s*=\s*([0-9.eE+-]+)'),
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
    r"""Returns ours/theirs for `code`, the factor its density must carry."""
    theirs, _ = extract(code)
    return OURS['matter']/theirs


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
