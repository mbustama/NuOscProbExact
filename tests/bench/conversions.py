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
adapters, and ``tests/test_bench_pipeline.py`` checks each factor twice: once by
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
#: Not a code: an analytic expression evaluated inside this harness.  It has
#: no source to read a constant out of, so it takes THIS library's -- which
#: is also what makes its point fair, since it is our own series judged
#: against our own reference in our own units.
ALIAS = {'Second-order expansion': 'NuOscProbExact'}


def _resolve(code):
    r"""Maps an alias to the code whose conventions it borrows."""
    return ALIAS.get(code, code)


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


#: Where each code's hbar c lives in its pinned source, and in what form.
#:
#: This exists because `hbar_c_cosine_scale` took 1.97327e-7 as a default
#: argument that no caller ever overrode, while its own docstring asserted
#: that all three compiled Earth codes carry that value.  Prob3++ does not:
#: `mosc.c` hard-codes ``LoEfac = 2.534``, four significant figures, implying
#: 1.9731650e-7.  The difference is 5.3e-5 relative -- three orders of
#: magnitude above the 4.2e-8 the manifest calls the old figure's deliberate
#: cap -- and it put a measured 3.4e-4 floor on Prob3++ that was pure
#: convention rather than algorithm.
#:
#: Three forms appear, so each is converted rather than assumed:
#:   'ev_m'   the value IS hbar c in eV m;
#:   'ev_km'  hbar c in eV km, as GLoBES stores it;
#:   'loefac' the combination ``1e-9/hbar_c*1e3/2``, as Prob3++ and nuCraft
#:            store it, from which hbar c is 1e-6/(2*value).
HBAR_C_SITES = {
    'NuFast-LBL':   ('nufast-lbl/NuFast_LBL.cpp', 'ev_m',
                     r'eVsqkm_to_GeV_over4\s*=\s*1e-9\s*/\s*([0-9.eE+-]+)'),
    'NuFast-Earth': ('nufast-earth/src/Oscillation.cpp', 'ev_m',
                     r'eVsqkm_to_GeV_over2\s*=\s*1\.e-9\s*/\s*([0-9.eE+-]+)'),
    'GLoBES':       ('globes-3.2.18/globes/globes.h', 'ev_km',
                     r'#else\s*#define GLB_EV_TO_KM_FACTOR\s+([0-9.eE+-]+)'),
    'Prob3++':      ('prob3/mosc.c', 'loefac',
                     r'const double LoEfac\s*=\s*([0-9.eE+-]+)'),
    'nuCraft':      ('nucraft/NuCraft.py', 'loefac',
                     r'exp\(-([0-9.]+)j/en'),
}

#: nuCraft's neutral-current entry, the sibling of the charged-current one in
#: SITES.  Extracted so `reference_conventions.never_absorbed` has something
#: to compute from: the manifest calls their ratio (0.50161 where isoscalar
#: is exactly 1/2) an error rather than a convention, and a reference that
#: absorbed it would forgive the defect.
NUCRAFT_NC_SITE = ('nucraft/NuCraft.py',
                   r'self\.A\s*=\s*array\(\[[0-9.eE+-]+\*self\.y\[0\],\s*0\.,\s*0\.\]'
                   r'\+\[([0-9.eE+-]+)\*')


def _read(rel, pattern, what):
    path = os.path.join(BUILD, rel)
    with open(path) as handle:
        found = re.search(pattern, handle.read())
    if not found:
        raise RuntimeError('%s: %s not found in %s' % (what, pattern, path))
    return float(found.group(1)), path


def hbar_c(code):
    r"""Returns `code`'s own hbar c in eV m, read from its pinned source.

    Never a default argument.  A per-code reference has to be built in the
    code's own units, and a code that rounds hbar c to four figures is
    solving a measurably different problem from one that does not.
    """
    code = _resolve(code)
    rel, form, pattern = HBAR_C_SITES[code]
    value, _ = _read(rel, pattern, code + ' hbar c')
    if form == 'ev_m':
        return value
    if form == 'ev_km':
        return value*1.0e3
    if form == 'loefac':
        return 1.0e-6/(2.0*value)
    raise RuntimeError('%s: unknown hbar c form %r' % (code, form))


def km_to_inv_ev(code):
    r"""Returns a kilometre in eV^-1 in `code`'s own convention."""
    code = _resolve(code)
    if code == 'nuSQuIDS':
        import nuSQuIDS as nsq                      # its constants are its own
        units = nsq.Const()
        return units.km*units.eV
    if code == 'NuOscProbExact':
        return OURS['km_to_inv_ev']
    return 1.0e3/hbar_c(code)


def nucraft_nc_constant():
    r"""Returns nuCraft's neutral-current matter entry from its own source."""
    rel, pattern = NUCRAFT_NC_SITE
    value, _ = _read(rel, pattern, 'nuCraft NC entry')
    return value


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


def matter_constant(code):
    r"""Returns `code`'s matter constant in the shared ``A = k rhoYe / 2e9``
    normalisation, read from that code's own source or its own header.

    Five codes carry a literal, so `SITES` reads it.  Two do not:

    * nuSQuIDS builds it at run time as
      ``sqrt(2) G_F N_A cm^-3`` (`src/nuSQuIDS.cpp`, ``HI_constants``, used in
      ``HI()`` as ``CC = HI_constants * rho * Ye``).  Every factor comes from
      its own ``Const()``, so this is its value and not a transcription of it.
      As a check that this reproduces the C++: the result lands within 1e-8 of
      NuFast-Earth's independently written literal, as it should, since both
      use the modern CODATA G_F and N_A.
    * this library's own is derived in `_ours_matter_constant`.

    Without this the reference builder had nothing to construct nuSQuIDS'
    potential from --- ``mass_defect('nuSQuIDS')`` raised ``KeyError``.
    """
    code = _resolve(code)
    if code == 'nuSQuIDS':
        import math
        import nuSQuIDS as nsq
        c = nsq.Const()
        # sqrt2 is a mathematical constant, not one of theirs to round.
        return math.sqrt(2.0)*c.GF*c.Na*c.cm**-3.0*2.0e9
    if code == 'NuOscProbExact':
        return OURS['matter']
    theirs, _ = extract(code)
    if code == 'GLoBES':
        # GLB_V_FACTOR is V per (g/cm^3) where the others carry A = 2 E V per
        # GeV; the bridge is 2 (from A's definition) times 1e9 (GeV to eV),
        # both exact unit factors rather than measured constants.
        theirs = 2.0e9*theirs
    return theirs


def nc_over_cc(code):
    r"""Returns the neutral-current to charged-current ratio `code` applies.

    Only meaningful at four flavors: among three active flavors the NC term
    is common to all of them and cancels out of any oscillation probability.

    For isoscalar matter with Y_e = 1/2 the exact value is -1/2.  nuCraft is
    the one code that does not reach it -- its two entries are independently
    rounded and their ratio is 0.50161 -- and per the manifest that is an
    error rather than a convention, so its reference is built with the exact
    ratio while keeping its constants for the units.  A reference that
    absorbed the rounding would forgive the defect.
    """
    if code == 'nuCraft':
        return -0.5          # deliberately NOT nucraft_nc_constant()/its CC
    return -0.5


def mass_defect(code):
    r"""Returns ours/theirs for `code`, the factor its density must carry.

    GLoBES's constant is in a different convention from the others'
    ``YerhoE2a``: GLB_V_FACTOR is the potential V in eV per (g/cm^3) of
    electron density, where the others carry :math:`A = 2 E V` per GeV.
    The bridge is :math:`2 \times 10^9` -- the 2 from the definition of A
    and the :math:`10^9` from GeV to eV, both exact unit factors rather
    than measured constants.
    """
    return OURS['matter']/matter_constant(code)


def globes_ne_mantle():
    r"""Returns the electron fraction GLoBES multiplies in by itself.

    ``glb_probability.c`` computes ``density * GLB_V_FACTOR *
    GLB_Ne_MANTLE``, so a density handed to GLoBES must divide our
    :math:`Y_e` by this rather than assume the two are equal.
    """
    value, _ = extract('GLoBES-Ne')
    return value


def hbar_c_cosine_scale(their_hbar_c_ev_m=None, code=None):
    r"""Returns the factor a chord's cosine carries to match a code's length.

    All three compiled Earth codes hard-code the same
    :math:`\hbar c = 1.97327 \times 10^{-7}` eV m, which makes a kilometre
    5.0677302143e9 eV^-1 against this library's 5.06773e9.  The chord is
    :math:`L = -2 R \cos\theta_z`, so the mismatch is absorbed by handing them
    a cosine shrunk by this factor rather than by changing the geometry.
    """
    if their_hbar_c_ev_m is None and code is None:
        raise RuntimeError(
            'hbar_c_cosine_scale needs a code: the three compiled Earth codes '
            'do NOT share one hbar c -- Prob3++ carries 2.534, implying '
            '1.9731650e-7 against the others\' 1.97327e-7.  Pass code=... so '
            'the value is read from that code\'s own pinned source.')
    theirs_km = (1.0e3/their_hbar_c_ev_m if their_hbar_c_ev_m is not None
                 else km_to_inv_ev(code))
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
    code = _resolve(code)
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
    sterile = p.get('sterile_3plus1')
    if sterile:
        out['dmsq41_ev2'] = sterile['dmsq41_ev2']
        for k in ('14', '24', '34'):
            out['s%ssq' % k] = sterile['s%ssq' % k]
            out['theta%s_rad' % k] = math.asin(math.sqrt(sterile['s%ssq' % k]))
    return out
