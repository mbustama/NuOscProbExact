# -*- coding: utf-8 -*-
r"""nuCraft r22, from the pinned tarball under ``.bench-build/src/nucraft``.

``CalcWeights`` takes a particle list but loops over it internally, so this
adapter is honest about what that buys: ``batches_energy = False`` with
``batch_symbol = "interface-only"``.  The whole (energy, zenith) grid is
still handed over in ONE call, so that no code in these figures is looped
from outside when its own interface would take the stack.

delta_CP enters nuCraft through the constructor's angle list -- there is no
setter -- so ``configure()`` rebuilds the ``NuCraft`` instance around the
hoisted ``EarthModel``.  That reconstruction is the per-step cost this
code's interface charges a delta_CP scan, and the AMORTIZED protocol is
defined to include it.

Conventions, all derived rather than typed:

* Angles.  nuCraft takes mixing ANGLES in degrees, not sin^2; they come
  from ``conversions.for_code('nuCraft')``, which derives
  theta = arcsin(sqrt(sin^2 theta)) from the manifest.  The CP phase rides
  on the (1,3) tuple, in degrees.
* Density.  Not rescaled.  nuCraft's CC matter entry (NuCraft.py:212) is
  absorbed into its own reference instead.
* Profile.  nuCraft's file route is Python 2 only, so this library's PREM
  goes in through ``EarthModel.profInt``, exactly as tests/prem_scan.py
  does; the inner-core radius comes from ``earth.PREM_BOUNDARIES``.

The knob is ``numPrec = 10**-knob`` (manifest domain 1e-2..1e-10, so knob
in 2..10); ``knob <= 0`` leaves nuCraft's own default.  Constant-density
problems are refused: nuCraft propagates through its Earth model only.
"""

import math
import os
import sys
import warnings

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(os.path.dirname(_HERE)))
sys.path.insert(0, os.path.join(_REPO, 'src'))
sys.path.insert(0, os.path.join(_REPO, '.bench-build', 'src', 'nucraft'))
sys.path.insert(0, os.path.dirname(_HERE))                # for conversions

import numpy as np                                             # noqa: E402

import conversions                                             # noqa: E402
import globaldefs as gd                                        # noqa: E402
import earth                                                   # noqa: E402

# Importing nuCraft installs a `warnings.showwarning` replacement with the
# Python 2 signature, which raises TypeError the moment nuCraft actually
# warns.  Keep the standard one, as tests/prem_scan.py does.
_SHOWWARNING = warnings.showwarning
import NuCraft as _nucraft_module                              # noqa: E402
from NuCraft import NuCraft, EarthModel                        # noqa: E402
warnings.showwarning = _SHOWWARNING

# r22 opens with `from numpy import *`, and numpy >= 2 exports `bool` again,
# so inside NuCraft.py the name `bool` is numpy.bool_ -- which makes its
# `assert type(vacuum) is bool` reject a genuine Python False.  Handing it
# an instance of ITS OWN idea of bool satisfies the check under any numpy,
# without patching the pinned source.
_NUCRAFT_FALSE = _nucraft_module.__dict__.get('bool', bool)(False)


class NuCraftAdapter(object):

    name = 'nuCraft'

    def capabilities(self):
        return {
            'batches_energy': False,
            'batches_zenith': False,
            'batch_symbol': ('interface-only -- CalcWeights takes a '
                             'particle list but loops over it internally'),
            'knob_name': 'log10_numPrec',
            'knob_domain': [2, 3, 4, 5, 6, 7, 8, 9, 10],
        }

    def setup(self, problem):
        p = problem
        if not p.costhz:
            raise SystemExit('nuCraft has no constant-density mode; it '
                             'propagates through its Earth model only')

        em = EarthModel(
            'prem', y=(p.ye, p.ye, p.ye),
            rICore=float(np.asarray(earth.PREM_BOUNDARIES, dtype=float)[0]))
        # nuCraft's file route is Python 2 only, so the profile goes in here.
        # No density rescaling: nuCraft's own matter constants are absorbed
        # into ITS reference (with the isoscalar ratio forced to exactly 1/2,
        # because its 0.50161 is an error rather than a convention).
        #
        # profInt is called from inside nuCraft's ODE right-hand side, once
        # per solver step.  A Python lambda around this library's scalar
        # density_prem -- which range-checks and searchsorts per call -- put
        # our own cost inside its inner loop, which is the same mistake that
        # inflated another code's published point.  nuCraft builds itself an
        # InterpolatedUnivariateSpline(k=1); it gets the identical object
        # type here, over THIS library's PREM, with the shell boundaries
        # carried as near-duplicate knots so no interpolated segment
        # straddles a density jump.
        em.profInt = self._profile_spline(p.ye)
        self._em = em

        c = conversions.for_code('nuCraft')
        self._deg = {k: math.degrees(c['theta%s_rad' % k])
                     for k in ('12', '13', '23')}
        # nuCraft's masses list is [m1, Dm21, Dm31]; the absolute scale m1
        # does not enter oscillations.
        self._masses = [1.0, p.dm21, p.dm31]

        # atmMode 0 with atmHeight 0 is surface to surface; CalcWeights
        # takes the ZENITH ANGLE in radians, not its cosine.  One row per
        # (energy, zenith) pair, 14 = nu_mu, energies in GeV.
        self._rows = [(14, float(e), float(np.arccos(cz)))
                      for cz in p.costhz for e in p.energies_gev]
        self._num_prec = 10.0**-p.knob if p.knob > 0 else None

        self.configure(p.dcp)

    @staticmethod
    def _profile_spline(ye, n_per_shell=200):
        r"""This library's PREM as the spline type nuCraft builds for itself."""
        from scipy.interpolate import InterpolatedUnivariateSpline
        edges = np.concatenate(([0.0], np.asarray(earth.PREM_BOUNDARIES,
                                                  dtype=float)))
        edges = edges[edges <= gd.EARTH_RADIUS]
        if edges[-1] < gd.EARTH_RADIUS:
            edges = np.concatenate((edges, [gd.EARTH_RADIUS]))
        eps = 1.0e-9                      # km, to carry the jump
        pieces = []
        for a, b in zip(edges[:-1], edges[1:]):
            inner = np.linspace(a, b, n_per_shell)
            inner[0], inner[-1] = a + eps, b - eps
            pieces.append(inner)
        r = np.unique(np.concatenate(([0.0], np.concatenate(pieces),
                                      [gd.EARTH_RADIUS])))
        rho = earth.density_prem(np.clip(r, 0.0, gd.EARTH_RADIUS))
        return InterpolatedUnivariateSpline(r, rho, k=1)

    def configure(self, dcp):
        angles = [(1, 2, self._deg['12']),
                  (1, 3, self._deg['13'], math.degrees(dcp)),
                  (2, 3, self._deg['23'])]
        self._nc = NuCraft(self._masses, angles, earthModel=self._em,
                           detectorDepth=0.0, atmHeight=0.0)

    def reset(self):
        r"""Nothing to reset: ``CalcWeights`` re-integrates every particle on
        every call, so each repetition is already cold."""

    def evaluate(self):
        kwargs = {'atmMode': 0, 'vacuum': _NUCRAFT_FALSE}
        if self._num_prec is not None:
            kwargs['numPrec'] = self._num_prec
        weights = self._nc.CalcWeights(self._rows, **kwargs)
        return float(sum(w[1] for w in weights))   # P(nu_mu -> nu_mu)


ADAPTER = NuCraftAdapter
