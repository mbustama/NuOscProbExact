# -*- coding: utf-8 -*-
r"""The second-order expansion in :math:`\alpha` and :math:`s_{13}`.

The one point in the constant-density figure that had never been measured.
It arrived in ``tests/speed_accuracy.json`` as ``0.2495 us`` at ``6.51e-3``
with no generator, no version and no formula recorded, and survived every
regeneration since because nothing here could reproduce it.  The paper names
it -- "the second-order expansion in alpha and s13 ... not a numerical error
but the truncation of the series" -- and cites Cervera et al., Freund, and
Akhmedov et al.  This is that formula, so the point becomes a measurement
like its neighbours instead of a number nobody can regenerate.

THE CHANNEL.  This expansion is the APPEARANCE formula.  It gives
:math:`P(\nu_\mu \to \nu_e)` and does not, in this form, give the other two
entries of the row, so :meth:`capabilities` declares one channel and the
accuracy sweep scores that one.  Inventing the disappearance entries from a
different expansion, or closing the row by unitarity the series does not
respect, would put numbers in a figure that no published formula returns --
which is the failure this whole comparison exists to avoid.  The figure's y
axis already reads :math:`\max|\Delta P_{\nu_\mu \to \nu_e}|`, so scoring the
appearance channel is what it claims to show.

CONVENTIONS.  This is not an external code and has no constants of its own:
it is an analytic expression evaluated in THIS library's conventions, with
the same matter potential and the same baseline in inverse eV as every other
series in that figure.  Its reference is therefore this library's, and the
comparison is like for like.

THE PLATEAU IS REAL.  Its error does not fall with any setting, because
there is no setting -- the residual is the truncated third order, not
arithmetic.  It is one of only two honest plateaus in these figures, beside
NuFast-LBL running out of double precision.
"""

import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_REPO, 'src'))
sys.path.insert(0, os.path.dirname(_HERE))

import conversions                                              # noqa: E402


class SecondOrderAdapter(object):
    r"""P(nu_mu -> nu_e) to second order in alpha and s13."""

    name = 'Second-order expansion'

    def capabilities(self):
        return {
            'batches_energy': True,
            'batches_zenith': False,
            'batch_symbol': 'numpy expression over the energy array',
            'knob_name': 'none',
            'knob_domain': [0],
            # One channel, and the sweep is told so rather than being handed
            # two fabricated ones.
            'channels': ['numu->nue'],
        }

    def setup(self, problem):
        p = problem
        self._scan_dmsq31 = getattr(p, 'scan', 'dcp') == 'dmsq31'
        self._e_ev = np.asarray(p.energies_gev, dtype=float)*1.0e9
        self._dm21, self._dm31 = p.dm21, p.dm31
        self._dcp = p.dcp

        self._th12 = math.asin(math.sqrt(p.s12sq))
        self._th13 = math.asin(math.sqrt(p.s13sq))
        self._th23 = math.asin(math.sqrt(p.s23sq))

        # V_CC in eV, in this library's normalisation -- the same one its own
        # series in this figure is measured against.
        k = conversions.matter_constant('NuOscProbExact')
        self._vcc = k*p.density*p.ye/2.0e9
        # Baseline in inverse eV, so that Delta is dimensionless.
        self._l_inv_ev = p.L_km*conversions.km_to_inv_ev('NuOscProbExact')
        self._probs = None

    def configure(self, v):
        if self._scan_dmsq31:
            self._dm31 = v
        else:
            self._dcp = v

    def reset(self):
        r"""Nothing is cached between calls: the expression is closed form."""
        self._probs = None

    def _appearance(self):
        r"""Cervera et al., second order in alpha and s13.

        P = s23^2 sin^2(2 th13) sin^2[(1-A)D] / (1-A)^2
          + alpha sin(2th13) sin(2th12) sin(2th23) cos(D + dcp)
              [sin(A D)/A] [sin((1-A)D)/(1-A)]
          + alpha^2 cos^2(th23) sin^2(2 th12) sin^2(A D)/A^2

        with D = dm31 L / 4E, A = 2 E V / dm31 and alpha = dm21/dm31.  The
        two removable singularities at A -> 0 and A -> 1 are carried by
        np.sinc, which is exactly sin(pi x)/(pi x) and finite at zero, rather
        than by an epsilon guard that would bias the very resonance region
        this expansion is used in.
        """
        e = self._e_ev
        dm31, dm21 = self._dm31, self._dm21
        alpha = dm21/dm31
        D = dm31*self._l_inv_ev/(4.0*e)
        A = 2.0*e*self._vcc/dm31

        s2_13 = math.sin(2.0*self._th13)
        s2_12 = math.sin(2.0*self._th12)
        s2_23 = math.sin(2.0*self._th23)
        s23sq = math.sin(self._th23)**2
        c23sq = math.cos(self._th23)**2

        # sin(x D)/x = D sinc(x D / pi); finite as x -> 0 by construction.
        sin_AD_over_A = D*np.sinc(A*D/math.pi)
        sin_1mAD_over_1mA = D*np.sinc((1.0 - A)*D/math.pi)

        term1 = s23sq*s2_13**2*sin_1mAD_over_1mA**2
        term2 = (alpha*s2_13*s2_12*s2_23*np.cos(D + self._dcp)
                 * sin_AD_over_A*sin_1mAD_over_1mA)
        term3 = alpha**2*c23sq*s2_12**2*sin_AD_over_A**2
        return term1 + term2 + term3

    def evaluate(self):
        self._probs = self._appearance()
        return float(np.sum(self._probs))

    def probabilities(self):
        r"""The appearance channel only; see the module docstring."""
        if self._probs is None:
            self._probs = self._appearance()
        return [float(x) for x in np.asarray(self._probs).ravel()]


ADAPTER = SecondOrderAdapter
