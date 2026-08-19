# -*- coding: utf-8 -*-
r"""NuOscProbExact -- this library, driven through its batched entry points.

Everything invariant under the delta_CP scan is built in ``setup()``: the
energy and zenith grids, the mass-defect-free matter potential (it is our
own convention), and the knob decode.  ``configure()`` rebuilds the vacuum
Hamiltonian, because delta_CP lives inside it -- that is the per-step cost a
fit actually pays with this library.  ``evaluate()`` makes ONE call for the
whole grid: ``earth.probabilities_3nu_earth`` batches over energies AND
zenith angles together, and the constant-density path stacks the
Hamiltonians and hands them to ``oscprob3nu.probabilities_3nu`` in one call.

The knob is one integer carrying both precision dials, documented in
``capabilities()``:

* ``knob > 0``  -- ``n_slabs_per_segment = knob`` (the slab dial);
* ``knob < 0``  -- ``rtol = 10**knob`` (the tolerance dial, reachable for
  knob in -3..-5; the slab search starts at 1 so the tolerance, not the
  start, picks the count);
* ``knob == 0`` -- the library default (``n_slabs_per_segment = 8``).

Oscillation parameters arrive through the Problem, which runner.py fills
from conversions.py; nothing physical is typed here.
"""

import math
import os
import sys

_REPO = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(_REPO, 'src'))

import numpy as np                                             # noqa: E402

import globaldefs as gd                                        # noqa: E402
import earth                                                   # noqa: E402
import hamiltonians3nu                                         # noqa: E402
import oscprob3nu                                              # noqa: E402

#: P(nu_mu -> nu_mu) in the flattened 3x3 the library returns.
_PMM = 4


class NuOscProbExact(object):

    name = 'NuOscProbExact'

    def capabilities(self):
        return {
            'batches_energy': True,
            'batches_zenith': True,
            'batch_symbol': ('earth.probabilities_3nu_earth(E array, costhz '
                             'array) / oscprob3nu.probabilities_3nu(H stack)'),
            'knob_name': 'n_slabs_per_segment(knob>0) | rtol=10^knob(knob<0)',
            # The slab dial runs the full 1..256.  The tolerance dial stops
            # at 1e-5, which is where it actually stops: on CHORD/12x1 the
            # slab product's error estimate bottoms out near 7e-8 at the
            # 1024-slab ceiling, so rtol=1e-6 and tighter raise rather than
            # return.  The manifest used to advertise 1e-3..1e-15; ten of
            # those thirteen settings crash, and declaring a domain a sweep
            # cannot survive is the mirror of the knob LBL-3 was about.
            # Reach the small residuals through the slab dial instead.
            'knob_domain': [1, 2, 4, 8, 16, 32, 64, 128, 256, -3, -4, -5],
        }

    def setup(self, problem):
        p = problem
        self._e_ev = np.asarray(p.energies_gev, dtype=float)*1.0e9
        self._costhz = (np.asarray(p.costhz, dtype=float)
                        if p.costhz else None)
        self._ye = p.ye
        self._sines = (math.sqrt(p.s12sq), math.sqrt(p.s23sq),
                       math.sqrt(p.s13sq))
        self._dm21, self._dm31 = p.dm21, p.dm31

        knob = int(p.knob)
        if knob > 0:
            self._n_slabs, self._rtol = knob, None
        elif knob < 0:
            self._n_slabs, self._rtol = 1, 10.0**knob
        else:
            self._n_slabs, self._rtol = 8, None    # the library default

        # Constant-density pieces, invariant under the scan: the potential
        # and the baseline in natural units, both from the library itself.
        self._vcc = earth.matter_potential(p.density, p.ye)
        self._L_inv_ev = p.L_km*gd.CONV_KM_TO_INV_EV

        self.configure(p.dcp)

    def configure(self, dcp):
        s12, s23, s13 = self._sines
        self._h_vac = np.asarray(
            hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
                s12, s23, s13, dcp, self._dm21, self._dm31))

    def reset(self):
        r"""Nothing to reset, and the one cache that survives is on the right
        side of the line.

        ``evaluate()`` rebuilds the matter Hamiltonians and the whole slab
        product on every call.  The chord geometry does NOT get rebuilt:
        ``earth._earth_slabs_cached`` is an ``lru_cache`` keyed on
        ``(costhz, n_slabs_per_segment)`` alone, so a second call reuses it
        --- verified by watching ``cache_info()`` register hits.

        That is deliberate and is the same side of the setup/request line
        every other code sits on: it is pure geometry, independent of the
        oscillation parameters, and it corresponds to GLoBES' chord
        decomposition, Prob3++'s profile file and NuFast-Earth's
        ``trajectories_calculated`` --- all built in ``setup`` and all left
        standing by ``reset``.  Nothing oscillation-dependent survives a
        repetition here, which is the property ``reset`` exists to enforce."""

    def evaluate(self):
        if self._costhz is not None:
            # A scalar angle takes the whole energy array down one chord;
            # arrays broadcast, so the oscillogram is asked for the way the
            # docstring says: energies and angles on different axes.
            if self._costhz.size == 1:
                energy, costhz = self._e_ev, float(self._costhz[0])
            else:
                energy, costhz = self._e_ev[None, :], self._costhz[:, None]
            probs = earth.probabilities_3nu_earth(
                self._h_vac, energy, costhz,
                n_slabs_per_segment=self._n_slabs,
                electron_fraction=self._ye, rtol=self._rtol)
            return float(np.sum(np.asarray(probs)[..., _PMM]))
        h = np.asarray(hamiltonians3nu.hamiltonian_3nu_matter(
            self._h_vac, self._e_ev, self._vcc))
        probs = oscprob3nu.probabilities_3nu(h, self._L_inv_ev)
        return float(np.sum(np.asarray(probs)[..., _PMM]))


ADAPTER = NuOscProbExact
