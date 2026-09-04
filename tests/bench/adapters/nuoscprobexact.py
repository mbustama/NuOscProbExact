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

#: The nu_mu row of the flattened 3x3 the library returns: numu->nue,
#: numu->numu, numu->nutau.  A row rather than one channel because the
#: constant-density figure plots appearance while the accuracy sweep scores
#: disappearance, and an artifact holding one cannot serve the other.
#: Ceiling handed to the tolerance search, overriding `slabs.N_SLABS_MAX`.
#: The shipped default of 1024 is a guardrail against a mis-stated tolerance,
#: not a limit of the discretisation; raising it here leaves every other
#: caller of the library on the shipped default.
_N_MAX_TOLERANCE = 65536

_NUMU_ROW = (3, 4, 5)
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
            # The tolerance dial used to stop at 1e-5 because rtol=1e-6 and
            # tighter RAISED -- but that was the 1024-slab ceiling, which is
            # only the DEFAULT of `n_max`, not a limit of the method.  The
            # constant's own docstring says to pass `n_max` to move it, and
            # `_N_MAX_TOLERANCE` below does.  The two dials then trace one
            # curve: measured on CHORD/12x1, rtol 3e-5 and 512 slabs both
            # return 2.1875e-07, rtol 1e-5 and 1024 both 5.4681e-08, and
            # rtol 1e-7 reaches 7.8942e-10 where it used to raise.
            'knob_domain': [1, 2, 4, 8, 16, 32, 64, 128, 256,
                            -3, -4, -5, -6, -7, -8],
        }

    def environment(self):
        r"""Returns what this library actually ran with, for the artifact.

        ADVERSARIAL.md's sixth attack is that THIS code might be
        under-represented, and its finding was that nothing asserted the fast
        path at run time.  A missing numba, a flipped ``USE_NUMBA`` or a
        moved batching threshold would have sent the work down the NumPy
        path and published that as ours, with no error and no hint in the
        artifact.  So the fast path is asserted here rather than assumed, and
        what was actually entered is recorded.
        """
        import fastkernels
        detail = {
            'have_numba': bool(fastkernels.HAVE_NUMBA),
            'use_numba': bool(getattr(fastkernels, 'USE_NUMBA', False)),
            'available': bool(fastkernels.available()),
            'kernel_entered': self._kernel,
            'batched_energy_and_zenith': self._costhz is not None,
        }
        try:
            import numba
            detail['numba_threads'] = int(numba.get_num_threads())
        except Exception:                                      # noqa: BLE001
            detail['numba_threads'] = None
        return detail

    def _assert_fast_path(self):
        r"""Refuses to measure this library on its slow path."""
        import fastkernels
        if not (fastkernels.HAVE_NUMBA and fastkernels.available()
                and getattr(fastkernels, 'USE_NUMBA', False)):
            raise SystemExit(
                'NuOscProbExact: the compiled kernels are not live '
                '(HAVE_NUMBA=%r, USE_NUMBA=%r, available=%r).  Measuring the '
                'NumPy path and publishing it as this library would '
                'under-represent it; install the fast extra and re-run.'
                % (fastkernels.HAVE_NUMBA,
                   getattr(fastkernels, 'USE_NUMBA', None),
                   fastkernels.available()))

    def setup(self, problem):
        p = problem
        # Which parameter the scan turns.  Objection Earth-2 named Dmsq31
        # as a realistic thing for a fit to move; a code that caches does
        # not invalidate the same work for it as for delta_CP.
        self._scan_dmsq31 = getattr(problem, 'scan', 'dcp') == 'dmsq31'

        self._assert_fast_path()
        # Which compiled kernel this problem reaches: the chord path for an
        # Earth crossing, the constant-density path otherwise.  Recorded in
        # every artifact so a silent fallback is visible afterwards.
        self._kernel = ('earth_chords_3nu_kernel' if p.costhz
                        else 'probabilities_3nu_kernel')
        self._e_ev = np.asarray(p.energies_gev, dtype=float)*1.0e9
        self._costhz = (np.asarray(p.costhz, dtype=float)
                        if p.costhz else None)
        self._ye = p.ye
        self._sines = (math.sqrt(p.s12sq), math.sqrt(p.s23sq),
                       math.sqrt(p.s13sq))
        self._dm21, self._dm31 = p.dm21, p.dm31
        self._dcp = p.dcp

        knob = int(p.knob)
        override = getattr(p, 'rtol_override', 0.0)
        if override:
            # An explicit tolerance, so the dial is not confined to powers
            # of ten.  n_slabs stays at the library default and the search
            # starts from it, exactly as for a knob-derived rtol.
            self._n_slabs, self._rtol = 1, float(override)
        elif knob > 0:
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

    def configure(self, v):
        s12, s23, s13 = self._sines
        dcp, dm31 = ((self._dcp, v) if self._scan_dmsq31
                     else (v, self._dm31))
        self._h_vac = np.asarray(
            hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
                s12, s23, s13, dcp, self._dm21, dm31))

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

    def probabilities(self):
        r"""Untimed.  The scored channel over the whole grid, in grid order.

        Grid order is zenith outer, energy inner, which is what the
        broadcast in ``_probs`` already produces: indexing energies on the
        last axis and angles on the one before it makes ``ravel`` walk
        exactly that way.
        """
        return [float(v) for v in self._probs().reshape(-1)]

    def _probs(self):
        r"""The nu_mu row for the whole grid: (..., 3), grid order then channel."""
        if self._costhz is not None:
            if self._costhz.size == 1:
                energy, costhz = self._e_ev, float(self._costhz[0])
            else:
                energy, costhz = self._e_ev[None, :], self._costhz[:, None]
            probs = earth.probabilities_3nu_earth(
                self._h_vac, energy, costhz,
                n_slabs_per_segment=self._n_slabs,
                electron_fraction=self._ye, rtol=self._rtol,
                n_max=_N_MAX_TOLERANCE)
            return np.asarray(probs)[..., _NUMU_ROW]
        h = np.asarray(hamiltonians3nu.hamiltonian_3nu_matter(
            self._h_vac, self._e_ev, self._vcc))
        return np.asarray(
            oscprob3nu.probabilities_3nu(h, self._L_inv_ev))[..., _NUMU_ROW]

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
                electron_fraction=self._ye, rtol=self._rtol,
                n_max=_N_MAX_TOLERANCE)
            return float(np.sum(np.asarray(probs)[..., _PMM]))
        h = np.asarray(hamiltonians3nu.hamiltonian_3nu_matter(
            self._h_vac, self._e_ev, self._vcc))
        probs = oscprob3nu.probabilities_3nu(h, self._L_inv_ev)
        return float(np.sum(np.asarray(probs)[..., _PMM]))


ADAPTER = NuOscProbExact
