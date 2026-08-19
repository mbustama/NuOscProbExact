# -*- coding: utf-8 -*-
r"""nuSQuIDS from the pinned wheel, in its multiple-energy mode.

Everything invariant under the delta_CP scan is built in ``setup()``: the
body (the PREM table for chord problems, a constant-density slab
otherwise), the track(s), the solver with its energy vector, every mixing
angle and mass splitting, and the solver tolerance.  ``configure()`` moves
delta_CP alone, through ``Set_CPPhase``.  ``evaluate()`` resets the initial
state and evolves the WHOLE energy stack in one solve -- the multiple-energy
constructor is the batched entry point, and driving it one energy at a time
would rebuild the solver's advantage away.  Zenith is not batched
(``batches_zenith = False``): with more than one zenith the tracks are
looped inside ``evaluate()``.

Conversions, all derived rather than typed:

* Density.  nuSQuIDS takes rho in g/cm^3 and computes electron density as
  ``rho * N_A * Y_e`` with ITS OWN Avogadro constant (exposed as
  ``nsq.Const().Na``); this library carries ``Y_e rho / m_bar``.  The handed
  density is scaled by the ratio -- the usual nuclear mass defect, here
  0.99209 -- computed from ``globaldefs`` and nuSQuIDS' own constant.
* Length.  nuSQuIDS' km in eV^-1 differs from this library's in the fifth
  decimal of hbar c; the chord is L = -2 R cos th_z, so the cosine is
  shrunk by the ratio, derived at run time from ``Const`` and
  ``globaldefs``.

The knob is the solver tolerance: ``rel_error = abs_error = 10**-knob``
(manifest domain 1e-3..1e-12, so knob in 3..12); ``knob <= 0`` leaves the
solver at its own defaults.
"""

import os
import sys

_REPO = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(_REPO, 'src'))

import numpy as np                                             # noqa: E402

import globaldefs as gd                                        # noqa: E402
import earth                                                   # noqa: E402
import nuSQuIDS as nsq                                         # noqa: E402


def _aligned_radii(n_per_shell):
    r"""Grid points that land exactly on every PREM shell boundary.

    The same table construction as tests/prem_scan.py: piecewise-linear
    interpolation nodes aligned with the discontinuities, holding each
    jump with a 1e-9 km offset.
    """
    edges = np.concatenate(([0.0], np.asarray(earth.PREM_BOUNDARIES,
                                              dtype=float)))
    edges = edges[edges <= gd.EARTH_RADIUS]
    if edges[-1] < gd.EARTH_RADIUS:
        edges = np.concatenate((edges, [gd.EARTH_RADIUS]))
    eps = 1.0e-9                              # km, to carry the jump
    pieces = []
    for a, b in zip(edges[:-1], edges[1:]):
        inner = np.linspace(a, b, n_per_shell)
        inner[0], inner[-1] = a + eps, b - eps
        pieces.append(inner)
    r = np.concatenate(([0.0], np.concatenate(pieces), [gd.EARTH_RADIUS]))
    return np.unique(r)


class NuSQuIDS(object):

    name = 'nuSQuIDS'

    def capabilities(self):
        return {
            'batches_energy': True,
            'batches_zenith': False,
            'batch_symbol': ('nuSQUIDS(E vector, ...) multiple-energy '
                             'constructor + EvolveState'),
            'knob_name': 'log10_solver_tolerance',
            'knob_domain': [3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        }

    def setup(self, problem):
        p = problem
        units = nsq.Const()
        e = np.asarray(p.energies_gev, dtype=float)*units.GeV

        # ours/theirs: electrons per gram from globaldefs against nuSQuIDS'
        # own Avogadro constant -- the nuclear mass defect, derived.
        electrons_per_gram = gd.CONV_G_TO_EV/((gd.MASS_PROTON
                                               + gd.MASS_NEUTRON)/2.0)
        density_scale = electrons_per_gram/units.Na

        if p.costhz:
            r = _aligned_radii(200)
            rho = earth.density_prem(np.clip(r, 0.0, gd.EARTH_RADIUS))
            body = nsq.EarthAtm((r/gd.EARTH_RADIUS).tolist(),
                                (rho*density_scale).tolist(),
                                np.full(r.size, p.ye).tolist())
            body.SetAtmosphereHeight(0.0)
            # nuSQuIDS' km against ours, absorbed into the chord's cosine.
            czscale = 1.0/((units.km*units.eV)/gd.CONV_KM_TO_INV_EV)
            self._tracks = [body.MakeTrackWithCosine(cz*czscale)
                            for cz in p.costhz]
        else:
            body = nsq.ConstantDensity(p.density*density_scale, p.ye)
            self._tracks = [nsq.ConstantDensity.Track(p.L_km*units.km)]

        self._single = e.size == 1
        if self._single:
            s = nsq.nuSQUIDS(3, nsq.NeutrinoType.neutrino)
            s.Set_E(float(e[0]))
            self._state = np.array([0.0, 1.0, 0.0])
        else:
            s = nsq.nuSQUIDS(e, 3, nsq.NeutrinoType.neutrino, False)
            self._state = np.zeros((e.size, 3))
            self._state[:, 1] = 1.0
        self._e = e

        s.Set_MixingAngle(0, 1, np.arcsin(np.sqrt(p.s12sq)))
        s.Set_MixingAngle(0, 2, np.arcsin(np.sqrt(p.s13sq)))
        s.Set_MixingAngle(1, 2, np.arcsin(np.sqrt(p.s23sq)))
        s.Set_SquareMassDifference(1, p.dm21)
        s.Set_SquareMassDifference(2, p.dm31)
        s.Set_CPPhase(0, 2, p.dcp)
        if p.knob > 0:
            s.Set_rel_error(10.0**-p.knob)
            s.Set_abs_error(10.0**-p.knob)
        s.Set_Body(body)
        self._s = s

    def configure(self, dcp):
        self._s.Set_CPPhase(0, 2, dcp)

    def evaluate(self):
        sink = 0.0
        for track in self._tracks:
            self._s.Set_Track(track)
            self._s.Set_initial_state(self._state, nsq.Basis.flavor)
            self._s.EvolveState()
            if self._single:
                sink += self._s.EvalFlavor(1)      # P(nu_mu -> nu_mu)
            else:
                sink += sum(self._s.EvalFlavor(1, float(ee))
                            for ee in self._e)
        return float(sink)


ADAPTER = NuSQuIDS
