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

Conventions.  Nothing is rescaled.  nuSQuIDS computes electron density as
``rho * N_A * Y_e`` with its own Avogadro constant where this library
carries ``Y_e rho / m_bar``, and its km differs from this library's in the
fifth decimal of hbar c.  Both are absorbed into nuSQuIDS' OWN 50-digit
reference -- ``conversions.km_to_inv_ev('nuSQuIDS')`` reads the second from
``Const()`` -- so what the solver is handed here is the honest problem and
its residual measures its algorithm.

The knob is the solver tolerance: ``rel_error = abs_error = 10**-knob``
(manifest domain 1e-3..1e-12, so knob in 3..12).  ``knob <= 0`` does NOT
leave the solver at its own defaults -- those fail on both Earth grids --
but sets an explicit 1e-7; see ``setup()``.
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
        self._constant_density_mode = False
        units = nsq.Const()
        e = np.asarray(p.energies_gev, dtype=float)*units.GeV

        # No density rescaling and no cosine rescaling.  nuSQuIDS' own
        # Avogadro constant and its own km are absorbed into ITS reference,
        # so what it is handed here is the honest physical problem.

        if p.costhz:
            r = _aligned_radii(200)
            rho = earth.density_prem(np.clip(r, 0.0, gd.EARTH_RADIUS))
            body = nsq.EarthAtm((r/gd.EARTH_RADIUS).tolist(),
                                rho.tolist(),
                                np.full(r.size, p.ye).tolist())
            body.SetAtmosphereHeight(0.0)
            self._tracks = [body.MakeTrackWithCosine(cz) for cz in p.costhz]
        else:
            body = nsq.ConstantDensity(p.density, p.ye)
            self._tracks = [nsq.ConstantDensity.Track(p.L_km*units.km)]
            # Constant density has a closed form and nuSQuIDS ships it:
            # with this set, EvolveState skips the ODE entirely and evolves
            # each node algebraically.  A seasoned user of this code would
            # never integrate an ODE through constant density, and running
            # it that way is the same omission objection LBL-3 was about,
            # aimed at a different code.  Verified against the ODE path on
            # CONST/60E: they agree to 7.2e-7, which is the ODE's own
            # tolerance.
            self._constant_density_mode = True

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
        # ALWAYS set the tolerance.  Leaving nuSQuIDS at its constructor
        # defaults is not a runnable configuration for this comparison: with
        # them the GSL solver fails outright ("Error in GSL ODE solver") on
        # BOTH Earth grids in the manifest -- at cosz = -1.0, the first point
        # of OSC/100x100, and at cosz = -0.9, the whole of CHORD/12x1.  Every
        # explicit tolerance tried, 1e-6 through 1e-10, succeeds.  That went
        # unnoticed because the cosine used to be rescaled by 1 - 1.4e-7 to
        # absorb this code's hbar c, which nudged the angle off the failing
        # value; removing that rescaling for the per-code references exposed
        # it.  knob = 0 therefore means an explicit mid-sweep tolerance, not
        # "whatever the constructor picked".
        tolerance = 10.0**-p.knob if p.knob > 0 else 1.0e-7
        s.Set_rel_error(tolerance)
        s.Set_abs_error(tolerance)
        self._tolerance = tolerance
        if getattr(self, '_constant_density_mode', False):
            s.Set_AllowConstantDensityOscillationOnlyEvolution(True)
        s.Set_Body(body)
        self._s = s

    def configure(self, dcp):
        self._s.Set_CPPhase(0, 2, dcp)

    def reset(self):
        r"""Nothing to reset: ``evaluate()`` calls ``Set_initial_state`` and
        ``EvolveState`` for every track, so each repetition re-solves from
        scratch and is already cold."""

    def evaluate(self):
        sink = 0.0
        for track in self._tracks:
            self._s.Set_Track(track)
            self._s.Set_initial_state(self._state, nsq.Basis.flavor)
            self._s.EvolveState()
            if self._single:
                sink += self._s.EvalFlavor(1)      # P(nu_mu -> nu_mu)
            else:
                # AtNode, not EvalFlavor(flavor, E): the grid energies ARE
                # this solver's nodes, so the interpolating overload would
                # charge nuSQuIDS an interpolation inside the timed region
                # and carry its interpolation error into the accuracy
                # number -- both against it, for nothing.
                sink += sum(self._s.EvalFlavorAtNode(1, i)
                            for i in range(self._e.size))
        return float(sink)


ADAPTER = NuSQuIDS
