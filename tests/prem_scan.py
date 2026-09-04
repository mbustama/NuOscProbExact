# -*- coding: utf-8 -*-
r"""Regenerates the frozen PREM speed-accuracy data used by the paper.

Companion to :mod:`tests.nusquids_scan`.  That one freezes a curve at
constant density; this one freezes the two *Earth* speed-accuracy planes,
at three flavors and at 3+1, so that the notebook can draw them without
building or importing any external code.

Run it, from the repository root, in an environment that has nuSQuIDS and
nuCraft on the path::

    python tests/prem_scan.py > tests/prem_speed_accuracy.json

It takes about seven minutes, nearly all of it in the arbitrary-precision
referee.

How this library is called
--------------------------
Since 1.12.0 the Earth routines take an array of energies, so the twelve
points are one batched call and not twelve scalar ones.  That is how the
question would actually be asked, and it is what makes the comparison
like for like: NuFast-Earth is handed the same twelve energies through
``Set_Spectra``, and the GLoBES and Prob3++ drivers loop over them inside
C.  Batching alone is worth between 4x and 10x here, on top of the
compiled evolution-operator kernel that 1.13.0 gave the slab path.

At four flavors :data:`oscprob4nu.ROOT_STRATEGY` selects how the latent
roots are found, and both routes are measured and frozen here even though
the figure now draws only the default.  On this problem they are not
an accuracy trade: they agree **bit for bit** at every slab count, because
the discretisation error swamps a root difference of 3.6e-17 against
3.9e-16.  The two curves therefore differ horizontally and not vertically,
and what the panel shows is the cost of the choice, not its consequence.

What a PREM figure can and cannot measure
-----------------------------------------
The constant-density figure measures how exactly a code evaluates the
probability for a *given* Hamiltonian, and bottoms out at round-off.  A
PREM figure cannot: the density varies continuously and every code,
including this one, approximates the profile.  So the y-axis here is
error against a *converged* PREM solution, and each code's position is
set by how finely it slices or integrates the profile rather than by the
accuracy of its probability formula.  This library therefore appears as
a curve, with ``n_slabs_per_segment`` as its dial, not as a horizontal
line.

The referee
-----------
The slab product converges at **second** order in the slab width -- the
density is sampled at slab midpoints, and the PREM shell boundaries are
cut exactly, so the leading error is O(h^2), not O(h).  Measured
convergence exponents are 2.000 at 3, 10 and 40 GeV.  The referee is
therefore the second-order Richardson extrapolation

    ref = (4*P(256) - P(128))/3

of a 30-digit ``mpmath`` slab product, *not* the first-order
``2*P(32) - P(16)``, which for this problem amplifies the error rather
than cancelling it.

Because that referee shares :func:`earth.earth_slabs` with the code it
judges, it is itself refereed: an adaptive DOP853 integration of
``dpsi/dx = -i H(x) psi`` through the continuous profile, with the shell
boundaries handed to the integrator as segment ends, agrees with it to
between 6.5e-13 and 2.3e-11 over the energies used here.  Two
discretisation families with nothing in common but the PREM definition
itself.

Conventions matched, and how each was established
-------------------------------------------------
Skip any of these and the comparison measures bookkeeping, not physics.

* **Atmosphere.**  ``nsq.EarthAtm``'s production height defaults to 22 km
  and nuCraft's to 20 km, which lengthen the chord at
  cos th_z = -0.9 by 24.4 km and by a comparable amount.  Both are set to
  zero, after which nuSQuIDS' baseline reproduces
  ``earth.distance_traveled_inside_earth(-0.9)`` = 11467.8 km exactly.

* **Length unit.**  nuSQuIDS' km is 5.0677307162e9 eV^-1 against this
  library's 5.06773e9, a relative excess of 1.413e-7.  An earlier revision
  absorbed that by shrinking the cos th_z it was handed, and that trick is
  wrong on a chord: the cosine sets the geometry, so the turning point
  moved 1.7 km shallower and every PREM crossing moved with it, and
  nuSQuIDS integrated the matter term through a measurably different
  Earth -- 5.1e-7 of probability at these energies, half of a 1.03e-6
  floor that was blamed on the discontinuities.  A km mismatch is a
  phase-per-length mismatch, so it is absorbed in the phases instead:
  every Delta m^2 nuSQuIDS receives is scaled by 1/1.000000141317, the
  density scale carries the same power for the matter term, and the
  cosine goes in raw.  In vacuum over this chord the two codes then agree
  to **1.3e-15**; with no correction at all they differ by 1.5e-6.

  nuCraft's km is 5.067732e9 eV^-1 -- the ``-2.533866j`` hard-coded in
  its right-hand side is -2j*(km in 1/eV)/2e9 -- an excess of 3.947e-7,
  and it is treated identically: its Delta m^2 are scaled, its density
  scale carries the matching power, and the cosine goes in raw.  In
  vacuum at 3+1 it then agrees with this library to **3.8e-15**.

* **Profile.**  Both external codes ship their own PREM tabulation --
  nuSQuIDS a 201-point table with ye = 0.4656/0.4957, nuCraft a 58-point
  linear spline -- and using them would compare Earth models rather than
  solvers.  Both are handed this library's PREM instead, at ye = 0.5.

  For nuSQuIDS this uses the documented ``(x, rho, ye)`` constructor, on a
  grid whose points land exactly on the PREM shell boundaries.  That
  matters: a uniform grid can only place a density discontinuity to within
  half a spacing, and the chord crosses ten of them, which holds nuSQuIDS
  at 3.4e-6 no matter how dense the table.  Boundary-aligned, its error
  follows the table -- 2.7e-9 at 200 points per shell, 2.8e-10 at the
  12800 used here -- because ``EarthAtm`` fits an Akima spline through
  the points and the spline, not the solver, is what converges.  (An
  earlier revision saw no node sensitivity at all, 9.1e-7 from 200 points
  onward; that was the two convention errors above swamping it.)

  nuCraft has no such constructor, and its documented file route calls the
  Python-2 builtin ``file()``, so it cannot run at all on Python 3.  The
  profile is therefore installed by replacing ``EarthModel.profInt``.

* **Matter potential.**  nuCraft's V_CC is larger than this library's by a
  constant factor.  Its source comment says the normalisation is
  sqrt(2) G_F rho/m_N with m_N = 0.939 GeV, but the constant it actually
  carries, 15.256e-5, is the atomic-mass-unit value: this library's
  equivalent constant is 15.14423e-5, and 15.14423/15.256 = 0.9926737.
  One power of its km excess rides on top, exactly as for nuSQuIDS,
  giving 0.9926733.  An earlier revision stopped at the bare constant
  ratio, 3.9e-7 high, and that -- with the same power missing from the
  mass splittings -- held nuCraft's three-flavor curve at 2.5e-6 however
  small its numPrec; with both carried, the same run follows its dial
  down to 3.6e-10, and at constant density the two codes agree to
  1.7e-11.  (The residual scan that pinned the GLoBES factor gave
  0.9926748 here; like nuSQuIDS' unity scan, it sat beside standing
  convention errors, so this factor too is derived, not scanned.)
  nuSQuIDS takes rho N_A Y_e in its own constants, so its factor
  is derived from its own ``Const``: its G_F and N_A match `globaldefs`
  bit for bit, but three powers of its hbar c live in its cm^-3 and a
  fourth rides on the km bridge, giving 0.9920927 where the bare mass
  defect built from this library's masses is 0.9920924.  That 4.24e-7 gap
  cost 5.3e-7 of probability here, and the residual scan that pinned the
  other codes' factors did not catch it: scanned with the old cosine
  trick's geometry error standing beside it, the minimum sat at 1.0000000
  with the floor intact.  So this factor is derived, not scanned.

Per-code traps
--------------
* nuSQuIDS once appeared to floor at ~1e-6 here regardless of tolerance,
  stepper (RK4/RKF45/RKCK/RK8PD/MSADAMS all identical) or h_max, and the
  floor was blamed on the PREM discontinuities.  Both halves of it were
  this harness's own bookkeeping: 5.3e-7 from a density scale built out
  of this library's masses alone (see **Matter potential**) and 5.1e-7
  from the cosine shrink (see **Length unit**), summing to the measured
  1.0314e-6.  With both conventions carried exactly, the same solver
  follows its tolerance dial down to 2.8e-10 at three flavors and
  1.3e-10 at 3+1 on the 12800-point table.
* nuCraft's ``CalcWeights`` takes the **zenith angle in radians**, though
  ``ConstructMixingMatrix`` takes its angles in degrees.  Passing degrees
  to both leaves a disagreement of 0.198 that no convention matching
  removes.
* nuCraft's default ``rICore`` is 1121.5 km against the 1221.5 km in its
  own PREM table, and ``EarthModel.__init__`` writes ``self.rE == rE``
  where it means ``=``, so the ``rE`` keyword is silently ignored.  Both
  are worth reporting upstream.

* **nuCraft cannot get below 2.8e-3 at 3+1, and it is not the solver.**
  Its sterile entry and its charged-current entry are set by two
  independently rounded constants, ``7.6525e-5*(1-y)`` and
  ``15.256e-5*y``, whose ratio is 0.5016 where the isoscalar value
  ``(1/2)(1-Y_e)/Y_e`` is exactly 0.5 -- 0.32% high.  Scaling the density
  cannot separate them, because it scales both, so no setting removes that
  floor.  Forcing the ratio to 0.5 by hand drops the same run to 2.8e-7,
  which is what its solver is actually worth.  The curve is drawn as
  released, floor and all, because that is what a user of nuCraft gets;
  the paper says where the floor comes from rather than leaving a reader
  to conclude the solver is poor.  Worth reporting upstream with the
  other two.

The three compiled codes
------------------------
NuFast-Earth, GLoBES and Prob3++ are C or C++ and cannot be driven from
here.  Their drivers, and the raw output each produced, are in
``tests/external_drivers``; this module only reads the ``.txt`` files, so
none of the three has to be installed to regenerate the figure.  The
README beside them carries the full account.  In brief: all three are
handed this library's PREM at Y_e = 0.5, on the four-major-layer,
n-sub-shell scheme each uses natively, from a header generated out of
``earth._PREM_COEFFS`` rather than transcribed.  All three carry
hbar c = 1.97327e-7 eV m, so each gets a cosine shrunk by 4.23e-8, which
takes them from ~1e-7 to 1.2e-14, 1.8e-14 and (for Prob3++) its own
floor in vacuum.  And all three carry a *different* rounding of the same
atomic-mass-unit constant -- 1.526493231029146e-4, 7.63247e-14 and
1.52588e-4 -- so each needs its own factor, 0.9920938, 0.992093 and
0.9924922; scanning returns 1.000000 for the first two.

All three are three-flavor only, so none of them appears at 3+1.
"""

import json
import os
import sys
import time
import warnings

import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp

sys.path.insert(0, 'src')
sys.path.insert(0, '/home/mbustamante/Research/AtmNu/nuCraft/trunk')

import globaldefs as gd                                       # noqa: E402
import earth                                                  # noqa: E402
import hamiltonians3nu                                        # noqa: E402
import hamiltonians4nu                                        # noqa: E402
import oscprob4nu                                             # noqa: E402

import nuSQuIDS as nsq                                        # noqa: E402

# Importing nuCraft installs a `warnings.showwarning` replacement with the
# Python 2 signature, which raises TypeError the moment nuCraft actually
# warns -- as it does when its unitarity check misses a tight `numPrec`.
# Keep the standard one so that the warning prints instead of killing the
# run.  This is the third Python 2 remnant in nuCraft r22, after the
# `file()` call in EarthModel and the `self.rE == rE` no-op.
_SHOWWARNING = warnings.showwarning
from NuCraft import NuCraft, EarthModel                       # noqa: E402
warnings.showwarning = _SHOWWARNING

UNITS = nsq.Const()

COSTHZ = -0.9
YE = 0.5
DPS = 30

# NuFit 4.0, normal ordering, as `globaldefs` carries it.
TH12 = np.arcsin(np.sqrt(0.310))
TH13 = np.arcsin(np.sqrt(2.240e-2))
TH23 = np.arcsin(np.sqrt(0.582))
DCP = 217.0/180.0*np.pi
DM21 = 7.39e-5
DM31 = 2.525e-3

# The 3+1 point of the paper's sterile oscillogram.
DM41 = 1.0
TH14 = np.arcsin(np.sqrt(0.10))
TH24 = np.arcsin(np.sqrt(0.10))
TH34 = 0.0

# The tolerance dial.  `_n_for_tolerance` starts its search at
# `n_slabs_per_segment` and only ever doubles upward, so leaving that at its
# default of 8 pinned the coarse end of this curve at the accuracy of eight
# slabs however loose the tolerance: the curve was a truncated copy of the
# slab sweep rather than an independent span of it.  Starting the search at
# 1 makes the two coextensive, and these tolerances then select 1, 2, 8, 16,
# 32, 64, 256, 512 and 1024 slabs per segment.
#
# The dial used to stop at 1e-5, and the reason given was that tighter
# tolerances "cannot be met within `slabs.N_SLABS_MAX`".  That ceiling is a
# guardrail, not a limit: it is the DEFAULT of the `n_max` argument, and the
# constant's own docstring says to pass `n_max` to move it.  Raising it here
# rather than shipping a larger default leaves every other caller alone.
# The tolerance dial and the slab dial are then the same curve -- measured,
# rtol 3e-5 and 512 slabs return an identical 2.1875e-07, as do 1e-5 and
# 1024 -- so the two are coextensive rather than one being a truncation of
# the other.
N_MAX_TOLERANCE = 65536
RTOLS = (3e0, 1e0, 3e-1, 3e-2, 1e-2, 1e-3, 1e-4, 3e-5, 1e-5,
         1e-6, 1e-7, 1e-8)

E_GEV_3NU = np.logspace(np.log10(3.0), np.log10(40.0), 12)
E_GEV_4NU = np.logspace(np.log10(300.0), np.log10(30000.0), 12)

# This library's normalisation constant, A = k rho Y_e / 2e9, derived from
# its own matter potential rather than transcribed.
_OUR_A_CONST = 2.0*1.0e9*earth.matter_potential(1.0, 0.5)/0.5

# nuSQuIDS' km in eV^-1 against ours: one power of hbar c, 1 + 1.413e-7.
_KM_NUSQUIDS_OVER_OURS = (UNITS.km*UNITS.eV)/gd.CONV_KM_TO_INV_EV

# The factor nuSQuIDS' density must carry so that the potential it builds,
# sqrt(2) G_F N_A rho Y_e / cm^3 in ITS constants, integrates to this
# library's V_CC phase over ITS km.  Its G_F and N_A are bit-identical to
# `globaldefs`', but its cm^-3 carries THREE powers of its hbar c, and a
# fourth power rides on the Delta m^2 bridge below.  Deriving this factor
# from this library's masses alone -- the bare nuclear mass defect,
# 0.9920924 -- misses those powers: it is 4.24e-7 low, which cost 5.3e-7
# of probability and was half of the 1.03e-6 floor this figure carried
# for a while.
DENSITY_SCALE_NUSQUIDS = (
    _OUR_A_CONST
    / (np.sqrt(2.0)*UNITS.GF*UNITS.Na*UNITS.cm**-3.0*2.0e9)
    / _KM_NUSQUIDS_OVER_OURS)

# The other half of that floor was the old length matching, which shrank
# the cosine handed to nuSQuIDS by the km ratio.  On a chord the cosine is
# geometry, not bookkeeping: shrinking it moved the turning point 1.7 km
# shallower and every PREM crossing with it, so nuSQuIDS propagated
# through a measurably different Earth -- 5.1e-7 of probability at these
# energies.  A km mismatch is a phase-per-length mismatch, so it is
# absorbed where it lives instead: every Delta m^2 handed to nuSQuIDS is
# scaled by this factor, the density scale above carries the same power
# for the matter term, and the cosine goes in raw.  In vacuum over this
# chord the two codes then agree to 1.3e-15.
DM_SCALE_NUSQUIDS = 1.0/_KM_NUSQUIDS_OVER_OURS

# nuCraft's km in eV^-1 against ours: the -2.533866j hard-coded in its
# right-hand side is -2j*(km in 1/eV)/2e9, so its km is 5.067732e9 against
# this library's 5.06773e9, an excess of 3.947e-7.  Same treatment as
# nuSQuIDS: every Delta m^2 it receives carries this factor, the density
# scale below carries the same power for the matter term, and the cosine
# goes in raw.
_KM_NUCRAFT_OVER_OURS = 2.533866*2.0e9/gd.CONV_KM_TO_INV_EV
DM_SCALE_NUCRAFT = 1.0/_KM_NUCRAFT_OVER_OURS

# The factor nuCraft's density must carry so that the CC potential it
# builds from its own constant, 15.256e-5 (NuCraft.py:212), integrates to
# this library's V_CC phase over ITS km.  Stopping at the bare constant
# ratio, 0.9926737, misses the km power -- the same class of error that
# cost nuSQuIDS 5.3e-7 -- and held nuCraft's three-flavor curve at 2.5e-6
# however small its numPrec.
DENSITY_SCALE_NUCRAFT = _OUR_A_CONST/15.256e-5*DM_SCALE_NUCRAFT

# nuCraft's sterile matter entry against its charged-current one, both
# from NuCraft.py:212: 7.6525e-5/15.256e-5 = 0.5016059, where its own
# formula -- one constant times (2*Y_e, 0, 0, 1-Y_e) -- makes the
# isoscalar value exactly 1/2.  A dimensionless internal inconsistency,
# 0.32% high, not a convention.  Its 3+1 referee is built with THIS ratio,
# so the plotted error measures its solver; see the module docstring.
NUCRAFT_NC_OVER_CC = 7.6525e-5/15.256e-5

H_VAC_3NU = np.asarray(
    hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF,
        gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF))
H_VAC_4NU = np.asarray(
    hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(
        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF,
        np.sin(TH14), np.sin(TH24), np.sin(TH34),
        gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF, DM41))


# --------------------------------------------------------------------------
# The referee
# --------------------------------------------------------------------------

def _slab_hamiltonians(energy_ev, n, n_flavors, nc_over_cc=None):
    """Per-slab Hamiltonians and widths, exactly as `earth` builds them.

    `nc_over_cc` overrides the sterile-to-CC potential ratio |V_NC|/V_CC;
    the default None takes `earth.matter_potential_nc`, which at
    Y_e = 1/2 is exactly one half.  nuCraft's referee passes its own
    0.5016059.  At three flavors the NC term is common to all active
    flavors and drops out, so the parameter has nothing to act on.
    """
    widths_km, densities = earth.earth_slabs(COSTHZ, n_slabs_per_segment=n)
    vcc = earth.matter_potential(densities, YE)
    if n_flavors == 3:
        h = hamiltonians3nu.hamiltonian_3nu_matter(H_VAC_3NU, energy_ev, vcc)
    else:
        if nc_over_cc is None:
            vnc = earth.matter_potential_nc(densities, electron_fraction=YE)
        else:
            vnc = -nc_over_cc*vcc
        h = hamiltonians4nu.hamiltonian_4nu_matter(
            H_VAC_4NU, energy_ev, vcc, vnc)
    return np.asarray(h), widths_km


def mp_probability(energy_ev, n, n_flavors, dps=DPS, nc_over_cc=None):
    r"""P(nu_mu -> nu_mu) from a `dps`-digit slab product."""
    mp.mp.dps = dps
    h, widths_km = _slab_hamiltonians(energy_ev, n, n_flavors, nc_over_cc)
    u = mp.eye(n_flavors)
    for hk, wk in zip(h, widths_km):
        m = mp.matrix(n_flavors, n_flavors)
        for i in range(n_flavors):
            for j in range(n_flavors):
                m[i, j] = mp.mpc(complex(hk[i, j]))
        # U = U_n ... U_1: the first slab crossed is the rightmost factor.
        u = mp.expm(-1j*m*(mp.mpf(float(wk))*mp.mpf(gd.CONV_KM_TO_INV_EV)))*u
    return float(abs(u[1, 1])**2.0)


#: Slab counts the Richardson referee is built from, per flavor count.
#:
#: Three flavors converge fast enough that 128 and 256 leave a referee that
#: agrees with the independent ODE integration to 2.3e-11, some three
#: decades below the smallest error plotted against it.
#:
#: 3+1 does not.  A 1 eV^2 splitting oscillates far faster along the same
#: chord, so the same pair leaves 1.2e-9 -- and every code in that panel now
#: reaches 1.2e-9 to 1.8e-9, which is to say they were all measuring the
#: referee rather than themselves.  Measured at 300 GeV against the ODE
#: integration: Richardson(128,256) 1.169e-09, Richardson(256,512)
#: 6.674e-11, Richardson(512,1024) 1.998e-12.  The last is three decades
#: clear of the codes and costs about 68 s per energy instead of 18 s.
_REFEREE_SLABS = {3: (128, 256), 4: (512, 1024)}


def reference(energy_ev, n_flavors, nc_over_cc=None):
    """Second-order Richardson extrapolation of the slab product."""
    n_coarse, n_fine = _REFEREE_SLABS[n_flavors]
    fine = mp_probability(energy_ev, n_fine, n_flavors, nc_over_cc=nc_over_cc)
    coarse = mp_probability(energy_ev, n_coarse, n_flavors,
                            nc_over_cc=nc_over_cc)
    return (4.0*fine - coarse)/3.0


def ode_reference(energy_ev, n_flavors, tol=1.0e-13, nc_over_cc=None):
    """The same limit, from a discretisation with nothing in common.

    Integrates dpsi/dx = -i H(x) psi through the *continuous* profile,
    restarting at every PREM shell boundary so that no step straddles a
    density discontinuity.
    """
    total_km = earth.distance_traveled_inside_earth(COSTHZ)
    edges = earth.prem_layer_edges_along_chord(COSTHZ)

    def rhs(s_km, y):
        r = earth.earth_radial_distance_from_depth(COSTHZ, s_km)
        rho = earth.density_prem(r)
        vcc = earth.matter_potential(rho, YE)
        if n_flavors == 3:
            h = hamiltonians3nu.hamiltonian_3nu_matter(
                H_VAC_3NU, energy_ev, vcc)
        else:
            if nc_over_cc is None:
                vnc = earth.matter_potential_nc(rho, electron_fraction=YE)
            else:
                vnc = -nc_over_cc*vcc
            h = hamiltonians4nu.hamiltonian_4nu_matter(
                H_VAC_4NU, energy_ev, vcc, vnc)
        psi = y[:n_flavors] + 1j*y[n_flavors:]
        d = -1j*np.asarray(h).dot(psi)*gd.CONV_KM_TO_INV_EV
        return np.concatenate((d.real, d.imag))

    psi = np.zeros(n_flavors, dtype=complex)
    psi[1] = 1.0
    bounds = np.concatenate(([0.0], edges, [total_km]))
    for a, b in zip(bounds[:-1], bounds[1:]):
        if b <= a:
            continue
        sol = solve_ivp(rhs, (a, b),
                        np.concatenate((psi.real, psi.imag)),
                        method='DOP853', rtol=tol, atol=tol)
        y = sol.y[:, -1]
        psi = y[:n_flavors] + 1j*y[n_flavors:]
    return float(abs(psi[1])**2.0)


# --------------------------------------------------------------------------
# The codes
# --------------------------------------------------------------------------

def ours(energies_ev, n, n_flavors):
    r"""P(nu_mu -> nu_mu) for the whole energy stack, in one call.

    Since 1.12.0 the Earth routines take an array of energies, so the
    twelve points are one batched call rather than twelve scalar ones.
    That is how a user would ask this question, and it is the like-for-like
    comparison: NuFast-Earth is handed the same twelve through
    ``Set_Spectra``, and the GLoBES and Prob3++ drivers loop over them
    inside C.  Batching is worth between 4x and 10x here on its own.
    """
    if n_flavors == 3:
        return np.asarray(earth.probabilities_3nu_earth(
            H_VAC_3NU, energies_ev, COSTHZ, n_slabs_per_segment=n,
            electron_fraction=YE))[..., 4].ravel()
    return np.asarray(earth.probabilities_4nu_earth(
        H_VAC_4NU, energies_ev, COSTHZ, n_slabs_per_segment=n,
        electron_fraction=YE))[..., 5].ravel()


def _aligned_radii(n_per_shell):
    """Grid points that land exactly on every PREM shell boundary."""
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


def nusquids_body(n_per_shell=12800):
    # 12800 for the same reason as in tests/bench/adapters/nusquids.py: at
    # 200 points per shell the Akima spline EarthAtm fits through the table,
    # not the solver, sets the floor -- 2.7e-9 here against 2.8e-10 at
    # 12800.  (An earlier revision measured no node sensitivity at all and
    # kept 200; that was the two convention errors described above
    # swamping it.)
    r = _aligned_radii(n_per_shell)
    rho = earth.density_prem(np.clip(r, 0.0, gd.EARTH_RADIUS))
    body = nsq.EarthAtm((r/gd.EARTH_RADIUS).tolist(),
                        (rho*DENSITY_SCALE_NUSQUIDS).tolist(),
                        np.full(r.size, YE).tolist())
    body.SetAtmosphereHeight(0.0)
    return body


def nusquids(energies_ev, body, tol, n_flavors):
    r"""P(nu_mu -> nu_mu) for the whole energy stack, in one solver.

    nuSQuIDS has a multiple-energy constructor, and it is used here for the
    same reason this library's batched call is: driving it one energy at a
    time, rebuilding the solver for each, costs it a factor of 1.4 to 1.8
    that it does not have to pay.  GLoBES and Prob3++ have no equivalent
    and are looped, which is not a handicap imposed on them but the only
    interface they expose.

    The multiple-energy mode shares its adaptive stepping across the stack,
    so it does not return quite what the one-at-a-time mode does -- 2.5e-2
    apart at tolerance 1e-3, 1.8e-8 at 1e-10, both measured under the old
    conventions.  The batched mode is the one scored here.
    """
    e = np.asarray(energies_ev, dtype=float)/1.0e9*UNITS.GeV
    s = nsq.nuSQUIDS(e, n_flavors, nsq.NeutrinoType.neutrino, False)
    s.Set_MixingAngle(0, 1, TH12)
    s.Set_MixingAngle(0, 2, TH13)
    s.Set_MixingAngle(1, 2, TH23)
    # Every Delta m^2 carries the km bridge, so that phase per handed km
    # equals phase per this library's km; the cosine goes in raw.
    s.Set_SquareMassDifference(1, DM21*DM_SCALE_NUSQUIDS)
    s.Set_SquareMassDifference(2, DM31*DM_SCALE_NUSQUIDS)
    s.Set_CPPhase(0, 2, DCP)
    if n_flavors == 4:
        s.Set_MixingAngle(0, 3, TH14)
        s.Set_MixingAngle(1, 3, TH24)
        s.Set_MixingAngle(2, 3, TH34)
        s.Set_SquareMassDifference(3, DM41*DM_SCALE_NUSQUIDS)
    s.Set_rel_error(tol)
    s.Set_abs_error(tol)
    s.Set_Body(body)
    s.Set_Track(body.MakeTrackWithCosine(COSTHZ))
    state = np.zeros((e.size, n_flavors))
    state[:, 1] = 1.0
    s.Set_initial_state(state, nsq.Basis.flavor)
    s.EvolveState()
    return np.array([s.EvalFlavor(1, ee) for ee in e])


def nucraft_instance(n_flavors):
    em = EarthModel('prem', y=(YE, YE, YE), rICore=1221.5)
    # nuCraft's file route is Python 2 only, so the profile goes in here.
    em.profInt = lambda r: (
        earth.density_prem(min(float(r), gd.EARTH_RADIUS))
        * DENSITY_SCALE_NUCRAFT)
    # Every Delta m^2 carries nuCraft's km bridge, exactly as nuSQuIDS';
    # the cosine goes in raw so the geometry is untouched.
    masses = [1.0, DM21*DM_SCALE_NUCRAFT, DM31*DM_SCALE_NUCRAFT]
    angles = [(1, 2, np.degrees(TH12)), (1, 3, np.degrees(TH13), 217.0),
              (2, 3, np.degrees(TH23))]
    if n_flavors == 4:
        masses.append(DM41*DM_SCALE_NUCRAFT)
        angles += [(1, 4, np.degrees(TH14)), (2, 4, np.degrees(TH24)),
                   (3, 4, np.degrees(TH34))]
    return NuCraft(masses, angles, earthModel=em,
                   detectorDepth=0.0, atmHeight=0.0)


def nucraft(energies_ev, instance, prec):
    r"""The same, for nuCraft, which takes a list of particles.

    Handing it the whole list rather than one at a time is worth nothing --
    1.00x, and bit-identical answers -- because CalcWeights loops over the
    list internally.  It is done this way regardless, so that no code in
    these figures is looped from outside when its own interface would take
    the stack.
    """
    # atmMode 0 with atmHeight 0 is surface to surface; zenith in radians.
    zenith = np.arccos(COSTHZ)
    rows = [(14, float(e)/1.0e9, zenith) for e in np.atleast_1d(energies_ev)]
    return np.array([w[1] for w in instance.CalcWeights(rows, atmMode=0,
                                                        numPrec=prec)])


# --------------------------------------------------------------------------
# The sweep
# --------------------------------------------------------------------------

def timed(call, energies, repeat=5):
    """Microseconds per probability, best of `repeat` passes.

    `call` takes one energy; the whole stack is looped over here.  Used
    for nuSQuIDS and nuCraft, whose interfaces are one point at a time.
    """
    best = np.inf
    for _ in range(repeat):
        start = time.perf_counter()
        for energy_ev in energies:
            call(energy_ev)
        best = min(best, (time.perf_counter()-start)/len(energies)*1.0e6)
    return best*1.0


def timed_batch(call, n_energies, repeat=7, min_block=0.05):
    """The same, for a `call` that takes the whole stack at once.

    The first pass is thrown away: with Numba present it pays for
    compiling the kernel specialisation, which a user pays once per
    session and not once per call.

    Since 1.13.0 a batched call over twelve energies can finish in a few
    hundred microseconds, which is too short to time one at a time: the
    first attempt at this produced a curve whose cost fell between 64 and
    128 slabs per segment, which is not something the code does.  So the
    call is repeated inside the timed block until the block lasts at least
    `min_block` seconds, the way `timeit` autoranges, and the best of
    `repeat` such blocks is taken.

    That also matters for the two four-flavor root strategies, whose gap
    is smaller than the run-to-run spread: timed in alternated pairs the
    ratio sits at 0.94-1.02 while individual pairs range from 0.40 to
    1.63.  Nothing here manufactures a difference that is not there.
    """
    call()
    inner = 1
    while True:
        start = time.perf_counter()
        for _ in range(inner):
            call()
        elapsed = time.perf_counter() - start
        if elapsed >= min_block or inner >= 1 << 20:
            break
        inner = max(inner + 1, int(inner*min_block/max(elapsed, 1e-9)*1.2))

    best = np.inf
    for _ in range(repeat):
        start = time.perf_counter()
        for _ in range(inner):
            call()
        elapsed = time.perf_counter() - start
        best = min(best, elapsed/inner/n_energies*1.0e6)
    return best*1.0


def timed_batch_interleaved(calls, n_energies, repeat=7, min_block=0.05):
    """Time several batched calls against each other, alternating.

    Timing one route to completion and then the next lets machine drift
    masquerade as a cost difference, and here it did: run in blocks, the
    eigensolver came out 1.4x slower than double-double at four flavors,
    while alternated pairs put the ratio at 0.94-1.02 with individual
    pairs ranging from 0.40 to 1.63.  The 1.4x was the ordering.

    So the routes are interleaved, and the direction of the sweep is
    reversed on alternate passes, which cancels a drift that is linear in
    time.  `oscprob4nu`'s own docstring times its two strategies this way
    and says why.
    """
    for call in calls:
        call()

    inner = 1
    while True:
        start = time.perf_counter()
        for _ in range(inner):
            calls[0]()
        elapsed = time.perf_counter() - start
        if elapsed >= min_block or inner >= 1 << 20:
            break
        inner = max(inner + 1, int(inner*min_block/max(elapsed, 1e-9)*1.2))

    best = [np.inf]*len(calls)
    for rep in range(repeat):
        order = range(len(calls)) if rep % 2 == 0 \
            else range(len(calls)-1, -1, -1)
        for i in order:
            start = time.perf_counter()
            for _ in range(inner):
                calls[i]()
            elapsed = time.perf_counter() - start
            best[i] = min(best[i], elapsed/inner/n_energies*1.0e6)
    return best


def series(name, points):
    return {'name': name, 'points': points}


def compiled_series(name, filename, ref):
    """A series read back from one of the compiled drivers.

    NuFast-Earth, Prob3++ and GLoBES are C or C++ and cannot be driven from
    here.  Their drivers live in ``tests/external_drivers`` together with the
    raw output each produced; see the README there for what every one of them
    had to be told before it was propagating through the same Earth as this
    code, and for the two checks each was put through.
    """
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'external_drivers', filename)
    points = []
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            n = int(fields[0])
            us = float(fields[1])
            p = np.array([float(v) for v in fields[2:]])
            points.append({
                'label': str(n),
                'n_shells_per_layer': n,
                'max_abs_error': float(np.max(np.abs(p-ref))),
                'us_per_probability': us})
    return series(name, points)


def build(n_flavors, energies_gev):
    energies = energies_gev*1.0e9
    ref = np.array([reference(e, n_flavors) for e in energies])
    ode = np.array([ode_reference(e, n_flavors) for e in energies])

    # The dial runs past where the other two codes flatten.  The referee is
    # built from n = 128 and 256, so those two points are also the standard
    # Richardson error estimator for themselves; what makes them honest is
    # the independent ODE check on the referee's absolute value, four
    # decades below the finest point plotted.
    # At four flavors the latent roots can come by either of two routes,
    # and they are shown separately because the choice is a real one.  On
    # this problem they are not an accuracy trade at all: the two agree
    # bit for bit at every slab count, because the discretisation error
    # swamps a root difference of 3.6e-17 against 3.9e-16.  What separates
    # them is cost, so the two curves differ horizontally and not
    # vertically.  At three flavors the question does not arise.
    if n_flavors == 3:
        variants = [('NuOscProbExact', None)]
    else:
        variants = [('NuOscProbExact (double-double)', 'double-double'),
                    ('NuOscProbExact (eigensolver)', 'eigensolver')]

    out = []
    def with_strategy(strategy, n):
        """One batched call, under the named root strategy."""
        def call():
            if strategy is not None:
                oscprob4nu.ROOT_STRATEGY = strategy
            return ours(energies, n, n_flavors)
        return call

    pts = [[] for _ in variants]
    for n in (1, 2, 4, 8, 16, 32, 64, 128, 256,
              512, 1024, 2048, 4096, 8192):
        calls = [with_strategy(strategy, n) for _, strategy in variants]
        # accuracy first: order cannot affect it
        errors = []
        for call in calls:
            errors.append(float(np.max(np.abs(call() - ref))))
        # then cost, with the routes alternated against each other
        times = timed_batch_interleaved(calls, len(energies))
        for i, (_, strategy) in enumerate(variants):
            pts[i].append({
                'label': str(n),
                'n_slabs_per_segment': n,
                'root_strategy': strategy,
                'max_abs_error': errors[i],
                'us_per_probability': times[i]})
    for i, (name, _) in enumerate(variants):
        out.append(series(name, pts[i]))

    # The other dial this library exposes.  Since 1.12.0 the Earth routines
    # take rtol/atol and choose the slab count themselves, doubling until
    # the answer settles, which is the same kind of knob nuSQuIDS and
    # nuCraft expose and a fairer thing to compare against them than a raw
    # slab count.  `return_n_slabs` reports what it picked.
    if n_flavors == 3:
        pts = []
        for rtol in RTOLS:
            p, n_used = earth.probabilities_3nu_earth(
                H_VAC_3NU, energies, COSTHZ, electron_fraction=YE,
                n_slabs_per_segment=1, rtol=rtol, n_max=N_MAX_TOLERANCE,
                return_n_slabs=True)
            p = np.asarray(p)[..., 4].ravel()
            pts.append({
                'label': '%.0e' % rtol,
                'rtol': rtol,
                'n_slabs_chosen': int(np.max(n_used)),
                'max_abs_error': float(np.max(np.abs(p-ref))),
                'us_per_probability': timed_batch(
                    lambda r=rtol: earth.probabilities_3nu_earth(
                        H_VAC_3NU, energies, COSTHZ, electron_fraction=YE,
                        n_slabs_per_segment=1, rtol=r,
                        n_max=N_MAX_TOLERANCE), len(energies))})
        out.append(series('NuOscProbExact (tolerance)', pts))
    else:
        pts = []
        for rtol in RTOLS:
            p, n_used = earth.probabilities_4nu_earth(
                H_VAC_4NU, energies, COSTHZ, electron_fraction=YE,
                n_slabs_per_segment=1, rtol=rtol, n_max=N_MAX_TOLERANCE,
                return_n_slabs=True)
            p = np.asarray(p)[..., 5].ravel()
            pts.append({
                'label': '%.0e' % rtol,
                'rtol': rtol,
                'n_slabs_chosen': int(np.max(n_used)),
                'max_abs_error': float(np.max(np.abs(p-ref))),
                'us_per_probability': timed_batch(
                    lambda r=rtol: earth.probabilities_4nu_earth(
                        H_VAC_4NU, energies, COSTHZ, electron_fraction=YE,
                        n_slabs_per_segment=1, rtol=r,
                        n_max=N_MAX_TOLERANCE), len(energies))})
        out.append(series('NuOscProbExact (tolerance)', pts))
    if n_flavors == 4:
        oscprob4nu.ROOT_STRATEGY = 'double-double'      # back to the default

    body = nusquids_body()
    pts = []
    for tol in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-12):
        p = nusquids(energies, body, tol, n_flavors)
        pts.append({
            'label': '%.0e' % tol,
            'tolerance': tol,
            'max_abs_error': float(np.max(np.abs(p-ref))),
            'us_per_probability': timed_batch(
                lambda t=tol: nusquids(energies, body, t, n_flavors),
                len(energies))})
    out.append(series('nuSQuIDS', pts))

    # nuCraft is judged, like every code here, against a referee in its
    # own conventions.  At three flavors the NC term is common to the
    # active flavors and drops out, so the shared referee already is that
    # referee.  At 3+1 it is not: nuCraft's sterile matter entry is an
    # independently rounded constant whose ratio to its CC entry is
    # 0.5016059 where its own formula makes the isoscalar value exactly
    # 1/2 (see the module docstring), and scored against an exact-ratio
    # referee it sits flat at 2.8e-3 however small numPrec is.  So its
    # 3+1 referee carries its own ratio, the plotted error measures its
    # solver -- 3.2e-4 down to 1.2e-9 across the dial -- and the 0.32%
    # inconsistency is reported in the caption, not in the curve.
    if n_flavors == 4:
        ref_nucraft = np.array([
            reference(e, n_flavors, nc_over_cc=NUCRAFT_NC_OVER_CC)
            for e in energies])
        ode_nucraft = np.array([
            ode_reference(e, n_flavors, nc_over_cc=NUCRAFT_NC_OVER_CC)
            for e in energies])
    else:
        ref_nucraft, ode_nucraft = ref, ode
    inst = nucraft_instance(n_flavors)
    pts = []
    for prec in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10):
        p = nucraft(energies, inst, prec)
        pts.append({
            'label': '%.0e' % prec,
            'num_prec': prec,
            'max_abs_error': float(np.max(np.abs(p-ref_nucraft))),
            'us_per_probability': timed_batch(
                lambda q=prec: nucraft(energies, inst, q), len(energies),
                repeat=3)})
    out.append(series('nuCraft', pts))

    # The three compiled codes are three-flavor only: NuFast-Earth takes
    # three angles and two splittings, and Prob3++ and GLoBES both
    # hard-wire double P[3][3].  GLoBES exposes glbRegisterProbabilityEngine
    # and Kopp's snu extension uses it to provide 3+n, but that is separate
    # software with its own conventions, not a flag.
    if n_flavors == 3:
        out.append(compiled_series('NuFast-Earth', 'nufast_earth_prem.txt',
                                   ref))
        out.append(compiled_series('GLoBES', 'globes_prem.txt', ref))
        out.append(compiled_series('Prob3++', 'prob3_prem.txt', ref))

    panel = {
        'energy_gev': energies_gev.tolist(),
        'reference': ref.tolist(),
        'reference_ode_crosscheck': ode.tolist(),
        'reference_vs_ode_max_abs': float(np.max(np.abs(ref-ode))),
        'series': out,
    }
    if n_flavors == 4:
        panel['reference_nucraft'] = ref_nucraft.tolist()
        panel['reference_nucraft_ode_crosscheck'] = ode_nucraft.tolist()
        panel['reference_nucraft_vs_ode_max_abs'] = float(
            np.max(np.abs(ref_nucraft-ode_nucraft)))
    return panel


def main():
    payload = {
        'generated_by': 'tests/prem_scan.py',
        'note': (
            'Speed-accuracy comparison through the Earth on the PREM '
            'profile at cos(theta_z) = -0.9, channel P(nu_mu -> nu_mu). '
            'Accuracy is max |P - P_ref| over twelve energies against a '
            'second-order Richardson extrapolation, (4*P(256)-P(128))/3, '
            'of a 30-digit mpmath slab product; that referee is itself '
            'checked against an adaptive DOP853 integration of the '
            'continuous profile. The y-axis is therefore error against a '
            'CONVERGED PREM solution, not an exact one, and each code sits '
            'where its treatment of the profile puts it rather than where '
            'the accuracy of its probability formula would. Every code is '
            'judged in its own conventions; for nuCraft at 3+1 that means '
            'a referee built with its own sterile-to-CC potential ratio, '
            '0.5016059 (reference_nucraft), so its curve measures its '
            'solver and the 0.32% inconsistency of that ratio is reported '
            'in the caption instead. See the module docstring of '
            'tests/prem_scan.py for every convention matched.'),
        'costhz': COSTHZ,
        'baseline_km': earth.distance_traveled_inside_earth(COSTHZ),
        'electron_fraction': YE,
        'mpmath_dps': DPS,
        'density_scale_nusquids': DENSITY_SCALE_NUSQUIDS,
        'density_scale_nucraft': DENSITY_SCALE_NUCRAFT,
        'dm2_scale_nusquids': DM_SCALE_NUSQUIDS,
        'dm2_scale_nucraft': DM_SCALE_NUCRAFT,
        'nucraft_nc_over_cc': NUCRAFT_NC_OVER_CC,
        'nusquids_version': '1.13.3',
        'nucraft_version': 'r22',
        'three_flavor': build(3, E_GEV_3NU),
        'sterile_3plus1': build(4, E_GEV_4NU),
    }
    print(json.dumps(payload, indent=1))


if __name__ == '__main__':
    main()
