# -*- coding: utf-8 -*-
r"""One 50-digit reference per code, each in that code's own convention.

A code's accuracy number is a statement about its algorithm only if the thing
it is differenced against solves *the problem that code was handed*.  The
previous comparison rescaled every code's density and baseline so that all of
them reproduced this library's potential and this library's kilometre, and
then differenced them against one shared reference.  That made this library
the only code receiving a problem it was designed for, and it capped the
figure at the size of the largest deliberate mismatch.

Here each code is handed honest physical inputs and gets its own reference,
built with:

* its own :math:`\hbar c`, from `conversions.km_to_inv_ev`;
* its own matter normalisation, from `conversions.matter_constant`;
* the same mixing parameters as every other code, in that code's own
  parameterisation, taken from the very doubles the code received;
* the same :math:`Y_e`;
* the profile named by the manifest's ``profile_basis`` for the figure.

What is easy to get wrong here is invisible: a reference in the wrong
convention makes a code look *better* than it is, and nothing downstream
notices.  So every quantity below is either read from `conversions` -- which
reads each code's pinned source -- or is a shared input, and none is typed.

Precision
---------
Two different things limit a reference, and only one of them is arithmetic.

In double precision the Richardson-extrapolated slab product floors near
``2e-12``: the Richardson sequence itself converges cleanly at fourth order,
but a second extrapolation pass gains nothing, because by then the error is
accumulated round-off over the ~1e4 matrix products rather than
discretisation.  No amount of extrapolation gets past that.  Hence mpmath.

The second limit is the inputs, and here the two profile bases differ:

* ``as-handed`` -- the double slab widths and densities the code received
  *are* the problem statement, so converting them to ``mpf`` is exact and is
  the right thing.  Nothing is lost.
* ``continuous`` -- the geometry is part of the answer, so the chord
  crossings and the PREM polynomial are evaluated in extended precision here
  rather than converted from double.

`self_test` checks the achieved precision rather than asserting it.
"""

import os
import sys

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), 'src'))

import conversions                                             # noqa: E402

#: Working precision.  The target is a reference good to 1e-18 or better, and
#: the codes it judges can reach round-off near 1e-16 now that no convention
#: mismatch floors them.  50 digits leaves thirty-odd to spare against the
#: accumulated round-off that limited the double-precision build.
DPS = 50

#: The nu_mu row: P(numu->nue), P(numu->numu), P(numu->nutau).  A row rather
#: than one channel because the constant-density figure plots appearance while
#: the accuracy sweep scores disappearance; a reference holding one cannot
#: serve the other, and the evolution operator gives all three for free.
# P(nu_beta -> nu_alpha) = |U[alpha][beta]|^2 for psi(L) = U psi(0), so the
# nu_mu ROW -- P(numu->nue), P(numu->numu), P(numu->nutau) -- is the second
# COLUMN of U, not its second row.  Taking the row instead agrees on the
# diagonal and is wrong by 1e-2 off it, which is precisely how it hid: the
# disappearance channel the sweep scores was correct and the two appearance
# channels were not.
NUMU_ROW = ((0, 1), (1, 1), (2, 1))

#: Kept for callers that want the single disappearance channel.
SCORED = (1, 1)


def _mpf(x):
    r"""Exact promotion of a double.  ``mp.mpf(float)`` loses nothing."""
    return mp.mpf(float(x))


def pmns(s12sq, s13sq, s23sq, dcp):
    r"""Returns the PMNS matrix in mpmath, built from the doubles given.

    The angles are derived here in extended precision from the same
    ``sin^2`` doubles every code received, rather than from a double-precision
    ``arcsin`` that was then promoted.
    """
    s12, s13, s23 = (mp.sqrt(_mpf(v)) for v in (s12sq, s13sq, s23sq))
    c12, c13, c23 = (mp.sqrt(1 - s*s) for s in (s12, s13, s23))
    e = mp.exp(1j*_mpf(dcp))
    ec = mp.conj(e)
    return mp.matrix([
        [c12*c13,                        s12*c13,                  s13*ec],
        [-s12*c23 - c12*s23*s13*e,  c12*c23 - s12*s23*s13*e,  s23*c13],
        [s12*s23 - c12*c23*s13*e,  -c12*s23 - s12*c23*s13*e,  c23*c13],
    ])


def mass_matrix(s12sq, s13sq, s23sq, dcp, dm21, dm31):
    r"""Returns ``U diag(0, dm21, dm31) U^dagger`` in mpmath, in eV^2."""
    u = pmns(s12sq, s13sq, s23sq, dcp)
    d = mp.diag([mp.mpf(0), _mpf(dm21), _mpf(dm31)])
    return u*d*u.transpose_conj()


def matter_potential(code, density_g_cm3, ye):
    r"""Returns V_CC in eV, in `code`'s own normalisation.

    Every external code writes the potential as ``A = k rhoYe / 2e9`` with
    ``rhoYe`` in g cm^-3, so ``V = k rhoYe / 2e9``.  ``k`` is that code's own
    constant, read from its own pinned source (or, for nuSQuIDS, built from
    its own ``Const``).
    """
    k = conversions.matter_constant(code)
    return _mpf(k)*_mpf(density_g_cm3)*_mpf(ye)/mp.mpf('2e9')


def _propagate(m2, potential_ev, energy_ev, length_inv_ev):
    r"""Returns the evolution operator across one constant-density slab."""
    h = m2/(2*_mpf(energy_ev))
    h[0, 0] += potential_ev
    return mp.expm(-1j*h*length_inv_ev)


def constant_density(code, energies_gev, l_km, density_g_cm3, ye, params):
    r"""Returns P(numu->numu) for a constant-density baseline, exactly.

    No discretisation enters: one slab, one matrix exponential.  The only
    error is arithmetic, and at `DPS` digits that is ~1e-50.  So a
    constant-density reference is exact for every practical purpose, and any
    residual a code shows against it is entirely the code's own.
    """
    mp.mp.dps = DPS
    m2 = mass_matrix(params['s12sq'], params['s13sq'], params['s23sq'],
                     params['dcp_rad'], params['dmsq21_ev2'],
                     params['dmsq31_ev2'])
    v = matter_potential(code, density_g_cm3, ye)
    length = _mpf(l_km)*_mpf(conversions.km_to_inv_ev(code))
    out = []
    for e_gev in energies_gev:
        u = _propagate(m2, v, _mpf(e_gev)*mp.mpf('1e9'), length)
        out.append([abs(u[c])**2 for c in NUMU_ROW])
    return out


def slab_product(code, energy_gev, widths_km, densities_g_cm3, ye, params,
                 m2=None):
    r"""Returns P(numu->numu) across a sequence of constant-density slabs.

    The first slab crossed is the rightmost factor, matching how `earth`
    composes them.
    """
    if m2 is None:
        m2 = mass_matrix(params['s12sq'], params['s13sq'], params['s23sq'],
                         params['dcp_rad'], params['dmsq21_ev2'],
                         params['dmsq31_ev2'])
    km = _mpf(conversions.km_to_inv_ev(code))
    energy_ev = _mpf(energy_gev)*mp.mpf('1e9')
    u = mp.eye(3)
    for width, rho in zip(widths_km, densities_g_cm3):
        v = matter_potential(code, rho, ye)
        u = _propagate(m2, v, energy_ev, _mpf(width)*km)*u
    return [abs(u[c])**2 for c in NUMU_ROW]


def richardson(coarse, fine):
    r"""Returns the second-order Richardson referee ``(4*fine - coarse)/3``.

    The slab product converges at **second** order --- measured ratio 4.00 per
    doubling --- so this is the right combination and ``2*fine - coarse`` is
    not.  The first-order form amplifies the error it is meant to cancel, and
    getting this wrong once already cost a session's work.
    """
    if isinstance(fine, (list, tuple)):
        return [(4*f - c)/3 for f, c in zip(fine, coarse)]
    return (4*fine - coarse)/3


def romberg(values_by_n):
    r"""Returns repeated Richardson extrapolation over a ladder of slab counts.

    The slab product's error expands in even powers of the slab width, so one
    extrapolation removes the ``h^2`` term, a second the ``h^4``, and so on.
    In double precision the second pass bought nothing --- successive
    differences stalled near ``2e-12`` at a ratio of 1.2 instead of 16 ---
    because accumulated round-off had already taken over.  At `DPS` digits
    that floor is gone, the measured first-pass ratio is 16.006, and the
    later passes gain as they should.  That is what makes a chord reference
    below ``1e-18`` reachable at a sane slab count rather than needing tens
    of thousands.

    `values_by_n` maps slab count to value, each count twice the last.  The
    return is the fully extrapolated corner and the triangle it came from,
    so the caller can see the convergence rather than trust it.
    """
    ns = sorted(values_by_n)
    row = [values_by_n[n] for n in ns]
    vector = isinstance(row[0], (list, tuple))
    table = [row]
    j = 1
    while len(row) > 1:
        factor = mp.mpf(4)**j
        if vector:
            row = [[(factor*b - a)/(factor - 1)
                    for a, b in zip(row[k], row[k + 1])]
                   for k in range(len(row) - 1)]
        else:
            row = [(factor*row[k + 1] - row[k])/(factor - 1)
                   for k in range(len(row) - 1)]
        table.append(row)
        j += 1
    return table[-1][0], table


def earth_chord_romberg(code, energy_gev, costhz, ye, params,
                        ladder=(32, 64, 128, 256, 512)):
    r"""Returns the Romberg-extrapolated chord probability and its error bar.

    The error estimate is the gap between the two most extrapolated corners,
    which is the honest thing to quote: it is measured on this problem rather
    than assumed from the convergence order.
    """
    import earth as _earth
    mp.mp.dps = DPS
    m2 = mass_matrix(params['s12sq'], params['s13sq'], params['s23sq'],
                     params['dcp_rad'], params['dmsq21_ev2'],
                     params['dmsq31_ev2'])
    values = {}
    for n in ladder:
        widths, densities = _earth.earth_slabs(costhz, n_slabs_per_segment=n)
        values[n] = slab_product(code, energy_gev, widths, densities, ye,
                                 params, m2=m2)
    corner, table = romberg(values)
    shorter, _ = romberg({n: values[n] for n in ladder[:-1]})
    return corner, abs(corner - shorter), table


# ---------------------------------------------------------------------------
# Chord geometry in extended precision.
#
# The Romberg corner stalled at 1e-17 with the geometry promoted from double.
# That is not arithmetic -- the arithmetic runs at DPS digits -- it is that
# extrapolating n -> infinity is a statement about the CONTINUOUS profile, so
# the crossing radii and the PREM polynomial are part of the answer and their
# double representation error (~1e-16 relative) breaks the smooth h^2
# expansion the extrapolation assumes.  Computed here in mpmath instead.
# ---------------------------------------------------------------------------

def _prem_coeffs_mp():
    r"""Returns this library's PREM polynomial coefficients as mpf."""
    import earth as _earth
    return [[_mpf(a) for a in row] for row in _earth._PREM_COEFFS]


def _prem_boundaries_mp():
    import earth as _earth
    return [_mpf(b) for b in _earth.PREM_BOUNDARIES]


def prem_density_mp(r, coeffs=None, bounds=None, radius=None):
    r"""Returns this library's PREM density at radius `r`, in mpmath.

    Horner in ``x = r/R``, with the shell selected the way ``density_prem``
    selects it (``searchsorted(B, r, side='left')``).
    """
    import globaldefs as gd
    coeffs = coeffs if coeffs is not None else _prem_coeffs_mp()
    bounds = bounds if bounds is not None else _prem_boundaries_mp()
    radius = radius if radius is not None else _mpf(gd.EARTH_RADIUS)
    if r > radius:
        r = radius
    k = 0
    while k < len(bounds) and bounds[k] < r:
        k += 1
    a = coeffs[k]
    x = r/radius
    return a[0] + x*(a[1] + x*(a[2] + x*a[3]))


def chord_edges_mp(costhz, boundaries, radius):
    r"""Returns the chord distances at which `boundaries` are crossed.

    Mirrors ``prem_layer_edges_along_chord`` formula for formula: with
    ``d = -2 R cos``, ``rmin^2 = R^2 (1 - cos^2)`` and
    ``u = -R cos -/+ sqrt(rb^2 - rmin^2)``, the crossing sits at ``d - u``.
    """
    cz = _mpf(costhz)
    d = -2*radius*cz
    rmin2 = radius*radius*(1 - cz*cz)
    out = set()
    for rb in boundaries:
        disc = rb*rb - rmin2
        if disc <= 0:
            continue
        root = mp.sqrt(disc)
        for u in (-radius*cz - root, -radius*cz + root):
            if 0 < u < d:
                out.add(d - u)
    return sorted(out)


def radius_at_mp(costhz, l, radius):
    r"""Returns the radius at chord distance `l`, mirroring `earth`'s form."""
    cz = _mpf(costhz)
    d = -2*radius*cz
    u = d - l
    r2 = radius*radius + u*u + 2*radius*u*cz
    return mp.sqrt(r2) if r2 > 0 else mp.mpf(0)


def chord_slabs_mp(costhz, boundaries, n_per_segment, density_fn,
                   piecewise_constant=False):
    r"""Returns ``(widths, densities)`` along a chord, all in mpmath.

    `boundaries` are the radii at which the density is discontinuous, so
    every one becomes a mandatory slab edge and no slab straddles a jump ---
    the same reason `earth_slabs` cuts there.

    When `piecewise_constant` is true the density is uniform between
    consecutive boundaries, so one slab per segment is not an approximation
    but the exact answer, and `n_per_segment` is ignored.  That is the case
    for the three compiled Earth codes, whose profile is a stack of uniform
    shells: their reference carries no discretisation error at all.
    """
    import globaldefs as gd
    # Set here rather than assumed: mpmath's default is 15 digits, so a
    # caller who reaches the geometry before anything else would silently
    # get a double-precision chord out of a function whose whole purpose is
    # not to be one.
    mp.mp.dps = DPS
    radius = _mpf(gd.EARTH_RADIUS)
    d = -2*radius*_mpf(costhz)
    edges = [mp.mpf(0)] + chord_edges_mp(costhz, boundaries, radius) + [d]

    widths, densities = [], []
    for start, end in zip(edges[:-1], edges[1:]):
        if end <= start:
            continue
        n = 1 if piecewise_constant else n_per_segment
        step = (end - start)/n
        for i in range(n):
            a = start + i*step
            widths.append(step)
            densities.append(density_fn(radius_at_mp(costhz, a + step/2,
                                                     radius)))
    return widths, densities


def squids_akima(x, y):
    r"""Returns SQuIDS' own Akima spline as an evaluator, ported exactly.

    Not scipy's.  ``Akima1DInterpolator`` and this one agree exactly AT the
    data nodes and differ by ~2e-8 between them, because the variants handle
    the end slopes and the weighting differently.  Measured against
    nuSQuIDS' own ``EarthAtm::density`` along a chord, substituting scipy's
    put a ~2e-9 relative error into the profile -- which would have floored
    every nuSQuIDS accuracy number near 1e-9, three orders above the 1e-12
    end of its own tolerance knob, while looking exactly like solver error.

    Ported from ``nuSQuIDS/src/tools.cpp``, ``AkimaSpline::initialize``.
    Evaluated in double because that is what nuSQuIDS evaluates it in: the
    spline it solves against IS a double-precision object, so reproducing
    its arithmetic is fidelity, not a loss.
    """
    x = [float(v) for v in x]
    y = [float(v) for v in y]
    n = len(x)
    if n < 3:
        raise ValueError('the Akima spline needs at least three points')

    def m(i):
        if i == -2:
            return 3*m(0) - 2*m(1)
        if i == -1:
            return 2*m(0) - m(1)
        if i == n - 1:
            return 2*m(n - 2) - m(n - 3)
        if i == n:
            return 3*m(n - 2) - 2*m(n - 3)
        return (y[i + 1] - y[i])/(x[i + 1] - x[i])

    segments = []
    l_ip1, ne_ip1 = -1.0, -1.0
    m_im2 = m_im1 = m_i = m_ip1 = None
    for i in range(n - 1):
        h_i = x[i + 1] - x[i]
        if i == 0:
            m_im2, m_im1, m_i, m_ip1 = m(-2), m(-1), m(0), m(1)
        m_ip2 = m(i + 2)

        if l_ip1 >= 0:                       # carried from the last iteration
            l_i, ne_i = l_ip1, ne_ip1
        else:
            l_i = abs(m_im2 - m_im1)
            ne_i = l_i + abs(m_i - m_ip1)
        t_ri = ((1 - l_i/ne_i)*m_im1 + (l_i/ne_i)*m_i) if ne_i > 0 else m_i

        l_ip1 = abs(m_im1 - m_i)
        ne_ip1 = l_ip1 + abs(m_ip1 - m_ip2)
        t_lip1 = (((1 - l_ip1/ne_ip1)*m_i + (l_ip1/ne_ip1)*m_ip1)
                  if ne_ip1 > 0 else m_i)

        segments.append((x[i], y[i], t_ri,
                         (3*m_i - 2*t_ri - t_lip1)/h_i,
                         (t_ri + t_lip1 - 2*m_i)/(h_i*h_i)))
        m_im2, m_im1, m_i, m_ip1 = m_im1, m_i, m_ip1, m_ip2

    import bisect
    starts = [seg[0] for seg in segments]

    def evaluate(value):
        v = float(value)
        k = max(0, bisect.bisect_right(starts, v) - 1)
        x0, a0, a1, a2, a3 = segments[k]
        t = v - x0
        return (((a3*t) + a2)*t + a1)*t + a0

    return evaluate


def profile_for(code, n_shells=256):
    r"""Returns the profile `code` was actually handed.

    This is the piece the manifest's ``profile_basis: as-handed`` names, and
    getting it wrong is the invisible error: a reference on the wrong profile
    measures a difference of Earth models and calls it an algorithm.

    Returns ``(boundaries, density_fn, piecewise_constant)``.
    """
    import earth as _earth
    import globaldefs as gd
    radius = _mpf(gd.EARTH_RADIUS)
    coeffs, bounds = _prem_coeffs_mp(), _prem_boundaries_mp()

    def continuous(r):
        return prem_density_mp(r, coeffs, bounds, radius)

    if code in ('NuFast-Earth', 'Prob3++', 'GLoBES'):
        # Four major layers, each cut into n equal sub-shells held at their
        # midpoint density -- exactly what nufast_earth.cpp's OurPREM,
        # prob3.cpp's write_profile and globes.cpp's chord() build.
        layers = [bounds[0], bounds[1], bounds[2], radius]
        shell_edges, shell_rho = [], []
        lo = mp.mpf(0)
        for hi in layers:
            for i in range(n_shells):
                shell_edges.append(lo + (i + 1)*(hi - lo)/n_shells)
                shell_rho.append(continuous(lo + (i + mp.mpf('0.5'))
                                            * (hi - lo)/n_shells))
            lo = hi

        def stepped(r):
            k = 0
            while k < len(shell_edges) - 1 and shell_edges[k] < r:
                k += 1
            return shell_rho[k]

        return shell_edges, stepped, True

    if code in ('nuSQuIDS', 'nuCraft'):
        # Both are handed a spline over the same aligned node set: nuCraft an
        # InterpolatedUnivariateSpline(k=1), nuSQuIDS an Akima built inside
        # EarthAtm.  Each code evaluates its spline in double, so the double
        # spline IS the problem it was handed and promoting it is faithful;
        # what must not happen is silently substituting the continuous PREM.
        import numpy as np
        from scipy.interpolate import InterpolatedUnivariateSpline
        edges = np.concatenate(([0.0], np.asarray(_earth.PREM_BOUNDARIES,
                                                  dtype=float)))
        edges = edges[edges <= gd.EARTH_RADIUS]
        if edges[-1] < gd.EARTH_RADIUS:
            edges = np.concatenate((edges, [gd.EARTH_RADIUS]))
        eps = 1.0e-9
        pieces = []
        for a, b in zip(edges[:-1], edges[1:]):
            inner = np.linspace(a, b, 200)
            inner[0], inner[-1] = a + eps, b - eps
            pieces.append(inner)
        nodes = np.unique(np.concatenate(([0.0], np.concatenate(pieces),
                                          [gd.EARTH_RADIUS])))
        rho = _earth.density_prem(np.clip(nodes, 0.0, gd.EARTH_RADIUS))
        if code == 'nuCraft':
            # The adapter hands nuCraft this very object, so the reference
            # and the code share one profile by construction.
            spline = InterpolatedUnivariateSpline(nodes, rho, k=1)
        else:
            # nuSQuIDS builds its OWN Akima from these nodes, so it has to
            # be reproduced rather than approximated with another library's.
            spline = squids_akima(nodes, rho)

        def splined(r):
            return _mpf(float(spline(min(float(r), gd.EARTH_RADIUS))))

        # EVERY spline knot is a mandatory slab edge, not just the PREM
        # boundaries.  A spline is only piecewise smooth: its density has a
        # kink at each knot, so a slab straddling one sees a profile with no
        # convergent h^2 expansion and the extrapolation stalls -- measured,
        # it stalled at 5.8e-10 for nuCraft and 2.8e-10 for nuSQuIDS, eight
        # orders short.  Cutting at the knots puts every slab inside one
        # polynomial piece, which is the same reason `earth_slabs` cuts at
        # PREM boundaries.  The segments are then ~3 km wide already, so the
        # ladder starts at one sub-slab per segment rather than thirty-two.
        return [_mpf(v) for v in nodes if 0.0 < v < gd.EARTH_RADIUS], \
            splined, False

    return bounds, continuous, False


def earth_chord_reference(code, energy_gev, costhz, ye, params,
                          n_shells=256, ladder=(32, 64, 128, 256, 512, 1024)):
    r"""Returns ``(probability, error_estimate)`` for `code` on a chord.

    Built on the continuous PREM in `code`'s own hbar c and matter
    normalisation.  The Romberg ladder extrapolates the slab count to
    infinity, so what is returned is the continuous-profile answer rather
    than any code's discretisation of it.
    """
    mp.mp.dps = DPS

    # CONTINUOUS for every code -- manifest reference_conventions.profile_basis.
    # There is one physical Earth and each code is scored by how far its answer
    # sits from it.  `profile_for` still builds the per-code as-handed profiles
    # and the fidelity audit still uses them, but no production reference is
    # built on one: scoring a code against its own approximation forgives the
    # discretisation, which measured as a spurious six-order advantage for
    # Prob3++ and GLoBES over this library on the same physics.
    import globaldefs as gd
    _coeffs, _bounds = _prem_coeffs_mp(), _prem_boundaries_mp()
    _radius = _mpf(gd.EARTH_RADIUS)
    boundaries, constant = _bounds, False

    def density_fn(r):
        return prem_density_mp(r, _coeffs, _bounds, _radius)

    m2 = mass_matrix(params['s12sq'], params['s13sq'], params['s23sq'],
                     params['dcp_rad'], params['dmsq21_ev2'],
                     params['dmsq31_ev2'])

    if constant:
        widths, densities = chord_slabs_mp(costhz, boundaries, 1, density_fn,
                                           piecewise_constant=True)
        return slab_product(code, energy_gev, widths, densities, ye, params,
                            m2=m2), float(mp.mpf(10)**(-DPS + 5))

    # A profile whose knots are already slab edges needs only a short ladder:
    # the segments are small and the density inside each is a single smooth
    # polynomial.  A continuous profile cut only at the nine PREM boundaries
    # has segments hundreds of km wide and needs the long one.
    #
    # FIVE rungs, not four, and the reason is the error estimate rather than
    # the value.  The estimate below is the standard Romberg one -- the gap
    # between the returned corner and the corner one extrapolation level
    # shallower -- so it reports the error of the corner it DISCARDS, which
    # is conservative and correct but only useful if that shallower corner is
    # itself good.  With four rungs the shallower corner is level 2, and its
    # error is ~3e-10 while the returned level-3 corner is at ~2e-13: the
    # estimate overstated the truth by a factor of 1700 and made two
    # references look eight orders worse than they were.  A fifth rung moves
    # the comparison to the level-3 corner and the quoted figure becomes
    # ~2e-13 (nuCraft) and ~5e-14 (nuSQuIDS), measured.
    #
    # Both are comfortably below nuSQuIDS' tightest tolerance setting, which
    # is what these references have to resolve.  A sixth rung would quote
    # ~1e-15 and costs six minutes an entry; nothing being judged needs it,
    # and the double-precision spline evaluation floors the whole
    # construction near 1e-17 in any case.
    if len(boundaries) > 64:
        ladder = (1, 2, 4, 8, 16)

    values = {}
    for n in ladder:
        widths, densities = chord_slabs_mp(costhz, boundaries, n, density_fn)
        values[n] = slab_product(code, energy_gev, widths, densities, ye,
                                 params, m2=m2)
    # The returned error is the standard Romberg estimate: how far the
    # corner moved when the last rung was added.  It bounds the error of the
    # shallower corner, so it is an upper bound on the returned one.
    corner, _ = romberg(values)
    shorter, _ = romberg({n: values[n] for n in ladder[:-1]})
    if isinstance(corner, (list, tuple)):
        return corner, max(abs(a - b) for a, b in zip(corner, shorter))
    return corner, abs(corner - shorter)


def earth_chord(code, energy_gev, costhz, ye, params, n_slabs=256,
                extrapolate=True):
    r"""Returns P(numu->numu) along a chord, Richardson-extrapolated.

    ``as-handed``: the slabs are the ones `earth` builds and the code was
    given, promoted exactly.  Their double-precision widths are the problem
    statement, not an approximation of it.
    """
    import earth as _earth
    mp.mp.dps = DPS
    m2 = mass_matrix(params['s12sq'], params['s13sq'], params['s23sq'],
                     params['dcp_rad'], params['dmsq21_ev2'],
                     params['dmsq31_ev2'])

    def at(n):
        widths, densities = _earth.earth_slabs(costhz, n_slabs_per_segment=n)
        return slab_product(code, energy_gev, widths, densities, ye, params,
                            m2=m2)

    if not extrapolate:
        return at(n_slabs)
    return richardson(at(n_slabs//2), at(n_slabs))


def self_test(code='NuOscProbExact', costhz=-0.9, energy_gev=10.0, ye=0.5):
    r"""Measures the precision actually achieved, rather than asserting it.

    Returns a dict.  Two independent checks:

    * the constant-density reference recomputed at two working precisions,
      which isolates arithmetic;
    * the Richardson sequence on the chord, whose successive differences
      should keep falling by ~16 -- in double they stalled near 2e-12, which
      is what said the limit was round-off rather than discretisation.
    """
    params = conversions.for_code(code)
    out = {'code': code}

    global DPS
    keep = DPS
    try:
        DPS = 30
        lo = constant_density(code, [energy_gev], 1300.0, 3.0, ye, params)[0]
        DPS = 60
        hi = constant_density(code, [energy_gev], 1300.0, 3.0, ye, params)[0]
        mp.mp.dps = 60
        out['constant_density_dps30_vs_dps60'] = float(abs(hi - lo))
    finally:
        DPS = keep

    mp.mp.dps = DPS
    m2 = mass_matrix(params['s12sq'], params['s13sq'], params['s23sq'],
                     params['dcp_rad'], params['dmsq21_ev2'],
                     params['dmsq31_ev2'])
    import earth as _earth
    raw = {}
    for n in (32, 64, 128, 256):
        widths, densities = _earth.earth_slabs(costhz, n_slabs_per_segment=n)
        raw[n] = slab_product(code, energy_gev, widths, densities, ye, params,
                              m2=m2)
    rich = {n: richardson(raw[n], raw[2*n]) for n in (32, 64, 128)}
    out['richardson_64_vs_128'] = float(abs(rich[64] - rich[128]))
    out['richardson_32_vs_64'] = float(abs(rich[32] - rich[64]))
    out['richardson_gain'] = (out['richardson_32_vs_64']
                              / out['richardson_64_vs_128']
                              if out['richardson_64_vs_128'] else float('inf'))
    return out


if __name__ == '__main__':
    import json
    print(json.dumps(self_test(), indent=2))
