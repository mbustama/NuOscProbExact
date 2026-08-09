# -*- coding: utf-8 -*-
r"""Tests of the Earth model: PREM, chord geometry, and probabilities.

The Earth routines are checked against things that are known
independently of this code: PREM's own published shell densities and
total mass, the geometry of a sphere, and --- for the probabilities ---
the constant-density result that a slabbed calculation must reproduce
when the slabs are all given the same density.
"""

import numpy as np
import pytest

import earth
import globaldefs as gd
import hamiltonians2nu
import hamiltonians3nu
import hamiltonians4nu
import oscprob3nu
import slabs

from conftest import ATOL


def h_vacuum_3nu():
    r"""Returns the energy-independent three-flavor vacuum Hamiltonian."""
    return hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, gd.DCP_NO_BF,
        gd.D21_NO_BF, gd.D31_NO_BF)


def h_vacuum_2nu():
    r"""Returns the energy-independent two-flavor vacuum Hamiltonian."""
    return hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(
        gd.S12_NO_BF, gd.D21_NO_BF)


def test_dms_to_decimal_handles_the_sign():
    r"""A negative degree does not have its minutes added back."""
    assert earth.dms_to_decimal(10, 30, 0) == pytest.approx(10.5)
    assert earth.dms_to_decimal(-10, 30, 0) == pytest.approx(-10.5)
    assert earth.dms_to_decimal(0, 0, 0) == pytest.approx(0.0)


def test_named_locations_round_trip():
    r"""Every predefined location is retrievable, however it is cased."""
    for name in earth.LOC_COORDS_DMS:
        lat, lon = earth.coordinates_of_named_location(name)
        assert len(lat) == 3 and len(lon) == 3
        assert -90.0 <= earth.dms_to_decimal(*lat) <= 90.0
        assert -180.0 <= earth.dms_to_decimal(*lon) <= 180.0

    assert (earth.coordinates_of_named_location('SOUTH POLE')
            == earth.coordinates_of_named_location('south_pole'))


def test_an_unknown_location_raises():
    r"""The error names the available locations rather than just failing."""
    with pytest.raises(ValueError, match='not a predefined location'):
        earth.coordinates_of_named_location('atlantis')


@pytest.mark.parametrize('r, expected', [
    (0.0, 13.0885),          # centre
    (1221.5, 12.7636),       # inner-core boundary, from the inner-core fit
    (6371.0, 1.0200),        # surface: ocean layer
])
def test_prem_density_at_known_radii(r, expected):
    r"""PREM reproduces its own published values at the shell edges."""
    assert earth.density_prem(r) == pytest.approx(expected, abs=1.e-3)


def test_prem_density_is_monotonic_below_the_lid():
    r"""The Earth does not get denser as you move away from the centre.

    True everywhere below 6151 km, which is where PREM's LID and low
    velocity zone begin.  That shell is deliberately excluded rather
    than tolerated: its density is :math:`2.6910 + 0.6924 x`, whose
    coefficient is *positive*, so PREM really does have the density
    rising with radius there.  A test that swept the whole Earth would
    fail on correct behaviour.
    """
    r = np.linspace(0.0, 6151.0, 20000)
    rho = earth.density_prem(r)

    assert np.all(np.diff(rho) <= 1.e-9)


def test_prem_density_rises_in_the_lid():
    r"""The one shell where density increases outward, checked
    explicitly so the exclusion above is a statement rather than a
    blind spot."""
    # Starting just inside the shell: the bins are right-closed, so r =
    # 6151.0 exactly still belongs to the denser shell below it.
    r = np.linspace(6151.001, 6346.6, 100)
    rho = earth.density_prem(r)

    assert np.all(np.diff(rho) > 0.0)
    # And the boundary itself is a genuine drop, not a rise
    assert earth.density_prem(6151.0) > earth.density_prem(6151.001)


def test_prem_integrates_to_the_mass_of_the_earth():
    r"""The strongest single check on the density profile.

    Integrating 4 pi r^2 rho(r) over the Earth must give its known mass,
    5.972e24 kg.  This catches a wrong coefficient, a mis-ordered shell,
    or an off-by-one in the boundary lookup, none of which the spot
    checks above would notice.
    """
    r = np.linspace(0.0, gd.EARTH_RADIUS, 200001)       # km
    rho = earth.density_prem(r)                          # g cm^-3

    # 4 pi r^2 rho dr, with r in km -> cm and rho in g cm^-3, giving grams.
    # The trapezoid is written out rather than taken from NumPy: `trapz`
    # was removed in NumPy 2.0 and `trapezoid` does not exist before it,
    # and the test matrix spans both.
    integrand = 4.0*np.pi*(r*1.e5)**2.0*rho
    x = r*1.e5
    mass_g = np.sum(0.5*(integrand[1:]+integrand[:-1])*np.diff(x))

    assert mass_g*1.e-3 == pytest.approx(5.972e24, rel=1.e-3)


def test_density_outside_the_earth_raises():
    r"""A radius beyond the surface is an error, not an extrapolation."""
    with pytest.raises(ValueError, match='cannot exceed'):
        earth.density_prem(gd.EARTH_RADIUS*1.01)

    with pytest.raises(ValueError, match='cannot be negative'):
        earth.density_prem(-1.0)

    # ... but a radius a hair outside, from floating-point arithmetic on a
    # chord endpoint, is clamped rather than rejected.
    assert earth.density_prem(gd.EARTH_RADIUS*(1.0 + 1.e-12)) > 0.0


def test_density_accepts_an_array():
    r"""Array input gives the same answers as scalar input."""
    r = np.array([0.0, 3000.0, 6000.0])
    assert np.allclose(earth.density_prem(r),
                       [earth.density_prem(float(x)) for x in r])


def test_matter_potential_reproduces_the_crust_constant():
    r"""The conversion agrees with globaldefs' own crust potential.

    globaldefs computes VCC_EARTH_CRUST from a 3 g/cm^3 crust by the
    same chain of constants; if this disagrees, one of the two is wrong.
    """
    assert earth.matter_potential(gd.DENSITY_MATTER_CRUST_G_PER_CM3) == \
        pytest.approx(gd.VCC_EARTH_CRUST, rel=1.e-12)


def test_matter_potential_is_linear_in_density():
    r"""Doubling the density doubles the potential."""
    assert earth.matter_potential(6.0) == \
        pytest.approx(2.0*earth.matter_potential(3.0), rel=1.e-12)
    assert np.allclose(earth.matter_potential(np.array([1.0, 2.0])),
                       [earth.matter_potential(1.0),
                        earth.matter_potential(2.0)])


@pytest.mark.parametrize('costhz, expected', [
    (-1.0, 2.0*gd.EARTH_RADIUS),
    (-0.5, gd.EARTH_RADIUS),
    (0.0, 0.0),
    (0.7, 0.0),
])
def test_chord_length_for_a_direction(costhz, expected):
    r"""The chord length is -2 R costhz, and zero going upward."""
    assert earth.distance_traveled_inside_earth(costhz) == \
        pytest.approx(expected)


def test_radius_along_a_chord_is_symmetric():
    r"""A chord is symmetric about its midpoint, and starts and ends on
    the surface."""
    costhz = -0.6
    d = earth.distance_traveled_inside_earth(costhz)
    l = np.linspace(0.0, d, 51)
    r = earth.earth_radial_distance_from_depth(costhz, l)

    assert r[0] == pytest.approx(gd.EARTH_RADIUS, abs=1.e-6)
    assert r[-1] == pytest.approx(gd.EARTH_RADIUS, abs=1.e-6)
    assert np.allclose(r, r[::-1], atol=1.e-6)
    # The closest approach is R sin(theta) = R sqrt(1 - costhz^2)
    assert r.min() == pytest.approx(
        gd.EARTH_RADIUS*np.sqrt(1.0-costhz**2.0), rel=1.e-3)


def test_radius_along_a_chord_rejects_bad_distances():
    r"""Distances off the end of the chord are errors."""
    with pytest.raises(ValueError, match='cannot exceed'):
        earth.earth_radial_distance_from_depth(-1.0, 1.e5)
    with pytest.raises(ValueError, match='cannot be negative'):
        earth.earth_radial_distance_from_depth(-1.0, -1.0)


def test_prem_edges_lie_on_boundaries():
    r"""Every crossing returned really is at a PREM boundary radius."""
    costhz = -0.85
    edges = earth.prem_layer_edges_along_chord(costhz)
    d = earth.distance_traveled_inside_earth(costhz)

    assert len(edges) > 0
    assert np.all((edges > 0.0) & (edges < d))
    assert np.all(np.diff(edges) > 0.0)

    r_at_edges = earth.earth_radial_distance_from_depth(costhz, edges)
    for r in r_at_edges:
        assert np.min(np.abs(earth.PREM_BOUNDARIES - r)) < 1.e-6


def test_a_shallow_chord_crosses_fewer_shells():
    r"""A chord that skims the surface never reaches the core."""
    deep = earth.prem_layer_edges_along_chord(-1.0)
    shallow = earth.prem_layer_edges_along_chord(-0.05)

    assert len(shallow) < len(deep)
    assert len(earth.prem_layer_edges_along_chord(0.3)) == 0


def test_earth_slabs_tile_the_chord():
    r"""The slabs cover the whole chord, exactly once."""
    costhz = -0.75
    widths, densities = earth.earth_slabs(costhz, n_slabs_per_segment=4)

    assert widths.shape == densities.shape
    assert np.all(widths > 0.0)
    assert widths.sum() == pytest.approx(
        earth.distance_traveled_inside_earth(costhz), rel=1.e-12)
    assert np.all(densities > 0.0)


def test_earth_slabs_refine_with_more_sub_slabs():
    r"""More sub-slabs means more slabs, still tiling the same chord."""
    coarse, _ = earth.earth_slabs(-0.9, n_slabs_per_segment=2)
    fine, _ = earth.earth_slabs(-0.9, n_slabs_per_segment=16)

    assert len(fine) == 8*len(coarse)
    assert fine.sum() == pytest.approx(coarse.sum(), rel=1.e-12)


def test_earth_slabs_rejects_upward_directions():
    r"""There is no path through the Earth for costhz >= 0."""
    with pytest.raises(ValueError, match='must be negative'):
        earth.earth_slabs(0.0)
    with pytest.raises(ValueError, match='at least 1'):
        earth.earth_slabs(-1.0, n_slabs_per_segment=0)


def test_costhz_between_two_locations():
    r"""The zenith angle of a chord matches its length."""
    lat1, lon1 = earth.coordinates_of_named_location('north_pole')
    lat2, lon2 = earth.coordinates_of_named_location('south_pole')

    chord = earth.chord_length_inside_earth(lat1, lon1, lat2, lon2)
    costhz = earth.costhz_between_points_on_surface(lat1, lon1, lat2, lon2)

    # Pole to pole is straight through the centre
    assert chord == pytest.approx(2.0*gd.EARTH_RADIUS, rel=1.e-9)
    assert costhz == pytest.approx(-1.0, rel=1.e-9)


def test_a_location_to_itself_has_no_chord():
    r"""Two coincident points give a degenerate trajectory, and say so."""
    with pytest.raises(ValueError, match='coincide'):
        earth.probabilities_3nu_between_locations(
            h_vacuum_3nu(), 1.e9, 'kamioka', 'kamioka')


def test_earth_probabilities_are_normalized():
    r"""Each initial flavor's probabilities sum to one."""
    prob = np.array(earth.probabilities_3nu_earth(
        h_vacuum_3nu(), 1.e9, -0.8, n_slabs_per_segment=4))

    assert len(prob) == 9
    assert np.all(prob >= 0.0)
    for start in (0, 3, 6):
        assert prob[start:start+3].sum() == pytest.approx(1.0, abs=ATOL)

    prob2 = np.array(earth.probabilities_2nu_earth(
        h_vacuum_2nu(), 1.e9, -0.8, n_slabs_per_segment=4))
    assert len(prob2) == 4
    assert prob2[0]+prob2[1] == pytest.approx(1.0, abs=ATOL)


def test_earth_probabilities_converge_with_refinement():
    r"""The answer settles as the slabs are refined, at second order.

    The discretisation is the only approximation in the calculation, so
    this is the check that it is under control.  It is deliberately made
    in the asymptotic regime, from 32 sub-slabs upward: at coarse
    resolution successive differences are not reliably monotone --- at
    1 GeV the step from 8 to 16 happens to be larger than the one from 4
    to 8, because the value at 8 lands close to the value at 4 by
    coincidence --- and a test anchored there would fail on correct
    behaviour.

    Sampling the density at the midpoint of each slab is second-order
    accurate, so each doubling should cut the error by about four.
    """
    args = (h_vacuum_3nu(), 1.e9, -0.8)
    p = {n: np.array(earth.probabilities_3nu_earth(*args,
                                                   n_slabs_per_segment=n))
         for n in (32, 64, 128, 256)}

    d64 = np.max(np.abs(p[64]-p[32]))
    d128 = np.max(np.abs(p[128]-p[64]))
    d256 = np.max(np.abs(p[256]-p[128]))

    assert d256 < d128 < d64
    # Second order means a ratio near four; allow a generous band, since
    # what is being checked is the order, not a precise constant.
    for coarse, fine in ((d64, d128), (d128, d256)):
        assert 2.5 < coarse/fine < 6.0


def test_a_uniform_earth_reproduces_the_constant_density_result():
    r"""With one density everywhere, slabbing must add nothing.

    Replacing PREM by a constant makes the Hamiltonian constant along
    the whole chord, so the slabbed answer has to equal the ordinary
    single-Hamiltonian one.  This is what confirms that the slab
    machinery, the unit conversions and the potential conversion are all
    wired together correctly, independently of PREM itself.
    """
    costhz = -0.7
    energy = 1.e9
    density = 5.0

    widths_km, _ = earth.earth_slabs(costhz, n_slabs_per_segment=3)
    vcc = earth.matter_potential(density)

    h_uniform = hamiltonians3nu.hamiltonian_3nu_matter(
        h_vacuum_3nu(), energy, vcc)
    stack = np.repeat(np.asarray(h_uniform)[None, ...], len(widths_km),
                      axis=0)

    slabbed = np.array(slabs.probabilities_3nu_slabs(
        stack, widths_km*gd.CONV_KM_TO_INV_EV))
    direct = np.array(oscprob3nu.probabilities_3nu(
        h_uniform, widths_km.sum()*gd.CONV_KM_TO_INV_EV))

    assert np.allclose(slabbed, direct, atol=ATOL)


def test_between_locations_matches_the_explicit_chord():
    r"""The named-location wrapper is the costhz calculation, no more."""
    lat1, lon1 = earth.coordinates_of_named_location('cern')
    lat2, lon2 = earth.coordinates_of_named_location('kamioka')
    costhz = earth.costhz_between_points_on_surface(lat1, lon1, lat2, lon2)

    h_vac = h_vacuum_3nu()
    via_names = np.array(earth.probabilities_3nu_between_locations(
        h_vac, 1.e9, 'cern', 'kamioka', n_slabs_per_segment=3))
    via_costhz = np.array(earth.probabilities_3nu_earth(
        h_vac, 1.e9, costhz, n_slabs_per_segment=3))

    assert np.allclose(via_names, via_costhz, atol=ATOL)

    h_vac2 = h_vacuum_2nu()
    assert np.allclose(
        np.array(earth.probabilities_2nu_between_locations(
            h_vac2, 1.e9, 'cern', 'kamioka', n_slabs_per_segment=3)),
        np.array(earth.probabilities_2nu_earth(
            h_vac2, 1.e9, costhz, n_slabs_per_segment=3)),
        atol=ATOL)


# --------------------------------------------------------------------------
# Four flavors through the Earth, which was not possible before 1.11.0
# --------------------------------------------------------------------------

def _vacuum_4nu(s14=0.0, s24=0.0, s34=0.0, dm41=1.0):
    r"""Returns a 3+1 energy-independent vacuum Hamiltonian."""
    return hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(
        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, s14, s24, s34,
        gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF, dm41)


def _vacuum_3nu():
    r"""Returns the matching three-flavor vacuum Hamiltonian."""
    return np.asarray(
        hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
            gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, gd.DCP_NO_BF,
            gd.D21_NO_BF, gd.D31_NO_BF))


@pytest.mark.parametrize('costhz', [-1.0, -0.8, -0.3, -0.05])
def test_a_decoupled_sterile_crosses_the_earth_as_three_flavors_do(costhz):
    r"""With the sterile angles off, the active block is the 3nu answer.

    The strongest check available for the four-flavor Earth path, and
    the reason it is worth having: it runs the whole PREM chord --- the
    boundary crossings, the per-slab densities, the composition of
    nineteen segments' worth of operators --- through machinery that had
    never been exercised at four flavors, and compares it against a
    module that has.  A wrong sterile matter entry, a wrong slab order
    or a wrong flavor ordering would all show here.
    """
    four = np.asarray(earth.probabilities_4nu_earth(
        _vacuum_4nu(), 1.0e10, costhz)).reshape(4, 4)
    three = np.asarray(earth.probabilities_3nu_earth(
        _vacuum_3nu(), 1.0e10, costhz)).reshape(3, 3)

    assert np.max(np.abs(four[:3, :3] - three)) < 1.0e-11
    assert np.max(np.abs(four[:3, 3])) == 0.0
    assert np.max(np.abs(four[3, :3])) == 0.0
    assert np.allclose(four.sum(axis=1), 1.0, atol=1.0e-9)


def test_the_sterile_state_feels_the_neutral_current_potential():
    r"""V_NC does not cancel once a sterile state is present.

    It is flavor-universal across the three active states, so at three
    flavors it is proportional to the identity and drops out.  A sterile
    state does not feel it, so removing it from all four leaves
    ``-V_NC`` on the sterile entry --- and that entry is what places the
    sterile matter resonance.  Dropping it would be invisible in vacuum,
    which is exactly why it is pinned here.
    """
    assert earth.matter_potential_nc(3.0) < 0.0
    assert np.isclose(earth.matter_potential_nc(3.0), gd.VNC_EARTH_CRUST)

    # Twice the density is twice the potential, and the ratio to V_CC is
    # fixed by the isoscalar assumption
    assert np.isclose(earth.matter_potential_nc(6.0),
                      2.0*earth.matter_potential_nc(3.0))
    assert np.isclose(-earth.matter_potential_nc(3.0)
                      / earth.matter_potential(3.0), 0.5)

    mixed = _vacuum_4nu(s14=np.sqrt(0.10), s24=np.sqrt(0.10))
    with_nc = np.asarray(earth.probabilities_4nu_earth(mixed, 1.0e10, -0.8))
    assert abs(sum(with_nc[0:4]) - 1.0) < 1.0e-9
    # The sterile channel is actually open, or this test proves nothing
    assert with_nc[3] > 1.0e-6


def test_4nu_earth_converges_as_the_slabs_are_refined():
    r"""Refining the discretisation settles, as at three flavors."""
    h_vacuum = _vacuum_4nu(s14=np.sqrt(0.10), s24=np.sqrt(0.10))
    values = [np.asarray(earth.probabilities_4nu_earth(
        h_vacuum, 1.0e10, -0.8, n_slabs_per_segment=n))[0]
        for n in (2, 4, 8, 16, 32)]

    differences = np.abs(np.diff(values))
    assert differences[-1] < differences[0]
    assert differences[-1] < 1.0e-3


def test_4nu_between_locations_matches_the_chord_it_names():
    r"""The named-pair wrapper is the chord calculation, not a new one."""
    h_vacuum = _vacuum_4nu(s14=np.sqrt(0.10))
    costhz = earth.costhz_between_points_on_surface(
        *earth.coordinates_of_named_location('fermilab'),
        *earth.coordinates_of_named_location('homestake'))

    named = earth.probabilities_4nu_between_locations(
        h_vacuum, 1.0e9, 'fermilab', 'homestake')
    direct = earth.probabilities_4nu_earth(h_vacuum, 1.0e9, costhz)

    assert np.allclose(named, direct, atol=0.0)


def test_a_coordinate_smaller_than_one_degree_south_is_expressible():
    r"""The sign can ride on the minutes when the degree part is zero.

    ``dms_to_decimal`` puts the sign on the degree part, which works
    until the degree part is zero: ``(-0, 30, 0)`` is ``(0, 30, 0)`` in
    Python, so 0 deg 30' South silently came back as 0 deg 30' North.
    No predefined location reaches that band, which is why it went
    unnoticed, but a user-supplied one can.
    """
    assert earth.dms_to_decimal(0, -30, 0) == -0.5
    assert earth.dms_to_decimal(0, 30, 0) == 0.5
    assert earth.dms_to_decimal(-5, 30, 0) == -5.5
    assert earth.dms_to_decimal(5, 30, 0) == 5.5

    # Every predefined location still round-trips
    for name in earth.LOC_COORDS_DMS:
        lat, lon = earth.coordinates_of_named_location(name)
        assert -90.0 <= earth.dms_to_decimal(*lat) <= 90.0
        assert -180.0 <= earth.dms_to_decimal(*lon) <= 180.0


def test_an_unknown_location_says_what_the_known_ones_are():
    r"""And does not bury that behind a chained KeyError."""
    with pytest.raises(ValueError, match='is not a predefined location'):
        earth.coordinates_of_named_location('atlantis')

    try:
        earth.coordinates_of_named_location('atlantis')
    except ValueError as error:
        assert error.__cause__ is None
        assert error.__suppress_context__ is True
        assert 'kamioka' in str(error)


def test_the_radius_along_a_chord_stays_real_at_the_endpoints():
    r"""Round-off at the endpoints must not surface as a nan.

    ``r2`` vanishes nowhere on a chord, but the expression for it can
    land a few ulp below zero at ``l = 0`` and ``l = d``.  It is clipped
    at zero rather than passed through ``abs``, which would also hide a
    genuinely negative value.
    """
    for costhz in (-1.0, -0.5, -1.0e-9):
        d = earth.distance_traveled_inside_earth(costhz)
        for l in (0.0, d, d/2.0):
            r = earth.earth_radial_distance_from_depth(costhz, l)
            assert np.isfinite(r)
            assert 0.0 <= r <= gd.EARTH_RADIUS*(1.0 + 1.0e-12)
        assert np.isclose(earth.earth_radial_distance_from_depth(costhz, 0.0),
                          gd.EARTH_RADIUS)
        assert np.isclose(earth.earth_radial_distance_from_depth(costhz, d),
                          gd.EARTH_RADIUS)


def test_the_neutron_fraction_can_be_given_explicitly():
    r"""A profile that is not isoscalar can say so.

    The default takes the neutron fraction as the complement of the
    electron fraction, which is the isoscalar assumption; a caller
    modelling something else supplies it directly.
    """
    isoscalar = earth.matter_potential_nc(3.0)

    assert np.isclose(earth.matter_potential_nc(3.0, neutron_fraction=0.5),
                      isoscalar)
    assert np.isclose(earth.matter_potential_nc(3.0, neutron_fraction=1.0),
                      2.0*isoscalar)
    # Zero neutrons, no neutral-current potential at all
    assert earth.matter_potential_nc(3.0, neutron_fraction=0.0) == 0.0


@pytest.mark.parametrize('costhz', [-1.5, -2.0, 1.5, -1.0e6, 42.0])
def test_a_costhz_outside_its_range_is_refused(costhz):
    r"""``costhz`` is a cosine, and a chord cannot exceed the diameter.

    The chord length is ``-2 R costhz``, an expression happy to accept
    anything: ``costhz = -1.5`` gave 19 113 km against an Earth diameter
    of 12 742 km, and that chord then acquired seventy-six
    plausible-looking slabs with densities spanning the whole PREM
    range.  Nothing downstream noticed.
    """
    with pytest.raises(ValueError, match=r'must lie in \[-1, 1\]'):
        earth.distance_traveled_inside_earth(costhz)

    if costhz < 0.0:
        with pytest.raises(ValueError, match=r'must lie in \[-1, 1\]'):
            earth.earth_slabs(costhz, 4)


def test_the_boundaries_of_costhz_are_still_allowed():
    r"""Inclusive: the diametric chord and the horizon are both real."""
    assert np.isclose(earth.distance_traveled_inside_earth(-1.0),
                      2.0*gd.EARTH_RADIUS)
    assert earth.distance_traveled_inside_earth(1.0) == 0.0
    assert earth.distance_traveled_inside_earth(0.0) == 0.0
    assert len(earth.earth_slabs(-1.0, 2)[0]) > 0


# Published baselines of real experiments, in km.  These are what the
# chord between the two sites is *for*, and they are the only external
# check on the haversine geometry: every other test here compares the
# geometry against itself.  The tolerance absorbs the spherical-Earth
# approximation and the fact that a published baseline is quoted between
# beam target and detector hall rather than between site centroids.
PUBLISHED_BASELINES = [
    ('cern', 'gran_sasso', 732.0),      # CNGS
    ('fermilab', 'homestake', 1300.0),  # DUNE
    ('tokai', 'kamioka', 295.0),        # T2K
    ('cern', 'pyhaasalmi', 2288.0),     # LBNO, as proposed
]


@pytest.mark.parametrize('site_1,site_2,published', PUBLISHED_BASELINES)
def test_named_baselines_match_the_published_ones(site_1, site_2, published):
    r"""The chord between two sites is the experiment's baseline.

    Nothing else in this module checks the haversine against anything
    outside the code: replacing ``cos(lat2)`` by ``sin(lat2)`` in it left
    the whole suite green, which means the geometry of every
    ``*_between_locations`` call rested on no external evidence at all.
    """
    chord = earth.chord_length_inside_earth(
        *earth.coordinates_of_named_location(site_1),
        *earth.coordinates_of_named_location(site_2))

    assert abs(chord - published)/published < 0.02, (
        '%s to %s came out %.1f km against a published %.1f km'
        % (site_1, site_2, chord, published))


def test_the_haversine_is_right_at_the_poles_too():
    r"""A pole is where a wrong latitude term shows up most.

    Mid-latitude pairs are forgiving --- ``cos`` and ``sin`` are close
    together near 45 degrees, which is where most of these sites are ---
    so the geometry is pinned somewhere it cannot hide: a chord to the
    South Pole, and the diameter between the two poles.
    """
    poles = earth.chord_length_inside_earth(
        *earth.coordinates_of_named_location('north_pole'),
        *earth.coordinates_of_named_location('south_pole'))
    assert np.isclose(poles, 2.0*gd.EARTH_RADIUS)

    # CERN sits at 46 deg 14' N, so its angular distance to the South
    # Pole is 136.23 deg and the chord is 2 R sin(half of that)
    cern_to_pole = earth.chord_length_inside_earth(
        *earth.coordinates_of_named_location('cern'),
        *earth.coordinates_of_named_location('south_pole'))
    expected = 2.0*gd.EARTH_RADIUS*np.sin(np.radians(90.0 + 46.2338)/2.0)
    assert np.isclose(cern_to_pole, expected, rtol=1.0e-4)

    # A site is zero from itself, whatever its latitude
    for name in ('north_pole', 'kamioka', 'south_pole'):
        coords = earth.coordinates_of_named_location(name)
        assert earth.chord_length_inside_earth(*coords, *coords) < 1.0e-9


@pytest.mark.parametrize('boundary', list(earth.PREM_BOUNDARIES))
def test_a_radius_exactly_on_a_shell_boundary_takes_the_inner_shell(boundary):
    r"""The PREM bins are right-closed, and the code says so.

    ``density_prem`` looks the shell up with ``side='left'``, whose
    comment explains that this reproduces bins of the form
    ``r <= 1221.5``.  Nothing tested it: switching to ``side='right'``
    left the suite green, because a slab midpoint never lands exactly on
    a boundary.  A caller evaluating the density *at* a boundary is
    entitled to the documented answer.
    """
    on = earth.density_prem(boundary)
    just_inside = earth.density_prem(boundary*(1.0 - 1.0e-12))
    just_outside = earth.density_prem(boundary*(1.0 + 1.0e-12))

    assert np.isclose(on, just_inside, rtol=1.0e-9)

    # Eight of the nine boundaries are genuine density discontinuities,
    # jumping by between 1.9% and 61%, so the two sides must differ
    # there --- without that, this test would pass on either convention
    # and prove nothing.  The ninth, 5771 km, is not a discontinuity:
    # PREM changes polynomial there but the two meet, agreeing to
    # 5.9e-8 relative, which is continuity to the precision the
    # published coefficients are quoted at.  It is a change of fit
    # inside the transition zone, not a change of rock.
    if boundary != 5771.0:
        assert not np.isclose(on, just_outside, rtol=1.0e-3)
    else:
        assert np.isclose(on, just_outside, rtol=1.0e-6)


# --- Antineutrinos through the Earth ----------------------------------
#
# Before 1.13.1 an antineutrino crossing had to be built by hand, as
# notebooks/13_antineutrinos.ipynb showed: conjugate the vacuum
# Hamiltonian, negate every potential, then compose the slabs.  The flag
# has to reproduce exactly that, and keep the batched PREM path.


def _vacuum_hamiltonian(n_flavors):
    r"""Returns the energy-independent vacuum Hamiltonian at any n."""
    if n_flavors == 2:
        return h_vacuum_2nu()
    if n_flavors == 3:
        return h_vacuum_3nu()
    return _vacuum_4nu(s14=np.sqrt(0.1), s24=np.sqrt(0.1))


def _antineutrino_by_hand(h_vac, energy, costhz, n_flavors,
                          n_slabs_per_segment=4):
    r"""The hand-built idiom the flag replaces."""
    widths_km, densities = earth.earth_slabs(costhz, n_slabs_per_segment)
    potentials = [-earth.matter_potential(densities)]
    if n_flavors == 4:
        potentials.append(-earth.matter_potential_nc(densities))
    builder = {2: hamiltonians2nu.hamiltonian_2nu_matter,
               3: hamiltonians3nu.hamiltonian_3nu_matter,
               4: hamiltonians4nu.hamiltonian_4nu_matter}[n_flavors]
    h = builder(np.conj(np.asarray(h_vac, dtype=complex)), energy,
                *potentials)
    return np.asarray({2: slabs.probabilities_2nu_slabs,
                       3: slabs.probabilities_3nu_slabs,
                       4: slabs.probabilities_4nu_slabs}[n_flavors](
                           h, widths_km*gd.CONV_KM_TO_INV_EV))


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
@pytest.mark.parametrize('costhz', [-1.0, -0.9, -0.3])
def test_the_antineutrino_flag_is_the_hand_built_construction(n_flavors,
                                                              costhz):
    r"""The flag reproduces conjugate-and-flip-both, exactly."""
    h_vac = _vacuum_hamiltonian(n_flavors)
    energy = 5.0e9

    routine = {2: earth.probabilities_2nu_earth,
               3: earth.probabilities_3nu_earth,
               4: earth.probabilities_4nu_earth}[n_flavors]
    got = np.asarray(routine(h_vac, energy, costhz, n_slabs_per_segment=4,
                             antineutrino=True))
    expected = _antineutrino_by_hand(h_vac, energy, costhz, n_flavors)

    assert np.allclose(got, expected, rtol=0.0, atol=1.0e-14)


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
def test_antineutrinos_differ_from_neutrinos_in_matter(n_flavors):
    r"""The flag is not a no-op: matter separates the two.

    The energy is chosen per flavor count to sit near that case's
    resonance: the two-flavor Hamiltonian here carries the solar
    splitting, whose MSW resonance in the Earth is around 100 MeV, not
    the few GeV where the atmospheric one sits.
    """
    h_vac = _vacuum_hamiltonian(n_flavors)
    energy = {2: 1.0e8, 3: 5.0e9, 4: 5.0e9}[n_flavors]

    routine = {2: earth.probabilities_2nu_earth,
               3: earth.probabilities_3nu_earth,
               4: earth.probabilities_4nu_earth}[n_flavors]
    nu = np.asarray(routine(h_vac, energy, -0.9, n_slabs_per_segment=4))
    bar = np.asarray(routine(h_vac, energy, -0.9, n_slabs_per_segment=4,
                             antineutrino=True))

    assert np.abs(nu - bar).max() > 1.0e-3


def test_conjugating_alone_is_not_enough():
    r"""Both operations are needed, which is the point of the flag.

    Conjugating the vacuum term without reversing the potential is the
    classic error; it must not agree with the real antineutrino answer.
    """
    h_vac = _vacuum_hamiltonian(3)
    energy = 5.0e9

    bar = np.asarray(earth.probabilities_3nu_earth(
        h_vac, energy, -0.9, n_slabs_per_segment=4, antineutrino=True))
    conjugate_only = np.asarray(earth.probabilities_3nu_earth(
        np.conj(np.asarray(h_vac, dtype=complex)), energy, -0.9,
        n_slabs_per_segment=4))

    assert np.abs(bar - conjugate_only).max() > 1.0e-3


def test_the_antineutrino_grid_is_batched_like_the_neutrino_one():
    r"""A whole oscillogram, one call, with the flag set."""
    h_vac = _vacuum_hamiltonian(3)
    energies = np.logspace(0.0, 1.0, 5)*1.0e9
    costhz = np.array([-0.95, -0.6, -0.2])

    grid = earth.probabilities_3nu_earth(
        h_vac, energies[:, None], costhz[None, :], n_slabs_per_segment=4,
        antineutrino=True)
    assert grid.shape == (5, 3, 9)

    for i, e in enumerate(energies):
        for j, c in enumerate(costhz):
            one = np.asarray(earth.probabilities_3nu_earth(
                h_vac, e, c, n_slabs_per_segment=4, antineutrino=True))
            assert np.allclose(grid[i, j], one, atol=1.0e-14)


def test_antineutrinos_between_locations_carry_the_flag():
    r"""The named-site wrappers pass it through."""
    h_vac = _vacuum_hamiltonian(3)
    energy = 5.0e9

    bar = np.asarray(earth.probabilities_3nu_between_locations(
        h_vac, energy, 'south_pole', 'kamioka', n_slabs_per_segment=4,
        antineutrino=True))
    costhz = earth._costhz_of_named_pair('south_pole', 'kamioka')
    expected = _antineutrino_by_hand(h_vac, energy, costhz, 3)

    assert np.allclose(bar, expected, rtol=0.0, atol=1.0e-14)


def test_antineutrinos_reach_the_compiled_four_flavor_kernel():
    r"""The batched kernel path, at four flavors, for antineutrinos.

    `_probabilities_earth_batch` reverses both potentials inside its
    compiled branch, and only the four-flavor kernel is handed the
    neutral-current one.  Reaching that line needs all three at once: a
    stack big enough for the kernel to be worth dispatching to, four
    flavors, and ``antineutrino=True``.
    """
    h_vac = _vacuum_hamiltonian(4)
    energies = np.logspace(2.0, 4.0, 64)*1.0e9
    costhz = -0.8

    bar = earth.probabilities_4nu_earth(h_vac, energies, costhz,
                                        n_slabs_per_segment=4,
                                        antineutrino=True)
    assert bar.shape == (64, 16)

    expected = np.stack([_antineutrino_by_hand(h_vac, e, costhz, 4)
                         for e in energies])
    assert np.allclose(bar, expected, rtol=0.0, atol=1.0e-12)

    nu = earth.probabilities_4nu_earth(h_vac, energies, costhz,
                                       n_slabs_per_segment=4)
    assert np.abs(nu - bar).max() > 1.0e-3
