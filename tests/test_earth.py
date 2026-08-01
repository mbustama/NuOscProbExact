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
        np.arcsin(gd.S12_NO_BF), gd.D21_NO_BF)


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
