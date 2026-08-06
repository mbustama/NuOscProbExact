# -*- coding: utf-8 -*-
r"""The Earth as a sequence of slabs: PREM, chord geometry, and probabilities.

A neutrino crossing the Earth sees a matter density that changes
continuously along its path, so the Hamiltonian is not constant and the
exact expansions of :mod:`oscprob2nu` and :mod:`oscprob3nu` do not apply
to the trajectory as a whole.  They apply to any piece of it over which
the density is taken to be constant, which is what this module builds:
the chord is cut into slabs, each slab is solved exactly by
:mod:`slabs`, and the operators are multiplied.

The density comes from the Preliminary Reference Earth Model (PREM)
:cite:`Dziewonski:1981xy`, a piecewise-polynomial fit to seismological
data, given as :math:`\rho(x)` with :math:`x = r/R_\oplus`.

Where the slabs are cut matters, and two different things are going on.

Between shells the density *jumps*, so the chord is first split at every
crossing of a PREM shell boundary: no amount of subdivision recovers a
discontinuity that straddles a slab.  This gives a set of **chord
segments**.  Note a segment is not a shell --- a chord enters and leaves
each shell it reaches, so a diametric chord has 19 segments across 10
shells.

Within a shell the density *varies smoothly*, since PREM gives it as a
polynomial in :math:`x = r/R_\oplus` rather than a constant, and a
segment can be long: crossing the mantle, the density changes by 21%
over a single 2200 km segment.  Each segment is therefore divided
further into ``n_slabs_per_segment`` equal sub-slabs, with the density
taken at the midpoint of each.

Midpoint sampling is second-order accurate, so the result converges to
the continuous answer as the sub-slabs are refined; the routines take
that number as an argument so a caller can watch it converge rather than
trust it.  Sampling at the midpoint rather than an end also matters for
the segment that straddles the closest approach: it enters and exits at
the same radius, so its two ends have identical density while the
interior differs from them by 2.5%.

Units follow the rest of the library: energies in eV, baselines in
eV\ :sup:`-1`, potentials in eV.  The exceptions are the geometry
routines, which work in km because that is how the Earth is described,
and `density_prem`, which returns g cm\ :sup:`-3` because that is how
PREM is stated.  `matter_potential` is the bridge between the two.

The Earth is treated as a sphere of radius `globaldefs.EARTH_RADIUS`.

Routine listings
----------------

    * dms_to_decimal - Degrees, minutes, seconds to decimal degrees
    * coordinates_of_named_location - Coordinates of a named site
    * density_prem - PREM density at a radius
    * matter_potential - Charged-current potential from a density
    * matter_potential_nc - Neutral-current potential, for a sterile state
    * distance_traveled_inside_earth - Chord length for a given costhz
    * earth_radial_distance_from_depth - Radius at a point on the chord
    * prem_layer_edges_along_chord - Where a chord crosses PREM shells
    * chord_length_inside_earth - Chord between two surface locations
    * costhz_between_points_on_surface - Its zenith angle
    * earth_slabs - Slab widths and densities along a chord
    * slabs_for_tolerance - Subdivision needed for a stated tolerance
    * probabilities_2nu_earth - Two-flavor probabilities across the Earth
    * probabilities_3nu_earth - Three-flavor probabilities across the Earth
    * probabilities_4nu_earth - Four-flavor probabilities across the Earth
    * probabilities_2nu_between_locations - Between two named sites
    * probabilities_3nu_between_locations - Between two named sites
    * probabilities_4nu_between_locations - Between two named sites
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

__all__ = ['LOC_COORDS_DMS', 'PREM_BOUNDARIES',
           # The chunking constants.  Each carries an autodoc docstring and
           # `MAX_CHUNK_BYTES` is named in the changelog as the knob to
           # retune, but none of them was listed here, so `automodule`
           # skipped all five and the documentation was written for a page
           # it never reached.
           'CHUNK_BYTES_FALLBACK', 'CHUNK_BYTES_MIN', 'CHUNK_BYTES_MAX',
           'MIN_CHUNK_ENERGIES', 'MAX_CHUNK_BYTES',
           'dms_to_decimal', 'coordinates_of_named_location',
           'density_prem', 'matter_potential', 'matter_potential_nc',
           'distance_traveled_inside_earth',
           'earth_radial_distance_from_depth',
           'prem_layer_edges_along_chord', 'chord_length_inside_earth',
           'costhz_between_points_on_surface', 'earth_slabs',
           'slabs_for_tolerance',
           'probabilities_2nu_earth', 'probabilities_3nu_earth',
           'probabilities_4nu_earth',
           'probabilities_2nu_between_locations',
           'probabilities_3nu_between_locations',
           'probabilities_4nu_between_locations']

import os
from functools import lru_cache
from typing import Optional, Tuple, Union

import numpy as np

import fastkernels
import globaldefs as gd
import hamiltonians2nu
import hamiltonians3nu
import hamiltonians4nu
import oscprob4nu
import slabs


LOC_COORDS_DMS = {
    'baikal':      {'lat': (51, 45, 54),    'lon': (104, 24, 54)},
    'cern':        {'lat': (46, 14, 1.80),  'lon': (6, 3, 11.40)},
    'desy':        {'lat': (53, 34, 19.79), 'lon': (9, 52, 27.59)},
    'ess':         {'lat': (55, 44, 6),     'lon': (13, 15, 5.04)},
    'fermilab':    {'lat': (41, 49, 55),    'lon': (-88, 15, 26)},
    'gran_sasso':  {'lat': (42, 25, 15.8),  'lon': (13, 30, 58.43)},
    'homestake':   {'lat': (44, 21, 5.76),  'lon': (-103, 45, 4.68)},
    'kamioka':     {'lat': (36, 25, 50.05), 'lon': (137, 18, 41.15)},
    'km3net_arca': {'lat': (36, 16, 0),     'lon': (16, 6, 0)},
    'km3net_orca': {'lat': (42, 48, 0),     'lon': (6, 2, 0)},
    'north_pole':  {'lat': (90, 0, 0),      'lon': (0, 0, 0)},
    'pyhaasalmi':  {'lat': (63, 39, 31),    'lon': (26, 2, 28)},
    'snolab':      {'lat': (46, 28, 18),    'lon': (-81, 11, 12)},
    'south_pole':  {'lat': (-90, 0, 0),     'lon': (0, 0, 0)},
    'tokai':       {'lat': (36, 27, 59),    'lon': (140, 36, 24)},
}
r"""dict: Predefined locations, in ISO 6709 convention.

North latitudes and East longitudes are positive; South and West are
negative.  Each entry gives ``lat`` and ``lon`` as (degree, minute,
second) tuples.  The same set of sites as the sibling Magnus package,
so a trajectory named in one can be reproduced in the other.
"""

PREM_BOUNDARIES = np.array([1221.5, 3480.0, 5701.0, 5771.0, 5971.0, 6151.0,
                            6346.6, 6356.0, 6368.0])
r"""numpy.ndarray: Outer radius of each PREM shell but the last.

The last shell ends at `globaldefs.EARTH_RADIUS`.  Units: [km].
"""

_PREM_COEFFS = np.array([
    [13.0885,  0.0,    -8.8381,  0.0],
    [12.5815, -1.2638, -3.6426, -5.5281],
    [7.9565,  -6.4761,  5.5283, -3.0807],
    [5.3197,  -1.4836,  0.0,     0.0],
    [11.2494, -8.0298,  0.0,     0.0],
    [7.1089,  -3.8045,  0.0,     0.0],
    [2.6910,   0.6924,  0.0,     0.0],
    [2.900,    0.0,     0.0,     0.0],
    [2.600,    0.0,     0.0,     0.0],
    [1.020,    0.0,     0.0,     0.0],
])


def dms_to_decimal(
    degrees: Union[int, float],
    minutes: Union[int, float],
    seconds: Union[int, float]
) -> float:
    r"""Returns a (degree, minute, second) coordinate in decimal degrees.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    degrees : int or float
        Degree part of the coordinate.  Carries the sign: a location at
        5 degrees South is ``(-5, ...)``.  For a coordinate between zero
        and one degree South or West, where the degree part cannot carry
        a sign, negate the minutes instead: 0 deg 30' S is
        ``(0, -30, 0)``.
    minutes : int or float
        Minute part of the coordinate.  Normally positive; a negative
        value flips the sign of the whole coordinate, which is the only
        way to express a southern or western coordinate smaller than one
        degree.
    seconds : int or float
        Second part of the coordinate, taken as positive.

    Returns
    -------
    float
        The coordinate in decimal degrees.

    Examples
    --------
    .. jupyter-execute::

        import earth

        print('%.6f' % earth.dms_to_decimal(36, 25, 50.05))
    """
    magnitude = abs(degrees) + abs(minutes)/60.0 + seconds/3600.0

    # The sign lives on the degree part, so a negative latitude must not
    # have its minutes and seconds added back the other way.  Between
    # zero and one degree the degree part is zero and cannot carry a
    # sign, so the minutes are allowed to carry it instead --- without
    # that, 0 deg 30' South is inexpressible, and silently comes back
    # North.
    negative = degrees < 0 or (degrees == 0 and minutes < 0)

    return -magnitude if negative else magnitude


def coordinates_of_named_location(
    loc_name: str
) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
    r"""Returns the coordinates of a predefined location.

    Looks ``loc_name`` up in `LOC_COORDS_DMS`, case-insensitively and
    treating spaces as underscores.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    loc_name : str
        Name of the location, e.g. ``'kamioka'`` or ``'south pole'``.

    Returns
    -------
    tuple of tuple of float
        The latitude and longitude, each as (degree, minute, second).

    Raises
    ------
    ValueError
        If the name is not one of the predefined locations.

    Examples
    --------
    .. jupyter-execute::

        import earth

        lat, lon = earth.coordinates_of_named_location('South Pole')
        print(lat, lon)
    """
    key = loc_name.lower().replace(' ', '_')
    try:
        entry = LOC_COORDS_DMS[key]
    except KeyError:
        # `from None` keeps the KeyError out of the traceback: it is an
        # implementation detail of the lookup, and chaining it only
        # buries the message that explains what to do.
        raise ValueError(
            'coordinates_of_named_location: %r is not a predefined '
            'location; the available names, in earth.LOC_COORDS_DMS, are: '
            '%s' % (loc_name, ', '.join(sorted(LOC_COORDS_DMS)))) from None

    return entry['lat'], entry['lon']


def density_prem(
    r: Union[int, float, list, np.ndarray],
    tol: float = 1.e-8
) -> Union[float, np.ndarray]:
    r"""Returns the matter density inside the Earth, according to PREM.

    Evaluates the Preliminary Reference Earth Model
    :cite:`Dziewonski:1981xy` at a radial distance measured from the
    centre of the Earth.  Accepts a single radius or an array of radii,
    evaluated in one vectorised pass.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    r : int, float, list or numpy.ndarray
        Radial distance(s) from the centre of the Earth, in units of km.
    tol : float, optional
        Relative tolerance by which a radius may exceed
        `globaldefs.EARTH_RADIUS` before a ValueError is raised.  Radii
        within the tolerance are clamped onto the surface, which is what
        makes a chord endpoint computed in floating point safe to pass
        in.  Default: 1e-8.

    Returns
    -------
    float or numpy.ndarray
        The matter density, in units of g cm\ :sup:`-3`.

    Raises
    ------
    ValueError
        If any radius exceeds `globaldefs.EARTH_RADIUS` by more than the
        relative tolerance, or if any radius is negative.

    Examples
    --------
    .. jupyter-execute::

        import earth

        print('%.4f' % earth.density_prem(0.0))
        print('%.4f' % earth.density_prem(6371.0))
    """
    scalar_input = (np.ndim(r) == 0)
    r = np.asarray(r, dtype=float)

    if np.any(r < 0.0):
        raise ValueError('density_prem: radial distance cannot be negative')

    x = r/gd.EARTH_RADIUS
    if np.any(x - 1.0 > tol):
        raise ValueError(
            'density_prem: radial distance cannot exceed '
            'globaldefs.EARTH_RADIUS = %s km by more than the relative '
            'tolerance tol = %s' % (gd.EARTH_RADIUS, tol))

    # Clamp radii within tolerance of the surface onto the surface
    r = np.minimum(r, gd.EARTH_RADIUS)
    x = np.minimum(x, 1.0)

    # Look up the PREM shell of each radius --- side='left' reproduces the
    # right-closed bins of the piecewise definition, e.g. r <= 1221.5 ---
    # and evaluate the density polynomial by Horner's rule.
    c = _PREM_COEFFS[np.searchsorted(PREM_BOUNDARIES, r, side='left')]
    density = c[..., 0] + x*(c[..., 1] + x*(c[..., 2] + x*c[..., 3]))

    return float(density) if scalar_input else density


def matter_potential(
    density: Union[int, float, list, np.ndarray],
    electron_fraction: float = gd.ELECTRON_FRACTION_EARTH_CRUST
) -> Union[float, np.ndarray]:
    r"""Returns the charged-current matter potential for a density.

    Returns :math:`V_{CC} = \sqrt{2} G_F n_e`, the potential that
    `hamiltonians3nu.hamiltonian_3nu_matter` and its two-flavor
    counterpart expect.  It is positive for neutrinos; pass its negative
    for antineutrinos.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    density : int, float, list or numpy.ndarray
        Matter density, in units of g cm\ :sup:`-3`.
    electron_fraction : float, optional
        Electrons per nucleon.  Default:
        `globaldefs.ELECTRON_FRACTION_EARTH_CRUST`, which is 0.5 and is
        a good approximation everywhere in the Earth.

    Returns
    -------
    float or numpy.ndarray
        The potential :math:`V_{CC}`, in units of eV.

    Examples
    --------
    .. jupyter-execute::

        import earth

        print('%.4e' % earth.matter_potential(3.0))
    """
    scalar_input = (np.ndim(density) == 0)
    density = np.asarray(density, dtype=float)

    # Electron number density in eV^3, by the same route as
    # globaldefs.NUM_DENSITY_E_EARTH_CRUST
    num_density_e = (density*gd.CONV_G_TO_EV
                     / ((gd.MASS_PROTON+gd.MASS_NEUTRON)/2.0)
                     * electron_fraction
                     / pow(gd.CONV_CM_TO_INV_EV, 3.0))
    potential = np.sqrt(2.0)*gd.GF*num_density_e

    return float(potential) if scalar_input else potential


def matter_potential_nc(
    density: Union[int, float, list, np.ndarray],
    neutron_fraction: Optional[float] = None,
    electron_fraction: float = gd.ELECTRON_FRACTION_EARTH_CRUST
) -> Union[float, np.ndarray]:
    r"""Returns the neutral-current matter potential for a density.

    Returns :math:`V_{NC} = -G_F n_n/\sqrt{2}`, which is **negative**
    for neutrinos.  It is felt equally by all three active flavors, so
    at two and three flavors it is proportional to the identity and
    drops out of the probabilities entirely --- which is why
    `matter_potential` alone serves them.

    It does not drop out once a sterile state is present, because the
    sterile state does not feel it.  Removing it from all four states
    costs only a global phase and leaves :math:`-V_{NC}` on the sterile
    entry; see :func:`hamiltonians4nu.hamiltonian_4nu_matter`.

    .. versionadded:: 1.11.0

    Parameters
    ----------
    density : int, float, list or numpy.ndarray
        Matter density, in units of g cm\ :sup:`-3`.
    neutron_fraction : float, optional
        Neutrons per nucleon.  Default: ``1 - electron_fraction``, the
        isoscalar value, since a nucleon is either a proton --- matched
        by an electron --- or a neutron.
    electron_fraction : float, optional
        Electrons per nucleon, used only to derive `neutron_fraction`
        when that is not given.  Default:
        `globaldefs.ELECTRON_FRACTION_EARTH_CRUST`.

    Returns
    -------
    float or numpy.ndarray
        The potential :math:`V_{NC}`, in units of eV.  Negative for
        neutrinos.

    Examples
    --------
    .. jupyter-execute::

        import earth

        print('%.4e' % earth.matter_potential_nc(3.0))
    """
    if neutron_fraction is None:
        neutron_fraction = 1.0 - electron_fraction

    scalar_input = (np.ndim(density) == 0)
    density = np.asarray(density, dtype=float)

    # Neutron number density in eV^3, by the same route as the electron
    # one in `matter_potential`
    num_density_n = (density*gd.CONV_G_TO_EV
                     / ((gd.MASS_PROTON+gd.MASS_NEUTRON)/2.0)
                     * neutron_fraction
                     / pow(gd.CONV_CM_TO_INV_EV, 3.0))
    potential = -gd.GF*num_density_n/np.sqrt(2.0)

    return float(potential) if scalar_input else potential


def _check_costhz(costhz: Union[int, float], caller: str) -> None:
    r"""Raises unless `costhz` is a cosine.

    Every geometry routine here funnels through
    `distance_traveled_inside_earth`, whose chord length is
    :math:`-2 R \cos\theta_z`.  That expression is happy to be handed a
    number outside :math:`[-1, 1]` and returns a chord longer than the
    Earth --- 19 113 km for ``costhz = -1.5``, against a diameter of
    12 742 km --- which then acquires a full set of plausible-looking
    slabs and densities.  Nothing further downstream notices, so it is
    caught here.

    Parameters
    ----------
    costhz : int or float
        Cosine of the zenith angle, which must lie in :math:`[-1, 1]`.
    caller : str
        Name of the calling routine, used in the error message.

    Returns
    -------
    None
        Nothing; the routine either returns or raises.

    Raises
    ------
    ValueError
        If ``costhz`` lies outside :math:`[-1, 1]`, or is not a number.

    .. versionadded:: 1.11.0
    """
    if not -1.0 <= costhz <= 1.0:
        raise ValueError(
            '%s: costhz is the cosine of the zenith angle and so must lie '
            'in [-1, 1]; got %r.  A value outside that range describes no '
            'direction, and would give a chord longer than the Earth.'
            % (caller, costhz))


def distance_traveled_inside_earth(costhz: Union[int, float]) -> float:
    r"""Returns the chord length through the Earth for a given direction.

    The neutrino is assumed to arrive at a detector on the surface, not
    underground, so the distance is zero for any *down*-going direction,
    ``costhz >= 0``, which reaches the detector from the sky without
    entering the Earth at all.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    costhz : int or float
        Cosine of the zenith angle of the neutrino direction.
        ``costhz = -1`` is straight up through the centre of the Earth.

    Returns
    -------
    float
        The chord length, in units of km.

    Raises
    ------
    ValueError
        If ``costhz`` lies outside :math:`[-1, 1]`, where it describes
        no direction.

    Examples
    --------
    .. jupyter-execute::

        import earth

        print('%.1f' % earth.distance_traveled_inside_earth(-1.0))
        print('%.1f' % earth.distance_traveled_inside_earth(0.5))
    """
    _check_costhz(costhz, 'distance_traveled_inside_earth')

    return 0.0 if costhz >= 0.0 else -2.0*gd.EARTH_RADIUS*costhz


def earth_radial_distance_from_depth(
    costhz: Union[int, float],
    l: Union[int, float, list, np.ndarray],
    tol: float = 1.e-8
) -> Union[float, np.ndarray]:
    r"""Returns the radius at a point along a chord through the Earth.

    A neutrino with direction ``costhz`` travels from ``l = 0`` at its
    point of entry to ``l =`` `distance_traveled_inside_earth`
    (``costhz``) at the detector.  This returns its distance from the
    centre of the Earth at a given ``l``.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    costhz : int or float
        Cosine of the zenith angle of the neutrino direction.
    l : int, float, list or numpy.ndarray
        Distance(s) along the chord from the point of entry, in units of
        km.
    tol : float, optional
        Absolute tolerance, in km, by which ``l`` may exceed the chord
        length before a ValueError is raised; values within the
        tolerance are clamped onto the exit point.  Default: 1e-8.

    Returns
    -------
    float or numpy.ndarray
        The radial distance from the centre of the Earth, in units of
        km.

    Raises
    ------
    ValueError
        If any ``l`` is negative or exceeds the chord length for this
        ``costhz`` by more than the tolerance.

    Examples
    --------
    .. jupyter-execute::

        import earth

        print('%.1f' % earth.earth_radial_distance_from_depth(-1.0, 6371.0))
    """
    scalar_input = (np.ndim(l) == 0)
    l = np.asarray(l, dtype=float)

    d = distance_traveled_inside_earth(costhz)

    if np.any(l < 0.0):
        raise ValueError('earth_radial_distance_from_depth: l cannot be '
                         'negative')
    if np.any(l - d > tol):
        raise ValueError(
            'earth_radial_distance_from_depth: l cannot exceed the distance '
            'traveled inside the Earth, %s km, for costhz = %s'
            % (d, costhz))

    # Clamp values within tolerance of the exit point onto the exit point
    l = np.minimum(l, d)

    u = d - l
    r2 = (gd.EARTH_RADIUS*gd.EARTH_RADIUS + u*u
          + 2.0*gd.EARTH_RADIUS*u*costhz)

    # r2 is a squared distance and cannot be negative for an `l` inside
    # the chord, which the clamping above guarantees; round-off can still
    # take it a few ulp below zero at the endpoints, where it vanishes.
    # Clipping at zero absorbs that without also hiding a genuinely
    # negative value the way `np.abs` would.
    r = np.sqrt(np.maximum(r2, 0.0))

    return float(r) if scalar_input else r


def prem_layer_edges_along_chord(costhz: Union[int, float]) -> np.ndarray:
    r"""Returns where a chord through the Earth crosses PREM boundaries.

    The density is discontinuous across a PREM shell boundary, so a slab
    that straddles one cannot represent it however finely the rest of
    the chord is divided.  These positions are therefore mandatory slab
    edges, and `earth_slabs` uses them as such.

    The crossings solve :math:`r(l) = r_b` for each boundary radius
    :math:`r_b`, a quadratic in :math:`l`: with :math:`u = d - l` and
    :math:`d = -2 R \cos\theta_z`,

    .. math::

       u^2 + 2 R \cos\theta_z\, u + \left(R^2 - r_b^2\right) = 0 .

    .. versionadded:: 1.8.0

    Parameters
    ----------
    costhz : int or float
        Cosine of the zenith angle of the neutrino direction.  Crossings
        exist only for ``costhz < 0``.

    Returns
    -------
    numpy.ndarray
        Sorted crossing positions along the chord, in units of km, each
        strictly inside ``(0, d)``.  Empty if the chord crosses no
        boundary.

    Examples
    --------
    .. jupyter-execute::

        import earth

        print(len(earth.prem_layer_edges_along_chord(-1.0)))
        print(len(earth.prem_layer_edges_along_chord(0.5)))
    """
    if costhz >= 0.0:
        return np.array([])

    radius = gd.EARTH_RADIUS
    d = -2.0*radius*costhz                            # chord length [km]
    rmin2 = radius*radius*(1.0 - costhz*costhz)       # closest approach^2

    crossings = []
    for rb in PREM_BOUNDARIES:
        disc = rb*rb - rmin2
        if disc <= 0.0:                               # never reaches this depth
            continue
        s = np.sqrt(disc)
        for u in (-radius*costhz - s, -radius*costhz + s):
            # The two roots are d/2 +/- s, and both lie strictly inside
            # (0, d) exactly when r_b < R, which holds for every entry of
            # PREM_BOUNDARIES (the outermost is 6368 km against a radius of
            # 6371).  The guard is therefore never false as things stand,
            # and is kept only so that adding a boundary at the surface
            # would degrade gracefully rather than emit a spurious edge.
            if 0.0 < u < d:                       # pragma: no branch
                crossings.append(d - u)

    # `np.unique` sorts, so the crossings need not be sorted first
    return np.unique(np.array(crossings))


def chord_length_inside_earth(
    lat1_dms: Tuple[float, float, float],
    lon1_dms: Tuple[float, float, float],
    lat2_dms: Tuple[float, float, float],
    lon2_dms: Tuple[float, float, float]
) -> float:
    r"""Returns the straight-line distance between two surface locations.

    Computes the chord --- the straight line through the Earth's
    interior, not the great-circle arc over its surface --- between two
    points on a spherical Earth, via the haversine formula for the
    central angle.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    lat1_dms : tuple of float
        Latitude of the first location, as (degree, minute, second).
    lon1_dms : tuple of float
        Longitude of the first location, as (degree, minute, second).
    lat2_dms : tuple of float
        Latitude of the second location, as (degree, minute, second).
    lon2_dms : tuple of float
        Longitude of the second location, as (degree, minute, second).

    Returns
    -------
    float
        The chord length, in units of km.

    Examples
    --------
    .. jupyter-execute::

        import earth

        lat1, lon1 = earth.coordinates_of_named_location('fermilab')
        lat2, lon2 = earth.coordinates_of_named_location('homestake')
        print('%.1f' % earth.chord_length_inside_earth(lat1, lon1, lat2, lon2))
    """
    lat1 = np.radians(dms_to_decimal(*lat1_dms))
    lon1 = np.radians(dms_to_decimal(*lon1_dms))
    lat2 = np.radians(dms_to_decimal(*lat2_dms))
    lon2 = np.radians(dms_to_decimal(*lon2_dms))

    # Haversine formula for the central angle
    a = (np.sin((lat2-lat1)/2.0)**2.0
         + np.cos(lat1)*np.cos(lat2)*np.sin((lon2-lon1)/2.0)**2.0)
    central_angle = 2.0*np.arctan2(np.sqrt(a), np.sqrt(1.0-a))

    return float(2.0*gd.EARTH_RADIUS*np.sin(central_angle/2.0))


def costhz_between_points_on_surface(
    lat1_dms: Tuple[float, float, float],
    lon1_dms: Tuple[float, float, float],
    lat2_dms: Tuple[float, float, float],
    lon2_dms: Tuple[float, float, float]
) -> float:
    r"""Returns the zenith angle of the chord between two locations.

    The cosine of the zenith angle at which a neutrino must travel to
    reach the second location from the first through the Earth's
    interior.  Both are assumed to be on the surface, so the result is
    never positive.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    lat1_dms : tuple of float
        Latitude of the first location, as (degree, minute, second).
    lon1_dms : tuple of float
        Longitude of the first location, as (degree, minute, second).
    lat2_dms : tuple of float
        Latitude of the second location, as (degree, minute, second).
    lon2_dms : tuple of float
        Longitude of the second location, as (degree, minute, second).

    Returns
    -------
    float
        Cosine of the zenith angle of the connecting chord.

    Examples
    --------
    .. jupyter-execute::

        import earth

        lat1, lon1 = earth.coordinates_of_named_location('cern')
        lat2, lon2 = earth.coordinates_of_named_location('gran_sasso')
        print('%.6f' % earth.costhz_between_points_on_surface(lat1, lon1,
                                                        lat2, lon2))
    """
    chord = chord_length_inside_earth(lat1_dms, lon1_dms, lat2_dms, lon2_dms)

    return -0.5*chord/gd.EARTH_RADIUS


@lru_cache(maxsize=256)
def _earth_slabs_cached(
    costhz: float,
    n_slabs_per_segment: int
) -> Tuple[np.ndarray, np.ndarray]:
    r"""Returns the cached slab widths and densities along a chord.

    The chord geometry depends on the zenith angle alone --- not on the
    neutrino energy, nor on the flavor count, nor on the Hamiltonian.  A
    scan over energy at fixed `costhz` therefore recomputed the identical
    PREM crossing for every point, which measured 176 of the 396
    microseconds an Earth probability took: the single largest cost in
    the call, and entirely redundant.

    The arrays handed back are marked read-only, so that the copy the
    public `earth_slabs` makes is the only writable one and an accidental
    write in here raises rather than silently poisoning the cache for
    every later caller.

    Validation lives in this function rather than in the wrapper because
    `_earth_hamiltonians` calls it directly, and a check that only the
    public path performs is a check with a way around it.
    """
    _check_costhz(costhz, 'earth_slabs')
    if costhz >= 0.0:
        raise ValueError(
            'earth_slabs: costhz must be negative for the neutrino to cross '
            'the Earth; got %s' % costhz)
    if n_slabs_per_segment < 1:
        raise ValueError('earth_slabs: n_slabs_per_segment must be at least 1')

    d = distance_traveled_inside_earth(costhz)
    edges = np.concatenate(([0.0], prem_layer_edges_along_chord(costhz), [d]))

    widths = []
    midpoints = []
    for start, end in zip(edges[:-1], edges[1:]):
        # `prem_layer_edges_along_chord` returns strictly increasing values
        # strictly inside (0, d), and np.unique has already removed exact
        # duplicates, so `edges` is strictly increasing and this cannot
        # fire.  It is kept against a crossing that rounds onto an endpoint,
        # which would otherwise produce a zero-width slab.
        if end <= start:                          # pragma: no cover
            continue
        sub = np.linspace(start, end, n_slabs_per_segment+1)
        widths.append(np.diff(sub))
        midpoints.append((sub[:-1]+sub[1:])/2.0)

    widths = np.concatenate(widths)
    midpoints = np.concatenate(midpoints)
    densities = density_prem(
        earth_radial_distance_from_depth(costhz, midpoints))

    # A chord through a spherically symmetric Earth meets every radius
    # twice, symmetrically about its closest approach, so both of these
    # are palindromes.  The densities come out exactly palindromic
    # already; the widths do not, because each segment is cut by its own
    # `linspace` and the two halves round differently --- by about
    # 1e-12 km on a 100 km slab.  Averaging each element with its mirror
    # makes both exact, since floating-point addition is commutative and
    # so the two ends of a pair get bitwise identical results.
    #
    # This is not housekeeping.  `_palindromic` decides on exact equality
    # whether a chord can be composed at half cost, and a difference in
    # the last bit is the difference between taking that path and not.
    widths = (widths + widths[::-1])/2.0
    densities = (densities + densities[::-1])/2.0

    widths.flags.writeable = False
    densities.flags.writeable = False

    return widths, densities


def earth_slabs(
    costhz: Union[int, float],
    n_slabs_per_segment: int = 8
) -> Tuple[np.ndarray, np.ndarray]:
    r"""Returns the slab widths and densities along a chord.

    Cuts the chord at every PREM shell boundary it crosses, divides each
    resulting segment into ``n_slabs_per_segment`` equal sub-slabs, and
    evaluates the density at the midpoint of each.  The boundary cuts
    are what make the discretisation converge quickly: they keep every
    slab inside a single shell, where the density is smooth.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    costhz : int or float
        Cosine of the zenith angle of the neutrino direction.  Must be
        negative, so that the neutrino crosses the Earth at all.
    n_slabs_per_segment : int, optional
        Number of equal sub-slabs per chord segment.  A segment runs
        between consecutive PREM boundary crossings; a chord crosses most
        shells twice, so there are more segments than shells.  Default: 8.

    Returns
    -------
    tuple of numpy.ndarray
        The slab widths, in units of km, and the density in each slab,
        in units of g cm\ :sup:`-3`, ordered along the trajectory.

    Raises
    ------
    ValueError
        If ``costhz >= 0``, so that there is no path through the Earth,
        or if ``n_slabs_per_segment`` is not positive.

    Examples
    --------
    .. jupyter-execute::

        import earth

        widths, densities = earth.earth_slabs(-1.0, n_slabs_per_segment=2)
        print(len(widths), '%.1f' % sum(widths))
    """
    widths, densities = _earth_slabs_cached(float(costhz),
                                            int(n_slabs_per_segment))

    # A copy, because the cached arrays are shared with every other
    # caller and this one is entitled to modify what it is given.
    return widths.copy(), densities.copy()


def _earth_hamiltonians(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float, list, np.ndarray],
    costhz: Union[int, float, list, np.ndarray],
    n_slabs_per_segment: int,
    electron_fraction: float,
    n_flavors: int
) -> Tuple[np.ndarray, np.ndarray]:
    r"""Returns the per-slab Hamiltonians and widths for a chord.

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent vacuum Hamiltonian.
    energy : int, float or numpy.ndarray
        Neutrino energy, in units of eV, or an array of energies, in
        which case one chord of Hamiltonians is built per energy.
    costhz : int or float
        Cosine of the zenith angle.
    n_slabs_per_segment : int
        Number of equal sub-slabs per chord segment.  A segment runs
        between consecutive PREM boundary crossings; a chord crosses most
        shells twice, so there are more segments than shells.
    electron_fraction : float
        Electrons per nucleon.
    n_flavors : int
        Number of neutrino flavors, 2, 3, or 4.

    Returns
    -------
    tuple of numpy.ndarray
        The Hamiltonians, of shape ``(n, n_flavors, n_flavors)`` for a
        scalar energy and ``(..., n, n_flavors, n_flavors)`` for an
        array of them, and the slab widths in units of eV\ :sup:`-1`.
    """
    # The cached arrays directly: this routine only reads them, and
    # the multiplication below allocates its own result.
    widths_km, densities = _earth_slabs_cached(float(costhz),
                                               int(n_slabs_per_segment))
    potentials = matter_potential(densities, electron_fraction)

    # The slab axis is the last one the potentials carry, so the energy
    # gains a trailing axis of its own to broadcast against it: a scalar
    # energy still yields one chord, and an array of energies yields one
    # chord each, through the same expression.  The matter potential
    # depends on the geometry alone, so a scan over energy builds it once
    # here rather than once per energy.
    energy = np.asarray(energy, dtype=float)[..., None]

    # The Hamiltonian builders take an array of potentials and return one
    # Hamiltonian per entry, so the whole chord is built in one call.
    if n_flavors == 2:
        h = hamiltonians2nu.hamiltonian_2nu_matter(
            h_vacuum_energy_independent, energy, potentials)
    elif n_flavors == 3:
        h = hamiltonians3nu.hamiltonian_3nu_matter(
            h_vacuum_energy_independent, energy, potentials)
    else:
        # A sterile state does not feel the neutral-current potential,
        # which therefore no longer cancels and has to be built too
        h = hamiltonians4nu.hamiltonian_4nu_matter(
            h_vacuum_energy_independent, energy, potentials,
            matter_potential_nc(densities,
                                electron_fraction=electron_fraction))

    return h, widths_km*gd.CONV_KM_TO_INV_EV


CHUNK_BYTES_FALLBACK = 16*1024*1024
r"""int: Module-level constant.

Chunk size assumed when the cache size cannot be read, in bytes.  A
middling last-level cache for a machine of the era; see
`MAX_CHUNK_BYTES`.

.. versionadded:: 1.12.0
"""

CHUNK_BYTES_MIN = 4*1024*1024
r"""int: Module-level constant.

Smallest chunk the detected cache size may produce, in bytes.  Below
roughly this the per-chunk overhead starts to cost more than the cache
residency buys.

.. versionadded:: 1.12.0
"""

CHUNK_BYTES_MAX = 64*1024*1024
r"""int: Module-level constant.

Largest chunk the detected cache size may produce, in bytes.  A server
with a very large last-level cache should not turn that into a very
large allocation.

.. versionadded:: 1.12.0
"""

MIN_CHUNK_ENERGIES = 32
r"""int: Module-level constant.

Fewest energies a chunk may hold, whatever the byte budget says.  A
four-flavor chord with a thousand slabs costs a quarter of a megabyte
per energy, and cutting that into chunks of one or two would spend more
time re-entering the kernel than it saved.

.. versionadded:: 1.12.0
"""


_SYSFS_CACHE = '/sys/devices/system/cpu/cpu0/cache'


def _cache_bytes_from_sysconf() -> Optional[int]:
    r"""Returns the largest cache size :func:`os.sysconf` reports, or None.

    Some POSIX builds carry ``SC_LEVEL*_CACHE_SIZE`` names and some do
    not; the ones that do not raise `ValueError` when asked, which is
    why each name is tried separately rather than in one block.

    Returns
    -------
    int or None
        The largest cache size in bytes, or None if none is reported.
    """
    largest = None
    for name in ('SC_LEVEL4_CACHE_SIZE', 'SC_LEVEL3_CACHE_SIZE',
                 'SC_LEVEL2_CACHE_SIZE'):
        try:
            size = os.sysconf(name)
        except (ValueError, OSError, AttributeError):
            continue
        if isinstance(size, int) and size > 0:
            largest = size if largest is None else max(largest, size)

    return largest


def _cache_bytes_from_sysfs() -> Optional[int]:
    r"""Returns the largest cache size ``sysfs`` reports, or None.

    Linux only, and the most reliable of the three where it exists.

    Returns
    -------
    int or None
        The largest cache size in bytes, or None if it cannot be read.
    """
    scale = {'K': 1024, 'M': 1024*1024, 'G': 1024*1024*1024}
    largest = None

    for entry in os.listdir(_SYSFS_CACHE):
        try:
            with open(os.path.join(_SYSFS_CACHE, entry, 'size')) as handle:
                text = handle.read().strip()
        except OSError:
            continue
        if not text:
            continue
        unit = scale.get(text[-1].upper(), 1)
        try:
            size = int(text[:-1] if unit > 1 else text)*unit
        except ValueError:
            continue
        if largest is None or size > largest:
            largest = size

    return largest


def _cache_bytes_from_sysctl() -> Optional[int]:
    r"""Returns the largest cache size ``sysctl`` reports, or None.

    macOS, through :mod:`ctypes` rather than the ``sysctl`` command:
    spawning a subprocess while a module is still being imported is a
    great deal more to go wrong than this is worth.  Apple Silicon
    reports no ``hw.l3cachesize``, so the per-performance-level L2 is
    asked for as well, that being the largest cache those machines have.

    Returns
    -------
    int or None
        The largest cache size in bytes, or None if none is reported.
    """
    import ctypes
    import ctypes.util

    path = ctypes.util.find_library('c')
    if path is None:
        return None
    libc = ctypes.CDLL(path)
    if not hasattr(libc, 'sysctlbyname'):
        return None

    largest = None
    for name in (b'hw.l3cachesize', b'hw.perflevel0.l2cachesize',
                 b'hw.l2cachesize'):
        value = ctypes.c_uint64(0)
        length = ctypes.c_size_t(ctypes.sizeof(value))
        # `pointer` rather than `byref`: the extra object it builds costs
        # nothing at import time, and it is writable from Python, so the
        # parsing below can be tested on a machine with no `sysctl` at
        # all rather than only asserted about.
        if libc.sysctlbyname(name, ctypes.pointer(value),
                             ctypes.pointer(length), None, 0) != 0:
            continue
        if value.value > 0:
            largest = (value.value if largest is None
                       else max(largest, value.value))

    return largest


def _last_level_cache_bytes() -> Optional[int]:
    r"""Returns the size of the largest CPU cache, in bytes, or None.

    There is no portable way to ask, so three ways are tried and the
    first that answers wins: :func:`os.sysconf`, which some POSIX builds
    carry; Linux's ``sysfs``; and macOS's ``sysctl``.  Windows is not
    probed --- doing it means untested :mod:`ctypes` calls into
    ``kernel32`` running at import time on machines this was never run
    on, which is a poor trade for a hint.

    Every probe is wrapped, and broadly.  This value only tunes how a
    long scan is cut into pieces: being wrong costs some speed and
    nothing else, so no failure of it should be able to stop
    :mod:`earth` from importing.  `CHUNK_BYTES_FALLBACK` is deliberately
    a plausible answer rather than a degenerate one.

    Returns
    -------
    int or None
        The largest cache size in bytes, or None if nothing answered.
    """
    for probe in (_cache_bytes_from_sysconf, _cache_bytes_from_sysfs,
                  _cache_bytes_from_sysctl):
        try:
            size = probe()
        except Exception:                       # noqa: BLE001
            continue
        if size:
            return size

    return None


MAX_CHUNK_BYTES = min(CHUNK_BYTES_MAX,
                      max(CHUNK_BYTES_MIN,
                          _last_level_cache_bytes() or CHUNK_BYTES_FALLBACK))
r"""int: Module-level constant.

Rough ceiling on the Hamiltonian stack an array of energies may build at
once, in bytes.  Longer scans are evaluated in chunks of that size.
Set it to retune; nothing caches the value.

There are two reasons to chunk, and the second is the one that sets the
number.

The first is memory.  The stack is proportional to the scan length: a
hundred thousand energies across a 120-slab chord is 1.6 GB at three
flavors and, counting the traceless copy the expansion needs, nearly
6 GB at four.  That is an ordinary oscillogram, not an abusive input.

The second is that **the batched kernel is memory-bound, not
compute-bound**, and this is what makes the chunk size worth choosing
rather than merely bounding.  The stack is written by the Hamiltonian
builder and then streamed by the kernel, which does little arithmetic
per byte; if it fits in the last-level cache the second pass is nearly
free, and if it does not, every slab is fetched from memory.  Measured
on one machine, cost per probability was 7.9 microseconds with an
8 MB working set and 16.4 with a 540 MB one --- a factor of two paid for
nothing but traffic.  Interleaved against a 64 MB chunk, a cache-sized
one was 1.3x to 1.8x quicker across both chord lengths and all three
flavor counts, and never slower.

So the default is the detected last-level cache, clamped into
``[CHUNK_BYTES_MIN, CHUNK_BYTES_MAX]``, falling back to
`CHUNK_BYTES_FALLBACK` where it cannot be read.  That is a *guess at the
right order of magnitude*, not a tuned constant: it was measured on one
12 MB machine, the optimum is broad, and any value near the cache beats
one far above it.  A machine whose cache is shared between busy cores
may do better with less.  It is a plain module attribute for that
reason.

.. versionadded:: 1.12.0
"""


def _probabilities_earth_batch(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: np.ndarray,
    costhz: Union[int, float, list, np.ndarray],
    n_slabs_per_segment: int,
    electron_fraction: float,
    n_flavors: int,
    caller: str
) -> np.ndarray:
    r"""Returns the probabilities for an array of energies, in chunks.

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent vacuum Hamiltonian.
    energy : numpy.ndarray
        Array of neutrino energies, in units of eV.
    costhz : int or float
        Cosine of the zenith angle.
    n_slabs_per_segment : int
        Number of equal sub-slabs per chord segment.
    electron_fraction : float
        Electrons per nucleon.
    n_flavors : int
        Number of neutrino flavors, 2, 3, or 4.
    caller : str
        Name of the calling routine, used in error messages.

    Returns
    -------
    numpy.ndarray
        The probabilities, of shape ``(..., n_flavors*n_flavors)``.
    """
    energy = np.asarray(energy, dtype=float)

    # The chord is the same for every energy, so one look at the geometry
    # says how much a single energy costs and therefore how many fit.
    widths_km, densities = _earth_slabs_cached(float(costhz),
                                               int(n_slabs_per_segment))
    per_energy = widths_km.shape[0]*n_flavors*n_flavors*16
    chunk = max(MIN_CHUNK_ENERGIES, MAX_CHUNK_BYTES//per_energy)

    flat = energy.reshape(-1)

    # The fused kernels build each slab's Hamiltonian as they go, so the
    # stack that the chunking below exists to bound is never allocated at
    # all: what they read is one potential and one width per slab, shared
    # by every energy in the scan.  Nothing about the physics differs ---
    # this is the same midpoint Hamiltonian, and agrees bit for bit.
    if fastkernels.worthwhile_slabs(n_flavors,
                                    flat.shape[0]*widths_km.shape[0]):
        widths = widths_km*gd.CONV_KM_TO_INV_EV
        potentials = matter_potential(densities, electron_fraction)
        if n_flavors == 2:
            u = fastkernels.earth_chords_2nu_kernel(
                h_vacuum_energy_independent, flat, potentials, widths)
        elif n_flavors == 3:
            u = fastkernels.earth_chords_3nu_kernel(
                h_vacuum_energy_independent, flat, potentials, widths)
        else:
            u = fastkernels.earth_chords_4nu_kernel(
                h_vacuum_energy_independent, flat, potentials,
                matter_potential_nc(densities,
                                    electron_fraction=electron_fraction),
                widths, oscprob4nu.POLISH_ROOTS)
        # P_ab = |U_ba|^2, initial flavor varying slowest
        p = np.abs(np.swapaxes(u, -1, -2))**2.0
        return p.reshape(energy.shape + (n_flavors*n_flavors,))

    if flat.shape[0] <= chunk:
        h, widths = _earth_hamiltonians(h_vacuum_energy_independent, flat,
                                        costhz, n_slabs_per_segment,
                                        electron_fraction, n_flavors)
        out = slabs._probabilities_slabs_batch(h, widths, n_flavors, caller)
    else:
        out = np.empty((flat.shape[0], n_flavors*n_flavors), dtype=float)
        for start in range(0, flat.shape[0], chunk):
            piece = flat[start:start+chunk]
            h, widths = _earth_hamiltonians(
                h_vacuum_energy_independent, piece, costhz,
                n_slabs_per_segment, electron_fraction, n_flavors)
            out[start:start+chunk] = slabs._probabilities_slabs_batch(
                h, widths, n_flavors, caller)

    # The caller's batch shape, with the probabilities as the last axis
    return out.reshape(energy.shape + (n_flavors*n_flavors,))


def slabs_for_tolerance(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float, list, np.ndarray],
    costhz: Union[int, float, list, np.ndarray],
    n_flavors: int = 3,
    rtol: Optional[float] = None,
    atol: Optional[float] = None,
    n_start: int = 8,
    n_max: int = slabs.N_SLABS_MAX,
    electron_fraction: float = gd.ELECTRON_FRACTION_EARTH_CRUST
) -> int:
    r"""Returns the subdivision an Earth crossing needs for a tolerance.

    The discretisation error of an Earth crossing is strongly
    energy-dependent --- at the default eight sub-slabs per segment it
    spans more than an order of magnitude between 3 and 40 GeV --- so a
    fixed ``n_slabs_per_segment`` does not give a fixed accuracy.  This
    returns the subdivision that does, for the direction and energies
    asked about, by refining until the measured error meets the
    tolerance.

    Every probability the call returns must meet the tolerance, and when
    an array of energies is given the answer covers all of them: the
    subdivision is set by the worst-converging entry, which is what
    makes one number safe to reuse across a scan.

    .. versionadded:: 1.12.0

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent vacuum Hamiltonian, of the flavor count given
        by `n_flavors`.
    energy : int, float or array_like
        Neutrino energy, in units of eV, or an array of energies, in
        which case the answer meets the tolerance at every one of them.
    costhz : int, float or array_like
        Cosine of the zenith angle of the neutrino direction.  Must be
        negative.  May be an array, which is how an oscillogram is
        asked for: index the energies and the angles on different axes,
        as ``energy[None, :]`` against ``costhz[:, None]``, and they
        broadcast into a grid.  Each distinct angle costs one pass, so a
        grid is far cheaper than a loop over its points.
    n_flavors : int, optional
        Number of neutrino flavors, 2, 3, or 4.  Default: 3.
    rtol : float, optional
        Relative tolerance, taken against the probability itself.
        Default: None.
    atol : float, optional
        Absolute tolerance.  Default: None.  At least one of the two
        must be given, and when both are, the threshold is
        ``atol + rtol*abs(P)``, the convention of `numpy.isclose`.
    n_start : int, optional
        Coarsest subdivision to try, and the smallest that can be
        returned.  Default: 8, the default of the probability routines.
    n_max : int, optional
        Largest subdivision to try.  Default: `slabs.N_SLABS_MAX`.
    electron_fraction : float, optional
        Electrons per nucleon.  Default:
        `globaldefs.ELECTRON_FRACTION_EARTH_CRUST`.

    Returns
    -------
    int
        The subdivision to pass as ``n_slabs_per_segment``.

    Raises
    ------
    ValueError
        If the tolerances are invalid or absent, if ``costhz >= 0``, or
        if the tolerance is not met by ``n_max``.

    Notes
    -----
    This costs several Earth crossings, the whole refinement adding up
    to roughly twice the evaluation at the subdivision it returns.  It
    is meant to be called **once** for a scan and its answer passed to
    the calls in the loop, not called per probability --- passing
    ``rtol`` to `probabilities_3nu_earth` for every point of a scan
    repeats this search at every point.

    Examples
    --------
    .. jupyter-execute::

        import earth
        import globaldefs as gd
        import hamiltonians3nu

        H = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
            gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, gd.DCP_NO_BF,
            gd.D21_NO_BF, gd.D31_NO_BF)
        print(earth.slabs_for_tolerance(H, 1.0e10, -0.8, atol=1.0e-5))
    """
    def evaluate(n: int) -> np.ndarray:
        return np.asarray(_probabilities_earth(
            h_vacuum_energy_independent, energy, costhz, n,
            electron_fraction, n_flavors, 'slabs_for_tolerance'),
            dtype=float)

    n, _ = slabs._n_for_tolerance(evaluate, rtol, atol, n_start, n_max,
                                  'slabs_for_tolerance')

    return n


def _probabilities_earth(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float, list, np.ndarray],
    costhz: Union[int, float, list, np.ndarray],
    n_slabs_per_segment: int,
    electron_fraction: float,
    n_flavors: int,
    caller: str
) -> Union[Tuple[float, ...], np.ndarray]:
    r"""Returns the probabilities for one subdivision, scalar or batched.

    The common body of the three public entry points, so that the
    tolerance search can reach the same code they do.

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent vacuum Hamiltonian.
    energy : int, float or array_like
        Neutrino energy, in units of eV, or an array of energies.
    costhz : int or float
        Cosine of the zenith angle.
    n_slabs_per_segment : int
        Number of equal sub-slabs per chord segment.
    electron_fraction : float
        Electrons per nucleon.
    n_flavors : int
        Number of neutrino flavors, 2, 3, or 4.
    caller : str
        Name of the calling routine, used in error messages.

    Returns
    -------
    tuple of float or numpy.ndarray
        The probabilities, as a tuple for a scalar energy and an array
        of shape ``(..., n_flavors*n_flavors)`` for an array of them.
    """
    # An array of angles is a grid, however the energies are shaped
    if np.ndim(costhz) != 0:
        return _probabilities_earth_grid(
            h_vacuum_energy_independent, energy, costhz, n_slabs_per_segment,
            electron_fraction, n_flavors, caller)

    if np.ndim(energy) == 0:
        h, widths = _earth_hamiltonians(h_vacuum_energy_independent, energy,
                                        costhz, n_slabs_per_segment,
                                        electron_fraction, n_flavors)
        if n_flavors == 2:
            return slabs.probabilities_2nu_slabs(h, widths)
        if n_flavors == 3:
            return slabs.probabilities_3nu_slabs(h, widths)
        return slabs.probabilities_4nu_slabs(h, widths)

    return _probabilities_earth_batch(h_vacuum_energy_independent, energy,
                                      costhz, n_slabs_per_segment,
                                      electron_fraction, n_flavors, caller)


def _probabilities_earth_grid(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float, list, np.ndarray],
    costhz: Union[list, np.ndarray],
    n_slabs_per_segment: int,
    electron_fraction: float,
    n_flavors: int,
    caller: str
) -> np.ndarray:
    r"""Returns the probabilities over a grid of energies and angles.

    An oscillogram, in other words.  The energies and angles broadcast
    against each other in the usual way, so a grid is asked for as
    ``probabilities_3nu_earth(h, energies[None, :], costhz[:, None])``
    --- the same idiom `oscprob3nu.probabilities_3nu` uses for a stack
    of Hamiltonians against a stack of baselines.

    The angles are handled one at a time rather than all at once, and
    that is not a compromise: the chord geometry is what changes with
    the angle, so two angles share neither their slab widths nor their
    number of slabs, and there is nothing for a single kernel call to
    share.  What *is* shared is every energy at a given angle, which is
    the axis the fused kernel already spreads over, so the work goes
    from one call per grid point to one call per distinct angle.

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent vacuum Hamiltonian.
    energy : int, float or array_like
        Neutrino energies, in units of eV.
    costhz : array_like
        Cosines of the zenith angle, broadcastable against `energy`.
    n_slabs_per_segment : int
        Number of equal sub-slabs per chord segment.
    electron_fraction : float
        Electrons per nucleon.
    n_flavors : int
        Number of neutrino flavors, 2, 3, or 4.
    caller : str
        Name of the calling routine, used in error messages.

    Returns
    -------
    numpy.ndarray
        The probabilities, of shape ``(..., n_flavors*n_flavors)``,
        on the broadcast grid.

    Raises
    ------
    ValueError
        If the energies and angles do not broadcast together.
    """
    try:
        energy_b, costhz_b = np.broadcast_arrays(
            np.asarray(energy, dtype=float), np.asarray(costhz, dtype=float))
    except ValueError:
        raise ValueError(
            '%s: energy of shape %s and costhz of shape %s do not broadcast '
            'together; for a grid, index them on different axes, as '
            'energy[None, :] and costhz[:, None]'
            % (caller, (np.shape(energy),), (np.shape(costhz),))) from None

    flat_energy = energy_b.reshape(-1)
    flat_costhz = costhz_b.reshape(-1)
    out = np.empty((flat_energy.shape[0], n_flavors*n_flavors), dtype=float)

    # One pass per distinct angle, each carrying all of that angle's
    # energies.  np.unique also means a grid that repeats an angle --- a
    # broadcast one always does --- pays for its geometry once.
    for angle in np.unique(flat_costhz):
        at_angle = flat_costhz == angle
        out[at_angle] = _probabilities_earth_batch(
            h_vacuum_energy_independent, flat_energy[at_angle], float(angle),
            n_slabs_per_segment, electron_fraction, n_flavors, caller)

    return out.reshape(energy_b.shape + (n_flavors*n_flavors,))


def _probabilities_earth_tol(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float, list, np.ndarray],
    costhz: Union[int, float, list, np.ndarray],
    n_slabs_per_segment: int,
    electron_fraction: float,
    n_flavors: int,
    rtol: Optional[float],
    atol: Optional[float],
    n_max: int,
    return_n_slabs: bool,
    caller: str
) -> Union[Tuple[float, ...], np.ndarray, tuple]:
    r"""Returns the probabilities, refining first if a tolerance is set.

    With neither tolerance given this is `_probabilities_earth` and
    nothing else, so the default path is exactly what it was before
    tolerances existed.

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent vacuum Hamiltonian.
    energy : int, float or array_like
        Neutrino energy, in units of eV, or an array of energies.
    costhz : int or float
        Cosine of the zenith angle.
    n_slabs_per_segment : int
        Subdivision to use, or the coarsest to try when refining.
    electron_fraction : float
        Electrons per nucleon.
    n_flavors : int
        Number of neutrino flavors, 2, 3, or 4.
    rtol : float or None
        Relative tolerance, or None.
    atol : float or None
        Absolute tolerance, or None.
    n_max : int
        Largest subdivision to try when refining.
    return_n_slabs : bool
        Whether to return the subdivision used alongside the answer.
    caller : str
        Name of the calling routine, used in error messages.

    Returns
    -------
    tuple of float, numpy.ndarray, or tuple
        The probabilities, paired with the subdivision used when
        `return_n_slabs` is set.
    """
    if rtol is None and atol is None:
        p = _probabilities_earth(h_vacuum_energy_independent, energy, costhz,
                                 n_slabs_per_segment, electron_fraction,
                                 n_flavors, caller)
        n = n_slabs_per_segment
    else:
        def evaluate(n_try: int) -> np.ndarray:
            return np.asarray(_probabilities_earth(
                h_vacuum_energy_independent, energy, costhz, n_try,
                electron_fraction, n_flavors, caller), dtype=float)

        # The search already evaluated the answer at the subdivision it
        # settled on, so there is nothing left to compute here
        n, p = slabs._n_for_tolerance(evaluate, rtol, atol,
                                      n_slabs_per_segment, n_max, caller)
        if np.ndim(energy) == 0:
            # A scalar energy returns a tuple however the answer was
            # reached, which the search's array does not preserve
            p = tuple(float(x) for x in p)

    return (p, n) if return_n_slabs else p


def probabilities_2nu_earth(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float, list, np.ndarray],
    costhz: Union[int, float, list, np.ndarray],
    n_slabs_per_segment: int = 8,
    electron_fraction: float = gd.ELECTRON_FRACTION_EARTH_CRUST,
    rtol: Optional[float] = None,
    atol: Optional[float] = None,
    n_max: int = slabs.N_SLABS_MAX,
    return_n_slabs: bool = False
) -> Union[Tuple[float, float, float, float], np.ndarray]:
    r"""Returns the two-flavor probabilities across the Earth.

    Builds the PREM slabs along the chord for the given direction and
    propagates through them exactly, slab by slab.

    .. versionadded:: 1.8.0

    .. versionchanged:: 1.12.0
       Accepts an array of energies, returning one row of probabilities
       per energy.  A scalar energy returns exactly what it returned
       before.

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent two-flavor vacuum Hamiltonian, as returned by
        `hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent`.
    energy : int, float or array_like
        Neutrino energy, in units of eV, or an array of energies.  The
        whole array crosses the same chord, so the geometry and the
        matter potentials are built once for the scan rather than once
        per energy.
    costhz : int, float or array_like
        Cosine of the zenith angle of the neutrino direction.  Must be
        negative.  May be an array, which is how an oscillogram is
        asked for: index the energies and the angles on different axes,
        as ``energy[None, :]`` against ``costhz[:, None]``, and they
        broadcast into a grid.  Each distinct angle costs one pass, so a
        grid is far cheaper than a loop over its points.
    n_slabs_per_segment : int, optional
        Number of equal sub-slabs per chord segment.  A segment runs
        between consecutive PREM boundary crossings; a chord crosses most
        shells twice, so there are more segments than shells.  Default: 8.
    electron_fraction : float, optional
        Electrons per nucleon.  Default:
        `globaldefs.ELECTRON_FRACTION_EARTH_CRUST`.
    rtol : float, optional
        Relative tolerance on every returned probability.  Default:
        None, meaning ``n_slabs_per_segment`` is used as given.  When
        either tolerance is set, the chord is refined until the measured
        discretisation error meets it, starting from
        ``n_slabs_per_segment``; see `slabs_for_tolerance`, which does
        the search and which is the cheaper way to run a scan, since it
        can be called once and its answer reused.
    atol : float, optional
        Absolute tolerance on every returned probability.  Default:
        None.  When both are given the threshold is
        ``atol + rtol*abs(P)``, the convention of `numpy.isclose`.
    n_max : int, optional
        Largest subdivision the refinement may try before giving up and
        raising.  Default: `slabs.N_SLABS_MAX`.
    return_n_slabs : bool, optional
        Whether to return the subdivision actually used alongside the
        probabilities.  Default: False.  Worth setting when a tolerance
        is in play, since a tight one can quietly cost a great deal of
        refinement.

    Returns
    -------
    tuple of float or numpy.ndarray
        The probabilities
        :math:`P_{ee}, P_{e\mu}, P_{\mu e}, P_{\mu\mu}` as a tuple for a
        scalar energy, or an array of shape ``(..., 4)`` in the same
        order for an array of energies.

    Raises
    ------
    ValueError
        If ``costhz >= 0`` or ``n_slabs_per_segment`` is not positive.
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np
        import earth
        import globaldefs as gd
        import hamiltonians2nu

        h_vac = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(
            gd.S12_NO_BF, gd.D21_NO_BF)

        # One energy, then a whole scan through the same call
        print('%.6f' % earth.probabilities_2nu_earth(h_vac, 1.0e9, -1.0)[1])

        energies = np.array([1.0, 5.0, 10.0])*1.0e9
        prob = earth.probabilities_2nu_earth(h_vac, energies, -1.0)
        print(prob.shape, np.round(prob[:, 1], 4))
    """
    return _probabilities_earth_tol(
        h_vacuum_energy_independent, energy, costhz, n_slabs_per_segment,
        electron_fraction, 2, rtol, atol, n_max, return_n_slabs,
        'probabilities_2nu_earth')


def probabilities_3nu_earth(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float, list, np.ndarray],
    costhz: Union[int, float, list, np.ndarray],
    n_slabs_per_segment: int = 8,
    electron_fraction: float = gd.ELECTRON_FRACTION_EARTH_CRUST,
    rtol: Optional[float] = None,
    atol: Optional[float] = None,
    n_max: int = slabs.N_SLABS_MAX,
    return_n_slabs: bool = False
) -> Union[Tuple[float, float, float, float, float, float, float, float,
                 float], np.ndarray]:
    r"""Returns the three-flavor probabilities across the Earth.

    Builds the PREM slabs along the chord for the given direction and
    propagates through them exactly, slab by slab.  Raising
    ``n_slabs_per_segment`` and watching the result settle is the way to
    confirm the discretisation is fine enough for the energy in
    question; the number needed grows as the oscillation length falls.

    .. versionadded:: 1.8.0

    .. versionchanged:: 1.12.0
       Accepts an array of energies, returning one row of probabilities
       per energy.  A scalar energy returns exactly what it returned
       before.

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent three-flavor vacuum Hamiltonian, as returned
        by `hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent`.
    energy : int, float or array_like
        Neutrino energy, in units of eV, or an array of energies.  The
        whole array crosses the same chord, so the geometry and the
        matter potentials are built once for the scan rather than once
        per energy, and the chords are composed in a single pass.
    costhz : int, float or array_like
        Cosine of the zenith angle of the neutrino direction.  Must be
        negative.  May be an array, which is how an oscillogram is
        asked for: index the energies and the angles on different axes,
        as ``energy[None, :]`` against ``costhz[:, None]``, and they
        broadcast into a grid.  Each distinct angle costs one pass, so a
        grid is far cheaper than a loop over its points.
    n_slabs_per_segment : int, optional
        Number of equal sub-slabs per chord segment.  A segment runs
        between consecutive PREM boundary crossings; a chord crosses most
        shells twice, so there are more segments than shells.  Default: 8.
    electron_fraction : float, optional
        Electrons per nucleon.  Default:
        `globaldefs.ELECTRON_FRACTION_EARTH_CRUST`.
    rtol : float, optional
        Relative tolerance on every returned probability.  Default:
        None, meaning ``n_slabs_per_segment`` is used as given.  When
        either tolerance is set, the chord is refined until the measured
        discretisation error meets it, starting from
        ``n_slabs_per_segment``; see `slabs_for_tolerance`, which does
        the search and which is the cheaper way to run a scan, since it
        can be called once and its answer reused.
    atol : float, optional
        Absolute tolerance on every returned probability.  Default:
        None.  When both are given the threshold is
        ``atol + rtol*abs(P)``, the convention of `numpy.isclose`.
    n_max : int, optional
        Largest subdivision the refinement may try before giving up and
        raising.  Default: `slabs.N_SLABS_MAX`.
    return_n_slabs : bool, optional
        Whether to return the subdivision actually used alongside the
        probabilities.  Default: False.  Worth setting when a tolerance
        is in play, since a tight one can quietly cost a great deal of
        refinement.

    Returns
    -------
    tuple of float or numpy.ndarray
        The nine probabilities, with the initial flavor varying slowest,
        as a tuple for a scalar energy or an array of shape ``(..., 9)``
        in the same order for an array of energies.

    Raises
    ------
    ValueError
        If ``costhz >= 0`` or ``n_slabs_per_segment`` is not positive.
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np
        import earth
        import globaldefs as gd
        import hamiltonians3nu

        h_vac = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
            gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, gd.DCP_NO_BF,
            gd.D21_NO_BF, gd.D31_NO_BF)

        # A core-crossing chord at one energy
        print('P_mue = %.6f'
              % earth.probabilities_3nu_earth(h_vac, 1.0e10, -1.0)[3])

        # An array of energies returns the whole scan from one call
        energies = np.array([1.0, 5.0, 10.0])*1.0e9
        prob = earth.probabilities_3nu_earth(h_vac, energies, -1.0)
        print(prob.shape, np.round(prob[:, 3], 4))
    """
    return _probabilities_earth_tol(
        h_vacuum_energy_independent, energy, costhz, n_slabs_per_segment,
        electron_fraction, 3, rtol, atol, n_max, return_n_slabs,
        'probabilities_3nu_earth')


def probabilities_2nu_between_locations(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float, list, np.ndarray],
    loc_name_1: str,
    loc_name_2: str,
    n_slabs_per_segment: int = 8,
    electron_fraction: float = gd.ELECTRON_FRACTION_EARTH_CRUST,
    rtol: Optional[float] = None,
    atol: Optional[float] = None,
    n_max: int = slabs.N_SLABS_MAX,
    return_n_slabs: bool = False
) -> Union[Tuple[float, float, float, float], np.ndarray]:
    r"""Returns the two-flavor probabilities between two named locations.

    Convenience wrapper: looks both locations up in `LOC_COORDS_DMS`,
    finds the zenith angle of the chord joining them, and evaluates
    `probabilities_2nu_earth` along it.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent two-flavor vacuum Hamiltonian.
    energy : int, float or array_like
        Neutrino energy, in units of eV, or an array of energies, which
        `probabilities_3nu_earth` and its siblings evaluate as a single
        scan across the shared chord.
    loc_name_1 : str
        Name of the source location, e.g. ``'fermilab'``.
    loc_name_2 : str
        Name of the detector location, e.g. ``'homestake'``.
    n_slabs_per_segment : int, optional
        Number of equal sub-slabs per chord segment.  A segment runs
        between consecutive PREM boundary crossings; a chord crosses most
        shells twice, so there are more segments than shells.  Default: 8.
    electron_fraction : float, optional
        Electrons per nucleon.  Default:
        `globaldefs.ELECTRON_FRACTION_EARTH_CRUST`.
    rtol : float, optional
        Relative tolerance on every returned probability.  Default:
        None, meaning ``n_slabs_per_segment`` is used as given.  When
        either tolerance is set, the chord is refined until the measured
        discretisation error meets it, starting from
        ``n_slabs_per_segment``; see `slabs_for_tolerance`, which does
        the search and which is the cheaper way to run a scan, since it
        can be called once and its answer reused.
    atol : float, optional
        Absolute tolerance on every returned probability.  Default:
        None.  When both are given the threshold is
        ``atol + rtol*abs(P)``, the convention of `numpy.isclose`.
    n_max : int, optional
        Largest subdivision the refinement may try before giving up and
        raising.  Default: `slabs.N_SLABS_MAX`.
    return_n_slabs : bool, optional
        Whether to return the subdivision actually used alongside the
        probabilities.  Default: False.  Worth setting when a tolerance
        is in play, since a tight one can quietly cost a great deal of
        refinement.

    Returns
    -------
    tuple of float or numpy.ndarray
        The probabilities
        :math:`P_{ee}, P_{e\mu}, P_{\mu e}, P_{\mu\mu}` as a tuple for
        a scalar energy, or an array of shape ``(..., 4)`` for an array
        of them.

    Raises
    ------
    ValueError
        If either name is not predefined, or if the two locations
        coincide, so that there is no chord between them.
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np
        import earth
        import globaldefs as gd
        import hamiltonians2nu

        h_vac = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(
            gd.S12_NO_BF, gd.D21_NO_BF)

        prob = earth.probabilities_2nu_between_locations(
            h_vac, 1.0e9, 'cern', 'gran_sasso')
        print('P_ee = %.6f' % prob[0])
    """
    costhz = _costhz_of_named_pair(loc_name_1, loc_name_2)

    return probabilities_2nu_earth(h_vacuum_energy_independent, energy,
                                   costhz, n_slabs_per_segment,
                                   electron_fraction, rtol, atol, n_max,
                                   return_n_slabs)


def probabilities_3nu_between_locations(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float, list, np.ndarray],
    loc_name_1: str,
    loc_name_2: str,
    n_slabs_per_segment: int = 8,
    electron_fraction: float = gd.ELECTRON_FRACTION_EARTH_CRUST,
    rtol: Optional[float] = None,
    atol: Optional[float] = None,
    n_max: int = slabs.N_SLABS_MAX,
    return_n_slabs: bool = False
) -> Union[Tuple[float, float, float, float, float, float, float, float,
                 float], np.ndarray]:
    r"""Returns the three-flavor probabilities between two named locations.

    Convenience wrapper: looks both locations up in `LOC_COORDS_DMS`,
    finds the zenith angle of the chord joining them, and evaluates
    `probabilities_3nu_earth` along it.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent three-flavor vacuum Hamiltonian.
    energy : int, float or array_like
        Neutrino energy, in units of eV, or an array of energies, which
        `probabilities_3nu_earth` and its siblings evaluate as a single
        scan across the shared chord.
    loc_name_1 : str
        Name of the source location, e.g. ``'cern'``.
    loc_name_2 : str
        Name of the detector location, e.g. ``'gran_sasso'``.
    n_slabs_per_segment : int, optional
        Number of equal sub-slabs per chord segment.  A segment runs
        between consecutive PREM boundary crossings; a chord crosses most
        shells twice, so there are more segments than shells.  Default: 8.
    electron_fraction : float, optional
        Electrons per nucleon.  Default:
        `globaldefs.ELECTRON_FRACTION_EARTH_CRUST`.
    rtol : float, optional
        Relative tolerance on every returned probability.  Default:
        None, meaning ``n_slabs_per_segment`` is used as given.  When
        either tolerance is set, the chord is refined until the measured
        discretisation error meets it, starting from
        ``n_slabs_per_segment``; see `slabs_for_tolerance`, which does
        the search and which is the cheaper way to run a scan, since it
        can be called once and its answer reused.
    atol : float, optional
        Absolute tolerance on every returned probability.  Default:
        None.  When both are given the threshold is
        ``atol + rtol*abs(P)``, the convention of `numpy.isclose`.
    n_max : int, optional
        Largest subdivision the refinement may try before giving up and
        raising.  Default: `slabs.N_SLABS_MAX`.
    return_n_slabs : bool, optional
        Whether to return the subdivision actually used alongside the
        probabilities.  Default: False.  Worth setting when a tolerance
        is in play, since a tight one can quietly cost a great deal of
        refinement.

    Returns
    -------
    tuple of float or numpy.ndarray
        The nine probabilities, with the initial flavor varying slowest,
        as a tuple for a scalar energy or an array of shape ``(..., 9)``
        for an array of them.

    Raises
    ------
    ValueError
        If either name is not predefined, or if the two locations
        coincide, so that there is no chord between them.
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np
        import earth
        import globaldefs as gd
        import hamiltonians3nu

        h_vac = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
            gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, gd.DCP_NO_BF,
            gd.D21_NO_BF, gd.D31_NO_BF)

        # The CNGS chord, looked up rather than given
        prob = earth.probabilities_3nu_between_locations(
            h_vac, 1.0e9, 'cern', 'gran_sasso')
        print('P_mue = %.6f' % prob[3])
    """
    costhz = _costhz_of_named_pair(loc_name_1, loc_name_2)

    return probabilities_3nu_earth(h_vacuum_energy_independent, energy,
                                   costhz, n_slabs_per_segment,
                                   electron_fraction, rtol, atol, n_max,
                                   return_n_slabs)


def _costhz_of_named_pair(loc_name_1: str, loc_name_2: str) -> float:
    r"""Returns the chord zenith angle between two named locations.

    Parameters
    ----------
    loc_name_1 : str
        Name of the first location.
    loc_name_2 : str
        Name of the second location.

    Returns
    -------
    float
        Cosine of the zenith angle of the connecting chord.

    Raises
    ------
    ValueError
        If either name is unknown, or if the two coincide.
    """
    lat1, lon1 = coordinates_of_named_location(loc_name_1)
    lat2, lon2 = coordinates_of_named_location(loc_name_2)
    costhz = costhz_between_points_on_surface(lat1, lon1, lat2, lon2)

    if costhz >= 0.0:
        raise ValueError(
            'the locations %r and %r coincide, so there is no chord between '
            'them' % (loc_name_1, loc_name_2))

    return costhz


def probabilities_4nu_earth(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float, list, np.ndarray],
    costhz: Union[int, float, list, np.ndarray],
    n_slabs_per_segment: int = 8,
    electron_fraction: float = gd.ELECTRON_FRACTION_EARTH_CRUST,
    rtol: Optional[float] = None,
    atol: Optional[float] = None,
    n_max: int = slabs.N_SLABS_MAX,
    return_n_slabs: bool = False
) -> Union[Tuple[float, ...], np.ndarray]:
    r"""Returns the four-flavor probabilities across the Earth.

    Builds the PREM slabs along the chord for the given direction and
    propagates through them exactly, slab by slab, exactly as at two and
    three flavors.  The one thing that is new is the potential: a
    sterile state does not feel the neutral-current interaction, so
    :math:`V_{NC}` no longer cancels between the flavors and is built
    per slab alongside :math:`V_{CC}`.  This is what puts the sterile
    matter resonance where it belongs; see
    :func:`hamiltonians4nu.hamiltonian_4nu_matter`.

    .. versionadded:: 1.11.0

    .. versionchanged:: 1.12.0
       Accepts an array of energies, returning one row of probabilities
       per energy.  A scalar energy returns exactly what it returned
       before.

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent four-flavor vacuum Hamiltonian, as returned
        by
        `hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent`.
    energy : int, float or array_like
        Neutrino energy, in units of eV, or an array of energies.  The
        whole array crosses the same chord, so the geometry and both
        matter potentials are built once for the scan rather than once
        per energy.
    costhz : int, float or array_like
        Cosine of the zenith angle of the neutrino direction.  Must be
        negative.  May be an array, which is how an oscillogram is
        asked for: index the energies and the angles on different axes,
        as ``energy[None, :]`` against ``costhz[:, None]``, and they
        broadcast into a grid.  Each distinct angle costs one pass, so a
        grid is far cheaper than a loop over its points.
    n_slabs_per_segment : int, optional
        Number of equal sub-slabs per chord segment.  A segment runs
        between consecutive PREM boundary crossings; a chord crosses most
        shells twice, so there are more segments than shells.  Default: 8.
    electron_fraction : float, optional
        Electrons per nucleon.  The neutron fraction is taken as its
        complement, which is what sets :math:`V_{NC}`.  Default:
        `globaldefs.ELECTRON_FRACTION_EARTH_CRUST`.
    rtol : float, optional
        Relative tolerance on every returned probability.  Default:
        None, meaning ``n_slabs_per_segment`` is used as given.  When
        either tolerance is set, the chord is refined until the measured
        discretisation error meets it, starting from
        ``n_slabs_per_segment``; see `slabs_for_tolerance`, which does
        the search and which is the cheaper way to run a scan, since it
        can be called once and its answer reused.
    atol : float, optional
        Absolute tolerance on every returned probability.  Default:
        None.  When both are given the threshold is
        ``atol + rtol*abs(P)``, the convention of `numpy.isclose`.
    n_max : int, optional
        Largest subdivision the refinement may try before giving up and
        raising.  Default: `slabs.N_SLABS_MAX`.
    return_n_slabs : bool, optional
        Whether to return the subdivision actually used alongside the
        probabilities.  Default: False.  Worth setting when a tolerance
        is in play, since a tight one can quietly cost a great deal of
        refinement.

    Returns
    -------
    tuple of float or numpy.ndarray
        The sixteen probabilities, with the initial flavor varying
        slowest, as a tuple for a scalar energy or an array of shape
        ``(..., 16)`` in the same order for an array of energies.  With
        the fourth state read as sterile, the flavor order is
        :math:`(\nu_e, \nu_\mu, \nu_\tau, \nu_s)`.

    Raises
    ------
    ValueError
        If ``costhz >= 0`` or ``n_slabs_per_segment`` is not positive.
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np
        import earth
        import globaldefs as gd
        import hamiltonians4nu

        # The sterile parameters are illustrative; `globaldefs` carries
        # no best fit for them.
        h_vac = hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(
            gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF,
            np.sqrt(0.10), np.sqrt(0.10), 0.0,
            gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF, 1.0)

        prob = earth.probabilities_4nu_earth(h_vac, 1.0e10, -1.0,
                                             n_slabs_per_segment=4)
        print('%d probabilities, P_mue = %.6f' % (len(prob), prob[4]))
    """
    return _probabilities_earth_tol(
        h_vacuum_energy_independent, energy, costhz, n_slabs_per_segment,
        electron_fraction, 4, rtol, atol, n_max, return_n_slabs,
        'probabilities_4nu_earth')


def probabilities_4nu_between_locations(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float, list, np.ndarray],
    loc_name_1: str,
    loc_name_2: str,
    n_slabs_per_segment: int = 8,
    electron_fraction: float = gd.ELECTRON_FRACTION_EARTH_CRUST,
    rtol: Optional[float] = None,
    atol: Optional[float] = None,
    n_max: int = slabs.N_SLABS_MAX,
    return_n_slabs: bool = False
) -> Union[Tuple[float, ...], np.ndarray]:
    r"""Returns the four-flavor probabilities between two named locations.

    Convenience wrapper: looks both locations up in `LOC_COORDS_DMS`,
    finds the zenith angle of the chord joining them, and evaluates
    `probabilities_4nu_earth` along it.

    .. versionadded:: 1.11.0

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent four-flavor vacuum Hamiltonian.
    energy : int, float or array_like
        Neutrino energy, in units of eV, or an array of energies, which
        `probabilities_3nu_earth` and its siblings evaluate as a single
        scan across the shared chord.
    loc_name_1 : str
        Name of the source location, e.g. ``'cern'``.
    loc_name_2 : str
        Name of the detector location, e.g. ``'gran_sasso'``.
    n_slabs_per_segment : int, optional
        Number of equal sub-slabs per chord segment.  A segment runs
        between consecutive PREM boundary crossings; a chord crosses most
        shells twice, so there are more segments than shells.  Default: 8.
    electron_fraction : float, optional
        Electrons per nucleon.  Default:
        `globaldefs.ELECTRON_FRACTION_EARTH_CRUST`.
    rtol : float, optional
        Relative tolerance on every returned probability.  Default:
        None, meaning ``n_slabs_per_segment`` is used as given.  When
        either tolerance is set, the chord is refined until the measured
        discretisation error meets it, starting from
        ``n_slabs_per_segment``; see `slabs_for_tolerance`, which does
        the search and which is the cheaper way to run a scan, since it
        can be called once and its answer reused.
    atol : float, optional
        Absolute tolerance on every returned probability.  Default:
        None.  When both are given the threshold is
        ``atol + rtol*abs(P)``, the convention of `numpy.isclose`.
    n_max : int, optional
        Largest subdivision the refinement may try before giving up and
        raising.  Default: `slabs.N_SLABS_MAX`.
    return_n_slabs : bool, optional
        Whether to return the subdivision actually used alongside the
        probabilities.  Default: False.  Worth setting when a tolerance
        is in play, since a tight one can quietly cost a great deal of
        refinement.

    Returns
    -------
    tuple of float or numpy.ndarray
        The sixteen probabilities, with the initial flavor varying
        slowest, as a tuple for a scalar energy or an array of shape
        ``(..., 16)`` for an array of them.

    Raises
    ------
    ValueError
        If either name is not predefined, or if the two locations
        coincide, so that there is no chord between them.
    
    Examples
    --------
    .. jupyter-execute::

        import numpy as np
        import earth
        import globaldefs as gd
        import hamiltonians4nu

        # The sterile parameters are illustrative; `globaldefs` carries
        # no best fit for them.
        h_vac = hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(
            gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF,
            np.sqrt(0.10), np.sqrt(0.10), 0.0,
            gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF, 1.0)

        prob = earth.probabilities_4nu_between_locations(
            h_vac, 1.0e9, 'fermilab', 'homestake', n_slabs_per_segment=4)
        print('P_mue = %.6f' % prob[4])
    """
    costhz = _costhz_of_named_pair(loc_name_1, loc_name_2)

    return probabilities_4nu_earth(h_vacuum_energy_independent, energy,
                                   costhz, n_slabs_per_segment,
                                   electron_fraction, rtol, atol, n_max,
                                   return_n_slabs)
