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
    * distance_traveled_inside_earth - Chord length for a given costhz
    * earth_radial_distance_from_depth - Radius at a point on the chord
    * prem_layer_edges_along_chord - Where a chord crosses PREM shells
    * chord_length_inside_earth - Chord between two surface locations
    * costhz_between_points_on_surface - Its zenith angle
    * earth_slabs - Slab widths and densities along a chord
    * probabilities_2nu_earth - Two-flavor probabilities across the Earth
    * probabilities_3nu_earth - Three-flavor probabilities across the Earth
    * probabilities_2nu_between_locations - Between two named sites
    * probabilities_3nu_between_locations - Between two named sites
"""

from __future__ import print_function

__version__ = "1.8"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

__all__ = ['LOC_COORDS_DMS', 'PREM_BOUNDARIES',
           'dms_to_decimal', 'coordinates_of_named_location',
           'density_prem', 'matter_potential',
           'distance_traveled_inside_earth',
           'earth_radial_distance_from_depth',
           'prem_layer_edges_along_chord', 'chord_length_inside_earth',
           'costhz_between_points_on_surface', 'earth_slabs',
           'probabilities_2nu_earth', 'probabilities_3nu_earth',
           'probabilities_2nu_between_locations',
           'probabilities_3nu_between_locations']

from typing import Optional, Tuple, Union

import numpy as np

import globaldefs as gd
import hamiltonians2nu
import hamiltonians3nu
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
        5 degrees South is ``(-5, ...)``.
    minutes : int or float
        Minute part of the coordinate, taken as positive.
    seconds : int or float
        Second part of the coordinate, taken as positive.

    Returns
    -------
    float
        The coordinate in decimal degrees.

    Examples
    --------
    >>> print('%.6f' % dms_to_decimal(36, 25, 50.05))
    36.430569
    """
    magnitude = abs(degrees) + minutes/60.0 + seconds/3600.0

    # The sign lives on the degree part, so a negative latitude must not
    # have its minutes and seconds added back the other way.
    return -magnitude if degrees < 0 else magnitude


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
    >>> lat, lon = coordinates_of_named_location('South Pole')
    >>> print(lat, lon)
    (-90, 0, 0) (0, 0, 0)
    """
    key = loc_name.lower().replace(' ', '_')
    try:
        entry = LOC_COORDS_DMS[key]
    except KeyError:
        raise ValueError(
            'coordinates_of_named_location: %r is not a predefined '
            'location; the available names, in earth.LOC_COORDS_DMS, are: '
            '%s' % (loc_name, ', '.join(sorted(LOC_COORDS_DMS))))

    return entry['lat'], entry['lon']


def density_prem(
    r: Union[int, float, list, np.ndarray],
    tol: Optional[float] = 1.e-8
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
    >>> print('%.4f' % density_prem(0.0))
    13.0885
    >>> print('%.4f' % density_prem(6371.0))
    1.0200
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
    electron_fraction: Optional[float] = gd.ELECTRON_FRACTION_EARTH_CRUST
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
    >>> print('%.4e' % matter_potential(3.0))
    1.1356e-13
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


def distance_traveled_inside_earth(costhz: Union[int, float]) -> float:
    r"""Returns the chord length through the Earth for a given direction.

    The neutrino is assumed to arrive at a detector on the surface, not
    underground, so the distance is zero for any upward-going direction,
    ``costhz >= 0``.

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

    Examples
    --------
    >>> print('%.1f' % distance_traveled_inside_earth(-1.0))
    12742.0
    >>> print('%.1f' % distance_traveled_inside_earth(0.5))
    0.0
    """
    return 0.0 if costhz >= 0.0 else -2.0*gd.EARTH_RADIUS*costhz


def earth_radial_distance_from_depth(
    costhz: Union[int, float],
    l: Union[int, float, list, np.ndarray],
    tol: Optional[float] = 1.e-8
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
    >>> print('%.1f' % earth_radial_distance_from_depth(-1.0, 6371.0))
    0.0
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
    r = np.sqrt(np.abs(r2))

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
    >>> print(len(prem_layer_edges_along_chord(-1.0)))
    18
    >>> print(len(prem_layer_edges_along_chord(0.5)))
    0
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

    return np.unique(np.array(sorted(crossings)))


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
    >>> lat1, lon1 = coordinates_of_named_location('fermilab')
    >>> lat2, lon2 = coordinates_of_named_location('homestake')
    >>> print('%.1f' % chord_length_inside_earth(lat1, lon1, lat2, lon2))
    1284.7
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
    >>> lat1, lon1 = coordinates_of_named_location('cern')
    >>> lat2, lon2 = coordinates_of_named_location('gran_sasso')
    >>> print('%.6f' % costhz_between_points_on_surface(lat1, lon1,
    ...                                                 lat2, lon2))
    -0.057179
    """
    chord = chord_length_inside_earth(lat1_dms, lon1_dms, lat2_dms, lon2_dms)

    return -0.5*chord/gd.EARTH_RADIUS


def earth_slabs(
    costhz: Union[int, float],
    n_slabs_per_segment: Optional[int] = 8
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
    >>> widths, densities = earth_slabs(-1.0, n_slabs_per_segment=2)
    >>> print(len(widths), '%.1f' % sum(widths))
    38 12742.0
    """
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

    return widths, densities


def _earth_hamiltonians(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float],
    costhz: Union[int, float],
    n_slabs_per_segment: int,
    electron_fraction: float,
    n_flavors: int
) -> Tuple[np.ndarray, np.ndarray]:
    r"""Returns the per-slab Hamiltonians and widths for a chord.

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent vacuum Hamiltonian.
    energy : int or float
        Neutrino energy, in units of eV.
    costhz : int or float
        Cosine of the zenith angle.
    n_slabs_per_segment : int
        Number of equal sub-slabs per chord segment.  A segment runs
        between consecutive PREM boundary crossings; a chord crosses most
        shells twice, so there are more segments than shells.
    electron_fraction : float
        Electrons per nucleon.
    n_flavors : int
        Number of neutrino flavors, 2 or 3.

    Returns
    -------
    tuple of numpy.ndarray
        The Hamiltonians, of shape ``(n, n_flavors, n_flavors)``, and
        the slab widths in units of eV\ :sup:`-1`.
    """
    widths_km, densities = earth_slabs(costhz, n_slabs_per_segment)
    potentials = matter_potential(densities, electron_fraction)

    # The Hamiltonian builders take an array of potentials and return one
    # Hamiltonian per entry, so the whole chord is built in one call.
    if n_flavors == 2:
        h = hamiltonians2nu.hamiltonian_2nu_matter(
            h_vacuum_energy_independent, energy, potentials)
    else:
        h = hamiltonians3nu.hamiltonian_3nu_matter(
            h_vacuum_energy_independent, energy, potentials)

    return h, widths_km*gd.CONV_KM_TO_INV_EV


def probabilities_2nu_earth(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float],
    costhz: Union[int, float],
    n_slabs_per_segment: Optional[int] = 8,
    electron_fraction: Optional[float] = gd.ELECTRON_FRACTION_EARTH_CRUST
) -> Tuple[float, float, float, float]:
    r"""Returns the two-flavor probabilities across the Earth.

    Builds the PREM slabs along the chord for the given direction and
    propagates through them exactly, slab by slab.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent two-flavor vacuum Hamiltonian, as returned by
        `hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent`.
    energy : int or float
        Neutrino energy, in units of eV.
    costhz : int or float
        Cosine of the zenith angle of the neutrino direction.  Must be
        negative.
    n_slabs_per_segment : int, optional
        Number of equal sub-slabs per chord segment.  A segment runs
        between consecutive PREM boundary crossings; a chord crosses most
        shells twice, so there are more segments than shells.  Default: 8.
    electron_fraction : float, optional
        Electrons per nucleon.  Default:
        `globaldefs.ELECTRON_FRACTION_EARTH_CRUST`.

    Returns
    -------
    tuple of float
        The probabilities
        :math:`P_{ee}, P_{e\mu}, P_{\mu e}, P_{\mu\mu}`.

    Raises
    ------
    ValueError
        If ``costhz >= 0`` or ``n_slabs_per_segment`` is not positive.
    """
    h, widths = _earth_hamiltonians(h_vacuum_energy_independent, energy,
                                    costhz, n_slabs_per_segment,
                                    electron_fraction, 2)

    return slabs.probabilities_2nu_slabs(h, widths)


def probabilities_3nu_earth(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float],
    costhz: Union[int, float],
    n_slabs_per_segment: Optional[int] = 8,
    electron_fraction: Optional[float] = gd.ELECTRON_FRACTION_EARTH_CRUST
) -> Tuple[float, float, float, float, float, float, float, float, float]:
    r"""Returns the three-flavor probabilities across the Earth.

    Builds the PREM slabs along the chord for the given direction and
    propagates through them exactly, slab by slab.  Raising
    ``n_slabs_per_segment`` and watching the result settle is the way to
    confirm the discretisation is fine enough for the energy in
    question; the number needed grows as the oscillation length falls.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent three-flavor vacuum Hamiltonian, as returned
        by `hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent`.
    energy : int or float
        Neutrino energy, in units of eV.
    costhz : int or float
        Cosine of the zenith angle of the neutrino direction.  Must be
        negative.
    n_slabs_per_segment : int, optional
        Number of equal sub-slabs per chord segment.  A segment runs
        between consecutive PREM boundary crossings; a chord crosses most
        shells twice, so there are more segments than shells.  Default: 8.
    electron_fraction : float, optional
        Electrons per nucleon.  Default:
        `globaldefs.ELECTRON_FRACTION_EARTH_CRUST`.

    Returns
    -------
    tuple of float
        The nine probabilities, with the initial flavor varying slowest.

    Raises
    ------
    ValueError
        If ``costhz >= 0`` or ``n_slabs_per_segment`` is not positive.
    """
    h, widths = _earth_hamiltonians(h_vacuum_energy_independent, energy,
                                    costhz, n_slabs_per_segment,
                                    electron_fraction, 3)

    return slabs.probabilities_3nu_slabs(h, widths)


def probabilities_2nu_between_locations(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float],
    loc_name_1: str,
    loc_name_2: str,
    n_slabs_per_segment: Optional[int] = 8,
    electron_fraction: Optional[float] = gd.ELECTRON_FRACTION_EARTH_CRUST
) -> Tuple[float, float, float, float]:
    r"""Returns the two-flavor probabilities between two named locations.

    Convenience wrapper: looks both locations up in `LOC_COORDS_DMS`,
    finds the zenith angle of the chord joining them, and evaluates
    `probabilities_2nu_earth` along it.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent two-flavor vacuum Hamiltonian.
    energy : int or float
        Neutrino energy, in units of eV.
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

    Returns
    -------
    tuple of float
        The probabilities
        :math:`P_{ee}, P_{e\mu}, P_{\mu e}, P_{\mu\mu}`.

    Raises
    ------
    ValueError
        If either name is not predefined, or if the two locations
        coincide, so that there is no chord between them.
    """
    costhz = _costhz_of_named_pair(loc_name_1, loc_name_2)

    return probabilities_2nu_earth(h_vacuum_energy_independent, energy,
                                   costhz, n_slabs_per_segment,
                                   electron_fraction)


def probabilities_3nu_between_locations(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float],
    loc_name_1: str,
    loc_name_2: str,
    n_slabs_per_segment: Optional[int] = 8,
    electron_fraction: Optional[float] = gd.ELECTRON_FRACTION_EARTH_CRUST
) -> Tuple[float, float, float, float, float, float, float, float, float]:
    r"""Returns the three-flavor probabilities between two named locations.

    Convenience wrapper: looks both locations up in `LOC_COORDS_DMS`,
    finds the zenith angle of the chord joining them, and evaluates
    `probabilities_3nu_earth` along it.

    .. versionadded:: 1.8.0

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent three-flavor vacuum Hamiltonian.
    energy : int or float
        Neutrino energy, in units of eV.
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

    Returns
    -------
    tuple of float
        The nine probabilities, with the initial flavor varying slowest.

    Raises
    ------
    ValueError
        If either name is not predefined, or if the two locations
        coincide, so that there is no chord between them.
    """
    costhz = _costhz_of_named_pair(loc_name_1, loc_name_2)

    return probabilities_3nu_earth(h_vacuum_energy_independent, energy,
                                   costhz, n_slabs_per_segment,
                                   electron_fraction)


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
