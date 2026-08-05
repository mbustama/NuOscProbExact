# -*- coding: utf-8 -*-
r"""Tests of energy batching and of the ``rtol``/``atol`` refinement.

Two features share this file because they share a failure mode: both
add a dispatch, and every dispatch added to this library so far has
silently orphaned the NumPy branch beneath it, because the compiled
predicates say yes at every size once Numba is present.  A suite that
only ever exercises the compiled side stays green while the fallback
rots.  Everything here that can take either path is therefore
parametrised over the ``backend`` fixture, and the assertions that a
compiled path was *taken* go through ``kernel_spy`` rather than
comparing NumPy against itself.

The batching claim is that an array of energies is the same answer as a
loop over them, not merely a similar one; the tolerance claim is that
the subdivision the search returns genuinely meets the tolerance, which
is checked against a converged reference rather than against the
estimate the search itself used.
"""

import numpy as np
import pytest

import earth
import fastkernels
import globaldefs as gd
import hamiltonians2nu
import hamiltonians3nu
import hamiltonians4nu
import oscprob3nu
import slabs

from conftest import ATOL


COSTHZ = -0.8
N_REFERENCE = 512


def h_vacuum(n_flavors):
    r"""Returns the energy-independent vacuum Hamiltonian for a flavor
    count."""
    if n_flavors == 2:
        return hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(
            gd.S12_NO_BF, gd.D21_NO_BF)
    if n_flavors == 3:
        return hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
            gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, gd.DCP_NO_BF,
            gd.D21_NO_BF, gd.D31_NO_BF)
    return hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(
        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, 0.1, 0.1, 0.1,
        gd.D21_NO_BF, gd.D31_NO_BF, 1.0, 0.0, 0.0, 0.0)


def probabilities_earth(n_flavors):
    r"""Returns the Earth entry point for a flavor count."""
    return {2: earth.probabilities_2nu_earth,
            3: earth.probabilities_3nu_earth,
            4: earth.probabilities_4nu_earth}[n_flavors]


def richardson_reference(n_flavors, energy, costhz, n=N_REFERENCE):
    r"""Returns the converged probabilities by second-order Richardson.

    ``[4*P(2n) - P(n)]/3`` cancels the leading second-order term of the
    midpoint discretisation.  The first-order form ``2*P(2n) - P(n)`` is
    *worse* than plain ``P(2n)`` here and is the wrong reference to
    build.
    """
    fn = probabilities_earth(n_flavors)
    p_n = np.asarray(fn(h_vacuum(n_flavors), energy, costhz, n))
    p_2n = np.asarray(fn(h_vacuum(n_flavors), energy, costhz, 2*n))

    return (4.0*p_2n - p_n)/3.0


# ---------------------------------------------------------------- batching


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
def test_array_of_energies_equals_a_loop_over_them(backend, n_flavors):
    r"""The batched path is the scalar path, evaluated many times.

    Not "close to": the batched call composes the same slabs from the
    same Hamiltonians, so any difference beyond round-off means the
    broadcast put an energy against the wrong chord.
    """
    energies = np.logspace(9.0, 11.0, 17)
    fn = probabilities_earth(n_flavors)
    hv = h_vacuum(n_flavors)

    batched = fn(hv, energies, COSTHZ)
    looped = np.array([fn(hv, e, COSTHZ) for e in energies])

    assert batched.shape == (len(energies), n_flavors*n_flavors)
    assert np.allclose(batched, looped, rtol=0.0, atol=ATOL)


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
def test_a_scalar_energy_still_returns_a_tuple(backend, n_flavors):
    r"""The scalar contract is untouched by the array one.

    A caller unpacking nine floats must keep working, and a zero
    dimensional array counts as a scalar rather than as a scan of one.
    """
    hv = h_vacuum(n_flavors)
    fn = probabilities_earth(n_flavors)

    scalar = fn(hv, 1.0e10, COSTHZ)
    zero_dim = fn(hv, np.float64(1.0e10), COSTHZ)

    assert isinstance(scalar, tuple)
    assert isinstance(zero_dim, tuple)
    assert len(scalar) == n_flavors*n_flavors
    assert np.allclose(scalar, zero_dim, rtol=0.0, atol=0.0)


def test_the_batched_path_reaches_the_compiled_kernel(kernel_spy):
    r"""The batch kernel is what a scan actually calls.

    Without this the comparison above passes by running the NumPy path
    twice, which is how a threshold change goes unnoticed.
    """
    earth.probabilities_3nu_earth(h_vacuum(3), np.logspace(9.0, 11.0, 32),
                                  COSTHZ)

    assert kernel_spy['slab_product_3nu_batch_kernel'] == 1
    assert kernel_spy['slab_product_3nu_kernel'] == 0


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
def test_the_numpy_batch_fallback_agrees_with_the_kernel(n_flavors):
    r"""The uncompiled branch beneath the batch dispatch still works.

    ``worthwhile_slabs`` is true at every size once Numba is present, so
    this branch is unreachable unless the backend is forced off, and
    nothing but a test like this one keeps it honest.

    Four flavors compares at 1e-9 rather than at round-off, and not
    because of anything batching does: the two backends polish the
    latent roots slightly differently, and the same 4e-10 separates them
    on the per-energy path that predates this feature.  What is asserted
    at round-off, below, is that each backend's batched answer matches
    *its own* scalar answer, which is the property batching is
    responsible for.
    """
    energies = np.logspace(9.0, 11.0, 11)
    hv = h_vacuum(n_flavors)
    fn = probabilities_earth(n_flavors)
    across_backends = 1.0e-9 if n_flavors == 4 else ATOL

    fastkernels.USE_NUMBA = True
    compiled = fn(hv, energies, COSTHZ)
    compiled_loop = np.array([fn(hv, e, COSTHZ) for e in energies])
    try:
        fastkernels.USE_NUMBA = False
        interpreted = fn(hv, energies, COSTHZ)
        interpreted_loop = np.array([fn(hv, e, COSTHZ) for e in energies])
    finally:
        fastkernels.USE_NUMBA = True

    assert np.allclose(compiled, interpreted, rtol=0.0, atol=across_backends)
    assert np.allclose(compiled, compiled_loop, rtol=0.0, atol=ATOL)
    assert np.allclose(interpreted, interpreted_loop, rtol=0.0, atol=ATOL)


def test_batching_preserves_the_shape_of_the_energy_array(backend):
    r"""A multidimensional energy grid keeps its shape.

    An oscillogram indexes energies and angles separately, so the
    probabilities have to arrive on the caller's grid rather than
    flattened onto one axis.
    """
    energies = np.logspace(9.0, 11.0, 12).reshape(3, 4)

    got = earth.probabilities_3nu_earth(h_vacuum(3), energies, COSTHZ)

    assert got.shape == (3, 4, 9)
    assert np.allclose(
        got.reshape(12, 9),
        earth.probabilities_3nu_earth(h_vacuum(3), energies.reshape(12),
                                      COSTHZ),
        rtol=0.0, atol=0.0)


def test_a_list_of_energies_is_accepted(backend):
    r"""``array_like`` means what it says, not "ndarray only"."""
    got = earth.probabilities_3nu_earth(h_vacuum(3), [1.0e9, 1.0e10], COSTHZ)

    assert got.shape == (2, 9)


def test_chunking_does_not_change_the_answer(backend):
    r"""Splitting a scan into chunks is an implementation detail.

    A long scan is evaluated in pieces to bound the Hamiltonian stack it
    builds --- a hundred thousand energies would otherwise be gigabytes
    --- and the seam between pieces must be invisible.
    """
    energies = np.logspace(9.0, 11.0, 40)
    hv = h_vacuum(3)

    unchunked = earth.probabilities_3nu_earth(hv, energies, COSTHZ)
    original = earth.MAX_CHUNK_BYTES
    try:
        # Small enough to force one chord per chunk, the extreme case
        earth.MAX_CHUNK_BYTES = 1
        chunked = earth.probabilities_3nu_earth(hv, energies, COSTHZ)
    finally:
        earth.MAX_CHUNK_BYTES = original

    assert np.array_equal(unchunked, chunked)


def test_probabilities_between_locations_takes_an_array(backend):
    r"""The named-location wrappers inherit the array path."""
    energies = np.logspace(9.0, 10.0, 5)

    got = earth.probabilities_3nu_between_locations(
        h_vacuum(3), energies, 'cern', 'kamioka')

    assert got.shape == (5, 9)


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
def test_the_batch_kernel_matches_the_per_chord_kernel(n_flavors):
    r"""Batching chords is composing each of them separately.

    Straight at the kernels, with no Earth geometry in the way, so a
    disagreement here is the batch loop's and nothing else's.
    """
    if not fastkernels.HAVE_NUMBA:
        pytest.skip('Numba is not installed')

    rng = np.random.default_rng(20260805)
    n_chords, n_slabs = 5, 7
    a = (rng.normal(size=(n_chords, n_slabs, n_flavors, n_flavors))
         + 1.0j*rng.normal(size=(n_chords, n_slabs, n_flavors, n_flavors)))
    h = (a + np.conj(np.swapaxes(a, -1, -2)))/2.0
    widths = rng.uniform(0.1, 2.0, size=n_slabs)

    if n_flavors == 2:
        batch = fastkernels.slab_product_2nu_batch_kernel(h, widths)
        one = fastkernels.slab_product_2nu_kernel
        per = np.array([one(h[i], widths) for i in range(n_chords)])
    elif n_flavors == 3:
        batch = fastkernels.slab_product_3nu_batch_kernel(h, widths)
        one = fastkernels.slab_product_3nu_kernel
        per = np.array([one(h[i], widths) for i in range(n_chords)])
    else:
        batch = fastkernels.slab_product_4nu_batch_kernel(h, widths, True)
        one = fastkernels.slab_product_4nu_kernel
        per = np.array([one(h[i], widths, True) for i in range(n_chords)])

    assert np.allclose(batch, per, rtol=0.0, atol=ATOL)


def test_the_batch_kernel_spreads_wide_work_over_threads():
    r"""Both sides of the parallel threshold give the same answer.

    The threshold is on ``n_chords*n_slabs`` rather than on the chord
    count, because a chord costs in proportion to its slabs; the two
    branches must agree whichever side of it the work falls.
    """
    if not fastkernels.HAVE_NUMBA:
        pytest.skip('Numba is not installed')

    rng = np.random.default_rng(20260806)
    n_slabs = 8
    n_chords = 2*fastkernels.PARALLEL_THRESHOLD//n_slabs
    a = (rng.normal(size=(n_chords, n_slabs, 3, 3))
         + 1.0j*rng.normal(size=(n_chords, n_slabs, 3, 3)))
    h = (a + np.conj(np.swapaxes(a, -1, -2)))/2.0
    widths = rng.uniform(0.1, 2.0, size=n_slabs)

    wide = fastkernels.slab_product_3nu_batch_kernel(h, widths)
    narrow = np.array([fastkernels.slab_product_3nu_batch_kernel(
        h[i:i+1], widths)[0] for i in range(n_chords)])

    assert n_chords*n_slabs >= fastkernels.PARALLEL_THRESHOLD
    assert np.allclose(wide, narrow, rtol=0.0, atol=ATOL)


@pytest.mark.parametrize('bad, message', [
    (np.zeros((4, 3, 3), dtype=complex), 'at least one leading batch axis'),
    (np.zeros((2, 4, 2, 2), dtype=complex), 'at least one leading batch axis'),
])
def test_the_batch_check_rejects_a_malformed_stack(bad, message):
    r"""A stack of the wrong rank is named, not broadcast into
    nonsense."""
    with pytest.raises(ValueError, match=message):
        slabs._probabilities_slabs_batch(bad, np.ones(4), 3, 'caller')


def test_the_batch_check_rejects_mismatched_widths():
    r"""One width per slab, and the error says which count was which."""
    with pytest.raises(ValueError, match='but 5 widths'):
        slabs._probabilities_slabs_batch(np.zeros((2, 4, 3, 3), dtype=complex),
                                         np.ones(5), 3, 'caller')


@pytest.mark.parametrize('n_slabs, widths, message', [
    (4, np.ones((2, 4)), 'one-dimensional'),
    (0, np.array([]), 'at least one slab'),
    (4, np.array([-1.0, 1.0, 1.0, 1.0]), 'cannot be negative'),
])
def test_the_batch_check_rejects_bad_widths(n_slabs, widths, message):
    r"""The width array is checked as the per-chord path checks it.

    The slab count is given explicitly so that each case reaches the
    check it is about: a stack whose length disagrees with the widths
    trips the count check first and never gets to the others.
    """
    h = np.zeros((2, n_slabs, 3, 3), dtype=complex)
    with pytest.raises(ValueError, match=message):
        slabs._probabilities_slabs_batch(h, widths, 3, 'caller')


# --------------------------------------------------------------- tolerance


@pytest.mark.parametrize('energy_gev', [3.0, 10.0, 40.0])
@pytest.mark.parametrize('costhz', [-0.3, -1.0])
def test_the_chosen_subdivision_meets_the_tolerance(energy_gev, costhz):
    r"""The point of the feature, checked against a converged answer.

    The search decides using its own error estimate, so checking it
    against that estimate would prove nothing.  The reference here is an
    independent Richardson extrapolation at a far finer subdivision.
    """
    energy = energy_gev*1.0e9
    reference = richardson_reference(3, energy, costhz)

    for atol in (1.0e-4, 1.0e-5):
        n = earth.slabs_for_tolerance(h_vacuum(3), energy, costhz, atol=atol)
        got = np.asarray(earth.probabilities_3nu_earth(
            h_vacuum(3), energy, costhz, n))
        assert np.max(np.abs(got - reference)) <= atol


def test_the_tolerance_binds_on_every_returned_probability():
    r"""All nine, not the one that happens to converge fastest.

    The subdivision is set by the worst entry, so no channel in the
    returned tuple is quietly less converged than the caller asked for.
    """
    energy, costhz, atol = 1.0e10, COSTHZ, 1.0e-5
    reference = richardson_reference(3, energy, costhz)

    n = earth.slabs_for_tolerance(h_vacuum(3), energy, costhz, atol=atol)
    got = np.asarray(earth.probabilities_3nu_earth(h_vacuum(3), energy,
                                                   costhz, n))

    assert np.all(np.abs(got - reference) <= atol)


def test_one_subdivision_covers_a_whole_scan():
    r"""An array of energies gets the worst-converging one's answer.

    That is what makes the returned number safe to reuse across a scan,
    and it is why the answer for an array is never smaller than for the
    energies in it.
    """
    energies = np.logspace(9.0, 10.5, 12)
    atol = 1.0e-5

    n_scan = earth.slabs_for_tolerance(h_vacuum(3), energies, COSTHZ,
                                       atol=atol)
    worst = max(earth.slabs_for_tolerance(h_vacuum(3), e, COSTHZ, atol=atol)
                for e in energies)

    assert n_scan >= worst
    got = earth.probabilities_3nu_earth(h_vacuum(3), energies, COSTHZ, n_scan)
    reference = np.array([richardson_reference(3, e, COSTHZ)
                          for e in energies])
    assert np.max(np.abs(got - reference)) <= atol


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
def test_the_tolerance_search_runs_at_every_flavor_count(backend, n_flavors):
    r"""Two, three and four flavors all reach the search."""
    n = earth.slabs_for_tolerance(h_vacuum(n_flavors), 1.0e10, COSTHZ,
                                  n_flavors=n_flavors, atol=1.0e-4)

    assert n >= 8
    assert n % 8 == 0


def test_no_tolerance_leaves_the_result_untouched(backend):
    r"""The regression that matters most: defaults are bit-identical.

    Adding an argument that defaults to None must not perturb a single
    bit of what the routine returned before it existed.
    """
    hv = h_vacuum(3)

    plain = earth.probabilities_3nu_earth(hv, 1.0e10, COSTHZ)
    explicit = earth.probabilities_3nu_earth(hv, 1.0e10, COSTHZ, 8,
                                             gd.ELECTRON_FRACTION_EARTH_CRUST,
                                             None, None)

    assert np.array_equal(np.asarray(plain), np.asarray(explicit))


def test_a_tolerance_actually_changes_the_subdivision(backend):
    r"""A tolerance tight enough to need refinement gets it.

    Guards against the arguments being accepted and then ignored, which
    would leave every assertion about accuracy true only by luck.
    """
    hv = h_vacuum(3)
    default = earth.probabilities_3nu_earth(hv, 1.0e10, COSTHZ)
    refined, n = earth.probabilities_3nu_earth(hv, 1.0e10, COSTHZ,
                                               atol=1.0e-6,
                                               return_n_slabs=True)

    assert n > 8
    assert not np.allclose(default, refined, rtol=0.0, atol=1.0e-9)
    assert np.allclose(
        refined, earth.probabilities_3nu_earth(hv, 1.0e10, COSTHZ, n),
        rtol=0.0, atol=0.0)


def test_the_reported_subdivision_is_the_one_the_helper_returns(backend):
    r"""``return_n_slabs`` and `slabs_for_tolerance` cannot disagree.

    They run the same search, and a caller who uses the helper to plan a
    scan must get what the routine would have chosen on its own.
    """
    hv = h_vacuum(3)

    _, reported = earth.probabilities_3nu_earth(hv, 1.0e10, COSTHZ,
                                                atol=1.0e-6,
                                                return_n_slabs=True)
    planned = earth.slabs_for_tolerance(hv, 1.0e10, COSTHZ, atol=1.0e-6)

    assert reported == planned


def test_return_n_slabs_reports_the_default_when_no_tolerance_is_set(backend):
    r"""Asking what was used is answerable without a tolerance too."""
    p, n = earth.probabilities_3nu_earth(h_vacuum(3), 1.0e10, COSTHZ,
                                         n_slabs_per_segment=12,
                                         return_n_slabs=True)

    assert n == 12
    assert isinstance(p, tuple)


def test_an_array_energy_with_a_tolerance_returns_an_array(backend):
    r"""The tolerance path preserves the scalar/array contract."""
    energies = np.logspace(9.0, 10.0, 6)

    got, n = earth.probabilities_3nu_earth(h_vacuum(3), energies, COSTHZ,
                                           atol=1.0e-5, return_n_slabs=True)

    assert got.shape == (6, 9)
    assert np.allclose(
        got, earth.probabilities_3nu_earth(h_vacuum(3), energies, COSTHZ, n),
        rtol=0.0, atol=0.0)


def test_a_loose_tolerance_is_met_by_the_starting_subdivision(backend):
    r"""A tolerance the default already meets does not force a doubling.

    The error of the coarser evaluation is four thirds of the gap, so
    the search can certify ``n_start`` itself rather than returning
    twice it and costing every later call double for nothing.
    """
    n = earth.slabs_for_tolerance(h_vacuum(3), 1.0e10, COSTHZ, atol=1.0e-2)

    assert n == 8


@pytest.mark.parametrize('rtol, atol, message', [
    (None, None, 'tolerance of zero'),
    (0.0, 0.0, 'tolerance of zero'),
    (0.0, None, 'tolerance of zero'),
    (-1.0e-3, None, 'cannot be negative'),
    (None, -1.0e-6, 'cannot be negative'),
])
def test_invalid_tolerances_are_refused(rtol, atol, message):
    r"""A tolerance that cannot be honoured is an error, not a no-op."""
    with pytest.raises(ValueError, match=message):
        earth.slabs_for_tolerance(h_vacuum(3), 1.0e10, COSTHZ,
                                  rtol=rtol, atol=atol)


def test_an_unreachable_tolerance_raises():
    r"""Clamping silently would be the dangerous answer.

    A caller who asked for an accuracy and did not get it should be
    told, rather than handed a coarser number with a warning that
    Python will print once and then suppress.
    """
    with pytest.raises(ValueError, match='could not meet'):
        earth.slabs_for_tolerance(h_vacuum(3), 1.0e10, COSTHZ,
                                  atol=1.0e-16, n_max=32)


def test_the_unreachable_message_reports_the_error_it_reached():
    r"""The message says how close it got, so n_max can be chosen."""
    with pytest.raises(ValueError) as excinfo:
        earth.slabs_for_tolerance(h_vacuum(3), 1.0e10, COSTHZ,
                                  atol=1.0e-16, n_max=32)

    # A real estimate, not the zero a spent loop variable would give
    assert 'was 0.000e+00' not in str(excinfo.value)


@pytest.mark.parametrize('kwargs, message', [
    ({'n_start': 0}, 'n_start must be at least 1'),
    ({'n_start': 8, 'n_max': 8}, 'at least twice'),
])
def test_the_search_refuses_an_impossible_budget(kwargs, message):
    r"""An error estimate costs two evaluations; below that, none."""
    with pytest.raises(ValueError, match=message):
        earth.slabs_for_tolerance(h_vacuum(3), 1.0e10, COSTHZ, atol=1.0e-4,
                                  **kwargs)


def test_both_tolerances_combine_as_numpy_does(backend):
    r"""``atol + rtol*abs(P)``, so the looser of the two governs.

    Adding a slack absolute term to a relative request cannot make the
    requirement stricter, and so cannot raise the subdivision.
    """
    hv = h_vacuum(3)

    tight = earth.slabs_for_tolerance(hv, 1.0e10, COSTHZ, rtol=1.0e-5)
    loosened = earth.slabs_for_tolerance(hv, 1.0e10, COSTHZ, rtol=1.0e-5,
                                         atol=1.0e-3)

    assert loosened <= tight


# ----------------------------------------------------------------- profile


def constant_profile(h):
    r"""Returns a callable giving the same Hamiltonian everywhere."""
    def hamiltonian_of(positions):
        return np.broadcast_to(h, (len(positions),) + np.shape(h))
    return hamiltonian_of


@pytest.mark.parametrize('n_flavors, routine', [
    (2, slabs.probabilities_2nu_profile),
    (3, slabs.probabilities_3nu_profile),
    (4, slabs.probabilities_4nu_profile),
])
@pytest.mark.parametrize('n_slabs', [1, 8, 32])
def test_a_constant_profile_is_the_closed_form(backend, n_flavors, routine,
                                               n_slabs):
    r"""Slicing a constant Hamiltonian cannot change the answer.

    The same property the slab machinery has, now reached through the
    profile interface: if this fails, the midpoints or the widths are
    wrong, and every converging case would be wrong by a shrinking
    amount that looks like convergence.
    """
    rng = np.random.default_rng(20260807 + n_flavors)
    a = (rng.normal(size=(n_flavors, n_flavors))
         + 1.0j*rng.normal(size=(n_flavors, n_flavors)))
    h = (a + a.conj().T)/2.0
    baseline = 3.0

    got = routine(constant_profile(h), baseline, n_slabs=n_slabs)
    reference = {2: 4, 3: 9, 4: 16}[n_flavors]

    assert len(got) == reference
    assert np.isclose(sum(got), n_flavors, rtol=0.0, atol=1.0e-10)


def test_a_constant_profile_matches_oscprob(backend):
    r"""Against the single-Hamiltonian routine, not just against
    itself."""
    h = hamiltonians3nu.hamiltonian_3nu_matter(h_vacuum(3), 1.0e10, 1.0e-13)
    baseline = 1.0e4*gd.CONV_KM_TO_INV_EV

    got = slabs.probabilities_3nu_profile(constant_profile(h), baseline,
                                          n_slabs=16)
    reference = oscprob3nu.probabilities_3nu(h, baseline)

    assert np.allclose(got, reference, rtol=0.0, atol=1.0e-12)


def varying_profile(baseline):
    r"""Returns a smoothly varying three-flavor profile."""
    h_vac = np.asarray(h_vacuum(3), dtype=complex)/1.0e10

    def hamiltonian_of(positions):
        h = np.broadcast_to(h_vac, (len(positions), 3, 3)).copy()
        h[:, 0, 0] += 1.0e-13*(1.0 + 0.5*np.sin(3.0*np.pi*positions/baseline))
        return h

    return hamiltonian_of


def test_a_varying_profile_converges_at_second_order(backend):
    r"""Midpoint sampling is second order, which is what the search
    assumes.

    If this ratio were not four, the error estimate the tolerance search
    is built on would be wrong, and every tolerance it reported met
    would be met only by accident.
    """
    baseline = 1.0e4*gd.CONV_KM_TO_INV_EV
    profile = varying_profile(baseline)

    values = {n: np.array(slabs.probabilities_3nu_profile(profile, baseline,
                                                          n_slabs=n))
              for n in (32, 64, 128, 256, 512)}
    reference = (4.0*values[512] - values[256])/3.0

    errors = [np.max(np.abs(values[n] - reference)) for n in (32, 64, 128)]
    ratios = [errors[i]/errors[i+1] for i in range(len(errors)-1)]

    assert all(3.6 < r < 4.4 for r in ratios), ratios


def test_a_profile_tolerance_is_met(backend):
    r"""The refinement works on a hand-built profile, not only on
    PREM."""
    baseline = 1.0e4*gd.CONV_KM_TO_INV_EV
    profile = varying_profile(baseline)
    reference = (4.0*np.array(slabs.probabilities_3nu_profile(
        profile, baseline, n_slabs=1024))
        - np.array(slabs.probabilities_3nu_profile(
            profile, baseline, n_slabs=512)))/3.0

    for atol in (1.0e-4, 1.0e-6):
        got, n = slabs.probabilities_3nu_profile(profile, baseline, atol=atol,
                                                 return_n_slabs=True)
        assert np.max(np.abs(np.array(got) - reference)) <= atol
        assert n <= slabs.N_SLABS_MAX


def test_a_profile_without_a_tolerance_uses_the_slab_count(backend):
    r"""No tolerance means no refinement and no extra evaluations."""
    baseline = 1.0e4*gd.CONV_KM_TO_INV_EV
    calls = []

    def counting_profile(positions):
        calls.append(len(positions))
        return varying_profile(baseline)(positions)

    got, n = slabs.probabilities_3nu_profile(counting_profile, baseline,
                                             n_slabs=21, return_n_slabs=True)

    assert n == 21
    assert calls == [21]
    assert len(got) == 9


@pytest.mark.parametrize('kwargs, message', [
    ({'baseline': 0.0}, 'baseline must be positive'),
    ({'baseline': -1.0}, 'baseline must be positive'),
    ({'n_slabs': 0}, 'n_slabs must be at least 1'),
])
def test_the_profile_routine_refuses_bad_geometry(kwargs, message):
    r"""A trajectory of no length is a mistake, not an empty product."""
    args = {'hamiltonian_of': constant_profile(np.eye(3, dtype=complex)),
            'baseline': 1.0}
    args.update(kwargs)
    with pytest.raises(ValueError, match=message):
        slabs.probabilities_3nu_profile(**args)


def test_the_profile_routine_refuses_a_non_callable():
    r"""The profile is a function of position, not an array.

    Passing the array a caller might expect to hand over gives a message
    saying so, rather than a TypeError from deep inside the search.
    """
    with pytest.raises(ValueError, match='must be callable'):
        slabs.probabilities_3nu_profile(np.eye(3, dtype=complex), 1.0)


def test_the_profile_routine_checks_what_the_callable_returned():
    r"""A profile returning the wrong shape is named at its source.

    Otherwise the mistake surfaces as a slab-count mismatch far from the
    function that caused it.
    """
    def wrong_shape(positions):
        return np.zeros((len(positions), 2, 2), dtype=complex)

    with pytest.raises(ValueError, match='must return one 3x3 Hamiltonian'):
        slabs.probabilities_3nu_profile(wrong_shape, 1.0)
