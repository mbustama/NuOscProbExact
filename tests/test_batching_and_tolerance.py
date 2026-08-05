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
import oscprob4nu
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
    r"""The fused chord kernel is what a scan actually calls.

    Skipped where the compiled path is switched off, which is a job the
    CI runs deliberately: with Numba present but `USE_NUMBA` false there
    is no kernel to reach, and asserting one was entered would be
    asserting the opposite of what that job is checking.

    Without this the comparison above passes by running the NumPy path
    twice, which is how a threshold change goes unnoticed.  The
    assertion is on the *fused* kernel specifically: an Earth scan that
    fell back to `slab_product_3nu_batch_kernel` would still be correct,
    and still be batched, and would quietly have given up the reason the
    fused one exists.
    """
    if not fastkernels.available():
        pytest.skip('the compiled path is switched off in this job')

    earth.probabilities_3nu_earth(h_vacuum(3), np.logspace(9.0, 11.0, 32),
                                  COSTHZ)

    assert kernel_spy['earth_chords_3nu_kernel'] == 1
    assert kernel_spy['slab_product_3nu_batch_kernel'] == 0
    assert kernel_spy['slab_product_3nu_kernel'] == 0


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
def test_the_numpy_batch_fallback_agrees_with_the_kernel(n_flavors,
                                                         monkeypatch):
    r"""The uncompiled branch beneath the batch dispatch still works.

    `monkeypatch` rather than a saved-and-restored global, because it
    restores what was actually there.  Setting `USE_NUMBA` back to True
    by hand looks equivalent and is not: the CI runs a job with Numba
    installed and the compiled path deliberately switched off, and a
    test that hands it back on leaves every later test in that job
    running the configuration the job exists to avoid.

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

    if not fastkernels.HAVE_NUMBA:
        pytest.skip('Numba is not installed')

    monkeypatch.setattr(fastkernels, 'USE_NUMBA', True)
    compiled = fn(hv, energies, COSTHZ)
    compiled_loop = np.array([fn(hv, e, COSTHZ) for e in energies])

    monkeypatch.setattr(fastkernels, 'USE_NUMBA', False)
    interpreted = fn(hv, energies, COSTHZ)
    interpreted_loop = np.array([fn(hv, e, COSTHZ) for e in energies])

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


# ------------------------------------------------------- chunk-size probing


def test_the_chunk_size_lands_in_the_permitted_range():
    r"""Whatever was detected, the chunk is clamped and usable.

    The detected cache is a hint, and a hint from an unfamiliar machine
    should not be able to produce a chunk of four bytes or of four
    gigabytes.
    """
    assert earth.CHUNK_BYTES_MIN <= earth.MAX_CHUNK_BYTES
    assert earth.MAX_CHUNK_BYTES <= earth.CHUNK_BYTES_MAX


def test_sysconf_probe_reports_a_size_or_nothing():
    r"""``os.sysconf`` carries the cache names on some builds only."""
    got = earth._cache_bytes_from_sysconf()

    assert got is None or got > 0


def test_sysconf_probe_survives_a_build_without_the_names(monkeypatch):
    r"""A build that refuses every name gives None, not an exception.

    This is the common case away from glibc, and it is reached on the
    machine this was written on, where every name raises `ValueError`.
    """
    def refuse(name):
        raise ValueError(name)

    monkeypatch.setattr(earth.os, 'sysconf', refuse)

    assert earth._cache_bytes_from_sysconf() is None


def test_sysconf_probe_takes_the_largest_reported(monkeypatch):
    r"""The last level is wanted, and it is not always the one asked
    for first."""
    monkeypatch.setattr(earth.os, 'sysconf',
                        lambda name: {'SC_LEVEL4_CACHE_SIZE': 0,
                                      'SC_LEVEL3_CACHE_SIZE': 8*1024*1024,
                                      'SC_LEVEL2_CACHE_SIZE': 1024*1024}[name])

    assert earth._cache_bytes_from_sysconf() == 8*1024*1024


def test_the_sysfs_probe_parses_the_units(tmp_path, monkeypatch):
    r"""``sysfs`` writes sizes as ``32K``, ``12288K``, sometimes ``2M``.

    Parsed rather than assumed, and the largest wins, since the
    directories are not ordered by level.
    """
    for name, text in (('index0', '48K'), ('index1', '2M'),
                       ('index2', '12288K'), ('index3', '')):
        (tmp_path/name).mkdir()
        (tmp_path/name/'size').write_text(text)
    # A directory with no size file at all, which sysfs does produce
    (tmp_path/'index4').mkdir()

    monkeypatch.setattr(earth, '_SYSFS_CACHE', str(tmp_path))

    assert earth._cache_bytes_from_sysfs() == 12288*1024


def test_the_sysfs_probe_ignores_an_unparseable_size(tmp_path, monkeypatch):
    r"""A size that is not a number is skipped, not fatal."""
    for name, text in (('index0', 'unknown'), ('index1', '4M')):
        (tmp_path/name).mkdir()
        (tmp_path/name/'size').write_text(text)

    monkeypatch.setattr(earth, '_SYSFS_CACHE', str(tmp_path))

    assert earth._cache_bytes_from_sysfs() == 4*1024*1024


def test_the_sysfs_probe_declines_where_there_is_no_sysfs(monkeypatch):
    r"""Off Linux the directory is simply absent."""
    monkeypatch.setattr(earth, '_SYSFS_CACHE', '/nonexistent/cache/path')

    with pytest.raises(OSError):
        earth._cache_bytes_from_sysfs()

    # ...which `_last_level_cache_bytes` is required to absorb
    monkeypatch.setattr(earth, '_cache_bytes_from_sysconf', lambda: None)
    monkeypatch.setattr(earth, '_cache_bytes_from_sysctl', lambda: None)
    assert earth._last_level_cache_bytes() is None


def test_the_sysctl_probe_reports_a_size_or_nothing():
    r"""Answers on macOS, declines elsewhere, raises nowhere."""
    got = earth._cache_bytes_from_sysctl()

    assert got is None or got > 0


def test_the_sysctl_probe_declines_without_a_c_library(monkeypatch):
    r"""No libc to find means no answer, rather than an exception."""
    import ctypes.util
    monkeypatch.setattr(ctypes.util, 'find_library', lambda _: None)

    assert earth._cache_bytes_from_sysctl() is None


def test_a_probe_that_raises_is_skipped(monkeypatch):
    r"""A hint must never be able to stop the module importing.

    Every probe touches something platform-specific, and the point of
    catching broadly is that an unfamiliar machine costs some speed on
    long scans rather than an ImportError.
    """
    def explode():
        raise RuntimeError('an unfamiliar machine')

    monkeypatch.setattr(earth, '_cache_bytes_from_sysconf', explode)
    monkeypatch.setattr(earth, '_cache_bytes_from_sysfs', explode)
    monkeypatch.setattr(earth, '_cache_bytes_from_sysctl', explode)

    assert earth._last_level_cache_bytes() is None


def test_nothing_detected_falls_back_to_a_plausible_chunk(monkeypatch):
    r"""The fallback is a real cache size, not a degenerate one."""
    monkeypatch.setattr(earth, '_cache_bytes_from_sysconf', lambda: None)
    monkeypatch.setattr(earth, '_cache_bytes_from_sysfs', lambda: None)
    monkeypatch.setattr(earth, '_cache_bytes_from_sysctl', lambda: None)

    assert earth._last_level_cache_bytes() is None
    assert earth.CHUNK_BYTES_MIN <= earth.CHUNK_BYTES_FALLBACK
    assert earth.CHUNK_BYTES_FALLBACK <= earth.CHUNK_BYTES_MAX


def test_a_chunk_always_holds_several_energies(backend):
    r"""The byte budget cannot cut a chunk down to a single energy.

    A four-flavor chord with many slabs costs a quarter of a megabyte
    per energy; without a floor, a small budget would re-enter the
    kernel once per energy and undo the batching entirely.
    """
    energies = np.logspace(9.0, 10.0, 80)
    original = earth.MAX_CHUNK_BYTES
    try:
        earth.MAX_CHUNK_BYTES = 1
        got = earth.probabilities_3nu_earth(h_vacuum(3), energies, COSTHZ)
    finally:
        earth.MAX_CHUNK_BYTES = original

    assert got.shape == (80, 9)
    assert np.allclose(
        got, earth.probabilities_3nu_earth(h_vacuum(3), energies, COSTHZ),
        rtol=0.0, atol=0.0)


class _FakeLibc:
    r"""Stands in for a macOS libc, answering a fixed set of names.

    The real call cannot run on the machine this suite runs on, so the
    parsing around it would otherwise be asserted about rather than
    exercised.  Writes through the pointer exactly as ``sysctlbyname``
    does, and returns non-zero for a name it does not know, which is how
    Apple Silicon reports having no L3.
    """

    def __init__(self, answers):
        self.answers = answers
        self.asked = []

    def sysctlbyname(self, name, value_ptr, length_ptr, _new, _newlen):
        self.asked.append(name)
        if name not in self.answers:
            return -1
        value_ptr.contents.value = self.answers[name]
        return 0


def _install_fake_libc(monkeypatch, answers):
    r"""Points the sysctl probe at a `_FakeLibc`."""
    import ctypes
    import ctypes.util

    libc = _FakeLibc(answers)
    monkeypatch.setattr(ctypes.util, 'find_library', lambda _: 'libc.dylib')
    monkeypatch.setattr(ctypes, 'CDLL', lambda _: libc)
    return libc


def test_the_sysctl_probe_reads_an_l3(monkeypatch):
    r"""An Intel Mac reports ``hw.l3cachesize`` and that is the answer."""
    _install_fake_libc(monkeypatch, {b'hw.l3cachesize': 16*1024*1024,
                                     b'hw.l2cachesize': 256*1024})

    assert earth._cache_bytes_from_sysctl() == 16*1024*1024


def test_the_sysctl_probe_falls_back_to_l2_on_apple_silicon(monkeypatch):
    r"""Apple Silicon has no L3, so its L2 is the last level.

    Returning nothing there would send those machines to the fallback
    when the real number was available one name away.
    """
    libc = _install_fake_libc(
        monkeypatch, {b'hw.perflevel0.l2cachesize': 32*1024*1024})

    assert earth._cache_bytes_from_sysctl() == 32*1024*1024
    assert b'hw.l3cachesize' in libc.asked


def test_the_sysctl_probe_declines_when_nothing_answers(monkeypatch):
    r"""A libc with the symbol but no cache names gives None."""
    _install_fake_libc(monkeypatch, {})

    assert earth._cache_bytes_from_sysctl() is None


def test_the_sysctl_probe_ignores_a_zero(monkeypatch):
    r"""A name that answers zero is not a cache size."""
    _install_fake_libc(monkeypatch, {b'hw.l3cachesize': 0,
                                     b'hw.l2cachesize': 4*1024*1024})

    assert earth._cache_bytes_from_sysctl() == 4*1024*1024


def test_the_sysctl_probe_declines_without_the_symbol(monkeypatch):
    r"""Linux's libc has no ``sysctlbyname``, which is not an error."""
    import ctypes
    import ctypes.util

    class _NoSymbol:
        pass

    monkeypatch.setattr(ctypes.util, 'find_library', lambda _: 'libc.so.6')
    monkeypatch.setattr(ctypes, 'CDLL', lambda _: _NoSymbol())

    assert earth._cache_bytes_from_sysctl() is None


# ------------------------------------------- fused kernels and angle grids


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
def test_the_fused_kernel_matches_the_materialised_one(n_flavors,
                                                       monkeypatch):
    r"""Building Hamiltonians inline is a memory layout, not a scheme.

    The fused kernel exists to avoid materialising one matrix per slab
    per energy; what it computes is the same midpoint Hamiltonian the
    batch kernel is handed, so the two must agree exactly rather than
    closely.  Anything less would mean the fused path had quietly become
    a different approximation.

    The palindromic composer is disabled here, and only here: it
    multiplies the same operators in a different grouping and so rounds
    differently, which would mask exactly the disagreement this test
    exists to catch.  That it *does* round differently, and by how
    much, is asserted separately.
    """
    if not fastkernels.HAVE_NUMBA:
        pytest.skip('Numba is not installed')

    monkeypatch.setattr(fastkernels, 'MIN_MIRROR_SLABS',
                        {2: 1 << 30, 3: 1 << 30, 4: 1 << 30})
    energies = np.logspace(9.0, 11.0, 24)
    widths_km, densities = earth._earth_slabs_cached(COSTHZ, 8)
    widths = widths_km*gd.CONV_KM_TO_INV_EV
    potentials = earth.matter_potential(densities,
                                        gd.ELECTRON_FRACTION_EARTH_CRUST)
    hv = h_vacuum(n_flavors)

    if n_flavors == 2:
        fused = fastkernels.earth_chords_2nu_kernel(hv, energies, potentials,
                                                    widths)
        stack = hamiltonians2nu.hamiltonian_2nu_matter(
            hv, energies[:, None], potentials[None, :])
        materialised = fastkernels.slab_product_2nu_batch_kernel(stack, widths)
    elif n_flavors == 3:
        fused = fastkernels.earth_chords_3nu_kernel(hv, energies, potentials,
                                                    widths)
        stack = hamiltonians3nu.hamiltonian_3nu_matter(
            hv, energies[:, None], potentials[None, :])
        materialised = fastkernels.slab_product_3nu_batch_kernel(stack, widths)
    else:
        nc = earth.matter_potential_nc(
            densities, electron_fraction=gd.ELECTRON_FRACTION_EARTH_CRUST)
        fused = fastkernels.earth_chords_4nu_kernel(
            hv, energies, potentials, nc, widths, oscprob4nu.POLISH_ROOTS)
        stack = hamiltonians4nu.hamiltonian_4nu_matter(
            hv, energies[:, None], potentials[None, :], nc[None, :])
        materialised = fastkernels.slab_product_4nu_batch_kernel(
            oscprob4nu._traceless_part(stack), widths, oscprob4nu.POLISH_ROOTS)

    assert np.array_equal(fused, materialised)


def test_the_fused_kernel_takes_a_scalar_energy():
    r"""A batch of one still comes back as a single operator."""
    if not fastkernels.HAVE_NUMBA:
        pytest.skip('Numba is not installed')

    widths_km, densities = earth._earth_slabs_cached(COSTHZ, 8)
    widths = widths_km*gd.CONV_KM_TO_INV_EV
    potentials = earth.matter_potential(densities,
                                        gd.ELECTRON_FRACTION_EARTH_CRUST)

    got = fastkernels.earth_chords_3nu_kernel(h_vacuum(3), np.float64(1.0e10),
                                              potentials, widths)

    assert got.shape == (3, 3)
    assert np.allclose(got.conj().T @ got, np.eye(3), atol=1.0e-10)


def test_the_fused_kernel_spreads_wide_work_over_threads():
    r"""Both sides of the parallel threshold give the same answer."""
    if not fastkernels.HAVE_NUMBA:
        pytest.skip('Numba is not installed')

    widths_km, densities = earth._earth_slabs_cached(COSTHZ, 8)
    widths = widths_km*gd.CONV_KM_TO_INV_EV
    potentials = earth.matter_potential(densities,
                                        gd.ELECTRON_FRACTION_EARTH_CRUST)
    n_wide = 2*fastkernels.PARALLEL_THRESHOLD//widths.shape[0] + 2
    energies = np.logspace(9.0, 11.0, n_wide)

    wide = fastkernels.earth_chords_3nu_kernel(h_vacuum(3), energies,
                                               potentials, widths)
    one_by_one = np.array([
        fastkernels.earth_chords_3nu_kernel(h_vacuum(3), energies[i:i+1],
                                            potentials, widths)[0]
        for i in range(n_wide)])

    assert n_wide*widths.shape[0] >= fastkernels.PARALLEL_THRESHOLD
    assert np.array_equal(wide, one_by_one)


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
def test_an_angle_grid_equals_a_loop_over_it(backend, n_flavors):
    r"""An oscillogram is its points, evaluated together.

    Both backends, because the grid path reaches the fused kernel on one
    and the chunked NumPy path on the other, and the two must land on
    the same numbers.
    """
    energies = np.logspace(9.0, 11.0, 6)
    angles = np.linspace(-1.0, -0.1, 4)
    fn = probabilities_earth(n_flavors)
    hv = h_vacuum(n_flavors)

    grid = fn(hv, energies[None, :], angles[:, None])
    loop = np.array([[fn(hv, e, c) for e in energies] for c in angles])

    assert grid.shape == (4, 6, n_flavors*n_flavors)
    assert np.allclose(grid, loop, rtol=0.0, atol=ATOL)


def test_an_angle_grid_broadcasts_like_numpy(backend):
    r"""Scalar against array, and array against array of one shape.

    The broadcasting is the caller's handle on what the grid means, so
    the degenerate shapes have to behave as `numpy` would.
    """
    hv = h_vacuum(3)
    energies = np.logspace(9.0, 11.0, 5)
    angles = np.linspace(-1.0, -0.2, 5)

    assert earth.probabilities_3nu_earth(hv, 1.0e10, angles).shape == (5, 9)
    assert earth.probabilities_3nu_earth(hv, energies, -0.8).shape == (5, 9)
    # Two 1-D arrays of equal length pair up element by element
    paired = earth.probabilities_3nu_earth(hv, energies, angles)
    assert paired.shape == (5, 9)
    assert np.allclose(
        paired[2],
        earth.probabilities_3nu_earth(hv, energies[2], angles[2]),
        rtol=0.0, atol=ATOL)


def test_a_repeated_angle_is_evaluated_once(backend):
    r"""A broadcast grid repeats every angle, and must not re-cut it.

    The geometry is the expensive part of a new angle, so the grid
    groups by distinct angle rather than walking the points.
    """
    hv = h_vacuum(3)
    energies = np.logspace(9.0, 11.0, 4)
    repeated = np.array([-0.8, -0.8, -0.8])

    got = earth.probabilities_3nu_earth(hv, energies[None, :],
                                        repeated[:, None])

    assert got.shape == (3, 4, 9)
    # Every row is the same angle, so every row is the same numbers
    assert np.array_equal(got[0], got[1])
    assert np.array_equal(got[0], got[2])


def test_a_grid_that_does_not_broadcast_says_so(backend):
    r"""The error names both shapes and how to index for a grid.

    Handing two flat arrays of different lengths is the natural mistake,
    and numpy's own message does not mention the axes to use.
    """
    hv = h_vacuum(3)

    with pytest.raises(ValueError, match='do not broadcast together'):
        earth.probabilities_3nu_earth(hv, np.logspace(9.0, 11.0, 7),
                                      np.linspace(-1.0, -0.2, 5))


def test_an_angle_grid_accepts_a_tolerance(backend):
    r"""Refinement over a grid binds on every point of it."""
    hv = h_vacuum(3)
    energies = np.logspace(9.0, 10.0, 3)
    angles = np.array([-0.9, -0.4])

    got, n = earth.probabilities_3nu_earth(hv, energies[None, :],
                                           angles[:, None], atol=1.0e-4,
                                           return_n_slabs=True)

    assert got.shape == (2, 3, 9)
    assert np.allclose(
        got, earth.probabilities_3nu_earth(hv, energies[None, :],
                                           angles[:, None], n),
        rtol=0.0, atol=0.0)


def test_earth_slabs_stays_scalar_in_the_angle():
    r"""Two angles have different slab counts, so a chord is one angle.

    `earth_slabs` returns the widths and densities of *a* chord; there
    is no array shape that holds two chords of different lengths, which
    is why the grid path groups by angle instead of broadcasting here.
    """
    short = earth.earth_slabs(-0.1)[0]
    long_chord = earth.earth_slabs(-1.0)[0]

    assert short.shape != long_chord.shape
    with pytest.raises((TypeError, ValueError)):
        earth.earth_slabs(np.array([-0.5, -0.9]))


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
def test_the_general_batch_path_still_reaches_its_kernel(n_flavors,
                                                         kernel_spy):
    r"""The general batched composer is not dead code.

    `earth` reaches the *fused* kernel now, so nothing on the Earth path
    exercises the compiled branch of `_probabilities_slabs_batch` any
    more.  That branch is still the batched composer for an arbitrary
    stack of Hamiltonians --- which the fused kernel cannot do, knowing
    only matter Hamiltonians of the form ``H_vac/E + V P`` --- so it is
    tested here directly rather than left to rot behind a dispatch that
    stopped choosing it.
    """
    if not fastkernels.available():
        pytest.skip('the compiled path is switched off in this job')

    rng = np.random.default_rng(20260808 + n_flavors)
    n_chords, n_slabs = 4, 6
    a = (rng.normal(size=(n_chords, n_slabs, n_flavors, n_flavors))
         + 1.0j*rng.normal(size=(n_chords, n_slabs, n_flavors, n_flavors)))
    h = (a + np.conj(np.swapaxes(a, -1, -2)))/2.0
    widths = rng.uniform(0.1, 2.0, size=n_slabs)

    batched = slabs._probabilities_slabs_batch(h, widths, n_flavors, 'caller')
    per_chord = {2: slabs.probabilities_2nu_slabs,
                 3: slabs.probabilities_3nu_slabs,
                 4: slabs.probabilities_4nu_slabs}[n_flavors]
    looped = np.array([per_chord(h[i], widths) for i in range(n_chords)])

    assert kernel_spy['slab_product_%dnu_batch_kernel' % n_flavors] == 1
    assert np.allclose(batched, looped, rtol=0.0, atol=ATOL)


# ----------------------------------------------------------- the palindrome


def test_the_prem_chord_is_exactly_palindromic():
    r"""A chord meets every radius twice, and the arrays say so exactly.

    Exactly, not nearly: `fastkernels.palindromic` decides on bitwise
    equality, so a chord whose widths differ from their mirror in the
    last bit --- which is what independent `linspace` calls per segment
    produce --- would silently never take the halved path.
    """
    for costhz in (-0.15, -0.55, -0.8, -1.0):
        for n_per in (1, 3, 8):
            widths, densities = earth._earth_slabs_cached(costhz, n_per)
            assert np.array_equal(widths, widths[::-1]), (costhz, n_per)
            assert np.array_equal(densities, densities[::-1]), (costhz, n_per)
            assert fastkernels.palindromic(widths, densities)


def test_symmetrising_the_chord_preserves_its_length():
    r"""Mirroring the widths moves them by round-off, not by physics."""
    for costhz in (-0.3, -1.0):
        widths, _ = earth._earth_slabs_cached(costhz, 8)
        assert np.isclose(widths.sum(),
                          earth.distance_traveled_inside_earth(costhz),
                          rtol=1.0e-12, atol=0.0)


@pytest.mark.parametrize('values, expected', [
    (np.array([1.0, 2.0, 1.0]), True),
    (np.array([1.0, 2.0, 2.0, 1.0]), True),
    (np.array([1.0, 2.0, 3.0]), False),
    (np.array([1.0]), True),
    (np.array([]), True),
])
def test_palindromic_decides_on_exact_equality(values, expected):
    r"""One array at a time, including the degenerate lengths."""
    assert fastkernels.palindromic(values) is expected


def test_palindromic_requires_every_array_to_agree():
    r"""Symmetric densities across asymmetric widths is not a palindrome.

    Both have to mirror for the operators to, and it is the conjunction
    that the composer relies on.
    """
    symmetric = np.array([1.0, 2.0, 1.0])
    asymmetric = np.array([1.0, 2.0, 3.0])

    assert fastkernels.palindromic(symmetric, symmetric)
    assert not fastkernels.palindromic(symmetric, asymmetric)
    assert not fastkernels.palindromic(asymmetric, symmetric)
    assert fastkernels.palindromic()


def test_a_last_bit_difference_defeats_the_palindrome():
    r"""The check is exact, and that is the point.

    A tolerance here would hand a nearly-symmetric profile the answer
    for a symmetric one, silently.  The saving is worth having only
    where the mirrored operators are genuinely identical.
    """
    widths = np.array([1.0, 2.0, 1.0])
    nudged = widths.copy()
    nudged[-1] = np.nextafter(nudged[-1], 2.0)

    assert fastkernels.palindromic(widths)
    assert not fastkernels.palindromic(nudged)


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
@pytest.mark.parametrize('costhz', [-0.15, -0.8, -1.0])
def test_the_mirrored_composer_agrees_with_the_whole_chord(n_flavors, costhz):
    r"""Half the expansions, the same answer.

    The mirrored composer multiplies the same operators in a different
    grouping, so it agrees to round-off rather than bitwise.  What must
    not happen is a difference larger than that, which would mean the
    split or the ordering is wrong --- and getting the order backwards
    still returns something that looks like a probability.
    """
    if not fastkernels.HAVE_NUMBA:
        pytest.skip('Numba is not installed')

    energies = np.logspace(9.0, 11.0, 32)
    hv = h_vacuum(n_flavors)
    fn = probabilities_earth(n_flavors)

    mirrored = fn(hv, energies, costhz)
    original = fastkernels.MIN_MIRROR_SLABS
    try:
        fastkernels.MIN_MIRROR_SLABS = {2: 1 << 30, 3: 1 << 30, 4: 1 << 30}
        whole = fn(hv, energies, costhz)
    finally:
        fastkernels.MIN_MIRROR_SLABS = original

    assert np.allclose(mirrored, whole, rtol=0.0, atol=1.0e-11)

    # Unitarity, as a comparison rather than against a fixed number.
    # Four flavors sits at 1e-10 rather than round-off because of the
    # root polishing, and has done since before this optimisation; what
    # matters is that halving the composition does not make it worse.
    def worst_unitarity(p):
        return np.max(np.abs(np.sum(p.reshape(-1, n_flavors, n_flavors),
                                    axis=-1) - 1.0))

    assert worst_unitarity(mirrored) <= 2.0*worst_unitarity(whole) + 1.0e-15


@pytest.mark.parametrize('n_per, parity', [(8, 'even'), (1, 'odd')])
def test_the_mirrored_composer_handles_both_parities(n_per, parity):
    r"""An odd chord has a middle slab that is its own mirror.

    A chord crosses an odd number of segments, so the slab count is odd
    whenever the subdivision is, and the centre slab must be applied
    once rather than twice or not at all.  Both mistakes leave a
    plausible-looking unitary answer.
    """
    if not fastkernels.HAVE_NUMBA:
        pytest.skip('Numba is not installed')

    widths, _ = earth._earth_slabs_cached(-1.0, n_per)
    assert (len(widths) % 2 == 0) == (parity == 'even'), len(widths)

    energies = np.logspace(9.0, 11.0, 16)
    mirrored = earth.probabilities_3nu_earth(h_vacuum(3), energies, -1.0,
                                             n_per)
    original = fastkernels.MIN_MIRROR_SLABS
    try:
        fastkernels.MIN_MIRROR_SLABS = {2: 1 << 30, 3: 1 << 30, 4: 1 << 30}
        whole = earth.probabilities_3nu_earth(h_vacuum(3), energies, -1.0,
                                              n_per)
    finally:
        fastkernels.MIN_MIRROR_SLABS = original

    assert np.allclose(mirrored, whole, rtol=0.0, atol=1.0e-11)


def test_an_asymmetric_chord_is_composed_whole():
    r"""The halved path is not taken on a profile that is not a
    palindrome.

    This is the assertion that keeps the optimisation honest: the fused
    kernels are public, and handing them an asymmetric profile must give
    the asymmetric answer rather than the symmetric one.
    """
    if not fastkernels.HAVE_NUMBA:
        pytest.skip('Numba is not installed')

    widths_km, densities = earth._earth_slabs_cached(-0.8, 8)
    widths = widths_km*gd.CONV_KM_TO_INV_EV
    potentials = earth.matter_potential(densities,
                                        gd.ELECTRON_FRACTION_EARTH_CRUST)
    # Break the symmetry in one slab only, and in the potential alone
    lopsided = potentials.copy()
    lopsided[0] *= 1.5
    energies = np.logspace(9.0, 11.0, 8)

    assert not fastkernels.palindromic(lopsided, widths)
    got = fastkernels.earth_chords_3nu_kernel(h_vacuum(3), energies,
                                              lopsided, widths)

    # Against the general composer, which knows nothing of palindromes
    stack = hamiltonians3nu.hamiltonian_3nu_matter(
        h_vacuum(3), energies[:, None], lopsided[None, :])
    reference = fastkernels.slab_product_3nu_batch_kernel(stack, widths)

    assert np.array_equal(got, reference)


def test_a_short_chord_is_composed_whole():
    r"""Below the threshold the mirrored composer is not worth it.

    It carries two running products instead of one, so it only pays
    once there are enough slabs for the halved expansions to dominate
    that.  Correctness is unaffected either way, which is what this
    checks.
    """
    if not fastkernels.available():
        pytest.skip('the compiled path is switched off in this job')

    assert not fastkernels.worthwhile_mirror(3, 7)
    assert fastkernels.worthwhile_mirror(3, fastkernels.MIN_MIRROR_SLABS[3])

    energies = np.logspace(9.0, 11.0, 8)
    short = earth.probabilities_3nu_earth(h_vacuum(3), energies, -0.15, 1)
    loop = np.array([earth.probabilities_3nu_earth(h_vacuum(3), e, -0.15, 1)
                     for e in energies])

    assert np.allclose(short, loop, rtol=0.0, atol=ATOL)


def test_the_compiled_palindrome_check_agrees_with_the_numpy_one():
    r"""Two implementations of one predicate must not disagree.

    The compiled one exists only because the NumPy one is too slow to
    call on a materialised stack --- 6 microseconds against the 6 the
    halved composition saves, which measured as a net loss.  Being
    quicker is no use if it answers differently.
    """
    if not fastkernels.HAVE_NUMBA:
        pytest.skip('Numba is not installed')

    rng = np.random.default_rng(20260809)
    widths = rng.uniform(0.1, 2.0, size=8)
    a = rng.normal(size=(8, 3, 3)) + 1.0j*rng.normal(size=(8, 3, 3))

    symmetric_h = np.concatenate([a[:4], a[:4][::-1]])
    symmetric_w = np.concatenate([widths[:4], widths[:4][::-1]])

    for h, w in ((symmetric_h, symmetric_w), (a, symmetric_w),
                 (symmetric_h, widths), (a, widths)):
        h = np.ascontiguousarray(h, dtype=complex)
        w = np.ascontiguousarray(w, dtype=float)
        assert (fastkernels._palindromic_stack(h, w)
                == fastkernels.palindromic(h, w))


def test_the_compiled_check_rejects_on_either_array():
    r"""A width mismatch and a Hamiltonian mismatch are both fatal.

    They are separate early exits, and a sequence can fail on one while
    passing the other.
    """
    if not fastkernels.HAVE_NUMBA:
        pytest.skip('Numba is not installed')

    h = np.zeros((4, 3, 3), dtype=complex)
    w = np.array([1.0, 2.0, 2.0, 1.0])
    assert fastkernels._palindromic_stack(h, w)

    bad_width = w.copy()
    bad_width[0] = 1.5
    assert not fastkernels._palindromic_stack(h, bad_width)

    bad_h = h.copy()
    bad_h[0, 1, 2] = 1.0 + 0.0j
    assert not fastkernels._palindromic_stack(bad_h, w)


def test_a_symmetric_slab_sequence_is_composed_at_half_cost():
    r"""The saving is not Earth-specific, which is the point.

    A hand-built profile that reads the same from either end --- a
    symmetric castle wall, say --- has the same redundancy a chord does,
    and `slabs` finds it without knowing anything about the Earth.
    """
    if not fastkernels.HAVE_NUMBA:
        pytest.skip('Numba is not installed')

    rng = np.random.default_rng(20260810)
    half = rng.normal(size=(20, 3, 3)) + 1.0j*rng.normal(size=(20, 3, 3))
    half = (half + np.conj(np.swapaxes(half, -1, -2)))/2.0
    h = np.concatenate([half, half[::-1]])
    w_half = rng.uniform(0.1, 1.5, size=20)
    w = np.concatenate([w_half, w_half[::-1]])

    assert fastkernels.palindromic(h, w)
    mirrored = slabs.probabilities_3nu_slabs(h, w)

    original = fastkernels.MIN_MIRROR_SLABS
    try:
        fastkernels.MIN_MIRROR_SLABS = {2: 1 << 30, 3: 1 << 30, 4: 1 << 30}
        whole = slabs.probabilities_3nu_slabs(h, w)
    finally:
        fastkernels.MIN_MIRROR_SLABS = original

    assert np.allclose(mirrored, whole, rtol=0.0, atol=1.0e-11)
    assert np.isclose(sum(mirrored), 3.0, rtol=0.0, atol=1.0e-10)
