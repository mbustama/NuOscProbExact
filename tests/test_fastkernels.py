# -*- coding: utf-8 -*-
r"""Tests of the optional Numba backend and of the fallback path.

:mod:`fastkernels` compiles the batched expansions when Numba is
installed.  The point of these tests is that the two paths are
interchangeable: whichever is in use, the answers must be the same to
round-off, and turning the compiled one off must change nothing but the
speed.

Every equivalence test here runs with ``fastkernels.USE_NUMBA`` forced
off as well as on, so the NumPy path is exercised even on a machine
where Numba is installed --- otherwise it would rot unnoticed.
"""

import numpy as np
import pytest
import scipy.linalg

import globaldefs as gd
import hamiltonians4nu

import fastkernels
import oscprob2nu
import oscprob3nu
import oscprob4nu
import slabs

from conftest import ATOL, as_nested_list, random_hermitian


needs_numba = pytest.mark.skipif(not fastkernels.HAVE_NUMBA,
                                 reason='Numba is not installed')


# The stiff 3+1 spectrum, which is the four-flavor case that separates a
# faithful kernel from a plausible one.  An eV-scale sterile puts four
# orders of magnitude between the eigenvalues and clusters three of them.
STIFF_BASELINE = 1300.0*gd.CONV_KM_TO_INV_EV
STIFF_31 = np.asarray(hamiltonians4nu.hamiltonian_4nu_matter(
    hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(
        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF,
        np.sqrt(0.10), np.sqrt(0.10), 0.0, gd.DCP_NO_BF,
        gd.D21_NO_BF, gd.D31_NO_BF, 1.0),
    1.0e9, gd.VCC_EARTH_CRUST, gd.VNC_EARTH_CRUST))


def stiff_reference():
    r"""Returns the sixteen probabilities for `STIFF_31`, from ``expm``.

    An independent computation, so that the two paths are compared
    against something outside both rather than only against each other.
    """
    traceless = STIFF_31 - np.trace(STIFF_31).real/4.0*np.eye(4)
    operator = scipy.linalg.expm(-1.j*traceless*STIFF_BASELINE)

    return (np.abs(operator.T)**2).reshape(16)


# --------------------------------------------------------------------------
# The switch itself
# --------------------------------------------------------------------------

def test_available_follows_the_switch(monkeypatch):
    r"""available() is the conjunction of the import and the switch."""
    monkeypatch.setattr(fastkernels, 'USE_NUMBA', False)
    assert fastkernels.available() is False
    monkeypatch.setattr(fastkernels, 'USE_NUMBA', True)
    assert fastkernels.available() is fastkernels.HAVE_NUMBA


def test_library_imports_without_numba(monkeypatch):
    r"""With Numba reported absent, the library still works.

    This is the state of an installation without the optional extra.
    """
    monkeypatch.setattr(fastkernels, 'HAVE_NUMBA', False)
    assert fastkernels.available() is False

    h_stack = np.stack([np.diag([1.0, 2.0, -3.0]).astype(complex)]*40)
    prob = oscprob3nu.probabilities_3nu(h_stack, 1.5)
    assert prob.shape == (40, 9)
    assert np.allclose(prob.reshape(-1, 3, 3).sum(axis=-1), 1.0, atol=ATOL)


# --------------------------------------------------------------------------
# The two paths agree
# --------------------------------------------------------------------------

@needs_numba
@pytest.mark.parametrize('n', [1, 5, 17, 64, 300, 1000])
def test_3nu_kernel_matches_numpy_path(rng, monkeypatch, n, kernel_spy):
    r"""Compiled and NumPy paths agree for stacks of every size.

    The sizes straddle both thresholds: `oscprob3nu.SMALL_BATCH`, below
    which the NumPy path defers to the scalar one, and
    `fastkernels.PARALLEL_THRESHOLD`, above which the kernel spreads
    over cores.
    """
    h_stack = np.stack([random_hermitian(rng, 3) for _ in range(n)])
    l_stack = rng.uniform(0.1, 50.0, n)

    monkeypatch.setattr(fastkernels, 'USE_NUMBA', True)
    with_numba = oscprob3nu.probabilities_3nu(h_stack, l_stack)
    assert kernel_spy['probabilities_3nu_kernel'] == 1

    monkeypatch.setattr(fastkernels, 'USE_NUMBA', False)
    without = oscprob3nu.probabilities_3nu(h_stack, l_stack)

    assert with_numba.shape == without.shape == (n, 9)
    assert np.allclose(with_numba, without, atol=1.0e-12)


@needs_numba
@pytest.mark.parametrize('n', [1, 17, 300, 1000])
def test_2nu_kernel_matches_numpy_path(rng, monkeypatch, n, kernel_spy):
    r"""Compiled and NumPy paths agree, for two flavors.

    `fastkernels.MIN_BATCH` keeps the kernel away from two-flavor stacks
    below fifty thousand elements, where NumPy is quicker.  Test stacks
    that size would be unwieldy, so the threshold is lowered here --- and
    the spy asserts the kernel really was reached, because without that
    this test silently compared the NumPy path against itself.
    """
    h_stack = np.stack([random_hermitian(rng, 2) for _ in range(n)])
    l_stack = rng.uniform(0.1, 50.0, n)

    monkeypatch.setattr(fastkernels, 'USE_NUMBA', True)
    monkeypatch.setattr(fastkernels, 'MIN_BATCH', {2: 1, 3: 1})
    with_numba = oscprob2nu.probabilities_2nu(h_stack, l_stack)
    assert kernel_spy['probabilities_2nu_kernel'] == 1

    monkeypatch.setattr(fastkernels, 'USE_NUMBA', False)
    without = oscprob2nu.probabilities_2nu(h_stack, l_stack)

    assert np.allclose(with_numba, without, atol=ATOL)


@needs_numba
@pytest.mark.parametrize('n', [1, 5, 17, 64, 300, 1000])
def test_4nu_kernel_matches_numpy_path(rng, monkeypatch, n, kernel_spy):
    r"""Compiled and NumPy paths agree, for four flavors.

    The tolerance is looser than the three-flavor one above, and that is
    a property of the expansion rather than of the backend: see
    `test_4nu_stiff_spectrum_is_as_close_as_it_can_be`.  For a generic
    spectrum the two paths still land within a few ulp of each other.
    """
    h_stack = np.stack([random_hermitian(rng, 4) for _ in range(n)])
    l_stack = rng.uniform(0.1, 50.0, n)

    monkeypatch.setattr(fastkernels, 'USE_NUMBA', True)
    with_numba = oscprob4nu.probabilities_4nu(h_stack, l_stack)
    assert kernel_spy['probabilities_4nu_kernel'] == 1

    monkeypatch.setattr(fastkernels, 'USE_NUMBA', False)
    without = oscprob4nu.probabilities_4nu(h_stack, l_stack)

    assert with_numba.shape == without.shape == (n, 16)
    assert np.allclose(with_numba, without, atol=1.0e-12)


@needs_numba
def test_4nu_kernel_polishes_its_roots(monkeypatch, kernel_spy):
    r"""The compiled path refines the latent roots, as the NumPy one does.

    ``POLISH_ROOTS`` is module-level state that an ``@njit`` function
    cannot read at call time, so it is threaded through as an argument.
    If that thread were ever cut the kernel would still return plausible
    numbers --- just five hundred times worse on a stiff spectrum, which
    is precisely the case the refinement exists for.  This pins both
    settings, so a kernel that silently ignored the flag would fail
    whichever way the flag was set.
    """
    monkeypatch.setattr(fastkernels, 'USE_NUMBA', True)
    h_stack = np.stack([STIFF_31]*4)
    l_stack = np.full(4, STIFF_BASELINE)

    monkeypatch.setattr(oscprob4nu, 'POLISH_ROOTS', True)
    polished = oscprob4nu.probabilities_4nu(h_stack, l_stack)
    monkeypatch.setattr(oscprob4nu, 'POLISH_ROOTS', False)
    unpolished = oscprob4nu.probabilities_4nu(h_stack, l_stack)
    assert kernel_spy['probabilities_4nu_kernel'] == 2

    reference = stiff_reference()
    close = np.max(np.abs(polished[0] - reference))
    far = np.max(np.abs(unpolished[0] - reference))

    assert close < 1.0e-9
    assert far > 10.0*close


@needs_numba
def test_evolution_operator_matches_between_the_paths(backend, rng):
    r"""The batched evolution operator agrees on both backends.

    The counterpart of the probability checks above, and the one that
    was missing: `slabs` and `earth` compose *operators* across adjacent
    slabs, so they never call `probabilities_3nu`.  Until
    `fastkernels.evolution_operator_3nu_kernel` existed they had no
    compiled path at all, and installing the optional extra bought an
    Earth crossing nothing.

    Running under the `backend` fixture is the point: with Numba present
    `worthwhile(3, size)` is true at every size, so the NumPy branch of
    `_evolution_operator_3nu_batch` is unreachable in that configuration
    and would otherwise be covered by nothing.
    """
    h_stack = np.stack([random_hermitian(rng, 3) for _ in range(64)])
    l_stack = rng.uniform(0.1, 5.0, size=64)

    u = np.asarray(oscprob3nu.evolution_operator_3nu(h_stack, l_stack))
    assert u.shape == (64, 3, 3)

    # Against the scalar path, element by element.
    for k in range(64):
        one = np.asarray(oscprob3nu.evolution_operator_3nu(
            as_nested_list(h_stack[k]), float(l_stack[k])))
        assert np.allclose(u[k], one, atol=ATOL)

    # And unitary, on whichever backend ran.
    assert np.allclose(np.einsum('nij,nkj->nik', u, u.conj()),
                       np.eye(3), atol=1.0e-10)


def test_evolution_operator_kernel_is_what_slabs_reaches(monkeypatch,
                                                         kernel_spy):
    r"""A slab composition enters the compiled operator kernel.

    Guards the wiring rather than the arithmetic: it is the dispatch in
    `_evolution_operator_3nu_batch` that was missing, and a kernel that
    exists but is never called is the bug this replaced.
    """
    if not fastkernels.HAVE_NUMBA:
        pytest.skip('needs the optional Numba backend')
    monkeypatch.setattr(fastkernels, 'USE_NUMBA', True)

    h = np.stack([np.diag([1.0, 2.0, -3.0]).astype(complex)]*12)
    widths = np.full(12, 0.3)
    slabs.probabilities_3nu_slabs(h, widths)

    assert kernel_spy['evolution_operator_3nu_kernel'] == 1, (
        'the slab path did not reach the compiled evolution-operator kernel')


def test_kernel_matches_the_scalar_path(rng, monkeypatch):
    r"""The compiled kernel agrees with the scalar routine element by
    element, which is the reference both batched paths answer to."""
    monkeypatch.setattr(fastkernels, 'USE_NUMBA', True)
    h_stack = np.stack([random_hermitian(rng, 3) for _ in range(40)])
    l_stack = rng.uniform(0.1, 50.0, 40)

    batched = oscprob3nu.probabilities_3nu(h_stack, l_stack)
    scalar = np.array([oscprob3nu.probabilities_3nu(as_nested_list(h), l)
                       for h, l in zip(h_stack, l_stack)])
    assert np.allclose(batched, scalar, atol=1.0e-12)


# --------------------------------------------------------------------------
# Everything the NumPy path guarantees, the kernel must too
# --------------------------------------------------------------------------

def test_broadcasting_patterns(backend, rng):
    r"""All three broadcasting patterns work on whichever backend."""
    # Deliberately larger than oscprob3nu.SMALL_BATCH, so that the array
    # path is reached rather than the short-stack shortcut
    h_stack = np.stack([random_hermitian(rng, 3) for _ in range(13)])
    l_stack = np.linspace(0.5, 4.0, 11)

    assert oscprob3nu.probabilities_3nu(h_stack, 2.0).shape == (13, 9)
    assert oscprob3nu.probabilities_3nu(
        as_nested_list(h_stack[0]), l_stack).shape == (11, 9)
    grid = oscprob3nu.probabilities_3nu(h_stack[:, None, :, :],
                                        l_stack[None, :])
    assert grid.shape == (13, 11, 9)

    for i in range(13):
        for j in range(11):
            expected = oscprob3nu.probabilities_3nu(
                as_nested_list(h_stack[i]), l_stack[j])
            assert np.allclose(grid[i, j], expected, atol=1.0e-12)


@needs_numba
def test_4nu_stiff_spectrum_is_as_close_as_it_can_be(monkeypatch):
    r"""On a stiff 3+1 spectrum the two paths part company at 1e-10.

    That is not a defect in either of them, and the number is not a
    fudge: the Newton-form reconstruction of :math:`U_4` cancels by a
    factor of about :math:`10^6` on this spectrum --- its four terms
    reach :math:`9.7 \times 10^5` and sum to a matrix of modulus one ---
    so a last-bit difference anywhere in that sum arrives at
    :math:`10^{-10}`, and the two paths reorder the arithmetic by
    construction.  Measured, the amplification predicts 2.3e-10 and the
    two paths differ by 2.4e-10.

    Asserting round-off agreement here would therefore be asserting
    something false about the mathematics.  What *can* be required is
    that both paths sit inside the accuracy
    :data:`oscprob4nu.POLISH_ROOTS` claims, against a reference outside
    both, and that neither is quietly the worse of the two by an order
    of magnitude --- which is what would show if the kernel had dropped
    the refinement or botched the quartic.
    """
    h_stack = np.stack([STIFF_31]*8)
    l_stack = np.full(8, STIFF_BASELINE)
    reference = stiff_reference()

    monkeypatch.setattr(fastkernels, 'USE_NUMBA', True)
    with_numba = oscprob4nu.probabilities_4nu(h_stack, l_stack)
    monkeypatch.setattr(fastkernels, 'USE_NUMBA', False)
    without = oscprob4nu.probabilities_4nu(h_stack, l_stack)

    assert np.allclose(with_numba, without, atol=1.0e-9)

    kernel_error = np.max(np.abs(with_numba[0] - reference))
    numpy_error = np.max(np.abs(without[0] - reference))
    assert kernel_error < 1.0e-9
    assert numpy_error < 1.0e-9
    assert kernel_error < 100.0*numpy_error


DEGENERATE_4NU = [
    [2.5, 2.5, 2.5, 2.5],             # proportional to the identity
    [0.0, 0.0, 0.0, 0.0],             # the zero Hamiltonian
    [2.0, 2.0, -2.0, -2.0],           # two degenerate pairs
    [3.0, -1.0, -1.0, -1.0],          # a triply degenerate root
    [1.0, 1.0, 2.0, -4.0],            # one degenerate pair
    [1.0, 1.0 + 1.e-13, 2.0, -4.0],   # nearly, but not exactly, degenerate
    [1.0e6, 1.0e6, 1.0e6, 1.0e6],     # degenerate and large
]


@pytest.mark.parametrize('baseline', [1.0, 37.0])
def test_4nu_degenerate_spectra_on_both_backends(backend, rng, baseline):
    r"""Repeated latent roots survive whichever backend evaluates them.

    The kernel carries its own copy of the confluent divided
    differences, so this is not covered by the NumPy path's tests.  Half
    of these spectra are the ones that made the rejected Vandermonde
    solve singular; a kernel that reintroduced it for speed would fail
    here rather than in production.
    """
    degenerate = [np.diag(spectrum).astype(complex)
                  for spectrum in DEGENERATE_4NU]
    ordinary = [random_hermitian(rng, 4) for _ in range(20)]
    h_stack = np.stack(degenerate + ordinary)
    l_stack = np.full(len(h_stack), baseline)

    prob = oscprob4nu.probabilities_4nu(h_stack, l_stack)

    assert np.all(np.isfinite(prob))
    assert np.allclose(prob.reshape(-1, 4, 4).sum(axis=-1), 1.0, atol=ATOL)

    for index, matrix in enumerate(degenerate):
        traceless = matrix - np.trace(matrix).real/4.0*np.eye(4)
        reference = np.abs(scipy.linalg.expm(
            -1.j*traceless*baseline).T)**2
        assert np.allclose(prob[index], reference.reshape(16), atol=1.0e-11)


@needs_numba
def test_4nu_nearly_degenerate_pairs_agree_between_the_paths(monkeypatch):
    r"""The two paths agree where the Newton step is refused.

    This is the case that made the guard necessary, and it is the one
    place where the kernel and the NumPy path could not simply be
    expected to track each other: the resolvent root that vanishes for a
    degenerate pair is computed to a different last bit by each, and a
    square root near zero turns that into the pair's entire separation.
    One path would then see two identical roots and the other two roots
    a hair apart, and the unguarded step treated those two situations
    completely differently --- the paths disagreed by 186 on quantities
    that cannot exceed one.

    With the step refused whenever it would cross a neighbour, both
    paths decline it together and the disagreement is round-off again.
    """
    monkeypatch.setattr(fastkernels, 'USE_NUMBA', True)
    rng = np.random.default_rng(20240202)

    matrices = []
    for _ in range(120):
        separation = 10.0**rng.uniform(-16.0, -6.0)
        base = rng.normal()
        spectrum = np.array([base, base + separation,
                             base + rng.uniform(0.5, 2.0), 0.0])
        spectrum[3] = -spectrum[:3].sum()
        a = rng.standard_normal((4, 4)) + 1.j*rng.standard_normal((4, 4))
        unitary, _ = np.linalg.qr(a)
        rotated = (unitary @ np.diag(spectrum).astype(complex)
                   @ unitary.conj().T)
        matrices.append(0.5*(rotated + rotated.conj().T))

    h_stack = np.stack(matrices)
    l_stack = rng.uniform(0.1, 50.0, len(matrices))

    with_numba = oscprob4nu.probabilities_4nu(h_stack, l_stack)
    monkeypatch.setattr(fastkernels, 'USE_NUMBA', False)
    without = oscprob4nu.probabilities_4nu(h_stack, l_stack)

    assert np.all(np.isfinite(with_numba))
    assert np.allclose(with_numba.reshape(-1, 4, 4).sum(axis=-1), 1.0,
                       atol=1.0e-9)
    assert np.allclose(with_numba, without, atol=1.0e-9)


def test_4nu_broadcasting_patterns(backend, rng):
    r"""All three broadcasting patterns work on whichever backend."""
    h_stack = np.stack([random_hermitian(rng, 4) for _ in range(13)])
    l_stack = np.linspace(0.5, 4.0, 11)

    assert oscprob4nu.probabilities_4nu(h_stack, 2.0).shape == (13, 16)
    assert oscprob4nu.probabilities_4nu(
        as_nested_list(h_stack[0]), l_stack).shape == (11, 16)
    grid = oscprob4nu.probabilities_4nu(h_stack[:, None, :, :],
                                        l_stack[None, :])
    assert grid.shape == (13, 11, 16)

    for i in range(13):
        for j in range(11):
            expected = oscprob4nu.probabilities_4nu(
                as_nested_list(h_stack[i]), l_stack[j])
            assert np.allclose(grid[i, j], expected, atol=1.0e-12)


def test_4nu_probabilities_are_normalized(backend, rng):
    r"""Normalization and positivity hold on whichever backend."""
    h_stack = np.stack([random_hermitian(rng, 4) for _ in range(400)])
    l_stack = rng.uniform(0.1, 500.0, 400)

    prob = oscprob4nu.probabilities_4nu(h_stack, l_stack).reshape(400, 4, 4)
    assert np.allclose(prob.sum(axis=-1), 1.0, atol=1.0e-12)
    assert np.allclose(prob.sum(axis=-2), 1.0, atol=1.0e-12)
    assert np.all(prob >= -ATOL) and np.all(prob <= 1.0+ATOL)


def test_4nu_empty_and_single_element_stacks(backend):
    r"""Degenerate stack sizes behave on whichever backend."""
    empty = oscprob4nu.probabilities_4nu(np.zeros((0, 4, 4), dtype=complex),
                                         np.zeros(0))
    assert empty.shape == (0, 16)

    one = oscprob4nu.probabilities_4nu(np.eye(4, dtype=complex)[None], 1.0)
    assert one.shape == (1, 16)
    assert np.allclose(one.reshape(4, 4), np.eye(4), atol=ATOL)


def test_4nu_zero_baseline(backend, rng):
    r"""L = 0 gives the identity on whichever backend."""
    h_stack = np.stack([random_hermitian(rng, 4) for _ in range(30)])
    prob = oscprob4nu.probabilities_4nu(h_stack, np.zeros(30))
    assert np.allclose(prob.reshape(30, 4, 4), np.eye(4), atol=ATOL)


def test_4nu_scalar_call_never_reaches_the_kernel(backend, rng, kernel_spy):
    r"""A single Hamiltonian and baseline keeps returning a tuple.

    Four flavors is the one count where the scalar entry point shares
    the batched machinery, so the dispatch has to exclude it explicitly
    rather than fall out of a separate code path.  Were it not excluded,
    a scalar call would return an array of sixteen and ``len(prob)``
    would still be sixteen --- the shape of the bug is that nothing
    obvious breaks.
    """
    matrix = as_nested_list(random_hermitian(rng, 4))

    prob = oscprob4nu.probabilities_4nu(matrix, 1.5)

    assert isinstance(prob, tuple)
    assert len(prob) == 16
    assert kernel_spy['probabilities_4nu_kernel'] == 0


DEGENERATE = [np.eye(3, dtype=complex),
              np.zeros((3, 3), dtype=complex),
              np.diag([1.0, 1.0, -2.0]).astype(complex),
              np.diag([2.0, -1.0, -1.0]).astype(complex),
              1.0e6*np.eye(3, dtype=complex)]


def test_degenerate_hamiltonians(backend, rng):
    r"""Degenerate spectra survive whichever backend evaluates them.

    The compiled kernel carries its own copy of the degeneracy branch,
    so this is not covered by the NumPy path's tests.
    """
    ordinary = [random_hermitian(rng, 3) for _ in range(20)]
    h_stack = np.stack(DEGENERATE + ordinary)
    l_stack = np.full(len(h_stack), 3.0)

    prob = oscprob3nu.probabilities_3nu(h_stack, l_stack)
    assert np.all(np.isfinite(prob))
    assert np.allclose(prob.reshape(-1, 3, 3).sum(axis=-1), 1.0, atol=1.0e-12)

    scalar = np.array([oscprob3nu.probabilities_3nu(as_nested_list(h), 3.0)
                       for h in h_stack])
    assert np.allclose(prob, scalar, atol=1.0e-12)


def test_probabilities_are_normalized(backend, rng):
    r"""Normalization and positivity hold on whichever backend."""
    h_stack = np.stack([random_hermitian(rng, 3) for _ in range(400)])
    l_stack = rng.uniform(0.1, 500.0, 400)

    prob = oscprob3nu.probabilities_3nu(h_stack, l_stack).reshape(400, 3, 3)
    assert np.allclose(prob.sum(axis=-1), 1.0, atol=1.0e-12)
    assert np.allclose(prob.sum(axis=-2), 1.0, atol=1.0e-12)
    assert np.all(prob >= -ATOL) and np.all(prob <= 1.0+ATOL)


def test_empty_and_single_element_stacks(backend):
    r"""Degenerate stack sizes behave on whichever backend."""
    empty = oscprob3nu.probabilities_3nu(np.zeros((0, 3, 3), dtype=complex),
                                         np.zeros(0))
    assert empty.shape == (0, 9)

    one = oscprob3nu.probabilities_3nu(np.eye(3, dtype=complex)[None], 1.0)
    assert one.shape == (1, 9)
    assert np.allclose(one.reshape(3, 3), np.eye(3), atol=ATOL)


def test_zero_baseline(backend, rng):
    r"""L = 0 gives the identity on whichever backend."""
    h_stack = np.stack([random_hermitian(rng, 3) for _ in range(30)])
    prob = oscprob3nu.probabilities_3nu(h_stack, np.zeros(30))
    assert np.allclose(prob.reshape(30, 3, 3), np.eye(3), atol=ATOL)


def test_small_batch_matches_large_batch(backend, rng):
    r"""The small-stack shortcut returns what the array path would.

    `oscprob3nu.SMALL_BATCH` sends short stacks through the scalar
    routine; padding the same Hamiltonians out past the threshold must
    not change the answers for the ones in common.
    """
    h_stack = np.stack([random_hermitian(rng, 3) for _ in range(64)])
    l_stack = rng.uniform(0.1, 20.0, 64)

    short = oscprob3nu.probabilities_3nu(h_stack[:5], l_stack[:5])
    long = oscprob3nu.probabilities_3nu(h_stack, l_stack)
    assert np.allclose(short, long[:5], atol=1.0e-12)


@needs_numba
def test_chi_4nu_matches_numpy_and_pivots_when_it_must():
    r"""The kernel's determinant is its own unit, and tested as one.

    It has to be.  The determinant feeds the Newton refinement, and the
    refinement's accuracy does not survive into the probabilities --- the
    reconstruction's own conditioning swamps it --- so no end-to-end test
    can tell a good determinant from a poor one.  Mutation testing
    confirmed that: disabling the partial pivoting entirely left the
    whole suite green.  Checking it here, against
    :func:`numpy.linalg.det`, is what makes the choice of algorithm
    defended by a test rather than only by a docstring.

    The second half is the case that pivoting exists for.  Evaluating
    ``chi`` at ``psi`` equal to a diagonal entry makes the leading entry
    of ``psi*1 - H~`` exactly zero, and elimination without pivoting
    divides by it.
    """
    rng = np.random.default_rng(31415)
    scratch = np.empty((4, 4), dtype=np.complex128)

    for _ in range(200):
        traceless = random_hermitian(rng, 4)
        traceless -= np.trace(traceless).real/4.0*np.eye(4)
        traceless = np.ascontiguousarray(traceless)
        psi = float(rng.uniform(-2.0, 2.0))

        ours = fastkernels._chi_4nu(traceless, psi, scratch)
        reference = np.linalg.det(psi*np.eye(4) - traceless).real
        assert np.isclose(ours, reference, rtol=1.0e-10, atol=1.0e-14)

    # A leading entry that vanishes at the evaluation point
    singular = np.array([[0.5, 1.0, 0.0, 0.0],
                         [1.0, -0.5, 0.0, 0.0],
                         [0.0, 0.0, 0.25, 0.75],
                         [0.0, 0.0, 0.75, -0.25]], dtype=complex)
    assert np.isclose(np.trace(singular).real, 0.0)

    ours = fastkernels._chi_4nu(singular, 0.5, scratch)
    reference = np.linalg.det(0.5*np.eye(4) - singular).real

    assert np.isfinite(ours)
    assert not np.isclose(reference, 0.0), (
        'the pivoting case must have a non-vanishing determinant, or a '
        'routine that gave up and returned zero would pass')
    assert np.isclose(ours, reference, rtol=1.0e-12)


@needs_numba
def test_kernels_are_cached_on_disk():
    r"""The kernels are declared with cache=True.

    Without it every interpreter start pays seconds of compilation,
    which would make the optional backend a poor trade.
    """
    from numba.core.caching import NullCache
    for name in ('_run_3nu_serial', '_run_3nu_parallel',
                 '_run_2nu_serial', '_run_2nu_parallel',
                 '_run_4nu_serial', '_run_4nu_parallel'):
        dispatcher = getattr(fastkernels, name)
        assert not isinstance(dispatcher._cache, NullCache), name


# --------------------------------------------------------------------------
# The kernel is only used where it is faster
# --------------------------------------------------------------------------

def test_worthwhile_respects_the_measured_thresholds(monkeypatch):
    r"""The kernel is declined for stacks it would not speed up.

    A backend that is sometimes slower than the path it replaces is
    worse than no backend.  The two-flavor expansion reduces to a square
    root and a sine per element, which NumPy already does about as well
    as compiled code can, so below `fastkernels.MIN_BATCH` the NumPy
    path is kept.
    """
    if not fastkernels.HAVE_NUMBA:
        pytest.skip('Numba is not installed')
    monkeypatch.setattr(fastkernels, 'USE_NUMBA', True)

    assert fastkernels.worthwhile(3, 1) is True
    assert fastkernels.worthwhile(2, 1) is False
    assert fastkernels.worthwhile(2, fastkernels.MIN_BATCH[2]-1) is False
    assert fastkernels.worthwhile(2, fastkernels.MIN_BATCH[2]) is True

    # Four flavors was measured to win at every size, like three, and
    # unlike three it has no short-stack shortcut to lose to
    assert fastkernels.worthwhile(4, 1) is True


def test_worthwhile_is_false_when_unavailable(monkeypatch):
    r"""With the backend off, no stack is worth compiling for."""
    monkeypatch.setattr(fastkernels, 'USE_NUMBA', False)
    assert fastkernels.worthwhile(3, 10**9) is False
    assert fastkernels.worthwhile(2, 10**9) is False


def test_declining_the_kernel_does_not_change_the_answers(rng, monkeypatch):
    r"""A stack below the two-flavor threshold gives the same numbers as
    one above it."""
    if not fastkernels.HAVE_NUMBA:
        pytest.skip('Numba is not installed')
    monkeypatch.setattr(fastkernels, 'USE_NUMBA', True)
    assert not fastkernels.worthwhile(2, 100)

    h_stack = np.stack([random_hermitian(rng, 2) for _ in range(100)])
    l_stack = rng.uniform(0.1, 20.0, 100)

    declined = oscprob2nu.probabilities_2nu(h_stack, l_stack)
    monkeypatch.setattr(fastkernels, 'MIN_BATCH', {2: 1, 3: 1})
    forced = oscprob2nu.probabilities_2nu(h_stack, l_stack)
    assert np.allclose(declined, forced, atol=ATOL)


def test_fixed_hamiltonian_scan_on_both_backends(backend, rng):
    r"""One Hamiltonian over many baselines, on whichever backend.

    On the NumPy path this is the branch that forms the sum over the
    latent roots as a small matrix product rather than three broadcasts,
    and it is also the shape whose axis alignment needed the padding in
    `_u_coefficients_3nu_batch`.  With Numba installed the dispatch
    sends it to the kernel instead, so without running this on both
    backends the NumPy branch would go untested in a default run.
    """
    h = random_hermitian(rng, 3)
    l_stack = np.linspace(0.1, 30.0, 40)

    batched = oscprob3nu.probabilities_3nu(as_nested_list(h), l_stack)
    scalar = np.array([oscprob3nu.probabilities_3nu(as_nested_list(h), l)
                       for l in l_stack])

    assert batched.shape == (40, 9)
    assert np.allclose(batched, scalar, atol=1.0e-12)


def test_fixed_hamiltonian_grid_on_both_backends(backend, rng):
    r"""The same for a grid, whose baselines carry an extra axis."""
    h = random_hermitian(rng, 3)
    l_stack = np.linspace(0.1, 30.0, 12)

    grid = oscprob3nu.probabilities_3nu(as_nested_list(h), l_stack[None, :])
    assert grid.shape == (1, 12, 9)
    for j in range(12):
        expected = oscprob3nu.probabilities_3nu(as_nested_list(h), l_stack[j])
        assert np.allclose(grid[0, j], expected, atol=1.0e-12)


def test_2nu_fixed_hamiltonian_scan_on_both_backends(backend, rng):
    r"""The two-flavor equivalent, above its own short-stack threshold."""
    h = random_hermitian(rng, 2)
    l_stack = np.linspace(0.1, 30.0, 40)

    batched = oscprob2nu.probabilities_2nu(as_nested_list(h), l_stack)
    scalar = np.array([oscprob2nu.probabilities_2nu(as_nested_list(h), l)
                       for l in l_stack])

    assert batched.shape == (40, 4)
    assert np.allclose(batched, scalar, atol=ATOL)
