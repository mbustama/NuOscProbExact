# -*- coding: utf-8 -*-
r"""Tests of the four-neutrino SU(4) expansion.

The independent computations these check against are the same ones used
at two and three flavors --- ``scipy.linalg.expm`` for the evolution
operator and ``numpy.linalg.eigh`` for the spectrum --- plus one that
only exists at four flavors: switching the three sterile mixing angles
off must reproduce :mod:`oscprob3nu` exactly, since a decoupled fourth
state cannot change the active block.
"""

import numpy as np
import pytest
import scipy.linalg

import globaldefs as gd
import hamiltonians3nu
import hamiltonians4nu
import oscprob3nu
import oscprob4nu


BASELINE = 1300.0*gd.CONV_KM_TO_INV_EV
ENERGY = 1.0e9
DM41 = 1.0


def random_hermitian(rng):
    r"""Returns a random Hermitian 4x4 matrix."""
    matrix = (rng.standard_normal((4, 4))
              + 1.j*rng.standard_normal((4, 4)))

    return 0.5*(matrix + matrix.conj().T)


def traceless(matrix):
    r"""Returns the traceless part of a 4x4 matrix."""
    return matrix - np.trace(matrix).real/4.0*np.eye(4)


def vacuum_4nu(s14=0.0, s24=0.0, s34=0.0):
    r"""Returns a 3+1 energy-independent vacuum Hamiltonian."""
    return hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(
        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, s14, s24, s34,
        gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF, DM41)


def test_generators_are_orthonormal_traceless_and_hermitian():
    r"""The fifteen generators obey the SU(4) normalization."""
    lam = oscprob4nu.LAMBDA_SU4

    assert lam.shape == (15, 4, 4)
    gram = np.einsum('aij,bji->ab', lam, lam).real/2.0
    assert np.allclose(gram, np.eye(15), atol=1.e-14)
    assert np.allclose(np.trace(lam, axis1=1, axis2=2), 0.0, atol=1.e-14)
    assert np.allclose(lam, np.conj(np.swapaxes(lam, 1, 2)), atol=0.0)


def test_su3_star_product_identity_fails_at_four_flavors():
    r"""``(h*h)*h = |h|^2 h/3`` is an accident of three flavors.

    The three-flavor expansion leans on it; the four-flavor one cannot,
    which is why ``((h*h)*h)_a`` enters as independent data.  If this
    ever started passing, the extra term would be redundant and the
    module would be carrying dead weight.
    """
    rng = np.random.default_rng(7)
    lam = oscprob4nu.LAMBDA_SU4

    coefficients = rng.standard_normal(15)
    matrix = np.einsum('a,aij->ij', coefficients, lam)
    invariant_2, _, _ = oscprob4nu.su4_invariants(matrix)

    cubed = matrix @ matrix @ matrix
    tower = (0.5*np.einsum('aij,ji->a', lam, cubed).real
             - 0.5*invariant_2*coefficients)

    deviation = (np.max(np.abs(tower - invariant_2*coefficients/3.0))
                 / np.max(np.abs(tower)))
    assert deviation > 0.1


def test_latent_roots_match_eigenvalues():
    r"""The quartic roots are the eigenvalues of the traceless part."""
    rng = np.random.default_rng(1)

    worst = 0.0
    for _ in range(400):
        matrix = random_hermitian(rng)
        roots = oscprob4nu.psi_roots_4nu(*oscprob4nu.su4_invariants(matrix))
        reference = np.sort(np.linalg.eigvalsh(traceless(matrix)))
        worst = max(worst, np.max(np.abs(np.asarray(roots) - reference)))

    assert worst < 1.e-11


def test_evolution_operator_matches_expm():
    r"""U_4 is the matrix exponential of the traceless Hamiltonian."""
    rng = np.random.default_rng(2)

    worst = 0.0
    for _ in range(400):
        matrix = random_hermitian(rng)
        baseline = rng.uniform(0.1, 50.0)
        operator = np.asarray(
            oscprob4nu.evolution_operator_4nu(matrix, baseline))
        reference = scipy.linalg.expm(-1.j*traceless(matrix)*baseline)
        worst = max(worst, np.max(np.abs(operator - reference)))

    assert worst < 1.e-12


def test_probabilities_are_unitary():
    r"""Each initial flavor's probabilities sum to one."""
    rng = np.random.default_rng(3)

    worst = 0.0
    for _ in range(400):
        matrix = random_hermitian(rng)
        probabilities = np.asarray(oscprob4nu.probabilities_4nu(
            matrix, rng.uniform(0.1, 50.0))).reshape(4, 4)
        worst = max(worst, np.max(np.abs(probabilities.sum(axis=1) - 1.0)))

    assert worst < 1.e-12


def test_su4_coefficients_reconstruct_the_operator():
    r"""``U_4 = u_0 + i u_a lambda^a`` reproduces the operator.

    This is the test of the closed-form ``u_a``, which are the
    four-flavor analogue of Eqs. (10)-(11) of the paper and the part of
    this module with no counterpart to copy from.
    """
    rng = np.random.default_rng(4)

    worst = 0.0
    for _ in range(200):
        matrix = random_hermitian(rng)
        baseline = rng.uniform(0.1, 20.0)
        coefficients = oscprob4nu.evolution_operator_4nu_u_coefficients(
            matrix, baseline)
        rebuilt = (coefficients[0]*np.eye(4)
                   + 1.j*np.einsum('a,aij->ij',
                                   np.asarray(coefficients[1:]),
                                   oscprob4nu.LAMBDA_SU4))
        operator = np.asarray(
            oscprob4nu.evolution_operator_4nu(matrix, baseline))
        worst = max(worst, np.max(np.abs(rebuilt - operator)))

    assert worst < 1.e-13


def test_probability_ordering_matches_the_other_flavors():
    r"""``P[4*alpha+beta]`` is ``|U[beta][alpha]|^2``."""
    matrix = np.diag([1.0, 2.0, 3.0, -6.0]).astype(complex)
    matrix[0, 1] = 0.3
    matrix[1, 0] = 0.3

    operator = np.asarray(oscprob4nu.evolution_operator_4nu(matrix, 1.0))
    probabilities = np.asarray(oscprob4nu.probabilities_4nu(matrix, 1.0))

    for alpha in range(4):
        for beta in range(4):
            assert np.isclose(probabilities[4*alpha + beta],
                              abs(operator[beta][alpha])**2)


@pytest.mark.parametrize('spectrum', [
    [2.5, 2.5, 2.5, 2.5],          # proportional to the identity
    [0.0, 0.0, 0.0, 0.0],          # the zero Hamiltonian
    [2.0, 2.0, -2.0, -2.0],        # two degenerate pairs
    [3.0, -1.0, -1.0, -1.0],       # a triply degenerate root
    [1.0, 1.0, 2.0, -4.0],         # one degenerate pair
    [1.0, 1.0 + 1.e-13, 2.0, -4.0],   # nearly, but not exactly, degenerate
])
@pytest.mark.parametrize('baseline', [1.0, 37.0])
def test_degenerate_spectra_are_exact(spectrum, baseline):
    r"""Repeated latent roots are handled without a special case.

    Interpolating in Newton form rather than solving a Vandermonde
    system is what makes this work: the Vandermonde matrix is singular
    the moment two roots coincide, and half of these spectra made it so.
    """
    matrix = np.diag(spectrum).astype(complex)

    operator = np.asarray(
        oscprob4nu.evolution_operator_4nu(matrix, baseline))
    reference = scipy.linalg.expm(-1.j*traceless(matrix)*baseline)

    assert np.all(np.isfinite(operator))
    assert np.max(np.abs(operator - reference)) < 1.e-11


def test_batched_agrees_with_scalar(backend):
    r"""A stack gives what the scalar path gives, to round-off.

    Run on both backends, because with Numba installed the batched call
    is the compiled kernel while the one-by-one calls stay on the NumPy
    path --- so this is the test that compares the two implementations
    element by element, and without the fixture it would only ever
    exercise whichever backend happened to be present.

    The bar was ``1e-15`` while both sides ran the same code, which is
    tighter than two implementations of the same expansion can be
    expected to agree; a few ulp on probabilities of order one is
    ``4e-15``.
    """
    rng = np.random.default_rng(5)
    stack = np.stack([random_hermitian(rng) for _ in range(40)])
    baselines = rng.uniform(1.0, 10.0, 40)

    batched = oscprob4nu.probabilities_4nu(stack, baselines)
    one_by_one = np.array([oscprob4nu.probabilities_4nu(matrix, baseline)
                           for matrix, baseline in zip(stack, baselines)])

    assert batched.shape == (40, 16)
    assert np.max(np.abs(batched - one_by_one)) < 1.e-13


def test_broadcasting_gives_an_oscillogram(backend):
    r"""Hamiltonians and baselines broadcast against each other."""
    rng = np.random.default_rng(6)
    stack = np.stack([random_hermitian(rng) for _ in range(12)])
    baselines = rng.uniform(1.0, 10.0, 5)

    probabilities = oscprob4nu.probabilities_4nu(stack[:, None],
                                                 baselines[None, :])

    assert probabilities.shape == (12, 5, 16)


def test_a_decoupled_sterile_state_reproduces_three_flavors():
    r"""With the sterile angles off, the active block is the 3nu one.

    The strongest check available: an independent module, written years
    earlier against a different algebra, computing the same numbers.
    """
    h_vacuum_4nu = vacuum_4nu()
    h_vacuum_3nu = np.asarray(
        hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
            gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, gd.DCP_NO_BF,
            gd.D21_NO_BF, gd.D31_NO_BF))

    four = np.asarray(oscprob4nu.probabilities_4nu(
        h_vacuum_4nu/ENERGY, BASELINE)).reshape(4, 4)
    three = np.asarray(oscprob3nu.probabilities_3nu(
        h_vacuum_3nu/ENERGY, BASELINE)).reshape(3, 3)

    assert np.max(np.abs(four[:3, :3] - three)) < 1.e-11
    assert np.max(np.abs(four[:3, 3])) < 1.e-14
    assert np.max(np.abs(four[3, :3])) < 1.e-14


def test_a_decoupled_sterile_state_reproduces_three_flavors_in_matter():
    r"""The same, with a matter potential and its sterile entry."""
    h_vacuum_4nu = vacuum_4nu()
    h_vacuum_3nu = np.asarray(
        hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
            gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, gd.DCP_NO_BF,
            gd.D21_NO_BF, gd.D31_NO_BF))

    four = np.asarray(oscprob4nu.probabilities_4nu(
        hamiltonians4nu.hamiltonian_4nu_matter(
            h_vacuum_4nu, ENERGY, gd.VCC_EARTH_CRUST, gd.VNC_EARTH_CRUST),
        BASELINE)).reshape(4, 4)
    three = np.asarray(oscprob3nu.probabilities_3nu(
        hamiltonians3nu.hamiltonian_3nu_matter(
            h_vacuum_3nu, ENERGY, gd.VCC_EARTH_CRUST),
        BASELINE)).reshape(3, 3)

    assert np.max(np.abs(four[:3, :3] - three)) < 1.e-11


def test_mixing_matrix_is_unitary_and_reduces_to_the_pmns_matrix():
    r"""The 3+1 matrix is unitary, and contains the PMNS one."""
    mixing = hamiltonians4nu.mixing_matrix_4nu(
        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF,
        np.sqrt(0.10), np.sqrt(0.10), np.sqrt(0.05), gd.DCP_NO_BF)
    assert np.max(np.abs(mixing.conj().T @ mixing - np.eye(4))) < 1.e-14

    decoupled = hamiltonians4nu.mixing_matrix_4nu(
        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, 0.0, 0.0, 0.0,
        gd.DCP_NO_BF)
    pmns = np.asarray(hamiltonians3nu.pmns_mixing_matrix(
        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, gd.DCP_NO_BF))
    assert np.max(np.abs(decoupled[:3, :3] - pmns)) < 1.e-14


def test_sterile_matter_entry_has_the_right_sign_and_size():
    r"""The sterile entry is ``-V_NC``, positive and half of ``V_CC``.

    Getting this wrong is the four-flavor analogue of the antineutrino
    sign trap: invisible in vacuum, and it moves the resonance.
    """
    h_vacuum = vacuum_4nu(s14=np.sqrt(0.10), s24=np.sqrt(0.10))
    h_matter = hamiltonians4nu.hamiltonian_4nu_matter(
        h_vacuum, ENERGY, gd.VCC_EARTH_CRUST, gd.VNC_EARTH_CRUST)

    added = h_matter - h_vacuum/ENERGY

    assert np.isclose(added[0, 0].real, gd.VCC_EARTH_CRUST)
    assert np.isclose(added[1, 1].real, 0.0)
    assert np.isclose(added[2, 2].real, 0.0)
    assert np.isclose(added[3, 3].real, -gd.VNC_EARTH_CRUST)
    assert gd.VNC_EARTH_CRUST < 0.0
    assert np.isclose(-gd.VNC_EARTH_CRUST/gd.VCC_EARTH_CRUST, 0.5)


def test_nsi_leaves_the_sterile_row_and_column_alone():
    r"""Sterile states have no standard interactions, so no NSI."""
    h_vacuum = vacuum_4nu(s14=np.sqrt(0.10))
    h_nsi = hamiltonians4nu.hamiltonian_4nu_nsi(
        h_vacuum, ENERGY, gd.VCC_EARTH_CRUST, gd.VNC_EARTH_CRUST, gd.EPS_3)
    h_matter = hamiltonians4nu.hamiltonian_4nu_matter(
        h_vacuum, ENERGY, gd.VCC_EARTH_CRUST, gd.VNC_EARTH_CRUST)

    difference = h_nsi - h_matter

    assert np.max(np.abs(difference[3, :])) < 1.e-30
    assert np.max(np.abs(difference[:, 3])) < 1.e-30
    assert np.max(np.abs(difference[:3, :3])) > 0.0


def test_root_polishing_is_what_makes_a_stiff_spectrum_accurate():
    r"""The refinement against the matrix earns its cost.

    A 3+1 spectrum with an eV-scale sterile spans four orders of
    magnitude, and the closed-form roots then carry only what the three
    invariants retain in double precision.  This pins the improvement so
    that removing the refinement cannot pass unnoticed.
    """
    h_matter = hamiltonians4nu.hamiltonian_4nu_matter(
        vacuum_4nu(s14=np.sqrt(0.10), s24=np.sqrt(0.10)),
        ENERGY, gd.VCC_EARTH_CRUST, gd.VNC_EARTH_CRUST)

    reference = np.abs(scipy.linalg.expm(
        -1.j*traceless(h_matter)*BASELINE).T)**2

    original = oscprob4nu.POLISH_ROOTS
    try:
        oscprob4nu.POLISH_ROOTS = True
        polished = np.max(np.abs(np.asarray(oscprob4nu.probabilities_4nu(
            h_matter, BASELINE)).reshape(4, 4) - reference))
        oscprob4nu.POLISH_ROOTS = False
        unpolished = np.max(np.abs(np.asarray(oscprob4nu.probabilities_4nu(
            h_matter, BASELINE)).reshape(4, 4) - reference))
    finally:
        oscprob4nu.POLISH_ROOTS = original

    assert polished < 1.e-9
    assert unpolished > 10.0*polished


def nearly_degenerate_pair(rng, separation):
    r"""Returns a Hermitian 4x4 with two eigenvalues that nearly coincide.

    Built by conjugating a chosen spectrum with a random unitary, so the
    degeneracy is not aligned with the basis and cannot be spotted by
    looking at the matrix.
    """
    base = rng.normal()
    spectrum = np.array([base, base + separation,
                         base + rng.uniform(0.5, 2.0), 0.0])
    spectrum[3] = -spectrum[:3].sum()

    matrix = rng.standard_normal((4, 4)) + 1.j*rng.standard_normal((4, 4))
    unitary, _ = np.linalg.qr(matrix)
    rotated = unitary @ np.diag(spectrum).astype(complex) @ unitary.conj().T

    return 0.5*(rotated + rotated.conj().T)


def test_a_nearly_degenerate_pair_does_not_derail_the_refinement():
    r"""The Newton step is refused where it would cross a neighbour.

    The step divides by a product of gaps, so a root that nearly
    coincides with another divides by almost nothing.  Before the step
    was guarded, a pair separated by one unit in the last place was
    thrown across the spectrum --- 0.8793 became 0.0180 --- and the
    probabilities that followed were not probabilities at all, reaching
    twenty-one and summing to sixty-nine.

    The guard is that a step may not carry a root more than halfway to
    its nearest neighbour.  Whether a pair lands on identical bits or on
    adjacent ones is decided by the last bit of a square root taken near
    zero, so testing the derivative against exactly zero, which is what
    this did before, was a knife edge: it gave different answers for a
    stack and a scalar call, and for the NumPy path and the kernel.

    Sweeping the separation across the range where it bites, rather than
    pinning one case, because the failure was found at 2e-8 and the
    knife edge is at 1e-16.

    What is asserted is what the guard actually promises, which is less
    than accuracy.  A nearly degenerate pair is not something the closed
    form resolves well in the first place: Euler's reduction recovers
    the pair's separation as the square root of a resolvent root that
    vanishes as the pair closes, and a square root near zero turns the
    :math:`\epsilon` on that root into :math:`\sqrt{\epsilon}` on the
    separation, some :math:`10^{-8}` relative.  The refinement cannot
    undo that and is not meant to.  What it must not do is make it
    worse, which unguarded it did by seven orders of magnitude.

    So: the refined roots are never worse than the unrefined ones, and
    stay bounded.  Case by case, since that is where the content is ---
    on the sweep below the guard refuses the step outright for about
    forty per cent of the spectra, and which of them happens to carry
    the worst error is a matter of the last bit and differs between
    machines.  The aggregate comparison is therefore ``<=`` and not
    ``<``; asserting the strict form passed here and failed on CI under
    two of the five Python versions, for exactly that reason.
    """
    rng = np.random.default_rng(24680)

    worst_refined = 0.0
    worst_unrefined = 0.0
    original = oscprob4nu.POLISH_ROOTS
    try:
        for _ in range(400):
            separation = 10.0**rng.uniform(-16.0, -6.0)
            matrix = nearly_degenerate_pair(rng, separation)

            traceless = oscprob4nu._traceless_part(matrix)
            reference = np.sort(np.linalg.eigvalsh(traceless))
            scale = np.max(np.abs(reference))

            oscprob4nu.POLISH_ROOTS = True
            refined = oscprob4nu._latent_roots(traceless[None])[0]
            oscprob4nu.POLISH_ROOTS = False
            unrefined = oscprob4nu._latent_roots(traceless[None])[0]

            refined_error = np.max(np.abs(refined - reference))/scale
            unrefined_error = np.max(np.abs(unrefined - reference))/scale

            assert refined_error <= 1.5*unrefined_error + 1.e-15, (
                'refining a pair separated by %.2e made the roots worse, '
                '%.2e against %.2e relative; the step is meant to be '
                'refused there, not amplified'
                % (separation, refined_error, unrefined_error))

            worst_refined = max(worst_refined, refined_error)
            worst_unrefined = max(worst_unrefined, unrefined_error)
    finally:
        oscprob4nu.POLISH_ROOTS = original

    # Unguarded, the worst of these reached 4.8 --- roots thrown clean
    # across a spectrum of order one
    assert worst_refined < 1.e-4, (
        'a nearly degenerate pair moved the latent roots by %.2e relative'
        % worst_refined)
    assert worst_refined <= worst_unrefined


def test_the_step_guard_leaves_ordinary_spectra_untouched():
    r"""The guard costs nothing where it is not needed.

    A guard that quietly declined useful refinements would show up as a
    loss of accuracy on generic spectra rather than as a failure, so it
    is pinned from the other side too: on ordinary random Hamiltonians
    the refined roots must still reach round-off, which they cannot do
    from the closed form alone.

    That the refinement *helps* is asserted where it is unambiguous, in
    `test_root_polishing_is_what_makes_a_stiff_spectrum_accurate`.  Here
    the aggregate is only required not to get worse, because on a
    generic spectrum the closed form is already near round-off and which
    element carries the worst error is decided by the last bit.
    """
    rng = np.random.default_rng(13579)

    worst_refined = 0.0
    worst_unrefined = 0.0
    original = oscprob4nu.POLISH_ROOTS
    try:
        for _ in range(200):
            matrix = random_hermitian(rng)
            traceless = oscprob4nu._traceless_part(matrix)
            reference = np.sort(np.linalg.eigvalsh(traceless))
            scale = np.max(np.abs(reference))

            oscprob4nu.POLISH_ROOTS = True
            refined = oscprob4nu._latent_roots(traceless[None])[0]
            oscprob4nu.POLISH_ROOTS = False
            unrefined = oscprob4nu._latent_roots(traceless[None])[0]

            worst_refined = max(worst_refined,
                                np.max(np.abs(refined - reference))/scale)
            worst_unrefined = max(worst_unrefined,
                                  np.max(np.abs(unrefined - reference))/scale)
    finally:
        oscprob4nu.POLISH_ROOTS = original

    assert worst_refined < 1.e-14
    assert worst_refined <= worst_unrefined


def test_roots_of_a_traceless_free_hamiltonian_all_vanish():
    r"""A Hamiltonian proportional to the identity has no evolution.

    All three invariants vanish, so the resolvent cubic has a triple
    root and its trigonometric form degenerates; the confluent branch
    has to return zeros rather than divide by zero.
    """
    roots = oscprob4nu.psi_roots_4nu(0.0, 0.0, 0.0)

    assert len(roots) == 4
    assert np.allclose(roots, 0.0, atol=1.e-300)


def test_evolution_operator_returns_an_array_for_a_stack():
    r"""A stacked call returns an array, a scalar one a nested list."""
    rng = np.random.default_rng(8)
    stack = np.stack([random_hermitian(rng) for _ in range(6)])

    batched = oscprob4nu.evolution_operator_4nu(stack, 3.0)
    single = oscprob4nu.evolution_operator_4nu(stack[0], 3.0)

    assert isinstance(batched, np.ndarray)
    assert batched.shape == (6, 4, 4)
    assert isinstance(single, list)
    assert np.max(np.abs(np.asarray(single) - batched[0])) < 1.e-15


def test_probabilities_accept_a_nested_list():
    r"""A plain nested list works, as at two and three flavors."""
    nested = [[1.0, 0.0, 0.0, 0.0],
              [0.0, 2.0, 0.0, 0.0],
              [0.0, 0.0, 3.0, 0.0],
              [0.0, 0.0, 0.0, -6.0]]

    probabilities = oscprob4nu.probabilities_4nu(nested, 1.0)

    assert len(probabilities) == 16
    assert np.isclose(sum(probabilities[0:4]), 1.0)
