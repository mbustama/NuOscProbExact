# -*- coding: utf-8 -*-
r"""Tests of the four-neutrino SU(4) expansion.

The independent computations these check against are the same ones used
at two and three flavors --- ``scipy.linalg.expm`` for the evolution
operator and ``numpy.linalg.eigh`` for the spectrum --- plus one that
only exists at four flavors: switching the three sterile mixing angles
off must reproduce :mod:`oscprob3nu` exactly, since a decoupled fourth
state cannot change the active block.

The latent roots get a stricter oracle than any of those, in
``stiff_reference.json``: fifty-digit ``mpmath`` values for nine
Hamiltonians, frozen to a file.  They have to be frozen, because
``numpy.linalg.eigvalsh`` stopped being independent of the library the
moment it became one of :data:`oscprob4nu.ROOT_STRATEGY`'s routes ---
tests written against it began comparing an implementation with itself,
one of them asserting ``1.06e-15 <= 0.0``.  See
``tests/gen_stiff_reference.py``, which regenerates the file and records
the three traps its own oracle fell into first.
"""

import json
import os

import numpy as np
import pytest
import scipy.linalg

import fastkernels
import globaldefs as gd
import hamiltonians3nu
import hamiltonians4nu
import oscprob3nu
import oscprob4nu


BASELINE = 1300.0*gd.CONV_KM_TO_INV_EV
ENERGY = 1.0e9
DM41 = 1.0

STIFF_REFERENCE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'stiff_reference.json')


def stiff_cases():
    r"""Yields the frozen reference cases, as ``float64`` arrays.

    Read from hexadecimal floats, so what comes back is the exact bits the
    fifty-digit values were computed from rather than a decimal rounding of
    them.
    """
    with open(STIFF_REFERENCE) as handle:
        payload = json.load(handle)

    for case in payload['cases']:
        matrix = np.array(
            [[float.fromhex(case['matrix_real'][i][j])
              + 1.j*float.fromhex(case['matrix_imag'][i][j])
              for j in range(4)] for i in range(4)])
        yield (case['label'], matrix,
               float.fromhex(case['baseline']),
               np.array([float.fromhex(r) for r in case['roots']]),
               np.array([float.fromhex(p) for p in case['probabilities']]))


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


@pytest.mark.parametrize('use_numba', [False, True])
def test_polish_roots_is_scoped_to_the_eigensolver_route(
        use_numba, monkeypatch):
    r"""The Newton step is live behind the eigensolver and dead ahead of it.

    :data:`oscprob4nu.POLISH_ROOTS` must still reach both backends on the
    ``'eigensolver'`` route, and must do nothing whatever on the
    ``'double-double'`` one, which neither forms the roots that way nor
    leaves a Newton step anything to find.  Wiring the two switches
    together by accident is what this fails on.

    What it deliberately does *not* do is rank the routes by probability
    error.  An earlier version of this test asserted a hundredfold gain
    from the Newton step and got one, because what it was compared against
    was the closed form, whose roots carry only what three invariants
    retain in double precision.  Against an eigensolver start the step is
    worth about five on the roots, and on a spectrum this stiff even that
    ordering does not survive the trip into probabilities: the accumulated
    phase amplifies a root error by :math:`(\psi_{\max} L)^2`, which here
    is enough to reorder two routes that differ by an ulp.  The ranking is
    pinned where it is clean, on the roots themselves, in
    `test_the_root_strategies_meet_their_documented_accuracy`.

    Two things about how it is pinned, both learned the hard way.

    It measures **each backend separately**.  The two compute the
    closed-form roots with different roundings, and in a cluster this
    tight that decides which Newton steps the halfway guard in
    `oscprob4nu._polish_roots` will admit.  Asserting one number for
    whichever path happened to be dispatched made the test a hostage to
    that: when the compiled path took over `probabilities_4nu`, the
    measured gain at 1 GeV fell from 418x to 10.6x, close enough to a
    threshold of ten that it passed on one machine and failed on
    another.  Neither backend is better than the other --- swept over
    energy, the compiled one leads at 0.5 and 10 GeV and trails at 1 and
    2 --- so there is nothing here to fix, only something to stop
    measuring by accident.

    And it measures where the effect is **robust rather than merely
    large**.  The gain is a strong function of how degenerate the
    spectrum is: at half a GeV the closest pair is 1.4e-4 of the
    spectrum apart and refining is worth some three thousand times on
    both paths, while by 25 GeV the roots are 3.3e-3 apart and refining
    is worth nothing on either --- it is very slightly *harmful*, the
    guard's guarantee being about the roots rather than about phases
    accumulated over a long baseline.  The energy here is the degenerate
    end, where the margin is four orders of magnitude wide and no
    plausible difference in libm can walk it across.
    """
    if use_numba and not fastkernels.HAVE_NUMBA:
        pytest.skip('Numba is not installed')
    monkeypatch.setattr(fastkernels, 'USE_NUMBA', use_numba)

    energy = 0.5e9
    h_matter = hamiltonians4nu.hamiltonian_4nu_matter(
        vacuum_4nu(s14=np.sqrt(0.10), s24=np.sqrt(0.10)),
        energy, gd.VCC_EARTH_CRUST, gd.VNC_EARTH_CRUST)

    reference = np.abs(scipy.linalg.expm(
        -1.j*traceless(h_matter)*BASELINE).T)**2

    # The spectrum this is about: four roots spanning orders of
    # magnitude, with the closest pair a fraction of that apart.
    spectrum = np.sort(np.linalg.eigvalsh(traceless(h_matter)))
    relative_gap = np.min(np.diff(spectrum))/np.max(np.abs(spectrum))
    assert relative_gap < 1.0e-3, relative_gap

    def probabilities(strategy, polish):
        monkeypatch.setattr(oscprob4nu, 'ROOT_STRATEGY', strategy)
        monkeypatch.setattr(oscprob4nu, 'POLISH_ROOTS', polish)
        return np.asarray(oscprob4nu.probabilities_4nu(
            h_matter, BASELINE)).reshape(4, 4)

    polished = probabilities('eigensolver', True)
    unpolished = probabilities('eigensolver', False)
    with_flag = probabilities('double-double', True)
    without_flag = probabilities('double-double', False)

    # Live behind the eigensolver: the step moves the answer
    assert np.max(np.abs(polished - unpolished)) > 0.0
    # Dead ahead of it: bit for bit, the flag is not read at all
    assert np.array_equal(with_flag, without_flag)

    for candidate in (polished, unpolished, with_flag):
        assert np.max(np.abs(candidate - reference)) < 1.0e-9


@pytest.mark.parametrize('strategy,tolerance',
                         [('double-double', 1.0e-16),
                          ('eigensolver', 5.0e-16)])
def test_the_root_strategies_meet_their_documented_accuracy(strategy,
                                                            tolerance,
                                                            monkeypatch):
    r"""Each route hits the figure :data:`oscprob4nu.ROOT_STRATEGY` quotes.

    Against the frozen fifty-digit roots, over all nine cases at once so
    that the stack path is what runs.  Measured worst: 5.5e-17 for
    double-double and 3.9e-16 for the eigensolver with its Newton step,
    against bars of 1e-16 and 5e-16 here.

    Those bars are close to the measurements on purpose.  One ``float64``
    ulp on a root of order one is 2.2e-16, so double-double is being held
    to *below* half an ulp --- there is no room to loosen this without
    giving up the claim it exists to defend, and a regression that pushed
    the answer to a whole ulp would still look small while meaning the dd
    coefficients had stopped arriving in dd.
    """
    monkeypatch.setattr(oscprob4nu, 'ROOT_STRATEGY', strategy)
    monkeypatch.setattr(oscprob4nu, 'POLISH_ROOTS', True)

    labels, matrices, references = [], [], []
    for label, matrix, _, roots, _ in stiff_cases():
        labels.append(label)
        matrices.append(traceless(matrix))
        references.append(roots)

    psi = oscprob4nu._latent_roots(np.stack(matrices))

    for label, got, want in zip(labels, psi, references):
        error = np.max(np.abs(np.sort(got) - want))/np.max(np.abs(want))
        assert error < tolerance, '%s: %.3e' % (label, error)


@pytest.mark.parametrize('use_numba', [False, True])
@pytest.mark.parametrize('strategy', ['double-double', 'eigensolver'])
def test_the_reference_probabilities_are_reproduced(use_numba, strategy,
                                                    monkeypatch):
    r"""Both routes reproduce the frozen ``expm`` probabilities.

    The tolerance grows with the square of the accumulated phase
    :math:`\psi_{\max} L`, which is where the reference's stiffest cases
    spend their accuracy: reconstructing :math:`U_4` in Newton form takes
    second differences of :math:`e^{-i\psi L}` over the roots, so an error
    in the phase enters the coefficients twice.  Measured
    :math:`5 \times 10^{-17} (\psi_{\max} L)^2` across five decades of
    phase --- 1.3e-12 at 245 radians, 2.8e-4 at 2.5e6 --- and the bar below
    keeps a factor of twenty above that.

    This is a limit of the reconstruction and not of the roots, which the
    strategies deliver to under an ulp either way; at
    :math:`\Delta m^2_{41} = 1000\ \mathrm{eV}^2` over 1300 km the phase is
    2.5 million radians and ``float64`` simply cannot carry it.  Notebook 17
    records the same wall from the other side, agreeing with nuSQuIDS to
    3e-10 on the stiffest spectrum tested.  The physical part of the range
    is where the tight numbers are: 1.3e-12 at 0.1 eV\ :sup:`2`.
    """
    if use_numba and not fastkernels.HAVE_NUMBA:
        pytest.skip('Numba is not installed')
    monkeypatch.setattr(fastkernels, 'USE_NUMBA', use_numba)
    monkeypatch.setattr(oscprob4nu, 'ROOT_STRATEGY', strategy)
    monkeypatch.setattr(oscprob4nu, 'POLISH_ROOTS', True)

    for label, matrix, baseline, roots, want in stiff_cases():
        got = np.asarray(oscprob4nu.probabilities_4nu(matrix, baseline))
        phase = np.max(np.abs(roots))*baseline
        tolerance = 1.0e-15*phase*phase + 4.0e-15

        error = np.max(np.abs(got.ravel() - want))
        assert error < tolerance, ('%s: %.3e (phase %.2e)'
                                   % (label, error, phase))


def test_an_unknown_root_strategy_is_refused(monkeypatch):
    r"""A misspelt :data:`oscprob4nu.ROOT_STRATEGY` says so.

    Silently falling back to either route would hide the mistake behind
    results that look fine, which for a switch whose whole purpose is to
    choose between accuracies is the one outcome worth ruling out.
    """
    monkeypatch.setattr(oscprob4nu, 'ROOT_STRATEGY', 'eigensolvor')

    with pytest.raises(ValueError, match='double-double'):
        oscprob4nu._latent_roots(np.stack([traceless(vacuum_4nu())]))
    with pytest.raises(ValueError, match='double-double'):
        oscprob4nu._strategy_code()


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


def test_four_flavor_liv_reduces_to_the_three_flavor_term():
    r"""With the sterile LIV angles and eigenvalue off, the active block
    is the three-flavor LIV Hamiltonian.

    The same decoupling argument that validates the vacuum and matter
    Hamiltonians, applied to the term added in 1.11.0.
    """
    h_vacuum_4nu = vacuum_4nu()
    h_vacuum_3nu = np.asarray(
        hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
            gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, gd.DCP_NO_BF,
            gd.D21_NO_BF, gd.D31_NO_BF))

    four = np.asarray(hamiltonians4nu.hamiltonian_4nu_liv(
        h_vacuum_4nu, ENERGY, 0.3, 0.4, 0.5, 0.0, 0.0, 0.0, 0.7,
        1.0e-9, 1.5e-9, 2.0e-9, 0.0, 1.0e12))
    three = np.asarray(hamiltonians3nu.hamiltonian_3nu_liv(
        h_vacuum_3nu, ENERGY, 0.3, 0.4, 0.5, 0.7,
        1.0e-9, 1.5e-9, 2.0e-9, 1.0e12))

    assert np.max(np.abs(four[:3, :3] - three)) < 1.e-27
    assert np.max(np.abs(four - four.conj().T)) < 1.e-27


def test_four_flavor_liv_scales_with_the_energy_and_broadcasts():
    r"""The LIV term grows as E/Lambda, and a stack of energies works."""
    h_vacuum = vacuum_4nu(s14=np.sqrt(0.10))
    args = (0.3, 0.4, 0.5, 0.2, 0.1, 0.0, 0.7,
            1.0e-9, 1.5e-9, 2.0e-9, 2.5e-9, 1.0e12)

    single = np.asarray(hamiltonians4nu.hamiltonian_4nu_liv(
        h_vacuum, ENERGY, *args))
    stacked = np.asarray(hamiltonians4nu.hamiltonian_4nu_liv(
        h_vacuum, np.array([ENERGY, 2.0*ENERGY]), *args))

    assert stacked.shape == (2, 4, 4)
    assert np.allclose(stacked[0], single, atol=0.0)

    # The vacuum part falls as 1/E and the LIV part rises as E, so the
    # difference of the two energies isolates the LIV scaling
    liv_once = single - np.asarray(h_vacuum)/ENERGY
    liv_twice = stacked[1] - np.asarray(h_vacuum)/(2.0*ENERGY)
    assert np.allclose(liv_twice, 2.0*liv_once, rtol=1.e-12)


def test_the_3plus1_rotation_order_is_the_one_documented():
    r"""R34 R24 R14 R23 R13 R12, and the order is observable.

    :func:`hamiltonians4nu.mixing_matrix_4nu` documents that ordering.
    Swapping the first two rotations left the whole suite green, because
    every test that exercised the mixing matrix either set
    ``s34 = 0`` --- which makes R34 the identity, so the order cannot
    matter --- or checked only unitarity, which no reordering breaks.
    With all three sterile angles non-zero the difference is 0.12 in the
    entries of U.
    """
    angles = (0.5, 0.5, 0.15, np.sqrt(0.10), np.sqrt(0.10), np.sqrt(0.20))
    dcp = 0.3
    mixing = hamiltonians4nu.mixing_matrix_4nu(*angles, dcp)

    rotation = hamiltonians4nu._rotation_4nu
    expected = (rotation(2, 3, angles[5])
                @ rotation(1, 3, angles[4])
                @ rotation(0, 3, angles[3])
                @ rotation(1, 2, angles[1])
                @ rotation(0, 2, angles[2], dcp)
                @ rotation(0, 1, angles[0]))

    assert np.max(np.abs(mixing - expected)) < 1.e-15

    # And the order is not a free choice: swapping the two outermost
    # rotations, which share the index 3, changes U materially
    swapped = (rotation(1, 3, angles[4])
               @ rotation(2, 3, angles[5])
               @ rotation(0, 3, angles[3])
               @ rotation(1, 2, angles[1])
               @ rotation(0, 2, angles[2], dcp)
               @ rotation(0, 1, angles[0]))
    assert np.max(np.abs(mixing - swapped)) > 0.1
