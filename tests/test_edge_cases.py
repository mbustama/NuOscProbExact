# -*- coding: utf-8 -*-
r"""Tests of numerically awkward inputs.

Degenerate Hamiltonians make the SU(3) invariant :math:`|h|^2` vanish,
which the latent-root formula divides by; near-degenerate ones can push
the argument of the arc cosine marginally outside :math:`[-1, 1]`
through round-off.  Neither may produce NaN or a non-unitary operator.
"""

import importlib
import pathlib

import numpy as np
import pytest

import hamiltonians2nu
import hamiltonians3nu
import oscprob2nu
import oscprob3nu

from conftest import ATOL, as_nested_list, random_hermitian


DEGENERATE_3NU = [
    ('proportional to the identity', np.eye(3, dtype=complex)),
    ('zero', np.zeros((3, 3), dtype=complex)),
    ('large multiple of the identity', 1.0e6*np.eye(3, dtype=complex)),
    ('two degenerate eigenvalues', np.diag([1.0, 1.0, -2.0]).astype(complex)),
    ('nearly degenerate', np.diag([1.0, 1.0 + 1.0e-14, -2.0]).astype(complex)),
]


@pytest.mark.parametrize('label,h_matrix', DEGENERATE_3NU,
                         ids=[c[0] for c in DEGENERATE_3NU])
@pytest.mark.parametrize('l', [0.0, 1.0, 1.0e3])
def test_3nu_evolution_operator_is_finite_for_degenerate_hamiltonians(
        label, h_matrix, l):
    r"""U_3 stays finite and unitary for degenerate Hamiltonians."""
    u = np.array(oscprob3nu.evolution_operator_3nu(
        as_nested_list(h_matrix), l))
    assert np.all(np.isfinite(u)), 'non-finite U for a %s Hamiltonian' % label
    assert np.allclose(u.conj().T @ u, np.eye(3), atol=1.0e-10)


@pytest.mark.parametrize('label,h_matrix', DEGENERATE_3NU,
                         ids=[c[0] for c in DEGENERATE_3NU])
def test_3nu_probabilities_are_finite_for_degenerate_hamiltonians(
        label, h_matrix):
    r"""The 3nu probabilities stay finite and normalized."""
    prob = np.array(oscprob3nu.probabilities_3nu(
        as_nested_list(h_matrix), 1.0)).reshape(3, 3)
    assert np.all(np.isfinite(prob))
    assert np.allclose(prob.sum(axis=1), np.ones(3), atol=1.0e-10)


def test_identity_hamiltonian_gives_no_oscillations():
    r"""A Hamiltonian proportional to the identity has no traceless
    part, so no flavor transitions occur."""
    for scale in [0.0, 1.0, 1.0e6]:
        prob = np.array(oscprob3nu.probabilities_3nu(
            as_nested_list(scale*np.eye(3, dtype=complex)), 10.0)
        ).reshape(3, 3)
        assert np.allclose(prob, np.eye(3), atol=1.0e-10)


@pytest.mark.parametrize('h_matrix', [
    np.eye(2, dtype=complex),
    np.zeros((2, 2), dtype=complex),
    1.0e6*np.eye(2, dtype=complex),
])
def test_2nu_evolution_operator_is_finite_for_degenerate_hamiltonians(
        h_matrix):
    r"""U_2 stays finite and unitary when |h| vanishes."""
    u = np.array(oscprob2nu.evolution_operator_2nu(
        as_nested_list(h_matrix), 1.0))
    assert np.all(np.isfinite(u))
    assert np.allclose(u.conj().T @ u, np.eye(2), atol=1.0e-10)


def test_2nu_probabilities_are_finite_for_degenerate_hamiltonians():
    r"""The 2nu probabilities stay finite when |h| vanishes."""
    prob = np.array(oscprob2nu.probabilities_2nu(
        as_nested_list(np.eye(2, dtype=complex)), 1.0))
    assert np.all(np.isfinite(prob))
    assert np.allclose(prob.reshape(2, 2), np.eye(2), atol=1.0e-10)


def test_psi_roots_are_real_for_hermitian_hamiltonians(rng):
    r"""The latent roots of a Hermitian Hamiltonian are real.

    Round-off can push the arc-cosine argument marginally outside
    [-1, 1]; the implementation must clip it rather than return complex
    roots, which would spoil unitarity.
    """
    from conftest import random_hermitian
    for _ in range(200):
        h_matrix = random_hermitian(rng, 3)
        h_coeffs = oscprob3nu.hamiltonian_3nu_coefficients(
            as_nested_list(h_matrix))
        h2, h3 = oscprob3nu.su3_invariants(h_coeffs)
        psi = np.array(oscprob3nu.psi_roots(h2, h3))
        assert np.allclose(np.imag(psi), 0.0, atol=1.0e-9)


@pytest.mark.parametrize('h2,h3', [
    (1.0, 2.0/3.0/np.sqrt(3.0)*3.0),   # arc-cosine argument at -1
    (1.0, -2.0/3.0/np.sqrt(3.0)*3.0),  # arc-cosine argument at +1
    (1.0, 10.0),                       # argument far outside [-1, 1]
    (1.0, -10.0),
])
def test_psi_roots_clip_the_arc_cosine_argument(h2, h3):
    r"""psi_roots returns real, finite roots even when its argument is
    pushed outside the domain of the arc cosine."""
    psi = np.array(oscprob3nu.psi_roots(h2, h3))
    assert np.all(np.isfinite(psi))
    assert np.allclose(np.imag(psi), 0.0, atol=ATOL)


def test_large_baseline_keeps_unitarity(rng):
    r"""Unitarity survives many oscillation cycles."""
    from conftest import random_hermitian
    for _ in range(20):
        h_matrix = random_hermitian(rng, 3)
        u = np.array(oscprob3nu.evolution_operator_3nu(
            as_nested_list(h_matrix), 1.0e6))
        assert np.allclose(u.conj().T @ u, np.eye(3), atol=1.0e-8)


# --------------------------------------------------------------------------
# Input the expansions cannot make sense of, which used to pass silently
# --------------------------------------------------------------------------

NOT_HERMITIAN = {
    2: [[1.0+0.j, 2.0+0.j], [0.0+0.j, 1.0+0.j]],
    3: [[1.0+0.j, 2.0+0.j, 0.j], [0.j, 1.0+0.j, 0.j], [0.j, 0.j, 1.0+0.j]],
    4: [[1.0+0.j, 2.0+0.j, 0.j, 0.j], [0.j, 1.0+0.j, 0.j, 0.j],
        [0.j, 0.j, 1.0+0.j, 0.j], [0.j, 0.j, 0.j, 1.0+0.j]],
}


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
@pytest.mark.parametrize('routine', ['probabilities', 'evolution_operator'])
def test_a_non_hermitian_hamiltonian_is_refused(n_flavors, routine):
    r"""A non-Hermitian Hamiltonian raises rather than returning numbers.

    This is the failure that had to be caught by checking, because
    nothing downstream reveals it: the expansion returns probabilities
    that still **sum to one**, so the sanity check a caller would
    actually apply cannot tell that the answer is meaningless.  Before
    the check existed, the matrix below returned 1.000000 for the first
    row and no warning of any kind.
    """
    module = importlib.import_module('oscprob%dnu' % n_flavors)
    call = getattr(module, '%s_%dnu' % (routine, n_flavors))

    with pytest.raises(ValueError, match='not Hermitian'):
        call(NOT_HERMITIAN[n_flavors], 1.0)


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
def test_a_non_hermitian_stack_is_refused(n_flavors):
    r"""The batched path checks too, not only the scalar one."""
    module = importlib.import_module('oscprob%dnu' % n_flavors)
    call = getattr(module, 'probabilities_%dnu' % n_flavors)
    stack = np.stack([np.array(NOT_HERMITIAN[n_flavors], dtype=complex)]*50)

    with pytest.raises(ValueError, match='not Hermitian'):
        call(stack, np.full(50, 1.0))


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
def test_the_hermiticity_check_can_be_switched_off(n_flavors):
    r"""The switch is the escape hatch for scans that cannot afford it.

    Checking costs a pass over the stack, which is the same order as the
    evaluation; `CHECK_HERMITICITY` exists so that a caller whose
    Hamiltonians are Hermitian by construction can decline it.
    """
    module = importlib.import_module('oscprob%dnu' % n_flavors)
    call = getattr(module, 'probabilities_%dnu' % n_flavors)

    original = module.CHECK_HERMITICITY
    try:
        module.CHECK_HERMITICITY = False
        prob = call(NOT_HERMITIAN[n_flavors], 1.0)
    finally:
        module.CHECK_HERMITICITY = original

    assert len(prob) == n_flavors*n_flavors


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
def test_hamiltonians_the_library_builds_pass_the_check(n_flavors, rng):
    r"""The tolerance is loose enough for matrices assembled in floating
    point, which is the whole difficulty of checking this at all."""
    module = importlib.import_module('oscprob%dnu' % n_flavors)
    call = getattr(module, 'probabilities_%dnu' % n_flavors)

    for _ in range(20):
        h_matrix = random_hermitian(rng, n_flavors)
        prob = np.asarray(call(as_nested_list(h_matrix), 1.5))
        assert np.isclose(prob[:n_flavors].sum(), 1.0, atol=ATOL)


SINE_OUT_OF_RANGE = [
    ('hamiltonians2nu', 'mixing_matrix_2nu', (1.5,)),
    ('hamiltonians3nu', 'pmns_mixing_matrix', (1.5, 0.1, 0.1, 0.0)),
    ('hamiltonians4nu', 'mixing_matrix_4nu',
     (1.5, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0)),
]


@pytest.mark.parametrize('module_name,routine,args', SINE_OUT_OF_RANGE)
def test_a_sine_outside_its_range_is_refused_the_same_way(module_name,
                                                          routine, args):
    r"""All three flavor counts refuse it, and say the same thing.

    They did not.  Two and three flavors took the cosine with
    :mod:`math` and raised ``math domain error``, which names neither
    the parameter nor the value; four flavors took it with
    :func:`numpy.sqrt`, which returns ``nan`` and let it run silently
    into the probabilities.  Same invalid input, three behaviours.
    """
    module = importlib.import_module(module_name)

    with pytest.raises(ValueError, match='sine of an angle'):
        getattr(module, routine)(*args)


def test_a_sine_of_exactly_one_is_still_allowed():
    r"""The boundary is inclusive: maximal mixing is a legal angle."""
    assert np.isclose(hamiltonians2nu.mixing_matrix_2nu(1.0)[0][0], 0.0)
    assert np.isclose(hamiltonians3nu.pmns_mixing_matrix(1.0, 0.1, 0.1,
                                                         0.0)[0][0].real, 0.0)


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
def test_a_complex_diagonal_entry_is_refused(n_flavors):
    r"""Hermiticity constrains the diagonal too, not only the pairs.

    A Hamiltonian whose off-diagonal entries are properly conjugate can
    still fail to be Hermitian, by carrying an imaginary part on the
    diagonal.  That is a separate branch of the check, and a separate
    way to get a wrong answer that still sums to one.
    """
    module = importlib.import_module('oscprob%dnu' % n_flavors)
    call = getattr(module, 'probabilities_%dnu' % n_flavors)

    h_matrix = np.eye(n_flavors, dtype=complex)
    h_matrix[0, 0] = 1.0 + 0.5j

    with pytest.raises(ValueError, match='diagonal entry'):
        call([[complex(z) for z in row] for row in h_matrix], 1.0)

    # And as a stack.  The two paths stopped sharing an implementation when
    # `_check_hermitian` gained a branch for a single matrix, and this
    # assertion is what noticed: the scalar case above had been the only
    # thing covering the batched diagonal refusal, so adding that branch
    # left it exercised by nothing while the suite stayed green.
    with pytest.raises(ValueError, match='diagonal entry'):
        call(np.stack([h_matrix]*4), 1.0)


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
def test_an_empty_stack_passes_the_check_and_returns_empty(n_flavors):
    r"""Nothing to check, and nothing to compute.

    An empty stack has no entries to compare, so the check returns
    early; the scale it would otherwise take is undefined on an empty
    array and would raise from inside NumPy.
    """
    module = importlib.import_module('oscprob%dnu' % n_flavors)
    call = getattr(module, 'probabilities_%dnu' % n_flavors)

    empty = np.zeros((0, n_flavors, n_flavors), dtype=complex)
    result = call(empty, np.zeros(0))

    assert result.shape == (0, n_flavors*n_flavors)


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
def test_every_entry_point_honours_the_switch_in_both_positions(n_flavors,
                                                                 rng):
    r"""Both settings, on every public route into the expansions.

    The check sits at seven entry points, and a switch that is only ever
    exercised in one position is half tested: the branch that skips the
    check is the one a production scan actually takes.  This runs each
    route with the switch on and off, over the scalar path, the batched
    NumPy path and the compiled kernel, and requires the answers to be
    identical --- checking must not change what is computed.
    """
    module = importlib.import_module('oscprob%dnu' % n_flavors)
    probabilities = getattr(module, 'probabilities_%dnu' % n_flavors)
    operator = getattr(module, 'evolution_operator_%dnu' % n_flavors)

    scalar = random_hermitian(rng, n_flavors)
    # Large enough to reach the compiled kernel where one is used
    stack = np.stack([random_hermitian(rng, n_flavors) for _ in range(64)])
    baselines = rng.uniform(0.1, 20.0, 64)

    original = module.CHECK_HERMITICITY
    results = {}
    try:
        for setting in (True, False):
            module.CHECK_HERMITICITY = setting
            results[setting] = (
                np.asarray(probabilities(as_nested_list(scalar), 1.5)),
                np.asarray(operator(as_nested_list(scalar), 1.5)),
                np.asarray(probabilities(stack, baselines)),
                np.asarray(operator(stack, baselines)),
            )
    finally:
        module.CHECK_HERMITICITY = original

    for checked, unchecked in zip(results[True], results[False]):
        assert np.array_equal(checked, unchecked)


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
@pytest.mark.parametrize('bad', [np.inf, -np.inf, np.nan])
def test_a_non_finite_entry_cannot_disable_the_check(n_flavors, bad):
    r"""A non-finite entry must not switch the check off.

    The tolerance is relative to the largest entry, so one infinity made
    the scale infinite, the tolerance infinite, and every comparison
    false --- a Hamiltonian that was both non-finite *and* non-Hermitian
    then passed a check whose entire purpose is to refuse the second.
    The matrix below is non-Hermitian in the (0, 1) pair as well as
    non-finite, so it must be refused on one ground or the other.
    """
    module = importlib.import_module('oscprob%dnu' % n_flavors)
    call = getattr(module, 'probabilities_%dnu' % n_flavors)

    h_matrix = np.eye(n_flavors, dtype=complex)
    h_matrix[0, 1] = 2.0          # not the conjugate of h[1, 0], which is 0
    h_matrix[0, 0] = bad

    with pytest.raises(ValueError):
        call([[complex(z) for z in row] for row in h_matrix], 1.0)


@pytest.mark.parametrize('n_flavors', [2, 3, 4])
@pytest.mark.parametrize('bad', [np.inf, -np.inf, np.nan])
def test_a_non_finite_entry_is_refused_on_both_paths(n_flavors, bad):
    r"""Non-finite alone is enough, and the two paths must agree.

    The test above makes its matrix non-finite *and* non-Hermitian, so it
    is refused either way and cannot tell which ground did it.  This one
    is Hermitian apart from being non-finite, so only the non-finite
    branch can refuse it --- and it checks the scalar and the batched
    path separately, because they do not share an implementation.

    Both of those gaps were real.  The scalar path added in 1.11.0 built
    its scale with `max`, which keeps its running value when compared
    against a nan, so a nan never reached `isfinite` and a Hamiltonian
    the batched path refuses came back with probabilities instead.  The
    whole suite passed.
    """
    module = importlib.import_module('oscprob%dnu' % n_flavors)
    call = getattr(module, 'probabilities_%dnu' % n_flavors)

    h_matrix = np.eye(n_flavors, dtype=complex)
    h_matrix[0, 0] = bad

    with pytest.raises(ValueError, match='non-finite'):
        call([[complex(z) for z in row] for row in h_matrix], 1.0)

    with pytest.raises(ValueError, match='non-finite'):
        call(np.stack([h_matrix]*4), 1.0)


HELPER_COPIES = [
    ('_check_hermitian',
     ['oscprob2nu', 'oscprob3nu', 'oscprob4nu']),
    ('_cos_from_sin',
     ['hamiltonians2nu', 'hamiltonians3nu', 'hamiltonians4nu']),
]


@pytest.mark.parametrize('name,modules', HELPER_COPIES)
def test_the_duplicated_helpers_have_not_drifted(name, modules):
    r"""Six copies of two helpers, kept identical by a test.

    :mod:`oscprob2nu` and :mod:`oscprob3nu` are documented as
    self-contained --- copying either into another project is a
    supported way to use this library --- so a shared module for these
    would break the property that makes that work.  Duplication is the
    deliberate cost, and this is what stops it becoming drift: the
    bodies must match character for character, apart from the flavor
    count and the module named in the error message.

    That matters most for `_check_hermitian`, where a fix applied to one
    copy and not the others would leave two flavor counts silently
    accepting what the third refuses --- which is the class of bug the
    check was added to remove.
    """
    import inspect
    import re

    def normalise(source):
        source = re.sub(r'range\((?:i\+1, )?\d\)',
                        lambda m: m.group(0).replace(m.group(0)[-2], 'N'),
                        source)
        source = re.sub(r'oscprob\dnu', 'oscprobNnu', source)
        source = re.sub(r'hamiltonians\dnu', 'hamiltoniansNnu', source)
        source = re.sub(r'\(\.\.\., \d, \d\)', '(..., N, N)', source)
        source = re.sub(r'the \w+ independent pairs', 'the N independent pairs',
                        source)
        return source

    bodies = {}
    for module_name in modules:
        module = importlib.import_module(module_name)
        bodies[module_name] = normalise(
            inspect.getsource(getattr(module, name)))

    reference = bodies[modules[0]]
    for module_name in modules[1:]:
        assert bodies[module_name] == reference, (
            '%s.%s has drifted from %s.%s; the copies are duplicated on '
            'purpose and must stay identical'
            % (module_name, name, modules[0], name))


def test_the_degeneracy_tolerance_is_not_free_to_widen():
    r"""``DEGENERACY_TOL`` separates two exact expressions, not two
    approximations.

    Below it the two-projector form is used, which is exact for a
    repeated root and an approximation otherwise; above it the general
    Lagrange form is used, which is exact for distinct roots and
    singular at a repeated one.  Widening the tolerance from 1e-12 to
    1e-2 left the whole suite green while taking a spectrum whose roots
    differ by a part in a thousand from 2e-14 against ``expm`` to
    7e-3 --- eleven orders of magnitude, unnoticed.
    """
    from scipy.linalg import expm

    h_matrix = np.diag([1.0, 1.0 + 3.0e-3, -2.0 - 3.0e-3]).astype(complex)
    traceless = h_matrix - np.trace(h_matrix).real/3.0*np.eye(3)
    reference = expm(-1.j*traceless*5.0)

    operator = np.asarray(oscprob3nu.evolution_operator_3nu(
        as_nested_list(h_matrix), 5.0))
    assert np.max(np.abs(operator - reference)) < 1.0e-11

    # The separation above sits four orders above the tolerance, so the
    # general form must be the one taken; pin the tolerance itself, since
    # that is what decides it
    assert oscprob3nu.DEGENERACY_TOL <= 1.0e-10


# --------------------------------------------------------------------------
# The documented "copy one file into your project" route
# --------------------------------------------------------------------------

@pytest.mark.parametrize('module_name', ['oscprob2nu', 'oscprob3nu',
                                          'oscprob4nu'])
def test_a_core_module_works_copied_out_on_its_own(module_name, tmp_path):
    r"""One file, copied out, with nothing else from this repository.

    ``README.md`` and ``installation.rst`` both call this a supported way
    to use the library.  It stopped being one in 1.6.0, when the optional
    compiled backend arrived and was imported unconditionally: a lone
    copy raised ``ImportError: No module named 'fastkernels'``.  It went
    unnoticed for six releases because any check run from inside the
    repository finds ``src/`` on the path and imports the real module.

    So this runs in a subprocess with the repository stripped from
    ``sys.path``, which is the only way to reproduce a user's situation
    from within the project.
    """
    import shutil
    import subprocess
    import sys
    import textwrap

    source = pathlib.Path(__file__).resolve().parents[1]/'src'/(
        module_name + '.py')
    shutil.copy(source, tmp_path/(module_name + '.py'))

    n_flavors = int(module_name[len('oscprob')])
    script = textwrap.dedent('''
        import sys
        sys.path = [p for p in sys.path
                    if 'NuOscProbExact' not in p and p not in ('', '.')]
        sys.path.insert(0, %r)

        import numpy as np
        import %s as module

        n = %d
        call = getattr(module, 'probabilities_%%dnu' %% n)

        # The backend is absent, so `worthwhile` must decline it rather
        # than raise, and the NumPy path must carry both a scalar call and
        # a stack large enough to have reached the kernel.
        assert module.fastkernels.worthwhile(n, 10**9) is False

        h = np.diag(np.arange(1.0, n+1.0)).astype(complex)
        h[0, 1] = 0.3 + 0.2j
        h[1, 0] = 0.3 - 0.2j

        scalar = call([[complex(z) for z in row] for row in h], 1.5)
        assert len(scalar) == n*n
        assert abs(sum(scalar[0:n]) - 1.0) < 1.0e-12

        stack = call(np.stack([h]*500), np.full(500, 1.5))
        assert stack.shape == (500, n*n)
        assert abs(stack[0].reshape(n, n).sum(axis=1)[0] - 1.0) < 1.0e-12

        print('ok')
    ''') % (str(tmp_path), module_name, n_flavors)

    result = subprocess.run([sys.executable, '-c', script],
                            capture_output=True, text=True)

    assert result.returncode == 0, (
        '%s does not work copied out on its own:\n%s'
        % (module_name, result.stderr[-1500:]))
    assert result.stdout.strip().endswith('ok')
