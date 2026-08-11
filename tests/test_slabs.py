# -*- coding: utf-8 -*-
r"""Tests of propagation across a sequence of adjacent slabs.

The claim these have to check is not that the slab machinery runs, but
that composing slabs is the same operation as propagating once.  Two
independent handles on that: splitting a single constant Hamiltonian
into many identical slabs must reproduce the one-slab answer exactly,
and a sequence of genuinely different slabs must reproduce the product
of matrix exponentials computed by ``scipy.linalg.expm``.
"""

import numpy as np
import pytest
from scipy.linalg import expm

import oscprob2nu
import oscprob3nu
import oscprob4nu
import slabs

from conftest import ATOL, traceless


def random_hermitian(rng, n):
    r"""Returns a random ``n``-by-``n`` complex Hermitian matrix."""
    a = rng.normal(size=(n, n)) + 1.0j*rng.normal(size=(n, n))
    return (a + a.conj().T)/2.0


@pytest.mark.parametrize('n_flavors, evolution, slab_evolution',
                         [(2, oscprob2nu.evolution_operator_2nu,
                           slabs.evolution_operator_2nu_slabs),
                          (3, oscprob3nu.evolution_operator_3nu,
                           slabs.evolution_operator_3nu_slabs)])
@pytest.mark.parametrize('n_slabs', [1, 2, 7, 40])
def test_splitting_one_slab_reproduces_it(n_flavors, evolution,
                                          slab_evolution, n_slabs):
    r"""Cutting a constant Hamiltonian into pieces changes nothing.

    This is the property that makes slab propagation exact rather than
    approximate: with the same Hamiltonian throughout, the number of
    cuts is arbitrary and the answer cannot depend on it.
    """
    rng = np.random.default_rng(20260801 + n_slabs)
    h = random_hermitian(rng, n_flavors)
    total = 1.7

    stack = np.repeat(h[None, ...], n_slabs, axis=0)
    widths = np.full(n_slabs, total/n_slabs)

    assert np.allclose(np.asarray(slab_evolution(stack, widths)),
                       np.asarray(evolution(h, total)), atol=ATOL)


@pytest.mark.parametrize('n_flavors, slab_evolution',
                         [(2, slabs.evolution_operator_2nu_slabs),
                          (3, slabs.evolution_operator_3nu_slabs)])
def test_agrees_with_a_product_of_matrix_exponentials(n_flavors,
                                                      slab_evolution):
    r"""Different slabs compose as ``expm`` says they do.

    The independent computation, and the one that pins the *order* of
    the product: reversing it changes the answer whenever the slab
    Hamiltonians do not commute, which random ones do not.

    Each factor uses the *traceless* part of its Hamiltonian, as the
    single-slab tests do, because the expansions return
    :math:`e^{-i H_0 L}` and drop the phase carried by the trace.
    """
    rng = np.random.default_rng(4242)
    n_slabs = 6
    stack = np.array([random_hermitian(rng, n_flavors)
                      for _ in range(n_slabs)])
    widths = rng.uniform(0.2, 1.5, size=n_slabs)

    expected = np.eye(n_flavors, dtype=complex)
    for h, w in zip(stack, widths):
        expected = expm(-1.0j*traceless(h)*w) @ expected

    assert np.allclose(np.asarray(slab_evolution(stack, widths)), expected,
                       atol=ATOL)


@pytest.mark.parametrize('n_flavors, slab_evolution',
                         [(2, slabs.evolution_operator_2nu_slabs),
                          (3, slabs.evolution_operator_3nu_slabs)])
def test_the_order_of_the_slabs_matters(n_flavors, slab_evolution):
    r"""Reversing non-commuting slabs changes the result.

    Guards against an implementation that happens to pass the tests
    above by composing in the wrong order and being saved by symmetry.
    """
    rng = np.random.default_rng(99)
    stack = np.array([random_hermitian(rng, n_flavors) for _ in range(3)])
    widths = np.array([0.4, 0.9, 1.3])

    forward = np.asarray(slab_evolution(stack, widths))
    backward = np.asarray(slab_evolution(stack[::-1], widths[::-1]))

    assert not np.allclose(forward, backward, atol=1.e-6)


@pytest.mark.parametrize('n_flavors, slab_evolution',
                         [(2, slabs.evolution_operator_2nu_slabs),
                          (3, slabs.evolution_operator_3nu_slabs)])
def test_the_evolution_operator_is_unitary(n_flavors, slab_evolution):
    r"""A product of unitary operators is unitary."""
    rng = np.random.default_rng(7)
    stack = np.array([random_hermitian(rng, n_flavors) for _ in range(5)])
    widths = rng.uniform(0.1, 2.0, size=5)

    u = np.asarray(slab_evolution(stack, widths))

    assert np.allclose(u.conj().T @ u, np.eye(n_flavors), atol=ATOL)


@pytest.mark.parametrize('probabilities, n_flavors, n_prob',
                         [(slabs.probabilities_2nu_slabs, 2, 4),
                          (slabs.probabilities_3nu_slabs, 3, 9)])
def test_probabilities_are_normalized(probabilities, n_flavors, n_prob):
    r"""Each initial flavor's probabilities sum to one."""
    rng = np.random.default_rng(11)
    stack = np.array([random_hermitian(rng, n_flavors) for _ in range(4)])
    widths = rng.uniform(0.1, 1.0, size=4)

    prob = np.array(probabilities(stack, widths))

    assert len(prob) == n_prob
    assert np.all(prob >= 0.0)
    # The initial flavor varies slowest, so each block of n_flavors is one
    # initial flavor's worth of transition probabilities.
    for start in range(0, n_prob, n_flavors):
        assert prob[start:start+n_flavors].sum() == pytest.approx(1.0,
                                                                  abs=ATOL)


@pytest.mark.parametrize('probabilities, evolution, n_flavors',
                         [(slabs.probabilities_2nu_slabs,
                           oscprob2nu.probabilities_2nu, 2),
                          (slabs.probabilities_3nu_slabs,
                           oscprob3nu.probabilities_3nu, 3)])
def test_probabilities_match_the_single_slab_routine(probabilities,
                                                     evolution, n_flavors):
    r"""One slab is the ordinary, non-slab computation."""
    rng = np.random.default_rng(3)
    h = random_hermitian(rng, n_flavors)

    assert np.allclose(np.array(probabilities(h[None, ...], [2.3])),
                       np.array(evolution(h, 2.3)), atol=ATOL)


@pytest.mark.parametrize('slab_evolution, n_flavors',
                         [(slabs.evolution_operator_2nu_slabs, 2),
                          (slabs.evolution_operator_3nu_slabs, 3)])
def test_a_zero_width_slab_is_the_identity(slab_evolution, n_flavors):
    r"""A slab of no width contributes nothing, rather than failing."""
    rng = np.random.default_rng(5)
    h = random_hermitian(rng, n_flavors)
    stack = np.array([h, random_hermitian(rng, n_flavors)])

    with_zero = np.asarray(slab_evolution(stack, [1.1, 0.0]))
    without = np.asarray(slab_evolution(h[None, ...], [1.1]))

    assert np.allclose(with_zero, without, atol=ATOL)


@pytest.mark.parametrize('slab_evolution, n_flavors',
                         [(slabs.evolution_operator_2nu_slabs, 2),
                          (slabs.evolution_operator_3nu_slabs, 3)])
def test_malformed_slab_sequences_raise(slab_evolution, n_flavors):
    r"""Every guard in `slabs._check_slabs` fires on its own case."""
    good = np.zeros((3, n_flavors, n_flavors), dtype=complex)

    with pytest.raises(ValueError, match='one width per slab'):
        slab_evolution(good, [1.0, 2.0])

    with pytest.raises(ValueError, match='at least one slab'):
        slab_evolution(np.zeros((0, n_flavors, n_flavors)), [])

    with pytest.raises(ValueError, match='cannot be negative'):
        slab_evolution(good, [1.0, -1.0, 1.0])

    with pytest.raises(ValueError, match='must have shape'):
        slab_evolution(np.zeros((3, 4, 4)), [1.0, 1.0, 1.0])

    with pytest.raises(ValueError, match='must have shape'):
        slab_evolution(np.zeros((n_flavors, n_flavors)), [1.0])

    with pytest.raises(ValueError, match='one-dimensional'):
        slab_evolution(good, [[1.0], [1.0], [1.0]])


# --------------------------------------------------------------------------
# Four flavors, which slabs could not compose until 1.11.0
# --------------------------------------------------------------------------

def test_4nu_slabs_reproduce_a_product_of_matrix_exponentials(rng):
    r"""Composing four-flavor slabs is propagating once, per slab.

    The same claim the two- and three-flavor tests make, against the
    same independent reference: ``scipy.linalg.expm`` of each traceless
    slab Hamiltonian, multiplied in the order the neutrino meets them.
    """
    for n_slabs in (1, 2, 5, 17):
        h_stack = np.stack([random_hermitian(rng, 4) for _ in range(n_slabs)])
        widths = rng.uniform(0.1, 3.0, n_slabs)

        ours = np.asarray(slabs.evolution_operator_4nu_slabs(h_stack, widths))

        reference = np.eye(4, dtype=complex)
        for h_matrix, width in zip(h_stack, widths):
            reference = expm(-1.j*traceless(h_matrix)*width) @ reference

        assert np.allclose(ours, reference, atol=1.0e-12)


def test_4nu_one_hamiltonian_split_into_slabs_is_the_same_propagation(rng):
    r"""Cutting one constant Hamiltonian into slabs changes nothing."""
    h_matrix = random_hermitian(rng, 4)
    total = 4.0

    once = np.asarray(oscprob4nu.evolution_operator_4nu(h_matrix, total))
    for n_slabs in (2, 7, 40):
        split = np.asarray(slabs.evolution_operator_4nu_slabs(
            np.stack([h_matrix]*n_slabs), np.full(n_slabs, total/n_slabs)))
        assert np.allclose(split, once, atol=1.0e-12)


def test_4nu_slab_probabilities_are_unitary(rng):
    r"""Each initial flavor's sixteen probabilities sum to one."""
    h_stack = np.stack([random_hermitian(rng, 4) for _ in range(6)])
    prob = np.asarray(slabs.probabilities_4nu_slabs(
        h_stack, rng.uniform(0.1, 2.0, 6))).reshape(4, 4)

    assert np.allclose(prob.sum(axis=1), 1.0, atol=ATOL)
    assert np.allclose(prob.sum(axis=0), 1.0, atol=ATOL)


def test_4nu_slabs_reject_a_malformed_sequence(rng):
    r"""The same validation the other flavor counts get."""
    h_stack = np.stack([random_hermitian(rng, 4) for _ in range(3)])

    with pytest.raises(ValueError, match='one width per slab'):
        slabs.evolution_operator_4nu_slabs(h_stack, [1.0, 2.0])
    with pytest.raises(ValueError, match=r'shape \(n, 4, 4\)'):
        slabs.evolution_operator_4nu_slabs(
            np.stack([random_hermitian(rng, 3) for _ in range(3)]), [1., 2., 3.])
    with pytest.raises(ValueError, match='cannot be negative'):
        slabs.evolution_operator_4nu_slabs(h_stack, [1.0, -1.0, 1.0])


# --- Batches of chords sharing one geometry ---------------------------
#
# The public routines took one chord at a time until 1.13.1.  What has to
# hold now is that a batch is the same answer as the loop it replaces, for
# every flavor count and every leading shape, and that a single chord still
# returns exactly the tuple it always did.


@pytest.mark.parametrize('n_flavors, probabilities, evolution',
                         [(2, slabs.probabilities_2nu_slabs,
                           slabs.evolution_operator_2nu_slabs),
                          (3, slabs.probabilities_3nu_slabs,
                           slabs.evolution_operator_3nu_slabs),
                          (4, slabs.probabilities_4nu_slabs,
                           slabs.evolution_operator_4nu_slabs)])
def test_a_batch_of_chords_is_the_loop_it_replaces(n_flavors, probabilities,
                                                   evolution, rng):
    r"""Batched and per-chord evaluation agree, at every flavor count.

    They agree to round-off rather than bit for bit: the two take
    different paths through the compiled backend, which is why this is
    an ``allclose`` and not an equality.
    """
    n_slabs, n_chords = 5, 6
    widths = rng.uniform(0.1, 2.0, n_slabs)
    stack = np.stack([np.stack([random_hermitian(rng, n_flavors)
                                for _ in range(n_slabs)])
                      for _ in range(n_chords)])

    batched = np.asarray(probabilities(stack, widths))
    looped = np.stack([np.asarray(probabilities(chord, widths))
                       for chord in stack])
    assert batched.shape == (n_chords, n_flavors*n_flavors)
    assert np.allclose(batched, looped, atol=ATOL)

    batched_u = np.asarray(evolution(stack, widths))
    looped_u = np.stack([np.asarray(evolution(chord, widths))
                         for chord in stack])
    assert batched_u.shape == (n_chords, n_flavors, n_flavors)
    assert np.allclose(batched_u, looped_u, atol=ATOL)


@pytest.mark.parametrize('n_flavors, probabilities',
                         [(2, slabs.probabilities_2nu_slabs),
                          (3, slabs.probabilities_3nu_slabs),
                          (4, slabs.probabilities_4nu_slabs)])
def test_one_chord_still_returns_the_tuple_it_always_did(n_flavors,
                                                         probabilities, rng):
    r"""The un-batched call is untouched: a tuple, not an array."""
    n_slabs = 4
    widths = rng.uniform(0.1, 2.0, n_slabs)
    chord = np.stack([random_hermitian(rng, n_flavors)
                      for _ in range(n_slabs)])

    prob = probabilities(chord, widths)
    assert isinstance(prob, tuple)
    assert len(prob) == n_flavors*n_flavors


@pytest.mark.parametrize('leading', [(2,), (2, 3), (1, 1, 4)])
def test_any_number_of_leading_axes_is_carried_through(leading, rng):
    r"""The batch axes are whatever the caller brought, not just one."""
    n_slabs = 3
    widths = rng.uniform(0.1, 2.0, n_slabs)
    stack = np.stack([random_hermitian(rng, 3)
                      for _ in range(int(np.prod(leading))*n_slabs)])
    stack = stack.reshape(leading + (n_slabs, 3, 3))

    prob = np.asarray(slabs.probabilities_3nu_slabs(stack, widths))
    assert prob.shape == leading + (9,)

    flat = np.asarray(slabs.probabilities_3nu_slabs(
        stack.reshape(-1, n_slabs, 3, 3), widths))
    assert np.allclose(prob.reshape(-1, 9), flat, atol=ATOL)


def test_batched_slab_probabilities_are_unitary(rng):
    r"""Unitarity holds chord by chord across a batch."""
    n_slabs, n_chords = 4, 5
    widths = rng.uniform(0.1, 2.0, n_slabs)
    stack = np.stack([np.stack([random_hermitian(rng, 3)
                                for _ in range(n_slabs)])
                      for _ in range(n_chords)])

    prob = np.asarray(slabs.probabilities_3nu_slabs(
        stack, widths)).reshape(n_chords, 3, 3)
    assert np.allclose(prob.sum(axis=2), 1.0, atol=ATOL)
    assert np.allclose(prob.sum(axis=1), 1.0, atol=ATOL)


def test_a_batch_is_validated_like_a_single_chord(rng):
    r"""The batch path rejects the same malformed input, and says so."""
    stack = np.stack([np.stack([random_hermitian(rng, 3) for _ in range(3)])
                      for _ in range(2)])

    with pytest.raises(ValueError, match='one width per slab'):
        slabs.probabilities_3nu_slabs(stack, [1.0, 2.0])
    with pytest.raises(ValueError, match='cannot be negative'):
        slabs.probabilities_3nu_slabs(stack, [1.0, -1.0, 1.0])
    with pytest.raises(ValueError, match=r'shape \(\.\.\., n, 3, 3\)'):
        slabs.probabilities_3nu_slabs(
            np.zeros((2, 3, 4, 4), dtype=complex), [1.0, 2.0, 3.0])
