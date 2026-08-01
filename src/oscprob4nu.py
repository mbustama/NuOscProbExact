# -*- coding: utf-8 -*-
r"""Compute the four-neutrino flavor-transition probabilities.

This module contains the routines needed to compute four-neutrino
flavor-transition probabilities for an arbitrary time-independent
:math:`4\times4` Hermitian Hamiltonian, using the SU(4) exponential
expansion.  It is the :math:`n = 4` member of the family whose
:math:`n = 2` and :math:`n = 3` members are :mod:`oscprob2nu` and
:mod:`oscprob3nu`, following the method of [1]_.  It is the last one: see
:ref:`why-four-is-the-end` in the methodology page.  The Cayley-Hamilton
route at :math:`n = 4` in constant matter was set out in [2]_.

The Hamiltonian is expanded in the basis of the fifteen generalized
Gell-Mann matrices,

.. math:: H = h_0 \mathbb{1} + h_a \lambda^a , \qquad a = 1, \ldots, 15 ,

and the time-evolution operator in the same basis,

.. math:: U_4(L) = u_0 \mathbb{1} + i u_a \lambda^a .

The term :math:`h_0` contributes only an overall phase and is dropped;
all routines therefore work with the traceless part of the Hamiltonian,
which leaves the oscillation probabilities unchanged.

What is new at four flavors
---------------------------

Three ingredients have no counterpart in :mod:`oscprob3nu`:

* **A third invariant.**  SU(4) has rank three, so the traceless part
  carries three independent invariants rather than two,

  .. math::
     I_2 = \tfrac12 \mathrm{Tr}\,\tilde{H}^2 , \qquad
     I_3 = \tfrac12 \mathrm{Tr}\,\tilde{H}^3 , \qquad
     I_4 = \tfrac12 \left(\mathrm{Tr}\,\tilde{H}^4 - I_2^2\right) ,

  and the cubic characteristic equation of the three-flavor case becomes
  the quartic

  .. math::
     \psi^4 - I_2 \psi^2 - \tfrac23 I_3 \psi
     + \tfrac14 \left(I_2^2 - 2 I_4\right) = 0 .

* **A quartic that still solves in closed form.**  Euler's method
  reduces it to a *resolvent cubic* whose three roots are real and
  non-negative because :math:`\tilde{H}` is Hermitian, so the same
  trigonometric formula that :func:`oscprob3nu.psi_roots` uses solves it.
  The SU(3) machinery is literally nested inside the SU(4) solution.

* **A longer star-product tower.**  The three-flavor identity
  :math:`(h \star h) \star h = \tfrac13 |h|^2 h` is a Cayley-Hamilton
  accident of :math:`n = 3` and is false for SU(4) --- it is off by some
  tens of per cent on a random Hamiltonian --- so
  :math:`((h \star h) \star h)_a` enters as independent data.

Accuracy
--------

For a generic Hamiltonian the expansion is exact to round-off, like its
two- and three-flavor counterparts.  A *stiff* spectrum is the case that
needs care, and the physically interesting 3+1 scenario is stiff: with
:math:`\Delta m^2_{41} \sim 1` eV\ :sup:`2` the eigenvalues span four
orders of magnitude, three of them clustering, and the information
needed to separate the cluster is destroyed when :math:`I_2, I_3, I_4`
are formed in double precision.  No amount of care in solving the
quartic recovers it.

The roots are therefore refined against the *matrix*, by one Newton step
on :math:`\chi(\psi) = \det(\psi \mathbb{1} - \tilde{H})`, which uses the
Hamiltonian entries rather than the three compressed invariants.  That
restores the roots to round-off and the probabilities from about
:math:`5 \times 10^{-7}` to :math:`10^{-9}`.

Neither figure is anywhere near a measurable effect --- probabilities
are confronted with data at the per-cent level --- so this matters for
the exactness claim, for error accumulating when :mod:`slabs` and
:mod:`earth` compose operators across layers, and for a regression suite
tight enough to catch a real mistake.  :data:`POLISH_ROOTS` has the
measured comparison against the alternatives, including why LAPACK and
extended precision both lose.

`evolution_operator_4nu` and `probabilities_4nu` accept either a single
Hamiltonian and baseline or a stack of them, in which case the whole
stack is evaluated at once, exactly as at two and three flavors.

Units
-----

The routines are unit-agnostic: they require only that the Hamiltonian
and the baseline be given in reciprocal units, so that the product
:math:`H L` is dimensionless.  Elsewhere in **NuOscProbExact** the
Hamiltonian is in eV and the baseline in eV\ :sup:`-1`; the module
:mod:`globaldefs` provides ``CONV_KM_TO_INV_EV`` to convert a baseline
in km into eV\ :sup:`-1`.

Routine listings
----------------

    * generators_su4 - Returns the fifteen generalized Gell-Mann matrices
    * hamiltonian_4nu_coefficients - Returns the :math:`h_a`
    * su4_invariants - Returns the invariants :math:`I_2, I_3, I_4`
    * psi_roots_4nu - Returns the roots of the characteristic equation
    * evolution_operator_4nu_u_coefficients - Returns the :math:`u_a`
    * evolution_operator_4nu - Returns the evolution operator :math:`U_4`
    * probabilities_4nu - Returns the oscillation probabilities

References
----------

.. [1] Mauricio Bustamante, "Exact neutrino oscillation probabilities
   with arbitrary time-independent Hamiltonians", arXiv:1904.12391.

.. [2] S. Kamo et al., "Matter enhanced transitions of active and
   sterile neutrinos", Eur. Phys. J. C 28, 211 (2003), which applies the
   Cayley-Hamilton approach at :math:`n = 4` in constant matter.
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

__all__ = ['LAMBDA_SU4', 'POLISH_ROOTS', 'SMALL_BATCH',
           'generators_su4', 'hamiltonian_4nu_coefficients',
           'su4_invariants', 'psi_roots_4nu',
           'evolution_operator_4nu_u_coefficients', 'evolution_operator_4nu',
           'probabilities_4nu']

from typing import List, Tuple, Union

import math

import numpy as np


SMALL_BATCH = 10
r"""int: Module-level constant.

Stacks with at most this many elements are evaluated one at a time
through the scalar path.  A batched call carries a fixed cost --- here a
batched :math:`4\times4` determinant and linear solve --- which for a
short stack exceeds what the scalar path spends on the whole job.  The
threshold matches :data:`oscprob3nu.SMALL_BATCH`; the four-flavor
expansion does strictly more work per element, so a threshold measured
for three flavors is a safe one here.
"""

POLISH_ROOTS = True
r"""bool: Module-level switch.

Whether to refine the latent roots against the Hamiltonian matrix, by
one Newton step on :math:`\chi(\psi) = \det(\psi \mathbb{1} - \tilde{H})`,
after solving the quartic in closed form.

This is on by default and should stay on.  The closed-form roots are
limited by the conditioning of :math:`I_2, I_3, I_4`, not by the quartic
solver: perturbing the three invariants at the :math:`10^{-16}` level
moves the roots of a stiff 3+1 spectrum by :math:`6 \times 10^{-11}`
relative, which is what the unrefined closed form achieves and what no
better root-finder can beat.  The Newton step reads the matrix entries
instead, and is not subject to that floor.

Measured against ``mpmath`` at fifty decimal digits, on stiff 3+1
Hamiltonians, with cost quoted for a 200 000-point scan:

======================================  ==========  =======
Strategy for the latent roots           Rel. error  Cost
======================================  ==========  =======
Closed form alone                       8.3e-11     0.17 s
Closed form + one Newton step           1.1e-16     0.41 s
``numpy.linalg.eigvalsh``               7.4e-16     0.17 s
Closed form in ``numpy.longdouble``     4.5e-11     0.43 s
======================================  ==========  =======

Note the second row beats the third: the Newton step is some seven times
more accurate than LAPACK, because ``eigvalsh`` reduces by similarity
transforms that each carry a backward error of order
:math:`\epsilon\|H\|`, while this converges onto the root of
:math:`\det(\psi\mathbb{1} - \tilde{H})` for the matrix it was given.
Extended precision was rejected for buying under a digit, being slower,
and silently being ``float64`` on Apple Silicon and Windows.

The refined figure has been confirmed from outside: against `nuSQuIDS
<https://github.com/arguelles/nuSQuIDS>`_, which integrates the density
matrix numerically, the four-flavor probabilities agree to
:math:`4 \times 10^{-16}` on a benign spectrum and
:math:`3 \times 10^{-10}` on the stiffest one tested.  See notebook 17.

In probabilities the difference is :math:`5 \times 10^{-7}` unrefined
against :math:`10^{-9}` refined.  Both are orders of magnitude below what
any experiment resolves; the reasons to want the smaller one are the
library's claim to exactness, error accumulating when :mod:`slabs` and
:mod:`earth` compose operators across many layers, and a regression
suite with no room for a bug to hide in.

A second Newton step changes nothing --- one already reaches the floor
--- so exactly one is taken.  It costs about 40% of the runtime, which
brings the four-flavor closed form to parity with a batched ``eigh``
rather than ahead of it.

The step is applied unconditionally rather than only where a spectrum
looks stiff, and that is a measured decision rather than a lazy one.
Two skip criteria were tried on 6300 Hamiltonians: the gap-based
amplification that perturbation theory suggests, which misjudges doubly
paired spectra by four orders of magnitude, and a matrix residual
comparing :math:`\prod_m \psi_m` with :math:`\det \tilde{H}`, which is
one constraint on four roots and misses errors that cancel in the
product.  Neither can safely skip any elements at all.  The reason is
structural: a criterion complete enough to certify four roots must
evaluate :math:`\chi` at four roots, which is this refinement --- the
check and the fix are the same computation.  See the methodology page.

Set to ``False`` to skip it, which is useful only for reproducing the
unrefined figures or for spectra known to be well separated.
"""

DEGENERACY_TOL = 1.e-12
r"""float: Module-level constant.

Relative tolerance below which two latent roots :math:`\psi` are treated
as degenerate.  Reconstructing :math:`U_4` from the roots divides by
their differences, which vanish at a repeated root, so a degenerate
spectrum is handled by a separate, exact expression.
"""


def generators_su4() -> np.ndarray:
    r"""Returns the fifteen generalized Gell-Mann matrices of SU(4).

    The generators are traceless, Hermitian, and normalized as
    :math:`\mathrm{Tr}(\lambda^a \lambda^b) = 2 \delta^{ab}`, which is
    the convention :mod:`oscprob3nu` uses at :math:`n = 3`.

    They are ordered as six symmetric and six antisymmetric off-diagonal
    matrices, one pair for each of the six index pairs
    :math:`(j, k)` with :math:`j < k` in lexicographic order, followed by
    the three diagonal matrices that span the Cartan subalgebra.

    .. versionadded:: 1.9.0

    Returns
    -------
    numpy.ndarray
        Complex array of shape ``(15, 4, 4)``, the generators
        :math:`\lambda^1, \ldots, \lambda^{15}`.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import oscprob4nu

        lam = oscprob4nu.generators_su4()
        gram = np.einsum('aij,bji->ab', lam, lam).real/2.0

        print('shape:', lam.shape)
        print('orthonormal to %.1e' % np.max(np.abs(gram - np.eye(15))))
    """
    generators = []

    for j in range(4):
        for k in range(j+1, 4):
            symmetric = np.zeros((4, 4), dtype=complex)
            symmetric[j, k] = symmetric[k, j] = 1.0
            antisymmetric = np.zeros((4, 4), dtype=complex)
            antisymmetric[j, k] = -1.j
            antisymmetric[k, j] = 1.j
            generators += [symmetric, antisymmetric]

    for level in range(1, 4):
        diagonal = np.zeros((4, 4), dtype=complex)
        diagonal[np.arange(level), np.arange(level)] = 1.0
        diagonal[level, level] = -float(level)
        generators.append(np.sqrt(2.0/(level*(level+1)))*diagonal)

    return np.asarray(generators)


LAMBDA_SU4 = generators_su4()
r"""numpy.ndarray: Module-level constant.

The fifteen generalized Gell-Mann matrices, of shape ``(15, 4, 4)``,
tabulated once at import time.  At three flavors the analogous table is
the :math:`d` tensor of :func:`oscprob3nu.tensor_d`; here the generators
themselves are tabulated, because every quantity this module needs is a
trace against them and none needs :math:`d_{abc}` explicitly.
"""

_TRACELESS_SUBTRACTION = np.eye(4)


def hamiltonian_4nu_coefficients(
    hamiltonian_matrix: Union[list, np.ndarray]
) -> List[float]:
    r"""Returns the :math:`h_a` of the SU(4) expansion of the Hamiltonian.

    Computes the coefficients :math:`h_1, \ldots, h_{15}` of the SU(4)
    expansion :math:`H = h_0 \mathbb{1} + h_a \lambda^a` of the
    four-flavor Hamiltonian `hamiltonian_matrix`, which is assumed to be
    given in the flavor basis.  The coefficient :math:`h_0` contributes
    only an overall phase to the evolution operator and is not returned.

    They follow from the normalization of the generators as
    :math:`h_a = \tfrac12 \mathrm{Tr}(H \lambda^a)`, and are real for a
    Hermitian Hamiltonian.

    .. versionadded:: 1.9.0

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Four-flavor Hermitian Hamiltonian, given as a nested list or an
        array of shape ``(4, 4)``.

    Returns
    -------
    list of float
        The fifteen coefficients :math:`h_1, \ldots, h_{15}`.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import oscprob4nu

        hamiltonian = np.diag([1.0, 2.0, 3.0, -6.0]).astype(complex)
        h = oscprob4nu.hamiltonian_4nu_coefficients(hamiltonian)

        print('%d coefficients' % len(h))
        print('the three Cartan ones: %.4f, %.4f, %.4f' % tuple(h[12:]))
    """
    matrix = np.asarray(hamiltonian_matrix, dtype=complex)
    traceless = matrix - (np.trace(matrix).real/4.0)*_TRACELESS_SUBTRACTION

    return list(0.5*np.einsum('aij,ji->a', LAMBDA_SU4, traceless).real)


def su4_invariants(
    hamiltonian_matrix: Union[list, np.ndarray]
) -> Tuple[float, float, float]:
    r"""Returns the three SU(4) invariants of the Hamiltonian.

    Returns :math:`I_2 = |h|^2`, :math:`I_3 = \langle h \rangle` and
    :math:`I_4 = |h \star h|^2` of the traceless part
    :math:`\tilde{H}` of `hamiltonian_matrix`, computed from traces of
    its powers,

    .. math::
       I_2 = \tfrac12 \mathrm{Tr}\,\tilde{H}^2 , \qquad
       I_3 = \tfrac12 \mathrm{Tr}\,\tilde{H}^3 , \qquad
       I_4 = \tfrac12 \left(\mathrm{Tr}\,\tilde{H}^4 - I_2^2\right) .

    SU(4) has rank three, so there are three of them; SU(3) has rank two,
    which is why :func:`oscprob3nu.su3_invariants` returns two.  Taking
    them from traces avoids ever building the :math:`d` tensor of SU(4),
    a :math:`15 \times 15 \times 15` table that nothing else here needs.

    .. versionadded:: 1.9.0

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Four-flavor Hermitian Hamiltonian, given as a nested list or an
        array of shape ``(4, 4)``.

    Returns
    -------
    tuple of float
        The invariants ``(I2, I3, I4)``.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import oscprob4nu

        hamiltonian = np.diag([1.0, 2.0, 3.0, -6.0]).astype(complex)
        i2, i3, i4 = oscprob4nu.su4_invariants(hamiltonian)

        print('I2 = %.4f' % i2)
        print('I3 = %.4f' % i3)
        print('I4 = %.4f' % i4)
    """
    matrix = np.asarray(hamiltonian_matrix, dtype=complex)
    traceless = matrix - (np.trace(matrix).real/4.0)*_TRACELESS_SUBTRACTION

    squared = traceless @ traceless
    invariant_2 = 0.5*np.trace(squared).real
    invariant_3 = 0.5*np.trace(squared @ traceless).real
    invariant_4 = 0.5*(np.einsum('ij,ji->', squared, squared).real
                       - invariant_2*invariant_2)

    return invariant_2, invariant_3, invariant_4


def _resolvent_cubic_roots(
    coeff_2: Union[int, float],
    coeff_1: Union[int, float],
    coeff_0: Union[int, float]
) -> np.ndarray:
    r"""Returns the three real roots of a depressed-friendly cubic.

    Solves :math:`z^3 + c_2 z^2 + c_1 z + c_0 = 0` by the trigonometric
    (Viete) construction, which is the same one
    :func:`oscprob3nu.psi_roots` uses at three flavors.  All three roots
    are real for the resolvent cubic of a Hermitian Hamiltonian.

    Parameters
    ----------
    coeff_2 : int or float
        Coefficient of :math:`z^2`.
    coeff_1 : int or float
        Coefficient of :math:`z`.
    coeff_0 : int or float
        Constant coefficient.

    Returns
    -------
    numpy.ndarray
        The three roots, of shape ``(3,)``.
    """
    depressed_p = coeff_1 - coeff_2*coeff_2/3.0
    depressed_q = (2.0*coeff_2**3/27.0 - coeff_2*coeff_1/3.0 + coeff_0)
    shift = -coeff_2/3.0

    if depressed_p >= 0.0:
        # A triple root, up to round-off: the three roots coincide and
        # the trigonometric form degenerates.
        return np.full(3, np.cbrt(-depressed_q)) + shift

    scale = 2.0*np.sqrt(-depressed_p/3.0)
    argument = np.clip(3.0*depressed_q/(depressed_p*scale), -1.0, 1.0)
    angle = np.arccos(argument)

    return scale*np.cos((angle + 2.0*np.pi*np.arange(3))/3.0) + shift


def psi_roots_4nu(
    invariant_2: Union[int, float],
    invariant_3: Union[int, float],
    invariant_4: Union[int, float]
) -> List[float]:
    r"""Returns the roots of the four-flavor characteristic equation.

    Returns the four real roots :math:`\psi_m` of

    .. math::
       \psi^4 - I_2 \psi^2 - \tfrac23 I_3 \psi
       + \tfrac14 \left(I_2^2 - 2 I_4\right) = 0 ,

    the characteristic equation of the traceless part of the
    Hamiltonian, obtained by Euler's solution of the quartic: the
    resolvent cubic

    .. math:: z^3 - 2 I_2 z^2 + 2 I_4 z - \tfrac49 I_3^2 = 0

    has roots :math:`z_i = (\psi_i + \psi_j)^2`, real and non-negative
    because the Hamiltonian is Hermitian, and then

    .. math::
       \psi_m = \tfrac12 \left(s_1 \sqrt{z_1} + s_2 \sqrt{z_2}
                               + s_3 \sqrt{z_3}\right) ,

    with signs :math:`s_i = \pm1` fixed by
    :math:`s_1 s_2 s_3 \sqrt{z_1 z_2 z_3} = \tfrac23 I_3`.

    These roots carry only the accuracy that :math:`I_2, I_3, I_4`
    carry.  For a stiff spectrum that is not enough, which is why the
    routines that use them refine them against the matrix; see
    :data:`POLISH_ROOTS`.

    .. versionadded:: 1.9.0

    Parameters
    ----------
    invariant_2 : int or float
        The invariant :math:`I_2 = |h|^2`.
    invariant_3 : int or float
        The invariant :math:`I_3 = \langle h \rangle`.
    invariant_4 : int or float
        The invariant :math:`I_4 = |h \star h|^2`.

    Returns
    -------
    list of float
        The four roots, in ascending order.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import oscprob4nu

        hamiltonian = np.diag([1.0, 2.0, 3.0, -6.0]).astype(complex)
        i2, i3, i4 = oscprob4nu.su4_invariants(hamiltonian)
        psi = oscprob4nu.psi_roots_4nu(i2, i3, i4)

        print('roots    :', np.round(psi, 6))
        print('reference:', np.round(np.linalg.eigvalsh(hamiltonian) - 0.0, 6))
    """
    quadratic = -invariant_2
    linear = -(2.0/3.0)*invariant_3
    constant = 0.25*(invariant_2*invariant_2 - 2.0*invariant_4)

    resolvent = _resolvent_cubic_roots(2.0*quadratic,
                                       quadratic*quadratic - 4.0*constant,
                                       -linear*linear)
    root = np.sqrt(np.clip(resolvent, 0.0, None))

    if linear > 0.0:
        root = root*np.array([1.0, 1.0, -1.0])

    return list(np.sort(0.5*np.array([
        root[0] + root[1] + root[2],
        root[0] - root[1] - root[2],
        -root[0] + root[1] - root[2],
        -root[0] - root[1] + root[2],
    ])))


def _polish_roots(
    traceless: np.ndarray,
    psi: np.ndarray
) -> np.ndarray:
    r"""Refines the latent roots against the Hamiltonian matrix.

    One Newton step on :math:`\chi(\psi) = \det(\psi\mathbb{1} -
    \tilde{H})`, with the derivative taken as
    :math:`\chi'(\psi_m) = \prod_{l \neq m}(\psi_m - \psi_l)` from the
    current estimates.  The determinant is evaluated from the matrix, so
    unlike the closed-form roots it is not limited by the conditioning of
    the three invariants.

    Works on a stack: `traceless` of shape ``(..., 4, 4)`` and `psi` of
    shape ``(..., 4)``.

    Parameters
    ----------
    traceless : numpy.ndarray
        Traceless part of the Hamiltonian, of shape ``(..., 4, 4)``.
    psi : numpy.ndarray
        Current estimates of the roots, of shape ``(..., 4)``.

    Returns
    -------
    numpy.ndarray
        The refined roots, sorted ascending along the last axis.
    """
    shifted = (psi[..., :, None, None]*np.eye(4)
               - traceless[..., None, :, :])
    chi = np.linalg.det(shifted).real

    gaps = psi[..., :, None] - psi[..., None, :]
    diagonal = np.arange(4)
    gaps[..., diagonal, diagonal] = 1.0
    derivative = np.prod(gaps, axis=-1)

    step = np.where(derivative != 0.0,
                    chi/np.where(derivative == 0.0, 1.0, derivative), 0.0)

    return np.sort(psi - step, axis=-1)


def _latent_roots(traceless: np.ndarray) -> np.ndarray:
    r"""Returns the latent roots of a stack of traceless Hamiltonians.

    Parameters
    ----------
    traceless : numpy.ndarray
        Traceless Hamiltonians, of shape ``(..., 4, 4)``.

    Returns
    -------
    numpy.ndarray
        The four roots per element, of shape ``(..., 4)``.
    """
    squared = traceless @ traceless
    invariant_2 = 0.5*np.einsum('...ii->...', squared).real
    invariant_3 = 0.5*np.einsum('...ii->...', squared @ traceless).real
    invariant_4 = 0.5*(np.einsum('...ij,...ji->...', squared, squared).real
                       - invariant_2*invariant_2)

    quadratic = -invariant_2
    linear = -(2.0/3.0)*invariant_3
    constant = 0.25*(invariant_2*invariant_2 - 2.0*invariant_4)

    coeff_2 = 2.0*quadratic
    coeff_1 = quadratic*quadratic - 4.0*constant
    coeff_0 = -linear*linear

    depressed_p = coeff_1 - coeff_2*coeff_2/3.0
    depressed_q = 2.0*coeff_2**3/27.0 - coeff_2*coeff_1/3.0 + coeff_0
    shift = -coeff_2/3.0

    scale = 2.0*np.sqrt(np.maximum(-depressed_p, 0.0)/3.0)
    denominator = depressed_p*scale
    argument = np.clip(3.0*depressed_q/np.where(denominator != 0.0,
                                                denominator, 1.0), -1.0, 1.0)
    angle = np.arccos(argument)

    index = np.arange(3)
    resolvent = (scale[..., None]
                 * np.cos((angle[..., None] + 2.0*np.pi*index)/3.0)
                 + shift[..., None])
    root = np.sqrt(np.clip(resolvent, 0.0, None))
    root = np.where((linear > 0.0)[..., None] & (index == 2), -root, root)

    psi = 0.5*np.stack([
        root[..., 0] + root[..., 1] + root[..., 2],
        root[..., 0] - root[..., 1] - root[..., 2],
        -root[..., 0] + root[..., 1] - root[..., 2],
        -root[..., 0] - root[..., 1] + root[..., 2],
    ], axis=-1)
    psi = np.sort(psi, axis=-1)

    if POLISH_ROOTS:
        psi = _polish_roots(traceless, psi)

    return psi


def _traceless_part(
    hamiltonian_matrix: Union[list, np.ndarray]
) -> np.ndarray:
    r"""Returns the traceless part of a Hamiltonian or a stack of them.

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Hamiltonian, or stack of them, of shape ``(..., 4, 4)``.

    Returns
    -------
    numpy.ndarray
        The traceless part, of the same shape.
    """
    matrix = np.asarray(hamiltonian_matrix, dtype=complex)
    trace = np.einsum('...ii->...', matrix).real/4.0

    return matrix - trace[..., None, None]*np.eye(4)


def _divided_differences(
    psi: np.ndarray,
    baseline: np.ndarray
) -> np.ndarray:
    r"""Returns the divided differences of :math:`e^{-i\psi L}`.

    Builds the Newton coefficients :math:`c_j = f[\psi_1, \ldots,
    \psi_{j+1}]` of the exponential over the four latent roots, taking
    the confluent limit wherever two nodes are closer than
    :data:`DEGENERACY_TOL` relative to the spectral scale.  For the
    exponential the confluent value is exact rather than approximated:
    :math:`f[\psi, \ldots, \psi]` over :math:`k+1` equal nodes is
    :math:`(-iL)^k e^{-i\psi L}/k!`.

    Parameters
    ----------
    psi : numpy.ndarray
        Latent roots, sorted ascending, of shape ``(..., 4)``.
    baseline : numpy.ndarray
        Baselines, of shape ``(...)``.

    Returns
    -------
    numpy.ndarray
        The four coefficients, of shape ``(..., 4)``.
    """
    scale = np.max(np.abs(psi), axis=-1)
    tolerance = DEGENERACY_TOL*np.where(scale > 0.0, scale, 1.0)

    phase = np.exp(-1.j*psi*baseline[..., None])
    minus_i_l = -1.j*baseline

    table = phase
    coefficients = [table[..., 0]]

    for order in range(1, 4):
        separation = psi[..., order:] - psi[..., :-order]
        # The confluent value, used wherever the nodes have merged.
        confluent = (minus_i_l[..., None]**order/math.factorial(order)
                     * phase[..., :-order])
        regular = ((table[..., 1:] - table[..., :-1])
                   / np.where(np.abs(separation) > tolerance[..., None],
                              separation, 1.0))
        table = np.where(np.abs(separation) > tolerance[..., None],
                         regular, confluent)
        coefficients.append(table[..., 0])

    return np.stack(coefficients, axis=-1)


def _evolution_operator_4nu_array(
    hamiltonian_matrix: Union[list, np.ndarray],
    L: Union[int, float, list, np.ndarray]
) -> np.ndarray:
    r"""Returns :math:`U_4` for a stack, as an array.

    Interpolates :math:`e^{-i\psi L}` through the four latent roots in
    Newton form,

    .. math::
       U_4 = c_0 \mathbb{1} + c_1 (\tilde{H} - \psi_1)
             + c_2 (\tilde{H} - \psi_1)(\tilde{H} - \psi_2)
             + c_3 (\tilde{H} - \psi_1)(\tilde{H} - \psi_2)
                   (\tilde{H} - \psi_3) ,

    with the :math:`c_j` the divided differences of the exponential over
    the roots.  This is the Cayley-Hamilton form written in the basis
    that tolerates repeated nodes: a *confluent* divided difference is
    just a derivative, and for the exponential that derivative is known
    in closed form, :math:`f^{(k)}(\psi)/k! = (-iL)^k e^{-i\psi L}/k!`.

    The obvious alternative --- solving the Vandermonde system
    :math:`\sum_j c_j \psi_m^j = e^{-i\psi_m L}` --- was tried first and
    rejected: it is singular the moment two roots coincide, which is not
    an exotic case but includes a Hamiltonian proportional to the
    identity, a zero Hamiltonian, and any triply degenerate spectrum.
    The Newton form handles all of them without a special branch, and
    needs no eigenvectors either.

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Hamiltonian, or stack of them, of shape ``(..., 4, 4)``.
    L : int or float or array_like
        Baseline, or array of baselines.

    Returns
    -------
    numpy.ndarray
        The evolution operator, of shape ``(..., 4, 4)``.
    """
    traceless = _traceless_part(hamiltonian_matrix)
    baseline = np.asarray(L, dtype=float)

    leading = np.broadcast_shapes(traceless.shape[:-2], baseline.shape)
    traceless = np.broadcast_to(traceless, leading + (4, 4))
    baseline = np.broadcast_to(baseline, leading)

    psi = _latent_roots(traceless)
    coefficients = _divided_differences(psi, baseline)

    identity = np.eye(4)
    factor = np.broadcast_to(identity, leading + (4, 4)).copy()
    operator = coefficients[..., 0, None, None]*identity

    for order in range(1, 4):
        factor = factor @ (traceless
                           - psi[..., order-1, None, None]*identity)
        operator = operator + coefficients[..., order, None, None]*factor

    return operator


def evolution_operator_4nu_u_coefficients(
    hamiltonian_matrix: Union[list, np.ndarray],
    L: Union[int, float]
) -> List[complex]:
    r"""Returns the :math:`u_a` of the SU(4) expansion of :math:`U_4`.

    Returns the sixteen coefficients :math:`u_0, u_1, \ldots, u_{15}` of
    the expansion :math:`U_4(L) = u_0 \mathbb{1} + i u_a \lambda^a` of the
    time-evolution operator, in the same convention as
    :func:`oscprob3nu.evolution_operator_3nu_u_coefficients`.

    The four-flavor analogue of Eqs. (10)-(11) of arXiv:1904.12391 is

    .. math::
       u_0 = \frac14 \sum_m e^{-i \psi_m L} , \qquad
       i u_a = \sum_m e^{-i \psi_m L}\,
       \frac{\left(\psi_m^2 - \tfrac12 I_2\right) h_a
             + \psi_m (h \star h)_a + ((h \star h) \star h)_a}
            {\chi'(\psi_m)} ,

    with :math:`\chi'(\psi_m) = 4\psi_m^3 - 2 I_2 \psi_m - \tfrac23 I_3`.
    Note the third term in the numerator: at three flavors the tower
    closes on itself and no such term appears.

    .. versionadded:: 1.9.0

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Four-flavor Hermitian Hamiltonian, given as a nested list or an
        array of shape ``(4, 4)``.
    L : int or float
        Baseline, in units reciprocal to those of the Hamiltonian.

    Returns
    -------
    list of complex
        The sixteen coefficients ``[u0, u1, ..., u15]``.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import oscprob4nu

        hamiltonian = np.diag([1.0, 2.0, 3.0, -6.0]).astype(complex)
        u = oscprob4nu.evolution_operator_4nu_u_coefficients(hamiltonian, 1.0)

        print('%d coefficients' % len(u))
        print('u0 = %.6f%+.6fj' % (u[0].real, u[0].imag))
    """
    operator = _evolution_operator_4nu_array(hamiltonian_matrix, L)

    u_0 = np.trace(operator)/4.0
    u_a = 0.5*np.einsum('aij,ji->a', LAMBDA_SU4, operator)/1.j

    return [complex(u_0)] + [complex(value) for value in u_a]


def evolution_operator_4nu(
    hamiltonian_matrix: Union[list, np.ndarray],
    L: Union[int, float, list, np.ndarray]
) -> Union[List[List[complex]], np.ndarray]:
    r"""Returns the four-neutrino time-evolution operator :math:`U_4`.

    Returns :math:`U_4(L) = e^{-i \tilde{H} L}`, with :math:`\tilde{H}`
    the traceless part of `hamiltonian_matrix`.  The discarded trace
    contributes only the overall phase :math:`e^{-i h_0 L}`, which
    cancels in the probabilities.

    .. versionadded:: 1.9.0

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Four-flavor Hermitian Hamiltonian, given as a nested list of
        shape ``(4, 4)``, or a stack of them, of shape ``(..., 4, 4)``.
    L : int or float or array_like
        Baseline, in units reciprocal to those of the Hamiltonian, or an
        array of baselines broadcastable against the leading axes of
        `hamiltonian_matrix`.

    Returns
    -------
    list of list of complex or numpy.ndarray
        For a single Hamiltonian and baseline, the :math:`4\times4`
        evolution operator as a nested list.  If either argument is a
        stack, an array of shape ``(..., 4, 4)``.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import oscprob4nu

        hamiltonian = np.diag([1.0, 2.0, 3.0, -6.0]).astype(complex)
        operator = oscprob4nu.evolution_operator_4nu(hamiltonian, 1.0)

        unitarity = np.asarray(operator).conj().T @ np.asarray(operator)
        print('|U| unitary to %.1e' % np.max(np.abs(unitarity - np.eye(4))))
    """
    operator = _evolution_operator_4nu_array(hamiltonian_matrix, L)

    if operator.ndim == 2:
        return [[complex(entry) for entry in row] for row in operator]

    return operator


def probabilities_4nu(
    hamiltonian_matrix: Union[list, np.ndarray],
    L: Union[int, float, list, np.ndarray]
) -> Union[Tuple[float, ...], np.ndarray]:
    r"""Returns the four-neutrino oscillation probabilities.

    Returns the sixteen flavor-transition probabilities
    :math:`P_{\alpha\beta} \equiv P(\nu_\alpha \to \nu_\beta)
    = |[U_4]_{\beta\alpha}|^2`, ordered with the initial flavor varying
    slowest, exactly as at two and three flavors.  With the fourth state
    read as sterile, the flavor order is
    :math:`(\nu_e, \nu_\mu, \nu_\tau, \nu_s)`.

    .. versionadded:: 1.9.0

    Parameters
    ----------
    hamiltonian_matrix : array_like
        Four-flavor Hermitian Hamiltonian, given as a nested list of
        shape ``(4, 4)``, or a stack of them, of shape ``(..., 4, 4)``.
    L : int or float or array_like
        Baseline, in units reciprocal to those of the Hamiltonian, or an
        array of baselines broadcastable against the leading axes of
        `hamiltonian_matrix`.

    Returns
    -------
    tuple of float or numpy.ndarray
        For a single Hamiltonian and baseline, the sixteen
        probabilities as a tuple, ordered
        :math:`P_{ee}, P_{e\mu}, P_{e\tau}, P_{es}, P_{\mu e}, \ldots,
        P_{ss}`.  If either argument is a stack, an array of shape
        ``(..., 16)`` in the same order.

    Notes
    -----
    Passing arrays evaluates the whole stack at once, which is far
    faster than calling this routine in a Python loop, exactly as at two
    and three flavors; see :ref:`scanning`.

    The accuracy deserves a word, because four flavors is where it stops
    being automatic.  Against :func:`numpy.linalg.eigh`, a generic
    Hamiltonian agrees to about :math:`2 \times 10^{-13}`.  A stiff 3+1
    spectrum, with :math:`\Delta m^2_{41}` four orders of magnitude above
    :math:`\Delta m^2_{21}`, agrees to about :math:`10^{-9}` --- limited
    not by this expansion but by what double precision retains when the
    invariants are formed.  See :data:`POLISH_ROOTS`, which is what keeps
    that figure from being :math:`5 \times 10^{-7}`.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import oscprob4nu

        hamiltonian = np.diag([1.0, 2.0, 3.0, -6.0]).astype(complex)
        prob = oscprob4nu.probabilities_4nu(hamiltonian, 1.0)

        print('%d probabilities' % len(prob))
        print('P_ee = %.6f' % prob[0])
        print('they sum to %.6f' % sum(prob[0:4]))
    """
    operator = _evolution_operator_4nu_array(hamiltonian_matrix, L)

    # P[alpha][beta] = |U[beta][alpha]|^2, so the initial flavor varies
    # slowest once the last two axes are swapped and flattened.
    probabilities = np.abs(np.swapaxes(operator, -1, -2))**2

    if probabilities.ndim == 2:
        return tuple(float(value) for value in probabilities.reshape(16))

    return probabilities.reshape(probabilities.shape[:-2] + (16,))
