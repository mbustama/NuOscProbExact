# -*- coding: utf-8 -*-
r"""Compute three-neutrino Hamiltonians for selected scenarios.

This module contains routines that build the three-neutrino Hamiltonian
for a number of standard scenarios --- oscillations in vacuum, in
matter of constant density, in matter with non-standard interactions
(NSI), and in a CPT-odd Lorentz invariance-violating (LIV) background
--- together with the textbook oscillation formula in vacuum, which
serves to validate the exact SU(3) computation performed by
:mod:`oscprob3nu` and described in [1]_.

The Hamiltonians built here are meant to be passed to
:func:`oscprob3nu.probabilities_3nu`.  They are examples: the exact
computation accepts *any* Hermitian :math:`3\times3` Hamiltonian.

The routines that take a neutrino energy also accept an *array* of
energies, and then return one Hamiltonian per energy, stacked along a
leading axis.  That stack is exactly what :func:`oscprob3nu.probabilities_3nu`
expects, so a whole energy scan is two calls rather than a loop.

Units
-----

Throughout this module,

===========================  ==================================
Quantity                     Units
===========================  ==================================
Mass-squared differences     eV\ :sup:`2`
Neutrino energy              eV
Baseline                     eV\ :sup:`-1`
Matter potential             eV
LIV eigenvalues and scale    eV
CP-violation phases          radian
===========================  ==================================

The routine
:func:`hamiltonian_3nu_vacuum_energy_independent` returns the
energy-*independent* part of the vacuum Hamiltonian, i.e. it has units
of eV\ :sup:`2` and must still be divided by the neutrino energy.  The
module :mod:`globaldefs` provides ``CONV_KM_TO_INV_EV`` to convert a
baseline in km into eV\ :sup:`-1`.

Sign convention
---------------

The vacuum Hamiltonian is
:math:`H_{\rm vac} = U M^2 U^\dagger / (2E)`, with
:math:`M^2 = \mathrm{diag}(0, \Delta m^2_{21}, \Delta m^2_{31})` and
:math:`U` the PMNS matrix, so that adding a positive matter potential to
the :math:`ee` entry describes neutrinos, not antineutrinos.

Routine listings
----------------

    * pmns_mixing_matrix - Returns the :math:`3\times3` PMNS matrix
    * hamiltonian_3nu_vacuum_energy_independent - Returns :math:`H_{\rm vac}` without the :math:`1/E`
    * delta - Kronecker delta
    * J - Product of four entries of the PMNS matrix
    * probabilities_3nu_vacuum_std - Vacuum probabilities, standard formula
    * hamiltonian_3nu_matter - Returns :math:`H_{\rm matter}`
    * hamiltonian_3nu_nsi - Returns :math:`H_{\rm NSI}`
    * hamiltonian_3nu_liv - Returns :math:`H_{\rm LIV}`

References
----------

.. [1] Mauricio Bustamante, "Exact neutrino oscillation probabilities
   with arbitrary time-independent Hamiltonians", arXiv:1904.12391.

Created: 2019/04/17 17:14
Last modified: 2026/07/31
"""

from __future__ import print_function

__version__ = "1.1"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

__all__ = ['pmns_mixing_matrix',
           'hamiltonian_3nu_vacuum_energy_independent', 'delta', 'J',
           'probabilities_3nu_vacuum_std', 'hamiltonian_3nu_matter',
           'hamiltonian_3nu_nsi', 'hamiltonian_3nu_liv']

import numpy as np


_EE_PROJECTOR = np.array([[1.0, 0.0, 0.0],
                          [0.0, 0.0, 0.0],
                          [0.0, 0.0, 0.0]], dtype=complex)
r"""numpy.ndarray: Module-level constant.

The matrix that the charged-current matter potential multiplies: it
selects the :math:`ee` entry, since only :math:`\nu_e` interacts through
charged currents with the electrons in matter.
"""


def pmns_mixing_matrix(s12, s23, s13, dCP):
    r"""Returns the :math:`3\times3` PMNS mixing matrix.

    Computes and returns the complex :math:`3\times3` PMNS mixing
    matrix, parametrized by the three rotation angles
    :math:`\theta_{12}`, :math:`\theta_{23}`, :math:`\theta_{13}`, and
    the CP-violation phase :math:`\delta_{\rm CP}`, in the standard
    PDG parametrization.

    Parameters
    ----------
    s12 : float
        :math:`\sin\theta_{12}`.
    s23 : float
        :math:`\sin\theta_{23}`.
    s13 : float
        :math:`\sin\theta_{13}`.
    dCP : float
        CP-violation phase :math:`\delta_{\rm CP}` [radian].

    Returns
    -------
    list of list of complex
        The :math:`3\times3` PMNS mixing matrix, as a nested list.

    Examples
    --------
    >>> U = pmns_mixing_matrix(0.55, 0.76, 0.15, 0.0)
    >>> print('%.6f  %.6f' % (U[0][0].real, U[0][1].real))
    0.825716  0.543777
    """
    c12 = np.sqrt(1.0-s12*s12)
    c23 = np.sqrt(1.0-s23*s23)
    c13 = np.sqrt(1.0-s13*s13)

    cdCP = np.cos(dCP)
    sdCP = np.sin(dCP)

    U00 = c12*c13
    U01 = s12*c13
    U02 = s13*complex(cdCP, -sdCP)
    U10 = -s12*c23 - c12*s23*s13*complex(cdCP, sdCP)
    U11 = c12*c23 - s12*s23*s13*complex(cdCP, sdCP)
    U12 = s23*c13
    U20 = s12*s23 - c12*c23*s13*complex(cdCP, sdCP)
    U21 = -c12*s23 - s12*c23*s13*complex(cdCP, sdCP)
    U22 = c23*c13

    return [[complex(U00), complex(U01), U02],
            [U10, U11, complex(U12)],
            [U20, U21, complex(U22)]]


def hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31,
                                              compute_matrix_multiplication=False):
    r"""Returns the three-neutrino Hamiltonian for vacuum oscillations.

    Computes and returns the energy-independent part of the complex
    :math:`3\times3` three-neutrino Hamiltonian for oscillations in
    vacuum, parametrized by three mixing angles ---
    :math:`\theta_{12}`, :math:`\theta_{23}`, :math:`\theta_{13}` ---
    one CP-violation phase --- :math:`\delta_{\rm CP}` --- and two
    mass-squared differences --- :math:`\Delta m^2_{21}`,
    :math:`\Delta m^2_{31}`.  The Hamiltonian is
    :math:`H = \frac{1}{2} U M^2 U^\dagger`, with :math:`U` the PMNS
    matrix and :math:`M^2 = \mathrm{diag}(0, \Delta m^2_{21},
    \Delta m^2_{31})` the mass matrix.  The multiplicative factor
    :math:`1/E` is *not* applied.

    Parameters
    ----------
    s12 : float
        :math:`\sin\theta_{12}`.
    s23 : float
        :math:`\sin\theta_{23}`.
    s13 : float
        :math:`\sin\theta_{13}`.
    dCP : float
        CP-violation phase :math:`\delta_{\rm CP}` [radian].
    D21 : float
        Mass-squared difference :math:`\Delta m^2_{21}`
        [eV\ :sup:`2`].
    D31 : float
        Mass-squared difference :math:`\Delta m^2_{31}`
        [eV\ :sup:`2`].
    compute_matrix_multiplication : bool, optional
        If ``False`` (default), use the pre-computed closed-form
        expressions; if ``True``, carry out the matrix multiplication
        :math:`U M^2 U^\dagger` explicitly.  Both give the same result;
        the option exists as a cross-check.

    Returns
    -------
    numpy.ndarray
        The :math:`3\times3` complex Hamiltonian [eV\ :sup:`2`], to be
        divided by the neutrino energy before use.

    Examples
    --------
    >>> H = hamiltonian_3nu_vacuum_energy_independent(0.55, 0.76, 0.15,
    ...                                               0.0, 7.4e-5, 2.5e-3)
    >>> print('%.6e' % H[0][0].real)
    3.906567e-05
    """
    c12 = np.sqrt(1.0-s12*s12)
    c23 = np.sqrt(1.0-s23*s23)
    c13 = np.sqrt(1.0-s13*s13)

    f = 1./2.

    if not compute_matrix_multiplication:

        # All Hij have units of [eV^2]
        H00 = c13*c13*D21*s12*s12 + D31*s13*s13
        H01 = c12*c13*c23*D21*s12 + \
            c13*(D31-D21*s12*s12)*s13*s23*complex(np.cos(dCP), -np.sin(dCP))
        H02 = c13*c23*(D31-D21*s12*s12)*s13*complex(np.cos(dCP),
                                                    -np.sin(dCP)) - \
            c12*c13*D21*s12*s23
        H10 = c12*c13*c23*D21*s12 + \
            c13*(D31-D21*s12*s12)*s13*s23*complex(np.cos(dCP), np.sin(dCP))
        H11 = c12*c12*c23*c23*D21 + \
            (c13*c13*D31 + D21*s12*s12*s13*s13)*s23*s23 - \
            2.0*c12*c23*D21*s12*s13*s23*np.cos(dCP)
        H12 = c13*c13*c23*D31*s23 + \
            (c23*s12*s13*complex(np.cos(dCP), -np.sin(dCP)) + c12*s23) * \
            (-c12*c23*D21 + D21*s12*s13*s23*complex(np.cos(dCP),
                                                    np.sin(dCP)))
        H20 = c13*c23*(D31-D21*s12*s12)*s13*complex(np.cos(dCP),
                                                    np.sin(dCP)) - \
            c12*c13*D21*s12*s23
        H21 = c13*c13*c23*D31*s23 - \
            D21*(c23*s12*s13*complex(np.cos(dCP), np.sin(dCP)) + c12*s23) * \
            (c12*c23 - s12*s13*s23*complex(np.cos(dCP), -np.sin(dCP)))
        H22 = c23*c23*(c13*c13*D31 + D21*s12*s12*s13*s13) + \
            c12*c12*D21*s23*s23 + \
            2.0*c12*c23*D21*s12*s13*s23*np.cos(dCP)

        H = np.array([[H00*f, H01*f, H02*f],
                      [H10*f, H11*f, H12*f],
                      [H20*f, H21*f, H22*f]], dtype=complex)

    else:

        # PMNS matrix
        U = np.array(pmns_mixing_matrix(s12, s23, s13, dCP))
        # Mass matrix
        M2 = np.array([[0.0, 0.0, 0.0], [0.0, D21, 0.0], [0.0, 0.0, D31]])
        # Hamiltonian
        H = (f*(U @ M2 @ U.conj().T)).astype(complex)

    return H


def delta(a, b):
    r"""Returns the Kronecker delta :math:`\delta_{ab}`.

    Parameters
    ----------
    a : int
        First index.
    b : int
        Second index.

    Returns
    -------
    int
        1 if ``a == b``, 0 otherwise.

    Examples
    --------
    >>> print(delta(0, 0), delta(0, 1))
    1 0
    """
    if (a == b):
        return 1
    else:
        return 0


def J(U, alpha, beta, k, j):
    r"""Returns the quartic product of PMNS matrix entries.

    Returns :math:`J = U_{\alpha k}^* U_{\beta k} U_{\alpha j}
    U_{\beta j}^*`, with :math:`U` the PMNS mixing matrix.  This
    product appears in the standard expression for the three-neutrino
    oscillation probability in vacuum; its imaginary part is the
    Jarlskog invariant, up to a sign.

    Parameters
    ----------
    U : array_like
        The :math:`3\times3` complex PMNS mixing matrix.
    alpha : int
        Index of the initial flavor (0: :math:`e`, 1: :math:`\mu`,
        2: :math:`\tau`).
    beta : int
        Index of the final flavor (0: :math:`e`, 1: :math:`\mu`,
        2: :math:`\tau`).
    k : int
        First index of the sum over mass eigenstates (0, 1, or 2).
    j : int
        Second index of the sum over mass eigenstates (0, 1, or 2).

    Returns
    -------
    complex
        The product :math:`U_{\alpha k}^* U_{\beta k} U_{\alpha j}
        U_{\beta j}^*`.

    Examples
    --------
    >>> U = pmns_mixing_matrix(0.55, 0.76, 0.15, 0.0)
    >>> print('%.6f' % J(U, 0, 1, 1, 0).real)
    -0.097579
    """
    return np.conj(U[alpha][k])*U[beta][k]*U[alpha][j]*np.conj(U[beta][j])


def probabilities_3nu_vacuum_std(U, D21, D31, energy, L):
    r"""Returns the 3nu vacuum probabilities, standard computation.

    Returns the probabilities for three-neutrino oscillations in vacuum,
    computed with the standard analytical expression

    .. math::
       P_{\alpha\beta} = \delta_{\alpha\beta}
       - 4 \sum_{k>j} \mathrm{Re}(J_{kj})
         \sin^2\left(\frac{\Delta m^2_{kj} L}{4E}\right)
       + 2 \sum_{k>j} \mathrm{Im}(J_{kj})
         \sin\left(\frac{\Delta m^2_{kj} L}{2E}\right) .

    This routine exists to validate the exact SU(3) computation in
    :mod:`oscprob3nu`; the two agree to round-off.

    Parameters
    ----------
    U : array_like
        The :math:`3\times3` complex PMNS mixing matrix, as returned by
        `pmns_mixing_matrix`.
    D21 : float
        Mass-squared difference :math:`\Delta m^2_{21}`
        [eV\ :sup:`2`].
    D31 : float
        Mass-squared difference :math:`\Delta m^2_{31}`
        [eV\ :sup:`2`].
    energy : float
        Neutrino energy [eV].
    L : float
        Baseline [eV\ :sup:`-1`].

    Returns
    -------
    list of float
        The nine probabilities ``[Pee, Pem, Pet, Pme, Pmm, Pmt, Pte,
        Ptm, Ptt]``, ordered with the initial flavor varying slowest.

    See Also
    --------
    oscprob3nu.probabilities_3nu : The exact SU(3) computation.

    Examples
    --------
    >>> U = pmns_mixing_matrix(0.55, 0.76, 0.15, 0.0)
    >>> prob = probabilities_3nu_vacuum_std(U, 7.4e-5, 2.5e-3, 1.0e9,
    ...                                     5.0e12)
    >>> print('%.6f  %.6f' % (prob[0], prob[1]))
    0.992787  0.001981
    """
    D32 = D31-D21

    # Oscillation phases, 2 * Dm2 * L / (4 * E)
    arg21 = D21*L/2.0/energy
    arg31 = D31*L/2.0/energy
    arg32 = D32*L/2.0/energy

    s21 = np.sin(arg21)
    s31 = np.sin(arg31)
    s32 = np.sin(arg32)
    ss21 = pow(np.sin(arg21/2.0), 2.0)
    ss31 = pow(np.sin(arg31/2.0), 2.0)
    ss32 = pow(np.sin(arg32/2.0), 2.0)

    # Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt
    prob = [delta(alpha, beta)
            - 4.0 * (J(U, alpha, beta, 1, 0).real*ss21
                     + J(U, alpha, beta, 2, 0).real*ss31
                     + J(U, alpha, beta, 2, 1).real*ss32)
            + 2.0 * (J(U, alpha, beta, 1, 0).imag*s21
                     + J(U, alpha, beta, 2, 0).imag*s31
                     + J(U, alpha, beta, 2, 1).imag*s32)
            for alpha in [0, 1, 2] for beta in [0, 1, 2]]

    return prob


def hamiltonian_3nu_matter(h_vacuum_energy_independent, energy, VCC):
    r"""Returns the three-neutrino Hamiltonian for matter oscillations.

    Computes and returns the :math:`3\times3` three-neutrino
    Hamiltonian for oscillations in matter of constant density, obtained
    by adding the charged-current matter potential to the :math:`ee`
    entry of the vacuum Hamiltonian.

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent part of the three-neutrino vacuum
        Hamiltonian [eV\ :sup:`2`], as returned by
        `hamiltonian_3nu_vacuum_energy_independent`.  It is not
        modified.
    energy : float or array_like
        Neutrino energy [eV], or an array of energies, in which case one
        Hamiltonian is returned per energy.
    VCC : float
        Potential due to charged-current interactions of
        :math:`\nu_e` with electrons [eV].  Positive for neutrinos,
        negative for antineutrinos.

    Returns
    -------
    numpy.ndarray
        The :math:`3\times3` complex Hamiltonian [eV], of shape
        ``(3, 3)`` for a scalar energy and ``(..., 3, 3)`` for an
        array of energies.

    Examples
    --------
    >>> H_vac = hamiltonian_3nu_vacuum_energy_independent(0.55, 0.76, 0.15,
    ...                                                   0.0, 7.4e-5,
    ...                                                   2.5e-3)
    >>> H = hamiltonian_3nu_matter(H_vac, 1.0e9, 1.0e-13)
    >>> print('%.6e' % H[0][0].real)
    1.390657e-13
    """
    h_vacuum = np.asarray(h_vacuum_energy_independent, dtype=complex)
    energy = np.asarray(energy, dtype=float)
    VCC = np.asarray(VCC, dtype=float)

    # Indexing the energy with two trailing axes lets a scalar energy
    # return a single 3x3 matrix and an array of energies return one
    # matrix per energy, through the same expression
    return h_vacuum/energy[..., None, None] \
        + VCC[..., None, None]*_EE_PROJECTOR


def hamiltonian_3nu_nsi(h_vacuum_energy_independent, energy, VCC, eps):
    r"""Returns the three-neutrino Hamiltonian for oscillations w/ NSI.

    Computes and returns the :math:`3\times3` three-neutrino
    Hamiltonian for oscillations with non-standard interactions (NSI)
    in matter of constant density.

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent part of the three-neutrino vacuum
        Hamiltonian [eV\ :sup:`2`], as returned by
        `hamiltonian_3nu_vacuum_energy_independent`.  It is not
        modified.
    energy : float or array_like
        Neutrino energy [eV], or an array of energies, in which case one
        Hamiltonian is returned per energy.
    VCC : float
        Potential due to charged-current interactions of
        :math:`\nu_e` with electrons [eV].
    eps : array_like
        The six NSI strength parameters ``[eps_ee, eps_em, eps_et,
        eps_mm, eps_mt, eps_tt]``, adimensional.  The diagonal
        parameters are real; the off-diagonal ones may be complex, and
        their complex conjugates are placed in the lower off-diagonal
        entries so that the Hamiltonian stays Hermitian.

    Returns
    -------
    numpy.ndarray
        The :math:`3\times3` complex Hamiltonian [eV], of shape
        ``(3, 3)`` for a scalar energy and ``(..., 3, 3)`` for an
        array of energies.

    Examples
    --------
    >>> H_vac = hamiltonian_3nu_vacuum_energy_independent(0.55, 0.76, 0.15,
    ...                                                   0.0, 7.4e-5,
    ...                                                   2.5e-3)
    >>> H = hamiltonian_3nu_nsi(H_vac, 1.0e9, 1.0e-13,
    ...                         [0.06, -0.06+0.03j, 0.0, 1.2, 0.0, 0.0])
    >>> print('%+.6e%+.6ej' % (H[0][1].real, H[0][1].imag))
    +1.445471e-13+3.000000e-15j
    """
    h_vacuum = np.asarray(h_vacuum_energy_independent, dtype=complex)
    energy = np.asarray(energy, dtype=float)
    VCC = np.asarray(VCC, dtype=float)

    eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt = eps

    # The matrix is complex so that complex off-diagonal parameters keep
    # their imaginary parts; a real one would silently discard them
    nsi = np.array([[1.0+eps_ee, eps_em, eps_et],
                    [np.conj(eps_em), eps_mm, eps_mt],
                    [np.conj(eps_et), np.conj(eps_mt), eps_tt]],
                   dtype=complex)

    return h_vacuum/energy[..., None, None] + VCC[..., None, None]*nsi


def hamiltonian_3nu_liv(h_vacuum_energy_independent, energy, sxi12, sxi23,
                        sxi13, dxiCP, b1, b2, b3, Lambda):
    r"""Returns the three-neutrino Hamiltonian for oscillations w/ LIV.

    Computes and returns the :math:`3\times3` three-neutrino
    Hamiltonian for oscillations in a CPT-odd Lorentz
    invariance-violating (LIV) background.  The LIV term is
    :math:`(E/\Lambda) R B_3 R^\dagger`, with
    :math:`B_3 = \mathrm{diag}(b_1, b_2, b_3)` and :math:`R` a
    PMNS-like matrix built from the mixing angles :math:`\xi_{ij}` and
    the phase :math:`\delta_{\xi,\rm CP}` that relate the eigenvectors
    of :math:`B_3` to the flavor states.

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent part of the three-neutrino vacuum
        Hamiltonian [eV\ :sup:`2`], as returned by
        `hamiltonian_3nu_vacuum_energy_independent`.  It is not
        modified.
    energy : float or array_like
        Neutrino energy [eV], or an array of energies, in which case one
        Hamiltonian is returned per energy.
    sxi12 : float
        :math:`\sin\xi_{12}`, with :math:`\xi_{12}` one of the mixing
        angles between the space of the eigenvectors of :math:`B_3` and
        the flavor states.
    sxi23 : float
        :math:`\sin\xi_{23}`, likewise.
    sxi13 : float
        :math:`\sin\xi_{13}`, likewise.
    dxiCP : float
        CP-violation phase of the LIV operator :math:`B_3` [radian].
    b1 : float
        Eigenvalue :math:`b_1` of the LIV operator :math:`B_3` [eV].
    b2 : float
        Eigenvalue :math:`b_2` of the LIV operator :math:`B_3` [eV].
    b3 : float
        Eigenvalue :math:`b_3` of the LIV operator :math:`B_3` [eV].
    Lambda : float
        Energy scale :math:`\Lambda` of the LIV operator :math:`B_3`
        [eV].

    Returns
    -------
    numpy.ndarray
        The :math:`3\times3` complex Hamiltonian [eV], of shape
        ``(3, 3)`` for a scalar energy and ``(..., 3, 3)`` for an
        array of energies.

    Examples
    --------
    >>> H_vac = hamiltonian_3nu_vacuum_energy_independent(0.55, 0.76, 0.15,
    ...                                                   0.0, 7.4e-5,
    ...                                                   2.5e-3)
    >>> H = hamiltonian_3nu_liv(H_vac, 1.0e9, 0.3, 0.4, 0.5, 0.7, 1.0e-9,
    ...                         1.5e-9, 2.0e-9, 1.0e12)
    >>> print('%.6e' % H[0][0].real)
    1.322816e-12
    """
    h_vacuum = np.asarray(h_vacuum_energy_independent, dtype=complex)
    energy = np.asarray(energy, dtype=float)

    # PMNS-like mixing matrix
    R = np.array(pmns_mixing_matrix(sxi12, sxi23, sxi13, dxiCP))
    # B matrix
    B = np.array([[b1, 0.0, 0.0], [0.0, b2, 0.0], [0.0, 0.0, b3]])
    # LIV term, R.B3.R^dagger
    liv = R @ B @ R.conj().T

    f = energy/Lambda

    return h_vacuum/energy[..., None, None] + f[..., None, None]*liv
