# -*- coding: utf-8 -*-
r"""Compute two-neutrino Hamiltonians for selected scenarios.

This module contains routines that build the two-neutrino Hamiltonian
for a number of standard scenarios --- oscillations in vacuum, in
matter of constant density, in matter with non-standard interactions
(NSI), and in a CPT-odd Lorentz invariance-violating (LIV) background
--- together with the textbook oscillation formulas for vacuum and
matter, which serve to validate the exact SU(2) computation performed by
:mod:`oscprob2nu` and described in [1]_.

The Hamiltonians built here are meant to be passed to
:func:`oscprob2nu.probabilities_2nu`.  They are examples: the exact
computation accepts *any* Hermitian :math:`2\times2` Hamiltonian.

The routines that take a neutrino energy also accept an *array* of
energies, and then return one Hamiltonian per energy, stacked along a
leading axis.  That stack is exactly what :func:`oscprob2nu.probabilities_2nu`
expects, so a whole energy scan is two calls rather than a loop.

Units
-----

Throughout this module,

===========================  ==================================
Quantity                     Units
===========================  ==================================
Mass-squared difference      eV\ :sup:`2`
Neutrino energy              eV
Baseline                     eV\ :sup:`-1`
Matter potential             eV
LIV eigenvalues and scale    eV
===========================  ==================================

The routine
:func:`hamiltonian_2nu_vacuum_energy_independent` returns the
energy-*independent* part of the vacuum Hamiltonian, i.e. it has units
of eV\ :sup:`2` and must still be divided by the neutrino energy.  The
module :mod:`globaldefs` provides ``CONV_KM_TO_INV_EV`` to convert a
baseline in km into eV\ :sup:`-1`.

Sign convention
---------------

The vacuum Hamiltonian is

.. math::
   H_{\rm vac} = \frac{\Delta m^2}{4E}
   \begin{pmatrix} -\cos 2\theta & \sin 2\theta \\
                    \sin 2\theta & \cos 2\theta \end{pmatrix} ,

i.e. the mass eigenstate with the larger mass-squared value is the
second one.  This sign matters: it fixes the sign of the matter
potential *relative* to the vacuum term, and hence whether the routines
describe neutrinos (as they do) or antineutrinos.  An overall sign flip
of the vacuum Hamiltonian alone is invisible in vacuum but moves the
Mikheyev-Smirnov-Wolfenstein resonance from neutrinos to antineutrinos.

Routine listings
----------------

    * mixing_matrix_2nu - Returns the :math:`2\times2` rotation matrix
    * hamiltonian_2nu_vacuum_energy_independent - Returns :math:`H_{\rm vac}` without the :math:`1/E`
    * probabilities_2nu_vacuum_std - Vacuum probabilities, standard formula
    * hamiltonian_2nu_matter - Returns :math:`H_{\rm matter}`
    * probabilities_2nu_matter_std - Matter probabilities, standard formula
    * hamiltonian_2nu_nsi - Returns :math:`H_{\rm NSI}`
    * hamiltonian_2nu_liv - Returns :math:`H_{\rm LIV}`

References
----------

.. [1] Mauricio Bustamante, "Exact neutrino oscillation probabilities
   with arbitrary time-independent Hamiltonians", arXiv:1904.12391.

Created: 2019/04/21 15:00
Last modified: 2026/07/31
"""

from __future__ import print_function

__version__ = "1.1"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

__all__ = ['mixing_matrix_2nu', 'hamiltonian_2nu_vacuum_energy_independent',
           'probabilities_2nu_vacuum_std', 'hamiltonian_2nu_matter',
           'probabilities_2nu_matter_std', 'hamiltonian_2nu_nsi',
           'hamiltonian_2nu_liv']

import numpy as np


_EE_PROJECTOR = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=complex)
r"""numpy.ndarray: Module-level constant.

The matrix that the charged-current matter potential multiplies: it
selects the :math:`ee` entry, since only :math:`\nu_e` interacts through
charged currents with the electrons in matter.
"""


def mixing_matrix_2nu(sth):
    r"""Returns the :math:`2\times2` rotation matrix.

    Computes and returns the real :math:`2\times2` rotation matrix
    parametrized by a single rotation angle :math:`\theta`.

    Parameters
    ----------
    sth : float
        :math:`\sin\theta`, with :math:`\theta` in the first quadrant,
        so that :math:`\cos\theta = \sqrt{1-\sin^2\theta} \geq 0`.

    Returns
    -------
    list of list of float
        The rotation matrix ``[[cth, sth], [-sth, cth]]``, with
        ``cth`` = :math:`\cos\theta`.

    Examples
    --------
    >>> R = mixing_matrix_2nu(0.6)
    >>> print('%.6f  %.6f' % (R[0][0], R[0][1]))
    0.800000  0.600000
    """
    cth = np.sqrt(1.0-sth*sth)

    return [[cth, sth], [-sth, cth]]


def hamiltonian_2nu_vacuum_energy_independent(sth, Dm2,
                                              compute_matrix_multiplication=False):
    r"""Returns the two-neutrino Hamiltonian for vacuum oscillations.

    Computes and returns the energy-independent part of the real
    :math:`2\times2` two-neutrino Hamiltonian for oscillations in
    vacuum, parametrized by a single mixing angle :math:`\theta` and a
    single mass-squared difference :math:`\Delta m^2`.  The Hamiltonian
    is :math:`H = \frac{1}{4} R M^2 R^T`, with :math:`R` the rotation
    matrix and :math:`M^2 = \mathrm{diag}(-\Delta m^2, \Delta m^2)` the
    traceless mass matrix.  The multiplicative factor :math:`1/E` is
    *not* applied.

    Parameters
    ----------
    sth : float
        :math:`\sin\theta`.
    Dm2 : float
        Mass-squared difference :math:`\Delta m^2` [eV\ :sup:`2`].
    compute_matrix_multiplication : bool, optional
        If ``False`` (default), use the pre-computed closed-form
        expressions; if ``True``, carry out the matrix multiplication
        :math:`R M^2 R^T` explicitly.  Both give the same result; the
        option exists as a cross-check.

    Returns
    -------
    numpy.ndarray
        The :math:`2\times2` Hamiltonian [eV\ :sup:`2`], to be divided
        by the neutrino energy before use.

    Notes
    -----
    See the module-level *Sign convention* section: the mass eigenstate
    with the larger mass-squared value is the second one, so that adding
    a positive matter potential to the :math:`ee` entry describes
    neutrinos.

    Examples
    --------
    >>> H = hamiltonian_2nu_vacuum_energy_independent(0.5, 1.0)
    >>> print('%.6f  %.6f' % (H[0][0].real, H[0][1].real))
    -0.125000  0.216506
    """
    # Trigonometric identities, rather than arcsin followed by cos and
    # sin, keep this consistent with mixing_matrix_2nu and avoid a
    # needless round trip through the angle itself.
    cth = np.sqrt(1.0-sth*sth)
    c2th = 1.0-2.0*sth*sth
    s2th = 2.0*sth*cth

    f = 1./4.

    if not compute_matrix_multiplication:

        H00 = -Dm2*c2th
        H01 = Dm2*s2th
        H10 = H01
        H11 = -H00

        H = np.array([[H00*f, H01*f], [H10*f, H11*f]], dtype=complex)

    else:

        # Rotation matrix
        R = np.array(mixing_matrix_2nu(sth))
        # Traceless mass matrix
        M2 = np.array([[-Dm2, 0.0], [0.0, Dm2]])
        # Hamiltonian
        H = (f*(R @ M2 @ R.T)).astype(complex)

    return H


def probabilities_2nu_vacuum_std(sth, Dm2, energy, L):
    r"""Returns the 2nu vacuum probabilities, standard computation.

    Returns the probabilities for two-neutrino oscillations in vacuum,
    computed with the standard analytical expression

    .. math::
       P_{e\mu} = \sin^2 2\theta \sin^2\left(\frac{\Delta m^2 L}{4E}\right).

    This routine exists to validate the exact SU(2) computation in
    :mod:`oscprob2nu`; the two agree to round-off.

    Parameters
    ----------
    sth : float
        :math:`\sin\theta`.
    Dm2 : float
        Mass-squared difference :math:`\Delta m^2` [eV\ :sup:`2`].
    energy : float
        Neutrino energy [eV].
    L : float
        Baseline [eV\ :sup:`-1`].

    Returns
    -------
    list of float
        The probabilities ``[Pee, Pem, Pme, Pmm]``.

    See Also
    --------
    oscprob2nu.probabilities_2nu : The exact SU(2) computation.

    Examples
    --------
    >>> prob = probabilities_2nu_vacuum_std(0.5, 2.5e-3, 1.0e9, 5.0e12)
    >>> print('%.6f  %.6f' % (prob[0], prob[1]))
    0.999794  0.000206
    """
    arg = Dm2*L/4.0/energy
    cth = np.sqrt(1.0-sth*sth)
    s2th = 2.0*sth*cth

    Pem = s2th*s2th * pow(np.sin(arg), 2.0)
    Pme = Pem
    Pee = 1.0-Pem
    Pmm = 1.0-Pme

    return [Pee, Pem, Pme, Pmm]


def hamiltonian_2nu_matter(h_vacuum_energy_independent, energy, VCC):
    r"""Returns the two-neutrino Hamiltonian for matter oscillations.

    Computes and returns the :math:`2\times2` two-neutrino Hamiltonian
    for oscillations in matter of constant density, obtained by adding
    the charged-current matter potential to the :math:`ee` entry of the
    vacuum Hamiltonian.

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent part of the two-neutrino vacuum Hamiltonian
        [eV\ :sup:`2`], as returned by
        `hamiltonian_2nu_vacuum_energy_independent`.  It is not
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
        The :math:`2\times2` Hamiltonian [eV], of shape ``(2, 2)`` for a
        scalar energy and ``(..., 2, 2)`` for an array of energies.

    Examples
    --------
    >>> H_vac = hamiltonian_2nu_vacuum_energy_independent(0.5, 2.5e-3)
    >>> H = hamiltonian_2nu_matter(H_vac, 1.0e9, 1.0e-13)
    >>> print('%.6e' % H[0][0].real)
    -2.125000e-13
    """
    h_vacuum = np.asarray(h_vacuum_energy_independent, dtype=complex)
    energy = np.asarray(energy, dtype=float)
    VCC = np.asarray(VCC, dtype=float)

    # Indexing the energy with two trailing axes lets a scalar energy
    # return a single 2x2 matrix and an array of energies return one
    # matrix per energy, through the same expression
    return h_vacuum/energy[..., None, None] \
        + VCC[..., None, None]*_EE_PROJECTOR


def probabilities_2nu_matter_std(sth, Dm2, VCC, energy, L):
    r"""Returns the 2nu matter probabilities, standard computation.

    Returns the probabilities for two-neutrino oscillations in matter of
    constant density, computed with the standard analytical expression
    in terms of the effective mixing angle :math:`\theta_m` and
    effective mass-squared difference :math:`\Delta m^2_m`,

    .. math::
       \sin^2 2\theta_m = \frac{\sin^2 2\theta}
                               {\sin^2 2\theta + (\cos 2\theta - x)^2} ,
       \quad
       x \equiv \frac{2 V_{\rm CC} E}{\Delta m^2} .

    This routine exists to validate the exact SU(2) computation in
    :mod:`oscprob2nu`; the two agree to round-off.

    Parameters
    ----------
    sth : float
        :math:`\sin\theta`.
    Dm2 : float
        Mass-squared difference :math:`\Delta m^2` [eV\ :sup:`2`].
    VCC : float
        Potential due to charged-current interactions of
        :math:`\nu_e` with electrons [eV].
    energy : float
        Neutrino energy [eV].
    L : float
        Baseline [eV\ :sup:`-1`].

    Returns
    -------
    list of float
        The probabilities ``[Pee, Pem, Pme, Pmm]``.

    Notes
    -----
    The resonance sits at :math:`x = \cos 2\theta`, which for
    :math:`\theta < \pi/4` lies at positive energy, i.e. in the
    neutrino channel.  Note that :math:`\cos 2\theta` is *signed*:
    computing it as :math:`\sqrt{1 - \sin^2 2\theta}` would lose the
    sign and misplace the resonance for :math:`\theta > \pi/4`.

    See Also
    --------
    oscprob2nu.probabilities_2nu : The exact SU(2) computation.

    Examples
    --------
    >>> prob = probabilities_2nu_matter_std(0.5, 2.5e-3, 1.0e-13, 1.0e9,
    ...                                     5.0e12)
    >>> print('%.6f  %.6f' % (prob[0], prob[1]))
    0.985595  0.014405
    """
    x = 2.0*VCC*energy/Dm2
    cth = np.sqrt(1.0-sth*sth)
    s2th = 2.0*sth*cth
    s2thsq = s2th*s2th
    # cos(2*theta) is signed; sqrt(1 - sin^2(2*theta)) would not be
    c2th = 1.0-2.0*sth*sth

    denominator = s2thsq+pow(c2th-x, 2.0)

    Dm2m = Dm2*np.sqrt(denominator)
    s2thmsq = s2thsq / denominator

    arg = Dm2m*L/4.0/energy

    Pem = s2thmsq * pow(np.sin(arg), 2.0)
    Pme = Pem
    Pee = 1.0-Pem
    Pmm = 1.0-Pme

    return [Pee, Pem, Pme, Pmm]


def hamiltonian_2nu_nsi(h_vacuum_energy_independent, energy, VCC, eps):
    r"""Returns the two-neutrino Hamiltonian for oscillations with NSI.

    Computes and returns the :math:`2\times2` two-neutrino Hamiltonian
    for oscillations with non-standard interactions (NSI) in matter of
    constant density.

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent part of the two-neutrino vacuum Hamiltonian
        [eV\ :sup:`2`], as returned by
        `hamiltonian_2nu_vacuum_energy_independent`.  It is not
        modified.
    energy : float or array_like
        Neutrino energy [eV], or an array of energies, in which case one
        Hamiltonian is returned per energy.
    VCC : float
        Potential due to charged-current interactions of
        :math:`\nu_e` with electrons [eV].
    eps : array_like
        The three NSI strength parameters
        ``[eps_ee, eps_em, eps_mm]``, adimensional.  The diagonal
        parameters ``eps_ee`` and ``eps_mm`` are real; the off-diagonal
        ``eps_em`` may be complex, and its complex conjugate is placed
        in the lower off-diagonal entry so that the Hamiltonian stays
        Hermitian.

    Returns
    -------
    numpy.ndarray
        The :math:`2\times2` complex Hamiltonian [eV], of shape
        ``(2, 2)`` for a scalar energy and ``(..., 2, 2)`` for an
        array of energies.

    Examples
    --------
    >>> H_vac = hamiltonian_2nu_vacuum_energy_independent(0.5, 2.5e-3)
    >>> H = hamiltonian_2nu_nsi(H_vac, 1.0e9, 1.0e-13, [0.06, -0.06+0.03j,
    ...                                                 1.2])
    >>> print('%+.6e%+.6ej' % (H[0][1].real, H[0][1].imag))
    +5.352659e-13+3.000000e-15j
    """
    h_vacuum = np.asarray(h_vacuum_energy_independent, dtype=complex)
    energy = np.asarray(energy, dtype=float)
    VCC = np.asarray(VCC, dtype=float)

    eps_ee, eps_em, eps_mm = eps

    # The matrix is complex so that a complex eps_em keeps its imaginary
    # part; a real one would silently discard it
    nsi = np.array([[1.0+eps_ee, eps_em],
                    [np.conj(eps_em), eps_mm]], dtype=complex)

    return h_vacuum/energy[..., None, None] + VCC[..., None, None]*nsi


def hamiltonian_2nu_liv(h_vacuum_energy_independent, energy, sxi, b1, b2,
                        Lambda):
    r"""Returns the two-neutrino Hamiltonian for oscillations with LIV.

    Computes and returns the :math:`2\times2` two-neutrino Hamiltonian
    for oscillations in a CPT-odd Lorentz invariance-violating (LIV)
    background.  The LIV term is :math:`(E/\Lambda) R B_2 R^T`, with
    :math:`B_2 = \mathrm{diag}(b_1, b_2)` and :math:`R` the rotation by
    the angle :math:`\xi` between the eigenvectors of :math:`B_2` and
    the flavor states.

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent part of the two-neutrino vacuum Hamiltonian
        [eV\ :sup:`2`], as returned by
        `hamiltonian_2nu_vacuum_energy_independent`.  It is not
        modified.
    energy : float or array_like
        Neutrino energy [eV], or an array of energies, in which case one
        Hamiltonian is returned per energy.
    sxi : float
        :math:`\sin\xi`, with :math:`\xi` the rotation angle between the
        space of the eigenvectors of :math:`B_2` and the flavor states.
    b1 : float
        Eigenvalue :math:`b_1` of the LIV operator :math:`B_2` [eV].
    b2 : float
        Eigenvalue :math:`b_2` of the LIV operator :math:`B_2` [eV].
    Lambda : float
        Energy scale :math:`\Lambda` of the LIV operator :math:`B_2`
        [eV].

    Returns
    -------
    numpy.ndarray
        The :math:`2\times2` complex Hamiltonian [eV], of shape
        ``(2, 2)`` for a scalar energy and ``(..., 2, 2)`` for an
        array of energies.

    Examples
    --------
    >>> H_vac = hamiltonian_2nu_vacuum_energy_independent(0.5, 2.5e-3)
    >>> H = hamiltonian_2nu_liv(H_vac, 1.0e9, 0.6, 1.0e-9, 2.0e-9, 1.0e12)
    >>> print('%.6e' % H[0][0].real)
    1.047500e-12
    """
    h_vacuum = np.asarray(h_vacuum_energy_independent, dtype=complex)
    energy = np.asarray(energy, dtype=float)

    cxi = np.sqrt(1.0-sxi*sxi)

    # R.B2.R^T, with R the rotation by xi and B2 = diag(b1, b2)
    liv = np.array([[b1*cxi*cxi + b2*sxi*sxi, (-b1+b2)*cxi*sxi],
                    [(-b1+b2)*cxi*sxi, b2*cxi*cxi + b1*sxi*sxi]],
                   dtype=complex)

    f = energy/Lambda

    return h_vacuum/energy[..., None, None] + f[..., None, None]*liv
