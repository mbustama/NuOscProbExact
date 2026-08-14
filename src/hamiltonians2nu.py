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
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

__all__ = ['mixing_matrix_2nu', 'hamiltonian_2nu_vacuum_energy_independent',
           'probabilities_2nu_vacuum_std', 'hamiltonian_2nu_matter',
           'probabilities_2nu_matter_std', 'hamiltonian_2nu_nsi',
           'hamiltonian_2nu_liv']

from typing import List, Union

import math

import numpy as np


_EE_PROJECTOR = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=complex)
r"""numpy.ndarray: Module-level constant.

The matrix that the charged-current matter potential multiplies: it
selects the :math:`ee` entry, since only :math:`\nu_e` interacts through
charged currents with the electrons in matter.
"""


ANGLE_CONVENTIONS = ('sin', 'sin2', 'rad', 'deg')
r"""tuple of str: Module-level constant.

The ways a mixing angle may be given to the routines here, as the
``angles`` argument they all take.  ``'sin'`` is the default and is what
every earlier version accepted.

The other three exist because the published numbers are not sines.
Global fits quote :math:`\sin^2\theta_{ij}` --- NuFit's :math:`0.310`,
:math:`0.582`, :math:`0.02240` --- whose square roots are what
``'sin'`` wants, and passing the published value under the default is a
silent error rather than a loud one: :math:`0.310` is a perfectly legal
sine, of a different angle.  ``'sin2'`` lets the published number be
typed as printed.  ``'rad'`` and ``'deg'`` take the angle itself.

A CP-violating phase is not a mixing angle and has no sine to pass, so
``'sin'`` and ``'sin2'`` leave it in radians; under ``'rad'`` and
``'deg'`` it follows the angles, which lets a whole parameter set be
given in the units it was published in.
"""


def _sine_from(
    value: Union[int, float],
    angles: str,
    name: str,
    caller: str
) -> float:
    r"""Returns :math:`\sin\theta` from an angle in any convention.

    Parameters
    ----------
    value : int or float
        The angle, expressed as `angles` says.
    angles : str
        One of `ANGLE_CONVENTIONS`.
    name : str
        Name of the parameter, used in the error message.
    caller : str
        Name of the calling routine, used in the error message.

    Returns
    -------
    float
        :math:`\sin\theta`.

    Raises
    ------
    ValueError
        If `angles` is not one of `ANGLE_CONVENTIONS`, or if `value` is
        outside the range that convention allows.
    """
    # First, and free, so that the default costs one comparison
    if angles == 'sin':
        return value

    if angles == 'sin2':
        if not 0.0 <= value <= 1.0:
            raise ValueError(
                "%s: with angles='sin2', %s is the square of a sine and so "
                'lies in [0, 1]; got %r' % (caller, name, value))
        # Non-negative by construction, which loses the sign a mixing
        # angle outside the first octant would carry.  Fits quote the
        # square precisely because that sign is conventional.
        return math.sqrt(value)

    if angles == 'rad':
        return math.sin(value)

    if angles == 'deg':
        return math.sin(math.radians(value))

    raise ValueError(
        '%s: angles must be one of %s; got %r'
        % (caller, ', '.join(repr(c) for c in ANGLE_CONVENTIONS), angles))


def _phase_from(
    value: Union[int, float],
    angles: str
) -> float:
    r"""Returns a CP-violating phase in radians, in any convention.

    A phase has no sine to pass --- :math:`\sin\delta` would lose the
    quadrant --- so it is unchanged under ``'sin'`` and ``'sin2'``, and
    follows the angles under ``'rad'`` and ``'deg'``.

    Parameters
    ----------
    value : int or float
        The phase, in radians unless `angles` is ``'deg'``.
    angles : str
        One of `ANGLE_CONVENTIONS`.

    Returns
    -------
    float
        The phase, in radians.
    """
    return math.radians(value) if angles == 'deg' else value


def _cos_from_sin(sth: Union[int, float], name: str, caller: str) -> float:
    r"""Returns :math:`\cos\theta` from :math:`\sin\theta`, checked.

    The mixing parameters throughout **NuOscProbExact** are sines of the
    angles, so every one of them has to lie in :math:`[-1, 1]`.  Taking
    the cosine as :math:`\sqrt{1 - \sin^2\theta}` without checking
    turns a value outside that range into whatever the square root does
    with a negative argument, which differed between the flavor counts:
    :mod:`math` raised ``math domain error``, naming neither the
    parameter nor the value, while :obj:`numpy.sqrt` returned ``nan``
    and let it propagate silently into the probabilities.

    Parameters
    ----------
    sth : int or float
        Sine of the angle.
    name : str
        Name of the parameter, used in the error message.
    caller : str
        Name of the calling routine, used in the error message.

    Returns
    -------
    float
        :math:`\cos\theta`, taken non-negative.

    Raises
    ------
    ValueError
        If ``sth`` does not lie in :math:`[-1, 1]`, or is not a number.

    .. versionadded:: 1.11.0
    """
    if not -1.0 <= sth <= 1.0:
        raise ValueError(
            '%s: %s must be the sine of an angle and so lie in [-1, 1]; '
            'got %r.  The mixing parameters are sines, not angles.'
            % (caller, name, sth))

    return math.sqrt(1.0 - sth*sth)


def mixing_matrix_2nu(sth: Union[int, float], angles: str = 'sin') -> List[List[float]]:
    r"""Returns the :math:`2\times2` rotation matrix.

    Computes and returns the real :math:`2\times2` rotation matrix
    parametrized by a single rotation angle :math:`\theta`.

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.1.0
       :math:`\cos\theta` is taken as
       :math:`\sqrt{1 - \sin^2\theta}` rather than through the angle
       itself.  The matrix is real, and is still returned as a nested
       list, as it always has been.

    .. versionchanged:: 1.4.0
       Faster, with identical results --- all 42 figures generated by
       ``run_testsuite.py`` are byte-for-byte those of 1.3.0. The scalar
       path stopped dispatching NumPy for single numbers:
       :func:`numpy.real`, :func:`numpy.imag`, :obj:`numpy.arccos`,
       :func:`numpy.clip` and :obj:`numpy.sqrt` on one number give way
       to attribute access and the :mod:`math` module.

    Parameters
    ----------
    sth : float
        :math:`\sin\theta`, with :math:`\theta` in the first quadrant,
        so that :math:`\cos\theta = \sqrt{1-\sin^2\theta} \geq 0`.

    angles : str, optional
        How the mixing angles are given: ``'sin'`` for their sines, the
        default and what earlier versions accepted; ``'sin2'`` for the
        squares of their sines, which is how global fits publish them;
        ``'rad'`` or ``'deg'`` for the angles themselves.  See
        `ANGLE_CONVENTIONS`.
    Returns
    -------
    list of list of float
        The rotation matrix ``[[cth, sth], [-sth, cth]]``, with
        ``cth`` = :math:`\cos\theta`.

    Examples
    --------
    .. jupyter-execute::

        import hamiltonians2nu

        R = hamiltonians2nu.mixing_matrix_2nu(0.6)
        print('%.6f  %.6f' % (R[0][0], R[0][1]))
    """
    sth = _sine_from(sth, angles, 'sth', 'mixing_matrix_2nu')
    cth = _cos_from_sin(sth, 'sth', 'mixing_matrix_2nu')

    return [[cth, sth], [-sth, cth]]


def hamiltonian_2nu_vacuum_energy_independent(
    sth: Union[int, float],
    Dm2: Union[int, float],
    compute_matrix_multiplication: bool = False,
    angles: str = 'sin'
) -> np.ndarray:
    r"""Returns the two-neutrino Hamiltonian for vacuum oscillations.

    Computes and returns the energy-independent part of the real
    :math:`2\times2` two-neutrino Hamiltonian for oscillations in
    vacuum, parametrized by a single mixing angle :math:`\theta` and a
    single mass-squared difference :math:`\Delta m^2`.  The Hamiltonian
    is :math:`H = \frac{1}{4} R M^2 R^T`, with :math:`R` the rotation
    matrix and :math:`M^2 = \mathrm{diag}(-\Delta m^2, \Delta m^2)` the
    traceless mass matrix.  The multiplicative factor :math:`1/E` is
    *not* applied.

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.1.0
       The sign convention was corrected.  The Hamiltonian was built
       from :math:`M^2 = \mathrm{diag}(\Delta m^2, -\Delta m^2)`, which
       yields the negative of the textbook Hamiltonian.  In vacuum this
       is invisible, but it reverses the sign of the matter potential
       relative to the vacuum term, so results in matter, with NSI, or
       with LIV were the antineutrino ones.  It also returns a complex
       :class:`numpy.ndarray` rather than a nested list.

    .. versionchanged:: 1.4.0
       Faster, with identical results --- all 42 figures generated by
       ``run_testsuite.py`` are byte-for-byte those of 1.3.0. The scalar
       path stopped dispatching NumPy for single numbers:
       :func:`numpy.real`, :func:`numpy.imag`, :obj:`numpy.arccos`,
       :func:`numpy.clip` and :obj:`numpy.sqrt` on one number give way
       to attribute access and the :mod:`math` module.

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

    angles : str, optional
        How the mixing angles are given: ``'sin'`` for their sines, the
        default and what earlier versions accepted; ``'sin2'`` for the
        squares of their sines, which is how global fits publish them;
        ``'rad'`` or ``'deg'`` for the angles themselves.  See
        `ANGLE_CONVENTIONS`.
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
    .. jupyter-execute::

        import hamiltonians2nu

        H = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(0.5, 1.0)
        print('%.6f  %.6f' % (H[0][0].real, H[0][1].real))
    """
    sth = _sine_from(sth, angles, 'sth', 'hamiltonian_2nu_vacuum_energy_independent')
    # Trigonometric identities, rather than arcsin followed by cos and
    # sin, keep this consistent with mixing_matrix_2nu and avoid a
    # needless round trip through the angle itself.
    cth = _cos_from_sin(sth, 'sth',
                        'hamiltonian_2nu_vacuum_energy_independent')
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


def probabilities_2nu_vacuum_std(
    sth: Union[int, float],
    Dm2: Union[int, float],
    energy: Union[int, float],
    L: Union[int, float],
    angles: str = 'sin'
) -> List[float]:
    r"""Returns the 2nu vacuum probabilities, standard computation.

    Returns the probabilities for two-neutrino oscillations in vacuum,
    computed with the standard analytical expression

    .. math::
       P_{e\mu} = \sin^2 2\theta \sin^2\left(\frac{\Delta m^2 L}{4E}\right).

    This routine exists to validate the exact SU(2) computation in
    :mod:`oscprob2nu`; the two agree to round-off.

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.1.0
       The signature changed: the energy is now given in eV and the
       baseline in :math:`\mathrm{eV}^{-1}`, like the rest of the
       library, rather than in GeV and km.  The rounded constants 1.27
       and 2.54 that folded in the old conversion overstated every phase
       by 0.242%.

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

    angles : str, optional
        How the mixing angles are given: ``'sin'`` for their sines, the
        default and what earlier versions accepted; ``'sin2'`` for the
        squares of their sines, which is how global fits publish them;
        ``'rad'`` or ``'deg'`` for the angles themselves.  See
        `ANGLE_CONVENTIONS`.
    Returns
    -------
    list of float
        The probabilities ``[Pee, Pem, Pme, Pmm]``.

    See Also
    --------
    oscprob2nu.probabilities_2nu : The exact SU(2) computation.

    Examples
    --------
    .. jupyter-execute::

        import hamiltonians2nu

        prob = hamiltonians2nu.probabilities_2nu_vacuum_std(0.5, 2.5e-3, 1.0e9, 5.0e12)
        print('%.6f  %.6f' % (prob[0], prob[1]))
    """
    sth = _sine_from(sth, angles, 'sth', 'probabilities_2nu_vacuum_std')
    arg = Dm2*L/4.0/energy
    cth = np.sqrt(1.0-sth*sth)
    s2th = 2.0*sth*cth

    Pem = s2th*s2th * pow(np.sin(arg), 2.0)
    Pme = Pem
    Pee = 1.0-Pem
    Pmm = 1.0-Pme

    return [Pee, Pem, Pme, Pmm]


def hamiltonian_2nu_matter(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float, list, np.ndarray],
    VCC: Union[int, float, list, np.ndarray]
) -> np.ndarray:
    r"""Returns the two-neutrino Hamiltonian for matter oscillations.

    Computes and returns the :math:`2\times2` two-neutrino Hamiltonian
    for oscillations in matter of constant density, obtained by adding
    the charged-current matter potential to the :math:`ee` entry of the
    vacuum Hamiltonian.

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.1.0
       Results changed, following the sign-convention correction in
       `hamiltonian_2nu_vacuum_energy_independent`: this routine
       previously returned the antineutrino Hamiltonian when asked for
       the neutrino one, placing the MSW resonance on the wrong side.
       Returns a complex :class:`numpy.ndarray`.

    .. versionchanged:: 1.3.0
       Accepts an array of energies, returning one Hamiltonian per
       energy stacked along a leading axis; the matter potential may be
       an array too.  A scalar energy still returns a single matrix, and
       the results are bit-for-bit what the equivalent loop produced.

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
    VCC : float or array_like
        Potential due to charged-current interactions of
        :math:`\nu_e` with electrons [eV].  Positive for neutrinos,
        negative for antineutrinos.  May be an array, to scan across a
        density profile alongside the energy.

    Returns
    -------
    numpy.ndarray
        The :math:`2\times2` Hamiltonian [eV], of shape ``(2, 2)`` for a
        scalar energy and ``(..., 2, 2)`` for an array of energies.

    Examples
    --------
    .. jupyter-execute::

        import hamiltonians2nu

        H_vac = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(0.5, 2.5e-3)
        H = hamiltonians2nu.hamiltonian_2nu_matter(H_vac, 1.0e9, 1.0e-13)
        print('%.6e' % H[0][0].real)
    """
    h_vacuum = np.asarray(h_vacuum_energy_independent, dtype=complex)
    energy = np.asarray(energy, dtype=float)
    VCC = np.asarray(VCC, dtype=float)

    # Indexing the energy with two trailing axes lets a scalar energy
    # return a single 2x2 matrix and an array of energies return one
    # matrix per energy, through the same expression
    return h_vacuum/energy[..., None, None] \
        + VCC[..., None, None]*_EE_PROJECTOR


def probabilities_2nu_matter_std(
    sth: Union[int, float],
    Dm2: Union[int, float],
    VCC: Union[int, float],
    energy: Union[int, float],
    L: Union[int, float],
    angles: str = 'sin'
) -> List[float]:
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

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.1.0
       The signature changed: the energy is now given in eV and the
       baseline in :math:`\mathrm{eV}^{-1}`, like the rest of the
       library, rather than in GeV and km.  The rounded constants 1.27
       and 2.54 that folded in the old conversion overstated every phase
       by 0.242%.  The sign of :math:`\cos 2\theta` is also kept, where
       it was previously discarded by computing it as :math:`\sqrt{1 -
       \sin^2 2\theta}`; for :math:`\theta > \pi/4` that put the matter
       resonance on the wrong side.

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

    angles : str, optional
        How the mixing angles are given: ``'sin'`` for their sines, the
        default and what earlier versions accepted; ``'sin2'`` for the
        squares of their sines, which is how global fits publish them;
        ``'rad'`` or ``'deg'`` for the angles themselves.  See
        `ANGLE_CONVENTIONS`.
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
    .. jupyter-execute::

        import hamiltonians2nu

        prob = hamiltonians2nu.probabilities_2nu_matter_std(0.5, 2.5e-3, 1.0e-13, 1.0e9,
                                            5.0e12)
        print('%.6f  %.6f' % (prob[0], prob[1]))
    """
    sth = _sine_from(sth, angles, 'sth', 'probabilities_2nu_matter_std')
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


def hamiltonian_2nu_nsi(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float, list, np.ndarray],
    VCC: Union[int, float, list, np.ndarray],
    eps: Union[list, np.ndarray]
) -> np.ndarray:
    r"""Returns the two-neutrino Hamiltonian for oscillations with NSI.

    Computes and returns the :math:`2\times2` two-neutrino Hamiltonian
    for oscillations with non-standard interactions (NSI) in matter of
    constant density.

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.1.0
       The imaginary part of a complex :math:`\epsilon_{e\mu}` is no
       longer discarded: the vacuum Hamiltonian was real, so the array
       was ``float64`` and the in-place addition truncated the value.
       Results also changed with the sign-convention correction
       described under `hamiltonian_2nu_matter`.

    .. versionchanged:: 1.3.0
       Accepts an array of energies, returning one Hamiltonian per
       energy stacked along a leading axis; the matter potential may be
       an array too.  A scalar energy still returns a single matrix, and
       the results are bit-for-bit what the equivalent loop produced.

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
    VCC : float or array_like
        Potential due to charged-current interactions of
        :math:`\nu_e` with electrons [eV].  May be an array, to scan
        across a density profile alongside the energy.
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
    .. jupyter-execute::

        import hamiltonians2nu

        H_vac = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(0.5, 2.5e-3)
        H = hamiltonians2nu.hamiltonian_2nu_nsi(H_vac, 1.0e9, 1.0e-13, [0.06, -0.06+0.03j,
                                                        1.2])
        print('%+.6e%+.6ej' % (H[0][1].real, H[0][1].imag))
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


def hamiltonian_2nu_liv(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float, list, np.ndarray],
    sxi: Union[int, float],
    b1: Union[int, float],
    b2: Union[int, float],
    Lambda: Union[int, float]
) -> np.ndarray:
    r"""Returns the two-neutrino Hamiltonian for oscillations with LIV.

    Computes and returns the :math:`2\times2` two-neutrino Hamiltonian
    for oscillations in a CPT-odd Lorentz invariance-violating (LIV)
    background.  The LIV term is :math:`(E/\Lambda) R B_2 R^T`, with
    :math:`B_2 = \mathrm{diag}(b_1, b_2)` and :math:`R` the rotation by
    the angle :math:`\xi` between the eigenvectors of :math:`B_2` and
    the flavor states.

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.1.0
       :math:`\cos\xi` was computed as ``sqrt(1 - sxi - sxi)`` rather
       than ``sqrt(1 - sxi*sxi)``.  For :math:`0 < \sin\xi < 1/2` the
       LIV term was not a rotation at all, and for :math:`\sin\xi \geq
       1/2` the whole Hamiltonian became NaN.  Results also changed with
       the sign-convention correction described under
       `hamiltonian_2nu_matter`.

    .. versionchanged:: 1.3.0
       Accepts an array of energies, returning one Hamiltonian per
       energy stacked along a leading axis.  The LIV term scales with
       the energy rather than being added at constant strength, so it is
       formed per entry.  A scalar energy still returns a single matrix,
       and the results are bit-for-bit what the equivalent loop
       produced.

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
    .. jupyter-execute::

        import hamiltonians2nu

        H_vac = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(0.5, 2.5e-3)
        H = hamiltonians2nu.hamiltonian_2nu_liv(H_vac, 1.0e9, 0.6, 1.0e-9, 2.0e-9, 1.0e12)
        print('%.6e' % H[0][0].real)
    """
    h_vacuum = np.asarray(h_vacuum_energy_independent, dtype=complex)
    energy = np.asarray(energy, dtype=float)

    cxi = _cos_from_sin(sxi, 'sxi', 'hamiltonian_2nu_liv')

    # R.B2.R^T, with R the rotation by xi and B2 = diag(b1, b2)
    liv = np.array([[b1*cxi*cxi + b2*sxi*sxi, (-b1+b2)*cxi*sxi],
                    [(-b1+b2)*cxi*sxi, b2*cxi*cxi + b1*sxi*sxi]],
                   dtype=complex)

    f = energy/Lambda

    return h_vacuum/energy[..., None, None] + f[..., None, None]*liv
