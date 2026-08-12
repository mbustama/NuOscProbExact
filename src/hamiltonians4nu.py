# -*- coding: utf-8 -*-
r"""Sample four-neutrino Hamiltonians, for 3+1 scenarios.

This module builds the :math:`4\times4` Hamiltonians that
:mod:`oscprob4nu` evaluates, in the same spirit as
:mod:`hamiltonians3nu` does at three flavors and as [1]_ describes: they
are *examples*, not limitations.  :func:`oscprob4nu.probabilities_4nu` takes any Hermitian
:math:`4\times4` matrix, so a scenario not built here is a matrix away.

The fourth state is written as sterile throughout, so the flavor basis
is :math:`(\nu_e, \nu_\mu, \nu_\tau, \nu_s)`, and the mixing matrix
carries three extra angles :math:`\theta_{14}, \theta_{24},
\theta_{34}` on top of the three standard ones.  Nothing in
:mod:`oscprob4nu` depends on that reading: a fourth *active* state, or
any other four-level system, works identically.

Why 3+1 is in scope here
------------------------

A 3+1 system is often described as "leaky" from the three-flavor point
of view, because probability disappears from the active block into the
sterile state.  That is a statement about the :math:`3\times3`
subsystem, not about the physics: over all four states the evolution is
closed and unitary, which is exactly the assumption the SU(4) expansion
needs.  Treating it at :math:`n = 4` therefore brings it back inside the
scope of an exact closed-form method.

The matter potential
--------------------

Active neutrinos feel the charged-current potential :math:`V_{CC}` (the
electron flavor only) and the flavor-universal neutral-current potential
:math:`V_{NC}` (all three); a sterile state feels neither.  A term
proportional to the identity contributes only a global phase, so
subtracting :math:`V_{NC}\mathbb{1}` from all four states costs nothing
and leaves

.. math::
   A_4 = \mathrm{diag}\left(V_{CC},\, 0,\, 0,\, -V_{NC}\right) ,

with :math:`-V_{NC} = +G_F n_n/\sqrt{2}` positive.  The sterile entry is
therefore *not* zero, and getting it wrong is the four-flavor analogue
of the antineutrino sign trap: the difference is invisible in vacuum and
sets the position of the matter resonance.

Routine listings
----------------

    * mixing_matrix_4nu - Returns the 3+1 mixing matrix
    * hamiltonian_4nu_vacuum_energy_independent - Vacuum Hamiltonian
    * hamiltonian_4nu_matter - Adds matter of constant density
    * hamiltonian_4nu_nsi - Adds non-standard interactions
    * hamiltonian_4nu_liv - Adds a Lorentz invariance-violating term

References
----------

.. [1] Mauricio Bustamante, "Exact neutrino oscillation probabilities
   with arbitrary time-independent Hamiltonians", arXiv:1904.12391.
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

__all__ = ['mixing_matrix_4nu',
           'hamiltonian_4nu_vacuum_energy_independent',
           'hamiltonian_4nu_matter', 'hamiltonian_4nu_nsi',
           'hamiltonian_4nu_liv']

from typing import Union

import math

import numpy as np


_EE_PROJECTOR_4NU = np.diag([1.0, 0.0, 0.0, 0.0])
_SS_PROJECTOR_4NU = np.diag([0.0, 0.0, 0.0, 1.0])


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


def _rotation_4nu(
    index_1: int,
    index_2: int,
    sth: Union[int, float],
    delta: Union[int, float] = 0.0
) -> np.ndarray:
    r"""Returns a complex rotation in the ``(index_1, index_2)`` plane.

    Parameters
    ----------
    index_1 : int
        First index of the rotation plane.
    index_2 : int
        Second index of the rotation plane.
    sth : int or float
        Sine of the rotation angle --- *not* the angle itself.
    delta : int or float, optional
        CP-violating phase attached to this rotation.  Default: 0.

    Returns
    -------
    numpy.ndarray
        Complex array of shape ``(4, 4)``.
    """
    cth = _cos_from_sin(sth, 'sth', 'mixing_matrix_4nu')
    rotation = np.eye(4, dtype=complex)
    rotation[index_1, index_1] = cth
    rotation[index_2, index_2] = cth
    rotation[index_1, index_2] = sth*np.exp(-1.j*delta)
    rotation[index_2, index_1] = -sth*np.exp(1.j*delta)

    return rotation


def mixing_matrix_4nu(
    s12: Union[int, float],
    s23: Union[int, float],
    s13: Union[int, float],
    s14: Union[int, float],
    s24: Union[int, float],
    s34: Union[int, float],
    dCP: Union[int, float],
    d14: Union[int, float] = 0.0,
    d24: Union[int, float] = 0.0,
    angles: str = 'sin'
) -> np.ndarray:
    r"""Returns the 3+1 lepton mixing matrix.

    Built in the common 3+1 ordering

    .. math::
       U = R_{34} R_{24}(\delta_{24}) R_{14}(\delta_{14})
           R_{23} R_{13}(\delta_{CP}) R_{12} ,

    which reduces to the standard PDG three-flavor matrix of
    :func:`hamiltonians3nu.pmns_mixing_matrix` in the upper-left block
    when the three new angles vanish.

    All six mixing parameters are **sines of the angles**, not the
    angles, matching the convention used throughout
    **NuOscProbExact**.

    .. versionadded:: 1.9.0

    Parameters
    ----------
    s12 : int or float
        Sine of :math:`\theta_{12}`.
    s23 : int or float
        Sine of :math:`\theta_{23}`.
    s13 : int or float
        Sine of :math:`\theta_{13}`.
    s14 : int or float
        Sine of :math:`\theta_{14}`.
    s24 : int or float
        Sine of :math:`\theta_{24}`.
    s34 : int or float
        Sine of :math:`\theta_{34}`.
    dCP : int or float
        Standard Dirac CP-violating phase, in radian.
    d14 : int or float, optional
        Extra phase on the 1-4 rotation, in radian.  Default: 0.
    d24 : int or float, optional
        Extra phase on the 2-4 rotation, in radian.  Default: 0.

    angles : str, optional
        How the mixing angles are given: ``'sin'`` for their sines, the
        default and what earlier versions accepted; ``'sin2'`` for the
        squares of their sines, which is how global fits publish them;
        ``'rad'`` or ``'deg'`` for the angles themselves.  Under ``'rad'`` and
        ``'deg'`` the phases follow; under the other two they stay in
        radians, a phase having no sine to pass.  See
        `ANGLE_CONVENTIONS`.
    Returns
    -------
    numpy.ndarray
        Complex array of shape ``(4, 4)``.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import hamiltonians4nu

        mixing = hamiltonians4nu.mixing_matrix_4nu(
            np.sqrt(0.310), np.sqrt(0.582), np.sqrt(2.240e-2),
            np.sqrt(0.10), np.sqrt(0.10), 0.0, 217.0/180.0*np.pi)

        unitarity = mixing.conj().T @ mixing
        print('unitary to %.1e' % np.max(np.abs(unitarity - np.eye(4))))
    """
    s12 = _sine_from(s12, angles, 's12', 'mixing_matrix_4nu')
    s23 = _sine_from(s23, angles, 's23', 'mixing_matrix_4nu')
    s13 = _sine_from(s13, angles, 's13', 'mixing_matrix_4nu')
    s14 = _sine_from(s14, angles, 's14', 'mixing_matrix_4nu')
    s24 = _sine_from(s24, angles, 's24', 'mixing_matrix_4nu')
    s34 = _sine_from(s34, angles, 's34', 'mixing_matrix_4nu')
    dCP = _phase_from(dCP, angles)
    d14 = _phase_from(d14, angles)
    d24 = _phase_from(d24, angles)
    return (_rotation_4nu(2, 3, s34)
            @ _rotation_4nu(1, 3, s24, d24)
            @ _rotation_4nu(0, 3, s14, d14)
            @ _rotation_4nu(1, 2, s23)
            @ _rotation_4nu(0, 2, s13, dCP)
            @ _rotation_4nu(0, 1, s12))


def hamiltonian_4nu_vacuum_energy_independent(
    s12: Union[int, float],
    s23: Union[int, float],
    s13: Union[int, float],
    s14: Union[int, float],
    s24: Union[int, float],
    s34: Union[int, float],
    dCP: Union[int, float],
    D21: Union[int, float],
    D31: Union[int, float],
    D41: Union[int, float],
    d14: Union[int, float] = 0.0,
    d24: Union[int, float] = 0.0,
    angles: str = 'sin'
) -> np.ndarray:
    r"""Returns the energy-independent four-neutrino vacuum Hamiltonian.

    Returns :math:`U M^2 U^\dagger / 2`, with
    :math:`M^2 = \mathrm{diag}(0, \Delta m^2_{21}, \Delta m^2_{31},
    \Delta m^2_{41})`, so that the vacuum Hamiltonian at energy
    :math:`E` is this matrix divided by :math:`E`.

    The factor :math:`1/E` is deliberately left out, so that the
    energy-independent part can be computed once and reused across an
    energy scan --- the same arrangement as
    :func:`hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent`.

    .. versionadded:: 1.9.0

    Parameters
    ----------
    s12 : int or float
        Sine of :math:`\theta_{12}`.
    s23 : int or float
        Sine of :math:`\theta_{23}`.
    s13 : int or float
        Sine of :math:`\theta_{13}`.
    s14 : int or float
        Sine of :math:`\theta_{14}`.
    s24 : int or float
        Sine of :math:`\theta_{24}`.
    s34 : int or float
        Sine of :math:`\theta_{34}`.
    dCP : int or float
        Standard Dirac CP-violating phase, in radian.
    D21 : int or float
        Mass-squared difference :math:`\Delta m^2_{21}`, in eV\ :sup:`2`.
    D31 : int or float
        Mass-squared difference :math:`\Delta m^2_{31}`, in eV\ :sup:`2`.
    D41 : int or float
        Mass-squared difference :math:`\Delta m^2_{41}`, in eV\ :sup:`2`.
    d14 : int or float, optional
        Extra phase on the 1-4 rotation, in radian.  Default: 0.
    d24 : int or float, optional
        Extra phase on the 2-4 rotation, in radian.  Default: 0.

    angles : str, optional
        How the mixing angles are given: ``'sin'`` for their sines, the
        default and what earlier versions accepted; ``'sin2'`` for the
        squares of their sines, which is how global fits publish them;
        ``'rad'`` or ``'deg'`` for the angles themselves.  Under ``'rad'`` and
        ``'deg'`` the phases follow; under the other two they stay in
        radians, a phase having no sine to pass.  See
        `ANGLE_CONVENTIONS`.
    Returns
    -------
    numpy.ndarray
        Complex array of shape ``(4, 4)``, in units of
        eV\ :sup:`2`.  Divide by the energy in eV to obtain a
        Hamiltonian in eV.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import globaldefs as gd
        import hamiltonians4nu
        import oscprob4nu

        h_vacuum = hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(
            gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF,
            np.sqrt(0.10), np.sqrt(0.10), 0.0,
            gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF, 1.0)

        prob = oscprob4nu.probabilities_4nu(
            h_vacuum/1.0e9, 1300.0*gd.CONV_KM_TO_INV_EV)

        print('P_mumu = %.6f' % prob[5])
        print('P_mus  = %.6f' % prob[7])
    """
    s12 = _sine_from(s12, angles, 's12', 'hamiltonian_4nu_vacuum_energy_independent')
    s23 = _sine_from(s23, angles, 's23', 'hamiltonian_4nu_vacuum_energy_independent')
    s13 = _sine_from(s13, angles, 's13', 'hamiltonian_4nu_vacuum_energy_independent')
    s14 = _sine_from(s14, angles, 's14', 'hamiltonian_4nu_vacuum_energy_independent')
    s24 = _sine_from(s24, angles, 's24', 'hamiltonian_4nu_vacuum_energy_independent')
    s34 = _sine_from(s34, angles, 's34', 'hamiltonian_4nu_vacuum_energy_independent')
    dCP = _phase_from(dCP, angles)
    d14 = _phase_from(d14, angles)
    d24 = _phase_from(d24, angles)
    mixing = mixing_matrix_4nu(s12, s23, s13, s14, s24, s34, dCP, d14, d24)
    masses = np.diag([0.0, D21, D31, D41]).astype(complex)

    return mixing @ masses @ mixing.conj().T / 2.0


def hamiltonian_4nu_matter(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float, list, np.ndarray],
    VCC: Union[int, float, list, np.ndarray],
    VNC: Union[int, float, list, np.ndarray]
) -> np.ndarray:
    r"""Returns the four-neutrino Hamiltonian in matter.

    Adds to the vacuum term the matter potential

    .. math::
       A_4 = \mathrm{diag}\left(V_{CC},\, 0,\, 0,\, -V_{NC}\right) ,

    which is what remains after the flavor-universal neutral-current
    potential of the three active states is removed as a global phase.
    Because the sterile state does not feel :math:`V_{NC}`, removing it
    leaves :math:`-V_{NC}` on the sterile entry rather than nothing.

    `VCC` is positive for neutrinos.  For antineutrinos, reverse the
    sign of **both** potentials and conjugate the vacuum term, exactly
    as at three flavors.

    .. versionadded:: 1.9.0

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent four-flavor vacuum Hamiltonian, of shape
        ``(4, 4)``, in eV\ :sup:`2`.
    energy : int or float or array_like
        Neutrino energy, in eV, or an array of energies.
    VCC : int or float or array_like
        Charged-current matter potential, in eV.  Positive for
        neutrinos.
    VNC : int or float or array_like
        Neutral-current matter potential, in eV.  Negative for
        neutrinos, equal to :math:`-G_F n_n/\sqrt{2}`; see
        `globaldefs.VNC_EARTH_CRUST`.

    Returns
    -------
    numpy.ndarray
        Complex array of shape ``(4, 4)`` for a scalar energy, or
        ``(..., 4, 4)`` for an array of energies, in eV.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import globaldefs as gd
        import hamiltonians4nu
        import oscprob4nu

        h_vacuum = hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(
            gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF,
            np.sqrt(0.10), np.sqrt(0.10), 0.0,
            gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF, 1.0)

        h_matter = hamiltonians4nu.hamiltonian_4nu_matter(
            h_vacuum, 1.0e9, gd.VCC_EARTH_CRUST, gd.VNC_EARTH_CRUST)

        prob = oscprob4nu.probabilities_4nu(
            h_matter, 1300.0*gd.CONV_KM_TO_INV_EV)

        print('P_ee   = %.6f' % prob[0])
        print('P_mumu = %.6f' % prob[5])
    """
    h_vacuum = np.asarray(h_vacuum_energy_independent, dtype=complex)
    energy = np.asarray(energy, dtype=float)
    VCC = np.asarray(VCC, dtype=float)
    VNC = np.asarray(VNC, dtype=float)

    return (h_vacuum/energy[..., None, None]
            + VCC[..., None, None]*_EE_PROJECTOR_4NU
            - VNC[..., None, None]*_SS_PROJECTOR_4NU)


def hamiltonian_4nu_nsi(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float, list, np.ndarray],
    VCC: Union[int, float, list, np.ndarray],
    VNC: Union[int, float, list, np.ndarray],
    eps: Union[list, np.ndarray]
) -> np.ndarray:
    r"""Returns the four-neutrino Hamiltonian with matter and NSI.

    Adds non-standard interactions to the *active* block only.  Sterile
    states have no standard-model interactions by construction, so they
    have no non-standard ones either: the sterile row and column of the
    NSI matrix are zero, and the sterile entry keeps the
    :math:`-V_{NC}` of :func:`hamiltonian_4nu_matter`.

    The `eps` parameters follow
    :func:`hamiltonians3nu.hamiltonian_3nu_nsi`: six of them, with the
    three off-diagonal ones allowed to be complex.

    .. versionadded:: 1.9.0

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent four-flavor vacuum Hamiltonian, of shape
        ``(4, 4)``, in eV\ :sup:`2`.
    energy : int or float or array_like
        Neutrino energy, in eV, or an array of energies.
    VCC : int or float or array_like
        Charged-current matter potential, in eV.
    VNC : int or float or array_like
        Neutral-current matter potential, in eV.
    eps : array_like
        The six NSI strength parameters
        ``[eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt]``, relative
        to `VCC`.  The off-diagonal ones may be complex.

    Returns
    -------
    numpy.ndarray
        Complex array of shape ``(4, 4)`` for a scalar energy, or
        ``(..., 4, 4)`` for an array of energies, in eV.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import globaldefs as gd
        import hamiltonians4nu
        import oscprob4nu

        h_vacuum = hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(
            gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF,
            np.sqrt(0.10), np.sqrt(0.10), 0.0,
            gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF, 1.0)

        h_nsi = hamiltonians4nu.hamiltonian_4nu_nsi(
            h_vacuum, 1.0e9, gd.VCC_EARTH_CRUST, gd.VNC_EARTH_CRUST,
            gd.EPS_3)

        prob = oscprob4nu.probabilities_4nu(
            h_nsi, 1300.0*gd.CONV_KM_TO_INV_EV)

        print('P_mue with NSI = %.6f' % prob[4])
    """
    h_vacuum = np.asarray(h_vacuum_energy_independent, dtype=complex)
    energy = np.asarray(energy, dtype=float)
    VCC = np.asarray(VCC, dtype=float)
    VNC = np.asarray(VNC, dtype=float)

    eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt = eps

    # Complex throughout, so that complex off-diagonal parameters keep
    # their imaginary parts; the sterile row and column stay zero.
    nsi = np.array(
        [[1.0+eps_ee, eps_em, eps_et, 0.0],
         [np.conj(eps_em), eps_mm, eps_mt, 0.0],
         [np.conj(eps_et), np.conj(eps_mt), eps_tt, 0.0],
         [0.0, 0.0, 0.0, 0.0]], dtype=complex)

    return (h_vacuum/energy[..., None, None]
            + VCC[..., None, None]*nsi
            - VNC[..., None, None]*_SS_PROJECTOR_4NU)


def hamiltonian_4nu_liv(
    h_vacuum_energy_independent: Union[list, np.ndarray],
    energy: Union[int, float, list, np.ndarray],
    sxi12: Union[int, float],
    sxi23: Union[int, float],
    sxi13: Union[int, float],
    sxi14: Union[int, float],
    sxi24: Union[int, float],
    sxi34: Union[int, float],
    dxiCP: Union[int, float],
    b1: Union[int, float],
    b2: Union[int, float],
    b3: Union[int, float],
    b4: Union[int, float],
    Lambda: Union[int, float],
    angles: str = 'sin'
) -> np.ndarray:
    r"""Returns the four-neutrino Hamiltonian for oscillations w/ LIV.

    The four-flavor counterpart of
    :func:`hamiltonians3nu.hamiltonian_3nu_liv`.  The LIV term is
    :math:`(E/\Lambda) R B_4 R^\dagger`, with
    :math:`B_4 = \mathrm{diag}(b_1, b_2, b_3, b_4)` and :math:`R` a
    mixing matrix of the same 3+1 form as
    :func:`mixing_matrix_4nu`, built from the angles :math:`\xi_{ij}`
    and the phase :math:`\delta_{\xi,\rm CP}` that relate the
    eigenvectors of :math:`B_4` to the flavor states.

    Nothing here privileges the fourth state: :math:`b_4` is an
    eigenvalue like the others, so a sterile neutrino may couple to the
    LIV background whether or not it couples to matter.  Setting the
    three new angles to zero recovers the three-flavor term in the
    active block, whatever :math:`b_4` is: with :math:`R` block
    diagonal, :math:`b_4` acts on the sterile state alone, and a term
    on a decoupled state is a phase.

    .. versionadded:: 1.11.0

    Parameters
    ----------
    h_vacuum_energy_independent : array_like
        Energy-independent four-flavor vacuum Hamiltonian, of shape
        ``(4, 4)``, in eV\ :sup:`2`.  It is not modified.
    energy : int or float or array_like
        Neutrino energy, in eV, or an array of energies.
    sxi12 : int or float
        Sine of :math:`\xi_{12}`.
    sxi23 : int or float
        Sine of :math:`\xi_{23}`.
    sxi13 : int or float
        Sine of :math:`\xi_{13}`.
    sxi14 : int or float
        Sine of :math:`\xi_{14}`.
    sxi24 : int or float
        Sine of :math:`\xi_{24}`.
    sxi34 : int or float
        Sine of :math:`\xi_{34}`.
    dxiCP : int or float
        CP-violation phase of the LIV operator, in radian.
    b1 : int or float
        Eigenvalue :math:`b_1` of the LIV operator :math:`B_4` [eV].
    b2 : int or float
        Eigenvalue :math:`b_2` of the LIV operator :math:`B_4` [eV].
    b3 : int or float
        Eigenvalue :math:`b_3` of the LIV operator :math:`B_4` [eV].
    b4 : int or float
        Eigenvalue :math:`b_4` of the LIV operator :math:`B_4` [eV].
    Lambda : int or float
        Energy scale :math:`\Lambda` of the LIV operator [eV].

    angles : str, optional
        How the mixing angles are given: ``'sin'`` for their sines, the
        default and what earlier versions accepted; ``'sin2'`` for the
        squares of their sines, which is how global fits publish them;
        ``'rad'`` or ``'deg'`` for the angles themselves.  Under ``'rad'`` and
        ``'deg'`` the phases follow; under the other two they stay in
        radians, a phase having no sine to pass.  See
        `ANGLE_CONVENTIONS`.
    Returns
    -------
    numpy.ndarray
        Complex array of shape ``(4, 4)`` for a scalar energy, or
        ``(..., 4, 4)`` for an array of energies, in eV.

    Examples
    --------
    .. jupyter-execute::

        import numpy as np

        import globaldefs as gd
        import hamiltonians4nu
        import oscprob4nu

        h_vacuum = hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(
            gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF,
            np.sqrt(0.10), np.sqrt(0.10), 0.0,
            gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF, 1.0)

        h_liv = hamiltonians4nu.hamiltonian_4nu_liv(
            h_vacuum, 1.0e9, 0.3, 0.4, 0.5, 0.0, 0.0, 0.0, 0.7,
            1.0e-9, 1.5e-9, 2.0e-9, 2.5e-9, 1.0e12)

        prob = oscprob4nu.probabilities_4nu(
            h_liv, 1300.0*gd.CONV_KM_TO_INV_EV)

        print('P_ee with LIV = %.6f' % prob[0])
    """
    sxi12 = _sine_from(sxi12, angles, 'sxi12', 'hamiltonian_4nu_liv')
    sxi23 = _sine_from(sxi23, angles, 'sxi23', 'hamiltonian_4nu_liv')
    sxi13 = _sine_from(sxi13, angles, 'sxi13', 'hamiltonian_4nu_liv')
    sxi14 = _sine_from(sxi14, angles, 'sxi14', 'hamiltonian_4nu_liv')
    sxi24 = _sine_from(sxi24, angles, 'sxi24', 'hamiltonian_4nu_liv')
    sxi34 = _sine_from(sxi34, angles, 'sxi34', 'hamiltonian_4nu_liv')
    dxiCP = _phase_from(dxiCP, angles)
    h_vacuum = np.asarray(h_vacuum_energy_independent, dtype=complex)
    energy = np.asarray(energy, dtype=float)

    rotation = mixing_matrix_4nu(sxi12, sxi23, sxi13, sxi14, sxi24, sxi34,
                                 dxiCP)
    operator = np.diag([b1, b2, b3, b4]).astype(complex)
    liv = rotation @ operator @ rotation.conj().T

    factor = energy/Lambda

    return (h_vacuum/energy[..., None, None]
            + factor[..., None, None]*liv)
