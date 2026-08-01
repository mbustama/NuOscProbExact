# -*- coding: utf-8 -*-
r"""Sample four-neutrino Hamiltonians, for 3+1 scenarios.

This module builds the :math:`4\times4` Hamiltonians that
:mod:`oscprob4nu` evaluates, in the same spirit as
:mod:`hamiltonians3nu` does at three flavors: they are *examples*, not
limitations.  :func:`oscprob4nu.probabilities_4nu` takes any Hermitian
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

References
----------

.. [1] Mauricio Bustamante, "Exact neutrino oscillation probabilities
   with arbitrary time-independent Hamiltonians", arXiv:1904.12391.
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

__all__ = ['mixing_matrix_4nu',
           'hamiltonian_4nu_vacuum_energy_independent',
           'hamiltonian_4nu_matter', 'hamiltonian_4nu_nsi']

from typing import Union

import numpy as np


_EE_PROJECTOR_4NU = np.diag([1.0, 0.0, 0.0, 0.0])
_SS_PROJECTOR_4NU = np.diag([0.0, 0.0, 0.0, 1.0])


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
    cth = np.sqrt(1.0 - sth*sth)
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
    d24: Union[int, float] = 0.0
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
    d24: Union[int, float] = 0.0
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
