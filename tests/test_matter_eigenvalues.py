# -*- coding: utf-8 -*-
r"""Checks the matter Hamiltonian against an independent analytic spectrum.

The frozen nuSQuIDS comparison in ``test_nusquids_comparison.py`` tests
the conventions against an outside code, but only at the seven
configurations that were frozen.  This tests the same thing at *any*
configuration, live, with no external dependency, by using a published
closed form for the eigenvalues.

Zaglauer and Schwarzer [1]_ give the exact eigenvalues of

.. math:: M^2_{\rm eff} = U M^2 U^\dagger + \mathrm{diag}(A, 0, 0)

as the roots of an explicit cubic in the *vacuum* parameters and
:math:`A = 2 E V_{CC}`, without ever forming the matrix.  So comparing
their roots against the spectrum of the Hamiltonian that
:func:`hamiltonians3nu.hamiltonian_3nu_matter` builds tests the
construction rather than the expansion: the PMNS parametrization, where
the matter potential is added, and the relative normalization of the two
terms.  Nothing else in the suite does that --- ``scipy.linalg.expm``
takes our Hamiltonian as given.

The coefficients here are written from the mixing angles directly,

.. math::
   |U_{e1}|^2 = c_{12}^2 c_{13}^2 , \quad
   |U_{e2}|^2 = s_{12}^2 c_{13}^2 , \quad
   |U_{e3}|^2 = s_{13}^2 ,

rather than read out of :func:`hamiltonians3nu.pmns_mixing_matrix`, so
that a mistake in our mixing matrix cannot cancel against itself.

Note what this does *not* reach: the Zaglauer-Schwarzer eigenvalues are
independent of :math:`\theta_{23}` and :math:`\delta_{CP}`, which is a
property of the matter term rather than a limitation of the
transcription --- only the electron row of :math:`U` enters.  Those two
parameters are covered by the vacuum reference formulas in
``test_reference_formulas.py`` and by the nuSQuIDS comparison.

References
----------

.. [1] H.W. Zaglauer and K.H. Schwarzer, "The mixing angles in matter
   for three generations of neutrinos and the MSW mechanism",
   Z. Phys. C 40, 273 (1988).
"""

import numpy as np
import pytest

import globaldefs as gd
import hamiltonians3nu


def zaglauer_schwarzer_eigenvalues(s12, s13, dm21, dm31, A):
    r"""Returns the exact eigenvalues of the matter-effective mass matrix.

    Parameters
    ----------
    s12 : float
        Sine of :math:`\theta_{12}`.
    s13 : float
        Sine of :math:`\theta_{13}`.
    dm21 : float
        :math:`\Delta m^2_{21}`, in eV\ :sup:`2`.
    dm31 : float
        :math:`\Delta m^2_{31}`, in eV\ :sup:`2`.
    A : float
        :math:`2 E V_{CC}`, in eV\ :sup:`2`.  Negative for antineutrinos.

    Returns
    -------
    numpy.ndarray
        The three eigenvalues, ascending, in eV\ :sup:`2`.
    """
    c13_squared = 1.0 - s13*s13
    ue1 = (1.0 - s12*s12)*c13_squared
    ue2 = s12*s12*c13_squared
    ue3 = s13*s13

    # lambda^3 - x lambda^2 + y lambda - z = 0
    x = dm21 + dm31 + A
    y = dm21*dm31 + A*(dm21*(1.0 - ue2) + dm31*(1.0 - ue3))
    z = A*dm21*dm31*ue1

    # Depressed form, then the trigonometric solution: all three roots are
    # real because the matrix is Hermitian.
    p = y - x*x/3.0
    q = -2.0*x**3/27.0 + x*y/3.0 - z

    scale = 2.0*np.sqrt(-p/3.0)
    angle = np.arccos(np.clip(3.0*q/(p*scale), -1.0, 1.0))

    return np.sort(scale*np.cos((angle + 2.0*np.pi*np.arange(3))/3.0)
                   + x/3.0)


def our_eigenvalues(dm31, energy, VCC):
    r"""Returns the spectrum of the Hamiltonian this library builds.

    Parameters
    ----------
    dm31 : float
        :math:`\Delta m^2_{31}`, in eV\ :sup:`2`.
    energy : float
        Neutrino energy, in eV.
    VCC : float
        Charged-current potential, in eV.

    Returns
    -------
    numpy.ndarray
        The three eigenvalues of :math:`2 E H`, ascending, in
        eV\ :sup:`2`.
    """
    vacuum = np.asarray(
        hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
            gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, gd.DCP_NO_BF,
            gd.D21_NO_BF, dm31))
    hamiltonian = np.asarray(
        hamiltonians3nu.hamiltonian_3nu_matter(vacuum, energy, VCC))

    return np.sort(np.linalg.eigvalsh(2.0*energy*hamiltonian))


def deviation(dm31, energy, VCC):
    r"""Returns the relative deviation between the two spectra."""
    ours = our_eigenvalues(dm31, energy, VCC)
    theirs = zaglauer_schwarzer_eigenvalues(
        gd.S12_NO_BF, gd.S13_NO_BF, gd.D21_NO_BF, dm31, 2.0*energy*VCC)

    return np.max(np.abs(ours - theirs))/np.max(np.abs(ours))


def potential(density):
    r"""Returns ``VCC`` at a given density, in eV."""
    return (gd.VCC_EARTH_CRUST*density
            / gd.DENSITY_MATTER_CRUST_G_PER_CM3)


@pytest.mark.parametrize('energy_gev', [0.3, 1.0, 2.5, 10.0, 100.0])
@pytest.mark.parametrize('density', [0.5, 3.0, 5.5, 13.1])
def test_matter_spectrum_matches_zaglauer_schwarzer(energy_gev, density):
    r"""The built Hamiltonian has the analytically predicted spectrum."""
    assert deviation(gd.D31_NO_BF, energy_gev*1.0e9,
                     potential(density)) < 1.0e-13


def test_it_holds_for_antineutrinos():
    r"""A negative potential is handled the same way.

    The sign of ``VCC`` is the whole difference between neutrinos and
    antineutrinos in matter, so an independent check of the spectrum at
    negative ``A`` is worth having on its own.
    """
    for energy_gev in (1.0, 5.0, 25.0):
        assert deviation(gd.D31_NO_BF, energy_gev*1.0e9,
                         -potential(3.0)) < 1.0e-13


def test_it_holds_for_the_inverted_ordering():
    r"""A negative mass-squared difference is handled the same way."""
    for energy_gev in (1.0, 5.0, 25.0):
        assert deviation(gd.D31_IO_BF, energy_gev*1.0e9,
                         potential(3.0)) < 1.0e-13


def test_it_holds_through_the_msw_resonance():
    r"""Across the resonance, where the eigenvalues come closest.

    The resonance sits near :math:`A = \Delta m^2_{31}\cos 2\theta_{13}`,
    which is where two eigenvalues approach each other and any error in
    the construction is most visible.
    """
    c2th13 = 1.0 - 2.0*gd.S13_NO_BF**2
    VCC = potential(3.0)
    resonant_energy = gd.D31_NO_BF*c2th13/(2.0*VCC)

    for factor in (0.5, 0.9, 1.0, 1.1, 2.0):
        assert deviation(gd.D31_NO_BF, resonant_energy*factor, VCC) < 1.0e-12


def test_zero_potential_recovers_the_vacuum_splittings():
    r"""At ``A = 0`` the eigenvalues are the vacuum ones."""
    eigenvalues = zaglauer_schwarzer_eigenvalues(
        gd.S12_NO_BF, gd.S13_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF, 0.0)

    assert np.allclose(eigenvalues,
                       np.sort([0.0, gd.D21_NO_BF, gd.D31_NO_BF]),
                       rtol=1.0e-12, atol=1.0e-18)
