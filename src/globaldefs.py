# -*- coding: utf-8 -*-
r"""Physical constants, unit-conversion factors, and default parameters.

This module contains the values of the physical constants, the
unit-conversion factors, and the default oscillation parameters used by
the sample-Hamiltonian and plotting modules of **NuOscProbExact**.

The core modules :mod:`oscprob2nu` and :mod:`oscprob3nu` do *not* need
these constants: they accept an arbitrary Hermitian Hamiltonian in
whatever units the user prefers.

Unless stated otherwise, quantities are expressed in natural units, in
which energies are in eV, mass-squared differences in eV\ :sup:`2`, and
baselines in eV\ :sup:`-1`.

Notes
-----
The lepton mixing parameters are the best-fit values from NuFit 4.0 [1]_,
including Super-Kamiokande atmospheric data, for both the normal
ordering (suffix ``_NO_BF``) and the inverted ordering (suffix
``_IO_BF``).  The non-standard-interaction strengths are taken from
[2]_.  These constants support the method of [3]_.

References
----------

.. [1] I. Esteban *et al.*, "Global analysis of three-flavour neutrino
   oscillations", JHEP 01, 106 (2019), arXiv:1811.05487 (NuFit 4.0).

.. [2] P. Coloma *et al.*, "Curtailing the dark side in non-standard
   neutrino interactions", arXiv:1805.04530.

.. [3] Mauricio Bustamante, "Exact neutrino oscillation probabilities
   with arbitrary time-independent Hamiltonians", arXiv:1904.12391.
"""

from __future__ import print_function

__version__ = "1.1"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

import numpy as np


###############################################################################
# Unit-conversion factors
###############################################################################

CONV_KM_TO_INV_EV = 5.06773e9
r"""float: Multiplicative conversion factor from km to eV\ :sup:`-1`.

Units: [eV\ :sup:`-1` km\ :sup:`-1`].
"""

CONV_CM_TO_INV_EV = CONV_KM_TO_INV_EV*1.e-5
r"""float: Multiplicative conversion factor from cm to eV\ :sup:`-1`.

Units: [eV\ :sup:`-1` cm\ :sup:`-1`].
"""

CONV_INV_EV_TO_CM = 1./CONV_CM_TO_INV_EV
r"""float: Multiplicative conversion factor from eV\ :sup:`-1` to cm.

Units: [eV cm].
"""

CONV_EV_TO_G = 1.783e-33
r"""float: Multiplicative conversion factor from eV to grams.

Converts a mass expressed in eV (i.e., eV/c\ :sup:`2`) into grams.
Units: [g eV\ :sup:`-1`].
"""

CONV_G_TO_EV = 1./CONV_EV_TO_G
r"""float: Multiplicative conversion factor from grams to eV.

Units: [eV g\ :sup:`-1`].
"""


###############################################################################
# Physical constants
###############################################################################

GF = 1.1663787e-23
r"""float: Fermi constant.

Units: [eV\ :sup:`-2`].
"""

MASS_ELECTRON = 0.5109989461e6
r"""float: Electron mass.

Units: [eV].
"""

MASS_PROTON = 938.272046e6
r"""float: Proton mass.

Units: [eV].
"""

MASS_NEUTRON = 939.565379e6
r"""float: Neutron mass.

Units: [eV].
"""


###############################################################################
# Earth geometry
###############################################################################

EARTH_RADIUS = 6371.0
r"""float: Mean radius of the Earth.

The IUGG mean radius, which is what the Preliminary Reference Earth
Model in :mod:`earth` is normalised against and what the chord geometry
there assumes.  The Earth is treated as a sphere throughout; the
equatorial and polar radii differ from this by about 0.3%, which is far
below the accuracy of any matter-density model.

Units: [km].
"""


###############################################################################
# Matter density in the Earth's crust
###############################################################################

ELECTRON_FRACTION_EARTH_CRUST = 0.5
r"""float: Electron fraction in the Earth's crust.

Units: [adimensional].
"""

DENSITY_MATTER_CRUST_G_PER_CM3 = 3.0
r"""float: Average matter density in the Earth's crust.

Units: [g cm\ :sup:`-3`].
"""

NUM_DENSITY_E_EARTH_CRUST = DENSITY_MATTER_CRUST_G_PER_CM3 * CONV_G_TO_EV \
                            / ((MASS_PROTON+MASS_NEUTRON)/2.0) \
                            * ELECTRON_FRACTION_EARTH_CRUST \
                            / pow(CONV_CM_TO_INV_EV, 3.0)
r"""float: Electron number density in the Earth's crust.

Units: [eV\ :sup:`3`].
"""

VCC_EARTH_CRUST = np.sqrt(2.0)*GF*NUM_DENSITY_E_EARTH_CRUST
r"""float: Charged-current matter potential in the Earth's crust.

Equal to :math:`\sqrt{2} G_F n_e`.  It is positive for neutrinos; use
its negative for antineutrinos.  Units: [eV].
"""


###############################################################################
# Lepton mixing parameters, normal ordering
#
# Best-fit values from NuFit 4.0, with Super-Kamiokande atmospheric data
###############################################################################

S12_NO_BF = np.sqrt(0.310)
r"""float: Lepton mixing angle :math:`\sin\theta_{12}`, normal ordering.

Units: [adimensional].
"""

S23_NO_BF = np.sqrt(0.582)
r"""float: Lepton mixing angle :math:`\sin\theta_{23}`, normal ordering.

Units: [adimensional].
"""

S13_NO_BF = np.sqrt(2.240e-2)
r"""float: Lepton mixing angle :math:`\sin\theta_{13}`, normal ordering.

Units: [adimensional].
"""

DCP_NO_BF = 217./180.*np.pi
r"""float: CP-violation phase :math:`\delta_{\rm CP}`, normal ordering.

Units: [radian].
"""

D21_NO_BF = 7.39e-5
r"""float: Mass-squared difference :math:`\Delta m^2_{21}`, normal
ordering.

Units: [eV\ :sup:`2`].
"""

D31_NO_BF = 2.525e-3
r"""float: Mass-squared difference :math:`\Delta m^2_{31}`, normal
ordering.

Units: [eV\ :sup:`2`].
"""


###############################################################################
# Lepton mixing parameters, inverted ordering
#
# Best-fit values from NuFit 4.0, with Super-Kamiokande atmospheric data
###############################################################################

S12_IO_BF = np.sqrt(0.310)
r"""float: Lepton mixing angle :math:`\sin\theta_{12}`, inverted
ordering.

Units: [adimensional].
"""

S23_IO_BF = np.sqrt(0.582)
r"""float: Lepton mixing angle :math:`\sin\theta_{23}`, inverted
ordering.

Units: [adimensional].
"""

S13_IO_BF = np.sqrt(2.263e-2)
r"""float: Lepton mixing angle :math:`\sin\theta_{13}`, inverted
ordering.

Units: [adimensional].
"""

DCP_IO_BF = 280./180.*np.pi
r"""float: CP-violation phase :math:`\delta_{\rm CP}`, inverted
ordering.

Units: [radian].
"""

D21_IO_BF = 7.39e-5
r"""float: Mass-squared difference :math:`\Delta m^2_{21}`, inverted
ordering.

Units: [eV\ :sup:`2`].
"""

D32_IO_BF = -2.512e-3
r"""float: Mass-squared difference :math:`\Delta m^2_{32}`, inverted
ordering.

Units: [eV\ :sup:`2`].
"""

D31_IO_BF = D32_IO_BF+D21_IO_BF
r"""float: Mass-squared difference :math:`\Delta m^2_{31}`, inverted
ordering.

Computed as :math:`\Delta m^2_{32} + \Delta m^2_{21}`.
Units: [eV\ :sup:`2`].
"""


###############################################################################
# Non-standard interaction (NSI) parameters
#
# Total NSI strengths computed using values of the u and d quark
# parameters compatible at 2 sigma with LMA+coherent, from [2]_
###############################################################################

EPS_EE = 0.06
r"""float: NSI strength parameter :math:`\epsilon_{ee}`.

Units: [adimensional].
"""

EPS_EM = -0.06
r"""float: NSI strength parameter :math:`\epsilon_{e\mu}`.

May in general be complex.  Units: [adimensional].
"""

EPS_ET = 0.0
r"""float: NSI strength parameter :math:`\epsilon_{e\tau}`.

May in general be complex.  Units: [adimensional].
"""

EPS_MM = 1.2
r"""float: NSI strength parameter :math:`\epsilon_{\mu\mu}`.

Units: [adimensional].
"""

EPS_MT = 0.0
r"""float: NSI strength parameter :math:`\epsilon_{\mu\tau}`.

May in general be complex.  Units: [adimensional].
"""

EPS_TT = 0.0
r"""float: NSI strength parameter :math:`\epsilon_{\tau\tau}`.

Units: [adimensional].
"""

EPS_2 = [EPS_EE, EPS_EM, EPS_MM]
r"""list of float: NSI strengths for two-neutrino oscillations.

Ordered as ``[eps_ee, eps_em, eps_mm]``, ready to be passed to
:func:`hamiltonians2nu.hamiltonian_2nu_nsi`.
Units: [adimensional].
"""

EPS_3 = [EPS_EE, EPS_EM, EPS_ET, EPS_MM, EPS_MT, EPS_TT]
r"""list of float: NSI strengths for three-neutrino oscillations.

Ordered as ``[eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt]``, ready
to be passed to :func:`hamiltonians3nu.hamiltonian_3nu_nsi`.
Units: [adimensional].
"""


###############################################################################
# Lorentz invariance-violating (LIV) parameters
#
# Compatible with the 90% C.L. upper limits on c^(4) from 1709.03434
###############################################################################

SXI12 = 0.0
r"""float: LIV mixing angle :math:`\sin\xi_{12}`.

Units: [adimensional].
"""

SXI23 = 0.0
r"""float: LIV mixing angle :math:`\sin\xi_{23}`.

Units: [adimensional].
"""

SXI13 = 0.0
r"""float: LIV mixing angle :math:`\sin\xi_{13}`.

Units: [adimensional].
"""

DXICP = 0.0
r"""float: LIV CP-violation phase :math:`\delta_{\xi,\rm CP}`.

Units: [radian].
"""

B1 = 1.e-9
r"""float: Eigenvalue :math:`b_1` of the LIV operator.

Units: [eV].
"""

B2 = 1.e-9
r"""float: Eigenvalue :math:`b_2` of the LIV operator.

Units: [eV].
"""

B3 = 2.e-9
r"""float: Eigenvalue :math:`b_3` of the LIV operator.

Units: [eV].
"""

LAMBDA = 1.e12
r"""float: Energy scale :math:`\Lambda` of the LIV operator.

Units: [eV].
"""
