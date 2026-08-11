# -*- coding: utf-8 -*-
r"""Physical constants, unit-conversion factors, and default parameters.

This module contains the values of the physical constants, the
unit-conversion factors, and the default oscillation parameters used by
the sample-Hamiltonian modules of **NuOscProbExact**, by :mod:`earth`,
and by the notebooks and worked examples.

The core modules :mod:`oscprob2nu`, :mod:`oscprob3nu` and
:mod:`oscprob4nu` do *not* need these constants: they accept an
arbitrary Hermitian Hamiltonian in whatever units the user prefers.

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

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

# Every constant this module defines, and nothing else.  The worked
# examples do `from globaldefs import *`, which without this also
# injects `np` into their namespace -- a name they neither asked for
# nor should shadow.  Every other module in src/ declares one.
__all__ = ['CONV_KM_TO_INV_EV', 'CONV_CM_TO_INV_EV', 'CONV_INV_EV_TO_CM',
           'CONV_EV_TO_G', 'CONV_G_TO_EV', 'GF', 'MASS_ELECTRON',
           'MASS_PROTON', 'MASS_NEUTRON', 'EARTH_RADIUS',
           'ELECTRON_FRACTION_EARTH_CRUST',
           'DENSITY_MATTER_CRUST_G_PER_CM3', 'NUM_DENSITY_E_EARTH_CRUST',
           'VCC_EARTH_CRUST', 'NEUTRON_FRACTION_EARTH_CRUST',
           'NUM_DENSITY_N_EARTH_CRUST', 'VNC_EARTH_CRUST',
           'ELECTRON_FRACTION_EARTH_CORE', 'ELECTRON_FRACTION_EARTH_MANTLE',
           'ELECTRON_FRACTION_EARTH_CRUST_LAYER',
           'ELECTRON_FRACTION_EARTH_OCEAN',
           'S12_NO_BF',
           'S23_NO_BF', 'S13_NO_BF', 'DCP_NO_BF', 'D21_NO_BF',
           'D31_NO_BF', 'S12_IO_BF', 'S23_IO_BF', 'S13_IO_BF',
           'DCP_IO_BF', 'D21_IO_BF', 'D32_IO_BF', 'D31_IO_BF', 'EPS_EE',
           'EPS_EM', 'EPS_ET', 'EPS_MM', 'EPS_MT', 'EPS_TT', 'EPS_2',
           'EPS_3', 'SXI12', 'SXI23', 'SXI13', 'DXICP', 'B1', 'B2', 'B3',
           'LAMBDA']

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

CONV_EV_TO_G = 1.78266192e-33
r"""float: Multiplicative conversion factor from eV to grams.

Converts a mass expressed in eV (i.e., eV/c\ :sup:`2`) into grams.
Units: [g eV\ :sup:`-1`].

.. versionchanged:: 1.11.0
   Given to the precision of the CODATA value, 1.78266192e-33, rather
   than rounded to 1.783e-33.  The rounded value was off by
   :math:`1.9 \times 10^{-4}` relative, three orders of magnitude worse
   than every other constant in this module, which sit between
   :math:`10^{-7}` and :math:`10^{-9}`.  It propagates through
   `NUM_DENSITY_E_EARTH_CRUST` into `VCC_EARTH_CRUST`, and through
   :func:`earth.matter_potential` into every Earth crossing, so the
   matter potential moved by that much.  Far below anything measurable,
   and still worth not carrying into a release.
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

NEUTRON_FRACTION_EARTH_CRUST = 1.0 - ELECTRON_FRACTION_EARTH_CRUST
r"""float: Neutrons per nucleon in the Earth's crust.

The crust is close to isoscalar, so with an electron fraction of one
half there is about one neutron per electron.  Units: [adimensional].
"""

NUM_DENSITY_N_EARTH_CRUST = NUM_DENSITY_E_EARTH_CRUST \
                            * NEUTRON_FRACTION_EARTH_CRUST \
                            / ELECTRON_FRACTION_EARTH_CRUST
r"""float: Neutron number density in the Earth's crust.

Units: [eV\ :sup:`3`].
"""


###############################################################################
# Electron fraction inside the Earth
###############################################################################

ELECTRON_FRACTION_EARTH_CORE = 0.4656
r"""float: Electron fraction in the Earth's core, :math:`r \leq 3480` km.

Iron, :math:`Z/A = 26/55.845`.  The core is neutron-richer than anything
above it, so this is the lowest of the four.

PREM\ [1]_ is a density model and carries no composition, so an electron
fraction has to come from somewhere else; these four values do, and each
is quoted with the material it belongs to.  Assuming one half everywhere
instead is exactly isoscalar matter, which no part of the Earth is.

None of the four is a default anywhere: `ELECTRON_FRACTION_EARTH_CRUST`,
one half throughout, remains what `earth` assumes unless asked
otherwise.  See `earth.electron_fraction_prem`.

Units: [adimensional].

References
----------
.. [1] A. M. Dziewonski and D. L. Anderson, "Preliminary reference Earth
   model", Phys. Earth Planet. Interiors 25, 297 (1981).
"""

ELECTRON_FRACTION_EARTH_MANTLE = 0.4957
r"""float: Electron fraction in the mantle, :math:`3480 < r \leq 6346.6` km.

Peridotite.  See `ELECTRON_FRACTION_EARTH_CORE` for why PREM does not
supply this.  Units: [adimensional].
"""

ELECTRON_FRACTION_EARTH_CRUST_LAYER = 0.4952
r"""float: Electron fraction in the crust, :math:`6346.6 < r \leq 6368` km.

Granitic, and within :math:`0.1\%` of the mantle's, so it is here for
explicitness rather than for effect.  Named apart from
`ELECTRON_FRACTION_EARTH_CRUST`, which is one half and is the uniform
value the code assumes by default.  Units: [adimensional].
"""

ELECTRON_FRACTION_EARTH_OCEAN = 0.5551
r"""float: Electron fraction in the ocean, :math:`r > 6368` km.

Seawater, :math:`Z/A = 10/18.015`.  Above one half, unlike every other
layer, because hydrogen has :math:`Z/A = 1`.

PREM's ocean is a global average, so it is a layer of the model rather
than of any particular baseline: a detector under rock has none.  For a
land chord, pass the crust's value here.  Units: [adimensional].
"""

VNC_EARTH_CRUST = -GF*NUM_DENSITY_N_EARTH_CRUST/np.sqrt(2.0)
r"""float: Neutral-current matter potential in the Earth's crust.

Equal to :math:`-G_F n_n/\sqrt{2}`, and **negative** for neutrinos.  It
is felt equally by all three active flavors, so at three flavors it is
proportional to the identity and drops out of the probabilities
entirely --- which is why :mod:`hamiltonians3nu` never needs it.

It matters as soon as a sterile state is present, because the sterile
state does not feel it: subtracting it from all four states, which costs
only a global phase, leaves :math:`-V_{NC}` on the sterile entry.  See
:func:`hamiltonians4nu.hamiltonian_4nu_matter`.

.. versionadded:: 1.9.0

Units: [eV].
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
# Total NSI strengths computed from the u and d quark parameters of [2]_.
#
# Read the numbers before reusing them: with EPS_EE = 0.06 and
# EPS_MM = 1.2, the combination that matter oscillations are actually
# sensitive to is eps_ee - eps_mm = -1.14, which is the LMA-D
# ("dark side") solution rather than the ordinary LMA one.  That is a
# deliberate choice of a large, visible effect for the worked examples
# and the notebooks, not a best fit: it makes the NSI curves differ from
# the standard ones by something a reader can see on a plot.  For a
# realistic study, substitute the constraint region you mean to explore.
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
