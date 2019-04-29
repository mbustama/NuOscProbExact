# -*- coding: utf-8 -*-
r"""Contains physical constants and unit-conversion constants.

This module contains contains values of physical constants and
unit-conversion factors used by the various modules of NuOscProbExact.
The core modules oscprob2nu.py and oscprob3nu.py do not require these
constants.

Created: 2019/04/17 17:03
Last modified: 2019/04/29 23:39
"""


__version__ = "1.0"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


from numpy import *
import numpy as np


CONV_KM_TO_INV_EV = 5.06773e9
r"""float: Module-level constant

Multiplicative conversion factor from km to eV^{-1}.
Units: [km^{-1} eV^{-1}].
"""

CONV_CM_TO_INV_EV = CONV_KM_TO_INV_EV*1.e-5
r"""float: Module-level constant

Multiplicative conversion factor from cm to eV^{-1}.
Units: [cm^{-1} eV^{-1}]
"""

CONV_INV_EV_TO_CM = 1./CONV_CM_TO_INV_EV
r"""float: Module-level constant

Multiplicative conversion factor from eV^{-1} to cm.
Units: [eV cm]
"""

CONV_EV_TO_G = 1.783e-33
r"""float: Module-level constant

Multiplicative conversion factor from eV^{-1} to grams.
Units: [g eV^{-1}]
"""

CONV_G_TO_EV = 1./CONV_EV_TO_G
r"""float: Module-level constant

Multiplicative conversion factor from grams to eV^{-1}.
Units: [eV g^{-1}]
"""

GF = 1.1663787e-23
r"""float: Module-level constant

Fermi constant.
Units: [eV^{-1}]
"""

MASS_ELECTRON = 0.5109989461e6
r"""float: Module-level constant

Electron mass.
Units: [eV]
"""

MASS_PROTON = 938.272046e6
r"""float: Module-level constant

Proton mass.
Units: [eV]
"""

MASS_NEUTRON = 939.565379e6
r"""float: Module-level constant

Neutron mass.
Units: [eV]
"""

ELECTRON_FRACTION_EARTH_CRUST = 0.5
r"""float: Module-level constant

Electron fraction in the Earth's crust.
Units: [Adimensional]
"""

DENSITY_MATTER_CRUST_G_PER_CM3 = 3.0
r"""float: Module-level constant

Average matter density in the Earth's crust.
Units: [g cm^{-3}]
"""

# NUM_DENSITY_E_EARTH_CRUST = DENSITY_MATTER_CRUST_G_PER_CM3 * CONV_G_TO_EV \
#                             / ((MASS_PROTON+MASS_NEUTRON)/2.0) \
#                             * ELECTRON_FRACTION_EARTH_CRUST \
#                             / pow(CONV_CM_TO_INV_EV, 3.0)
NUM_DENSITY_E_EARTH_CRUST = DENSITY_MATTER_CRUST_G_PER_CM3 * CONV_G_TO_EV \
                            / ((MASS_PROTON+MASS_NEUTRON)/2.0) \
                            * ELECTRON_FRACTION_EARTH_CRUST \
                            / pow(CONV_CM_TO_INV_EV, 3.0)
r"""float: Module-level constant

Electron number density in the Earth's crust
Units: [eV^3]
"""

VCC_EARTH_CRUST = sqrt(2.0)*GF*NUM_DENSITY_E_EARTH_CRUST
r"""float: Module-level constant

Charged-current matter potential in the Earth's crust.
Units: [eV]
"""

S12_NO_BF = sqrt(0.310)
r"""float: Module-level constant

Lepton mixing angle sin(theta_12), best fit from NuFit 4.0, assuming
normal ordering with SK atmospheric data.
Units: [Adimensional]
"""

S23_NO_BF = sqrt(0.582)
r"""float: Module-level constant

Lepton mixing angle sin(theta_23), best fit from NuFit 4.0, assuming
normal ordering with SK atmospheric data.
Units: [Adimensional]
"""

S13_NO_BF = sqrt(2.240e-2)
r"""float: Module-level constant

Lepton mixing angle sin(theta_13), best fit from NuFit 4.0, assuming
normal ordering with SK atmospheric data.
Units: [Adimensional]
"""

DCP_NO_BF = 217./180.*np.pi
r"""float: Module-level constant

Lepton CP-violation phase delta_CP, best fit from NuFit 4.0, assuming
normal ordering with SK atmospheric data.
Units: [radian]
"""

D21_NO_BF = 7.39e-5
r"""float: Module-level constant

Mass-squared difference Delta m^2_21, best fit from NuFit 4.0, assuming
normal ordering with SK atmospheric data.
Units: [eV^2]
"""

D31_NO_BF = 2.525e-3
r"""float: Module-level constant

Mass-squared difference Delta m^2_31, best fit from NuFit 4.0, assuming
normal ordering with SK atmospheric data.
Units: [eV^2]
"""

S12_IO_BF = sqrt(0.310)
r"""float: Module-level constant

Lepton mixing angle sin(theta_12), best fit from NuFit 4.0, assuming
inverted ordering with SK atmospheric data.
Units: [Adimensional]
"""

S23_IO_BF = sqrt(0.582)
r"""float: Module-level constant

Lepton mixing angle sin(theta_23), best fit from NuFit 4.0, assuming
inverted ordering with SK atmospheric data.
Units: [Adimensional]
"""

S13_IO_BF = sqrt(2.263e-2)
r"""float: Module-level constant

Lepton mixing angle sin(theta_13), best fit from NuFit 4.0, assuming
inverted ordering with SK atmospheric data.
Units: [Adimensional]
"""

DCP_IO_BF = 280./180.*np.pi
r"""float: Module-level constant

Lepton CP-violation phase delta_CP, best fit from NuFit 4.0, assuming
inverted ordering with SK atmospheric data.
Units: [radian]
"""

D21_IO_BF = 7.39e-5
r"""float: Module-level constant

Mass-squared difference Delta m^2_21, best fit from NuFit 4.0, assuming
normal ordering with SK atmospheric data.
Units: [eV^2]
"""

D32_IO_BF = -2.512e-3
r"""float: Module-level constant

Mass-squared difference Delta m^2_32, best fit from NuFit 4.0, assuming
normal ordering with SK atmospheric data.
Units: [eV^2]
"""

D31_IO_BF = D32_IO_BF+D21_IO_BF
r"""float: Module-level constant

Mass-squared difference Delta m^2_31, best fit from NuFit 4.0, assuming
inverted ordering with SK atmospheric data.
Units: [eV^2]
"""

EPS_EE = 0.06
r"""float: Module-level constant

Total NSI strength parameter eps_ee computed using values of the u and d
quark parameters compatible at 2sigma with LMA+coherent from 1805.04530.
Units: [Adimensional]
"""

EPS_EM = -0.06
r"""float: Module-level constant

Total NSI strength parameter eps_em computed using values of the u and d
quark parameters compatible at 2sigma with LMA+coherent from 1805.04530.
Units: [Adimensional]
"""

EPS_ET = 0.0
r"""float: Module-level constant

Total NSI strength parameter eps_et computed using values of the u and d
quark parameters compatible at 2sigma with LMA+coherent from 1805.04530.
Units: [Adimensional]
"""

EPS_MM = 1.2
r"""float: Module-level constant

Total NSI strength parameter eps_mm computed using values of the u and d
quark parameters compatible at 2sigma with LMA+coherent from 1805.04530.
Units: [Adimensional]
"""

EPS_MT = 0.0
r"""float: Module-level constant

Total NSI strength parameter eps_mt computed using values of the u and d
quark parameters compatible at 2sigma with LMA+coherent from 1805.04530.
Units: [Adimensional]
"""

EPS_TT = 0.0
r"""float: Module-level constant

Total NSI strength parameter eps_tt computed using values of the u and d
quark parameters compatible at 2sigma with LMA+coherent from 1805.04530.
Units: [Adimensional]
"""

EPS_2 = [EPS_EE, EPS_EM, EPS_MM]
r"""float: Module-level constant

Vector of total NSI strength parameters for two-neutrino oscillations.
Used in oscprob2nu_plot.py.
Units: [Adimensional]
"""

EPS_3 = [EPS_EE, EPS_EM, EPS_ET, EPS_MM, EPS_MT, EPS_TT]
r"""float: Module-level constant

Vector of total NSI strength parameters for three-neutrino oscillations.
Used in oscprob3nu_plot.py.
Units: [Adimensional]
"""

# LIV parameters
# Compatible with 90% C.L. upper limits on c^(4) from 1709.03434
SXI12 = 0.0
r"""float: Module-level constant

LIV lepton mixing angle sin(xi_12).
Units: [Adimensional]
"""

SXI23 = 0.0
r"""float: Module-level constant

LIV lepton mixing angle sin(xi_23).
Units: [Adimensional]
"""

SXI13 = 0.0
r"""float: Module-level constant

LIV lepton mixing angle sin(xi_13).
Units: [Adimensional]
"""

DXICP = 0.0
r"""float: Module-level constant

LIV CP-violation phase.
Units: [radian]
"""

B1 = 1.e-9
r"""float: Module-level constant

LIV eigenvalue b_1.
Units: [eV]
"""

B2 = 1.e-9
r"""float: Module-level constant

LIV eigenvalue b_2.
Units: [eV]
"""
B3 = 2.e-9
r"""float: Module-level constant

LIV eigenvalue b_3.
Units: [eV]
"""

LAMBDA = 1.e12 # [eV]
r"""float: Module-level constant

LIV energy scale Lambda.
Units: [eV]
"""
