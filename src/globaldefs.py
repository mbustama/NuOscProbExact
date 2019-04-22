#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = "0.1"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


"""
globaldefs.py:
    Contains global constants

Created: 2019/04/17 17:03
Last modified: 2019/04/17 17:03
"""


from numpy import *
import numpy as np


# Conversion factor from km to eV^{-1}
CONV_KM_TO_INV_EV = 1./1.23984e-9 # [km^{-1} eV]

# Conversion factor from eV^{-1} to cm
CONV_INV_EV_TO_CM = 1.23984e-4 # [cm^{-1} / eV^{-1}]
CONV_CM_TO_INV_EV = 1./CONV_INV_EV_TO_CM # [eV^{-1} cm^{-1}]

# Conversion factor from eV to g
CONV_EV_TO_G = 1.783e-33 # [g eV^{-1}]
CONV_G_TO_EV = 1./CONV_EV_TO_G # [eV g^{-1}]

# Fermi constant
GF = 1.1663787e-23 # [eV^-2]

# Masses
MASS_ELECTRON = 0.5109989461e6 # [eV]
MASS_PROTON = 938.272046e6 # [eV]
MASS_NEUTRON = 939.565379e6 # [eV]

# Number density of electrons in the Earth's crust
ELECTRON_FRACTION_EARTH_CRUST = 0.5
DENSITY_MATTER_CRUST_G_PER_CM3 = 3.0 # [g cm^{-3}]
NUM_DENSITY_E_EARTH_CRUST = DENSITY_MATTER_CRUST_G_PER_CM3 * CONV_G_TO_EV \
                            / ((MASS_PROTON+MASS_NEUTRON)/2.0) \
                            * ELECTRON_FRACTION_EARTH_CRUST \
                            / pow(CONV_CM_TO_INV_EV, 3.0) # [eV^3]

# Matter potential in the Earth's crust
VCC_EARTH_CRUST = sqrt(2.0)*GF*NUM_DENSITY_E_EARTH_CRUST # [eV]

# Lepton mixing parameters
# Best-fit values of mixing parameters, normal ordering
# From NuFit 4.0 with SK atmospheric data
S12_BF = sqrt(0.310)
S23_BF = sqrt(0.582)
S13_BF = sqrt(2.240e-2)
DCP_BF = 217./180.*np.pi # [rad]
D21_BF = 7.39e-5 # [eV^2]
D31_BF = 2.525e-3 # [eV^2]

# NSI parameters
# Compatible with the 2sigma LMA+COHERENT ranges from 1805.04530
EPS_EE = 0.06
EPS_EM = -0.06
EPS_ET = 0.0#-0.6
EPS_MM = 0.0#0.6
EPS_MT = -0.06
EPS_TT = 0.0#0.6
EPS_2 = [EPS_EE, EPS_EM, EPS_MM]
EPS_3 = [EPS_EE, EPS_EM, EPS_ET, EPS_MM, EPS_MT, EPS_TT]

# LIV parameters
# Compatible with 90% C.L. upper limits on c^(4) from 1709.03434
SXI12 = 0.0
SXI23 = 0.0
SXI13 = 0.0
DXICP = 0.0
B1 = 1.e-9 #0.0
B2 = 1.e-9 #1.e-16
B3 = 5.e-9 #0.0
LAMBDA = 1.e12 # [eV]

