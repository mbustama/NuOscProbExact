"""import oscprob3nu
import hamiltonians3nu
from globaldefs import *

h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_vacuum_energy_independent(  S12_BF, S23_BF,
                                                                                S13_BF, DCP_BF,
                                                                                D21_BF, D31_BF)

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]


Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = oscprob3nu.probabilities_3nu( \
                                                    np.multiply(1./energy, h_vacuum_energy_indep),
                                                    baseline*CONV_KM_TO_INV_EV)

print("Pee = %6.5f, Pem = %6.5f, Pet = %6.5f" % (Pee, Pem, Pet))
print("Pme = %6.5f, Pmm = %6.5f, Pmt = %6.5f" % (Pme, Pmm, Pmt))
print("Pte = %6.5f, Ptm = %6.5f, Ptt = %6.5f" % (Pte, Ptm, Ptt))
"""


"""
import oscprob3nu
import hamiltonians3nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_vacuum_energy_independent(  S12_BF, S23_BF,
                                                                                S13_BF, DCP_BF,
                                                                                D21_BF, D31_BF)
h_vacuum = np.multiply(1./energy, h_vacuum_energy_indep)
h_matter = hamiltonians3nu.hamiltonian_matter(h_vacuum, VCC_EARTH_CRUST)

Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = oscprob3nu.probabilities_3nu( \
                                                h_matter, baseline*CONV_KM_TO_INV_EV)

print("Pee = %6.5f, Pem = %6.5f, Pet = %6.5f" % (Pee, Pem, Pet))
print("Pme = %6.5f, Pmm = %6.5f, Pmt = %6.5f" % (Pme, Pmm, Pmt))
print("Pte = %6.5f, Ptm = %6.5f, Ptt = %6.5f" % (Pte, Ptm, Ptt))
"""

import oscprob3nu
import hamiltonians3nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(  S12_BF, S23_BF,
                                                                                    S13_BF, DCP_BF,
                                                                                    D21_BF, D31_BF)
h_vacuum = np.multiply(1./energy, h_vacuum_energy_indep)

h_coeffs = oscprob3nu.hamiltonian_3nu_coefficients(h_vacuum)
u_coeffs = oscprob3nu.evolution_operator_3nu_u_coefficients(h_vacuum, baseline*CONV_KM_TO_INV_EV)
evol_op = oscprob3nu.evolution_operator_3nu(h_vacuum, baseline*CONV_KM_TO_INV_EV)

print('h_coeffs = ', h_coeffs)
print()
print('u_coeffs = ', u_coeffs)
print()
print('evol_op = ', evol_op)
