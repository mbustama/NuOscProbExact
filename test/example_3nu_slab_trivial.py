# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 13:19:13 2023

@author: Gianfranco
"""

import sys
sys.path.append('../src')

import oscprob3nu_slab

hamiltonian_1 =[
                [1.0+0.0j, 1.0+2.0j, 3.0+0.0j],
                [1.0-2.0j, 3.0+0.0j, 1.0+2.0j],
                [3.0+0.0j, 1.0-2.0j, 1.0+0.0j]
]

hamiltonian_2 =[
                [2.0+0.0j, 2.0+4.0j, 6.0+0.0j],
                [2.0-4.0j, 6.0+0.0j, 2.0+4.0j],
                [6.0+0.0j, 2.0-4.0j, 2.0+0.0j]
]

hamiltonian_matrices = [hamiltonian_1, hamiltonian_2]

slab_start = [0.0, 0.5]

width_final_slab = 1

baseline = 1

Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = \
    oscprob3nu_slab.probabilities_3nu_slabs(hamiltonian_matrices, slab_start, width_final_slab, baseline)

print("Pee = %6.5f, Pem = %6.5f, Pet = %6.5f" % (Pee, Pem, Pet))
print("Pme = %6.5f, Pmm = %6.5f, Pmt = %6.5f" % (Pme, Pmm, Pmt))
print("Pte = %6.5f, Ptm = %6.5f, Ptt = %6.5f" % (Pte, Ptm, Ptt))




