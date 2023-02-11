# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 16:15:05 2022

@author: Gianfranco Colone
"""

import numpy as np
import oscprob3nu 
import hamiltonians3nu
from globaldefs import *

unit_matrix_3=[
    [1,0,0],
    [0,1,0],    
    [0,0,1]
]

def evolution_operator_3nu_slabs(hamiltonian_matrices, slab_start, 
    width_final_slab, L):

    try:
        assert hamiltonian_matrices is not None
        assert slab_start is not None
        assert width_final_slab is not None
        assert L is not None
    except AssertionError:
        print(evolution_operator_3nu_slabs.__name__ \
            + ": One or more of the input arguments is None")

    try:
        assert L >= slab_start[0]
    except AssertionError:
        print(evolution_operator_3nu_slabs.__name__ \
            + ": Requested L is smaller than start of first slab")

    try:
        assert L <= slab_start[len(slab_start)-1]+width_final_slab
    except AssertionError:
        print(evolution_operator_3nu_slabs.__name__ \
            + ": Requested L is larger than end of last slab")
            
    slab_start = np.append(slab_start, slab_start[-1] + width_final_slab)
    
    # Slabs up to L
    slab_start = [x for x in slab_start if x < L]
    num_slabs = len(slab_start)

    # Find the end coordinates of the slabs
    slab_end = [slab_start[i+1] if (i < num_slabs-1) else L
        for i in range(num_slabs)]

    # Multiply the evolution operator of each slab
    U = unit_matrix_3
    for i in range(num_slabs):
        slab_width = slab_end[i]-slab_start[i]
        Unew = oscprob3nu.evolution_operator_3nu(hamiltonian_matrices[i], 
                slab_width)
        U = np.matmul(Unew, U)

    return U

    
def probabilities_3nu_slabs(hamiltonian_matrices, slabs_initial, 
    width_final_slab, L):
    
    U = evolution_operator_3nu_slabs(hamiltonian_matrices, slabs_initial, 
        width_final_slab, L)

    Pee = abs(U[0][0])**2.
    Pem = abs(U[1][0])**2.
    Pet = abs(U[2][0])**2.
    Pme = abs(U[0][1])**2.
    Pmm = abs(U[1][1])**2.
    Pmt = abs(U[2][1])**2.
    Pte = abs(U[0][2])**2.
    Ptm = abs(U[1][2])**2.
    Ptt = abs(U[2][2])**2.

    return Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt
