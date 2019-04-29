# -*- coding: utf-8 -*-
r"""Run the trivial 3nu example shown in README.md.

Runs the trivial three-neutrino example shown in README.md

References
----------

.. [1] Mauricio Bustamante, "Exact neutrino oscillation probabilities:
   a fast general-purpose computation method for two and three neutrino
   flavors", arXiv:1904.XXXXX.

Created: 2019/04/29 23:22
Last modified: 2019/04/29 23:22
"""


from __future__ import print_function

__version__ = "1.0"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import sys
sys.path.append('../src')

import oscprob3nu

hamiltonian = [
                [1.0+0.0j, 0.0+2.0j, 0.0-1.0j],
                [0.0-2.0j, 3.0+0.0j, 3.0+0.0j],
                [0.0+1.0j, 3.0-0.0j, 5.0+0.0j]
]

L = 1.0

Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = \
    oscprob3nu.probabilities_3nu(hamiltonian, L)

print("Pee = %6.5f, Pem = %6.5f, Pet = %6.5f" % (Pee, Pem, Pet))
print("Pme = %6.5f, Pmm = %6.5f, Pmt = %6.5f" % (Pme, Pmm, Pmt))
print("Pte = %6.5f, Ptm = %6.5f, Ptt = %6.5f" % (Pte, Ptm, Ptt))
