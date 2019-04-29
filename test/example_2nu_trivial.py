# -*- coding: utf-8 -*-
r"""Run the trivial 2nu example shown in README.md.

Runs the trivial two-neutrino example shown in README.md

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

import oscprob2nu


hamiltonian = [
                [1.0+0.0j, 1.0+2.0j],
                [1.0-2.0j, 3.0+0.0j]
]

L = 1.0

Pee, Pem, Pme, Pmm = oscprob2nu.probabilities_2nu(hamiltonian, L)

print("Pee = %6.5f, Pem = %6.5f" % (Pee, Pem))
print("Pme = %6.5f, Pmm = %6.5f" % (Pme, Pmm))
