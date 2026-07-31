.. NuOscProbExact documentation master file

NuOscProbExact: Exact Neutrino Oscillation Probabilities
=========================================================

.. image:: https://img.shields.io/badge/arXiv-1904.12391-orange.svg
   :target: https://arxiv.org/abs/1904.12391
   :alt: arXiv:1904.12391

.. image:: https://zenodo.org/badge/182178323.svg
   :target: https://zenodo.org/badge/latestdoi/182178323
   :alt: DOI

.. image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License: MIT

.. image:: https://img.shields.io/badge/python-3.9+-blue.svg
   :target: https://www.python.org/downloads/
   :alt: Python 3.9+

.. important::
   **Important links:**

   * `GitHub repository <https://github.com/mbustama/NuOscProbExact>`_
   * `The paper <https://arxiv.org/abs/1904.12391>`_ (arXiv:1904.12391)
   * :doc:`changelog`

**NuOscProbExact** computes exact two-flavor and three-flavor neutrino
oscillation probabilities for arbitrary time-independent Hamiltonians.  It is
a Python implementation of the method of Ohlsson and Snellman
:cite:`Ohlsson:1999xb`, revisited in :cite:`Bustamante:2019ggq`.

The method expands the Hamiltonian and the time-evolution operator in the
bases of SU(2) and SU(3) matrices, which yields concise, analytical, and
**exact** closed-form expressions for the probabilities.

.. _what-exact-means:

What "exact" means here
------------------------

There is no approximation beyond floating-point round-off.  The evolution
operator produced by the SU(2) and SU(3) expansions is the matrix exponential
:math:`e^{-i H_0 L}`, and the test suite checks that it is, against an
independent computation:

.. list-table::
   :header-rows: 1
   :widths: 60 40

   * - Property
     - Measured agreement
   * - :math:`U_3` vs ``scipy.linalg.expm``, 200 random Hermitian Hamiltonians
     - 8e-15
   * - :math:`U_2` vs ``scipy.linalg.expm``, 200 random Hermitian Hamiltonians
     - 1.3e-15
   * - Unitarity, :math:`U^\dagger U - \mathbb{1}`
     - 5e-15
   * - Hard-coded :math:`d_{ijk}` vs :math:`\frac{1}{4}\mathrm{Tr}(\{\lambda_i,\lambda_j\}\lambda_k)`, all 512 entries
     - 2e-16
   * - Exact result vs the standard vacuum oscillation formula
     - 4e-19

When is NuOscProbExact a good fit?
-----------------------------------

#. **Your Hamiltonian is time-independent.**  The method assumes it, and in
   exchange gives a closed form rather than a numerical integration.  For a
   Hamiltonian that varies along the trajectory, propagate slab by slab or use
   a method built for it.

#. **You want an arbitrary Hamiltonian, not a fixed scenario.**  The core
   routines take any Hermitian :math:`2\times2` or :math:`3\times3` matrix.
   Non-standard interactions, Lorentz-invariance violation, sterile-like
   perturbations and matter effects are all just entries in that matrix; the
   bundled Hamiltonians in :mod:`hamiltonians2nu` and :mod:`hamiltonians3nu`
   are examples, not limitations.

#. **You need the evolution operator, not only the probabilities.**
   :func:`oscprob3nu.evolution_operator_3nu` returns :math:`U_3(L)` itself, so
   it can be composed across segments or used to propagate a density matrix.

Performance
-----------

A single probability takes about 13 microseconds for three flavors and 1.3
for two.  Scans --- a curve versus baseline or energy, an oscillogram over
both --- are what the code mostly does, and two things make those much
faster without changing any answer.

**Pass arrays instead of looping.**  Every core routine takes a stack of
Hamiltonians, an array of baselines, or both, and evaluates them in one
call.  That is worth roughly 25 to 60 times the equivalent Python loop, and
needs no extra dependency: the expensive part of the expansion, the
characteristic equation whose roots give the oscillation phases, depends on
the Hamiltonian alone, so a scan over baselines solves it once instead of
once per point.  See :ref:`scanning` for the three patterns.

**Install Numba, if the scans are large.**  With
``pip install "nuoscprobexact[fast]"`` the batched paths run as compiled
loops spread over the available cores; without it the NumPy path is used and
the results are the same to round-off.  On 2000-point scans, against the
Python loop:

.. list-table::
   :header-rows: 1
   :widths: 40 20 20 20

   * - Scan
     - Loop
     - Arrays
     - Arrays + Numba
   * - Three flavors, versus baseline
     - 30 ms
     - 1.1 ms
     - **0.24 ms**
   * - Three flavors, versus energy
     - 41 ms
     - 1.6 ms
     - **0.43 ms**
   * - Oscillogram, 100 x 100
     - 195 ms
     - 4.4 ms
     - **0.78 ms**
   * - Two flavors, versus baseline
     - 4.4 ms
     - **0.07 ms**
     - not used

Repeated runs on one machine vary by up to half again, so read these as
orders of magnitude.

The last row is not an omission.  The backend is used only where it has been
measured to win: for three flavors that is every stack size, but the
two-flavor expansion reduces to a square root and a sine per element, which
NumPy already does about as well as compiled code can.  Below fifty thousand
elements the NumPy path is kept; above it the kernel leads by about 1.3 to
1.8 times.  The library chooses without being asked.

:doc:`methodology` explains where the time goes, and what was tried and
rejected.

Getting started
---------------

.. code-block:: shell

   git clone https://github.com/mbustama/NuOscProbExact.git
   cd NuOscProbExact
   pip install -e .

.. code-block:: python

   import numpy as np
   import oscprob3nu
   import hamiltonians3nu
   from globaldefs import *

   h_vacuum_energy_indep = \
       hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
           S12_NO_BF, S23_NO_BF, S13_NO_BF, DCP_NO_BF, D21_NO_BF, D31_NO_BF)
   h_vacuum = np.multiply(1./1.e9, h_vacuum_energy_indep)

   Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = \
       oscprob3nu.probabilities_3nu(h_vacuum, 1.3e3*CONV_KM_TO_INV_EV)

See :doc:`installation` and :doc:`quickstart` for the longer version.

.. toctree::
   :maxdepth: 2
   :caption: User Guide:

   installation
   quickstart
   methodology
   functions
   references
   changelog

Citing
------

If you use **NuOscProbExact** in your work, please cite
:cite:`Bustamante:2019ggq`.  The BibTeX entry is on the
:doc:`references` page, and INSPIRE keeps
`an up-to-date record <https://inspirehep.net/literature/1731803>`_.

License
-------

**NuOscProbExact** is released under the `MIT License
<https://opensource.org/licenses/MIT>`_.  The full text ships with the source,
as ``LICENSE`` in the repository root.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
