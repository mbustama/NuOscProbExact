.. NuOscProbExact documentation master file

NuOscProbExact: Exact Neutrino Oscillation Probabilities
=========================================================

.. image:: https://github.com/mbustama/NuOscProbExact/actions/workflows/tests.yml/badge.svg
   :target: https://github.com/mbustama/NuOscProbExact/actions/workflows/tests.yml
   :alt: tests

.. image:: https://github.com/mbustama/NuOscProbExact/actions/workflows/lint.yml/badge.svg
   :target: https://github.com/mbustama/NuOscProbExact/actions/workflows/lint.yml
   :alt: Code Quality

.. image:: https://codecov.io/gh/mbustama/NuOscProbExact/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/mbustama/NuOscProbExact
   :alt: codecov

.. image:: https://img.shields.io/badge/docs-GitHub%20Pages-blue.svg
   :target: https://mbustama.github.io/NuOscProbExact/
   :alt: Documentation

.. image:: https://img.shields.io/pypi/v/nuoscprobexact.svg
   :target: https://pypi.org/project/nuoscprobexact/
   :alt: PyPI

.. image:: https://pepy.tech/badge/nuoscprobexact
   :target: https://pepy.tech/project/nuoscprobexact
   :alt: Downloads

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

.. image:: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json
   :target: https://github.com/astral-sh/ruff
   :alt: Code style: ruff

.. important::
   **Important links:**

   * :doc:`What it can compute, with code <recipes>`
   * `GitHub repository <https://github.com/mbustama/NuOscProbExact>`_
   * `The paper <https://arxiv.org/abs/1904.12391>`_ (arXiv:1904.12391)
   * `The notebooks <https://github.com/mbustama/NuOscProbExact/tree/main/notebooks>`_
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

#. **Your Hamiltonian is constant, or piecewise constant.**  The method
   assumes a Hamiltonian that does not change, and in exchange gives a closed
   form rather than a numerical integration.  A trajectory made of pieces is
   handled by :mod:`slabs`, which solves each piece exactly and multiplies the
   operators, and :mod:`earth` builds those pieces from the Preliminary
   Reference Earth Model.

   A profile that varies *smoothly* over an oscillation length is the case to
   avoid: it can be slabbed, but the step size is then set by the oscillation
   rather than by the density, which for the Sun means of order
   :math:`10^4` slabs per resonance crossing.  Use a Magnus-type method
   there.

#. **You want an arbitrary Hamiltonian, not a fixed scenario.**  The core
   routines take any Hermitian :math:`2\times2` or :math:`3\times3` matrix.
   Non-standard interactions, Lorentz-invariance violation, sterile-like
   perturbations and matter effects are all just entries in that matrix; the
   bundled Hamiltonians in :mod:`hamiltonians2nu` and :mod:`hamiltonians3nu`
   are examples, not limitations.

#. **You need the evolution operator, not only the probabilities.**
   :func:`oscprob3nu.evolution_operator_3nu` returns :math:`U_3(L)` itself, so
   it can be composed across segments or used to propagate a density matrix.

#. **You are scanning, not evaluating one point.**  Every core routine takes a
   stack of Hamiltonians, an array of baselines, or both broadcast against
   each other, and returns the whole scan in one call.

What it is not
--------------

* **Not a solver for continuously varying Hamiltonians.**  See the first point
  above, and :doc:`recipes` for where the boundary lies in practice.
* **Not a four-flavor code.**  The SU(2) and SU(3) expansions are specific to
  two and three flavors; a sterile fourth is outside what they cover.
* **Not a flux, cross-section or detector code.**  It computes oscillation
  probabilities and stops there.
* **Not a fitting framework.**  There is no likelihood machinery; the
  probabilities are meant to be handed to whatever does that.

Performance
-----------

A single probability takes about 8 microseconds for three flavors and 1
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
     - 38 ms
     - 1.8 ms
     - **0.31 ms**
   * - Three flavors, versus energy
     - 34 ms
     - 1.5 ms
     - **0.20 ms**
   * - Oscillogram, 100 x 100
     - 197 ms
     - 5.3 ms
     - **0.85 ms**
   * - Two flavors, versus baseline
     - 6.9 ms
     - **0.07 ms**
     - not used

Best of seven runs, interleaved, on one machine.  Repeated runs vary by tens
of per cent, so read these as orders of magnitude rather than constants ---
and if the exact figures matter to you, `notebook 09
<https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/09_performance.ipynb>`_
measures them on the machine that runs it.

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

   pip install nuoscprobexact

or, for a clone with the notebooks, the worked examples and the test suite:

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
   recipes
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
