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

**NuOscProbExact** computes exact two-, three- and four-flavor neutrino
oscillation probabilities for arbitrary time-independent Hamiltonians.  It is
a Python implementation of the method of Ohlsson and Snellman
:cite:`Ohlsson:1999xb`, revisited in :cite:`Bustamante:2019ggq`.

The method expands the Hamiltonian and the time-evolution operator in the
bases of SU(2), SU(3) and SU(4) matrices, which yields concise, analytical,
and **exact** closed-form expressions for the probabilities.

Four flavors brings 3+1 sterile scenarios into scope as *closed* four-state
systems, rather than as a leak out of the three-flavor block.  It is also
where the method ends: :ref:`why-four-is-the-end` explains why there is no
five-flavor counterpart, and why that is a fact about algebra rather than a
gap in the code.

.. _what-exact-means:

What "exact" means here
------------------------

There is no approximation beyond floating-point round-off.  The evolution
operator produced by the SU(2), SU(3) and SU(4) expansions is the matrix
exponential :math:`e^{-i H_0 L}`, and the test suite checks that it is,
against an independent computation:

.. list-table::
   :header-rows: 1
   :widths: 60 40

   * - Property
     - Measured agreement
   * - :math:`U_3` vs ``scipy.linalg.expm``, 100 random Hermitian Hamiltonians
     - 8e-15
   * - :math:`U_2` vs ``scipy.linalg.expm``, 100 random Hermitian Hamiltonians
     - 1.6e-15
   * - :math:`U_4` vs ``scipy.linalg.expm``, 400 random Hermitian Hamiltonians
     - 2.9e-14
   * - Unitarity, :math:`U^\dagger U - \mathbb{1}`
     - 5e-15
   * - Hard-coded :math:`d_{ijk}` vs :math:`\frac{1}{4}\mathrm{Tr}(\{\lambda_i,\lambda_j\}\lambda_k)`, all 512 entries
     - 2e-16
   * - Exact result vs the standard vacuum oscillation formula
     - 7e-16 (two flavors), 1e-14 (three)
   * - Four flavors with the sterile angles off, vs :mod:`oscprob3nu`
     - 3e-14, in vacuum and in matter alike
   * - Three flavors vs `nuSQuIDS <https://github.com/arguelles/nuSQuIDS>`_,
       an independent external code
     - 2e-15
   * - Four flavors vs nuSQuIDS
     - 4e-16 to 3e-10
   * - Matter spectrum vs the Zaglauer-Schwarzer closed form
     - 7e-16

Each figure is the worst case the corresponding test actually reaches, on
the conditions that test uses --- not the tolerance it asserts against,
which is looser.  ``tests/test_documented_figures.py`` guards the
performance numbers on this page; these accuracy ones are reproduced by
running the suite.

One line of context for the four-flavor row.  A *stiff* spectrum --- a 3+1
scenario with an eV-scale :math:`\Delta m^2_{41}` --- reaches about
:math:`10^{-9}` rather than :math:`10^{-14}` in probabilities, limited by
what the SU(4) invariants retain rather than by the expansion itself.  Since
1.13.0 the invariants are formed in double-double arithmetic, which takes the
*roots* to :math:`3.6 \times 10^{-17}`, under a fifth of a ``float64`` ulp; the
probability figure improves less, because rebuilding :math:`U_4` takes second
differences of :math:`e^{-i\psi L}` and so amplifies whatever root error
remains by the square of the accumulated phase.  Both figures are far below
anything an experiment can resolve, so this is a statement about the
*exactness claim* and about error accumulating over composed slabs, not about
whether the probabilities are good enough to use.  :ref:`stiff-spectra` gives
the mechanism, what the code does about it, and what the alternatives
measured.

The last three rows are the ones an internal suite cannot supply.  Everything
above them compares this library against itself or against a formula
transcribed from the same papers, so a convention that were wrong
*consistently* --- a mixing-matrix ordering, a sign, a unit --- would pass all
of it.  The external checks would not, and they cover the antineutrino rule
and the mass ordering explicitly.  `Notebook 17
<https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/17_cross_checks.ipynb>`_
works through both, including the two conventions that have to be matched
first and the one residual that turns out to be ours rather than a
disagreement.

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
   routines take any Hermitian :math:`2\times2`, :math:`3\times3` or
   :math:`4\times4` matrix.  Non-standard interactions, Lorentz-invariance
   violation, sterile states and matter effects are all just entries in that
   matrix; the bundled Hamiltonians in :mod:`hamiltonians2nu`,
   :mod:`hamiltonians3nu` and :mod:`hamiltonians4nu` are examples, not
   limitations.

#. **You need the evolution operator, not only the probabilities.**
   :func:`oscprob3nu.evolution_operator_3nu` returns :math:`U_3(L)` itself, so
   it can be composed across segments or used to propagate a density matrix.

#. **You are scanning, not evaluating one point.**  Every core routine takes a
   stack of Hamiltonians, an array of baselines, or both broadcast against
   each other, and returns the whole scan in one call.

What it is not
--------------

* **Not a solver for continuously varying Hamiltonians.**  See
  :ref:`use-magnus-instead` just below, which says plainly when to reach for
  a different tool.
* **Not a five-flavor code.**  The expansions run to SU(4) and stop there,
  because the closed form does: see :ref:`why-four-is-the-end`.  Four flavors
  covers 3+1, which is the case people actually ask for.
* **Not a flux, cross-section or detector code.**  It computes oscillation
  probabilities and stops there.
* **Not a fitting framework.**  There is no likelihood machinery; the
  probabilities are meant to be handed to whatever does that.

.. _use-magnus-instead:

When to use Magnus instead
--------------------------

**NuOscProbExact** assumes the Hamiltonian is constant, or piecewise
constant.  Everything it is good at follows from that, and so does the one
case where it is the wrong tool.

Reach for **Magnus** instead when **the Hamiltonian varies continuously and
appreciably over an oscillation length**.
A smoothly varying profile can always be approximated by slabs, but the step
size is then set by the oscillation rather than by the density, and the slab
count grows until the calculation is neither exact nor quick.  Concretely:

.. list-table::
   :header-rows: 1
   :widths: 46 27 27

   * - Situation
     - Use this
     - Because
   * - Constant density
     - **NuOscProbExact**
     - One closed form, no integration
   * - Piecewise constant, tens of layers --- the Earth through PREM
     - **NuOscProbExact** (:mod:`slabs`, :mod:`earth`)
     - Each layer solved exactly, operators multiplied
   * - Smoothly varying, slow against the oscillation
     - Either
     - Slabbing converges quickly
   * - Smoothly varying, fast against the oscillation --- the Sun, adiabatic MSW
     - **Magnus**
     - Slabbing needs :math:`\sim 10^4` steps per resonance crossing
   * - Genuinely open systems: decay, decoherence
     - Neither
     - Needs a Lindblad solver, not a unitary one

`Notebook 14
<https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/14_solar_and_adiabatic_msw.ipynb>`_
works the solar case through and shows exactly where the wall is, rather than
asserting it.

Performance
-----------

A single probability takes about 8 microseconds for three flavors and 1
for two.  Scans --- a curve versus baseline or energy, an oscillogram over
both --- are what the code mostly does, and two things make those much
faster without changing any answer.

**Pass arrays instead of looping.**  Every core routine takes a stack of
Hamiltonians, an array of baselines, or both, and evaluates them in one
call.  That is worth roughly 20 to 90 times the equivalent Python loop, and
needs no extra dependency: the expensive part of the expansion, the
characteristic equation whose roots give the oscillation phases, depends on
the Hamiltonian alone, so a scan over baselines solves it once instead of
once per point.  See :ref:`scanning` for the three patterns.

**Numba, which comes with the package.**  Since 1.13.0 it is a base
dependency rather than an extra, so the batched paths run as compiled loops
spread over the available cores with nothing further to install.  Where it is
unavailable, or switched off with :data:`fastkernels.USE_NUMBA`, the NumPy
path is used and the results are the same to round-off.  On 2000-point scans,
against the Python loop:

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

**Layered matter has its own kernel.**  :mod:`slabs` and :mod:`earth`
compose evolution operators rather than computing probabilities, so they
could not use the backend at all before 1.12.0 --- it offered probability
kernels only, and an Earth crossing quietly ran the NumPy path however the
extra was installed.  They now compose a whole trajectory in one compiled
pass, never materialising the stack of operators and never dispatching a
product per slab, which is worth ~12x, ~12x and ~9x on a single Earth
crossing at two, three and four flavors.  The threshold there is one slab
rather than the fifty thousand above: the comparison is against a Python
loop, so the kernel is ahead at every length.

**An Earth scan is one call, not a loop.**  Both the energy and the zenith
angle may be arrays, and they broadcast, so an oscillogram is
``probabilities_3nu_earth(h, energies[None, :], costhz[:, None])`` rather
than a loop over its points.  The geometry, the matter potentials and the
slab widths depend on the angle alone, so a scan builds them once instead
of once per energy; the chord kernels then build each slab's Hamiltonian
as they go, so the stack of operators that made the batched path
memory-bound is never allocated at all.  Over 2000 energies:

.. list-table::
   :header-rows: 1
   :widths: 34 22 22 22

   * - Earth scan, 2000 energies
     - NumPy loop
     - Compiled loop
     - Compiled array
   * - Two flavors
     - 542 ms
     - 44 ms
     - **1.2 ms**
   * - Three flavors
     - 895 ms
     - 77 ms
     - **7.5 ms**
   * - Four flavors
     - 2224 ms
     - 272 ms
     - **25 ms**

That is ~38x, ~10x and ~11x over the compiled loop, and ~460x, ~120x and
~87x over the NumPy one.  A 100 x 100 oscillogram at three flavors takes
42 ms against 371 ms for the loop, about 4 microseconds a point.

**A chord is a palindrome.**  A neutrino crossing a spherically symmetric
Earth meets every radius on the way in and again on the way out, so slab
:math:`j` and slab :math:`n-1-j` carry the same Hamiltonian, the same
width and therefore the same operator.  Composing from the centre outwards
computes each of them once instead of twice, which is worth ~1.4x, ~1.5x
and ~1.8x of the figures above.  Nothing about it is specific to PREM:
:func:`fastkernels.palindromic` decides, and :mod:`slabs` asks it of any
sequence handed to it, so a symmetric castle wall is composed at the same
discount.  :data:`fastkernels.USE_PALINDROME` switches it off, which is
how to ask for the plain left-to-right product when a comparison needs
one.

**Ask for an accuracy, not a slab count.**  ``n_slabs_per_segment`` fixes
the discretisation rather than the error, and the two are not the same:
the error is strongly energy- and angle-dependent.  Pass ``rtol`` or
``atol`` to any :mod:`earth` entry point and the chord is refined until
the measured error meets it, on every returned probability;
:func:`earth.slabs_for_tolerance` does the search on its own, which is the
cheap way to run a scan.  See :doc:`recipes`.

**One cost runs the other way.**  Every entry point verifies that the
Hamiltonian is Hermitian, because one that is not returns probabilities that
still sum to one and so betrays nothing.  Validating a stack is a pass over
it --- the same order of work as evaluating it, and the compiled kernel has
made evaluating it fast --- so it costs 1.3 to 1.8 times on a 2000-point scan
and 3.2 to 5.7 times on a 200 000-point one.  It is on by default anyway;
where the Hamiltonians come from a construction already trusted, set
:data:`oscprob3nu.CHECK_HERMITICITY` to ``False``, and likewise on the other
two modules.

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
   notebooks
   methodology
   functions
   references
   changelog

Author
------

**NuOscProbExact** was written by Mauricio Bustamante
(mbustamante@gmail.com).  Bug reports and questions are best raised as
`GitHub issues <https://github.com/mbustama/NuOscProbExact/issues>`_, which
leave a public record others can find.

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
