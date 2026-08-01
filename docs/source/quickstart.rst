Quickstart
==========

Every routine in **NuOscProbExact** takes a Hamiltonian as a plain nested
list or NumPy array and a baseline, and returns probabilities.  There is no
object to construct and no state to configure.

Units
-----

The two core modules are unit-agnostic: they require only that the
Hamiltonian and the baseline be given in reciprocal units, so that the
product :math:`H L` is dimensionless.

Everywhere else in **NuOscProbExact** the convention is:

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Quantity
     - Units
   * - Hamiltonian
     - eV
   * - Mass-squared differences
     - eV\ :sup:`2`
   * - Neutrino energy
     - eV
   * - Baseline
     - eV\ :sup:`-1`
   * - Matter potential
     - eV
   * - Mixing angles, CP phases
     - radian

:mod:`globaldefs` provides ``CONV_KM_TO_INV_EV`` to convert a baseline in km
into eV\ :sup:`-1`.

An arbitrary Hamiltonian
------------------------

The shortest possible use: hand the code a Hermitian matrix and a baseline.

.. jupyter-execute::

   import oscprob3nu

   hamiltonian = [
       [1.0+0.0j, 0.0+2.0j, 0.0-1.0j],
       [0.0-2.0j, 3.0+0.0j, 3.0+0.0j],
       [0.0+1.0j, 3.0-0.0j, -5.0+0.0j],
   ]

   Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = \
       oscprob3nu.probabilities_3nu(hamiltonian, 1.0)

   print('Pee = %.5f, Pem = %.5f, Pet = %.5f' % (Pee, Pem, Pet))

The probabilities are ordered with the initial flavor varying slowest, so
``Pem`` is :math:`P(\nu_e \to \nu_\mu)`.

The Hamiltonian must be Hermitian.  Its trace is discarded, since it
contributes only an overall phase that cancels in the probabilities.

Oscillations in vacuum
----------------------

.. jupyter-execute::

   import numpy as np
   import oscprob3nu
   import hamiltonians3nu
   from globaldefs import *

   energy = 1.e9      # [eV]
   baseline = 1.3e3   # [km]

   h_vacuum_energy_indep = \
       hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
           S12_NO_BF, S23_NO_BF, S13_NO_BF, DCP_NO_BF, D21_NO_BF, D31_NO_BF)
   h_vacuum = np.multiply(1./energy, h_vacuum_energy_indep)

   Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = \
       oscprob3nu.probabilities_3nu(h_vacuum, baseline*CONV_KM_TO_INV_EV)

   print('Pee = %.5f, Pem = %.5f, Pet = %.5f' % (Pee, Pem, Pet))
   print('Pme = %.5f, Pmm = %.5f, Pmt = %.5f' % (Pme, Pmm, Pmt))
   print('Pte = %.5f, Ptm = %.5f, Ptt = %.5f' % (Pte, Ptm, Ptt))

``hamiltonian_3nu_vacuum_energy_independent`` deliberately leaves out the
factor :math:`1/E`, so that the energy-independent part can be computed once
and reused across an energy scan.

Matter, NSI, and LIV
--------------------

The remaining sample Hamiltonians all take that same energy-independent
vacuum term and add to it:

.. jupyter-execute::

   # Matter of constant density
   h_matter = hamiltonians3nu.hamiltonian_3nu_matter(
       h_vacuum_energy_indep, energy, VCC_EARTH_CRUST)

   # Non-standard interactions
   h_nsi = hamiltonians3nu.hamiltonian_3nu_nsi(
       h_vacuum_energy_indep, energy, VCC_EARTH_CRUST, EPS_3)

   # Lorentz invariance-violating background
   h_liv = hamiltonians3nu.hamiltonian_3nu_liv(
       h_vacuum_energy_indep, energy, SXI12, SXI23, SXI13, DXICP,
       B1, B2, B3, LAMBDA)

The matter potential ``VCC`` is positive for neutrinos.  Pass its negative
for antineutrinos; see :ref:`sign-convention` for why that is the whole
difference.

Two flavors
-----------

:mod:`oscprob2nu` mirrors :mod:`oscprob3nu` throughout:

.. jupyter-execute::

   import oscprob2nu
   import hamiltonians2nu

   # `sth` is sin(theta), not the angle itself
   h2_vacuum_energy_indep = \
       hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(
           S23_NO_BF, D31_NO_BF)
   h2_vacuum = np.multiply(1./energy, h2_vacuum_energy_indep)

   Pee, Pem, Pme, Pmm = \
       oscprob2nu.probabilities_2nu(h2_vacuum, baseline*CONV_KM_TO_INV_EV)

   print('Pee = %.5f, Pem = %.5f' % (Pee, Pem))

The evolution operator
----------------------

If you need the evolution operator itself --- to compose it across segments,
or to propagate a density matrix --- ask for it directly:

.. jupyter-execute::

   U3 = oscprob3nu.evolution_operator_3nu(h_vacuum,
                                          baseline*CONV_KM_TO_INV_EV)

   print('|U3[0][0]| = %.6f' % abs(U3[0][0]))

and the expansion coefficients, if you want those:

.. jupyter-execute::

   h_coeffs = oscprob3nu.hamiltonian_3nu_coefficients(h_vacuum)
   u_coeffs = oscprob3nu.evolution_operator_3nu_u_coefficients(
       h_vacuum, baseline*CONV_KM_TO_INV_EV)

   print('h has %d coefficients' % len(h_coeffs))

Because :math:`H` is time-independent, the evolution operator obeys the group
property :math:`U(L_1+L_2) = U(L_2) U(L_1)`, which is what makes composing
across segments of constant density legitimate.

.. _scanning:

Scanning: pass arrays
---------------------

Every routine above accepts a *stack* of Hamiltonians, a *stack* of
baselines, or both, and evaluates the whole stack at once.  This is
between one and two orders of magnitude faster than the same calls in a
Python loop, and it is the recommended way to produce a curve or a grid.

Versus baseline --- one Hamiltonian, many baselines:

.. jupyter-execute::

   baselines = np.linspace(1.0, 1.3e4, 2000)*CONV_KM_TO_INV_EV
   prob = oscprob3nu.probabilities_3nu(h_vacuum, baselines)

   Pem = prob[:, 1]           # same ordering as the scalar return
   print('shape:', prob.shape)

The characteristic equation depends only on the Hamiltonian, so it is
solved once for the whole scan rather than once per baseline.

Versus energy --- many Hamiltonians, one baseline.  Since
:math:`H \propto 1/E`, build the stack by dividing the
energy-independent term by an array of energies:

.. jupyter-execute::

   energies = np.logspace(-1, 1, 2000)*1.e9                  # [eV]
   h_stack = h_vacuum_energy_indep/energies[:, None, None]
   prob = oscprob3nu.probabilities_3nu(h_stack, baseline*CONV_KM_TO_INV_EV)

   print('shape:', prob.shape)

The sample Hamiltonians do this for you: pass an array of energies and
they return one Hamiltonian per energy, in exactly that shape.

.. jupyter-execute::

   h_stack = hamiltonians3nu.hamiltonian_3nu_matter(
       h_vacuum_energy_indep, energies, VCC_EARTH_CRUST)
   prob = oscprob3nu.probabilities_3nu(h_stack, baseline*CONV_KM_TO_INV_EV)

   print('shape:', prob.shape)

so a whole scan in matter, with NSI, or with LIV is two calls and no
Python loop.  The matter potential may be an array too, for a scan
across a density profile alongside the energy.

An oscillogram is the outer combination of the two: give the
Hamiltonians and the baselines separate axes and let them broadcast.

.. jupyter-execute::

   prob = oscprob3nu.probabilities_3nu(h_stack[:, None, :, :],
                                       baselines[:200][None, :])

   print('shape:', prob.shape)   # (energies, baselines, 9)

A single Hamiltonian and a scalar baseline still return a plain tuple,
exactly as before, so existing code is unaffected.

Going faster still
------------------

If `Numba <https://numba.pydata.org>`_ is installed, the batched paths are
evaluated by compiled kernels instead of NumPy, which is worth roughly
1.5x to 15x on large stacks, depending on the number of flavors:

.. code-block:: shell

   pip install "nuoscprobexact[fast]"

Nothing in your code changes; :mod:`fastkernels` is picked up
automatically, and the answers are the same to round-off.  If it is not
installed, the NumPy path is used and everything works as before.

It is used only where it is faster.  For three flavors that is every stack
size, by between two and sixteen times; for two flavors the NumPy path is
already lean enough to win below about fifty thousand elements, so it is
kept there.  :func:`fastkernels.worthwhile` makes that choice from measured
thresholds, so installing the extra can only help.

The first call compiles, which takes a few seconds.  The kernels are
cached on disk, so later runs start in milliseconds.  To force the NumPy
path --- to compare the two, say --- set

.. code-block:: python

   import fastkernels
   fastkernels.USE_NUMBA = False

The scalar path is deliberately left uncompiled: a single probability
takes about sixteen microseconds, which is not worth a compilation
pause.  Short *stacks* are also evaluated one element at a time, since
below about ten elements the array machinery costs more than it saves.

More examples
-------------

The ``examples/`` directory holds a runnable script for each of the cases
above.  They are the same examples given in the `README
<https://github.com/mbustama/NuOscProbExact/blob/main/README.md>`_, with the
output you should expect.

.. note::

   This directory was called ``test/`` in version 1.0.0 of the code, and is
   named that way in version 2 of `the paper
   <https://arxiv.org/abs/1904.12391>`_.  It was renamed to ``examples/`` to
   stop it being confused with ``tests/``, which holds the regression suite.
