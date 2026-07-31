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

.. code-block:: python

   import oscprob3nu

   hamiltonian = [
       [1.0+0.0j, 0.0+2.0j, 0.0-1.0j],
       [0.0-2.0j, 3.0+0.0j, 3.0+0.0j],
       [0.0+1.0j, 3.0-0.0j, -5.0+0.0j],
   ]

   Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = \
       oscprob3nu.probabilities_3nu(hamiltonian, 1.0)

The probabilities are ordered with the initial flavor varying slowest, so
``Pem`` is :math:`P(\nu_e \to \nu_\mu)`.

The Hamiltonian must be Hermitian.  Its trace is discarded, since it
contributes only an overall phase that cancels in the probabilities.

Oscillations in vacuum
----------------------

.. code-block:: python

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

which gives

.. code-block:: text

   Pee = 0.92768, Pem = 0.01432, Pet = 0.05800
   Pme = 0.04023, Pmm = 0.37887, Pmt = 0.58090
   Pte = 0.03210, Ptm = 0.60680, Ptt = 0.36110

``hamiltonian_3nu_vacuum_energy_independent`` deliberately leaves out the
factor :math:`1/E`, so that the energy-independent part can be computed once
and reused across an energy scan.

Matter, NSI, and LIV
--------------------

The remaining sample Hamiltonians all take that same energy-independent
vacuum term and add to it:

.. code-block:: python

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

.. code-block:: python

   import oscprob2nu
   import hamiltonians2nu

   h_vacuum_energy_indep = \
       hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(
           S23_NO_BF, D31_NO_BF)
   h_vacuum = np.multiply(1./energy, h_vacuum_energy_indep)

   Pee, Pem, Pme, Pmm = \
       oscprob2nu.probabilities_2nu(h_vacuum, baseline*CONV_KM_TO_INV_EV)

The evolution operator
----------------------

If you need the evolution operator itself --- to compose it across segments,
or to propagate a density matrix --- ask for it directly:

.. code-block:: python

   U3 = oscprob3nu.evolution_operator_3nu(h_vacuum,
                                          baseline*CONV_KM_TO_INV_EV)

and the expansion coefficients, if you want those:

.. code-block:: python

   h_coeffs = oscprob3nu.hamiltonian_3nu_coefficients(h_vacuum)
   u_coeffs = oscprob3nu.evolution_operator_3nu_u_coefficients(
       h_vacuum, baseline*CONV_KM_TO_INV_EV)

Because :math:`H` is time-independent, the evolution operator obeys the group
property :math:`U(L_1+L_2) = U(L_2) U(L_1)`, which is what makes composing
across segments of constant density legitimate.

Scanning: pass arrays
---------------------

Every routine above accepts a *stack* of Hamiltonians, a *stack* of
baselines, or both, and evaluates the whole stack at once.  This is
between one and two orders of magnitude faster than the same calls in a
Python loop, and it is the recommended way to produce a curve or a grid.

Versus baseline --- one Hamiltonian, many baselines:

.. code-block:: python

   baselines = np.linspace(1.0, 1.3e4, 2000)*CONV_KM_TO_INV_EV
   prob = oscprob3nu.probabilities_3nu(h_vacuum, baselines)
   prob.shape        # (2000, 9)
   Pem = prob[:, 1]  # same ordering as the scalar return

The characteristic equation depends only on the Hamiltonian, so it is
solved once for the whole scan rather than once per baseline.

Versus energy --- many Hamiltonians, one baseline.  Since
:math:`H \propto 1/E`, build the stack by dividing the
energy-independent term by an array of energies:

.. code-block:: python

   energies = np.logspace(-1, 1, 2000)*1.e9                  # [eV]
   h_stack = h_vacuum_energy_indep/energies[:, None, None]
   prob = oscprob3nu.probabilities_3nu(h_stack, baseline*CONV_KM_TO_INV_EV)

An oscillogram is the outer combination of the two: give the
Hamiltonians and the baselines separate axes and let them broadcast.

.. code-block:: python

   prob = oscprob3nu.probabilities_3nu(h_stack[:, None, :, :],
                                       baselines[None, :])
   prob.shape        # (len(energies), len(baselines), 9)

A single Hamiltonian and a scalar baseline still return a plain tuple,
exactly as before, so existing code is unaffected.

More examples
-------------

The ``test/`` directory holds a runnable script for each of the cases above.
They are the same examples given in the `README
<https://github.com/mbustama/NuOscProbExact/blob/main/README.md>`_, with the
output you should expect.
