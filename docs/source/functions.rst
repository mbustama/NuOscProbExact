API reference
=============

Generated from the docstrings.  Every ``Examples`` block below is executable
and is run as a doctest by the regression suite, so the numbers shown are
what the code actually returns.

Core modules
------------

The two core modules are self-contained: they need only ``numpy``, and they
accept *any* Hermitian Hamiltonian.  Copying either one into your own project
is a supported way to use **NuOscProbExact**.

Both are unit-agnostic: they require only that the Hamiltonian and the
baseline be given in reciprocal units, so that :math:`H L` is dimensionless.

oscprob2nu
^^^^^^^^^^

.. automodule:: oscprob2nu
   :members:
   :undoc-members:

oscprob3nu
^^^^^^^^^^

.. automodule:: oscprob3nu
   :members:
   :undoc-members:

Sample Hamiltonians
-------------------

These build the Hamiltonian for a number of standard scenarios.  They are
examples, not a limit on what can be computed: the core routines accept any
Hermitian matrix.

Energies are in eV, mass-squared differences in eV\ :sup:`2`, baselines in
eV\ :sup:`-1`, and potentials in eV.  See :ref:`sign-convention` for the
sign the vacuum Hamiltonians adopt, and why it matters once a matter
potential is added.

hamiltonians2nu
^^^^^^^^^^^^^^^

.. automodule:: hamiltonians2nu
   :members:
   :undoc-members:

hamiltonians3nu
^^^^^^^^^^^^^^^

.. automodule:: hamiltonians3nu
   :members:
   :undoc-members:

Optional compiled backend
-------------------------

.. automodule:: fastkernels
   :members:
   :undoc-members:

Constants
---------

Physical constants, unit-conversion factors, and the NuFit 4.0 best-fit
oscillation parameters for both mass orderings.  The core modules do not need
any of these.

globaldefs
^^^^^^^^^^

.. automodule:: globaldefs
   :members:
   :undoc-members:
