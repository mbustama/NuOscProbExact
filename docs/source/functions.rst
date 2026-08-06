API reference
=============

Generated from the docstrings.  Every ``Examples`` block below is **executed
when this page is built**, so the results shown are what the code returns
rather than numbers written beside it.  The regression suite runs the same
blocks on every supported Python.

Core modules
------------

The three core modules are self-contained: they need only ``numpy``, and they
accept *any* Hermitian Hamiltonian.  Copying one into your own project is a
supported way to use **NuOscProbExact**.

All three are unit-agnostic: they require only that the Hamiltonian and the
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

oscprob4nu
^^^^^^^^^^

The :math:`n = 4` member, and the last one that exists in closed form; see
:ref:`why-four-is-the-end`.  Its accuracy on stiff 3+1 spectra deserves a
look before use: :ref:`stiff-spectra`.

.. automodule:: oscprob4nu
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

hamiltonians4nu
^^^^^^^^^^^^^^^

.. automodule:: hamiltonians4nu
   :members:
   :undoc-members:

Piecewise-constant matter
-------------------------

The exact expansions assume a Hamiltonian that does not change along the
trajectory.  A neutrino crossing the Earth does not have one, so its path is
cut into slabs, each solved exactly and the results multiplied.  Within a slab
nothing is approximated; the only approximation is how finely a continuously
varying profile is sliced, and that is an argument the caller controls.

slabs
^^^^^

.. automodule:: slabs
   :members:
   :undoc-members:

earth
^^^^^

.. automodule:: earth
   :members:
   :undoc-members:

The compiled backend
--------------------

.. automodule:: fastkernels
   :members:
   :undoc-members:

Constants
---------

Physical constants, unit-conversion factors, and the NuFit 4.0 best-fit
oscillation parameters :cite:`Esteban:2018azc` for both mass orderings.  The
core modules do not need any of these.  The sample non-standard-interaction
strengths are deliberately large, so that the worked examples show a visible
effect; the combination matter oscillations are sensitive to puts them in the
LMA-D region :cite:`Coloma:2019mbs` rather than anywhere near a fit.

globaldefs
^^^^^^^^^^

.. automodule:: globaldefs
   :members:
   :undoc-members:
