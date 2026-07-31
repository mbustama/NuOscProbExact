Installation
============

**NuOscProbExact** is written entirely in Python, so there is nothing to
compile or link.

Requirements
------------

The dependencies are deliberately minimal, and split by what you actually
want to do:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - To do this
     - You need
   * - Compute probabilities (:mod:`oscprob2nu`, :mod:`oscprob3nu`)
     - ``numpy``
   * - Use the bundled sample Hamiltonians
     - ``numpy``
   * - Generate the figures (``run_testsuite.py``, ``test/``)
     - ``numpy``, ``matplotlib``
   * - Run the regression suite (``tests/``)
     - ``numpy``, ``pytest``, ``scipy``
   * - Speed up large scans (optional)
     - ``numba``
   * - Build this documentation
     - the contents of ``docs/requirements.txt``

``scipy`` is used only by the test suite, to check the evolution operator
against an independent matrix exponential.  The library itself never imports
it.

**NuOscProbExact** requires Python 3.7 or newer.  It uses the ``@``
matrix-multiplication operator, so it will not run under Python 2.

Installing
----------

Clone the repository:

.. code-block:: shell

   git clone https://github.com/mbustama/NuOscProbExact.git
   cd NuOscProbExact

Then, optionally, install it so the modules are importable from anywhere:

.. code-block:: shell

   pip install -e .

This puts ``oscprob2nu``, ``oscprob3nu``, ``hamiltonians2nu``,
``hamiltonians3nu`` and ``globaldefs`` on your Python path under exactly those
names.  Skipping this step is fine: the bundled examples add ``../src`` to
``sys.path`` themselves, and you can do the same, or simply copy the two core
modules into your own project.

Optional extras install the dependencies for each task:

.. code-block:: shell

   pip install -e ".[fast]"    # numba, for the compiled batched kernels
   pip install -e ".[plots]"   # matplotlib, for the figures
   pip install -e ".[test]"    # pytest and scipy, for the regression suite
   pip install -e ".[docs]"    # sphinx and friends, for this documentation

Checking the installation
-------------------------

Run the regression suite:

.. code-block:: shell

   pytest

All 132 tests should pass.  They check the SU(2) and SU(3) machinery against
independent computations --- unitarity of the evolution operator, agreement
with ``scipy.linalg.expm``, agreement with the standard oscillation formulas,
and the sign conventions of the sample Hamiltonians --- and run every example
embedded in the docstrings.

Then run one of the worked examples:

.. code-block:: shell

   cd test
   python example_3nu_vacuum.py

and, if you have ``matplotlib``, the figure suite:

.. code-block:: shell

   python run_testsuite.py

which writes 42 figures to ``fig/``.

Building the documentation
--------------------------

.. code-block:: shell

   pip install -r docs/requirements.txt
   cd docs
   make html

The result lands in ``docs/build/html``; open ``index.html``.

File tree
---------

.. code-block:: text

   NuOscProbExact/
   ├── .gitignore                       # Build, cache, and generated-output artefacts
   ├── CHANGELOG.md                     # Notable changes, rendered as a docs page
   ├── LICENSE                          # MIT license
   ├── README.md                        # The file that you are reading
   ├── pyproject.toml                   # Packaging metadata and pytest configuration
   ├── run_testsuite.py                 # Generates the full suite of test figures
   ├── docs/                            # Sphinx documentation
   │   ├── Makefile                     # `make html` on Linux and macOS
   │   ├── make.bat                     # `make html` on Windows
   │   ├── requirements.txt             # Documentation-only dependencies
   │   └── source/
   │       ├── conf.py                  # Sphinx configuration
   │       ├── index.rst                # Landing page
   │       ├── installation.rst         # Requirements, installation, file tree
   │       ├── quickstart.rst           # Shortest path to a probability
   │       ├── methodology.rst          # The SU(2) and SU(3) expansions
   │       ├── functions.rst            # API reference, from the docstrings
   │       ├── references.rst           # Bibliography
   │       ├── refs.bib                 # BibTeX entries for the bibliography
   │       └── changelog.rst            # Includes the root CHANGELOG.md
   ├── fig/                             # Figures written by run_testsuite.py (initially empty)
   ├── img/                             # Pre-computed figures shown in README.md
   │   ├── prob_3nu_vacuum_vs_baseline_ee_em_et.png
   │   └── prob_3nu_vacuum_vs_energy_ee_em_et.png
   ├── src/                             # The library
   │   ├── oscprob2nu.py                # Two-flavor probabilities, SU(2) expansion
   │   ├── oscprob3nu.py                # Three-flavor probabilities, SU(3) expansion
   │   ├── hamiltonians2nu.py           # Example two-flavor Hamiltonians
   │   ├── hamiltonians3nu.py           # Example three-flavor Hamiltonians
   │   ├── globaldefs.py                # Physical constants and unit conversions
   │   └── fastkernels.py               # Optional Numba kernels, with a NumPy fallback
   ├── test/                            # Worked examples and figure generators
   │   ├── example_2nu_trivial.py       # Two-flavor, arbitrary Hamiltonian
   │   ├── example_2nu_vacuum.py        # Two-flavor, oscillations in vacuum
   │   ├── example_2nu_vacuum_coeffs.py # Two-flavor, expansion coefficients
   │   ├── example_3nu_trivial.py       # Three-flavor, arbitrary Hamiltonian
   │   ├── example_3nu_vacuum.py        # Three-flavor, oscillations in vacuum
   │   ├── example_3nu_vacuum_coeffs.py # Three-flavor, expansion coefficients
   │   ├── example_3nu_matter.py        # Three-flavor, oscillations in matter
   │   ├── example_3nu_nsi.py           # Three-flavor, matter with NSI
   │   ├── example_3nu_liv.py           # Three-flavor, LIV background
   │   ├── oscprob2nu_plot.py           # Two-flavor probability figures
   │   ├── oscprob3nu_plot.py           # Three-flavor probability figures
   │   ├── oscprob2nu_plotpaper.py      # The two-flavor figure from the paper
   │   └── oscprob3nu_plotpaper.py      # The three-flavor figure from the paper
   └── tests/                           # Regression suite, run with pytest
       ├── conftest.py                  # Shared fixtures and path setup
       ├── test_su3_algebra.py          # d tensor, star product, SU(3) invariants
       ├── test_evolution_operator.py   # U against an independent matrix exponential
       ├── test_probabilities.py        # Normalization, positivity, P = |U|^2
       ├── test_hamiltonians.py         # Sample Hamiltonians and sign conventions
       ├── test_reference_formulas.py   # Exact result against the standard formulas
       ├── test_edge_cases.py           # Degenerate and near-degenerate Hamiltonians
       ├── test_docstrings.py           # Runs the examples embedded in the docstrings
       ├── test_vectorized.py           # The batched path, against the scalar one
       ├── test_vectorized_hamiltonians.py  # Hamiltonians built for an array of energies
       ├── test_annotations.py          # Annotations, and their agreement with the docs
       ├── test_fastkernels.py          # Both backends, against each other
       └── test_file_tree.py            # Keeps this tree in step with the repository
