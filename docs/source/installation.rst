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
   * - Run the notebooks (``notebooks/``)
     - ``numpy``, ``matplotlib``, Jupyter
   * - Run the regression suite (``tests/``)
     - ``numpy``, ``pytest``, ``scipy``
   * - Speed up large scans (optional)
     - ``numba``
   * - Build this documentation
     - the contents of ``docs/requirements.txt``

``scipy`` is used only by the test suite, to check the evolution operator
against an independent matrix exponential.  The library itself never imports
it.

**NuOscProbExact** requires Python 3.9 or newer, and every release is tested
on 3.9, 3.10, 3.11, 3.12, and 3.13.  The floor comes from
:func:`numpy.broadcast_shapes`, which the batched paths use and which arrived
in NumPy 1.20; 3.9 is also the oldest version for which the optional ``numba``
backend still has a wheel.

Installing
----------

From PyPI (recommended)
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

   pip install nuoscprobexact

That is the whole installation.  The only required dependency is ``numpy``.

The optional extras add what each task needs, and can be combined:

.. code-block:: shell

   pip install "nuoscprobexact[fast]"       # numba, for the compiled kernels
   pip install "nuoscprobexact[notebooks]"  # Jupyter, matplotlib and scipy
   pip install "nuoscprobexact[test]"       # pytest and scipy
   pip install "nuoscprobexact[docs]"       # Sphinx and friends

This puts ``oscprob2nu``, ``oscprob3nu``, ``hamiltonians2nu``,
``hamiltonians3nu``, ``globaldefs``, ``fastkernels``, ``slabs`` and ``earth``
on your Python path under exactly those names --- the same names the paper and
the worked examples use.

From GitHub
^^^^^^^^^^^

Install from a clone if you want the notebooks, the worked examples from the
paper, the regression suite, or a version that is not yet released:

.. code-block:: shell

   git clone https://github.com/mbustama/NuOscProbExact.git
   cd NuOscProbExact
   pip install -e .

``-e`` installs in editable mode, so edits to ``src/`` take effect without
reinstalling.  The extras work the same way, for example
``pip install -e ".[fast,test]"``.

Without installing anything
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The two core modules are self-contained --- they need only ``numpy`` and
``cmath`` --- so copying ``src/oscprob2nu.py`` or ``src/oscprob3nu.py`` into
your own project is a supported way to use **NuOscProbExact**.  Adding
``src/`` to the path works too, and is what the bundled examples do:

.. code-block:: python

   import sys
   sys.path.append('/path/to/NuOscProbExact/src')

   import oscprob3nu

Checking the installation
-------------------------

Run the regression suite:

.. code-block:: shell

   pytest

Every test should either pass or be skipped, and skips are expected rather
than a sign of trouble: the only tests that skip are those for the optional
Numba backend, which stand down when ``numba`` is not installed --- the
default.  Install the ``fast`` extra and the whole suite runs.  It is a few
hundred tests and takes a few seconds.

The suite checks the SU(2) and SU(3) machinery against independent
computations: unitarity of the evolution operator, agreement with
``scipy.linalg.expm``, agreement with the standard oscillation formulas, the
``d`` tensor and star product of the SU(3) algebra, and the sign conventions
of the sample Hamiltonians.  Beyond that it exercises the batched evaluation
paths against the scalar ones, the two backends against each other, degenerate
and near-degenerate Hamiltonians, the energies and baselines the library is
actually used at, and the agreement between the type annotations and the
docstrings --- and it runs every example embedded in the docstrings as a
doctest.

Then run one of the worked examples:

.. code-block:: shell

   cd test
   python example_3nu_vacuum.py

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
   ├── .github/                         # Continuous integration (GitHub Actions)
   │   └── workflows/
   │       ├── tests.yml                # The suite: five Pythons, all three backends
   │       ├── lint.yml                 # ruff, and the docs build under -W
   │       ├── pages.yml                # Builds and deploys the docs to GitHub Pages
   │       └── publish.yml              # Publishes to PyPI on a GitHub Release
   ├── .gitignore                       # Build, cache, and generated-output artefacts
   ├── CHANGELOG.md                     # Notable changes, rendered as a docs page
   ├── LICENSE                          # MIT license
   ├── README.md                        # The file that you are reading
   ├── pyproject.toml                   # Packaging metadata and pytest configuration
   ├── docs/                            # Sphinx documentation
   │   ├── Makefile                     # `make html` on Linux and macOS
   │   ├── make.bat                     # `make html` on Windows
   │   ├── requirements.txt             # Documentation-only dependencies
   │   └── source/
   │       ├── conf.py                  # Sphinx configuration
   │       ├── index.rst                # Landing page
   │       ├── installation.rst         # Requirements, installation, file tree
   │       ├── quickstart.rst           # Shortest path to a probability
   │       ├── recipes.rst              # Numerical recipes, with pre-generated figures
   │       ├── methodology.rst          # The SU(2) and SU(3) expansions
   │       ├── functions.rst            # API reference, from the docstrings
   │       ├── references.rst           # Bibliography
   │       ├── refs.bib                 # BibTeX entries for the bibliography
   │       ├── changelog.rst            # Includes the root CHANGELOG.md
   │       └── _static/
   │           └── nuoscprobexact_logo.png
   ├── img/                             # Pre-computed figures shown in README.md
   │   ├── prob_3nu_vacuum_vs_baseline_ee_em_et.png
   │   ├── prob_3nu_vacuum_vs_energy_ee_em_et.png
   │   └── gallery/                     # Figures lifted from the notebooks, shown in README.md
   │       ├── gallery_biprobability.png
   │       ├── gallery_earth.png
   │       ├── gallery_matter.png
   │       ├── gallery_ordering.png
   │       ├── gallery_oscillogram.png
   │       ├── gallery_prem.png
   │       ├── gallery_profiles.png
   │       └── gallery_vacuum.png
   ├── notebooks/                       # Worked examples, with their figures stored inline
   │   ├── 01_basics.ipynb              # Units, one probability, and broadcasting
   │   ├── 02_vacuum_oscillations.ipynb # Against baseline and against energy
   │   ├── 03_matter_nsi_liv.ipynb      # Constant-density matter, NSI, and LIV
   │   ├── 04_oscillogram.ipynb         # Energy-baseline maps in one call
   │   ├── 05_biprobability.ipynb       # CP ellipses, in vacuum and in matter
   │   ├── 06_earth_and_prem.ipynb      # PREM, chord geometry, and slabs
   │   ├── 07_earth_probabilities.ipynb # Through the Earth, and between sites
   │   ├── 08_unusual_density_profiles.ipynb  # Castle-wall and other hand-built profiles
   │   ├── 09_performance.ipynb         # Looping vs broadcasting, and the backend
   │   ├── 10_paper_figures.ipynb       # The two figures from arXiv:1904.12391
   │   ├── 11_exact_vs_approximations.ipynb  # Where the textbook formulas break down
   │   ├── 12_ordering_and_octant.ipynb # Normal vs inverted, and the 23 octant
   │   ├── 13_antineutrinos.ipynb       # Conjugate and flip, and two ways to slip
   │   ├── 14_solar_and_adiabatic_msw.ipynb  # The MSW resonance, and the cost wall
   │   ├── 15_numerical_edge_cases.ipynb  # Degeneracies, and what does not go NaN
   │   └── make_notebooks.py            # Generates and executes all of the above
   ├── src/                             # The library
   │   ├── oscprob2nu.py                # Two-flavor probabilities, SU(2) expansion
   │   ├── oscprob3nu.py                # Three-flavor probabilities, SU(3) expansion
   │   ├── hamiltonians2nu.py           # Example two-flavor Hamiltonians
   │   ├── hamiltonians3nu.py           # Example three-flavor Hamiltonians
   │   ├── globaldefs.py                # Physical constants and unit conversions
   │   ├── fastkernels.py               # Optional Numba kernels, with a NumPy fallback
   │   ├── slabs.py                     # Propagation across adjacent slabs
   │   └── earth.py                     # PREM, chord geometry, and Earth crossings
   ├── test/                            # The worked examples from the paper
   │   ├── example_2nu_trivial.py       # Two-flavor, arbitrary Hamiltonian
   │   ├── example_2nu_vacuum.py        # Two-flavor, oscillations in vacuum
   │   ├── example_2nu_vacuum_coeffs.py # Two-flavor, expansion coefficients
   │   ├── example_3nu_trivial.py       # Three-flavor, arbitrary Hamiltonian
   │   ├── example_3nu_vacuum.py        # Three-flavor, oscillations in vacuum
   │   ├── example_3nu_vacuum_coeffs.py # Three-flavor, expansion coefficients
   │   ├── example_3nu_matter.py        # Three-flavor, oscillations in matter
   │   ├── example_3nu_nsi.py           # Three-flavor, matter with NSI
   │   └── example_3nu_liv.py           # Three-flavor, LIV background
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
       ├── test_physical_scales.py      # Both backends at the scales actually used
       ├── test_slabs.py                # Slab composition, against expm
       ├── test_earth.py                # PREM, geometry, and Earth probabilities
       └── test_file_tree.py            # Keeps this tree in step with the repository
