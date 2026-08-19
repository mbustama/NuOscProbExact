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
   * - Compute probabilities (:mod:`oscprob2nu`, :mod:`oscprob3nu`,
       :mod:`oscprob4nu`)
     - ``numpy``
   * - Use the bundled sample Hamiltonians
     - ``numpy``
   * - Run the notebooks (``notebooks/``)
     - ``numpy``, ``matplotlib``, Jupyter
   * - Run the regression suite (``tests/``)
     - ``numpy``, ``pytest``, ``scipy``, ``coverage``
   * - Speed up large scans (installed by default)
     - ``numba``
   * - Build this documentation
     - the contents of ``docs/requirements.txt``

``scipy`` is used only by the test suite, to check the evolution operator
against an independent matrix exponential.  The library itself never imports
it.

.. important::

   ``numba`` is a base dependency as of 1.13.0, not an optional extra, so the
   compiled path is what a plain install gets.  It brings a ceiling with it:
   ``numba`` declares an upper bound on ``numpy`` (0.66.0 requires
   ``numpy<2.5``), and installing into an environment holding a newer ``numpy``
   than that will downgrade it.  The bound is ``numba``'s own and cannot be
   avoided from this side.  If the newest ``numpy`` matters more than the speed,
   install with ``--no-deps`` and add ``numpy`` yourself, or leave ``numba``
   installed and set :data:`fastkernels.USE_NUMBA` to ``False``.

**NuOscProbExact** requires Python 3.9 or newer, and every release is tested
on 3.9, 3.10, 3.11, 3.12, and 3.13.  The floor comes from
:func:`numpy.broadcast_shapes`, which the batched paths use and which arrived
in NumPy 1.20; 3.9 is also the oldest version for which ``numba`` still has a
wheel, which is why it is not excluded by a version marker --- every Python
from 3.9 to 3.14 resolves a working ``numba``.

Installing
----------

From PyPI (recommended)
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

   pip install nuoscprobexact

That is the whole installation, and it includes the compiled backend:
``numpy`` and ``numba`` both arrive, and large scans run on compiled kernels
without anything further being asked for.

The optional extras add what each task needs, and can be combined:

.. code-block:: shell

   pip install "nuoscprobexact[fast]"       # numba; a base dependency, so a no-op
   pip install "nuoscprobexact[notebooks]"  # Jupyter, matplotlib and scipy
   pip install "nuoscprobexact[test]"       # pytest, scipy and coverage
   pip install "nuoscprobexact[docs]"       # Sphinx and friends

This puts ``oscprob2nu``, ``oscprob3nu``, ``oscprob4nu``, ``hamiltonians2nu``,
``hamiltonians3nu``, ``hamiltonians4nu``, ``globaldefs``, ``fastkernels``,
``slabs`` and ``earth`` on your Python path under exactly those names --- the
same names the paper and the worked examples use.

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

The three core modules are self-contained --- they need only ``numpy`` and
the standard library --- so copying ``src/oscprob2nu.py``,
``src/oscprob3nu.py`` or ``src/oscprob4nu.py`` into your own project is a
supported way to use **NuOscProbExact**.  Each imports :mod:`fastkernels`
if it is there and does without it if it is not, so a lone copy works and
simply runs the NumPy path; a test copies each of the three out and
exercises it that way.  Adding
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
than a sign of trouble: the only tests that skip are those for the Numba
backend, which stand down when ``numba`` is not installed.  Since ``numba``
is a base dependency, that now happens only if it has been uninstalled
deliberately, so the usual result is the whole suite running.  It is several
hundred tests and takes a few seconds.

The suite checks the SU(2), SU(3) and SU(4) machinery against independent
computations: unitarity of the evolution operator, agreement with
``scipy.linalg.expm``, agreement with the standard oscillation formulas, the
``d`` tensor and star product of the SU(3) algebra, and the sign conventions
of the sample Hamiltonians.  Beyond that it exercises the batched evaluation
paths against the scalar ones, the two backends against each other, degenerate
and near-degenerate Hamiltonians, the energies and baselines the library is
actually used at, and the agreement between the type annotations and the
docstrings --- and it runs every example embedded in the docstrings.

Then run one of the worked examples:

.. code-block:: shell

   cd examples
   python example_3nu_vacuum.py

.. note::

   This directory was called ``test/`` in version 1.0.0 of the code, and is
   named that way in version 2 of `the paper
   <https://arxiv.org/abs/1904.12391>`_.  It was renamed to ``examples/`` to
   stop it being confused with ``tests/``, which holds the regression suite.

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
   ├── RELEASE_NOTES.md                 # What changed since the published version
   ├── LICENSE                          # MIT license
   ├── README.md                        # Project overview and worked examples
   ├── pyproject.toml                   # Packaging metadata and pytest configuration
   ├── MANIFEST.in                      # Keeps the test suite out of the source distribution
   ├── examples/                        # Runnable scripts, one per scenario, linked from README.md
   │   ├── example_2nu_trivial.py       # Two-flavor, arbitrary Hamiltonian
   │   ├── example_2nu_vacuum.py        # Two-flavor, oscillations in vacuum
   │   ├── example_2nu_vacuum_coeffs.py # Two-flavor, expansion coefficients
   │   ├── example_3nu_trivial.py       # Three-flavor, arbitrary Hamiltonian
   │   ├── example_3nu_vacuum.py        # Three-flavor, oscillations in vacuum
   │   ├── example_3nu_vacuum_coeffs.py # Three-flavor, expansion coefficients
   │   ├── example_3nu_matter.py        # Three-flavor, oscillations in matter
   │   ├── example_3nu_nsi.py           # Three-flavor, matter with NSI
   │   └── example_3nu_liv.py           # Three-flavor, LIV background
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
   │       ├── notebooks.rst            # The twenty tutorial notebooks, in reading order
   │       ├── methodology.rst          # The SU(2), SU(3) and SU(4) expansions
   │       ├── functions.rst            # API reference, from the docstrings
   │       ├── references.rst           # Bibliography
   │       ├── refs.bib                 # BibTeX entries for the bibliography
   │       ├── changelog.rst            # Includes the root CHANGELOG.md
   │       └── _static/
   │           ├── nuoscprobexact_logo.png
   │           └── slabs_composition.svg  # How slabs compose, drawn for quickstart.rst
   ├── img/                             # Figures from earlier versions of README.md
   │   ├── anim_cp.gif                  # The demonstration reel, one scene per file
   │   ├── anim_earth.gif
   │   ├── anim_slabs.gif
   │   ├── anim_sterile.gif
   │   ├── prob_3nu_vacuum_vs_baseline_ee_em_et.png
   │   ├── prob_3nu_vacuum_vs_energy_ee_em_et.png
   │   └── gallery/                     # Figures lifted from the notebooks, shown in README.md
   │       ├── gallery_anim_cp.png
   │       ├── gallery_anim_earth.png
   │       ├── gallery_anim_slabs.png
   │       ├── gallery_anim_sterile.png
   │       ├── gallery_biprobability.png
   │       ├── gallery_earth.png
   │       ├── gallery_long_range.png
   │       ├── gallery_matter.png
   │       ├── gallery_ordering.png
   │       ├── gallery_oscillogram.png
   │       ├── gallery_prem.png
   │       ├── gallery_profiles.png
   │       ├── gallery_sterile.png
   │       ├── gallery_sterile_earth.png
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
   │   ├── 16_four_neutrinos.ipynb      # A 3+1 sterile state, through the SU(4) expansion
   │   ├── 17_cross_checks.ipynb        # Corroboration from nuSQuIDS and Zaglauer-Schwarzer
   │   ├── 18_evolution_operator.ipynb  # The operator, and the SU(n) coefficients
   │   ├── 19_animations.ipynb          # Four animated scenes, as stills; the reel they came from
   │   ├── 20_arbitrary_hamiltonian.ipynb  # A long-range force, carried through three profiles
   │   ├── make_notebooks.py            # Generates and executes all of the above
   │   └── coastlines_crude.json        # Coastlines for the Earth cutaway, from GSHHS
   ├── src/                             # The library
   │   ├── oscprob2nu.py                # Two-flavor probabilities, SU(2) expansion
   │   ├── oscprob3nu.py                # Three-flavor probabilities, SU(3) expansion
   │   ├── oscprob4nu.py                # Four-flavor probabilities, SU(4) expansion
   │   ├── hamiltonians2nu.py           # Example two-flavor Hamiltonians
   │   ├── hamiltonians3nu.py           # Example three-flavor Hamiltonians
   │   ├── hamiltonians4nu.py           # Example four-flavor (3+1) Hamiltonians
   │   ├── globaldefs.py                # Physical constants and unit conversions
   │   ├── fastkernels.py               # Numba kernels, with a NumPy fallback
   │   ├── slabs.py                     # Propagation across adjacent slabs
   │   └── earth.py                     # PREM, chord geometry, and Earth crossings
   ├── tools/                           # Scripts that are not part of the package
   │   └── make_demo_video.py           # Joins and shrinks the clips notebook 19 renders
   ├── resources/                       # Things the repository carries but does not install
   │   └── paper/                       # Source of the paper that documents this library
   │       ├── README.md                # How to build it, and what each file is
   │       ├── main.tex                 # The paper; BibTeX, no inlined bibliography
   │       ├── make_versions.py         # Derives the clean and the diffed versions
   │       ├── baseline_cpc_v1.tex      # The published version, for the diff
   │       ├── refs.bib                 # The bibliography, read by BibTeX at build
   │       ├── elsarticle.cls
   │       ├── elsarticle-num.bst
   │       └── figs/                    # Every figure, fourteen of them from notebook 10
   │           ├── validation.pdf
   │           ├── prob_2nu_vs_energy_compare.pdf
   │           ├── prob_3nu_vs_energy_compare.pdf
   │           ├── slabs_composition.pdf
   │           ├── earth_oscillogram.pdf
   │           ├── density_arrangement.pdf
   │           ├── sterile_earth_oscillogram.pdf
   │           ├── architecture.pdf
   │           ├── performance.pdf
   │           ├── exact_vs_approximations.pdf
   │           ├── speed_accuracy.pdf
   │           ├── prem_speed_accuracy.pdf
   │           ├── prem_speed_accuracy_3plus1.pdf
   │           ├── lri_earth.pdf
   │           └── Notes_on_SU_n__probability_relations.pdf
   └── tests/                           # Regression suite, run with pytest
       ├── conftest.py                  # Shared fixtures and path setup
       ├── bench/                       # Fair-comparison benchmark pipeline
       │   ├── manifest.json            # Pinned versions, build profiles, thread policy, capabilities
       │   ├── OBJECTIONS.md            # Each objection, the measurement answering it, and its test
       │   ├── requirements.lock        # Exact Python versions the benchmark venv installs
       │   ├── build.sh                 # Clones, hash-verifies and builds all seven at their pins
       │   ├── bench.hpp                # The C++ harness: owns main(), every clock, the statistics
       │   ├── machine.py               # Environment capture, the canary, and the rejection rule
       │   ├── conversions.py           # The one place a physical conversion factor is computed
       │   ├── gen_conversions.py       # Emits conversions.h so no adapter carries a physical literal
       │   └── adapters/                # One adapter per code, physics only
       │       └── nufast_earth.cpp     # NuFast-Earth, batched over energy and zenith, invariants hoisted
       ├── test_bench_pipeline.py       # The benchmark pipeline's fairness invariants
       ├── gen_stiff_reference.py       # Regenerates the fifty-digit four-flavor oracle
       ├── stiff_reference.json         # That oracle, frozen: nine Hamiltonians in hexadecimal floats
       ├── test_su3_algebra.py          # d tensor, star product, SU(3) invariants
       ├── test_oscprob4nu.py           # SU(4) algebra, quartic roots, 3+1 physics
       ├── test_evolution_operator.py   # U against an independent matrix exponential
       ├── test_probabilities.py        # Normalization, positivity, P = |U|^2
       ├── test_hamiltonians.py         # Sample Hamiltonians and sign conventions
       ├── test_reference_formulas.py   # Exact result against the standard formulas
       ├── test_matter_eigenvalues.py   # Matter spectrum, against Zaglauer-Schwarzer
       ├── test_edge_cases.py           # Degenerate and near-degenerate Hamiltonians
       ├── test_docstrings.py           # Runs the examples embedded in the docstrings
       ├── test_vectorized.py           # The batched path, against the scalar one
       ├── test_vectorized_hamiltonians.py  # Hamiltonians built for an array of energies
       ├── test_annotations.py          # Annotations, and their agreement with the docs
       ├── test_fastkernels.py          # Both backends, against each other
       ├── test_physical_scales.py      # Both backends at the scales actually used
       ├── test_slabs.py                # Slab composition, against expm
       ├── test_earth.py                # PREM, geometry, and Earth probabilities
       ├── test_batching_and_tolerance.py  # Energy batching, and the rtol/atol refinement
       ├── bit_capture.py               # Exact-bit capture, for refactors meant to be invisible
       ├── bit_compare.py               # Compares two captures, in ulps
       ├── test_documented_figures.py   # Keeps the quoted performance figures agreeing
       ├── test_paper.py                # Keeps the paper agreeing with the repository
       ├── test_version_consistency.py  # Keeps the version agreeing wherever it is implied
       ├── test_nusquids_comparison.py  # Against nuSQuIDS, an independent external code
       ├── nusquids_reference.py        # Regenerates the frozen nuSQuIDS reference data
       ├── nusquids_reference.json      # Those reference values, with their provenance
       ├── nusquids_scan.py             # Regenerates the frozen nuSQuIDS energy scan
       ├── nusquids_scan.json           # That scan, for the paper's figures
       ├── nufast_scan.json             # NuFast-LBL at two Newton settings
       ├── speed_accuracy.json          # The six-code constant-density speed-accuracy plane
       ├── timing_other_codes.json      # Timings behind the performance figure
       ├── prem_scan.py                 # Regenerates the two Earth speed-accuracy planes
       ├── prem_speed_accuracy.json     # Those planes, at three flavors and at 3+1
       ├── const_density_scan.json      # The exact reference for the comparison figure
       ├── external_drivers/            # Drivers for the codes that cannot be called from Python
       │   ├── README.md                # Every convention each one had to be told, and why
       │   ├── gen_prem_header.py       # Emits this library's PREM as a C header
       │   ├── nufast_drv.cpp           # NuFast-Earth on the PREM chord
       │   ├── nufast_earth_prem.txt
       │   ├── globes_drv.c             # GLoBES on the PREM chord
       │   ├── globes_prem.txt
       │   ├── prob3_drv.cpp            # Prob3++ on the PREM chord
       │   ├── prob3_prem.txt
       │   ├── perf_nusquids.py         # nuSQuIDS on the constant-density scan
       │   ├── globes_perf.c            # GLoBES on the same scan
       │   ├── prob3_perf.cpp           # Prob3++ on the same scan
       │   ├── nufast_lbl_perf.cpp      # NuFast-LBL at two Newton settings
       │   ├── const_scan.py            # The 50-digit reference for the comparison figure
       │   ├── globes_scan.c            # GLoBES on that comparison
       │   └── prob3_scan.cpp           # Prob3++ on that comparison
       └── test_file_tree.py            # Keeps this tree in step with the repository
