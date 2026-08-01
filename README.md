[![tests](https://github.com/mbustama/NuOscProbExact/actions/workflows/tests.yml/badge.svg)](https://github.com/mbustama/NuOscProbExact/actions/workflows/tests.yml)
[![Code Quality](https://github.com/mbustama/NuOscProbExact/actions/workflows/lint.yml/badge.svg)](https://github.com/mbustama/NuOscProbExact/actions/workflows/lint.yml)
[![codecov](https://codecov.io/gh/mbustama/NuOscProbExact/branch/main/graph/badge.svg)](https://codecov.io/gh/mbustama/NuOscProbExact)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue.svg)](https://mbustama.github.io/NuOscProbExact/)
[![PyPI](https://img.shields.io/pypi/v/nuoscprobexact.svg)](https://pypi.org/project/nuoscprobexact/)
[![Downloads](https://pepy.tech/badge/nuoscprobexact)](https://pepy.tech/project/nuoscprobexact)
[![arXiv](https://img.shields.io/badge/arXiv-1904.12391-orange.svg)](https://arxiv.org/abs/1904.12391)
[![DOI](https://zenodo.org/badge/182178323.svg)](https://zenodo.org/badge/latestdoi/182178323)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Code style: ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

# NuOscProbExact
Code to compute exact two- and three-neutrino oscillation probabilities using SU(2) and SU(3) expansions

> **Note:** The oscillation probabilities are computed exactly, with no approximation beyond floating-point round-off.  A regression test suite lives in `tests/` and can be run with `pytest`.

## What you can compute

Every figure below is produced by a notebook in [`notebooks/`](notebooks/), and the link under each one goes to the code that drew it.  The documentation collects the same material, with runnable snippets, on its [numerical recipes](https://mbustama.github.io/NuOscProbExact/recipes.html) page.

| | |
|:--:|:--:|
| <img src="img/gallery/gallery_vacuum.png" width="380"/><br/>**Oscillation probabilities** against baseline or energy, for two or three flavors.<br/>[notebook 02](notebooks/02_vacuum_oscillations.ipynb) | <img src="img/gallery/gallery_matter.png" width="380"/><br/>**Matter, NSI and Lorentz-invariance violation** — each just a different Hermitian matrix.<br/>[notebook 03](notebooks/03_matter_nsi_liv.ipynb) |
| <img src="img/gallery/gallery_oscillogram.png" width="380"/><br/>**Oscillograms** over energy and baseline: 57 600 probabilities in a single call.<br/>[notebook 04](notebooks/04_oscillogram.ipynb) | <img src="img/gallery/gallery_biprobability.png" width="380"/><br/>**CP violation**, as bi-probability ellipses in vacuum and in matter.<br/>[notebook 05](notebooks/05_biprobability.ipynb) |
| <img src="img/gallery/gallery_prem.png" width="380"/><br/>**The Earth's density**, from the Preliminary Reference Earth Model.<br/>[notebook 06](notebooks/06_earth_and_prem.ipynb) | <img src="img/gallery/gallery_earth.png" width="380"/><br/>**Neutrinos through the Earth**, in energy and zenith angle, or between two named sites.<br/>[notebook 07](notebooks/07_earth_probabilities.ipynb) |
| <img src="img/gallery/gallery_profiles.png" width="380"/><br/>**Arbitrary matter profiles** — castle walls and worse, exactly.<br/>[notebook 08](notebooks/08_unusual_density_profiles.ipynb) | <img src="img/gallery/gallery_ordering.png" width="380"/><br/>**Mass ordering and the θ₂₃ octant**, separated by matter through the Earth.<br/>[notebook 12](notebooks/12_ordering_and_octant.ipynb) |

## Contents

1. [What you can compute](#what-you-can-compute)

2. [What is NuOscProbExact?](#what-is-nuoscprobexact)

3. [Requirements](#requirements)

4. [Installation](#installation)

5. [Performance](#performance)

6. [Usage and examples](#usage-and-examples)
   1. [Basics](#basics)
   2. [Trivial example](#trivial-example)
   3. [Oscillations in vacuum: fixed energy and baseline](#oscillations-in-vacuum-fixed-energy-and-baseline)
   4. [Three-neutrino oscillations in vacuum: fixed energy, varying baseline](#three-neutrino-oscillations-in-vacuum-fixed-energy-varying-baseline)
   5. [Three-neutrino oscillations in vacuum: fixed baseline, varying energy](#three-neutrino-oscillations-in-vacuum-fixed-baseline-varying-energy)
   6. [Three-neutrino oscillations in matter](#three-neutrino-oscillations-in-matter)
   7. [Three-neutrino oscillations in matter with non-standard interactions (NSI)](#three-neutrino-oscillations-in-matter-with-non-standard-interactions-nsi)
   8. [Three-neutrino oscillations in a Lorentz invariance-violating (LIV) background](#three-neutrino-oscillations-in-a-lorentz-invariance-violating-liv-background)
   9. [Arbitrary Hamiltonians](#arbitrary-hamiltonians)

7. [Notebooks](#notebooks)

8. [Documentation and help](#documentation-and-help)

9. [Citing](#citing)


## What is NuOscProbExact?

**NuOscProbExact** is a Python implementation of the method developed by [Ohlsson & Snellman](https://arxiv.org/abs/hep-ph/9910546) to compute exact two-flavor and three-flavor neutrino oscillation probabilities for arbitrary time-independent Hamiltonians.  The method was revisited and the code presented in the paper *NuOscProbExact: a general-purpose code to compute exact two-flavor and three-flavor neutrino oscillation probabilities* ([arXiv:1904.12391](http://arxiv.org/abs/1904.12391)), by Mauricio Bustamante.

The method relies on expansions of the Hamiltonian and time-evolution operators in terms of SU(2) and SU(3) matrices in order to obtain concise, analytical, and exact expressions for the probabilities, that are also easy to implement and evaluate.  For details of the method, see the paper above.

### What it does

* **Exact probabilities for any Hermitian 2×2 or 3×3 Hamiltonian.**  There is no approximation beyond floating-point round-off.  Oscillations in vacuum, in matter, with non-standard interactions, in a Lorentz invariance-violating background and with sterile-like perturbations are not special cases in the code — each is a different matrix handed to the same routine.
* **The evolution operator itself**, not only the probabilities, so it can be composed across segments or used to propagate a density matrix.
* **Whole scans in one call.**  Every core routine accepts a stack of Hamiltonians, an array of baselines, or both broadcast against each other, which is tens of times faster than the equivalent Python loop and gives identical results.
* **Piecewise-constant matter.**  `slabs` propagates across a sequence of adjacent slabs of arbitrary width and density, solving each exactly and multiplying the operators.
* **The Earth.**  `earth` builds those slabs from the Preliminary Reference Earth Model, and computes probabilities along a given zenith angle or between two of fifteen predefined locations.
* **An optional compiled backend.**  With `numba` installed, large batched calls run on compiled kernels; without it the NumPy path is used and the answers are the same to round-off.

### What it does not do

* **Hamiltonians that vary continuously along the trajectory.**  The method assumes a Hamiltonian that does not change, and in exchange gives a closed form rather than a numerical integration.  A smoothly varying profile has to be approximated by slabs, and that is only practical while the profile varies slowly compared with the oscillation length.  It works well for the Earth; it is the wrong tool for the Sun, where the step size would be set by the oscillation rather than by the density, needing of order 10<sup>4</sup> slabs per crossing.  [Notebook 14](notebooks/14_solar_and_adiabatic_msw.ipynb) works that case through and shows where the wall is; for smoothly varying profiles use a Magnus-type method such as [Magnus](https://github.com/mbustama/Magnus).
* **More than three flavors.**  The SU(2) and SU(3) expansions are specific to two and three flavors.  A fourth, sterile flavor is outside what the closed forms cover.
* **Neutrino production, cross sections, fluxes or detector response.**  This computes oscillation probabilities and nothing downstream of them.
* **Fitting or statistics.**  There is no likelihood machinery here; the probabilities are meant to be handed to whatever does that.

**NuOscProbExact** was developed by Mauricio Bustamante.  If you use it in your work, please follow the directions on [Citing](#citing).


## Requirements

**NuOscProbExact** is fully written in Python 3.  It uses standard modules that are available, sometimes by default, as part of most Python installations, either stand-alone or via Anaconda.  Where a row names an extra, that is the one to install; the rest need nothing beyond `numpy`.  The commands are under [Installation](#installation) below, and are not repeated here.

| To do this | You need | Extra |
|---|---|---|
| Compute probabilities (`oscprob2nu.py`, `oscprob3nu.py`) | `numpy`, `cmath` | — |
| Use the bundled sample Hamiltonians (`hamiltonians2nu.py`, `hamiltonians3nu.py`) | `numpy`, `cmath`, `copy` | — |
| Propagate through layered matter or the Earth (`slabs.py`, `earth.py`) | `numpy` | — |
| Go faster on large scans *(optional)* | `numba` | `fast` |
| Run the notebooks (`notebooks/`) | `matplotlib`, Jupyter | `notebooks` |
| Run the regression suite (`tests/`) | `pytest`, `scipy` | `test` |
| Build the documentation | Sphinx and friends | `docs` |

Only `numpy` is ever required.  `scipy` is used by the test suite alone, to cross-check the evolution operator against an independent matrix exponential; the library itself never imports it.  `numba` is entirely optional — it is worth roughly 1.5x to 15x on large scans, depending on their size and the number of flavors, and without it the NumPy path is used and the results are identical to round-off.


## Installation

**NuOscProbExact** is pure Python: there is nothing to compile or link.

> **Python version:** The code requires Python 3.9 or newer, and every release is tested on 3.9, 3.10, 3.11, 3.12, and 3.13.  The floor comes from `numpy.broadcast_shapes`, which the batched paths use and which arrived in NumPy 1.20; 3.9 is also the oldest version for which the optional `numba` backend still has a wheel.

### From PyPI (recommended)

```shell
pip install nuoscprobexact
```

That is the whole installation.  The only required dependency is `numpy`.

The optional extras add what each task needs, and can be combined:

```shell
pip install "nuoscprobexact[fast]"       # numba, for the compiled batched kernels
pip install "nuoscprobexact[notebooks]"  # Jupyter, matplotlib and scipy, for notebooks/
pip install "nuoscprobexact[test]"       # pytest and scipy, to run the regression suite
pip install "nuoscprobexact[docs]"       # Sphinx and friends, to build the documentation
```

Then, in your own code:

```python
import numpy as np
import oscprob3nu
import hamiltonians3nu
import globaldefs as gd

h_vacuum = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
    gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, gd.DCP_NO_BF,
    gd.D21_NO_BF, gd.D31_NO_BF)

prob = oscprob3nu.probabilities_3nu(
    np.asarray(h_vacuum)/1.e9, 1300.0*gd.CONV_KM_TO_INV_EV)
```

The modules are installed under their bare names --- `oscprob2nu`, `oscprob3nu`, `hamiltonians2nu`, `hamiltonians3nu`, `globaldefs`, `fastkernels`, `slabs`, `earth` --- which is the same way the paper and the worked examples refer to them.

### From GitHub

Install from a clone if you want the notebooks, the worked examples from the paper, the regression suite, or a version that is not yet released:

```shell
git clone https://github.com/mbustama/NuOscProbExact.git
cd NuOscProbExact
pip install -e .
```

`-e` installs in editable mode, so edits to `src/` take effect without reinstalling.  The extras work the same way, for example `pip install -e ".[fast,test]"`.

A clone gives you the following file structure:

```text
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
├── examples/                        # Runnable scripts, the ones the README walks through
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
```

### Without installing anything

The two core modules are self-contained --- they need only `numpy` and `cmath` --- so copying `src/oscprob2nu.py` or `src/oscprob3nu.py` into your own project is a supported way to use **NuOscProbExact**.  Adding `src/` to the path works too, and is what the bundled examples do:

```python
import sys
sys.path.append('/path/to/NuOscProbExact/src')

import oscprob3nu
```

### Checking the installation

**Run the worked examples.**
   Inside the directory `examples/`, we provide several example files to get you started.  We also elaborate on these examples later in this README, and show the output thay you should expect from them.  To run any of the examples, just execute, *e.g.*,
   ```shell
   python example_2nu_trivial.py
   ```
   Inspecting the example files and reading their description below will help you to learn how to use **NuOscProbExact** in your own project.

   > **Renamed:** this directory was called `test/` in version 1.0.0 of the code, and is named that way in version 2 of [the paper](https://arxiv.org/abs/1904.12391).  It became `examples/` to stop it being confused with `tests/`, which holds the regression suite.

**Run the regression tests.**
   ```shell
   cd /home/MyProjects/NuOscProbExact
   pytest
   ```
   These check the SU(2) and SU(3) machinery against independent computations --- unitarity of the evolution operator, agreement with `scipy.linalg.expm`, agreement with the standard oscillation formulas, and the sign conventions of the sample Hamiltonians --- and run every example embedded in the docstrings.

**Open the notebooks.**
   ```shell
   cd /home/MyProjects/NuOscProbExact
   pip install -e ".[notebooks]"
   jupyter lab notebooks/
   ```
   Fifteen worked notebooks, numbered in reading order, covering the probabilities against baseline and against energy, matter and new physics, oscillograms, bi-probability plots, the Earth, arbitrary matter profiles, performance, the paper's own figures, the textbook approximations, mass ordering and the octant, antineutrinos, solar neutrinos, and numerical edge cases.  They carry their figures inline, so they can also just be read on GitHub.

## Performance

The probabilities are computed from a closed form, so a single one is quick — about **8 µs** for three flavors and **1 µs** for two.  Most real use, though, is a *scan*: a curve versus baseline or energy, or an oscillogram over both.  Two things make those much faster, and neither changes the answers.

### 1. Pass arrays instead of looping

Every core routine accepts a **stack** of Hamiltonians, an array of baselines, or both, and evaluates the whole thing in one call:

```python
# instead of this
prob = [oscprob3nu.probabilities_3nu(h_vacuum, l) for l in baselines]

# do this
prob = oscprob3nu.probabilities_3nu(h_vacuum, baselines)     # (N, 9)
```

The sample Hamiltonians take an array of energies too, so a scan in matter is two calls and no Python loop:

```python
h_stack = hamiltonians3nu.hamiltonian_3nu_matter(h_vacuum_energy_indep,
                                                 energies, VCC_EARTH_CRUST)
prob = oscprob3nu.probabilities_3nu(h_stack, baseline)
```

This is the single biggest win — **roughly 20–90×** — and it needs no extra dependency.  It works because the expansion's expensive part, the characteristic equation whose roots give the oscillation phases, depends on the Hamiltonian alone: a scan over baselines solves it *once* rather than once per point.

### 2. Install Numba, if the scans are large

```shell
pip install "nuoscprobexact[fast]"
```

Nothing in your code changes.  If [Numba](https://numba.pydata.org) is importable, the batched paths run as compiled machine-code loops spread over your cores instead of as a chain of NumPy array operations; if it is not, the NumPy path is used and the results are the same to round-off.

Measured on 2000-point scans, against the equivalent Python loop:

| Scan | loop | arrays | arrays + Numba |
|---|---|---|---|
| Three-flavor, vs. baseline | 38 ms | 1.8 ms (~21×) | 0.31 ms (**~120×**) |
| Three-flavor, vs. energy | 34 ms | 1.5 ms (~23×) | 0.20 ms (**~170×**) |
| Three-flavor oscillogram, 100×100 | 197 ms | 5.3 ms (~37×) | 0.85 ms (**~230×**) |
| Two-flavor, vs. baseline | 6.9 ms | 0.07 ms (~93×) | *not used — see below* |

Best of seven runs, interleaved, on one machine.  These are indicative, not precise: repeated runs vary by tens of per cent, so treat them as orders of magnitude.  [Notebook 09](notebooks/09_performance.ipynb) measures the same comparison on whatever machine runs it, which is the number to trust.

**The backend is not used where it would not help.**  For three flavors it wins at every stack size, by between two and sixteen times.  For two flavors it does not: that expansion reduces to a square root and a sine per element, which NumPy already does about as well as compiled code can, and the kernel additionally has to materialise the Hamiltonian stack.  Below fifty thousand elements the NumPy path is quicker, so it is kept; above, the kernel leads by about 1.3–1.8×.  The thresholds are measured, and the library picks whichever is faster without you doing anything.

Two costs, so the trade is visible: importing Numba takes about 140 ms against 65 ms for NumPy alone, and the first call compiles, which takes a few seconds.  The kernels are cached on disk, so later runs start in milliseconds.  This is why it is an optional extra and not a dependency.

### What you do not have to think about

* **Short stacks.** Below about ten elements the array machinery costs more than it saves, so those are evaluated one at a time automatically.
* **The scalar path.** It is deliberately left uncompiled: 8 µs is not worth a compilation pause on a first call.
* **Turning Numba off.** `fastkernels.USE_NUMBA = False` forces the NumPy path, which is how the test suite checks that the two agree.

One thing that *is* worth doing by hand: build the energy-independent part of the vacuum Hamiltonian once, outside any scan, since it does not depend on the energy.  The bundled examples all do this.

## Usage and examples

There are only two core modules: `oscprobn2nu.py` and `oscprob3nu.py`.  Each one is stand-alone (except for the dependencies described [above](#requirements)).  To use either module in your code, copy it to your project's working directory, or add their location to the paths where your environment looks for modules, *e.g.*,
```python
import sys

sys.path.append('../src')
```

In the examples below, we focus mostly on `oscprob3nu`, but what we show applies to `oscprob2nu` as well.


### Basics

Most of the time, you will be only interested in computing oscillation probabilities, not in the intermediate steps of the method.  The functions to compute and return the probabilities are `probabilities_3nu`, for the three-neutrino case, and `probabilities_2nu`, for the two-neutrino case.  Below, we show how to use them.

#### Three-neutrino oscillations

The function to compute three-neutrino probabilities is `probabilities_3nu` in the module `oscprob3nu`.  It takes as input parameters the `hamiltonian`, in the form of a 3x3 Hermitian matrix, and the baseline `L`.

This function returns the list of probabilities `Pee` (nu_e --> nu_e), `Pem` (nu_e --> nu_mu), `Pet` (nu_e --> nu_tau), `Pme` (nu_mu --> nu_e), `Pmm` (nu_mu --> nu_mu), `Pmt` (nu_mu --> nu_tau), `Pte` (nu_tau --> nu_e), `Ptm` (nu_tau --> nu_mu), and `Ptt` (nu_tau --> nu_tau).

To use it, call
```python
import oscprob3nu

hamiltonian = [[H11, H12, H13], [H21, H22, H23], [H31, H32, H33]]
Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = oscprob3nu.probabilities_3nu(hamiltonian, L)
```

#### Two-neutrino oscillations

The function to compute two-neutrino probabilities is `probabilities_2nu` in the module `oscprob2nu`.  It takes as input parameters the `hamiltonian`, in the form of a 2x2 Hermitian matrix, and the baseline `L`.

This function returns the list of probabilities `Pee` (nu_e --> nu_e), `Pem` (nu_e --> nu_mu), `Pme` (nu_mu --> nu_e), and `Pmm` (nu_mu --> nu_mu).  (These probabilities could also be `Pmm`, `Pmt`, `Pmt`, and `Ptt` instead, depending on what Hamiltonian you pass to `probabilities_2nu`.)

To use it, call
```python
import oscprob2nu

hamiltonian = [[H11, H12], [H21, H22]]
Pee, Pem, Pme, Pmm = oscprob2nu.probabilities_2nu(hamiltonian, L)
```

> **Important:** If you feed the code a non-Hermitian matrix, it will output nonsensical results

> **About the units:** The code in the modules `oscprob3nu` and `osprob2nu` does not assume units for any of the model parameters, so you need to make sure that you pass values with the correct units.  The module `globaldefs` contains conversion factors which might come in handy for this.


### Trivial example

#### Three-neutrino oscillations

As a first, trivial example, we pass an arbitrary Hamiltonian and baseline to `probabilities_3nu`:
```python
# Find this example in NuOscProbExact/examples/example_3nu_trivial.py

import oscprob3nu

hamiltonian = [
                [1.0+0.0j, 0.0+2.0j, 0.0-1.0j],
                [0.0-2.0j, 3.0+0.0j, 3.0+0.0j],
                [0.0+1.0j, 3.0-0.0j, 5.0+0.0j]
]

L = 1.0

Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = oscprob3nu.probabilities_3nu(hamiltonian, L)

print("Pee = %6.5f, Pem = %6.5f, Pet = %6.5f" % (Pee, Pem, Pet))
print("Pme = %6.5f, Pmm = %6.5f, Pmt = %6.5f" % (Pme, Pmm, Pmt))
print("Pte = %6.5f, Ptm = %6.5f, Ptt = %6.5f" % (Pte, Ptm, Ptt))
```
This returns
```shell
Pee = 0.34273, Pem = 0.41369, Pet = 0.24358
Pme = 0.41369, Pmm = 0.00485, Pmt = 0.58146
Pte = 0.24358, Ptm = 0.58146, Ptt = 0.17497
```

As expected, `Pme + Pmm + Pmt = 1`, and `Pte + Ptm + Ptt = 1`.

#### Two-neutrino oscillations

In this case, we use `probabilities_2nu`:
```python
# Find this example in NuOscProbExact/examples/example_2nu_trivial.py

import oscprob2nu

hamiltonian = [
                [1.0+0.0j, 1.0+2.0j],
                [1.0-2.0j, 3.0+0.0j]
]

L = 1.0

Pee, Pem, Pme, Pmm = oscprob2nu.probabilities_2nu(hamiltonian, L)

print("Pee = %6.5f, Pem = %6.5f" % (Pee, Pem))
print("Pme = %6.5f, Pmm = %6.5f" % (Pme, Pmm))
```
This returns
```shell
Pee = 0.66063, Pem = 0.33937
Pme = 0.33937, Pmm = 0.66063
```

As expected, `Pem == Pme` and `Pee + Pem = 1`.


### Oscillations in vacuum: fixed energy and baseline

#### Three-neutrino oscillations

Now we compute the three-neutrino oscillation probabilities in vacuum.  To do this, we can use the routine
```python
hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31)
```
that is provided in the `hamiltonians3nu` module.  It returns the 3x3 Hamiltonian for oscillations in vacuum.  The input parameters `s12`, `s23`, `s13`, `dCP`, `D21`, and `D31` are, respectively, sin(theta_12), sin(theta_23), sin(theta_13), delta_CP, Delta m_21^2, and Delta m_31^2.  For this example, we set them to their current best-fit values, which we pull from `globaldefs` (inspect that file for more information about these values).

> **Important:** The function `hamiltonian_3nu_vacuum_energy_independent` returns the Hamiltonian in vacuum **without** multiplying it by the *1/E* prefactor, where *E* is the neutrino energy.  It was done in this way so that, if we wish to compute the probabilities at different energies, we need to compute `hamiltonian_3nu_vacuum_energy_independent` only once, and then multiply it by a varying *1/E* prefactor.

```python
# Find this example in NuOscProbExact/examples/example_3nu_vacuum.py

import numpy as np

import oscprob3nu
import hamiltonians3nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

# Use the NuFit 4.0 best-fit values of the mixing parameters pulled from globaldefs
# NO means "normal ordering"; change NO to IO if you want to use inverted ordering
h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent( \
                                                                                S12_NO_BF, S23_NO_BF,
                                                                                S13_NO_BF, DCP_NO_BF,
                                                                                D21_NO_BF, D31_NO_BF)
h_vacuum = np.multiply(1./energy, h_vacuum_energy_indep)

# CONV_KM_TO_INV_EV is pulled from globaldefs; it converts km to eV^{-1}
Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = oscprob3nu.probabilities_3nu( h_vacuum,
                                                                            baseline*CONV_KM_TO_INV_EV)

print("Pee = %6.5f, Pem = %6.5f, Pet = %6.5f" % (Pee, Pem, Pet))
print("Pme = %6.5f, Pmm = %6.5f, Pmt = %6.5f" % (Pme, Pmm, Pmt))
print("Pte = %6.5f, Ptm = %6.5f, Ptt = %6.5f" % (Pte, Ptm, Ptt))
````
This returns
```shell
Pee = 0.92768, Pem = 0.01432, Pet = 0.05800
Pme = 0.04023, Pmm = 0.37887, Pmt = 0.58090
Pte = 0.03210, Ptm = 0.60680, Ptt = 0.36110
```

> **Computing anti-neutrino probabilities**: All of the examples shown in this README (and in the files inside the `examples/` directory) are for neutrinos, not anti-neutrinos.  If you wish to compute probabilities for anti-neutrinos, a simple way to do this is to pass `-dCP` instead of `dCP` to `hamiltonian_3nu_vacuum_energy_independent` (or to `hamiltonian_2nu_vacuum_energy_independent`).

> **About `globaldefs`**: This module contains physical constants and unit-conversion constants that are used in the examples and that you can use in your code.

Sometimes, you might be interested also in returning the coefficients `h1`, ..., `h8` of the expansion of the Hamiltonian in terms of Gell-Mann matrices (Table II in the paper), the coefficients `u0`, ..., `u8` of the SU(3) expansion of the associated time-evolution operator (Eqs. (13) and (14) in the paper), or the time-evolution operator `evol_operator` itself, as a 3x3 matrix (Eq. (15) in the paper).  See the paper [arXiv:1904.12391](http://arxiv.org/abs/1904.12391) for details on these quantities.

The module `oscprob3nu` has functions to do this:
```python
# Find this example in NuOscProbExact/examples/example_3nu_vacuum_coeffs.py

import numpy as np

import oscprob3nu
import hamiltonians3nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent( \
                                                                                S12_NO_BF, S23_NO_BF,
                                                                                S13_NO_BF, DCP_NO_BF,
                                                                                D21_NO_BF, D31_NO_BF)
h_vacuum = np.multiply(1./energy, h_vacuum_energy_indep)

h1, h2, h3, h4, h5, h6, h7, h8 = oscprob3nu.hamiltonian_3nu_coefficients(h_vacuum)
print('h1: {:.4e}'.format(h1))
print('h2: {:.4e}'.format(h2))
print('h3: {:.4e}'.format(h3))
print('h4: {:.4e}'.format(h4))
print('h5: {:.4e}'.format(h5))
print('h6: {:.4e}'.format(h6))
print('h7: {:.4e}'.format(h7))
print('h8: {:.4e}'.format(h8))
print()

u0, u1, u2, u3, u4, u5, u6, u7, u8 = oscprob3nu.evolution_operator_3nu_u_coefficients(  \
                                                                            h_vacuum,
                                                                            baseline*CONV_KM_TO_INV_EV)
print('u0: {:.4f}'.format(u0))
print('u1: {:.4f}'.format(u1))
print('u2: {:.4f}'.format(u2))
print('u3: {:.4f}'.format(u3))
print('u4: {:.4f}'.format(u4))
print('u5: {:.4f}'.format(u5))
print('u6: {:.4f}'.format(u6))
print('u7: {:.4f}'.format(u7))
print('u8: {:.4f}'.format(u8))
print()

evol_operator = oscprob3nu.evolution_operator_3nu(h_vacuum, baseline*CONV_KM_TO_INV_EV)
print('U3 = ')
with np.printoptions(precision=3, suppress=True):
    print(np.array(evol_operator))
```
This returns
```shell
h1: -1.0187e-13
h2: -8.4997e-14
h3: -3.4583e-13
h4: -1.0848e-13
h5: -7.2033e-14
h6: 5.9597e-13
h7: 1.5392e-15
h8: -8.2865e-14

u0: -0.3794+0.5072j
u1: 0.0318+0.1167j
u2: -0.0257+0.1095j
u3: -0.1270+0.4507j
u4: -0.1066+0.1569j
u5: -0.0217+0.0928j
u6: 0.1383-0.7580j
u7: 0.0093-0.0040j
u8: -0.0323+0.1084j

[[-0.893+0.362j -0.142+0.141j -0.179-0.014j]
 [-0.091-0.078j  0.009+0.615j  0.767+0.134j]
 [-0.135-0.199j  0.749+0.142j -0.254+0.545j]]
```



#### Two-neutrino oscillations

To compute the two-neutrino oscillation probabilities in vacuum, we can use the routine
```python
hamiltonian_2nu_vacuum_energy_independent(sth, Dm2)
```
that is provided in the `hamiltonians2nu` module.  The input parameters `sth`, and `Dm2` are, respectively, sin(theta), and Delta m^2.  For this example, we set them to current best-fit values for atmospheric neutrinos.

```python
# Find this example in NuOscProbExact/examples/example_2nu_vacuum.py

import numpy as np

import oscprob2nu
import hamiltonians2nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_energy_indep = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(S23_NO_BF, D31_NO_BF)
h_vacuum = np.multiply(1./energy, h_vacuum_energy_indep)

Pee, Pem, Pme, Pmm = oscprob2nu.probabilities_2nu(h_vacuum, baseline*CONV_KM_TO_INV_EV)

print("Pee = %6.5f, Pem = %6.5f" % (Pee, Pem))
print("Pme = %6.5f, Pmm = %6.5f" % (Pme, Pmm))
````
This returns
```shell
Pee = 0.29595, Pem = 0.70405
Pme = 0.70405, Pmm = 0.29595
```

Like in the three-neutrino case, we can also return the coefficients `h1`, `h2`, `h3` of the expansion of the Hamiltonian in terms of Pauli matrices (Table I in the paper), or the time-evolution operator `evol_operator` itself, as a 2x2 matrix (Eq. (5) in the paper).
```python
# Find this example in NuOscProbExact/examples/example_2nu_vacuum_coeffs.py

import numpy as np

import oscprob2nu
import hamiltonians2nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_energy_indep = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(S23_NO_BF, D31_NO_BF)
h_vacuum = np.multiply(1./energy, h_vacuum_energy_indep)

h1, h2, h3 = oscprob2nu.hamiltonian_2nu_coefficients(h_vacuum)
print('h1: {:.4e}'.format(h1))
print('h2: {:.4e}'.format(h2))
print('h3: {:.4e}'.format(h3))
print()

evol_operator = oscprob2nu.evolution_operator_2nu(h_vacuum, baseline*CONV_KM_TO_INV_EV)
print('U2 = ')
with np.printoptions(precision=3, suppress=True):
    print(np.array(evol_operator))
```
This returns
```shell
h1: 6.2270e-13
h2: -0.0000e+00
h3: 1.0352e-13

U2 =
[[-0.526+0.139j  0.   +0.839j]
 [ 0.   +0.839j -0.526-0.139j]]
```


### Three-neutrino oscillations in vacuum: fixed energy, varying baseline

Now we fix the energy at, say, 10 MeV, and vary the baseline between 1 and 500 km.  We use a fine grid in `L` so that the oscillations are clearly rendered.

```python
import numpy as np

import oscprob3nu
import hamiltonians3nu
from globaldefs import *

energy = 1.e7     # Neutrino energy [eV]

# Baselines, L
log10_l_min = 0.0  # log10 [km]
log10_l_max = 3.0  # log10 [km]
log10_l_npts = 1000
log10_l_val = np.linspace(log10_l_min, log10_l_max, log10_l_npts)  # [km]
l_val = [10.**x for x in log10_l_val]


h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent( \
                                                                                S12_NO_BF, S23_NO_BF,
                                                                                S13_NO_BF, DCP_NO_BF,
                                                                                D21_NO_BF, D31_NO_BF)
h_vacuum = np.multiply(1./energy, h_vacuum_energy_indep)

# Each element of prob: [Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt]
prob = [oscprob3nu.probabilities_3nu(h_vacuum, CONV_KM_TO_INV_EV*l) for l in l_val]
prob_ee = [x[0] for x in prob]  # Pee
prob_em = [x[1] for x in prob]  # Pem
prob_et = [x[2] for x in prob]  # Pet
```

To plot the data:
```python
from pylab import *
from matplotlib import *
import matplotlib as mpl

fig = plt.figure(figsize=[9,9])
ax = fig.add_subplot(1,1,1)

ax.plot(l_val, prob_ee, label=r'$P_{\nu_e \to \nu_e}$', color='C0', zorder=1)
ax.plot(l_val, prob_em, label=r'$P_{\nu_e \to \nu_\mu}$', color='C1', zorder=1)
ax.plot(l_val, prob_et, label=r'$P_{\nu_e \to \nu_\tau}$', color='C2', zorder=1)

ax.set_xlabel(r'Baseline $L$ [km]', fontsize=25)
ax.set_ylabel(r'Three-neutrino probability', fontsize=25)
ax.legend(loc='center left', frameon=False)
ax.set_xlim([10.**log10_l_min, 10.**log10_l_max])
ax.set_xscale('log')
ax.set_ylim([0.0, 1.0])

plt.show()
````

The same curve, and the equivalents for oscillations in matter, with non-standard interactions, and in a Lorentz invariance-violating background, are drawn in [notebook 02](notebooks/02_vacuum_oscillations.ipynb) and [notebook 03](notebooks/03_matter_nsi_liv.ipynb), which store their figures inline:

<img align="middle" class="center" src="https://github.com/mbustama/NuOscProbExact/blob/main/img/prob_3nu_vacuum_vs_baseline_ee_em_et.png" width="400"/>

The default mixing parameters used throughout come from the `globaldefs` module, which carries the NuFit best fit for both mass orderings.  For more on the scenarios themselves, see the paper [arXiv:1904.12391](http://arxiv.org/abs/1904.12391).


### Three-neutrino oscillations in vacuum: fixed baseline, varying energy

Now we fix the baseline at, say, 1300 km, and vary the energy between 100 MeV and 10 GeV.

```python
import numpy as np

import oscprob3nu
import hamiltonians3nu
from globaldefs import *

baseline = 1.3e3                       # Baseline [km]
baseline = baseline*CONV_KM_TO_INV_EV  # [eV^{-1}]

# Neutrino energies
log10_energy_min = -1.0 # [GeV]
log10_energy_max = 1.0  # [GeV]
log10_energy_npts = 200
log10_energy = np.linspace( log10_energy_min,
                            log10_energy_max,
                            log10_energy_npts)
energy = [10.**x for x in log10_energy] # [GeV]

h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent( \
                                                                                S12_NO_BF, S23_NO_BF,
                                                                                S13_NO_BF, DCP_NO_BF,
                                                                                D21_NO_BF, D31_NO_BF)

# Each element of prob: [Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt]
prob = [oscprob3nu.probabilities_3nu(np.multiply(1./x/1.e9, h_vacuum_energy_indep), baseline) \
        for x in energy]
prob_ee = [x[0] for x in prob]  # Pee
prob_em = [x[1] for x in prob]  # Pem
prob_et = [x[2] for x in prob]  # Pet
```

To plot the data:
```python
from pylab import *
from matplotlib import *
import matplotlib as mpl

fig = plt.figure(figsize=[9,9])
ax = fig.add_subplot(1,1,1)

ax.plot(energy, prob_ee, label=r'$P_{\nu_e \to \nu_e}$', color='C0', zorder=1)
ax.plot(energy, prob_em, label=r'$P_{\nu_e \to \nu_\mu}$', color='C1', zorder=1)
ax.plot(energy, prob_et, label=r'$P_{\nu_e \to \nu_\tau}$', color='C2', zorder=1)

ax.set_xlabel(r'Neutrino energy $E$ [GeV]', fontsize=25)
ax.set_ylabel(r'Three-neutrino probability', fontsize=25)
ax.legend(loc='center right', frameon=False)
ax.set_xlim([10.**log10_energy_min, 10.**log10_energy_max])
ax.set_xscale('log')
ax.set_ylim([0.0, 1.0])

plt.show()
````

The same curve is drawn in [notebook 02](notebooks/02_vacuum_oscillations.ipynb), which stores its figures inline:

<img align="middle" class="center" src="https://github.com/mbustama/NuOscProbExact/blob/main/img/prob_3nu_vacuum_vs_energy_ee_em_et.png" width="400"/>

For more on the oscillation scenarios themselves, see the paper [arXiv:1904.12391](http://arxiv.org/abs/1904.12391).


### Three-neutrino oscillations in vacuum: fixed baseline, varying energy

Now we fix the baseline at, say, 1300 km, and vary the energy between 100 MeV and 10 GeV.

```python
import numpy as np

import oscprob3nu
import hamiltonians3nu
from globaldefs import *

baseline = 1.3e3                       # Baseline [km]
baseline = baseline*CONV_KM_TO_INV_EV  # [eV^{-1}]

# Neutrino energies
log10_energy_min = -1.0 # [GeV]
log10_energy_max = 1.0  # [GeV]
log10_energy_npts = 200
log10_energy = np.linspace( log10_energy_min,
                            log10_energy_max,
                            log10_energy_npts)
energy = [10.**x for x in log10_energy] # [GeV]

h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent( \
                                                                                S12_NO_BF, S23_NO_BF,
                                                                                S13_NO_BF, DCP_NO_BF,
                                                                                D21_NO_BF, D31_NO_BF)

# Each element of prob: [Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt]
prob = [oscprob3nu.probabilities_3nu(np.multiply(1./x/1.e9, h_vacuum_energy_indep), baseline) \
        for x in energy]
prob_ee = [x[0] for x in prob]  # Pee
prob_em = [x[1] for x in prob]  # Pem
prob_et = [x[2] for x in prob]  # Pet
```

To plot the data:
```python
from pylab import *
from matplotlib import *
import matplotlib as mpl

fig = plt.figure(figsize=[9,9])
ax = fig.add_subplot(1,1,1)

ax.plot(energy, prob_ee, label=r'$P_{\nu_e \to \nu_e}$', color='C0', zorder=1)
ax.plot(energy, prob_em, label=r'$P_{\nu_e \to \nu_\mu}$', color='C1', zorder=1)
ax.plot(energy, prob_et, label=r'$P_{\nu_e \to \nu_\tau}$', color='C2', zorder=1)

ax.set_xlabel(r'Neutrino energy $E$ [GeV]', fontsize=25)
ax.set_ylabel(r'Three-neutrino probability', fontsize=25)
ax.legend(loc='center right', frameon=False)
ax.set_xlim([10.**log10_energy_min, 10.**log10_energy_max])
ax.set_xscale('log')
ax.set_ylim([0.0, 1.0])

plt.show()
````

The same curve is drawn in [notebook 02](notebooks/02_vacuum_oscillations.ipynb):


<img align="middle" class="center" src="https://github.com/mbustama/NuOscProbExact/blob/main/img/prob_3nu_vacuum_vs_energy_ee_em_et.png" width="400"/>



### Three-neutrino oscillations in matter

For oscillation in matter, we proceed in an analogous way as for oscillations in vacuum.  To compute the Hamiltonian in matter, we can use the routine `hamiltonian_matter` in the module `hamiltonians3nu`.  First, we need to compute the energy-independent `h_vacuum_energy_independent`, and then pass it to `hamiltonian_3nu_matter`, together with the `energy` and the neutrino-electron charged-current potential `VCC`, with V_CC = sqrt(2.0) * G_F * n_e.  This routine is called as:
```python
hamiltonian_3nu_matter(h_vacuum_energy_independent, energy, VCC)
```

In the example below, we set the matter potential to `VCC_EARTH_CRUST`, which is computed using the averaage electron density of the crust of the Earth (3 g cm^{-3}), and is read from `globaldefs`.
```python
# Find this example in NuOscProbExact/examples/example_3nu_matter.py

import oscprob3nu
import hamiltonians3nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent( \
                                                                                S12_NO_BF, S23_NO_BF,
                                                                                S13_NO_BF, DCP_NO_BF,
                                                                                D21_NO_BF, D31_NO_BF)

# Units of VCC_EARTH_CRUST: [eV]
h_matter = hamiltonians3nu.hamiltonian_3nu_matter(h_vacuum_energy_indep, energy, VCC_EARTH_CRUST)

Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = oscprob3nu.probabilities_3nu( h_matter,
                                                                            baseline*CONV_KM_TO_INV_EV)

print("Pee = %6.5f, Pem = %6.5f, Pet = %6.5f" % (Pee, Pem, Pet))
print("Pme = %6.5f, Pmm = %6.5f, Pmt = %6.5f" % (Pme, Pmm, Pmt))
print("Pte = %6.5f, Ptm = %6.5f, Ptt = %6.5f" % (Pte, Ptm, Ptt))
````
This returns
```shell
Pee = 0.95262, Pem = 0.00623, Pet = 0.04115
Pme = 0.02590, Pmm = 0.37644, Pmt = 0.59766
Pte = 0.02148, Ptm = 0.61733, Ptt = 0.36119
```

### Three-neutrino oscillations in matter with non-standard interactions (NSI)

For oscillation in matter with NSI, we can use the routine `hamiltonian_3nu_nsi` in the module `hamiltonians3nu`.  First, we need to compute `h_vacuum_energy_independent`, and then pass it to `hamiltonian_nsi`, together with `energy`, `VCC`, and a vector `eps` containing the NSI strength parameters, *i.e.*,
```python
eps = [eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt]
```
This routine is called as
```python
hamiltonian_3nu_nsi(h_vacuum_energy_independent, energy, VCC, eps)
```

In the example below, we set `eps` to its default value `EPS_3` pulled from `globaldefs`:
```python
# Find this example in NuOscProbExact/examples/example_3nu_nsi.py

import oscprob3nu
import hamiltonians3nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent( \
                                                                                S12_NO_BF, S23_NO_BF,
                                                                                S13_NO_BF, DCP_NO_BF,
                                                                                D21_NO_BF, D31_NO_BF)

# EPS_3 is the 3x3 matrix of NSI strength parameters, read from globaldefs; see that file for the values
h_nsi = hamiltonians3nu.hamiltonian_3nu_nsi(h_vacuum_energy_indep, energy, VCC_EARTH_CRUST, EPS_3)

Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = oscprob3nu.probabilities_3nu( h_nsi,
                                                                            baseline*CONV_KM_TO_INV_EV)

print("Pee = %6.5f, Pem = %6.5f, Pet = %6.5f" % (Pee, Pem, Pet))
print("Pme = %6.5f, Pmm = %6.5f, Pmt = %6.5f" % (Pme, Pmm, Pmt))
print("Pte = %6.5f, Ptm = %6.5f, Ptt = %6.5f" % (Pte, Ptm, Ptt))
````
This returns
```shell
Pee = 0.92494, Pem = 0.01758, Pet = 0.05749
Pme = 0.03652, Pmm = 0.32524, Pmt = 0.63824
Pte = 0.03855, Ptm = 0.65718, Ptt = 0.30427
```


### Three-neutrino oscillations in a Lorentz invariance-violating (LIV) background

For oscillation LIV, we can use the routine `hamiltonian_3nu_liv` in the module `hamiltonians3nu`.  As before, first, we need to compute `h_vacuum_energy_independent`, and then pass it to `hamiltonian_3nu_liv`, together with `energy` and the following LIV parameters: `sxi12` (sin(xi_12)), `sxi23` (sin(xi_23)), `sxi13` (sin(xi_13)), `dxiCP` (new CP-violation phase), `b1` (first eigenvalue of the LIV operator), `b2` (second eigenvalue), `b3` (third eigenvalue), and `Lambda` (energy scale of LIV).

This routine is called as
```python
hamiltonian_3nu_liv(h_vacuum_energy_independent, energy, sxi12, sxi23, sxi13, dxiCP, b1, b2, b3, Lambda)
```

In the example below, we set the LIV parameters to their default values pulled from `globaldefs`:
```python
# Find this example in NuOscProbExact/examples/example_3nu_liv.py

import oscprob3nu
import hamiltonians3nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(  S12_BF, S23_BF,
                                                                                    S13_BF, DCP_BF,
                                                                                    D21_BF, D31_BF)

# The values of the LIV parameters (SXI12, SXI23, SXI13, DXICP, B1, B2, B3, LAMBDA) are read
# from globaldefs
h_liv = hamiltonians3nu.hamiltonian_3nu_liv(h_vacuum_energy_indep, energy,
                                            SXI12, SXI23, SXI13, DXICP,
                                            B1, B2, B3, LAMBDA)

Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = oscprob3nu.probabilities_3nu( h_liv,
                                                                            baseline*CONV_KM_TO_INV_EV)

print("Pee = %6.5f, Pem = %6.5f, Pet = %6.5f" % (Pee, Pem, Pet))
print("Pme = %6.5f, Pmm = %6.5f, Pmt = %6.5f" % (Pme, Pmm, Pmt))
print("Pte = %6.5f, Ptm = %6.5f, Ptt = %6.5f" % (Pte, Ptm, Ptt))
````
This returns
```shell
Pee = 0.92721, Pem = 0.05299, Pet = 0.01980
Pme = 0.05609, Pmm = 0.25288, Pmt = 0.69103
Pte = 0.01670, Ptm = 0.69412, Ptt = 0.28917
```

### Arbitrary Hamiltonians

Of course, you can supply your custom Hamiltonian and compute the associated oscillation probabilities; see [Trivial example](#trivial-example) above.  Usually, you will want to add an extra term from your preferred model to the vacuum Hamiltonian.  To do that, take a cue from the examples above.

In the following example, the function `hamiltonian_mymodel`, supplied by you, should return a 3x3 matrix:
```python
import oscprob3nu
import hamiltonians3nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(  S12_BF, S23_BF,
                                                                                    S13_BF, DCP_BF,
                                                                                    D21_BF, D31_BF)
h_vacuum = np.multiply(1./energy, h_vacuum_energy_indep)
h_mymodel = h_vacuum + hamiltonian_mymodel(mymodel_parameters)

Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = oscprob3nu.probabilities_3nu( h_mymodel,
                                                                            baseline*CONV_KM_TO_INV_EV)

```
Though we do not show it here, `hamiltonian_mymodel` could also depend on `energy`.  The code for two-neutrino oscillations is analogous, but `hamiltonian_mymodel` should return a 2x2 matrix instead.


## Notebooks

Fifteen worked notebooks live in [`notebooks/`](notebooks/), numbered in reading order.  They carry their figures inline, so they render on GitHub without being run:

| Notebook | What it covers |
|---|---|
| [01 Basics](notebooks/01_basics.ipynb) | Units, one probability, and why to pass arrays rather than loop |
| [02 Vacuum oscillations](notebooks/02_vacuum_oscillations.ipynb) | Against baseline and against energy, checked against the textbook formula |
| [03 Matter, NSI, LIV](notebooks/03_matter_nsi_liv.ipynb) | Constant-density matter and two kinds of new physics |
| [04 Oscillograms](notebooks/04_oscillogram.ipynb) | A 240x240 energy-baseline map in a single call |
| [05 Bi-probability](notebooks/05_biprobability.ipynb) | CP ellipses, in vacuum and in matter |
| [06 The Earth and PREM](notebooks/06_earth_and_prem.ipynb) | The density profile, chord geometry, slabs, and their convergence |
| [07 Through the Earth](notebooks/07_earth_probabilities.ipynb) | Zenith-angle scans, an Earth oscillogram, and real baselines |
| [08 Unusual density profiles](notebooks/08_unusual_density_profiles.ipynb) | Castle-wall and serrated profiles, and parametric enhancement |
| [09 Performance](notebooks/09_performance.ipynb) | Looping versus broadcasting, and the compiled backend, measured live |
| [10 The paper's figures](notebooks/10_paper_figures.ipynb) | The two figures from [arXiv:1904.12391](https://arxiv.org/abs/1904.12391) |
| [11 Exact vs approximations](notebooks/11_exact_vs_approximations.ipynb) | Where the familiar formulas agree, and where they do not |
| [12 Ordering and octant](notebooks/12_ordering_and_octant.ipynb) | Normal against inverted, and the θ₂₃ octant degeneracy |
| [13 Antineutrinos](notebooks/13_antineutrinos.ipynb) | Conjugate *and* flip the potential — and two ways to get it wrong |
| [14 Solar and the MSW resonance](notebooks/14_solar_and_adiabatic_msw.ipynb) | The adiabatic resonance, validated — and why slabs are the wrong tool for it |
| [15 Numerical edge cases](notebooks/15_numerical_edge_cases.ipynb) | Degenerate spectra, and what returns a number instead of NaN |

Run them with `pip install -e ".[notebooks]"` and `jupyter lab notebooks/`.  Every one of them is executed by CI, so an example that stops working fails the build.


## Documentation and help

All of the modules provided in **NuOscProbExact** have been documented using Python docstrings, written in [numpydoc](https://numpydoc.readthedocs.io/) format so that they can be rendered directly by [Sphinx](https://www.sphinx-doc.org/) with the `sphinx.ext.napoleon` and `numpydoc` extensions.  They are human-readable by opening the source `.py` files.  Alternatively, they can be printed from within an interactive Python session.

Every `Examples` block in the docstrings is executed when the documentation is built, so the results shown on the [API page](https://mbustama.github.io/NuOscProbExact/functions.html) are produced by the code rather than pasted beside it, and cannot drift.  The regression suite runs the same blocks on every supported Python (`tests/test_docstrings.py`), which the documentation build --- one job, one interpreter --- would not catch.

A full Sphinx project lives in `docs/`.  Build it with
```shell
pip install -r docs/requirements.txt
cd docs && make html
```
and open `docs/build/html/index.html`.  It contains an installation guide, a quickstart, a description of the [method](docs/source/methodology.rst) and its sign conventions, the API reference generated from the docstrings, a bibliography, and the changelog.

Notable changes between versions are recorded in [CHANGELOG.md](CHANGELOG.md), which the documentation renders as its own page, so there is a single source of truth.

To view the documentation of a module from within an interactive Python session, run, *e.g.*,
```python
import oscprob3nu

print(oscprob3nu.__doc__)
```
This will print to screen a description of what the module does (in this example, `oscprob3nu`) and a list of the functions that it contains, including a description of each.

To view the documentation of a particular function from within an interactive Python session, run, *e.g.*,
```python
import oscprob3nu

help(oscprob3nu.hamiltonian_3nu_coefficients)
```
This will print to screen a description of what the function does (in the example above, `oscprob3nu.hamiltonian_3nu_coefficients`), a list and description of its input parameters, and a description of the values that it returns.


## Citing

If you use **NuOscProbExact** in your work, we ask you that you please cite the following paper: Mauricio Bustamante, *NuOscProbExact: a general-purpose code to compute exact two-flavor and three-flavor neutrino oscillation probabilities* ([arXiv:1904.12391](http://arxiv.org/abs/1904.12391)).

If you are citing **NuOscProbExact** in a document that will be uploaded to the arXiv, please consider using the LaTeX or BibTeX entries provided by INSPIRE ([link here](http://inspirehep.net/record/1731803/export/hx)):
```
@article{Bustamante:2019ggq,
      author         = "Bustamante, Mauricio",
      title          = "{NuOscProbExact: a general-purpose code to compute
                        exact two-flavor and three-flavor neutrino
                        oscillation probabilities}",
      year           = "2019",
      eprint         = "1904.12391",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      SLACcitation   = "%%CITATION = ARXIV:1904.12391;%%"
}
```





