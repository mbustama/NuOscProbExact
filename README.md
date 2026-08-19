[![tests](https://github.com/mbustama/NuOscProbExact/actions/workflows/tests.yml/badge.svg)](https://github.com/mbustama/NuOscProbExact/actions/workflows/tests.yml)
[![Code Quality](https://github.com/mbustama/NuOscProbExact/actions/workflows/lint.yml/badge.svg)](https://github.com/mbustama/NuOscProbExact/actions/workflows/lint.yml)
[![codecov](https://codecov.io/gh/mbustama/NuOscProbExact/branch/main/graph/badge.svg)](https://codecov.io/gh/mbustama/NuOscProbExact)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue.svg)](https://mbustama.github.io/NuOscProbExact/)
[![PyPI](https://img.shields.io/pypi/v/nuoscprobexact.svg)](https://pypi.org/project/nuoscprobexact/)
[![Downloads](https://pepy.tech/badge/nuoscprobexact)](https://pepy.tech/project/nuoscprobexact)
[![arXiv](https://img.shields.io/badge/arXiv-1904.12391-orange.svg)](https://arxiv.org/abs/1904.12391)
[![DOI](https://zenodo.org/badge/182178323.svg)](https://zenodo.org/badge/latestdoi/182178323)
[![INSPIRE](https://img.shields.io/badge/INSPIRE-cited%20by-003a6c.svg)](https://inspirehep.net/literature?sort=mostrecent&size=250&page=1&q=refersto%3Arecid%3A1731803&ui-citation-summary=true)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Code style: ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

# NuOscProbExact
Code to compute exact two-, three- and four-neutrino oscillation probabilities using SU(2), SU(3) and SU(4) expansions

> **Note:** The oscillation probabilities are computed exactly, with no approximation beyond floating-point round-off.  A regression test suite lives in `tests/` and can be run with `pytest`, and the results are cross-checked against [nuSQuIDS](https://github.com/arguelles/nuSQuIDS), an independent external code — see [notebook 17](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/17_cross_checks.ipynb).

## What you can compute

Every figure below is produced by a notebook in [`notebooks/`](https://github.com/mbustama/NuOscProbExact/tree/main/notebooks/), and the link under each one goes to the code that drew it.  The documentation collects the same material, with runnable snippets, on its [numerical recipes](https://mbustama.github.io/NuOscProbExact/recipes.html) page.

| | |
|:--:|:--:|
| <img src="https://raw.githubusercontent.com/mbustama/NuOscProbExact/main/img/gallery/gallery_vacuum.png" width="380"/><br/>**Oscillation probabilities** against baseline or energy, for two, three or four flavors.<br/>[notebook 02](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/02_vacuum_oscillations.ipynb) | <img src="https://raw.githubusercontent.com/mbustama/NuOscProbExact/main/img/gallery/gallery_matter.png" width="380"/><br/>**Matter, NSI and Lorentz-invariance violation** — each just a different Hermitian matrix.<br/>[notebook 03](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/03_matter_nsi_liv.ipynb) |
| <img src="https://raw.githubusercontent.com/mbustama/NuOscProbExact/main/img/gallery/gallery_oscillogram.png" width="380"/><br/>**Oscillograms** over energy and baseline: 57 600 probabilities in a single call.<br/>[notebook 04](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/04_oscillogram.ipynb) | <img src="https://raw.githubusercontent.com/mbustama/NuOscProbExact/main/img/gallery/gallery_biprobability.png" width="380"/><br/>**CP violation**, as bi-probability ellipses in vacuum and in matter.<br/>[notebook 05](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/05_biprobability.ipynb) |
| <img src="https://raw.githubusercontent.com/mbustama/NuOscProbExact/main/img/gallery/gallery_prem.png" width="380"/><br/>**The Earth's density**, from the Preliminary Reference Earth Model.<br/>[notebook 06](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/06_earth_and_prem.ipynb) | <img src="https://raw.githubusercontent.com/mbustama/NuOscProbExact/main/img/gallery/gallery_earth.png" width="380"/><br/>**Neutrinos through the Earth**, in energy and zenith angle, or between two named sites.<br/>[notebook 07](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/07_earth_probabilities.ipynb) |
| <img src="https://raw.githubusercontent.com/mbustama/NuOscProbExact/main/img/gallery/gallery_profiles.png" width="380"/><br/>**Arbitrary matter profiles** — castle walls and worse, exactly.<br/>[notebook 08](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/08_unusual_density_profiles.ipynb) | <img src="https://raw.githubusercontent.com/mbustama/NuOscProbExact/main/img/gallery/gallery_ordering.png" width="380"/><br/>**Mass ordering and the θ₂₃ octant**, separated by matter through the Earth.<br/>[notebook 12](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/12_ordering_and_octant.ipynb) |
| <img src="https://raw.githubusercontent.com/mbustama/NuOscProbExact/main/img/gallery/gallery_sterile.png" width="380"/><br/>**Four flavors: a 3+1 sterile state**, resolved at a short baseline.<br/>[notebook 16](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/16_four_neutrinos.ipynb) | <img src="https://raw.githubusercontent.com/mbustama/NuOscProbExact/main/img/gallery/gallery_sterile_earth.png" width="380"/><br/>**The sterile matter resonance** through the Earth, in energy and zenith angle.<br/>[notebook 16](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/16_four_neutrinos.ipynb) |
| <img src="https://raw.githubusercontent.com/mbustama/NuOscProbExact/main/img/gallery/gallery_long_range.png" width="380"/><br/>**A Hamiltonian of your own**, carried through a varying profile — here a long-range force through the Earth.<br/>[notebook 20](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/20_arbitrary_hamiltonian.ipynb) | |

## Contents

1. [What you can compute](#what-you-can-compute)

2. [What is NuOscProbExact?](#what-is-nuoscprobexact)
   1. [What it does](#what-it-does)
   2. [What it does not do](#what-it-does-not-do)
   3. [When to use Magnus instead](#when-to-use-magnus-instead)

3. [Requirements](#requirements)

4. [Installation](#installation)

5. [Performance](#performance)

6. [Usage and examples](#usage-and-examples)
   1. [Basics](#basics)
   2. [A first probability](#a-first-probability)
   3. [Whole scans in one call](#whole-scans-in-one-call)
   4. [Four flavors: a 3+1 sterile state](#four-flavors-a-31-sterile-state)
   5. [Arbitrary Hamiltonians](#arbitrary-hamiltonians)
   6. [Where the rest is](#where-the-rest-is)

7. [Notebooks](#notebooks)

8. [Documentation and help](#documentation-and-help)

9. [Citing](#citing)

10. [License](#license)


## What is NuOscProbExact?

**NuOscProbExact** is a Python implementation of the method developed by [Ohlsson & Snellman](https://arxiv.org/abs/hep-ph/9910546) to compute exact neutrino oscillation probabilities for arbitrary time-independent Hamiltonians.  The method was revisited and the code presented in the paper *NuOscProbExact: a general-purpose code to compute exact two-flavor and three-flavor neutrino oscillation probabilities* ([arXiv:1904.12391](https://arxiv.org/abs/1904.12391)), by Mauricio Bustamante.

The paper covers two and three flavors; the code has since been extended to **four**, through the SU(4) algebra, which brings 3+1 sterile scenarios into scope.  Four is where the closed form ends — see [why](https://mbustama.github.io/NuOscProbExact/methodology.html#why-the-method-stops-at-four).

The method relies on expansions of the Hamiltonian and time-evolution operators in terms of SU(2), SU(3) and SU(4) matrices in order to obtain concise, analytical, and exact expressions for the probabilities, that are also easy to implement and evaluate.  For details of the method, see the paper above; the four-flavor extension is documented in [the methodology page](https://mbustama.github.io/NuOscProbExact/methodology.html).

### What it does

* **Exact probabilities for any Hermitian 2×2, 3×3 or 4×4 Hamiltonian.**  There is no approximation beyond floating-point round-off, and the probabilities agree with the independently written [nuSQuIDS](https://github.com/arguelles/nuSQuIDS) to round-off once conventions are matched ([notebook 17](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/17_cross_checks.ipynb)).  Oscillations in vacuum, in matter, with non-standard interactions, in a Lorentz invariance-violating background and with sterile states are not special cases in the code — each is a different matrix handed to the same routine.
* **Four flavors, for 3+1 sterile scenarios.**  `oscprob4nu` carries the same closed-form treatment to SU(4), which is the last place it reaches: at five flavors the eigenvalues stop being expressible in radicals, and that is a theorem rather than a missing feature.  A 3+1 system is closed and unitary over all four states, so it sits squarely inside the method's assumptions rather than "leaking" out of a three-flavor block.
* **The evolution operator itself**, not only the probabilities, so it can be composed across segments or used to propagate a density matrix.
* **Whole scans in one call.**  Every core routine accepts a stack of Hamiltonians, an array of baselines, or both broadcast against each other, which is tens of times faster than the equivalent Python loop and gives identical results.
* **Piecewise-constant matter.**  `slabs` propagates across a sequence of adjacent slabs of arbitrary width and density, solving each exactly and multiplying the operators, at two, three or four flavors.
* **The Earth.**  `earth` builds those slabs from the Preliminary Reference Earth Model, and computes probabilities along a given zenith angle or between two of fifteen predefined locations.  A 3+1 crossing is included: the sterile state does not feel the neutral-current potential, so that potential stops cancelling and `earth` builds it per slab.  PREM fixes the density but not the composition, so the electron fraction is the caller's: it defaults to one half throughout, and `earth.electron_fraction_prem` supplies the value of each layer instead — core, mantle, crust and ocean.  That default will change to the layered values in a future release.
* **A compiled backend, installed by default.**  Large batched calls run on compiled kernels, and so does the whole slab composition behind `slabs` and `earth` — an Earth crossing is 7–14× quicker depending on the flavor count.  `numba` comes with the package as of 1.13.0; where it is unavailable or switched off, the NumPy path is used and the answers are the same to round-off.

### What it does not do

* **Hamiltonians that vary continuously along the trajectory.**  See [When to use Magnus instead](#when-to-use-magnus-instead) below — this is the one case where a different tool is the right answer, and it is worth knowing before you start.
* **More than four flavors.**  The expansions run to SU(4) and stop, because the closed form does: solving for the eigenvalues means solving the characteristic polynomial in radicals, and at degree five Abel–Ruffini says that cannot be done.  Four flavors covers 3+1, which is the case people actually ask for.
* **Neutrino production, cross sections, fluxes or detector response.**  This computes oscillation probabilities and nothing downstream of them.
* **Fitting or statistics.**  There is no likelihood machinery here; the probabilities are meant to be handed to whatever does that.

### When to use Magnus instead

**NuOscProbExact** assumes the Hamiltonian is constant, or piecewise constant.  Everything it is good at follows from that — and so does the one case where it is the wrong tool.

Use **Magnus** instead when **the Hamiltonian varies continuously and appreciably over an oscillation length**.  A smoothly varying profile can always be approximated by slabs, but then the step size is set by the oscillation rather than by the density, and the slab count grows until the calculation is neither exact nor quick.

| Situation | Use | Because |
|---|---|---|
| Constant density | **NuOscProbExact** | One closed form, no integration |
| Piecewise constant, tens of layers — the Earth through PREM | **NuOscProbExact** (`slabs`, `earth`) | Each layer solved exactly, operators multiplied |
| Smoothly varying, slow against the oscillation | either | Slabbing converges quickly |
| Smoothly varying, fast against the oscillation — the Sun, adiabatic MSW | **Magnus** | Slabbing needs ~10<sup>4</sup> steps per resonance crossing |
| Open systems: decay, decoherence | neither | Needs a Lindblad solver, not a unitary one |

[Notebook 14](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/14_solar_and_adiabatic_msw.ipynb) works the solar case through and shows exactly where the wall is, rather than asserting it.

**NuOscProbExact** was developed by Mauricio Bustamante.  If you use it in your work, please follow the directions on [Citing](#citing).


## Requirements

**NuOscProbExact** is fully written in Python 3.  It uses standard modules that are available, sometimes by default, as part of most Python installations, either stand-alone or via Anaconda.  `numpy` and `numba` are installed with the package, so the first four rows need nothing at all; where a row names an extra, that is the one to install.  The commands are under [Installation](#installation) below, and are not repeated here.

| To do this | You need | Extra |
|---|---|---|
| Compute probabilities (`oscprob2nu.py`, `oscprob3nu.py`, `oscprob4nu.py`) | `numpy`, `cmath` | — |
| Use the bundled sample Hamiltonians (`hamiltonians2nu.py`, `hamiltonians3nu.py`, `hamiltonians4nu.py`) | `numpy`, `cmath`, `copy` | — |
| Propagate through layered matter or the Earth (`slabs.py`, `earth.py`) | `numpy` | — |
| Go faster on large scans | `numba` | — *(installed by default)* |
| Run the notebooks (`notebooks/`) | `matplotlib`, Jupyter | `notebooks` |
| Run the regression suite (`tests/`) | `pytest`, `scipy`, `coverage`, `pytest-cov` | `test` |
| Build the documentation | Sphinx and friends | `docs` |

`numpy` and `numba` are installed for you; nothing else is needed to compute a probability.  `scipy` is used by the test suite alone, to cross-check the evolution operator against an independent matrix exponential; the library itself never imports it.

`numba` became a real dependency in 1.13.0 rather than an optional extra, because it is worth roughly 1.5x to 20x on large scans and there is no reason for that to depend on reading far enough down this page.  The NumPy path remains, is still tested on every push, and is what runs if `numba` is unavailable or switched off — the results are identical to round-off either way.  **One consequence to know about:** `numba` declares its own ceiling on `numpy` (0.66.0 requires `numpy<2.5`), so installing into an environment that holds a newer `numpy` than that will downgrade it.  If you need the newest `numpy` more than you need the speed, install with `--no-deps` and add `numpy` yourself, or set `fastkernels.USE_NUMBA = False` and ignore the compiled path.


## Installation

**NuOscProbExact** is pure Python: there is nothing to compile or link.

> **Python version:** The code requires Python 3.9 or newer, and every release is tested on 3.9, 3.10, 3.11, 3.12, and 3.13.  The floor comes from `numpy.broadcast_shapes`, which the batched paths use and which arrived in NumPy 1.20; 3.9 is also the oldest version for which the optional `numba` backend still has a wheel.

### From PyPI (recommended)

```shell
pip install nuoscprobexact
```

That is the whole installation, and it includes the compiled backend: `numpy` and `numba` both arrive, and large scans run on compiled kernels without anything further being asked for.

The optional extras add what each task needs, and can be combined:

```shell
pip install "nuoscprobexact[fast]"       # numba; now a base dependency, so this is a no-op
pip install "nuoscprobexact[notebooks]"  # Jupyter, matplotlib and scipy, for notebooks/
pip install "nuoscprobexact[test]"       # pytest, scipy and coverage, to run the suite
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

The modules are installed under their bare names --- `oscprob2nu`, `oscprob3nu`, `oscprob4nu`, `hamiltonians2nu`, `hamiltonians3nu`, `hamiltonians4nu`, `globaldefs`, `fastkernels`, `slabs`, `earth` --- which is the same way the paper and the worked examples refer to them.

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
    │   ├── ADVERSARIAL.md           # What the adversarial check must try to break
    │   ├── OBJECTIONS.md            # Each objection, the measurement answering it, and its test
    │   ├── requirements.lock        # Exact Python versions the benchmark venv installs
    │   ├── build.sh                 # Clones, hash-verifies and builds all seven at their pins
    │   ├── bench.hpp                # The C++ harness: owns main(), every clock, the statistics
    │   ├── machine.py               # Environment capture, the canary, and the rejection rule
    │   ├── conversions.py           # The one place a physical conversion factor is computed
    │   ├── gen_conversions.py       # Emits conversions.h so no adapter carries a physical literal
    │   ├── adapters/                # One adapter per code, physics only
    │   │   ├── nufast_earth.cpp     # NuFast-Earth, batched over energy and zenith, invariants hoisted
    │   │   ├── nufast_lbl.cpp       # NuFast-LBL, one batched call over the whole energy vector
    │   │   ├── prob3.cpp            # Prob3++, looped in C++ because that is its interface
    │   │   ├── globes.cpp           # GLoBES, looped, chord decomposition hoisted
    │   │   ├── nuoscprobexact.py    # This library, batched
    │   │   ├── nusquids.py          # nuSQuIDS through its multiple-energy constructor
    │   │   └── nucraft.py           # nuCraft, whose batching is interface-only
    │   └── runner.py                # The Python harness; the only place Python timing happens
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
```

### Without installing anything

The three core modules are self-contained --- they need only `numpy` and the standard library --- so copying `src/oscprob2nu.py`, `src/oscprob3nu.py` or `src/oscprob4nu.py` into your own project is a supported way to use **NuOscProbExact**.  Each imports `fastkernels` if it is available and does without it if it is not, so a lone copy works and simply runs the NumPy path; a test copies each of the three out and exercises it that way.  Adding `src/` to the path works too, and is what the bundled examples do:

```python
import sys
sys.path.append('/path/to/NuOscProbExact/src')

import oscprob3nu
```

### Checking the installation

**Run the worked examples.**
   Inside the directory `examples/`, we provide several example files to get you started.  Each is runnable as it stands and prints the probabilities it computes; [Usage and examples](#usage-and-examples) below walks through the first of them.  To run any of the examples, just execute, *e.g.*,
   ```shell
   python example_2nu_trivial.py
   ```
   Inspecting the example files and reading their description below will help you to learn how to use **NuOscProbExact** in your own project.

   > **Renamed:** this directory was called `test/` in version 1.0.0 of the code, and is named that way in version 2 of [the paper](https://arxiv.org/abs/1904.12391).  It became `examples/` to stop it being confused with `tests/`, which holds the regression suite.

**Run the regression tests.**
   ```shell
   cd /path/to/NuOscProbExact
   pytest
   ```
   These check the SU(2), SU(3) and SU(4) machinery against independent computations --- unitarity of the evolution operator, agreement with `scipy.linalg.expm`, agreement with the standard oscillation formulas, and the sign conventions of the sample Hamiltonians --- and run every example embedded in the docstrings.

**Open the notebooks.**
   ```shell
   cd /path/to/NuOscProbExact
   pip install -e ".[notebooks]"
   jupyter lab notebooks/
   ```
   Twenty worked notebooks, numbered in reading order, covering the probabilities against baseline and against energy, matter and new physics, oscillograms, bi-probability plots, the Earth, arbitrary matter profiles, performance, the paper's own figures, the textbook approximations, mass ordering and the octant, antineutrinos, solar neutrinos, numerical edge cases, four-neutrino 3+1 scenarios, cross-checks with other public codes, the evolution operator itself, animated parameter sweeps, and a Hamiltonian of your own carried through a varying profile.  They carry their figures inline, so they can also just be read on GitHub.

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

### 2. Numba, which you already have

Nothing to install and nothing in your code to change: [Numba](https://numba.pydata.org) comes with the package as of 1.13.0.  The batched paths run as compiled machine-code loops spread over your cores instead of as a chain of NumPy array operations.  Where Numba is unavailable, or switched off with `fastkernels.USE_NUMBA = False`, the NumPy path is used and the results are the same to round-off.

Measured on 2000-point scans, against the equivalent Python loop:

| Scan | loop | arrays | arrays + Numba |
|---|---|---|---|
| Three-flavor, vs. baseline | 38 ms | 1.8 ms (~21×) | 0.31 ms (**~120×**) |
| Three-flavor, vs. energy | 34 ms | 1.5 ms (~23×) | 0.20 ms (**~170×**) |
| Three-flavor oscillogram, 100×100 | 197 ms | 5.3 ms (~37×) | 0.85 ms (**~230×**) |
| Two-flavor, vs. baseline | 6.9 ms | 0.07 ms (~99×) | *not used — see below* |

Best of seven runs, interleaved, on one machine.  These are indicative, not precise: repeated runs vary by tens of per cent, so treat them as orders of magnitude.  [Notebook 09](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/09_performance.ipynb) measures the same comparison on whatever machine runs it, which is the number to trust.

**Checking the input costs more than the arithmetic on a large scan.**  Every entry point verifies that the Hamiltonian is Hermitian, because one that is not returns probabilities that still sum to one — so nothing downstream reveals the mistake.  Validating a stack is a pass over it, the same order of work as evaluating it: 1.3× to 1.8× at two thousand points, and 3.2× to 5.7× at two hundred thousand, where the compiled kernel has made the evaluation fast enough that the check dominates.  If your Hamiltonians come from a construction you already trust — everything `hamiltonians2nu`, `hamiltonians3nu` and `hamiltonians4nu` build is Hermitian to round-off — decline it with `oscprob3nu.CHECK_HERMITICITY = False`, and likewise on the other two modules.

**Layered matter gets its own kernel.**  `slabs` and `earth` compose evolution operators rather than computing probabilities, so until 1.12.0 they could not use the backend at all and quietly ran the NumPy path.  They now compose the whole trajectory in one compiled pass — never materialising the stack, never dispatching a product per slab — which is worth ~12×, ~12× and ~9× on a single Earth crossing at two, three and four flavors.  The threshold there is one slab, not `MIN_BATCH`: the comparison is against a Python loop, so the kernel is ahead at every length.

**An Earth scan is one call.**  Both the energy and the zenith angle may be arrays, and they broadcast, so an oscillogram is `probabilities_3nu_earth(h, energies[None, :], costhz[:, None])` rather than a loop over its points.  The geometry and the matter potentials depend on the angle alone, so a scan builds them once instead of once per energy; the chord kernels then build each slab's Hamiltonian as they go, so the stack of operators that made the batched path memory-bound is never allocated.  Over two thousand energies:

| Earth scan, 2000 energies | NumPy loop | Compiled loop | Compiled array |
|---|---|---|---|
| Two flavors | 542 ms | 44 ms | **1.2 ms** |
| Three flavors | 895 ms | 77 ms | **7.5 ms** |
| Four flavors | 2224 ms | 272 ms | **25 ms** |

That is ~38×, ~10× and ~11× over the compiled loop, and ~460×, ~120× and ~87× over the NumPy one.  A 100 × 100 oscillogram at three flavors takes 42 ms against 371 ms for the loop, about 4 µs a point.

**A chord is a palindrome.**  A neutrino crossing a spherically symmetric Earth meets every radius on the way in and again on the way out, so slab *j* and slab *n−1−j* carry the same Hamiltonian, the same width, and therefore the same operator.  Composing from the centre outwards computes each of them once instead of twice — worth ~1.4×, ~1.5× and ~1.8× of the figures above.  Nothing about it is specific to PREM: `fastkernels.palindromic` decides, and `slabs` asks it of any sequence, so a symmetric castle wall is composed at the same discount.  `fastkernels.USE_PALINDROME = False` asks for the plain left-to-right product instead.

**Ask for an accuracy, not a slab count.**  `n_slabs_per_segment` fixes the discretisation rather than the error, and the two are not the same thing: the error is strongly energy- and angle-dependent.  Pass `rtol` or `atol` to any `earth` entry point and the chord is refined until the measured error meets it, on *every* returned probability, raising rather than silently returning something coarser when it cannot.  `earth.slabs_for_tolerance` runs the same search on its own, which is the cheap way to set one subdivision for a whole scan.  The same machinery serves profiles that are not the Earth through `slabs.probabilities_{2,3,4}nu_profile`, which take the profile as a callable — scaled by the baseline rather than by the last sampled position, which moves as the refinement doubles.

**The backend is not used where it would not help.**  For three flavors it wins at every stack size, by between two and sixteen times.  For two flavors it does not: that expansion reduces to a square root and a sine per element, which NumPy already does about as well as compiled code can, and the kernel additionally has to materialise the Hamiltonian stack.  Below fifty thousand elements the NumPy path is quicker, so it is kept; above, the kernel leads by about 1.3–1.8×.  The thresholds are measured, and the library picks whichever is faster without you doing anything.

Two costs, so the trade is visible: importing Numba takes about 140 ms against 65 ms for NumPy alone, and the first call into a kernel compiles, which takes a few seconds.  The kernels are cached on disk, so later runs start in milliseconds.  Neither cost falls on a single scalar probability, which never enters a kernel at all — it is the closed form, and it returns in microseconds whether or not Numba is present.

That is why this was an optional extra until 1.13.0, and the argument that changed: the costs are paid once and are invisible next to a 12× Earth crossing, while an extra is paid attention to only by readers who get this far.

### What you do not have to think about

* **Short stacks.** Below thirteen elements at three flavors, and twelve at two, the array machinery costs more than it saves, so those are evaluated one at a time automatically.
* **The scalar path.** It is deliberately left uncompiled: 8 µs is not worth a compilation pause on a first call.
* **Turning Numba off.** `fastkernels.USE_NUMBA = False` forces the NumPy path, which is how the test suite checks that the two agree.

One thing that *is* worth doing by hand: build the energy-independent part of the vacuum Hamiltonian once, outside any scan, since it does not depend on the energy.  The bundled examples all do this.

## Usage and examples

There are three core modules, one per flavor count: `oscprob2nu.py`, `oscprob3nu.py` and `oscprob4nu.py`.  Each is stand-alone apart from the dependencies described [above](#requirements).  Install the package, or add `src/` to the path, which is what the bundled examples do:

```python
import sys

sys.path.append('../src')
```

What follows is the short version: what the functions take and return, and four examples that between them cover a single probability, a whole scan, four flavors, and your own Hamiltonian.  Everything else --- matter, non-standard interactions, Lorentz-invariance violation, oscillograms, the Earth, the expansion coefficients --- lives in the runnable scripts in [`examples/`](https://github.com/mbustama/NuOscProbExact/tree/main/examples/) and in the [notebooks](#notebooks), which store their figures inline and are executed by CI.  It is not repeated here, so there is one copy of each to keep correct.


### Basics

Most of the time you want probabilities, not the intermediate steps.  The routine is `probabilities_Nnu` in `oscprobNnu`, and it takes a Hermitian matrix and a baseline:

| Flavors | Call | Returns |
|---|---|---|
| 2 | `oscprob2nu.probabilities_2nu(h, L)` | 4 values: `Pee, Pem, Pme, Pmm` |
| 3 | `oscprob3nu.probabilities_3nu(h, L)` | 9 values: `Pee, Pem, Pet, Pme, ..., Ptt` |
| 4 | `oscprob4nu.probabilities_4nu(h, L)` | 16 values: `Pee, Pem, Pet, Pes, Pme, ..., Pss` |

In every case the initial flavor varies slowest, so `P[n*alpha + beta]` is P(nu_alpha -> nu_beta) for `n` flavors.  The two-flavor labels could equally be `Pmm, Pmt, Ptm, Ptt` --- which pair of flavors they describe is set by the Hamiltonian you pass, not by the code.

The evolution operator itself is available too, as `evolution_operator_Nnu(h, L)`, if you want to compose it across segments or propagate a density matrix rather than read off probabilities.

> **Important:** The Hamiltonian must be Hermitian, and every entry point checks that it is: one that is not raises `ValueError` rather than returning numbers.  The check is there because the numbers it would otherwise return still sum to one, so nothing downstream would reveal the mistake.  [Performance](#performance) gives what the check costs and how to decline it where your Hamiltonians are Hermitian by construction.

> **About the units:** These modules assume no units for any of the model parameters, so you need to pass values with consistent ones --- all that is required is that `H*L` be dimensionless.  The module `globaldefs` provides physical constants and conversion factors, including `CONV_KM_TO_INV_EV`, which converts a baseline in km to eV^{-1}.


### A first probability

Three-flavor oscillations in vacuum, at a fixed energy and baseline.  `hamiltonian_3nu_vacuum_energy_independent` returns the vacuum Hamiltonian **without** the *1/E* prefactor, so that a scan over energies computes it once and divides by a varying *E*:

```python
# Find this example in NuOscProbExact/examples/example_3nu_vacuum.py

import numpy as np

import oscprob3nu
import hamiltonians3nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

# NuFit best-fit mixing parameters, pulled from globaldefs.  NO means
# "normal ordering"; change NO to IO for inverted ordering
h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
    S12_NO_BF, S23_NO_BF, S13_NO_BF, DCP_NO_BF, D21_NO_BF, D31_NO_BF)
h_vacuum = np.asarray(h_vacuum_energy_indep)/energy

Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = oscprob3nu.probabilities_3nu(
    h_vacuum, baseline*CONV_KM_TO_INV_EV)

print("Pee = %6.5f, Pem = %6.5f, Pet = %6.5f" % (Pee, Pem, Pet))
print("Pme = %6.5f, Pmm = %6.5f, Pmt = %6.5f" % (Pme, Pmm, Pmt))
print("Pte = %6.5f, Ptm = %6.5f, Ptt = %6.5f" % (Pte, Ptm, Ptt))
```

This returns

```shell
Pee = 0.92768, Pem = 0.01432, Pet = 0.05800
Pme = 0.04023, Pmm = 0.37887, Pmt = 0.58090
Pte = 0.03210, Ptm = 0.60680, Ptt = 0.36110
```

Each row sums to one, as it must.

> **Antineutrinos:** pass `-dCP` instead of `dCP`, and flip the sign of the matter potential.  [Notebook 13](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/13_antineutrinos.ipynb) works through both, and the two ways to get it wrong.


### Whole scans in one call

Do not call the routine in a Python loop.  Every core routine accepts a **stack** of Hamiltonians, an array of baselines, or both broadcast against each other, and evaluates the lot in one call --- which is [tens of times faster](#performance) and gives identical results:

```python
import numpy as np

import oscprob3nu
import hamiltonians3nu
from globaldefs import *

baseline = 1.3e3*CONV_KM_TO_INV_EV       # [eV^{-1}]
energies = np.logspace(-1.0, 1.0, 200)*1.e9   # 0.1 to 10 GeV [eV]

h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
    S12_NO_BF, S23_NO_BF, S13_NO_BF, DCP_NO_BF, D21_NO_BF, D31_NO_BF)

# One Hamiltonian per energy, stacked along a leading axis
h_stack = np.asarray(h_vacuum_energy_indep)/energies[:, None, None]

prob = oscprob3nu.probabilities_3nu(h_stack, baseline)   # shape (200, 9)
prob_ee, prob_em, prob_et = prob[:, 0], prob[:, 1], prob[:, 2]

print("prob.shape =", prob.shape)
print("P_ee at %5.2f GeV = %.5f" % (energies[0]/1.e9, prob_ee[0]))
print("P_ee at %5.2f GeV = %.5f" % (energies[-1]/1.e9, prob_ee[-1]))
```

This returns

```shell
prob.shape = (200, 9)
P_ee at  0.10 GeV = 0.24693
P_ee at 10.00 GeV = 0.98582
```

The same works for a scan over baselines, or for both at once to build an oscillogram.  The sample Hamiltonians in matter accept an array of energies directly, so a matter scan is two calls and no loop.  [Notebook 02](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/02_vacuum_oscillations.ipynb) plots these curves, [notebook 04](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/04_oscillogram.ipynb) builds an oscillogram, and [notebook 09](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/09_performance.ipynb) measures what the broadcasting buys on whatever machine runs it.


### Four flavors: a 3+1 sterile state

`oscprob4nu` works exactly the same way, with a 4x4 Hamiltonian and sixteen probabilities.  With the fourth state read as sterile, the flavor order is (nu_e, nu_mu, nu_tau, nu_s):

```python
import numpy as np

import oscprob4nu
import hamiltonians4nu
from globaldefs import *

# Three extra mixing angles and one extra mass-squared splitting,
# here Dm41^2 = 1 eV^2
h_vacuum_energy_indep = hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(
    S12_NO_BF, S23_NO_BF, S13_NO_BF,
    np.sqrt(0.10), np.sqrt(0.10), 0.0,
    DCP_NO_BF, D21_NO_BF, D31_NO_BF, 1.0)

prob = oscprob4nu.probabilities_4nu(
    np.asarray(h_vacuum_energy_indep)/1.e9, 1.3e3*CONV_KM_TO_INV_EV)

print("%d probabilities" % len(prob))
print("Pee = %6.5f, Pem = %6.5f, Pet = %6.5f, Pes = %6.5f" % tuple(prob[0:4]))
print("they sum to %.5f" % sum(prob[0:4]))
```

```shell
16 probabilities
Pee = 0.76700, Pem = 0.00149, Pet = 0.05220, Pes = 0.17931
they sum to 1.00000
```

[Notebook 16](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/16_four_neutrinos.ipynb) works a 3+1 scenario through properly --- the sterile entry in the matter potential, a short-baseline scan, the sterile matter resonance through the Earth --- and explains why four flavors is where the closed form ends.


### Arbitrary Hamiltonians

Nothing above is a special case in the code: vacuum, matter, non-standard interactions and Lorentz-invariance violation are each just a different Hermitian matrix handed to the same routine.  So your own model is too.  Usually you will want to add a term to the vacuum Hamiltonian, where `hamiltonian_mymodel` is yours to write and returns a 3x3 matrix:

```python
import numpy as np

import oscprob3nu
import hamiltonians3nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
    S12_NO_BF, S23_NO_BF, S13_NO_BF, DCP_NO_BF, D21_NO_BF, D31_NO_BF)
h_vacuum = np.asarray(h_vacuum_energy_indep)/energy

h_mymodel = h_vacuum + hamiltonian_mymodel(mymodel_parameters)

prob = oscprob3nu.probabilities_3nu(h_mymodel, baseline*CONV_KM_TO_INV_EV)
```

`hamiltonian_mymodel` may depend on the energy too.  For two flavors it returns a 2x2 matrix instead, and for four, a 4x4 one.  Passing an arbitrary matrix directly, with no vacuum term at all, works exactly as you would expect --- see [`examples/example_3nu_trivial.py`](https://github.com/mbustama/NuOscProbExact/blob/main/examples/example_3nu_trivial.py).


### Where the rest is

Each of these is a runnable script; none of them is transcribed into this file, so there is a single copy to keep correct.

| Scenario | Script | Notebook |
|---|---|---|
| Arbitrary Hamiltonian, 2 and 3 flavors | [`example_2nu_trivial.py`](https://github.com/mbustama/NuOscProbExact/blob/main/examples/example_2nu_trivial.py), [`example_3nu_trivial.py`](https://github.com/mbustama/NuOscProbExact/blob/main/examples/example_3nu_trivial.py) | [01](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/01_basics.ipynb) |
| Vacuum, 2 and 3 flavors | [`example_2nu_vacuum.py`](https://github.com/mbustama/NuOscProbExact/blob/main/examples/example_2nu_vacuum.py), [`example_3nu_vacuum.py`](https://github.com/mbustama/NuOscProbExact/blob/main/examples/example_3nu_vacuum.py) | [02](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/02_vacuum_oscillations.ipynb) |
| Constant-density matter | [`example_3nu_matter.py`](https://github.com/mbustama/NuOscProbExact/blob/main/examples/example_3nu_matter.py) | [03](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/03_matter_nsi_liv.ipynb) |
| Matter with non-standard interactions | [`example_3nu_nsi.py`](https://github.com/mbustama/NuOscProbExact/blob/main/examples/example_3nu_nsi.py) | [03](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/03_matter_nsi_liv.ipynb) |
| Lorentz-invariance violation | [`example_3nu_liv.py`](https://github.com/mbustama/NuOscProbExact/blob/main/examples/example_3nu_liv.py) | [03](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/03_matter_nsi_liv.ipynb) |
| SU(2) and SU(3) expansion coefficients, and the evolution operator | [`example_2nu_vacuum_coeffs.py`](https://github.com/mbustama/NuOscProbExact/blob/main/examples/example_2nu_vacuum_coeffs.py), [`example_3nu_vacuum_coeffs.py`](https://github.com/mbustama/NuOscProbExact/blob/main/examples/example_3nu_vacuum_coeffs.py) | [18](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/18_evolution_operator.ipynb) |
| Layered matter, and the Earth through PREM | — | [06](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/06_earth_and_prem.ipynb), [07](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/07_earth_probabilities.ipynb) |

The [numerical recipes](https://mbustama.github.io/NuOscProbExact/recipes.html) page collects the same material as runnable snippets, the [tutorial notebooks](https://mbustama.github.io/NuOscProbExact/notebooks.html) page indexes the twenty in reading order, and the [API reference](https://mbustama.github.io/NuOscProbExact/functions.html) documents every routine, with examples that are executed when the documentation is built rather than pasted beside it.


## Notebooks

Twenty worked notebooks live in [`notebooks/`](https://github.com/mbustama/NuOscProbExact/tree/main/notebooks/), numbered in reading order.  They carry their figures inline, so they render on GitHub without being run.  The documentation lists them with the same descriptions, grouped by what they are for, on its [tutorial notebooks](https://mbustama.github.io/NuOscProbExact/notebooks.html) page:

| Notebook | What it covers |
|---|---|
| [01 Basics](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/01_basics.ipynb) | Units, one probability, and why to pass arrays rather than loop |
| [02 Vacuum oscillations](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/02_vacuum_oscillations.ipynb) | Against baseline and against energy, checked against the textbook formula |
| [03 Matter, NSI, LIV](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/03_matter_nsi_liv.ipynb) | Constant-density matter and two kinds of new physics |
| [04 Oscillograms](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/04_oscillogram.ipynb) | A 240x240 energy-baseline map in a single call |
| [05 Bi-probability](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/05_biprobability.ipynb) | CP ellipses, in vacuum and in matter |
| [06 The Earth and PREM](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/06_earth_and_prem.ipynb) | The density profile, chord geometry, slabs, and their convergence |
| [07 Through the Earth](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/07_earth_probabilities.ipynb) | Zenith-angle scans, an Earth oscillogram, and real baselines |
| [08 Unusual density profiles](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/08_unusual_density_profiles.ipynb) | Castle-wall and serrated profiles, and parametric enhancement |
| [09 Performance](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/09_performance.ipynb) | Looping versus broadcasting, and the compiled backend, measured live |
| [10 The paper's figures](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/10_paper_figures.ipynb) | The two figures from [arXiv:1904.12391](https://arxiv.org/abs/1904.12391) |
| [11 Exact vs approximations](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/11_exact_vs_approximations.ipynb) | Where the familiar formulas agree, and where they do not |
| [12 Ordering and octant](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/12_ordering_and_octant.ipynb) | Normal against inverted, and the θ₂₃ octant degeneracy |
| [13 Antineutrinos](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/13_antineutrinos.ipynb) | Conjugate *and* flip the potential — and two ways to get it wrong |
| [14 Solar and the MSW resonance](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/14_solar_and_adiabatic_msw.ipynb) | The adiabatic resonance, validated — and why slabs are the wrong tool for it |
| [15 Numerical edge cases](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/15_numerical_edge_cases.ipynb) | Degenerate spectra, and what returns a number instead of NaN |
| [16 Four neutrinos](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/16_four_neutrinos.ipynb) | A 3+1 sterile state through SU(4), and why the method stops at four |
| [17 Cross-checks with other codes](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/17_cross_checks.ipynb) | Corroboration from nuSQuIDS and from Zaglauer–Schwarzer, and the conventions that had to be matched |
| [18 The evolution operator](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/18_evolution_operator.ipynb) | The operator itself, the group property, and the SU(*n*) coefficients underneath |
| [19 Animated scenes](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/19_animations.ipynb) | Four parameter sweeps drawn as stills, and the reel they came from |
| [20 An arbitrary Hamiltonian](https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/20_arbitrary_hamiltonian.ipynb) | A Hamiltonian of your own carried through three varying profiles, and how to put one through the Earth |

Run them with `pip install -e ".[notebooks]"` and `jupyter lab notebooks/`.  Every one of them is executed by CI, so an example that stops working fails the build.


## Documentation and help

All of the modules provided in **NuOscProbExact** have been documented using Python docstrings, written in [numpydoc](https://numpydoc.readthedocs.io/) format so that they can be rendered directly by [Sphinx](https://www.sphinx-doc.org/) with the `numpydoc` extension.  They are human-readable by opening the source `.py` files.  Alternatively, they can be printed from within an interactive Python session.

Every `Examples` block in the docstrings is executed when the documentation is built, so the results shown on the [API page](https://mbustama.github.io/NuOscProbExact/functions.html) are produced by the code rather than pasted beside it, and cannot drift.  The regression suite runs the same blocks on every supported Python (`tests/test_docstrings.py`), which the documentation build --- one job, one interpreter --- would not catch.

A full Sphinx project lives in `docs/`.  Build it with
```shell
pip install -r docs/requirements.txt
cd docs && make html
```
and open `docs/build/html/index.html`.  It contains an installation guide, a quickstart, a description of the [method](https://github.com/mbustama/NuOscProbExact/blob/main/docs/source/methodology.rst) and its sign conventions, the API reference generated from the docstrings, a bibliography, and the changelog.

Notable changes between versions are recorded in [CHANGELOG.md](https://github.com/mbustama/NuOscProbExact/blob/main/CHANGELOG.md), which the documentation renders as its own page, so there is a single source of truth.

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

If you use **NuOscProbExact** in your work, we ask you that you please cite the following paper: Mauricio Bustamante, *NuOscProbExact: a general-purpose code to compute exact two-flavor and three-flavor neutrino oscillation probabilities* ([arXiv:1904.12391](https://arxiv.org/abs/1904.12391)).  The papers that cite it are listed [on INSPIRE](https://inspirehep.net/literature?sort=mostrecent&size=250&page=1&q=refersto%3Arecid%3A1731803&ui-citation-summary=true).

If you are citing **NuOscProbExact** in a document that will be uploaded to the arXiv, please consider using the LaTeX or BibTeX entries provided by INSPIRE ([link here](https://inspirehep.net/literature/1731803)):
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

The paper covers two and three flavors, which is what it was written about.  To cite *the software* — a particular version of it, including the four-flavor extension that came after the paper — use the Zenodo DOI badge at the top of this file, which resolves to the most recently archived release.  Zenodo mints a DOI per GitHub Release, so citing a specific version means citing the DOI archived for it.


## License

**NuOscProbExact** is released under the [MIT License](https://opensource.org/licenses/MIT).  The full text ships with the source, as [`LICENSE`](https://github.com/mbustama/NuOscProbExact/blob/main/LICENSE) in the repository root.

