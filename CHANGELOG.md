# Changelog

All notable changes to **NuOscProbExact** are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project uses [Semantic Versioning](https://semver.org/).

## [1.2.0] - 2026-07-31

Performance.  The exact SU(2) and SU(3) expansions now evaluate a whole stack
of Hamiltonians and baselines in one pass, which is what a scan over energy or
baseline, or an oscillogram, actually needs.

Nothing about the method changes, and nothing about the existing interface
changes: a single Hamiltonian with a scalar baseline returns exactly what it
returned before, at the same speed.

### Added

- `evolution_operator_2nu`, `evolution_operator_3nu`, `probabilities_2nu` and
  `probabilities_3nu` accept a stack of Hamiltonians of shape `(..., n, n)`, an
  array of baselines, or both, broadcast against each other.  The result is an
  array with the broadcast leading axes: `(..., 9)` for the three-flavor
  probabilities, `(..., 3, 3)` for the evolution operator.  Measured against
  the equivalent Python loop over 2000 points:

  | Scan | Speedup |
  |---|---|
  | Versus baseline: one Hamiltonian, many baselines | 80-100x |
  | Versus energy: many Hamiltonians, one baseline | 20-30x |
  | Oscillogram, 100 energies x 100 baselines | ~80x |
  | Two flavors, versus baseline | ~50x |

  The two scans differ because the latent roots depend on the Hamiltonian
  alone: scanning one Hamiltonian over many baselines solves the
  characteristic equation once and then evaluates only the phases, whereas an
  energy scan changes the Hamiltonian at every point.  The vectorised path
  agrees with the scalar one to 1.2e-13, and with `scipy.linalg.expm` to
  3.2e-14.

  For context, the vectorised expansion is now about as fast as diagonalising
  with LAPACK --- one `eigh` plus phases for a baseline scan, batched `eigh`
  for an energy scan.  That is the honest comparison to draw: the SU(3) route's
  advantage is that it is a closed form, not that it outruns `eigh`.  What the
  vectorisation removes is the disadvantage it used to carry.

- 24 tests covering the vectorised path (`tests/test_vectorized.py`): agreement
  with the scalar path element by element, all three broadcasting patterns,
  unitarity and normalization on the batched path, empty and length-one stacks,
  and degenerate Hamiltonians mixed into an otherwise ordinary stack.

- A section on scanning in `docs/source/quickstart.rst`, and a revised cost
  section in `docs/source/methodology.rst` carrying the measurements above.

### Changed

- Degenerate spectra are still handled exactly on the vectorised path, but the
  branch cannot be taken elementwise inside a vectorised expression.  The
  general formula is evaluated everywhere with vanishing denominators replaced
  by one, and the affected elements are then recomputed individually with the
  scalar routine.  Degeneracy is measure-zero among floating-point
  Hamiltonians, so that fallback loop is empty in essentially every real use;
  the tests exercise it by constructing degenerate cases deliberately.

- The scalar path is unchanged in speed.  Dispatching between the scalar and
  vectorised paths is done by inspecting the argument types rather than by
  calling `numpy.ndim`, which would convert a nested list to an array on every
  call: the first implementation did exactly that and cost the scalar path 47%,
  which is why the dispatch is written the way it is.  Measured best-of-five,
  interleaved against the previous release, the overhead is now -1%.

- The scalar three-flavor routine is refactored so that its core takes the
  SU(3) coefficients and invariants rather than the Hamiltonian matrix, since
  the vectorised path reuses it for degenerate elements.

## [1.1.0] - 2026-07-31

An audit of the code released alongside
[arXiv:1904.12391](https://arxiv.org/abs/1904.12391), covering the mathematics,
the implementation, and the documentation.

The exact SU(2) and SU(3) machinery at the heart of the method was found to be
correct: the evolution operators agree with `scipy.linalg.expm` to 8e-15 over
200 random Hermitian Hamiltonians, they are unitary to 5e-15, and all 512
entries of the hard-coded `d` tensor reproduce
d<sub>ijk</sub> = ¼Tr({λ<sub>i</sub>,λ<sub>j</sub>}λ<sub>k</sub>) to 2e-16. The
defects were in the layers around it: a shortcut expression for the two-flavor
probability, the sign convention of the two-flavor vacuum Hamiltonian, and
essentially all of the validation and documentation infrastructure.

**Two of the fixes below change published numbers.** If you have results from a
previous version, the ones to re-check are two-flavor probabilities for any
Hamiltonian with a complex off-diagonal entry, and *any* two-flavor computation
in matter, with NSI, or with LIV. Three-flavor results are unaffected.

### Added

- A regression test suite in `tests/`, run with `pytest`: 132 tests covering the
  SU(3) algebra, both evolution operators, the probabilities, the sample
  Hamiltonians, the standard oscillation formulas, degenerate Hamiltonians, and
  the docstring examples. Previously nothing in the repository checked a
  computed probability against an independent calculation --- `test/` holds
  worked examples and plot generators, but not a single assertion. Every defect
  fixed in this release has a test that fails without the fix.
- Handling of degenerate Hamiltonians in the SU(3) expansion. The coefficients
  u<sub>k</sub> come from Lagrange interpolation over the three latent roots,
  which divides by 3ψ<sub>m</sub>² − |h|², the derivative of the characteristic
  polynomial; that factor vanishes exactly at a repeated root. Two cases are now
  handled explicitly and exactly: |h|² = 0, where the Hamiltonian is
  proportional to the identity and U₃ = 𝟙, and a doubly degenerate root, where
  the confluent limit of the same expansion collapses the spectral decomposition
  onto a single projector. `diag(1, 1, −2)` now agrees with `scipy.linalg.expm`
  to 3.5e-15 where it previously returned NaN. The SU(2) counterpart takes the
  limit sin(|h|L)/|h| → L.
- `pyproject.toml`, so the package can be installed with `pip install -e .`.
  It keeps the flat module layout that the examples and the paper both assume,
  so `import oscprob3nu` works unchanged whether or not you install it. Optional
  extras separate matplotlib (figures), pytest and scipy (tests), and sphinx and
  numpydoc (documentation).
- A `Raises` contract on `tensor_d`, which now validates its indices.
- This changelog, rendered into the documentation by `docs/source/changelog.rst`.

### Changed

- The standard oscillation formulas `probabilities_2nu_vacuum_std`,
  `probabilities_2nu_matter_std` and `probabilities_3nu_vacuum_std` now take the
  neutrino energy in eV and the baseline in eV⁻¹, like the rest of the library.
  **This changes the signature of three routines.** They previously took GeV and
  km and folded the conversion into the rounded constants 1.27 and 2.54. Since
  these routines exist solely to validate the exact computation, the rounding
  mattered: the exact prefactor implied by this repository's own
  `CONV_KM_TO_INV_EV` is 1.266933…, so 1.27 overstates every oscillation phase
  by 0.242%, and near an oscillation minimum that becomes a large relative error
  in the probability. At L = 1000 km, E = 1 GeV the reference formula returned
  0.004125 against the exact 0.003204 --- a 29% discrepancy in the quantity
  meant to confirm the exact result. Agreement is now 4e-19.
- `hamiltonian_2nu_coefficients` and `hamiltonian_3nu_coefficients` return real
  floats. The coefficients h<sub>k</sub> of a Hermitian Hamiltonian are real by
  construction, but the routines returned a mixture: h₁ and h₂ as floats, h₃ and
  h₈ as complex, because they were built from arithmetic on the complex diagonal
  entries. The docstrings claimed the coefficients "are complex numbers, in
  general", which is not true for a Hermitian Hamiltonian.
- The sample Hamiltonians are returned as complex `numpy.ndarray` throughout,
  rather than a mixture of nested lists and real arrays.
- The `d` tensor is tabulated once at import time as a dense 8×8×8 array, and
  the star product is computed once per call instead of once per (k, m) pair. A
  single `probabilities_3nu` call used to evaluate `tensor_d` 2048 times --- 512
  from the 8³ loop for ⟨h⟩, and 1536 because `star(k, h_coeffs)` sat inside the
  loop over m, rebuilding a baseline-independent vector three times for each of
  eight k. None remain on the hot path, and an evaluation goes from roughly 700
  to roughly 40 microseconds. The algebra is unchanged, to 1e-14.
  *(Corrected in 1.2.0: the commit message for this change, and the description
  of pull request #2, quote "757 to 308 microseconds". That measurement was
  taken without a warm-up and is wrong; re-measured carefully, best of five
  runs, the figures are ~700 µs before and ~40 µs after --- a factor of ~17
  rather than ~2.5.)*
- The library docstrings are rewritten in numpydoc format, so Sphinx can render
  them through `sphinx.ext.napoleon` and the `numpydoc` extension; a trial
  autodoc build completes with no warnings. Each module gains `__all__`, an
  explicit units section, and --- for the vacuum Hamiltonians --- a statement of
  the sign convention and why it matters.
- Every worked example in the docstrings is now executable and is run as a
  doctest by `tests/test_docstrings.py`, so the numbers quoted in the
  documentation cannot drift from what the code returns.
- The README's quoted output for the two-flavor trivial example changes from
  `Pee = 0.93213, Pem = 0.06787` to `Pee = 0.66063, Pem = 0.33937`. That example
  uses H₁₂ = 1 + 2i, so it was itself an instance of the missing h₂ contribution
  below: the README was documenting the bug. The two-flavor coefficient example
  changes sign, following the vacuum-Hamiltonian correction. All three-flavor
  outputs were already correct and are unchanged.
- `np.matrix.transpose` is replaced by the `@` operator and `.T` / `.conj().T`.
  The three call sites transposed an ndarray through an unbound method of
  `np.matrix`, which has carried a `PendingDeprecationWarning` since NumPy 1.15.

### Fixed

- **`probabilities_2nu` dropped the h₂ contribution.** It computed
  P<sub>eμ</sub> = |h₁|²/|h|² sin²(|h|L), but the transition probability is
  |U<sub>μe</sub>|² = u₁² + u₂², so the |h₂|² term is also required. Since
  h₂ = −Im(H₁₂) vanishes whenever the off-diagonal entry is real, oscillations in
  vacuum and in matter of constant density were unaffected and the error went
  unnoticed. It affected every Hamiltonian with a complex off-diagonal entry:
  NSI with a complex ε<sub>eμ</sub>, and any CP-violating two-flavor scenario.
  The function's own docstring example is such a case, documenting
  P<sub>eμ</sub> = 0.495179 while the code returned exactly 0.
- **The two-flavor vacuum Hamiltonian used the opposite sign convention.** It was
  built from M² = diag(Δm², −Δm²), assigning the larger mass-squared value to the
  first mass eigenstate, which yields the negative of the textbook Hamiltonian.
  In vacuum this is invisible --- for a real Hamiltonian the probabilities are
  invariant under H → −H --- but it stops being invisible the moment
  `hamiltonian_2nu_matter` adds the matter potential to the *ee* entry, because
  that flips the sign of the potential relative to the vacuum term. The result
  satisfied, exactly, P[H<sub>code</sub> + V] = P[H<sub>textbook</sub> − V]: every
  two-flavor computation in matter, with NSI, or with LIV returned the
  **antineutrino** probability when asked for the neutrino one, and the
  Mikheyev-Smirnov-Wolfenstein resonance sat on the wrong side. For θ₁₂ at a
  3000 km baseline it appeared at 0.011 GeV instead of 0.162 GeV. The
  three-flavor Hamiltonian already used the correct convention.
- cos ξ in `hamiltonian_2nu_liv` was computed as `sqrt(1 - sxi - sxi)` instead of
  `sqrt(1 - sxi*sxi)`. For 0 < sin ξ < ½ the LIV term is not a rotation of
  diag(b₁, b₂) at all --- at sin ξ = 0.3 it gives sin² + cos² = 0.49 --- and for
  sin ξ ≥ ½ the square root turns negative and the whole Hamiltonian becomes NaN.
  `globaldefs` sets `SXI12 = 0`, the one value at which the typo is harmless,
  which is why every bundled example survived it.
- `hamiltonian_2nu_nsi` silently discarded the imaginary part of a complex
  ε<sub>eμ</sub>. The two-flavor vacuum Hamiltonian is real, so the array was
  `float64` and the in-place addition truncated the value, emitting only a
  `ComplexWarning`.
- `probabilities_2nu_matter_std` computed cos 2θ as sqrt(1 − sin²2θ), discarding
  its sign. For θ > π/4 --- which includes the best-fit θ₂₃, whose cos 2θ is
  −0.164 --- the matter resonance was placed on the wrong side.
- `psi_roots` evaluated `pow(h2, -1.5)`, which is infinite when the traceless part
  of the Hamiltonian vanishes; every latent root, the evolution operator, and all
  nine probabilities came back NaN. It also passed the arc-cosine argument
  unclipped, so round-off could push it outside [−1, 1] and `cmath.acos` would
  quietly return a complex angle, giving complex roots and a non-unitary
  evolution operator.
- `tensor_d` dispatched through a chain of `elif` branches with no final `else`,
  so any index outside 0–7 fell off the end and the function returned `None`. The
  failure surfaced far from its cause, as a `TypeError` in whatever arithmetic
  consumed the result.
- **`run_testsuite.py`, the entry point documented in the README, crashed on its
  first plot** with `NameError: name 'pylab' is not defined`. The four plotting
  modules call `pylab.savefig`, but import pylab with `from pylab import *`,
  which binds the names inside pylab and never `pylab` itself. pyplot is now
  imported explicitly. The suite runs end to end and writes 42 figures.
- `plot_probability_3nu_vacuum_vs_l_std` called `probabilities_3nu_std`, which is
  defined nowhere in the repository, so the routine raised `NameError` on every
  invocation. Its baseline conversion was also inverted, dividing by
  `CONV_KM_TO_INV_EV` where it should multiply.
- The commented-out cross-checks against the standard formulas in the two
  paper-figure scripts referenced `D21_BF` and `D31_BF`, which do not exist, and
  converted their baselines the wrong way. They are corrected but left commented,
  since enabling them changes published figures; the cross-check itself now lives
  in `tests/test_reference_formulas.py`.
- Documentation errors: every worked example in the docstrings was stale
  (`hamiltonian_2nu_coefficients` documented `[2j, 0.0, -1]` and returns
  `[0.0, -2.0, -1]`; `hamiltonian_3nu_coefficients` documented h₈ = −1.732 and
  returns 4.041; `evolution_operator_2nu` printed a 3×3 matrix copied from the
  three-flavor module; all nine probabilities of `probabilities_3nu` were wrong).
  The placeholder `arXiv:1904.XXXXX` survived in eleven files. The routine listing
  of `hamiltonians3nu` named the two-flavor routines. `globaldefs` documented the
  Fermi constant in eV⁻¹ rather than eV⁻², described `CONV_EV_TO_G` as converting
  eV⁻¹ to grams, and labelled several inverted-ordering constants as normal
  ordering. `tensor_d` cited a reference that lived in the module docstring, which
  numpydoc scopes per docstring, leaving the citation dangling.

### Removed

- `from numpy import *` from the five library modules, where it shadowed the
  builtin `sum`, `abs` and `min`, and supplied the deprecated `np.matrix` under
  the bare name `matrix`.
- The duplicated `import cmath as cmath` alongside `import cmath`, the unused
  `import oscprob3nu` in `hamiltonians2nu`, and the now-unnecessary
  `copy.deepcopy` in the six Hamiltonian builders, which preceded an
  `np.multiply` that already allocated a new array.
- The module-level `sigma_1`, `sigma_2`, `sigma_3`, `sigma`, `identity` and
  `base` lists in `oscprob2nu`, which no routine referenced.
- The stale "In the works: optimizing the code to run faster, using Cython" note
  and the Python 2 compatibility note from the README. The code requires Python
  3.7 or newer and uses the `@` operator.

## [1.0.0] - 2019-04-30

Initial release, accompanying
*NuOscProbExact: a general-purpose code to compute exact two-flavor and
three-flavor neutrino oscillation probabilities*
([arXiv:1904.12391](https://arxiv.org/abs/1904.12391)).

### Added

- `oscprob2nu`: exact two-flavor oscillation probabilities for an arbitrary
  time-independent Hermitian 2×2 Hamiltonian, via the SU(2) exponential
  expansion.
- `oscprob3nu`: exact three-flavor oscillation probabilities for an arbitrary
  time-independent Hermitian 3×3 Hamiltonian, via the SU(3) exponential
  expansion, using the `d` tensor of the SU(3) algebra and the latent roots of
  the characteristic equation.
- `hamiltonians2nu` and `hamiltonians3nu`: sample Hamiltonians for oscillations
  in vacuum, in matter of constant density, in matter with non-standard
  interactions, and in a CPT-odd Lorentz invariance-violating background,
  together with the standard oscillation formulas.
- `globaldefs`: physical constants, unit-conversion factors, and the NuFit 4.0
  best-fit oscillation parameters for both mass orderings.
- Worked examples in `test/`, and `run_testsuite.py` to generate the figures.
