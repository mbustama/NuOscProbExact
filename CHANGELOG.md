# Changelog

All notable changes to **NuOscProbExact** are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed

- **The Earth speed-accuracy figures are re-measured against the compiled
  slab path.**  Their timings were taken before the evolution-operator
  kernel existed, when the Earth routines still evaluated one energy at a
  time and composed the slabs in a serial Python loop.  Both figures now
  call the batched path: at three flavors the same accuracy that cost
  9010 us costs 119 us, and this library moves from the right-hand edge of
  the plane to the lower-left corner, ahead of GLoBES and NuFast-Earth at
  every accuracy the six codes share.  The accuracies are unchanged to the
  last digit, and the arbitrary-precision referee moved by 7e-18, which is
  the check that the speed-up left the physics alone.

- **The 3+1 panel now carries one curve per `oscprob4nu.ROOT_STRATEGY`.**
  On this problem the two routes are indistinguishable: bit-identical
  probabilities at every slab count, and costs within 0.4% of each other.
  Two things had to be fixed before that could be said. Timing the routes
  in blocks made the eigensolver look 1.4x slower, which was the ordering
  and not the code; they are now interleaved with the sweep direction
  reversed on alternate passes. And single-shot timing of a call that can
  finish in a few hundred microseconds produced a three-flavor curve whose
  cost *fell* between 64 and 128 slabs per segment, so the timer now
  autoranges the way `timeit` does.

### Changed

- **The figures are recoloured and rearranged.**  `NuOscProbExact` is red
  across all three, having swapped with `nuCraft`, and carries the same open
  circle everywhere so that one code is followed by eye between the
  constant-density and Earth planes.  The Earth panels' y-axis label now fits
  on one line.  The four-flavor panel draws one curve rather than two: both
  root strategies are still measured and frozen, but they agree to the last
  bit and cost the same to within 0.4%, so two labels sat on one line.

- **`performance.pdf` gains `GLoBES` and `Prob3++`, and `nuSQuIDS` reaches
  30 000 energies.**  It also drops to two curves for this library, the
  compiled kernel and one point at a time, with the plain array path between
  them left out; and the legend moves above the panel in two columns to make
  room.  The new codes are flat lines, which is the point: none of the four
  external codes has a batched entry point over energies, and `nuSQuIDS`'
  multiple-energy mode amortises building the solver rather than the
  integration.  `nuCraft` is not there because it has no constant-density
  mode to time.

- **`nuCraft` now appears in the 3+1 Earth panel**, flat at 2.8e-3.  It was
  left out on the grounds that publishing that floor would misrepresent its
  solver; drawn with a caption that says where the floor comes from -- two
  independently rounded constants whose ratio is 0.5016 against an exact
  0.5 -- it misrepresents nothing, and it is what a user of released nuCraft
  gets.

## [1.13.0] - 2026-08-06

### Changed

- **`numba` is a dependency, not an optional extra.**  A plain
  `pip install nuoscprobexact` now gets the compiled path.  It was in the
  `fast` extra, which meant the 12x on an Earth crossing was collected by
  whoever read far enough down the README — and this release is largely about
  making the Earth quick, so leaving the speed opt-in wasted most of it.

  Nothing in the code changed.  `HAVE_NUMBA`, `USE_NUMBA` and `worthwhile()`
  already preferred the compiled path whenever it existed, so this is
  packaging alone.  `fast = ["numba"]` is kept as a no-op, because
  `pip install "nuoscprobexact[fast]"` is what every earlier instruction said
  and an extra that vanishes turns those into errors.

  **The cost is numba's numpy ceiling, which we now inherit.**  numba 0.66.0
  requires `numpy<2.5`, so installing into an environment holding numpy 2.5.1
  downgrades it to 2.4.6 — verified in a clean virtual environment, not
  inferred.  Nor is it a passing inconvenience: across 0.60 to 0.66 the
  ceiling excluded the then-current numpy either on the day numba shipped
  (0.61.0, 0.66.0) or within three months, while numpy releases a minor
  roughly every six.  The pin is numba's own and no spelling here avoids it.
  The escape hatches are `--no-deps` plus numpy by hand, or leaving numba
  installed and setting `fastkernels.USE_NUMBA = False`.

  **Deliberately unmarked.**  A `python_version` marker was drafted and
  dropped: resolving the dependency set against every target from 3.9 to 3.14
  showed all six succeed — 3.9 gets numba 0.60.0 with numpy 2.0.2, the rest
  0.66.0 — so a marker would only have withheld the compiled path from
  someone who could have had it.  The premise for the marker, a Python where
  numba cannot be installed, does not currently exist.

  Two consequences worth knowing.  Extras can add but never subtract, so the
  bare install is no longer a way to *avoid* numba; that is what the two
  escape hatches above are for.  And the `Without Numba` CI job now
  uninstalls numba explicitly rather than merely omitting the extra — without
  that its own assertion fails, reporting a break that is not one.  It is
  described there as the fallback rather than the default, which is what it
  has become.

- **The four-flavor latent roots now come from invariants formed in
  double-double arithmetic.**  Worst relative error over the nine
  Hamiltonians in the new `tests/stiff_reference.json` goes from 3.9e-16 to
  **3.6e-17** — under a fifth of a `float64` ulp — for roughly a fifth more
  time on a full `probabilities_4nu` through the compiled kernel.  The error is
  exact; the cost is a range, landing at 1.18x, 1.20x and 1.25x by three
  measurement methods with a per-pair spread of 0.76x to 1.72x, so it is quoted
  as one.

  1.12.0 escaped the invariants' conditioning by refining against the matrix
  with a Newton step; this removes the conditioning instead.  `I_2, I_3, I_4`
  compress a 4x4 matrix into three numbers, and the amplification from those
  coefficients to the roots measures 2.3e9 — so a `float64` coefficient at
  1e-16 becomes a root at 1e-7, which is what the closed form alone scores
  (2.2e-07) and what no better root-finder can beat.  Carrying each
  coefficient as a pair of `float64` limbs, some 32 digits, leaves the same
  amplification acting on 1e-32 and landing at 1e-23; the roots are then
  limited by rounding the answer to `float64` rather than by the algebra.
  One Aberth sweep, also in double-double, takes the quartic there.

  Three things keep the cost to 25%, and the first two were mistakes before
  they were optimisations.  The residual-trace shift is removed from the
  eigensolver's start **in double-double**; doing it in `float64` discards
  the low limb the exact traceless-ing just computed, puts the start an ulp
  out, and costs a second Aberth sweep to recover — a second sweep that was
  briefly mistaken for something the iteration needed, since one sweep did
  then measure 3.9e-16.  Aberth's step is evaluated as `chi*P/(chi'*P -
  chi*S)`, for `P` and `S` the product and second symmetric function of the
  three gaps, rather than as `(chi/chi')/(1 - (chi/chi')*sum 1/gap)`: the
  same expression with one double-double division per root instead of five,
  and a dd division is three `float64` divisions where a dd multiply is
  none.  And `1/3` is baked as a dd constant instead of divided for, which
  also made one reference case exact that was not.  Together these halved
  the double-double overhead — roughly 50 ms to roughly 25 ms per 100 000
  points — and improved the worst root from 5.5e-17 to 3.6e-17.

  Exploiting the exactly-zero imaginary diagonal of `H~` inside `H~^2` was
  tried too, and measured no faster — the branches cost what the skipped
  multiplications saved. What remains is close to arithmetic floor: `H~^2` in
  double-double is about 3000 flops per element, which at 100 000 elements
  over eight cores accounts for the ~12 ms it takes.

  **`oscprob4nu.ROOT_STRATEGY` selects it,** defaulting to
  `'double-double'`, with `'eigensolver'` for the LAPACK route.
  `oscprob4nu.POLISH_ROOTS` keeps its meaning and is read only on the latter
  — the default route reaches the floor without it and leaves a Newton step
  nothing to find.  Both routes exist on both backends and agree to an ulp,
  so a result never depends on whether a compiler was present.

  Three things were measured and rejected, and are recorded because each
  looks like an obvious improvement.  Mirroring `H~^2` across its diagonal
  saves a quarter of the work and quietly costs `I_3` and `I_4` 1e-17,
  because a Hamiltonian built as `U M^2 U†` is Hermitian only to rounding
  and the mirror discards the asymmetry — Hermitizing exactly in
  double-double first makes the saving sound, and that is what ships.
  Starting from Euler's closed form instead of the eigensolver is twice as
  fast and exact on every stiff case, and still fails, because on a cluster
  Aberth converges *linearly* at ratio one half: from 1e-7 five sweeps reach
  3.8e-9 where thirty are needed.  Durand-Kerner needs one division per root
  against Aberth's five and is non-monotone in the sweep count — 3.9e-16,
  9.7e-17, 1.9e-16 — so it cannot be a default.

  **`fastkernels`' four-flavor kernel entry points take `strategy: int`
  where they took `polish: bool`.**  A compiled kernel cannot read a Python
  global at call time, so the two switches travel as one integer: 0 for
  double-double, 1 for the eigensolver with the Newton step, 2 for the
  eigensolver alone.  Callers going through `oscprob4nu` are unaffected.

  Two limits worth stating plainly.  On the NumPy fallback the ratio is of
  order 1.5-2x rather than 1.2x, the double-double primitives being
  elementwise array operations with nothing to amortise them over.  And this
  buys *roots*, not probabilities everywhere: reconstructing `U_4` in Newton
  form takes second differences of `exp(-i psi L)`, so a root error enters
  the coefficients twice and the probability error grows as
  `(psi_max L)^2` — measured at 5e-17 times that across five decades of
  phase.  At Delta m^2_41 = 1000 eV^2 over 1300 km the phase is 2.5 million
  radians and both routes land at 2e-4 together.  The physical range is
  where the tight figures are: 1.3e-12 at 0.1 eV^2.

### Added

- **A runnable example on the eight public routines that had none.**  The six
  `earth` entry points --- `probabilities_{2,3,4}nu_earth` and the three
  `_between_locations` wrappers --- are the most-called functions in the
  library and were the only ones in that module without one, and
  `probabilities_{2,4}nu_profile` were without one where their three-flavor
  sibling had it.  The `earth` examples show the array-of-energies form
  alongside the scalar one, since that is the addition of this release cycle
  most likely to be missed.  Every block is executed by Sphinx on the API
  page and again by `tests/test_docstrings.py` on all five Pythons.

- **A frozen fifty-digit oracle for the four-flavor roots,**
  `tests/stiff_reference.json`, with `tests/gen_stiff_reference.py` to
  regenerate it.  Nine Hamiltonians — the 3+1 family from 0.1 to 1000 eV^2,
  two random Hermitian, and two with a pair separated by 1e-8 and 1e-16 —
  stored as hexadecimal floats so a reader gets the exact bits the reference
  was computed from.

  It has to be frozen because `numpy.linalg.eigvalsh` stopped being an
  independent oracle the moment it became one of `ROOT_STRATEGY`'s routes:
  tests written against it began comparing an implementation with itself,
  one of them asserting `1.06e-15 <= 0.0`.  The generator uses `mpmath.eig`
  for the roots and `mpmath.expm` for the probabilities, neither sharing code
  with the library nor going through the invariants, and checks the stored
  roots against `det(psi*I - H~)` taken straight from the matrix.

  That check is there because the generator got it wrong first, and its
  docstring records all three traps.  The quartic built from the invariants
  is not always the matrix's characteristic polynomial — `float64`
  traceless-ing leaves a residual trace while that quartic pins its cubic
  coefficient to exactly zero, and the two disagree by up to 3.8e-7 wherever
  the residue is nonzero, agreeing to 2e-17 wherever it is not, which is what
  makes it easy to miss.  Two solvers agreeing is not evidence: the wrong
  reference was corroborated by `mpmath.eig` and led to `eigvalsh` being
  blamed for residuals that were the oracle's.  And probability conservation
  is not an accuracy measure, since a row of `|U|^2` sums to one for any
  unitary `U` however wrong its eigenvalues.

### Fixed


- **The documented way to write a varying profile was quietly wrong.**  The
  example in `slabs.probabilities_3nu_profile` scaled its potential by
  `x[-1]`, the last position it was handed.  Those positions are the
  *midpoints of the current refinement*, so they move every time the slab
  count doubles: dividing by the last one makes the profile itself a function
  of `n_slabs`.  Measured on a Hamiltonian with any mixing at all, that drops
  the convergence from fourth order to **first** — 1.03, 1.02, 1.01 over
  successive doublings against 4.00 for the same profile scaled by the
  baseline — and `atol=1e-8` then exhausts `n_max=1024` and raises where the
  corrected form meets it at **16 slabs**.

  It was invisible in the shipped example because that example's Hamiltonian
  was diagonal at every point, so there was no mixing, `P_ee` was exactly
  1.000000, and the refinement never had to do anything.  A reader adapting it
  to a Hamiltonian that oscillates got a `ValueError` and no clue why.  The
  example now carries an off-diagonal entry, so it exercises the refinement it
  is demonstrating, and the pitfall is named in the example, in
  `recipes.rst` and in `README.md`.

- **Five documented constants reached no page.**  `MAX_CHUNK_BYTES`,
  `CHUNK_BYTES_FALLBACK`, `CHUNK_BYTES_MIN`, `CHUNK_BYTES_MAX` and
  `MIN_CHUNK_ENERGIES` each carry an autodoc docstring, and `MAX_CHUNK_BYTES`
  is named in this changelog as the knob to retune — but none was in
  `earth.__all__`, so `automodule` skipped all five.  The documentation had
  been written for a page it never appeared on: zero occurrences in the built
  API reference, against eleven for a sibling that was listed.

- **`README.md`'s requirements table contradicted its own introduction.**  The
  lead-in said the rows without an extra "need nothing beyond `numpy`",
  directly above a `numba` row marked as installed by default.


- **The copied-out-module test stripped `sys.path` by name, not by path.**
  `test_a_core_module_works_copied_out_on_its_own` removed every entry whose
  *text* contained `NuOscProbExact`, which is not the set it meant: it also
  removes unrelated directories carrying the project's name, and the
  commonest of those is a virtual environment created inside the checkout.
  Its `site-packages` holds NumPy, so the subprocess lost NumPy and the test
  reported `oscprob2nu does not work copied out on its own:
  ModuleNotFoundError: No module named 'numpy'` — the library blamed for a
  path collision, on all three modules at once.  A `.venv` in the project
  root, which is a common layout, failed three tests on a clean checkout.

  Now the repository root, `src` and `tests` come off by `os.path.realpath`
  comparison.  A `startswith` test on the root, the obvious alternative,
  reproduces the bug exactly, since such an environment lives under it.  The
  assertion that `worthwhile()` returns False still proves `src` really left
  the path, so the test cannot quietly become vacuous.

  It also now fails loudly on its own account: if NumPy is unimportable once
  the stripping is done, the subprocess exits saying `TEST FAULT` and dumping
  `sys.path`, instead of letting it look like the library.  Probed by
  reinstating the old strip, which now names its own cause.

## [1.12.0] - 2026-08-02

### Added

- **`probabilities_{2,3,4}nu_earth` take an array of energies.**  They took a
  scalar, so an energy scan was N separate Python calls, each rebuilding the
  matter potentials, re-validating slabs the library had just built itself,
  and launching its own kernel.  Passing an array now evaluates the whole
  scan in one pass and returns `(..., 4)`, `(..., 9)` or `(..., 16)`; a
  scalar energy still returns a tuple, bit-for-bit what it returned before.
  The three `_between_locations` wrappers inherit it.

  Two things made this worth more than the bookkeeping it saves.  The matter
  potential depends on the geometry alone, so a scan at fixed zenith angle
  was recomputing the identical array once per energy — measured at 11% of a
  hundred-energy scan, alongside 9% spent re-validating.  And the chords are
  independent of one another, which the per-chord kernel could not exploit
  because the product *along* a chord does not commute: the new
  `slab_product_{2,3,4}nu_batch_kernel` spreads them over the cores instead,
  above a threshold on `n_chords*n_slabs` rather than on the chord count,
  since a chord costs in proportion to its slabs.

  A hundred-energy scan is **6.7x, 4.8x and 3.5x** quicker at two, three and
  four flavors than the loop it replaces; a thousand-energy scan is 9.6x,
  4.6x and 1.9x.  Against the uncompiled per-energy loop, a thousand-energy
  scan is 121x, 49x and 17x.

  Two limits worth stating.  Four flavors gains least and falls off with
  length because the Hamiltonian stack, and the traceless copy the expansion
  needs, dominate: this is memory, not arithmetic.  And the stack is
  proportional to the scan — a hundred thousand energies across a 120-slab
  chord would be 1.6 GB at three flavors and nearly 6 GB at four, which is an
  ordinary oscillogram rather than an abusive input — so long scans are
  evaluated in chunks of `earth.MAX_CHUNK_BYTES`, bit-identically and at no
  measurable cost.

  That chunk turned out to matter for a second reason, and it is the more
  interesting one: **the batched kernel is memory-bound, not compute-bound.**
  The stack is written by the Hamiltonian builder and then streamed by the
  kernel, which does little arithmetic per byte.  Cost per probability
  measured 7.9 µs with an 8 MB working set and 16.4 µs with a 540 MB one — a
  factor of two paid for traffic alone — and thread scaling flattens at eight
  threads and regresses at twelve, which is the same story from the other
  side.  `MAX_CHUNK_BYTES` therefore defaults to the **detected last-level
  cache** rather than to a fixed 64 MB, worth 1.3x to 1.8x on long scans
  across both chord lengths and all three flavor counts, and never slower.

  Detection tries `os.sysconf`, Linux's `sysfs` and macOS's `sysctl` (through
  `ctypes`, not a subprocess), and falls back to `CHUNK_BYTES_FALLBACK` where
  none answers — Windows, notably, which is not probed because doing it means
  untested `kernel32` calls at import time.  Every probe is wrapped broadly
  and deliberately: the value only tunes how a scan is cut up, so no failure
  of it should be able to stop the module importing.  It is a plain module
  attribute, and the default is a guess at the right order of magnitude
  rather than a tuned constant — the optimum is broad, and any value near the
  cache beats one far above it.

  What is *not* done: building the per-slab Hamiltonians inside the kernel.
  That is the largest remaining item, worth about 19% of a scan, and it was
  left alone deliberately because it would couple `fastkernels` to the
  matter-Hamiltonian convention that lives in `hamiltonians*nu`.  Batching
  captured the larger share without that.

- **`tests/bit_capture.py` and `tests/bit_compare.py`**, the harness for a
  change that is meant to leave every number alone.  Neither is collected by
  pytest; they are run in pairs, before and after, and the comparison
  demands `numpy.array_equal` rather than `allclose`, reporting any
  disagreement in ulps as well as absolutely.

  No golden file is committed with them, deliberately: `acos`, `cos` and
  `sin` come from the platform's libm and are not guaranteed bit-identical
  between implementations, so a stored baseline would fail somewhere in CI
  while saying nothing about the change under test.  Two captures taken on
  the same machine in the same session remove that variable.

  The harness was validated before being trusted, by mutating a kernel
  rather than by inspection: re-associating one expression in `_entries_3nu`
  from `(x*SQRT3)/6.0` to `x*(SQRT3/6.0)` moved 133 of 870 arrays, by up to
  17000 ulps — and moved *only* the three-flavor compiled arrays, leaving
  the other 725 untouched.  It captures evolution operators as well as
  probabilities for the same reason: `|U|^2` discards phase, and on that
  mutant the operators showed 17000 ulps where the probabilities showed
  1800.

  Its first use was to **decline** a change rather than certify one, which
  is worth recording.  A set of loop-invariant hoists and temporary
  removals in the per-slab arithmetic — bit-identical in principle, and
  estimated at six to ten per cent — turned out to be unmeasurable here:
  the run-to-run spread of the *same build* on this machine is 27% to 56%,
  and an untouched two-flavor control showed an apparent 0.906x "effect".
  Every candidate sat inside that.  They are left undone rather than
  committed on a claim this hardware cannot check; the earlier wins in this
  release are all far above that floor and stand.  The harness is committed
  so that whoever revisits them on a quieter machine does not rebuild it.

- **A chord is a palindrome, and half its operators were being computed
  twice.**  A neutrino crossing a spherically symmetric Earth meets every
  radius on the way in and again on the way out, so slab `j` and slab
  `n-1-j` carry the same density, the same width, and therefore the *same*
  evolution operator.  Every Earth probability was computing each of them
  twice.

  Writing `U = U_{n-1} ... U_0` and splitting at the centre gives
  `U = (U_0 U_1 ... U_{m-1})(U_{m-1} ... U_0) = A B` for even `n`, and
  `U = A U_m B` for odd `n`, with the middle slab its own mirror.  Both
  products accumulate in one pass over the first half, so each expansion is
  computed once and used twice.  Only the expansions halve — the matrix
  products still number `n`, since two running products are carried instead
  of one — and the expansion is about two thirds of a slab's cost, which is
  why this is worth about 1.5x rather than 2x.

  Measured end to end, interleaved: **1.34x to 1.50x at two flavors, 1.55x
  to 1.72x at three, and 1.78x to 1.80x at four**, the last being largest
  because `_operator_4nu` is the heaviest thing a slab does.  A three-flavor
  scan now costs 2.65 microseconds per probability, against 4.8 before and
  about 38 when this release opened.  An oscillogram is 1.47x, at 4.3
  microseconds per point; a single Earth call is 1.19x.

  **It is not restricted to PREM.**  `fastkernels.palindromic` decides, and
  `slabs` asks it of any slab sequence it is handed, so a symmetric castle
  wall or any hand-built profile that reads the same from either end is
  composed at the same discount.  The Earth is simply the case that always
  qualifies.

  Three things this cost, all worth recording.

  The test is **exact equality, never a tolerance**.  The saving relies on
  the mirrored operators being identical, which follows from identical
  inputs and from nothing weaker; a tolerance would hand a
  nearly-symmetric profile the answer for a symmetric one, silently, which
  is the one thing an optimisation must not do.

  For that test to succeed, `earth_slabs` now emits **exactly** palindromic
  widths.  It did not before: each segment was cut by its own `linspace`
  and the two halves rounded differently, by about 1e-12 km on a 100 km
  slab.  Averaging each width with its mirror makes both bitwise equal.
  This release's results are therefore **not** bit-identical to 1.11.0 on
  any Earth call, which is stated here rather than discovered later.
  Measured over 4640 values across four angles, four subdivisions and five
  energies, the shift is **6e-15 at two flavors, 7e-14 at three, and 4e-10
  at four**.

  The four-flavor figure is much the largest and is not this change's
  doing: that path polishes quartic latent roots, and its conditioning
  already separates the two *backends* by 4.0e-10 on code this release
  never touched.  A 1e-12 km perturbation lands in the same place.  For
  scale, the discretisation error at the default subdivision is 1e-4 to
  1e-5 — five orders of magnitude larger than even the four-flavor shift,
  and eleven larger than the three-flavor one.

  And the check has to be compiled.  Comparing a materialised `(n, 3, 3)`
  complex stack against its reverse through NumPy costs about 6
  microseconds — reversed views, temporaries, and a full pass whatever the
  answer — against the 6 microseconds the halved composition saves on a
  120-slab chord.  Measured, that turned a 1.55x win into a 0.98x loss.
  `_palindromic_stack` walks the first half and returns at the first
  disagreement, and is 8x cheaper.

  `fastkernels.USE_PALINDROME` switches it, and is True by default.  It is
  not a correctness switch --- the two orderings agree to a few times
  1e-15 --- but the plain left-to-right product is what a comparison needs
  when the question is whether a discrepancy is the palindrome's doing.

  Composition ordering is regrouped, so results move by 2e-15 to 5e-15
  against composing the whole chord.  Unitarity is unchanged — measured
  identical at three flavors and at four, marginally better for the
  mirrored composer at two of three angles.

- **The Earth kernels build their own Hamiltonians.**  `earth_chords_{2,3,4}nu_kernel`
  take the energy-independent vacuum Hamiltonian, the per-slab potentials and
  the widths, and construct each slab's matter Hamiltonian in registers.  The
  batch kernels they replace on this path were handed a stack — one matrix per
  slab per energy, 17 KB per chord, written by the builder and then streamed
  back — and that stack was what made a scan memory-bound.  What the fused
  kernels read instead is one potential and one width per slab, shared by every
  energy in the scan.

  Worth being plain about what this is and is not: **it is a memory layout, not
  a scheme.**  The Hamiltonian is the same midpoint Hamiltonian
  `hamiltonians3nu.hamiltonian_3nu_matter` builds, and the fused output is
  bit-for-bit the batch kernel's at all three flavor counts, which is what the
  tests assert.  The cost is that `fastkernels` now knows the
  matter-Hamiltonian convention, which is a real coupling and was weighed
  rather than assumed.

  On top of batching: **1.9x to 2.4x at three and four flavors**, more at two.
  Per probability, a three-flavor scan went from about 8.5 to 4.8 microseconds.
  The general batched composer in `slabs` stays — it takes an arbitrary stack
  of Hamiltonians, which the fused kernels cannot, knowing only
  ``H_vac/E + V P`` — and is now tested directly rather than through `earth`.

- **An array of zenith angles, which is to say oscillograms.**  The energies
  and angles broadcast against each other, so a grid is asked for as
  ``probabilities_3nu_earth(h, energies[None, :], costhz[:, None])`` — the
  idiom `oscprob3nu.probabilities_3nu` already uses for Hamiltonians against
  baselines.  Measured **8x** against the loop over grid points it replaces,
  at 4.8 microseconds per point for a 60 x 200 grid.

  The angles are evaluated one at a time, and that is not a shortcut: the chord
  geometry is what the angle changes, so two angles share neither their slab
  widths nor even how many slabs they have, and there is nothing for a single
  kernel call to share.  What they do share is every energy at a given angle,
  which is the axis the fused kernel already spreads over — so the work drops
  from one call per grid point to one call per distinct angle.  A broadcast
  grid repeats each angle many times, and `numpy.unique` means it is cut once.

  `earth_slabs` deliberately stays scalar in the angle: it returns *a* chord's
  widths and densities, and there is no array shape that holds two chords of
  different lengths.

- **`rtol` and `atol` on the Earth probability routines**, and
  `earth.slabs_for_tolerance`, which is what they call.  The discretisation
  error is strongly energy- and angle-dependent — at the default eight
  sub-slabs per segment it spans more than an order of magnitude between 3
  and 40 GeV — so a fixed `n_slabs_per_segment` does not give a fixed
  accuracy.  Either tolerance being set refines the chord until the measured
  error meets it; both unset leaves the result bit-identical to before, which
  is the regression that matters most here.

  The threshold is `atol + rtol*abs(P)`, the convention of `numpy.isclose`,
  and it binds on **every** returned probability, so the subdivision is set
  by the worst-converging entry and no channel in the tuple is quietly
  less converged than was asked for.  For an array of energies one
  subdivision covers the whole scan, chosen from the worst energy in it.

  The search doubles the subdivision until the error estimate passes, rather
  than solving the second-order law for the required `n` in one jump.  It
  costs about twice a single evaluation at the subdivision it returns —
  measured at 2.2x, the geometric series being dominated by its last term —
  and in exchange every value it returns has had its error *measured* rather
  than extrapolated, which is what makes it able to say when a tolerance
  cannot be met.  It raises then, rather than clamping and warning: a
  warning Python prints once and thereafter suppresses is a poor way to tell
  someone their stated accuracy was not honoured.

  Because the error of the coarser evaluation is four thirds of the gap, a
  tolerance the default already meets is met by `n_slabs_per_segment` itself
  rather than by twice it, so a loose tolerance costs nothing afterwards.
  `return_n_slabs=True` reports what was chosen; a tight tolerance can
  quietly cost a great deal of refinement and should be discoverable.

  Verified against second-order Richardson references at three zenith angles
  and three energies: the subdivision returned meets the tolerance in every
  case, and it ranges from 8 to 256 across them.  The helper is meant to be
  called **once** per scan and its answer reused — passing `rtol` to every
  point of a loop repeats the search at every point.

- **`slabs.probabilities_{2,3,4}nu_profile`**, the general counterpart of the
  Earth routines: they take the profile as a callable rather than knowing
  about PREM, so the same tolerance machinery serves a hand-built density
  profile, a castle wall, or a solar model.  `rtol` and `atol` do not appear
  in `oscprob{2,3,4}nu`, and deliberately: those compute an exact closed form
  for a single Hamiltonian, with no discretisation to refine, so a tolerance
  there would be a parameter that could not be honoured.

  Equal slabs sampled at their midpoints, which is second order and so
  refines by the law the search assumes — verified at a ratio of 4.00.
  Where a profile is *discontinuous* this is the wrong tool, since no
  refinement recovers a jump that straddles a slab; the docstring says so and
  says to split the trajectory at the discontinuities instead, which is what
  `earth` does with the PREM shells.

- **A compiled path for `slabs` and `earth`, which had none.**  The backend
  offered probability kernels only, so composing evolution operators across
  slabs — which is what layered matter and the Earth do, and they never call
  `probabilities_3nu` — ran the NumPy path however the `fast` extra was
  installed.  Measured before the fix: `probabilities_3nu_earth` took 2.11 ms
  with Numba and 2.09 ms without, entered a kernel zero times, and never even
  consulted `worthwhile()`.

  Six new kernels: `evolution_operator_{2,3,4}nu_kernel` and
  `slab_product_{2,3,4}nu_kernel`, with `MIN_SLAB_BATCH` and
  `worthwhile_slabs` governing the second set.  An Earth crossing is now
  **13.9x, 9.6x and 6.6x** quicker at two, three and four flavors, unchanged
  to 1e-14 or better.

  Three things are worth recording about how that number was reached.

  The operators were never the bottleneck.  Compiling them alone bought
  1.33x.  At 10 GeV and cos θz = −0.8, 120 slabs, the 396 µs divided as
  176 µs of PREM chord geometry, 126 µs composing the product in a Python
  loop, 55 µs of evolution operators and 10 µs building the Hamiltonians.

  `earth_slabs` depends on the zenith angle alone — not the energy, not the
  flavor count, not the Hamiltonian — yet an energy scan recomputed the
  identical chord for every point.  It is now cached; a hundred-energy scan
  reports 899 hits and one miss.  The public function still returns arrays
  the caller owns, because it copies; the cached originals are read-only so
  a stray write raises rather than poisoning the cache.

  And the slab product needed its own threshold.  `MIN_BATCH` weighs
  compiled arithmetic against NumPy arithmetic, which is why two flavors
  sits at fifty thousand; the slab product instead replaces a *Python loop
  of dispatched matrix products*, so it is ahead at every length.  Measured
  over 1 to 256 slabs it leads by 137x, 225x and 187x at one slab and 59x,
  13x and 7x at 256 — the margin narrowing with length, the reverse of the
  probability kernels, because the fixed cost being avoided is the
  caller's.  Reusing `MIN_BATCH` would have left two flavors on the NumPy
  path for every Earth crossing.

- **`fastkernels.evolution_operator_3nu_kernel`**, and the dispatch in
  `oscprob3nu` that reaches it.  The compiled backend had no
  evolution-operator kernel at all, only probability ones — so `slabs` and
  `earth`, which compose *operators* across adjacent slabs and never call
  `probabilities_3nu`, ran the NumPy path however the backend was
  configured.  Installing the `fast` extra bought an Earth crossing exactly
  nothing: measured before the fix, `probabilities_3nu_earth` took 2.11 ms
  with Numba and 2.09 ms without, and the kernel was entered zero times —
  `worthwhile()` was never even consulted.

  The arithmetic was already there.  `_one_3nu` built all nine entries of
  `U_3` and squared them on the next line, so the operator kernel is that
  computation stopping one step earlier; the shared part is factored into
  `_entries_3nu` and inlined into both, leaving the probability path
  unchanged.  Agreement with the NumPy path is 8.9e-15 over 200 random
  Hamiltonians plus the degenerate branches, which match exactly.

  **This is three flavors only.**  Four is the same shape of change —
  `_one_4nu` already materialises the operator in scratch — but two is not:
  `_one_2nu` never forms `U` at all, going straight to the closed-form
  `P_eμ`, so a two-flavor operator kernel is new code rather than a
  refactor.

  Worth stating plainly: **the win is smaller than the gap suggested.**  An
  Earth call is 1.33x faster, not an order of magnitude, because computing
  the operators was never the bottleneck.  At 10 GeV and cos θz = −0.8, 120
  slabs, the 396 µs breaks down as 176 µs of PREM chord geometry (44%),
  126 µs composing the product in a Python loop (32%), 55 µs of evolution
  operators (14%, now compiled) and 10 µs building the Hamiltonians.  The
  remaining levers are the composition loop, which could be folded into the
  kernel, and `earth_slabs`, which is purely geometric and recomputed for
  every energy in a scan even though it depends only on the zenith angle.

- **`tests/prem_scan.py` and `tests/prem_speed_accuracy.json`**, the frozen
  data behind the paper's two Earth speed–accuracy figures, at three flavors
  and at 3+1.  The module docstring carries the full provenance: the referee,
  the check on the referee, and every convention that had to be matched
  before six codes were propagating through the same Earth as this one.

- **`tests/external_drivers/`**, the C and C++ drivers for `NuFast-Earth`,
  `GLoBES` and `Prob3++`, with the raw output of each and a README recording
  what every one of them had to be told.  The previous round of this
  comparison built its drivers in a session scratchpad and lost them.
  `gen_prem_header.py` emits this library's PREM as a C header from
  `earth._PREM_COEFFS`, rather than anyone transcribing forty coefficients;
  it reproduces `earth.density_prem` to 5e-14 over 6372 radii.

- **Notebook 20, `20_arbitrary_hamiltonian.ipynb`**, the case the other
  nineteen leave out.  `03` puts a non-standard Hamiltonian through matter of
  constant density; `08` and `14` put the standard one through profiles that
  vary; `01` evaluates an arbitrary Hermitian matrix once.  None of them
  carries a user's own Hamiltonian through a profile that changes under it,
  which is what `probabilities_{2,3,4}nu_profile` was added for.

  The worked case is a long-range force from a gauged `L_e - L_mu` symmetry
  (arXiv:1808.02042), whose potential is a Yukawa integral over the electrons
  of the whole body rather than a reading of the local density.  It is carried
  through three profiles: an invented body, an exponential solar-like one, and
  the Earth.

  It also documents the one thing about `earth` that is not obvious from the
  API.  `probabilities_3nu_earth` builds `H = H_vac/E + V_CC*P_ee` itself, so a
  Hamiltonian not of that form cannot go through it at all; the way through is
  `earth_slabs` to `matter_potential` to your own construction to
  `probabilities_3nu_slabs`, and with the extra term switched off that path
  reproduces `probabilities_3nu_earth` to **zero difference** over sixteen
  (angle, energy) pairs, because it is the same arithmetic.  The same recipe is
  now on the numerical recipes page.

  Three checks run inside the notebook rather than being asserted by it: the
  potential against the closed form for a uniform ball, which holds at any
  mediator mass, to 1e-9 at `mR = 1`; second-order convergence of the
  midpoint sampling, measured at 2.01; and the energy-averaged solar `P_ee`
  against the analytic adiabatic MSW result with the `cos^4(theta13)`
  correction, agreeing to the sampling floor of the energy band.  Two traps
  are recorded in the code where they were hit: the exterior part of the
  Yukawa integral has to be accumulated inward from the surface, since taking
  it as (total - interior) loses the whole tail to cancellation once
  `exp(-mR)` is small; and every PREM shell boundary needs a quadrature node
  carrying the density on each side, or the potential converges at first order
  instead of second.

### Fixed

- **The documented notebook count had rotted, and nothing was guarding it.**
  `README.md` said "Eighteen worked notebooks" in two places and the
  generator's own docstring said "fifteen", while nineteen were shipping; the
  README's own table of notebooks stopped at 18, so `19_animations` was absent
  from the one list meant to be complete.  All corrected to twenty, the table
  completed, and `test_documented_figures.py` now derives the count from
  `notebooks/` and fails on any stale spelling of it.

  The guard was probed at four sites before being trusted, on the evidence of
  the last round that two of four new guards did not work when first written.
  It caught one stale phrase that had been missed by hand, and one of the
  four perturbations did *not* trip it — "all twenty share a setup cell" has
  no noun for the pattern to anchor to — so that sentence now names the noun
  rather than the pattern being loosened into false positives.

- **`RELEASE_NOTES.md` still described 1.11.0.**  The file the tree calls
  "what changed since the published version" was the v1.11.0 announcement in
  full — every section framed "since 1.0.0" — while 1.11.0 *is* the published
  version, so it described the past rather than the delta.  Rewritten for
  1.12.0, leading with the one thing an upgrader needs: Earth results are no
  longer bit-identical, by 6e-15, 7e-14 and 4e-10, against a discretisation
  error five orders of magnitude larger.

  It had also gone stale in the way everything here does when unguarded,
  claiming eighteen notebooks and 596 tests.  Because it now quotes the Earth
  crossing and palindrome figures, it joins the documents those two guards
  check, so it can no longer drift away from them; the test count is stated as
  a bound rather than a figure, since that is the one number in it certain to
  grow.

- **The Earth performance figures had drifted, and nothing was guarding them.**
  `13.9x, 9.6x and 6.6x` was stated in `README.md`, `index.rst` and
  `fastkernels`, agreed with itself in all three, and was wrong in two by the
  time this release shipped — the palindrome had made the compiled side
  quicker, so three and four flavors were being understated.  Re-measured, one
  Earth crossing is ~12x, ~12x and ~9x.

  `test_documented_figures.py` covered the constant-density tables and not
  these, which is exactly the gap its own docstring warns about: "a guard that
  covers part of a table is an invitation to the rest of it."  It now pins the
  Earth crossing speedups, the 2000-energy scan table, what the palindrome is
  worth, and that every document claiming the palindrome also names the switch
  that turns it off.  Each new guard was checked by perturbing a document and
  confirming it failed — two of the four did not, at first: a bare `1.5`
  matched the unrelated `1.5x to 20x` span, and `USE_PALINDROME` matched any
  misspelling containing it.

  The narrative was stale as well as the numbers.  `index.rst`, `quickstart.rst`
  and `README.md` described the slab kernel and none of what followed it —
  array energies and angles, the fused chord kernels, the palindrome, `rtol`
  and `atol`.  All three now carry the 2000-energy scan table, the oscillogram
  figure, and the tolerance machinery, with `recipes.rst` holding the worked
  examples.

- **`test_root_polishing_is_what_makes_a_stiff_spectrum_accurate` measured
  whichever backend happened to be dispatched.**  When the compiled path
  took over `probabilities_4nu` earlier in this release, the gain it
  measured at 1 GeV fell from 418x to 10.6x — close enough to its threshold
  of ten that it passed on one machine and failed on another.  It now runs
  per backend.

  Investigating that turned up no bug to fix, which is the useful part.
  The two paths round the closed-form roots differently, and in a cluster
  where the closest pair is 2.6e-4 of the spectrum apart, that decides
  which Newton steps the halfway guard in `_polish_roots` admits.  Neither
  path is better: swept over energy, the compiled one leads at 0.5 and
  10 GeV and trails at 1 and 2.

  The larger point is that the gain from refining is a function of how
  degenerate the spectrum is, not a constant to assert.  At 0.5 GeV, with
  the closest pair 1.4e-4 apart, refining is worth some three thousand
  times on both paths; by 25 GeV, at 3.3e-3, it is worth nothing on either
  and is very slightly harmful — the guard's guarantee being about the
  roots, not about phases accumulated over a long baseline.  The test now
  measures at the degenerate end, where the margin is four orders of
  magnitude wide, rather than where it happened to be written.

- **The convergence order of the slab product was being quoted as first.**
  It is second: the density is sampled at slab *midpoints* and the PREM shell
  boundaries are cut exactly, so the leading error is O(h^2).  Measured
  exponents are 2.000 at 3, 10 and 40 GeV.  The consequence is practical — a
  Richardson reference built for first order, `2*P(32) - P(16)`, is worse
  than plain `P(32)`, while `[4*P(256) - P(128)]/3` reaches 1e-11 and is
  confirmed there by an independent integration of the continuous profile.

- **Three external codes carry three different roundings of the same
  constant.**  `NuFast-Earth`'s 1.526493231029146e-4, `GLoBES`'
  7.63247e-14 and `Prob3++`'s 1.52588e-4 are all sqrt(2) G_F over the atomic
  mass unit, and no two agree, so each needs its own factor against this
  library's 1.514423e-4 — 0.9920938, 0.992093 and 0.9924922.  Scanning for
  the residual minimum returns 1.000000 for the first two.

- **The star-product deviation in `methodology.rst` quoted a sample
  extremum.**  "Between 30% and 230%" is the largest and smallest of one
  finite draw and grows with the sample; at two thousand Hamiltonians the
  same measurement gives 27% to 263%.  The median, 56%, is stable and is what
  the text now leans on, with the central 90% of draws quoted beside it.

## [1.11.0] - 2026-08-02

**The pre-publishing audit of `src/`.**  Two passes over all ten modules,
7378 lines, checking the code, the comments and the docstrings against each
other and against what the code actually does when run.  A minor release:
it adds public API — four-flavor propagation through layered matter and the
Earth, a Lorentz-violating four-flavor term, and a Hermiticity check — and
fixes six docstrings that stated the opposite of the code.

### Added

- **Four flavors through layered matter and the Earth.**
  `slabs.evolution_operator_4nu_slabs`, `slabs.probabilities_4nu_slabs`,
  `earth.probabilities_4nu_earth` and
  `earth.probabilities_4nu_between_locations`.  Everything they need already
  existed — `hamiltonians4nu.hamiltonian_4nu_matter` accepts an array of
  potentials, which is the shape `earth` feeds it — so a 3+1 user could
  compute a probability but could not propagate one through the Earth, which
  is the case the sterile matter resonance lives in.

- **`earth.matter_potential_nc`**, the neutral-current potential.  It is
  flavor-universal across the three active states, so at two and three
  flavors it is proportional to the identity and drops out, which is why
  `earth` never needed it.  A sterile state does not feel it, so it stops
  cancelling and leaves `-V_NC` on the sterile entry — the entry that places
  the resonance.  It reproduces `globaldefs.VNC_EARTH_CRUST` exactly.

- **`hamiltonians4nu.hamiltonian_4nu_liv`**, the missing counterpart to the
  three-flavor LIV term.  `b4` is an eigenvalue like the others, so a sterile
  state may couple to the LIV background whether or not it couples to matter.

- **`CHECK_HERMITICITY`**, in each of `oscprob2nu`, `oscprob3nu` and
  `oscprob4nu`.  A non-Hermitian Hamiltonian was accepted silently, and
  nothing downstream revealed it: the probabilities returned still **sum to
  one**, so the check a caller would actually apply could not tell the answer
  was meaningless.  All seven public entry points now verify it.

  It is not free, and the cost is stated rather than buried: validating a
  stack is a pass over it, the same order as evaluating it, and the compiled
  kernel has made evaluating it fast.  Interleaved, best of fifteen, the
  check costs 1.3x–1.8x on a 2000-point scan and 3.2x–5.7x on a 200 000-point
  one.  It defaults to `True` anyway, because a library that silently returns
  meaningless numbers costs its user more than the check does; set it to
  `False` for scans whose Hamiltonians are Hermitian by construction.

- **44 tests**, covering the four-flavor slab and Earth paths, the
  Hermiticity refusal on every entry point and both backends, the mixing-angle
  and `costhz` domain checks, the `earth` edge cases below, and a guard that
  the two helpers duplicated across three modules each do not drift apart —
  `oscprob2nu` and `oscprob3nu` are documented as self-contained, so a shared
  module would break the property that makes copying them work, and
  duplication is the deliberate cost.

### Fixed

- **Every image and link in `README.md` would have been broken on PyPI.**
  The README *is* the PyPI long description, and PyPI resolves relative
  URLs against `pypi.org` — so the ten-image gallery the file opens with,
  and thirty-two links to notebooks, examples, `LICENSE` and the changelog,
  would every one have 404ed on the page that `pip install` sends people
  to.  They are now absolute, images through `raw.githubusercontent.com`
  and the rest through `blob`/`tree`, which changes nothing about how
  GitHub renders them.  The twenty-seven in-page anchors are left alone.

  Worth knowing why nothing caught this: `twine check` passes, and passed
  before the fix.  It validates that the markup *parses* as PyPI will
  parse it, not that anything it points at resolves.  This was found by
  building the distribution and reading the `METADATA` out of the wheel.

- **A refusal branch lost its only test.** Giving `_check_hermitian` a
  path for a single matrix took the scalar and batched routes off a shared
  implementation, and the test for an imaginary entry on the diagonal
  passes a nested list — so it stopped covering the batched version of
  that refusal, which was then exercised by nothing.  Coverage fell from
  100% to 99.62%, comfortably above the floor and therefore green.  The
  test now checks both paths, as its non-finite counterpart already does
  for the same reason.

- **`CHECK_HERMITICITY` made a single probability nine times dearer, and
  nothing measured it.**  The check was written for stacks, as a dozen
  whole-array reductions, and that is the right shape at two hundred
  thousand elements and pure overhead at one: 59 µs to validate a single
  3x3, against 0.63 µs per element on a stack of two thousand.  A scalar
  three-flavor probability went from 8 µs to 143 µs — and
  `probabilities_3nu` paid it **twice**, being the only one of the three
  that reached its answer through the public `evolution_operator_3nu`,
  which validates the same matrix again.

  Two consequences, neither visible from the test suite, which checks that
  the documented figures agree with *each other* and never against the
  code.  The 8 µs quoted in five places became wrong by a factor of nine.
  And `SMALL_BATCH`, which exists to send short stacks down the cheaper
  scalar path, inverted: at its own threshold of ten the scalar path was
  **17.7x slower** than batching, so the optimisation had become a
  pessimisation at every size.

  `_check_hermitian` now takes a separate path for a single matrix,
  comparing the entries as Python complex numbers, and `probabilities_3nu`
  assembles the operator through a private helper rather than through the
  public routine.  Measured after: a scalar three-flavor probability costs
  17.5 µs against 143, a two-flavor one 4.3 against 34, and the check
  itself 3.6 µs against 59.  The check now costs 1.35x on the scalar path
  rather than 9.4x, in line with the 1.3x–1.8x it costs on a stack.

  The thresholds were then re-measured, interleaved, best of seven, since
  the numbers they were set from predate the check: the crossover is
  thirteen elements at three flavors and twelve at two, so `SMALL_BATCH`
  goes from 10 to 12 and from 6 to 11.  Note these govern the NumPy path
  only — with the compiled backend installed, `fastkernels.MIN_BATCH` is
  one at three and four flavors, so the kernel takes every stack before
  `SMALL_BATCH` is consulted.

  A later pass caught what that re-measurement broke around it.
  `oscprob4nu.SMALL_BATCH` was documented as matching the three-flavor
  threshold, which stopped being true the moment the three-flavor one
  moved; its docstring also described a scalar path that four flavors does
  not have, the quickstart having said so for two releases.  Nothing reads
  that constant, and it now says as much.  `methodology.rst` explained the
  two thresholds as differing because two flavors amortises its overhead
  sooner — an explanation for a gap of four elements that survived the gap
  becoming one, and one that contradicts the `oscprob2nu` docstring beside
  it: they are close *despite* the difference in work per element, because
  what is amortised is the array machinery's fixed cost.  And
  `_check_hermitian`'s own docstring described only the stack strategy,
  never the single-matrix branch that is now half of it.

- **A nan could pass the scalar Hermiticity check.**  The first version of
  the fast path above built its scale with `max`, which keeps its running
  value when the comparison is against a nan — every comparison with one
  being false — so a nan never reached `isfinite`, and a Hamiltonian the
  batched path refuses came back with probabilities instead.  The two paths
  disagreeing about what is valid is worse than the cost the fast path was
  written to remove.

  The whole suite passed.  The existing guard makes its matrix non-finite
  *and* non-Hermitian, so it is refused on the second ground and cannot
  tell whether the first works; and it never exercised the batched path.
  Both gaps are now covered by a test that is Hermitian apart from being
  non-finite and runs each path separately, confirmed to fail without the
  guard.  Non-finite entries are now caught per entry, before `max`.

- **Six docstring claims that the code contradicts.**  `pmns_mixing_matrix`
  and `mixing_matrix_2nu` each said 1.1.0 made them return a
  `numpy.ndarray` "rather than a nested list"; both return a nested list, and
  always have.  `hamiltonian_2nu_liv` and `hamiltonian_3nu_liv` each said
  1.3.0 let "the matter potential be an array too"; neither has a matter
  potential.  `psi_roots` said 1.5.0 formed the arc-cosine argument "as a
  division rather than a power of -1.5", where the line reads
  `pow(h2, -1.5)`.  `distance_traveled_inside_earth` called `costhz >= 0`
  "upward-going" two lines above calling `costhz = -1` "straight up".

- **`sin(theta)` outside [-1, 1] behaved three different ways**: two and
  three flavors raised `math domain error`, naming neither the parameter nor
  the value; four flavors returned `nan` and let it run silently into the
  probabilities.  One helper now serves all three.

- **`CONV_EV_TO_G`** was `1.783e-33` against CODATA's `1.78266192e-33` — off
  by `1.9e-4` relative, three orders of magnitude worse than every other
  constant in `globaldefs`.  It reaches `VCC_EARTH_CRUST` and every Earth
  crossing.  The nuSQuIDS cross-check is unmoved by it, at 3.9e-16 and
  3.0e-10 as before.

- **`earth.dms_to_decimal` could not express a coordinate between 0 and −1
  degrees.**  The sign rides on the degree part, and `-0` is `0`, so 0°30′S
  came back North.  No predefined location reaches that band; a user-supplied
  one can, so the minutes may now carry the sign.

- **`earth.earth_radial_distance_from_depth`** took `sqrt(abs(r2))`, which
  absorbs endpoint round-off and would equally hide a genuinely negative
  radicand.  It clips at zero.

- **`earth.coordinates_of_named_location`** raised its helpful message from
  inside `except KeyError` without `from None`, burying it under a chained
  traceback.

- **`probabilities_4nu` and friends carried no `versionchanged` for 1.10.1**,
  although that release changed their results for a nearly degenerate
  spectrum.

- **A non-finite entry disabled the Hermiticity check.**  Found by a later
  pass over the check itself: the tolerance is relative to the largest
  entry, so a single infinity made the scale infinite, the tolerance
  infinite, and every comparison false — a Hamiltonian that was both
  non-finite *and* non-Hermitian passed a check whose purpose is to refuse
  the second.  Caught for the cost of one `isfinite`, since the scale is
  computed anyway.

- **`earth`'s `costhz` was not required to be a cosine.**  The chord length
  is `-2 R costhz`, which accepts anything: `costhz = -1.5` gave a chord of
  19 113 km against an Earth diameter of 12 742 km, and that chord then
  acquired seventy-six plausible slabs spanning the whole PREM density
  range.  Nothing downstream noticed.  The geometry routines now refuse a
  value outside [-1, 1], inclusive at both ends.

### Changed

- **The NSI defaults are described accurately.**  Their comment claimed
  compatibility "at 2 sigma with LMA+coherent"; with `EPS_EE = 0.06` and
  `EPS_MM = 1.2` the combination matter oscillations are sensitive to is
  `eps_ee - eps_mm = -1.14`, which is the LMA-D solution.  The values are
  **unchanged** — they are deliberately large so the worked examples show a
  visible effect — and the comment now says so and warns against reusing them
  as a fit.

- **The version is written in one place.**  `docs/source/conf.py` derives
  `release` from `pyproject.toml` and `version` from that, rather than
  restating both by hand — the short form was never a fact in the first
  place, only a derivation someone had to remember to redo.  What cannot be
  derived is guarded instead, in `tests/test_version_consistency.py`: the
  changelog must lead with the declared version, its headings must descend
  without duplicates down to 1.0.0, every `versionadded` directive must name
  a released version and none may name an unreleased one, and `src/` must
  not restate the version at all.  Six of those seven checks existed only as
  scripts run by hand during this audit.

- **14 parameters annotated `Optional[T]` with non-None defaults** now say
  `T`; three private functions gained the `Parameters` sections every other
  private function in `src/` has; `globaldefs` no longer describes itself as
  serving "plotting modules", which no longer exist, nor omits `oscprob4nu`
  from the modules that do not need it.

- **Notebook 16** said, in a code comment, that `earth` does the full PREM
  profile — which it could not do at four flavors.  It now measures the
  averaged and PREM answers against each other.

- **`README.md`'s list of installed module names** omitted `oscprob4nu` and
  `hamiltonians4nu`, which have been installed since 1.9.0.  The `slabs` and
  `earth` bullets now say which flavor counts they cover, and both the
  README and the quickstart document `CHECK_HERMITICITY` and what it costs —
  a user-facing switch worth 3.2x to 5.7x on a large scan was described
  nowhere but its own docstring, where `USE_NUMBA` has been in the quickstart
  since 1.6.0.

- **A core module copied out on its own works again.**  `README.md` and
  `installation.rst` both call copying `src/oscprob3nu.py` into your own
  project "a supported way to use **NuOscProbExact**".  It stopped being one
  in 1.6.0, when the optional compiled backend arrived and was imported
  unconditionally: a lone copy raised `ImportError: No module named
  'fastkernels'`.  Six releases, unnoticed, because any check run from inside
  the repository finds `src/` on the path and imports the real module.  The
  import is now guarded, an absent backend answered the same way switching it
  off is, and a test copies each of the three modules out into a subprocess
  with the repository stripped from `sys.path`.

- **The accuracy table on the landing page** claimed "200 random Hermitian
  Hamiltonians" where the fixture provides 100, "2000" where the test does
  400, `4e-19` for an agreement that cannot be smaller than one ulp of a
  probability (measured: 7e-16 and 1e-14), and `7e-12` where the measurement
  is 3e-14 — that last appears to be the test's *assertion threshold*
  recorded as if it were the result.

- **Two documents gave different speedups for the same four scans.**
  `methodology.rst` said ~30x, ~25x, ~40x, ~70x; the timings tabulated in
  `index.rst` and `README.md` give 21x, 23x, 37x, 99x.  The ratios are now
  derived from those timings and guarded together.

- **The SU(3) star-product identity described as "37% off" at n=4** — one
  draw quoted as characteristic.  Over two hundred random Hamiltonians the
  deviation has a median of 56% and a range of 30% to 230%.

- **Four flavors reached the documentation.**  `quickstart.rst` had a "Two
  flavors" section mirroring three and nothing for four; `recipes.rst`, which
  the landing page links as "What it can compute, with code", mentioned it
  once in passing.  Both now carry executed examples, including the sterile
  matter entry and a PREM crossing.  `CHECK_HERMITICITY` reached `index.rst`
  and `methodology.rst`, where the measured cost of the check and the two
  attempts to reduce it are now recorded.

- **A diagram of how slabs compose**, in the quickstart section above,
  drawn as SVG rather than generated: a neutrino entering as `nu_alpha` and
  leaving as `nu_beta`, four slabs of differing width and density, the
  Hamiltonian and exactly-solved `U_k` each one contributes,
  and two dashed ties that cross to show the ordering — the slab crossed
  first is the rightmost factor in the product.  That reversal is the part
  of the API most easily got wrong, and it is the one thing prose states
  and a picture shows.  Hand-written SVG keeps the labels as real text,
  adds no plotting code to a page whose point is the shortest path to a
  probability, and costs the build nothing.  It carries a `title` and
  `desc` for screen readers and an explicit white background, so that
  opening the file on its own in a dark-mode viewer does not leave the
  slate text invisible.

- **`slabs` reached the quickstart too.**  The page named it once, in a
  clause inside the four-flavor section, so a reader could finish it without
  learning that the library handles piecewise-constant matter at all — the
  answer to "what about the Earth?", which is the first question the
  constant-density example provokes.  A section now follows the evolution
  operator, where the group property `U(L_1+L_2) = U(L_2) U(L_1)` was already
  stated and then left unused, and shows a three-layer profile and a uniform
  one cut into four slabs, which returns the same probability to 6e-16.
  `recipes.rst` was deliberately left alone: its "An arbitrary matter
  profile" already *is* that recipe, down to the castle wall that changes the
  answer at 0.44 GeV while the mean density does not, and a second entry
  would only duplicate it.

- **The documentation's inert code is now run by a test.**  Snippets shown as
  `jupyter-execute` are executed by the documentation build and cannot rot
  silently; snippets shown as `code-block` are rendered and never run.  The
  landing page's *Getting started* example — the first code a reader sees —
  was one of the latter.  So were the two switch snippets in
  `quickstart.rst`, where a renamed attribute would render perfectly and do
  nothing, which is the worst way for a documented escape hatch to fail.

- **The documentation now says who wrote it.**  `grep Bustamante docs/source`
  returned exactly one hit, the `author` field inside the BibTeX entry on the
  references page — so the docs site, the one artifact a reader reaches
  without seeing `README.md` or `pyproject.toml`, credited its author only
  incidentally and gave no way to reach him.  `index.rst` gains an **Author**
  section next to Citing and License, mirroring the line `README.md` has
  carried all along, with the address already declared in `pyproject.toml`
  and in every module's `__email__`, and a pointer to GitHub issues as the
  route that leaves a public record.  It is deliberately *not* folded into
  Citing: the name is already in the BibTeX entry directly below, and no
  citation format carries an e-mail address.

- **The PyPI project page gained the links it was missing.**  `Homepage`
  and `Paper` were the only two, so the deployed documentation and the
  release history — both of which exist — were reachable from PyPI only
  through the repository.  `Documentation`, `Changelog` and `Issues` are
  now declared, along with the `Development Status` and `Operating System`
  classifiers that the page filters on.  Verified by reading them back out
  of a built wheel rather than by trusting the source.

- **A two-pass audit of every documentation file and the README**, run
  before publishing, checking each claim against the code rather than against
  the neighbouring prose.  Twenty-three findings, of which these are the ones
  a reader would have acted on:

  - `README.md` warned, in a call-out box, that a non-Hermitian matrix
    "will output nonsensical results".  It has raised `ValueError` since
    `CHECK_HERMITICITY` landed earlier in this release — the box described
    exactly the behaviour this release removed, and contradicted the
    Performance section sixty lines below it.
  - `installation.rst` omitted `oscprob4nu` and `hamiltonians4nu` from the
    modules `pip install` puts on the path, and both Requirements tables
    omitted four flavors entirely.  This is the same omission recorded as
    fixed above: it was fixed in `README.md` and left in the sibling
    document that tells the same reader the same thing.
  - The link to **Magnus** 404s, from both `index.rst` and `README.md`,
    because that repository is private.  The recommendation stays; the
    hyperlink goes until it resolves.  Found by `sphinx linkcheck`, which
    is worth a place in the release checklist: 29 external URLs, one real
    failure.
  - `README.md` annotated its own benchmark row "~93x" where the two
    timings on that line give 99, which is also what `methodology.rst`
    quotes for the same scan.
  - `fastkernels` and `methodology.rst` published different tables under
    the same heading: ~9x against ~13x for one row, ~1.5x against ~1.4x for
    another, and a different fifth scan in each.  The module is where they
    are measured, so the page now quotes it.
  - `recipes.rst` said the castle-wall profile gives "nearly three times"
    the appearance probability of a uniform one; the block above it prints
    0.0104 against 0.0028, which is 3.7.

  Two guards are widened as a result, because in each case the drift
  happened underneath a test that covered part of the same table: the
  kernel-speedup check ran on the two four-flavor rows and now runs on all
  six, and the README's inline ratios, which nothing checked at all, are now
  derived from the timings beside them.  Both new checks were confirmed to
  fail on the values they replace.

- **The API reference documented the inert copy of `SMALL_BATCH`.**
  `automodule` honours `__all__`, and only `oscprob4nu` exported it — the
  one flavor count where nothing reads it.  The two that govern dispatch on
  every call, in `oscprob2nu` and `oscprob3nu`, appeared nowhere and were
  cross-referenced twice from `methodology.rst` as targets that did not
  exist.  Exporting them also takes Sphinx in nitpicky mode from 23
  unresolved references to **zero**; the remainder were `numpy.sqrt` and
  `numpy.arccos` given as `:func:`, which NumPy's inventory does not
  register as functions because they are ufuncs.  The ordinary build under
  `-W` is silent about every one of these.

- **The README gained a License section** — the landing page had one and it
  did not — a pointer to the Zenodo DOI as the way to cite the software
  rather than the paper, which covers only two and three flavors, and a
  table of contents that lists its sections in the order they appear.  Three
  `refs.bib` entries that were rendered but never cited — PREM, NuFit 4.0
  and the LMA-D analysis — are now cited where the prose already relies on
  them.

- **Smaller corrections**: "the two core modules" in three documents where
  there are three; "every routine above accepts a stack" where the
  coefficient routines raise `TypeError`; "returns probabilities" of routines
  that return an operator; the short-stack shortcut described as universal
  when four flavors has none.

### Known limits

- **The electron number density uses the free-nucleon mean mass**,
  `(m_p + m_n)/2`, rather than the atomic mass unit.  Nucleons in nuclei are
  bound, so this undercounts `n_e`, and hence `V_CC`, by about 0.85%.  It is
  applied consistently in `globaldefs` and `earth.matter_potential`, so it is
  a modelling choice rather than an inconsistency, and it is left alone
  deliberately.

## [1.10.1] - 2026-08-02

**A Newton step that could throw a latent root across the spectrum.**  A
correctness fix in `oscprob4nu._polish_roots`, present since 1.9.0 and
affecting the NumPy path and the compiled kernel alike.  Found by mutation
testing the four-flavor kernel, not caused by it.

### Fixed

- **The root refinement is refused where it would cross a neighbour.**  The
  Newton step divides by `chi'(psi_m) = prod_l!=m (psi_m - psi_l)`, a product
  of gaps, and was guarded only by `derivative != 0.0`.  That is the right
  guard with its threshold on a knife edge: a pair separated by one unit in
  the last place gives a derivative of order 1e-16 and a step of order one.
  Observed, a root at `0.8793` was refined to `0.0180`, and the sixteen
  numbers that followed were not probabilities — reaching 21.7 and summing
  to 69.

  Whether a nearly degenerate pair lands on identical bits or on adjacent
  ones is decided by the last bit of a square root taken near zero, so the
  old test gave different answers for a stack and a scalar call, and for the
  NumPy path and the kernel — which is how the two disagreed by 186 on
  quantities that cannot exceed one.

  The guard is the standard one for polishing polynomial roots: a step for a
  simple root may not carry it more than halfway to its nearest neighbour,
  and a step that wants to is evidence the root belongs to a cluster, where
  refining against a nearly singular `chi'` only destroys what the closed
  form had.

  On synthetic spectra with pair separations drawn between 1e-16 and 1e-6
  relative, roots wrong by more than 1e-6 relative went from **10.2% to
  none**, and the worst error from 4.8 to 2.2e-6.  Kernel-against-NumPy on
  the same population went from 186 to 5.6e-11.

- **What the fix does not do**, stated because the difference matters: it
  does not make a nearly degenerate pair accurate.  Euler's reduction
  recovers the pair's separation as `sqrt(z)` for a resolvent root `z` that
  vanishes as the pair closes, so the epsilon on `z` becomes `sqrt(epsilon)`
  on the separation, of order 1e-8 relative.  No Newton step against `chi`
  recovers that.  The guarantee is the weaker and correct one — refining
  never leaves the roots worse than the closed form left them — and a test
  asserts that ordering case by case, which is where the content is.  The
  guard refuses the step outright for about forty per cent of that sweep,
  so which spectrum carries the worst error is decided by the last bit;
  the aggregate comparison is `<=` rather than `<` for that reason, a
  distinction CI found under two of the five Python versions after the
  strict form passed locally.

### Unchanged

- **Two and three flavors**, verified byte-for-byte: 94 records covering both
  backends, scalar calls, stacks straddling every dispatch threshold, scans,
  grids, evolution operators, `earth` and `slabs`, hashing identically to
  1.9.0.
- **Every physical configuration tested.**  A 400 000-point scan through the
  sterile matter resonance at twelve times crust density holds a smallest
  relative eigenvalue gap of 7.7e-4, four orders above where the guard fires;
  both paths stay unitary to 1.6e-13 across it, exactly as before.  Over
  three thousand ordinary random Hermitian spectra the guard changed nothing,
  bit for bit, and on the stiff 3+1 spectrum the refined error stays at
  5.5e-16.

## [1.10.0] - 2026-08-02

**A compiled kernel for four neutrinos.**  `fastkernels` covered two and
three flavors; `oscprob4nu` was pure NumPy whether or not Numba was
installed.  It is now the third member of the family, and it is the one that
gains the most — 18x to 19x on large stacks, against 15x at three flavors.

A minor rather than a patch release, matching 1.6.0 and 1.9.0: it adds
public API, `fastkernels.probabilities_4nu_kernel`.  No existing module
changes behaviour, and nothing changes at all without Numba installed.

### Added

- **`fastkernels.probabilities_4nu_kernel`**, and the compiled machinery
  behind it: `_one_4nu`, which transcribes the whole four-flavor chain for a
  single element — traceless part, the three invariants from traces of
  powers, the quartic by Euler's reduction, the Newton refinement of the
  roots against the matrix, the divided differences of the exponential over
  them, and the Newton-form reconstruction of `U_4` — plus serial and
  parallel runners over a stack.

  Two decisions in it are worth recording.

  - **The kernel refines the roots, and `POLISH_ROOTS` is threaded through
    as an argument** rather than read from module state, which an `@njit`
    function cannot do at call time without recompiling.  `_run` grew an
    `extra` tuple for it, which two and three flavors leave empty.  A test
    pins both settings, so a kernel that ignored the flag would fail
    whichever way it was set.
  - **`chi(psi) = det(psi·1 − H~)` is evaluated by Gaussian elimination with
    partial pivoting**, written out for a 4×4 in a caller-supplied scratch
    buffer — the same factorization LAPACK performs for `numpy.linalg.det`
    on the NumPy path, with no allocation and no call.  It beats
    `numpy.linalg.det` under Numba by 3.5x.

    A Laplace expansion in the six 2×2 minors of the first two rows was
    written first, is **5.9x cheaper still**, and was rejected.  The
    determinant is evaluated *at a root*, where it is meant to vanish: on a
    stiff 3+1 spectrum the true value sits seventeen orders of magnitude
    below the products being summed, so an expansion that cancels them only
    at the end has nothing left, while elimination cancels while the entries
    are still full precision.  Measured against `mpmath` at sixty digits, on
    the clustered roots where `chi'` is 6e-35, the expansion was a thousand
    times the less accurate, refining those roots to 4.4e-15 relative
    against 5.5e-16; on a spectrum whose cluster is 1e-3 wide the gap is
    54x.  `POLISH_ROOTS` tabulates 1.1e-16, and a backend quietly delivering
    forty times that whenever an optional dependency is installed would make
    that table false.

    **What the measurement did not show** is any of this reaching the
    probabilities, and no test here distinguishes the two.  Below `psi·L ~ 1`
    both sit on the one-ulp floor; above it the reconstruction cancels by
    ~1e6 and swamps them, and which scores better is then noise — on the
    stiff spectrum at 1300 km the *rejected* expansion won, 3.1e-11 against
    1.5e-10.  The case for elimination is fidelity to the roots the NumPy
    path computes, not a demonstrated gain in the numbers handed back.  It
    costs about 40% of the kernel's serial runtime.

- **`fastkernels.MIN_BATCH[4] = 1`.**  Measured the way the other two
  thresholds were — alternating the paths through
  `oscprob4nu.probabilities_4nu`, best of nine rounds each.  The kernel
  leads by 15x at a single element, falls to 5x just below
  `PARALLEL_THRESHOLD` where it is still single-threaded, and settles at 18x
  once the threads are in use.  It is never behind, so the threshold is one.

  Note that `worthwhile` already fell back to `MIN_BATCH.get(n, 1)`, so four
  flavors returned `True` at every size before the entry existed.  Adding
  the dispatch without a measured threshold would have enabled the kernel
  everywhere by accident rather than by measurement; the entry is explicit
  so that it is on the record.

- **25 tests**, in `tests/test_fastkernels.py` and one parameterization in
  `tests/test_oscprob4nu.py`.  The `backend` fixture moved to `conftest.py`,
  since the four-flavor module's own suite now needs it too: with a kernel
  installed, every batched assertion there was otherwise being made about
  whichever backend happened to be present.

- **A unit test for `_chi_4nu`**, against `numpy.linalg.det` directly,
  including a matrix whose leading entry vanishes at the evaluation point
  and so requires the pivoting.  Mutation testing found that disabling the
  pivoting left the whole suite green, which follows from the determinant's
  accuracy not surviving into the probabilities: no end-to-end test can see
  it, so it is checked where it is observable.

- **A guard on the new figures**, in `tests/test_documented_figures.py`.
  The four-flavor speedups are stated in `fastkernels` and again in
  `methodology.rst`, and the one-sentence span in `README.md` and
  `quickstart.rst` has to bracket them — which is the shape of every drift
  that module already exists to catch.

### Changed

- **`tests/test_oscprob4nu.py::test_batched_agrees_with_scalar` asks for
  round-off rather than exactness.**  Its bar was `1e-15`, which held only
  while both sides ran the same NumPy code; the batched call is now the
  compiled kernel and the one-by-one calls stay on the NumPy path, and two
  implementations of the same expansion agree to a few ulp — 4.2e-15 on
  probabilities of order one.  It runs on both backends now, which makes it
  the element-by-element comparison of the two.

- **Notebook 09's four-flavor section** no longer says there is no compiled
  kernel for four flavors, and measures the ratio twice: 9.3x on the NumPy
  path, 4.9x with both compiled.  The narrowing is not the algebra getting
  cheaper but the fixed cost of driving it as forty whole-array passes,
  which the kernel does not pay.

### Known limits

- **A stiff 3+1 spectrum makes the two paths differ by ~1e-10**, and this is
  a property of the expansion rather than of either path.  The Newton-form
  reconstruction of `U_4` cancels by ~1e6 there — its four terms reach
  9.7e5 and sum to a matrix of modulus one — so a last-bit difference
  anywhere in that sum arrives at 1e-10.  The amplification predicts 2.3e-10
  and the paths differ by 2.4e-10.  Both stay inside the 1e-9 that
  `POLISH_ROOTS` documents; the equivalence test asserts that, and says why
  it cannot assert round-off.

- **`_polish_roots` is unsafe when two latent roots very nearly coincide,
  on both paths and on `main` before this release.**  Found by mutation
  testing this kernel, not caused by it.  The Newton step divides by
  `chi'(psi_m) = prod_l!=m (psi_m - psi_l)`, guarded only by
  `derivative != 0.0`, which is a knife-edge test: a pair separated by one
  ulp gives a derivative of ~1e-16 and a step that throws the root across
  the spectrum.  On synthetic spectra with pair separations drawn between
  1e-16 and 1e-6 relative, the NumPy path's `_latent_roots` returns roots
  wrong by more than 1e-6 relative for ~10% of them, on `main` unchanged.

  The kernel inherits this faithfully — handed the kernel's roots, the
  NumPy `_polish_roots` returns bit-identically wrong values — but it
  changes *which* inputs trigger it, because the resolvent root that
  vanishes for a degenerate pair is computed to a different last bit and
  `sqrt` near zero amplifies that into the pair's whole separation.

  Not reached by any physical configuration tested: a 400 000-point scan
  through the sterile matter resonance at twelve times crust density holds
  a smallest relative eigenvalue gap of 7.7e-4, and both paths stay unitary
  to 1.6e-13 across it.  A real fix belongs in `_polish_roots` and is its
  own change.

- **A single four-flavor probability still costs ~262 µs**, against 13.5 µs
  at three, because `oscprob4nu` has no scalar closed form and a scalar call
  runs the array machinery on a stack of one.  The kernel does not change
  this: the dispatch excludes scalar calls, which have to keep returning a
  tuple of sixteen.  `oscprob4nu.SMALL_BATCH` is documented as sending short
  stacks through "the scalar path" and is in fact read by nothing.  Left
  alone deliberately, as its own decision.

## [1.9.0] - 2026-08-01

**Four-neutrino oscillations.**  The Ohlsson–Snellman method is carried to
`n = 4` through the SU(4) algebra, which brings 3+1 sterile scenarios into
scope — and which is the last `n` where the method exists in closed form at
all.

A minor rather than a patch release, matching 1.8.0, which added `slabs` and
`earth`: this adds public API rather than correcting anything.  No existing
module changes behaviour.

### Added

- **`oscprob4nu`**, the four-flavor expansion.  Everything structural carries
  over from three flavors — expand in the fifteen generalized Gell-Mann
  matrices, factor out the global phase, interpolate the exponential over the
  latent roots, read the probabilities off `|U|²`.  Three ingredients are new:

  - **A third invariant.**  SU(4) has rank three, so the traceless part
    carries `I2`, `I3` and `I4`, and the cubic characteristic equation becomes
    a quartic.  All three come from traces of powers, which avoids ever
    building the SU(4) `d` tensor — a 15×15×15 table nothing else would need.
  - **A quartic that still solves in closed form.**  Euler's reduction gives a
    resolvent cubic whose roots are real and non-negative because the
    Hamiltonian is Hermitian, so the same trigonometric construction that
    `oscprob3nu.psi_roots` uses solves it.  The SU(3) machinery is literally
    nested inside the SU(4) solution.
  - **A longer star-product tower.**  The three-flavor identity
    `(h*h)*h = |h|² h/3` is a Cayley–Hamilton accident of `n = 3` and is false
    at `n = 4` — 37% off on a random Hamiltonian — so `((h*h)*h)_a` enters the
    `u_a` as independent data.  A test asserts it stays false, since if it
    ever held the extra term would be dead weight.

- **`hamiltonians4nu`**, sample 3+1 Hamiltonians: the 4×4 mixing matrix with
  three extra angles, the energy-independent vacuum term, matter, and matter
  with NSI.  As at three flavors these are examples, not limitations.

- **`globaldefs.VNC_EARTH_CRUST`.**  The neutral-current potential is
  flavor-universal across the active states, so at three flavors it is
  proportional to the identity and drops out entirely — which is why it has
  never been needed.  With a sterile state it does not drop out: removing it
  from all four costs only a global phase and leaves `-V_NC` on the sterile
  entry, positive and half of `V_CC` in the crust.  Getting that entry wrong is
  the four-flavor analogue of the antineutrino sign trap — invisible in vacuum,
  and it moves the resonance — so a test pins its sign and size.

- **Notebook 16**, working a 3+1 scenario through: the sixteen probabilities,
  the sterile entry in the matter potential, a short-baseline energy scan, the
  sterile matter resonance through the Earth, the three new algebraic
  ingredients evaluated live, the accuracy comparison, and the decoupling
  check against `oscprob3nu`.  Two gallery figures come from it.

- **`tests/test_oscprob4nu.py`**, 30 tests.  The strongest is the last kind:
  switching the three sterile angles off must reproduce `oscprob3nu` exactly,
  which it does to 7e-12 in vacuum and in matter — an independent module,
  written years earlier against a different algebra, computing the same
  numbers.

### Changed

- **The documentation no longer says the library is limited to two and three
  flavors.**  There were more than a dozen such claims, including a bullet in
  "What it is not" reading *"Not a four-flavor code … a sterile fourth is
  outside what they cover"*.  The paper's title is left exactly as published
  in all four places it appears; the extension is described as an extension.

- **When to use Magnus is now stated up front**, in both `README.md` and the
  documentation landing page, as a table rather than a remark buried in "What
  it does not do".  The rule is one line — reach for Magnus when the
  Hamiltonian varies appreciably over an oscillation length — and the table
  gives five cases with the reason for each, including the one neither tool
  handles, which is open-system evolution needing Lindblad.

- **`methodology.rst` gains a four-flavor section**, including *why the method
  stops at four*: the eigenvalues come from solving the characteristic
  polynomial in radicals, and Abel–Ruffini says degree five has no such
  solution, `S₅` not being soluble where `S₂`, `S₃` and `S₄` are.  A theorem,
  not a gap.  The philosophy survives even where the closed form does not —
  numerically computed eigenvalues feed the same machinery at any `n`; it
  simply stops being a closed form, which is this library's reason to exist.

### Notes on accuracy

Worth reading before using four flavors in anger, and worth keeping in
proportion: **none of what follows is near a measurable effect**, since
probabilities meet data at the per-cent level.

A generic four-flavor Hamiltonian agrees with `scipy.linalg.expm` to 3.7e-14.
A *stiff* 3+1 spectrum — an eV-scale `Δm²₄₁`, so eigenvalues spanning four
orders of magnitude with three clustered — is different: forming `I2`, `I3`
and `I4` in double precision compresses the 4×4 matrix into three numbers and
loses what separates the cluster.  Perturbing the three invariants at the
1e-16 level moves the roots by 6e-11 relative, which is ordinary
ill-conditioning of polynomial roots against their coefficients, so no better
root-finder helps.  Deflating the quartic to a cubic first was tried, and
does not.

The roots are therefore refined against the *matrix*, by one Newton step on
`χ(ψ) = det(ψ𝟙 - H̃)`.  Measured against `mpmath` at fifty decimal digits:

| Strategy for the roots | Relative error | Cost, 200k points |
|---|---|---|
| Closed form alone | 8.3e-11 | 0.17 s |
| **Closed form + one Newton step** | **1.1e-16** | 0.41 s |
| `numpy.linalg.eigvalsh` | 7.4e-16 | 0.17 s |
| Closed form in `numpy.longdouble` | 4.5e-11 | 0.43 s |

The Newton step is about seven times *more accurate than LAPACK*, since
`eigvalsh` reduces by similarity transforms each carrying a backward error of
order `ε‖H‖`.  Extended precision was rejected for buying under a digit, for
being slower, and for silently being `float64` on Apple Silicon and Windows.
`eigvalsh` was rejected for being less accurate and for meaning the module
would take its eigenvalues from LAPACK, which is the one thing this library
exists not to do.

In probabilities that is 5e-7 unrefined against 1e-9 refined.  The reasons to
want the smaller number are the exactness claim, error accumulating when
`slabs` and `earth` compose operators across many layers, and a regression
suite with no room for a bug to hide in.  It costs about 40% of the runtime,
bringing the four-flavor closed form to parity with a batched `eigh` rather
than ahead of it.  `oscprob4nu.POLISH_ROOTS` records the trade and switches
it off.

**Applying the refinement selectively was measured and rejected.**  Skipping
it where the spectrum is not stiff — the way `SMALL_BATCH` and `MIN_BATCH`
dispatch on a threshold — would recover that 40%.  Two criteria were tried on
6300 Hamiltonians: the gap-based amplification perturbation theory suggests,
which misjudges doubly degenerate pairs by four orders of magnitude, and a
matrix residual comparing `prod(ψ_m)` with `det(H̃)`, which is one constraint
on four roots and misses errors that cancel in the product.  Neither can
safely skip a single element.  The reason is structural: a criterion complete
enough to certify four roots must evaluate `χ` at four roots, which *is* the
refinement — the check and the fix are the same computation.  And since a 3+1
scan is stiff at every point, even a working criterion would skip nothing on
the workload that motivates the module.

### Fixed

- **Degenerate spectra at four flavors.**  The first implementation
  reconstructed `U₄` by solving a Vandermonde system for the Cayley–Hamilton
  coefficients, which is singular the moment two latent roots coincide — not
  an exotic case, but one that includes a Hamiltonian proportional to the
  identity, a zero Hamiltonian, and any triply degenerate spectrum, all of
  which raised `LinAlgError`.  It is replaced by Newton interpolation with
  divided differences, where a repeated node is a derivative and for the
  exponential that derivative is known exactly.  No special branch, no
  tolerance-dependent switch between formulas, and six degenerate spectra now
  agree with `expm` to 1e-12.

## [1.8.6] - 2026-08-01

A pre-publication audit, before the first upload to PyPI.  Nothing here
changes a probability; it is documentation that had drifted from the code,
metadata that had gone stale, and one extra that could not do what it
promised.

Two of the entries below deserve reading together, because they are the same
decision applied twice.  The file tree and the performance figures were both
kept in more than one place, both drifted, and both were caught by hand
rather than by the suite.  Neither copy is deleted — `README.md` is the PyPI
long description and has to stand alone — so instead the duplication is made
safe: the tree is now generated from one table, and the figures are now
checked against one table.  Duplication that cannot drift is not the same
problem as duplication that can.

### Fixed

- **Three recipes on the "Numerical recipes" page stopped short of what they
  promised.**  The page opens by saying each recipe is "a few lines and a
  figure", that the code shown is the same calculation as the notebook linked
  beside it, and that anything short enough to run is executed at build time
  and shows its output.

  "Between two places on the Earth" said "the chord between two named sites,
  and the probability along it", then printed only chords.  The probability was
  the point, and `earth.probabilities_3nu_between_locations` — the reason the
  named-location table exists at all — appeared nowhere on the page.  It is now
  called, and the section links notebook 07, which it had also been missing
  while every other recipe carried a link.

  "An oscillogram" was a static `code-block` whose last statement assigned
  `grid` and stopped, so it produced no output even in principle.  Nothing
  justified the exception: the 240×240 grid evaluates in 0.14 s.  It now runs
  and prints its shape and range.

  "Mass ordering and the octant" promised two open questions in its title and
  printed two mass-squared differences, with no probability anywhere and no
  mention of the octant, beneath a figure captioned with a claim the code did
  not support.  Both are now computed, following notebook 12.  The octant block
  is evaluated at 5 GeV rather than at the 2.5 GeV used for the ordering: at
  2.5 GeV the two octants give `P_mumu` of 0.019 against 0.0067, which is the
  deep disappearance minimum rather than an error, but would have read as
  contradicting the octant degeneracy it illustrates.  At 5 GeV `P_mumu` is
  0.4816 against 0.4768 — one per cent apart — while `P_mue` is 0.0245 against
  0.0299, twenty-two per cent apart.

- **The `docs` extra could not build the documentation.**  `pip install -e
  '.[docs]'` installed four of the ten packages `cd docs && make html` needs,
  and failed on the first missing Sphinx extension.  `conf.py` loads
  `sphinx-copybutton`, `sphinxcontrib-bibtex` and `jupyter-sphinx`, sets
  `html_theme` to `sphinx_rtd_theme`, and executes the narrative blocks in a
  real kernel, which needs `ipykernel`, `matplotlib` and `scipy`.  Nothing
  caught this because nothing uses the extra — both documentation jobs install
  `docs/requirements.txt`, which has always been complete.  The two lists are
  now the same set, bar `numpy`, which the extra omits as already a hard
  dependency.

- **`methodology.rst` contradicted itself on the scalar timing.**  The page
  states it twice, and the 1.8.4 reconciliation updated only the second
  occurrence, so one page has since said "about thirteen microseconds" and
  "about sixteen microseconds" a hundred lines apart.

- **`fastkernels`'s "Routine listings" omitted `available` and `worthwhile`**,
  two of its four public functions, and `worthwhile` is the one that decides
  whether the compiled path is taken at all.  The other seven modules were
  checked the same way, by parsing each listing against the public functions
  defined in the file; this was the only gap, in both directions.

- **`index.rst` contradicted its own benchmark table.**  It described the gain
  from passing arrays as "roughly 25 to 60 times" three lines above a table
  whose rows give 21×, 23×, 37× and 99×.  `README.md`, which carries the same
  table, had it right at "20–90×".  The page now says 20 to 90 as well, and a
  test derives the span from the table rather than trusting the prose.

### Added

- **`tests/test_documented_figures.py`, so the repeated performance figures
  cannot drift apart again.**  The same measurements appear in `README.md`,
  `index.rst`, `quickstart.rst`, `methodology.rst` and `fastkernels`, and they
  have to: the README is the PyPI long description and must stand alone.  The
  figures are now held once in the test and checked against every document
  that quotes them — the scalar timing in digits or in words, whichever that
  file uses; the twelve numbers of the 2000-point table, compared as the
  strings the documents print, so that "0.20 ms" drifting to "0.2 ms" is
  caught too; the quoted speedup span against the table's own ratios; and an
  explicit refusal of the superseded 13 and 16, so a stale copy-paste cannot
  reintroduce one quietly.  Each failure names the file that disagrees and the
  value it should carry.

### Changed

- **The scalar timings are re-measured: about 8 microseconds for three flavors
  and 1 for two, against the 16 and 2 quoted since 1.8.4.**  The figure
  appeared in five files and seven places, all agreeing with each other and
  none reproducing.  Measured after a warm-up loop, fastest of fifteen rounds
  of two thousand calls, timing overhead subtracted, repeated three times end
  to end: 8.18, 7.99 and 8.12 µs for three flavors, 1.10, 1.11 and 1.15 for
  two.  The spread within a case is under four per cent, far tighter than the
  gap to the figure replaced, so this is not the run-to-run variation the
  surrounding prose already warns about.

  The four batched speedups measured at the same 1.8.4 sitting are **not**
  re-checked here.  They are ratios between two paths measured together, which
  is more robust than an absolute time, and notebook 09 measures them on
  whatever machine runs it.

- **The file tree is generated rather than maintained by hand in two places.**
  It appears in `README.md` and in `installation.rst`, and both were edited by
  hand.  The old tests compared the two and checked that every tracked file's
  *basename* appeared in one of them, which left three ways to be wrong: the
  comparison was between stripped lines, so indentation could drift unnoticed;
  matching on basenames meant a file listed under the wrong directory still
  passed; and adding a file meant editing two documents, in the right place,
  with the box-drawing characters aligned.

  Both are now rendered from one ordered table in `tests/test_file_tree.py`,
  which pairs each path with the comment beside it.  Order and comments stay
  written down, since neither can be derived from the filesystem; membership
  is derived, and `test_tree_matches_git` requires the file entries to be
  exactly `git ls-files` in both directions.  `python tests/test_file_tree.py
  --write` updates both documents; with no arguments it reports whether they
  are current.  The renderer was written against the existing tree, and
  reproduces all 94 lines of both copies byte for byte, so the change altered
  no rendered output at all.

- **`docs/requirements.txt` is now just `-e .[docs]`.**  Two lists of
  documentation dependencies, which this same release had just finished making
  agree — and the reason they disagreed for so long is that nothing installed
  the extra, so nothing noticed.  There is now one list, in `pyproject.toml`.
  The editable install that requirements.txt performs makes the separate
  `pip install -e .` step redundant, so the two steps are merged into one in
  `lint.yml` and `pages.yml`, keeping the note about *why* it must be editable:
  jupyter-sphinx executes the narrative blocks in a real kernel, a separate
  process that `conf.py`'s `sys.path.insert` never reaches.

  Verified in a fresh virtualenv: `pip install -r docs/requirements.txt` from
  the repository root installs the package plus all ten documentation
  packages, `oscprob3nu.__file__` resolves into the working tree, and Sphinx
  then builds clean under `-W` with no separate package install.

- **`README.md` no longer lists the optional extras twice.**  "Requirements"
  gave two of them as install commands and "Installation", twenty lines below,
  gave all four again.  Requirements now says what each task needs and which
  extra provides it, as a table; Installation owns the commands.  No fact is
  dropped, and the `docs` row is new — Requirements had never mentioned it.

### Removed

- **The module `__version__` strings.**  Eight modules carried one by hand:
  five said `"1.1"`, `fastkernels` said `"1.6"`, `slabs` and `earth` said
  `"1.8"`, against a project at 1.8.5.  No release commit had ever touched
  them, so all eight had drifted — `oscprob3nu` claimed `"1.1"` while holding
  the vectorised batch path, the degeneracy handling and the Numba dispatch,
  none of which existed at 1.1.

  Syncing them would add eight more places every release must update, which is
  what produced the drift.  Deriving them from `importlib.metadata` would raise
  `PackageNotFoundError` for a copy used off `sys.path` without being
  installed, which is what the scripts in `examples/` do and what the paper
  tells readers to do.  So they are deleted.  Nothing in the repository read
  them: the documentation takes its version from `docs/source/conf.py` and the
  package metadata from `pyproject.toml`.  `__author__` and `__email__` are
  kept — still true, and they do not drift.

- **`from __future__ import print_function`, from the eight modules.**  A
  Python 2 shim in a package whose floor is 3.9, and a no-op on every
  interpreter this project has ever supported.  The nine scripts in `examples/`
  keep theirs: they are published alongside
  [arXiv:1904.12391](https://arxiv.org/abs/1904.12391), and the repository
  already treats their period character as deliberate — the star imports they
  open with are ruff-exempted in `pyproject.toml` for that reason.

## [1.8.5] - 2026-08-01

### Changed

- **`test/` is now `examples/`.**  Having `test/` and `tests/` side by side was
  a standing invitation to open the wrong one, and tab-completion could not
  tell them apart.  The new name says what the directory holds: nine runnable
  scripts, the ones `README.md` walks through line by line.

  Nothing functional depended on the old name.  `pytest` never collected from
  it (`testpaths = ["tests"]`), so the clash was cosmetic rather than real, and
  the scripts locate `src/` relative to themselves, so they run from the new
  location unchanged.  `git mv` keeps their history.

  The directory is called `test/` in version 1.0.0 of the code and in version 2
  of [arXiv:1904.12391](https://arxiv.org/abs/1904.12391).  `README.md` and the
  installation and quickstart pages say so, for anyone arriving from the paper.
  The layout had already moved away from what the paper describes — 1.8.1
  removed `run_testsuite.py` and the four plotting modules, and `src/` has
  gained three modules since — so this is one more difference in a list, not a
  new kind of divergence.

  Entries below this one still say `test/`.  They are dated records of what was
  true when written, and are left alone.

## [1.8.4] - 2026-08-01

### Changed

- **Docstring examples are executed when the documentation is built.**  Every
  `Examples` block is now a `.. jupyter-execute::` directive, following the
  pattern the sibling Magnus package uses: self-contained, importing what it
  needs, so it can be copied straight out of the page.  The API reference no
  longer shows `>>>` prompts with results pasted beside them — 42 blocks, 0
  prompts left.

  Each converted block was checked against the output its doctest documented,
  and all 42 reproduce it exactly, so nothing about what the examples *do*
  changed.

  `tests/test_docstrings.py` is repurposed rather than left hollow.  It now
  extracts and executes the same blocks on every supported Python — the
  documentation is built by one job on one interpreter, and an example that
  works on 3.12 and not on 3.9 would otherwise reach the published page
  unnoticed.  A second test refuses any `>>>` example that creeps back, since
  the run-test would not see it.

- **The performance figures are re-measured and reconciled.**  Three places
  quoted different numbers for the same four benchmarks: `fastkernels` said
  ~20x, ~15x, ~5x and ~4x; the 1.6.0 changelog entry said 12.5×, 15.0×, 2.3×
  and 3.7×; and neither matched what the code does now, which is ~15x, ~9x,
  ~3.5x and ~1.5x.

  The live claims — in `fastkernels`, the README, `index.rst`, `quickstart.rst`
  and `methodology.rst` — now carry the same measured figures, taken best of
  seven with the two paths interleaved, and say plainly that they move by tens
  of per cent between runs.  The 1.6.0 entry below is left as it was written:
  it records what was measured then, and rewriting it would falsify the record
  rather than correct it.

  The scalar timings drifted too, and are corrected: about 16 µs for three
  flavors and 2 µs for two, against the 13 µs and 1.3 µs previously quoted.

  Notebook 09 measures the same comparison on whatever machine runs it, and the
  documentation now points at it as the figure to trust.

### Fixed

- Two `test/` examples were referenced in `README.md` under names that do not
  exist — `example_2nu_vacuum_coefficients.py` and its three-flavor
  counterpart, where the files are `..._coeffs.py`.

- Nothing executed the worked examples in `test/`, which is why the broken
  references went unnoticed; a step in `lint.yml` now runs all nine.  They are
  the examples the paper refers to and the ones `README.md` walks through, so
  they are kept rather than removed, but they are now checked.

- `README.md` still described the notebooks as nine; there are fifteen.

## [1.8.3] - 2026-08-01

### Changed

- Installation instructions in `README.md` and `docs/source/installation.rst`
  now lead with `pip install nuoscprobexact` and are split three ways: from
  PyPI, from a clone of GitHub, and without installing anything at all.

  Cloning was previously the only route described, which made the ordinary
  case --- somebody who wants to use the library rather than work on it ---
  read the longest set of instructions.  The clone is still documented, and
  is still what you want for the notebooks, the worked examples from the
  paper, the regression suite, or an unreleased version.

  The third route is kept because it is genuinely supported: the two core
  modules need only `numpy` and `cmath`, so copying one into a project of
  your own works, and is what the paper assumes.

  **`pip install nuoscprobexact` does not work until the first release is
  published.**  The instructions describe the intended state; publishing is
  what makes them true.

## [1.8.2] - 2026-08-01

A pass over the README and the documentation, so that both say what the
library does now rather than what it did in 2019.

### Added

- A **What you can compute** gallery at the top of `README.md`: eight figures,
  each linking to the notebook that drew it.  The figures are extracted from
  the executed notebooks by `notebooks/make_notebooks.py`, so there is one
  piece of code behind each picture rather than a second copy that can drift.

- `docs/source/recipes.rst`, a numerical-recipes page collecting the same
  material with runnable snippets.  Eight of its examples are executed when the
  documentation is built, so the numbers on the page are produced by the code
  being documented.

- A statement of scope in both `README.md` and the documentation landing page —
  what the library does, and what it does not.  The second list is the one that
  was missing: no continuously varying Hamiltonians, no fourth flavor, no
  fluxes or cross sections, no fitting.

- `notebooks/make_notebooks.py`, the generator that produces and executes the
  fifteen notebooks and extracts the gallery figures.  It was scaffolding
  outside the repository until now.

### Changed

- The narrative documentation pages execute their examples at build time
  through `jupyter-sphinx`, rather than quoting output written by hand beside
  them.  `quickstart.rst` had one such block; its numbers were correct, and are
  now generated.

  Docstrings keep their `>>>` examples deliberately.  Those are run as doctests
  by the regression suite on every supported Python, and they are what `help()`
  shows in a terminal; converting them would have moved the check to a single
  interpreter and made the interactive reading worse.

- The scope description in `README.md` and `index.rst` now covers `slabs`,
  `earth`, PREM and the named locations, none of which existed when that text
  was written.

### Removed

- `Created:` and `Last modified:` lines from every module docstring.  They are
  version control's job, they had already drifted, and Sphinx rendered them
  into the API reference, where each appeared eight times.

## [1.8.1] - 2026-08-01

Seven worked notebooks, a logo, and `fig/` out of version control.  No change
to the library.

### Added

- `notebooks/`, fifteen Jupyter notebooks numbered in reading order: the
  basics, vacuum oscillations, matter and NSI and LIV, oscillograms,
  bi-probability plots, the Earth and PREM, probabilities through the Earth,
  unusual matter profiles, performance, the paper's own figures, exact versus
  the textbook approximations, mass ordering and the octant, antineutrinos,
  solar neutrinos, and numerical edge cases.  They carry their figures inline, so they render
  on GitHub without being run.

  Two of them go beyond what the old figure suite covered.  The unusual-profile
  notebook builds castle-wall, serrated and shuffled profiles by hand and holds
  their mean density fixed, so that any difference in the probabilities is due
  to the arrangement of the matter alone — including the parametric enhancement
  that a periodic profile produces.  The performance notebook measures, on the
  machine that runs it, the cost of looping against broadcasting and of the
  NumPy path against the compiled kernel.

  The solar notebook is as much a warning as a demonstration.  It validates
  the slab machinery against the analytic adiabatic MSW result — averaged over
  a narrow energy band, the two agree to 0.0009 — and then shows why the
  approach is impractical there anyway: a single averaged point costs 800 000
  slab evaluations, and a realistic calculation multiplies that by the
  production region, the spectrum and a third flavor.  It points at the
  Magnus package for profiles that vary smoothly over an oscillation length.

- A `notebooks` extra — Jupyter, matplotlib and scipy — and a `notebooks` job
  in `lint.yml` that executes every
  one of them.  A notebook is documentation that claims to work, and stored
  outputs make that claim persuasive without making it true — they were
  correct whenever the notebook was last run, which may predate the change
  that broke it.  The job also refuses a notebook stripped of its outputs,
  which would execute cleanly while showing a reader nothing.

- The project logo, wired in as `html_logo`, at the top of the documentation
  sidebar.

### Removed

- `run_testsuite.py` and the four plotting modules it drove
  (`test/oscprob2nu_plot.py`, `test/oscprob3nu_plot.py`, and the two
  `*_plotpaper.py`).  The notebooks cover the same figures and show them
  without anyone having to run anything or go looking in a directory
  afterwards, and the plotting modules had no other caller.  The worked
  examples in `test/`, which the paper refers to, are untouched.

  Two code blocks in `README.md` went with them.  Both called a module named
  `oscprob3nu_tests`, which has not existed under that name for years, so they
  had been broken well before this release; they now point at the notebook that
  draws the same curve.

- `fig/`, which is no longer tracked or written.  It held one committed
  `.gitignore` whose only job was to keep an empty directory alive.

- The `plots` extra, which existed for the plotting modules and had nothing
  left to install for.  `matplotlib` is still available through the new
  `notebooks` extra.

## [1.8.0] - 2026-08-01

Piecewise-constant matter.  The exact expansions assume a Hamiltonian that
does not change along the trajectory; a neutrino crossing the Earth does not
have one.  Two new modules close that gap without giving up exactness: the
path is cut into slabs, each is solved exactly, and the operators are
multiplied.  Within a slab nothing is approximated.

Nothing in the existing modules changed behaviour.

### Added

- `slabs`, which propagates across a sequence of adjacent slabs of arbitrary
  width and Hamiltonian, for two and three flavors:
  `evolution_operator_2nu_slabs`, `evolution_operator_3nu_slabs`,
  `probabilities_2nu_slabs`, `probabilities_3nu_slabs`.

  The per-slab operators are evaluated in one batched call, so `n` slabs cost
  one vectorised evaluation plus `n-1` matrix products rather than `n`
  separate evaluations.

  Note that each slab drops the phase carried by the trace of its
  Hamiltonian, as the single-slab routines do.  Their product therefore
  differs from the product of full matrix exponentials by one overall phase,
  which leaves every probability unchanged but matters to anyone comparing
  the operator itself against `scipy.linalg.expm`.

- `earth`, which builds those slabs from the Preliminary Reference Earth
  Model: `density_prem`, `matter_potential`, the chord geometry
  (`distance_traveled_inside_earth`, `earth_radial_distance_from_depth`,
  `prem_layer_edges_along_chord`, `chord_length_inside_earth`,
  `costhz_between_points_on_surface`), `earth_slabs`, and the probabilities
  across the Earth for a given zenith angle (`probabilities_2nu_earth`,
  `probabilities_3nu_earth`) or between two named locations
  (`probabilities_2nu_between_locations`,
  `probabilities_3nu_between_locations`).

  Two things decide where the slabs are cut.  *Between* shells the density
  jumps, so the chord is first split at every PREM boundary crossing — no
  amount of subdivision recovers a discontinuity that straddles a slab.
  That gives a set of chord segments, which are not the same as shells: a
  chord enters and leaves each shell it reaches, so a diametric chord has 19
  segments across 10 shells.  *Within* a shell the density varies smoothly,
  by as much as 21% over a single 2200 km mantle segment, so each segment is
  divided further into `n_slabs_per_segment` equal sub-slabs with the density
  taken at the midpoint.

  Midpoint sampling converges at second order — measured, not assumed: past
  32 sub-slabs per segment each doubling cuts the error by about four.

  The 15 predefined locations are the same set as the sibling Magnus
  package, so a trajectory named in one can be reproduced in the other.

- `globaldefs.EARTH_RADIUS`.

- 68 tests for the two modules.  The ones that matter are the independent
  checks: slab composition against a product of `scipy.linalg.expm` factors,
  splitting one Hamiltonian into many slabs reproducing the single-slab
  answer, PREM integrating to the mass of the Earth to 0.02%,
  `matter_potential` reproducing `globaldefs.VCC_EARTH_CRUST` from the crust
  density, and a uniform-density Earth reproducing the ordinary
  constant-Hamiltonian result.

## [1.7.0] - 2026-07-31

Continuous integration, an automatic coverage gate, and linting.  Nothing in
the library changed: no module in `src/` was touched except for a coverage
pragma, the 279 tests are the same 279 tests, and every probability this code
computes is the one it computed before.  What changed is that all of it is now
checked on every push instead of whenever someone remembered.

The one thing users will notice is the supported Python range.

### Added

- Four GitHub Actions workflows, in `.github/workflows/`.

  | Workflow | What it does |
  |---|---|
  | `tests.yml` | The suite on Python 3.9–3.13, plus a job for each of the three backend configurations, plus coverage |
  | `lint.yml` | `ruff check`, and a Sphinx build under `-W --keep-going` |
  | `pages.yml` | Builds and deploys the documentation to GitHub Pages on push to `main` |
  | `publish.yml` | Builds and publishes to PyPI on a published GitHub Release |

  The three backend configurations are separate jobs because they are separate
  claims: that the compiled kernels work, that the NumPy path gives the same
  answers across the whole suite, and that a plain `pip install` with no
  `numba` still works.  Each of the two backend jobs asserts its own
  preconditions, so neither can quietly become a duplicate of another and keep
  reporting success.

  `publish.yml` can also rehearse a release against TestPyPI, run manually from
  the Actions tab.  Worth doing because the real upload cannot be retried: PyPI
  refuses a second upload of a version permanently, even after the file is
  deleted.  The event picks the index — a published Release goes to PyPI, a
  manual run goes to TestPyPI — so neither can reach the other, and a rehearsal
  is stamped `.devN` so it can be repeated as often as needed.

- An automatic coverage gate, configured in `[tool.coverage]` in
  `pyproject.toml` so that a local `pytest --cov` measures and gates exactly as
  CI does.  Branch coverage is on, and the floor is 98%.

  The `@njit` kernels in `fastkernels` are excluded, because Numba compiles
  them to machine code that `coverage.py` cannot trace at all — they were
  reported as 118 missing lines that no test could ever close, despite being
  exercised against the NumPy path on every run.  Excluding the decorator
  rather than omitting the file keeps `available`, `worthwhile`, `_run` and
  both public wrappers measured, which is where a dispatch bug would hide.
  Coverage of what can be traced is 100%.

- A `[tool.ruff]` configuration and `ruff check` in CI.  The rule selection is
  pinned rather than left to ruff's default, which moves between releases.
  Three codebase-wide conventions are exempted with the reason recorded at
  each: `l` for the baseline, one-line guards in the `d_ijk` lookup table, and
  the star imports in the paper's worked examples.

- `pytest-cov` in the `test` extra, and version classifiers for Python 3.9
  through 3.13.

- Seven badges in `README.md` and `docs/source/index.rst`: tests, Code Quality,
  codecov, Documentation, PyPI, Downloads, and Code style: ruff.  The PyPI and
  Downloads badges will not resolve until the first release is published.

- A `.. versionadded::` directive on all 32 public functions, so the API
  reference says when each one entered the library.  The answer is mostly
  1.0.0: 28 of them have been there since the first release in 2019, and the
  only later additions are the four that came with the Numba backend in 1.6.0.
  That the public surface has been stable for seven years is worth being able
  to see at a glance.

- 57 `.. versionchanged::` directives across 25 of those functions, recording
  what changed and when, for every release that changed how a function behaves
  or how fast it runs.

  33 of them are changes a caller can observe: a changed signature, a newly
  accepted input type, a changed return type, or a changed returned value.
  Most are 1.1.0, where the audit corrected results; the rest are the batched
  interface in 1.2.0 and the batched Hamiltonian builders in 1.3.0.

  The other 24 are the performance releases — 1.4.0, 1.5.0 and 1.6.0 — which
  rewrote much of the library without moving a single number.  Each says so:
  the 1.4.0 directives note that all 42 figures are byte-for-byte those of
  1.3.0, and the 1.5.0 ones that the probabilities agree with 1.4.0 to 1.6e-13.

  Speedups are given as ratios, not absolute timings.  The 1.4.0 and 1.5.0
  tables in this file were measured in separate sessions and do not chain —
  1.4.0 ends at 14.5 µs where 1.5.0 begins at 40.4 µs — so a pair of absolute
  figures in a docstring would invite a comparison that is not valid.  The
  ratio is what each release actually established.

### Changed

- **`requires-python` is now `>=3.9`, raised from `>=3.7`.**  The old floor was
  never tested and its stated justification — the `@` operator — has held since
  Python 3.5, so it did not pin 3.7 or anything else.  3.9 is what the code
  actually needs: the batched paths call `numpy.broadcast_shapes`, which
  arrived in NumPy 1.20, and 3.9 is the oldest interpreter for which the
  optional `numba` backend still has a wheel.  All five supported versions are
  now in the CI matrix.

### Removed

- Eight unused imports, in three of the worked examples in `test/` and two
  modules in `tests/`.  These were what running a linter for the first time
  actually found.

## [1.6.0] - 2026-07-31

An optional Numba backend for the batched paths, and a shortcut for very short
stacks.  The library's only required dependency is still NumPy, and the
results are unchanged either way: all 42 figures generated by
`run_testsuite.py` remain byte-for-byte identical, and the two backends agree
to round-off.

### Added

- `fastkernels`, which tries to import Numba at module scope.  If it is
  present the batched expansions are compiled into fused machine-code loops
  spread over the available cores; if it is not, `HAVE_NUMBA` is `False`,
  nothing else in the module is defined, and the NumPy path is used exactly as
  before.  Install with `pip install "nuoscprobexact[fast]"`.

  Measured against the NumPy path:

  | Stack | Speedup |
  |---|---|
  | 200 000 energies, three flavors | **12.5×** |
  | 20 000 energies, three flavors | **15.0×** |
  | 200 000 baselines, two flavors | 3.7× |
  | 100 × 100 oscillogram | 2.3× |

  The gain comes from replacing roughly fifteen passes over N-element arrays,
  each writing a temporary the next reads back, with one loop that keeps its
  intermediates in registers.

  The kernels are declared `cache=True`, so the few seconds of first-call
  compilation are paid once per machine; later runs load from disk in
  milliseconds. Stacks below `PARALLEL_THRESHOLD` run single-threaded, since
  waking the thread pool costs more than it saves there.

  Two things stay on the NumPy path deliberately. The **scalar** routine is not
  compiled: one probability takes about thirteen microseconds, which is not
  worth a compilation pause on a first call. And `fastkernels.USE_NUMBA =
  False` forces the NumPy path at any time.

- `fastkernels.worthwhile`, which declines the kernel for stacks it would not
  speed up. This was not in the first draft, and the first draft was wrong for
  it: with the backend installed, a two-flavor baseline scan of 2000 points ran
  **2.6× slower** than without. That expansion reduces to a square root and a
  sine per element, which NumPy already does about as well as compiled code
  can, and the kernel additionally has to materialise the Hamiltonian stack —
  for a fixed-Hamiltonian scan, the same matrix repeated, costing 2.5 ms to
  copy at 200 000 points.

  The thresholds in `MIN_BATCH` were measured by alternating the two paths and
  taking the best of nine rounds each: three flavors wins at every size, two
  flavors from about fifty thousand elements. A backend that is sometimes
  slower than the path it replaces is worse than no backend.

- Tests for both backends in `tests/test_fastkernels.py`. Every equivalence
  check runs with the compiled path forced off as well as on, so the NumPy
  path is exercised even where Numba is installed rather than rotting
  unnoticed; the whole suite passes both ways, and in a third configuration
  where Numba is genuinely unimportable. The kernel carries its own copy of
  the degeneracy branch, so that is checked separately against the scalar
  path, and agrees exactly.

- A `kernel_spy` fixture, which counts entries into the compiled kernels so a
  test that means to exercise them can assert that it did.

  This exists because one did not. Raising the two-flavor dispatch threshold
  to fifty thousand elements meant the two-flavor equivalence tests, which run
  stacks of at most a thousand, stopped reaching the kernel altogether: they
  were comparing the NumPy path against itself, and passing. Instrumenting the
  suite showed the two-flavor kernel entered **once** in a whole run against
  thirty-nine entries for the three-flavor one. A test that passes for the
  wrong reason is worse than a missing one, because it reads as coverage.

- `tests/test_physical_scales.py`, which compares the evaluation paths at the
  scales the library is actually used at. Everything else compares them on
  random Hermitian matrices with entries of order one and baselines of order
  ten; a real Hamiltonian has entries around 10⁻¹³ eV and a baseline around
  10¹³ eV⁻¹, and agreement in one regime does not imply agreement in the
  other. These run the bundled vacuum, matter, NSI and LIV Hamiltonians at
  NuFit parameters over energies from 10 MeV to 100 GeV and baselines from 1
  km to 10⁷ km, on both backends, including an oscillogram and a
  near-degenerate spectrum split by one part in 10¹⁴. Agreement is 5e-15.

- Tests for three branches that measurement showed were never executed: the
  fixed-Hamiltonian matrix product and its axis padding, which live on the
  NumPy path and which a default run sends to the kernel instead; `psi_roots`
  returning zeros for a vanishing invariant; and the expansion core
  recomputing the star product when not given one.

  Line *and branch* coverage of the four core modules is now **100% on both
  backends** — 170 branches in `oscprob3nu` alone — where line coverage was
  98% and 99%. The compiled kernels cannot be traced by `coverage.py`, so the
  spy counts entries into them instead.

  `fastkernels` is also brought under the docstring and annotation guards that
  every other module was already subject to, which immediately caught an
  unannotated helper in it.

- A **Performance** section in the README, and one on the documentation
  landing page, built around the two things a user can act on: pass arrays
  instead of looping, which is the larger win and needs no extra dependency,
  and install the optional extra if the scans are large. `methodology.rst`
  gains sections on the short-stack shortcut and on the compiled backend,
  including why it is conditional. All of them state the two-flavor caveat
  rather than quoting only the flattering rows.

### Changed

- Stacks of at most `SMALL_BATCH` elements are evaluated one at a time through
  the scalar path. A batched call carries a couple of hundred microseconds of
  fixed cost whatever its length, so for a handful of points it spends more on
  the machinery than the scalar path spends on the whole job: a single
  Hamiltonian passed as a stack of one cost 211 µs through the array path
  against 66 through the scalar one.

  The thresholds are measured, and differ between the expansions because the
  two-flavor one does much less work per element and so amortises sooner:
  **eleven** elements for three flavors, **seven** for two. An earlier draft
  used sixteen for both, from a measurement that prepared the nested lists
  outside the timed region; the real path has to convert them, which is most
  of the per-element cost.

  That conversion is worth making explicitly: the scalar path is quicker on
  Python complex numbers than on NumPy scalars, by more than `.tolist()`
  costs.

## [1.5.0] - 2026-07-31

A deep pass over the whole computation, looking for work that was being done
and thrown away.  Nothing about the results changes: all 42 figures generated
by `run_testsuite.py` remain byte-for-byte identical, and the probabilities
agree with 1.4.0 to 1.6e-13 across every code path.

Measured best-of-seven, interleaved against 1.4.0:

| Benchmark | 1.4.0 | 1.5.0 | Speedup |
|---|---|---|---|
| Energy scan, 2000 points | 4.91 ms | 1.42 ms | **3.5x** |
| Two-flavor scan, 2000 points | 0.23 ms | 0.05 ms | **4.6x** |
| Oscillogram, 100 x 100 | 6.61 ms | 3.65 ms | **1.8x** |
| Baseline scan, 2000 points | 1.12 ms | 0.70 ms | **1.6x** |
| Three-flavor probability, scalar | 40.4 µs | 12.5 µs | **3.2x** |
| Two-flavor probability, scalar | 4.95 µs | 1.51 µs | **3.3x** |

Cumulatively, against the 1.1.0 audit release --- which had no batched
interface, so the comparison is a Python loop against a single call --- a
2000-point energy scan with the Hamiltonians built goes from 93 ms to 1.9 ms,
a factor of 48.  The same loop written today, still one point at a time, takes
39 ms: the scalar path accounts for 2.4x of that and vectorising for the rest.

### Changed

- The batched star product no longer contracts the dense 8x8x8 `d` tensor
  through `np.einsum`, which with no path plan walks the whole table for every
  element: for a 2000-point energy scan that single line was 70% of the total.
  `_star_all`, the sparse expansion the scalar path already used, vectorises
  unchanged and is 14x quicker there. (`optimize=True` was tried first: 1.1x.)

- `probabilities_2nu` and `probabilities_3nu` no longer build the evolution
  operator, square it, and then transpose and reshape the result into the
  order they return. They form the entries and square them as they go.

  For two flavors the saving is larger still. The coefficients of a Hermitian
  2x2 Hamiltonian are real, so |U_ee|² = u₀²+u₃² and |U_μe|² = |U_eμ|² =
  u₁²+u₂² — two distinct numbers, the second the complement of the first by
  unitarity. Neither the operator nor the coefficients are needed.

- The batched coefficients are laid out with the component index **first**, so
  that each h_k is contiguous. Every downstream step works one component at a
  time, and reading those from a strided `(..., 8)` view costs about a third
  more. On an oscillogram: forming u_k, 746 → 480 µs; the nine entries and
  their moduli, 2175 → 1598 µs. The two-flavor module gets the same layout, so
  the two mirror each other as their documentation claims.

- The coefficients u₀ and u_k are returned separately rather than concatenated
  into a `(..., 9)` array that every caller took apart again.

- Around the latent roots: √(|h|²) is taken once rather than three times; the
  arc-cosine argument is a division rather than a power of −1.5 (26.9 → 3.8
  µs); and the degeneracy test uses two calls to `np.minimum` instead of
  stacking three gap arrays to reduce along the axis it just created (106 → 13
  µs).

- The scalar exponentials use `cmath.rect(1, t)` rather than
  `cmath.exp(1j*t)` — the same value by construction, a third quicker.

### Added

- Tests for two-flavor broadcasting: one Hamiltonian against many baselines,
  and a grid. Neither had coverage, and the first is the shape that the layout
  change initially broke. Also: stacks with more than one leading axis, and
  shapes that do not broadcast, which must raise rather than quietly produce
  the wrong thing.

Note that the *ratio* between a vectorised scan and the equivalent Python loop
has narrowed across releases, from about 80x to about 30x for a baseline scan,
even though the vectorised scan itself keeps getting quicker: the loop got
quicker too.  `docs/source/methodology.rst` now carries measured ratios rather
than ones carried forward from an earlier release.

### Not changed, having been measured

Two plausible-looking optimisations were tested and rejected:

- `|z|²` as `z.real² + z.imag²` or `(z·z̄).real` beats `np.abs(z)**2` below
  about 10⁵ elements but is up to 1.6x **worse** above it, where the extra
  temporary stops fitting in cache. `np.abs(z)**2` is kept.

- Computing the batched complex exponential as `cos + i·sin` is slower than
  `np.exp` at every size that matters here.

## [1.4.0] - 2026-07-31

Type annotations throughout, and a pass over the arithmetic that was spending
more time in NumPy's dispatch machinery than in the calculation.

Nothing about the results changes: all 42 figures generated by
`run_testsuite.py` are byte-for-byte identical to those from 1.3.0.

### Added

- Type annotations on every parameter and return value, including the private
  helpers, in the style used by Magnus.  Sphinx renders them into the
  signature, so the API page now states the polymorphism of the core routines
  rather than leaving it to the prose.

- `tests/test_annotations.py`, which keeps the annotations and the docstrings
  in step: everything is annotated, every annotation resolves to a real type,
  every public parameter is documented, and an annotation that admits an array
  paired with a docstring that promises a scalar fails the suite --- in either
  direction.  That bidirectional check immediately caught a fix that had landed
  on the wrong function; a one-directional audit had passed it.

- `oscprob3nu._star_all`, the sparse expansion of
  :math:`(h \star h)_i = d_{ijk} h_j h_k` written out.  The dense
  8x8x8 table is kept for the batched path, where the array machinery pays for
  itself, but for eight numbers it spends several microseconds of dispatch on a
  few dozen multiplications.  The explicit form is six times quicker and agrees
  to 1.1e-16; a test checks it against `tensor_d` term by term.

### Changed

- The scalar paths no longer dispatch NumPy for single numbers: `np.real` and
  `np.imag` on one complex number, `np.arccos`, `np.clip` and `np.sqrt` on one
  float, all give way to attribute access and the `math` module.

- The sum over the three latent roots is rewritten using the fact that two of
  its factors do not depend on k,

      u_k = i [ (sum_m w_m psi_m) h_k - (sum_m w_m) (h*h)_k ] ,

  which forms those combinations once instead of inside the loop over the eight
  k, and removes an `(N, 3, 8)` intermediate --- eight times the size of the
  result --- from the batched path.

- The star product was computed twice per call, once to form the invariant
  ⟨h⟩ and again inside the expansion.  It is now computed once and passed on.

- `abs(z)**2` becomes `_abs2(z)`, which skips the square root that the square
  immediately undoes.

- The closed-form three-flavor vacuum Hamiltonian rebuilt the CP phase fifteen
  times across its nine entries, and recomputed two products in five places.
  Hoisting them shortens the expressions materially and takes the routine from
  75.0 to 32.1 microseconds.

Measured best-of-seven, interleaved against 1.3.0:

| Benchmark | 1.3.0 | 1.4.0 | Speedup |
|---|---|---|---|
| Three-flavor probability, scalar | 46.4 µs | 14.5 µs | **3.2x** |
| Two-flavor probability, scalar | 5.7 µs | 1.7 µs | **3.3x** |
| Energy scan, 2000 points | 5.75 ms | 4.99 ms | 1.15x |
| Oscillogram, 100 x 100 | 8.05 ms | 6.70 ms | 1.20x |
| Baseline scan, 2000 points | 1.33 ms | 1.29 ms | 1.03x |
| Vacuum Hamiltonian | 75.0 µs | 32.1 µs | 1.9x |

The scalar paths gain most, because they were the ones paying dispatch
overhead on every operation; the batched paths were already NumPy-bound, and
gain only where an intermediate array was removed.  The baseline scan needed
care in the other direction: for a single Hamiltonian the old expression was a
3-by-8 matrix product that BLAS does better than three broadcasts, and keeping
it as such is what took that case from a 6% regression back to parity.

## [1.3.0] - 2026-07-31

Completes the vectorisation started in 1.2.0.  The expansions could already
take a stack of Hamiltonians, but the routines that *build* those Hamiltonians
could not, so an energy scan still had to loop before the fast path ever saw
it.  Now the whole scan is two calls.

### Added

- `hamiltonian_2nu_matter`, `hamiltonian_2nu_nsi`, `hamiltonian_2nu_liv` and
  their three-flavor counterparts accept an array of energies and return one
  Hamiltonian per energy, stacked along a leading axis --- which is the shape
  `probabilities_3nu` consumes.  A scalar energy still returns a single
  `n`-by-`n` matrix.  The matter potential may also be an array, for a scan
  across a density profile alongside the energy.

  The results are bit-for-bit identical to the previous loop, not merely close:
  the arithmetic per element is unchanged.  A 200-point three-flavor energy
  scan in matter, Hamiltonians included, goes from 24.7 ms to 0.71 ms.

- 9 tests for the batched builders (`tests/test_vectorized_hamiltonians.py`),
  including the end-to-end scan that the plotting modules perform, an
  energy-by-baseline grid, and an array-valued matter potential.

### Changed

- Each builder is now a single expression --- the vacuum term divided by the
  energy, plus a constant matrix scaled by the potential --- rather than a
  sequence of entry-by-entry additions.  Indexing the energy with two trailing
  axes lets the same expression serve a scalar energy and an array of them.
  The matter potential multiplies a named module-level projector, and the NSI
  and LIV terms are written as the matrices they are.

- The four plotting modules use the batched interface: eighteen list
  comprehensions become single calls, and the nine `[x[k] for x in prob]`
  column extractions become `prob[:, k]`.  All 42 figures were generated both
  ways and are byte-for-byte identical.

  This makes the *computation* in those modules about 35x faster but the suite
  only about 6% faster --- 24.9 s to 23.3 s --- because matplotlib dominates
  the wall time. The case for the change is consistency with the library's own
  recommended usage, not speed.

- The default branch is renamed from `master` to `main`, and the documentation
  links that named it are updated.  GitHub redirects the old URLs, but a
  redirect is not a reason to keep serving stale links from our own docs.

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
