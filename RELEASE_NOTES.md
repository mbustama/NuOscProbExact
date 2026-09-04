This is everything since the published **v1.11.0**, which is versions 1.12.0, 1.13.0 and
1.13.1 together — 1.12.0 was declared but never published, so upgrading from 1.11.0 brings
all three.  What follows describes 1.12.0 and 1.13.0; 1.13.1 adds a batch of chords to the
slab routines, the electron-fraction values below, a single call for antineutrino
oscillograms, and a mean nucleon mass that follows `Y_e`.  The changelog has each in full.

It is mostly about the **Earth**. It had no compiled path at all, took one energy at a time,
and asked for a slab count where it should have asked for an accuracy; all three are fixed
here. Two changes are not about the Earth and are the ones most likely to affect you
anyway: **`numba` is now a dependency**, so a plain install gets the compiled path — and
brings numba's ceiling on `numpy` with it — and the **four-flavor latent roots are now
computed in double-double arithmetic**, which takes their worst error to under a fifth of a
`float64` ulp. All three are covered below.

The full history is in
[CHANGELOG.md](https://github.com/mbustama/NuOscProbExact/blob/main/CHANGELOG.md), which
records how every figure here was measured.

## Upgrading from 1.11.0: Earth results move in the last digits

**Results are not bit-identical to 1.11.0 on any Earth call.** This is stated here rather
than left to be discovered.

The chord palindrome below requires that a chord's slab widths be *exactly* symmetric. They
were not: each segment was cut by its own `linspace` and the two halves rounded differently,
by about 1e-12 km on a 100 km slab. Averaging each width with its mirror makes them bitwise
equal, and shifts the probabilities. Measured over 4640 values across four angles, four
subdivisions and five energies, the shift is **6e-15 at two flavors, 7e-14 at three, and
4e-10 at four**.

The four-flavor figure is much the largest, and it is the conditioning of the quartic latent
roots rather than the widths: that same conditioning already separated the two backends by
4.0e-10 on code the width change never touched. It is also the thing this release goes on to
fix — see the double-double roots below.

For scale, the discretisation error at the default subdivision is **1e-4 to 1e-5** — five
orders of magnitude larger than even the four-flavor shift, and eleven larger than the
three-flavor one. A scalar energy still returns a tuple, bit-for-bit what it returned before.

**Four-flavor results also move, everywhere, not only in the Earth.** The latent roots now
come from a different route, so any `probabilities_4nu` may change in its last digits.
Measured against what 1.11.0 computed, on 3+1 spectra at 0.5, 1 and 10 GeV and
Δm²₄₁ of 0.1 and 1 eV²: **bit-identical in two of those six, and at most 1.2e-10** in the
rest. The new values are the accurate ones, by a factor of about seven.

## What is new since 1.11.0

**The compiled backend is no longer optional.** `numba` moved from the `fast` extra into the
base dependencies, so a plain `pip install nuoscprobexact` gets the compiled path. It had
been an extra, which meant the 12× on an Earth crossing reached whoever read far enough down
the README.

The cost, stated because it will otherwise be discovered by surprise: `numba` declares its
own ceiling on `numpy` — 0.66.0 requires `numpy<2.5` — and **this library inherits it**.
Installing into an environment holding numpy 2.5.1 downgrades numpy to 2.4.6. That is not a
corner case: across numba 0.60 to 0.66, the ceiling excluded the then-current numpy either on
the day numba shipped or within three months, and numpy releases a minor roughly every six.
No spelling on our side avoids it, because the pin is numba's.

If the newest numpy matters more than the speed, install with `--no-deps` and add numpy
yourself, or keep numba installed and set `fastkernels.USE_NUMBA = False`. The NumPy path is
still there, still tested on every push, and still agrees to round-off.

**Four-flavor roots in double-double arithmetic.** Worst relative error over nine reference
Hamiltonians goes from 3.9e-16 to **3.6e-17** — under a fifth of a `float64` ulp — for
roughly a fifth more time on a full `probabilities_4nu` through the compiled kernel. The
error is exact; the cost is deliberately a range, since timed in alternated pairs it lands
anywhere from 1.18× to 1.25× depending on method, with individual pairs from 0.76× to 1.72×.

The problem was never the quartic solver. The three invariants `I₂, I₃, I₄` compress a 4×4
matrix into three numbers, and the amplification from those coefficients to the roots
measures **2.3e9** — so a coefficient rounded at 1e-16 becomes a root wrong at 1e-7, which is
what the closed form alone scores and what no better root-finder can beat. 1.12.0 got around
it by refining against the matrix; this removes it, by carrying each coefficient as a pair of
`float64` limbs, some 32 digits, so the same amplification acts on 1e-32 and lands at 1e-23.
The roots are then limited by rounding the answer to `float64` rather than by the algebra.

`oscprob4nu.ROOT_STRATEGY` chooses, defaulting to `'double-double'`, with `'eigensolver'` for
the LAPACK route; `POLISH_ROOTS` keeps its meaning and is read only on the latter. Both routes
run on both backends and agree to an ulp, so a result never depends on whether a compiler was
present. If you call `fastkernels`' four-flavor kernels directly, their `polish: bool`
argument is now `strategy: int`.

Two limits worth stating. On the NumPy fallback the cost ratio is of order 1.5–2× rather
than 1.2×. And this buys **roots**, not probabilities everywhere: rebuilding `U₄` takes second
differences of `exp(-iψL)`, so a root error enters twice and the probability error grows as
`(ψ_max L)²`. At Δm²₄₁ = 1000 eV² over 1300 km the accumulated phase is 2.5 million radians
and every route lands at 2e-4 together — that is `float64`, not the roots. The physical range
is where the tight figures are: **1.3e-12** at 0.1 eV².

**The Earth, compiled.** The backend offered probability kernels only, so composing evolution
operators across slabs — which is what layered matter and the Earth do — ran the NumPy path
however the `fast` extra was installed. Measured before the fix, `probabilities_3nu_earth`
took 2.11 ms with `numba` and 2.09 ms without, and entered a kernel zero times. The new
`earth_chords_{2,3,4}nu_kernel` build each slab's Hamiltonian in registers rather than
materialising a stack, and one Earth crossing is now **~12×, ~12× and ~9×** quicker at two,
three and four flavors.

**Whole Earth scans in one call.** `probabilities_{2,3,4}nu_earth` take an **array of
energies**; they took a scalar, so a scan was N separate calls, each rebuilding the matter
potentials and re-validating slabs the library had just built itself. A hundred-energy scan
is **6.7×, 4.8× and 3.5×** quicker than the loop it replaces; against the uncompiled
per-energy loop, a thousand-energy scan is **121×, 49× and 17×**. The three
`_between_locations` wrappers inherit it.

**Oscillograms, from an array of zenith angles.** Energies and angles broadcast against each
other, so a grid is asked for as
`probabilities_3nu_earth(h, energies[None, :], costhz[:, None])` — the idiom
`oscprob3nu.probabilities_3nu` already uses for Hamiltonians against baselines. Measured
**8×** against the loop over grid points, at 4.8 µs per point for a 60 × 200 grid. The angles
are still evaluated one at a time, and that is not a shortcut: the chord geometry is what the
angle changes, so two angles share neither their slab widths nor even how many slabs they
have.

**The chord palindrome.** A neutrino crossing a spherically symmetric Earth meets every radius
on the way in and again on the way out, so slab *j* and slab *n−1−j* carry the same
Hamiltonian, the same width, and therefore the same operator. Half of them were being computed
twice. Composing from the centre outwards computes each once, worth **~1.4×, ~1.5× and ~1.8×**
on top of the figures above. Nothing about it is specific to PREM — a symmetric castle wall is
composed at the same discount — and `fastkernels.USE_PALINDROME = False` asks for the plain
left-to-right product instead.

**A tolerance instead of a slab count.** `rtol` and `atol` on every Earth probability routine,
plus `earth.slabs_for_tolerance` to size the subdivision once for a whole scan. This matters
because the discretisation error is strongly energy-dependent — at the default eight sub-slabs
per segment it spans more than an order of magnitude between 3 and 40 GeV — so a fixed
`n_slabs_per_segment` does not give a fixed accuracy. Both tolerances unset leaves the result
bit-identical to before.

**Any profile, not just the Earth.** `slabs.probabilities_{2,3,4}nu_profile` take the profile
as a callable rather than knowing about PREM, so the same tolerance machinery serves a
hand-built density profile, a castle wall or a solar model. Equal slabs sampled at their
midpoints, which is second order and so refines by the law the search assumes — verified at a
ratio of 4.00. Where a profile is *discontinuous* this is the wrong tool, since no refinement
recovers a jump that straddles a slab; split the trajectory at the discontinuities instead,
which is what `earth` does with the PREM shells. The tolerances are deliberately absent from
`oscprob{2,3,4}nu`, which compute an exact closed form with no discretisation to refine.

**Long scans stay inside memory.** The Hamiltonian stack is proportional to the scan, so a
long one is chunked. The chunking turned out to matter for a second reason: the batched kernel
is **memory-bound, not compute-bound** — 7.9 µs per probability with an 8 MB working set
against 16.4 µs with a 540 MB one, a factor of two paid for traffic alone.
`earth.MAX_CHUNK_BYTES` therefore defaults to the detected last-level cache rather than to a
fixed size, worth 1.3× to 1.8× on long scans.

**A twentieth notebook.** `20_arbitrary_hamiltonian.ipynb` carries a Hamiltonian of your own
through three varying profiles — an invented body, an exponential solar-like one, and the
Earth — using a long-range force from a gauged L<sub>e</sub> − L<sub>μ</sub> symmetry as the
worked case. It also documents the one thing about `earth` that the API does not reveal:
`probabilities_3nu_earth` builds *H* = *H*<sub>vac</sub>/*E* + *V*<sub>CC</sub>*P*<sub>ee</sub>
itself, so a Hamiltonian not of that form cannot go through it at all. The way through is
`earth_slabs` → `matter_potential` → your own construction → `probabilities_3nu_slabs`, and
with the extra term switched off that path reproduces `probabilities_3nu_earth` to zero
difference, because it is the same arithmetic.

**Corroboration against five other codes.** `tests/prem_scan.py` and the frozen
`tests/prem_speed_accuracy.json` carry the Earth speed–accuracy planes against nuSQuIDS,
nuCraft, GLoBES, Prob3++ and NuFast, and `tests/external_drivers/` ships the C and C++ drivers
themselves, with the raw output of each and a README recording what every one of them had to
be told. `gen_prem_header.py` emits this library's PREM as a C header rather than anyone
transcribing forty coefficients, reproducing `earth.density_prem` to 5e-14 over 6372 radii.

**A harness for changes that must not move a number.** `tests/bit_capture.py` and
`tests/bit_compare.py` run in pairs and report in ulps. No golden file is committed,
deliberately: libm is not bit-identical across platforms. Its first real use was to **decline**
a change.

## Also fixed

- **The Earth performance figures had drifted, and nothing was guarding them.** `13.9x, 9.6x
  and 6.6x` was stated in three documents, agreed with itself in all three, and was wrong in
  two — the palindrome had made the compiled side quicker, so three and four flavors were being
  understated. Re-measured, and now guarded by tests that fail if one document disagrees.
- **The documented notebook count had rotted**, and the README's table of notebooks stopped at
  18. Corrected, completed, and guarded by a test that derives the count from `notebooks/`.
- **The convergence order of the slab product was being quoted as first.** It is second, which
  is what the midpoint sampling gives.

## The Earth's electron fraction: the values ship now, the default changes later

**The default `electron_fraction` will change from 0.5 to per-layer values.** The values and
the machinery to use them are in 1.13.1; only the default is deferred, and it is announced
here so that it is not discovered later.

PREM is a density model and carries no composition, so `Y_e` has to be supplied. Every `earth`
routine takes one half throughout, which is exactly isoscalar matter — and no layer of the
Earth is. The values that belong to the material of each layer are now in the code, and
`earth.electron_fraction_prem` returns them by radius: **0.4656** in the core (iron),
**0.4957** in the mantle (peridotite), **0.4952** in the crust (granitic), and **0.5551** in
the ocean (seawater, above one half because hydrogen has `Z/A = 1`). Every routine's
`electron_fraction` argument now accepts a callable of radius, so passing
`earth.electron_fraction_prem` itself is all it takes -- and that is the form that works
alongside `rtol` and `atol`, which choose the slab count for you.

**Nothing changes yet.** One half remains the default, so every result in this release is what
it was. The change is deferred because it is not small: over 0.3–30 GeV the median change in
`P(νμ→νe)` is **67% through the diameter, 3% through the mantle alone, and 0.6% on a shallow
chord**, and flipping it also means re-measuring the comparisons against external codes, which
are handed a matched potential.

**What to do about it.** If you want the layered values now, build the array as above. If you
want today's numbers to survive the change, pass `electron_fraction=0.5` explicitly; that will
keep working and will mean the same thing before and after.

## The paper now travels with the code

`resources/paper/` carries the source of the Computer Physics Communications article that
documents this library: `main.tex`, all fifteen figures, `elsarticle.cls` and
`elsarticle-num.bst` so it compiles without the Elsevier bundle, and `make_versions.py`,
which derives the marked-up diff against the published version. The figures in it are drawn
by `notebooks/10_paper_figures.ipynb` from data frozen in `tests/`, so a claim in the text
can be traced to the measurement behind it without leaving the repository.

**A clone gets it; `pip` does not.** Checked rather than assumed, on this release: the wheel
holds 15 entries and the source distribution 24, all modules or metadata, and neither
matches `resources`, `paper`, `.tex` or `.pdf`. Nothing about what you install has changed.

`tests/external_drivers/` holds the C and C++ drivers for the codes that cannot be called
from Python — GLoBES, Prob3++, NuFast — with the raw output of each and a README recording
every convention that had to be matched first: three different roundings of the atomic mass
unit among them, and a zenith angle in radians where the documentation says degrees.

The paper's arbitrary-Hamiltonian example is now a **long-range interaction** — a gauged
$L_e - L_\mu$ under which the Earth's electrons source a Yukawa potential, non-local and so
not expressible as a density profile — with the figure for it drawn stand-alone in notebook
10. It replaces a template that called an undefined `hamiltonian_mymodel`, which could
assert the library's generality but not show it.

## A note on the numbers

Every speedup above was measured on one laptop on one day, and figures like these move by tens
of per cent between machines. They are quoted because a claim with no number attached cannot
be checked, not because the third digit means anything. The palindrome is the one figure a
reader can switch off and re-measure directly, via `fastkernels.USE_PALINDROME`.

## Install

```shell
pip install nuoscprobexact
```

Python 3.9 or newer, tested on 3.9 through 3.13. `numpy` and `numba` are both installed, so
the compiled backend is what a plain install gets — see the note above on numba's numpy
ceiling. `pip install "nuoscprobexact[fast]"` still works and is now a no-op.

Twenty worked notebooks, and a regression suite of 797 tests at 100% coverage
against a 98% floor, cross-checked against [nuSQuIDS](https://github.com/arguelles/nuSQuIDS)
and against the Zaglauer–Schwarzer closed form for the matter spectrum.

## Links

- [Documentation](https://mbustama.github.io/NuOscProbExact/)
- [Quickstart](https://mbustama.github.io/NuOscProbExact/quickstart.html)
- [What it can compute, with code](https://mbustama.github.io/NuOscProbExact/recipes.html)
- [The paper](https://arxiv.org/abs/1904.12391)
