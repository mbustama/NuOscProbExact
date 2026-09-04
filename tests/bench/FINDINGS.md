# What the adversarial check broke

`ADVERSARIAL.md` named six guarantees and asked for each to be falsified rather
than reviewed. Five broke. This file records what broke, where, and how each
claim was checked, because a finding that cannot be re-checked is an opinion.

The check itself read source and ran the test suite. It started no measurement,
compiled nothing, and left the repository untouched; the one throwaway adapter
it wrote lived in a scratch copy of this directory. Every claim below was
verified against the pinned sources under `.bench-build/src`, not against the
comments describing them. Where a claim rests on reasoning rather than on
source, it says so in the line marked *inference*.

## The short list, worst first

| # | what | where | state |
|---|---|---|---|
| 1 | `OurPREM` leaves three base-class fields uninitialised | `adapters/nufast_earth.cpp:34` | source-confirmed |
| 2 | THROUGHPUT measures a cached copy for NuFast-Earth | `bench.hpp` `detail::throughput` | source-confirmed |
| 3 | All twelve invariants pass a deliberately unfair adapter | `../test_bench_pipeline.py` | reproduced |
| 4 | Prob3++ does not use the hbar c the conversion assumes | `conversions.py:137` | source-confirmed |
| 5 | `.bench-build/bin/` predates the contract it was built against | timestamps | source-confirmed |
| 6 | NuFast-Earth answers a constant-density request with a chord | `adapters/nufast_earth.cpp:106` | source-confirmed |
| 7 | The harness has no field for `n_layers` or `shells_total` | `bench.hpp` `main` | source-confirmed |
| 8 | nuSQuIDS is read out through its interpolating call | `adapters/nusquids.py:143` | source-confirmed |
| 9 | `OurPREM` hard-codes the electron fraction | `adapters/nufast_earth.cpp:42` | source-confirmed |
| 10 | The reference builder's inputs are incomplete | `manifest.json`, `conversions.py` | inference |

---

## 1. `OurPREM` never initialises the fields the engine reads

`NuFast::Earth_Density` (`.bench-build/src/nufast-earth/include/Earth.h:8`)
carries four fields a subclass is expected to fill:

```
int                 n_discontinuities;
std::vector<double> discontinuities;
double              r_E = 6371;
bool                constant_shells;
```

Only `r_E` has a default member initialiser. Every upstream subclass sets the
other three in its constructor --- `PREM_NUniformLayer`, `PREM_Full`,
`PREM_Four`, `PREM_Prob3`, `Constant` and both `NDiscontinuityLayer` variants,
at `Earth.cpp:44`, `:90`, `:108`, `:121`, `:132`, `:145`, `:156`, `:163`,
`:182`, `:191`, `:206`, `:210`, `:234` and `:284`.

`OurPREM` sets none of them. The word `discontinuities` appears exactly once in
`adapters/nufast_earth.cpp`, in a comment on line 15. `Set_Earth`
(`NuFastEarth.cpp`) stores the pointer and validates nothing. So
`new OurPREM(...)` hands the engine an object whose `n_discontinuities` and
`constant_shells` are indeterminate and whose `discontinuities` vector is
empty.

What the engine then does with them:

* `Geometry.cpp:47` loops `for (i = n_discontinuities - 1; i >= 0; i--)` and
  dereferences `discontinuities[i]` --- an out-of-bounds read on an empty
  vector, under a garbage bound.
* If `constant_shells` reads true, `Calculate_Eigens` evaluates
  `rhoYe(discontinuities[j] - 1e-8)` on the same empty vector, and
  `Calculate_Internal_Amplitudes` uses the result to index `eigens_constant[i]`.
* If it reads false, the engine takes the `eigens_varying` branch: one
  eigendecomposition per (energy, cosz, layer), with **no caching across zenith
  angles**. That caching is the whole subject of objection Earth-3, and
  `OSC/100x100` exists to expose it.
* With `n_discontinuities` at zero, `Mean_Density` collapses each half-chord to
  a single segment averaged on a 10 km step, discarding the radial structure
  entirely.

*Inference*: this is the most likely source of the constant
0.05-per-probability offset between `bench_nufast_earth` and `prob3`/`globes`
that `README.md` records as a loose thread. A defect that survives refinement
and ignores the precision knob has the right shape. Confirming it needs one
untimed checksum run before and after the fix.

**Fix.** `OurPREM` must publish its 1024 sub-shell boundaries as
`discontinuities`, ascending, with the surface as the last entry, set
`n_discontinuities` to their count, and set `constant_shells = true` --- which
is also the branch that gives NuFast-Earth its zenith caching. Shell index `j`
means the region below `discontinuities[j]`, which is the convention
`Calculate_Eigens` and `Calculate_Internal_Amplitudes` agree on.

## 2. THROUGHPUT measures a memcpy for NuFast-Earth

`bench::detail::throughput` calls `driver::configure(p.dcp)` exactly once,
before any clock, and then calls only `driver::evaluate()`.

NuFast-Earth's `evaluate()` is `g_engine->Get_Probabilities()`, which is
`if (not probabilities_calculated) Calculate_Probabilities(); return probabilities;`.
`probabilities_calculated` is cleared only by a setter, and the only setter the
adapter calls after setup is `Set_delta`, from `configure`. So the first
repetition computes and every repetition after it returns the cached
`vector<vector<Matrix3r>>` by value.

The autorange loop makes this worse rather than louder: it grows `reps` until a
block reaches 0.05 s, so it will keep multiplying copies until the copying alone
fills the block, and the reported microseconds per point falls accordingly.

NuFast-Earth is the only adapter affected. NuFast-LBL calls a free function,
Prob3++ calls `SetMNS`/`propagate` per point, GLoBES calls
`glb_probability_matrix` per energy, nuSQuIDS re-runs `Set_initial_state` and
`EvolveState`, nuCraft re-runs `CalcWeights`, and this library recomputes. It is
also the code whose author raised the objections, which is the worst place for a
fabricated number.

AMORTIZED is unaffected: it calls `configure` once per scan step, and
`Set_delta` clears the flag.

**Fix.** Either `throughput` calls `configure` before each `evaluate`, which
changes what the protocol means for every code, or the NuFast-Earth adapter
invalidates in `evaluate()`. The second is narrower and keeps the protocol
definition the author's.

## 3. The invariants pass an adapter that solves nothing

A throwaway adapter, written to fail as many invariants as possible, passed all
twelve. It:

* computes the whole answer in `setup()` and returns a stored scalar from
  `evaluate()`, so nothing timed does any work;
* ignores `dcp` in `configure()`, so the AMORTIZED scan moves nothing;
* claims `batches_energy` and `batches_zenith` with a fabricated
  `batch_symbol` while doing neither;
* reaches the matter constant as `float('1.5' '2588e-4')`. `_code_only()`
  strips string literals before the literal check runs, so the grep cannot see
  it.

Naming itself `NuOscProbExact` made `test_capabilities_claims_match_the_registry`
engage rather than skip, and it passed that too.

The invariants guard provenance --- who owns a clock, where a constant came
from, whether a version is pinned --- and none of them observes behaviour.
Three that would:

* `evaluate()`'s result must change when `configure()` changes `dcp`;
* `evaluate()`'s cost must scale with `problem.points()`;
* two knob values from the declared domain must give different answers.

## 4. Prob3++ does not use hbar c = 1.97327e-7 eV m

`conversions.hbar_c_cosine_scale` takes `their_hbar_c_ev_m=1.97327e-7` as a
default argument, no caller overrides it, and the generated `conversions.h`
prints "All three compiled Earth codes hard-code hbar c = 1.97327e-7 eV m".
Read from the pinned sources:

| code | site | value or implied hbar c |
|---|---|---|
| NuFast-LBL | `NuFast_LBL.cpp:11` | `1e-9 / 1.97327e-7 * 1e3 / 4` |
| NuFast-Earth | `src/Oscillation.cpp:14` | `1e-9 / 1.97327e-7 * 1e3 / 2` |
| GLoBES | `globes/globes.h:83` | `GLB_EV_TO_KM_FACTOR 1.97327e-10` |
| nuCraft | `NuCraft.py:550` | `1.266933` |
| **Prob3++** | **`mosc.c:203`, `mosc.c:489`** | **`2.534`, implying 1.9731650e-7** |

Prob3++ rounds to four significant figures. Its kilometre is 5.068000e9 eV^-1
against 5.067730e9 --- a relative difference of **5.3e-5**, three orders of
magnitude above the 4.2e-8 the manifest describes as the deliberate cap on the
old six-way figure. `prob3.cpp` applies `OUR_COSZ_HBARC_SCALE`, built from
NuFast's value, so 5.3e-5 of baseline error survives. On `CHORD/12x1`
(L about 11 468 km, 3 to 40 GeV) that is roughly 6e-4 radians of phase, which
puts a floor near 1e-4 on Prob3++'s accuracy that is entirely convention.

The deeper problem is that `conversions.py` extracts every matter constant from
pinned source and extracts hbar c for no code at all. A per-code reference has
to absorb each code's own hbar c; there is nothing to read it from.

**Fix.** Add hbar c to `SITES`, one entry per code, and derive the scale per
code rather than from one default argument.

## 5. The built binaries predate the contract

Every file in `.bench-build/bin/` is stamped 13:28. `bench.hpp` is 13:32 ---
when `probabilities()` and the untimed accuracy protocol were added ---
`conversions.h` is 13:39, and `nufast_earth.cpp` and `nufast_lbl.cpp` are 14:03.

The `_accuracy` binaries exist and will run. Invoking one with
`--protocol accuracy` reaches a harness that had no accuracy branch, so it falls
through to a timing protocol and prints a time. That is exactly the confusion
the untimed accuracy path was written to make impossible, available now from a
file that looks built.

Separately, `speed_accuracy.json` and `timing_other_codes.json` still carry no
`generated_by`, and `nufast_scan.json`'s names a citation rather than a path.
OWN-2 is unchanged.

## 6. NuFast-Earth answers a constant-density request with a chord

`adapters/nufast_earth.cpp` always calls `Set_Earth` and `Set_Spectra`, and when
`p.costhz` is empty it substitutes `g_z = {-0.9}`. A run on `CONST/60E`
therefore propagates through the Earth at one zenith angle and stamps the
artifact `"grid": "CONST/60E"`.

`Set_E_Spectra` and `Set_Trajectory` --- the single-trajectory mode added in
v1.1.0, which is the whole subject of objection Earth-4 --- are never called.
Earth-4 cannot be answered with what is built, and the artifact that would be
written looks correct.

## 7. Nothing in the artifact can carry the shell count

`test_shell_count_semantics_recorded` needs `knob.n_layers` and `shells_total`.
`bench.hpp`'s JSON emits one knob field, `"knob": {"<knob_name>": <int>}`, and
NuFast-Earth's `knob_name` is `eigenvalue_precision`. There is no field for
either quantity, and `n_per_layer` is hard-coded at 256, so objection Earth-1
cannot be answered as the table promises.

This is the same gap `README.md` records as "`n_layers` cannot actually be
swept", seen from the artifact side rather than the `Problem` side.

## 8. nuSQuIDS is read out through its interpolating call

`adapters/nusquids.py:143` uses `EvalFlavor(1, float(ee))`. The pinned wheel
also exposes `EvalFlavorAtNode(flavor, ie)`, documented "Flavor content at the
node (no interpolation)". The grid energies are the solver's own nodes, built
from the same vector handed to the multiple-energy constructor.

So the adapter pays an interpolation inside the timed region and carries
interpolation error into the accuracy number. Both cost nuSQuIDS, on both axes,
for nothing.

## 9. `OurPREM` hard-codes the electron fraction

`adapters/nufast_earth.cpp:42` is
`rho_[L * n_ + i] = our_prem_rho(r) * 0.5 * scale;`, and `setup()` never reads
`p.ye`. It is the only adapter that ignores it --- `nufast_lbl.cpp:67`,
`prob3.cpp:100` and `globes.cpp:124` all use it. The two agree today only
because `Problem::ye` defaults to 0.5, and a bare `0.5` is invisible to
`test_no_adapter_types_a_physical_constant`.

## 10. The reference builder's inputs are incomplete

*Inference throughout; the builder does not exist yet, so this is about what it
would be built from.*

`manifest.json` names two things absorbed into each code's reference: hbar c and
the code's own matter normalisation. A third differs per code and is not
mentioned --- the radial discretisation of PREM:

* the three compiled Earth adapters cut four major layers into n midpoint
  shells;
* `nusquids.py` builds a piecewise-linear table with 200 nodes per PREM shell;
* this library integrates the full nine-boundary PREM.

If a reference absorbs only the first two, each residual mixes algorithm error
with a profile difference, and at 1024 shells the profile term is not small.
The manifest has to say whether a code's reference is built on the continuous
PREM --- in which case the residual includes its discretisation, which is a
defensible choice --- or on the discretised profile it was actually handed, in
which case it measures the algorithm alone. Either is fine; leaving it unstated
is how a figure ends up meaning neither.

`reference_conventions.never_absorbed` also has nothing to build from yet:
`conversions.mass_defect('nuCraft')` extracts only the charged-current entry
(`15.256e-5`), and the sterile entry (`7.6525e-5`, giving the 0.50161 ratio the
manifest calls an error rather than a convention) is extracted nowhere.

---

## What did not break

Recorded because evidence of absence is worth as much as the findings, and
because a check that reports only failures cannot be told from one that invents
them.

**Attack 6, NuOscProbExact under-representation.** Re-measured by tracing kernel
entry with no clock read. Every row of `ADVERSARIAL.md`'s table holds:

| check | measured |
|---|---|
| `HAVE_NUMBA` / `USE_NUMBA` / `available()` | true / true / true |
| numba threads | 12 |
| `CONST/60E` enters | `probabilities_3nu_kernel` |
| `CHORD/12x1` enters | `earth_chords_3nu_kernel` |
| batching | energies and zenith together |
| `_PMM = 4` | correct: `probabilities_3nu` returns `(N, 9)` row-major, and each row triple sums to 1.000000000000 |
| knob branches | 0 gives 8 slabs, 8 gives 8 slabs, -3 to -15 give one slab and `rtol = 10**knob`; all run |

The gap `ADVERSARIAL.md` names is the real one: nothing asserts any of this
during a run.

**The two claims `README.md` says must not be re-derived wrongly.** Both hold.
`.bench-build/src/nufast-lbl/NuFast_LBL.cpp` carries `YerhoE2a = 1.52588e-4`,
the same rounding Prob3++ uses, and the generated
`OUR_PREM_MASS_DEFECT_NUFAST_LBL` and `OUR_PREM_MASS_DEFECT_PROB3` are equal to
every digit, as they must be. The committed adapter passes the whole energy
vector in one call.

**Everything else in the timed regions.** GLoBES recomputes its mixing matrix in
`configure`, which the protocol is defined to include; Prob3++ keeps `SetMNS`
inside `evaluate`'s loop, where the interface puts it; nuCraft rebuilds its
instance in `configure`, which is the cost its interface charges a scan; the
GLoBES chord decomposition is O(N^2) but sits in `setup`, outside every clock.
`P[1][1]`, `GetProb(2,2)`, `EvalFlavor(1)` and `_PMM = 4` all name
numu to numu, agreeing with `manifest.json`'s `scored_channel`.

---

# What was fixed, and what fixing it turned up

Applied and verified on 2026-08-19, after the check above. Every number below
came from the untimed `accuracy` protocol or a correctness probe; no timing run
has happened yet.

## Fixed

| # | fix | evidence |
|---|---|---|
| 1 | `OurPREM` now sets `discontinuities`, `n_discontinuities` and `constant_shells = true` | NuFast-Earth's deviation from this library fell from **max 0.63** to **max 1.6e-5**, and it now sits 1.4e-5 from GLoBES instead of 0.63 |
| 2 | `evaluate()` calls `Set_delta` before `Get_Probabilities()` | THROUGHPUT can no longer time a cached copy |
| 4 | `conversions.hbar_c` / `km_to_inv_ev` read each code's own hbar c from its pinned source; the single default argument is gone and a bare call now raises | five distinct values across the seven codes, Prob3++ the outlier at 5.3e-5 |
| 8 | nuSQuIDS read out with `EvalFlavorAtNode` | see the correction below |
| 9 | `OurPREM` takes `p.ye`, and the PREM layer radii come from `our_prem.h` rather than four typed literals | no adapter now types a boundary |
| -- | all four compiled adapters implement `probabilities()`, so `build.sh` links again and the untimed accuracy path exists | seven binaries built; `.bench-build/bin/` no longer predates the contract, closing 5 |
| -- | every absorbed rescaling removed from every adapter --- both mass defects and the cosine scale, in `nufast_lbl.cpp`, `nufast_earth.cpp`, `prob3.cpp`, `globes.cpp`, `nusquids.py`, `nucraft.py` | each code is now handed honest physical inputs, as `reference_conventions` requires |
| 3 | `../test_bench_behaviour.py` adds three behavioural invariants: the scan must move something, `evaluate()` must consume the grid it is normalised by, and a knob's endpoints must differ | each fails on the throwaway adapter that passed all twelve provenance invariants |

The Prob3++ prediction is worth recording because it was made before it was
measured. From `2.534` alone the expected floor was "near 1e-4"; measured
against this library on `CHORD/12x1`, Prob3++ sat **3.4e-4** away while GLoBES
--- same profile, same cosine, different hbar c only --- sat **1.4e-6**.

## A correction to finding 8

`EvalFlavorAtNode` and `EvalFlavor(flavor, E)` returned **identical** values on
`CHORD/12x1`, to every bit. So the interpolating overload was not costing
nuSQuIDS accuracy on that grid, and the finding as written overstated it. It
still costs an interpolated lookup per point inside the timed region, which is
why the change stands, but it is a speed fix and not an accuracy one.

## Two new findings, both surfaced by fixing the others

### 11. nuSQuIDS cannot run either Earth grid at its default tolerance

Removing the cosine rescaling --- which the per-code reference design requires
--- made nuSQuIDS fail outright:

```
SQUIDS::Evolve: Error in GSL ODE solver (failure)
```

Not a tiny-grid artifact: it fails at 2, 3, 5 and 12 energies alike. Bisected
to the cosine rather than the density, then to the solver tolerance:

| cosz | default tolerance | explicit 1e-6 / 1e-8 / 1e-10 |
|---|---|---|
| -1.0 | **fails** | all succeed |
| -0.9 | **fails** | all succeed |
| -0.95, -0.85, -0.8, -0.5, -0.1 | succeeds | succeeds |

`Set_h_max` and an explicit GSL stepper do not help; only setting
`rel_error`/`abs_error` does. **cosz = -1.0 is the first point of
`OSC/100x100`, and cosz = -0.9 is the whole of `CHORD/12x1`** --- so at the
harness's default knob of 0, which is what `if p.knob > 0` left at nuSQuIDS'
own defaults, nuSQuIDS could not have produced either Earth figure.

This went unseen because the old cosine rescaling multiplied the angle by
1 - 1.4e-7, which happens to land off the failing value. The adapter was
working by luck, and the luck was a factor we were about to remove for
unrelated and correct reasons.

**Fixed**: the tolerance is now always set explicitly, with knob = 0 meaning
1e-7 rather than "whatever the constructor picked". Both grids now run,
including cosz = -1.0.

### 12. This library's rtol knob was advertised four decades past where it works

`manifest.json` declared `rtol: 1e-3..1e-15` and the adapter declared a knob
domain down to -15. Measured on `CHORD/12x1`:

| rtol | result |
|---|---|
| 1e-3, 1e-4, 1e-5 | returns |
| 1e-6 through 1e-15 | **raises** --- the slab product's error estimate bottoms out near 7e-8 at the 1024-slab ceiling |

Ten of the thirteen declared tolerance settings crash. A sweep over the
declared domain would have died partway through, and declaring a domain that
cannot be swept is the mirror image of the omission objection LBL-3 was about
--- there we failed to try a setting a code exposed, here we advertised
settings this library does not reach.

**Fixed**: the declared domain is now the reachable band, in both the adapter
and the manifest, with the reason recorded in both. The small residuals are
reached through the slab dial, which does run its full 1..256.

## Still open

* The per-code reference builder (`README.md` item 3). Nothing is differenced
  against anything until it exists, and after the de-scaling above it is the
  only thing that can validate any code's accuracy number.
* Finding 6, NuFast-Earth's `single_trajectory_mode`, so objection Earth-4 can
  be answered.
* Finding 7, `n_layers` and `shells_total` in the artifact, for Earth-1.
* Finding 10, the manifest saying whether a reference is built on the
  continuous PREM or on the discretised profile the code was handed.
* The Python-side accuracy path: `runner.py` still accepts only
  `amortized|throughput`, and no Python adapter has a `probabilities()`.

---

# Corrections from the independent second pass

An independent review was run over the same tree, focused on how each code is
used. Its claims were re-verified here against the pinned sources before being
recorded. Three correct my own.

## Correction A: fix 2 was partial --- THROUGHPUT is still not comparable

`Set_delta` clears only `probabilities_calculated` (`NuFastEarth.cpp:103-111`).
So adding it to `evaluate()` stopped the memcpy, but the eigens and the
internal amplitudes --- the amplitude product over roughly 2000 chord segments
--- still run once at warm-up and are reused by every timed repetition. The
per-repetition cost is now one `Inner_Amplitude_to_Probability` rotation per
grid point.

Every other adapter redoes its full pipeline per repetition: NuFast-LBL
recomputes everything, Prob3++ re-propagates, GLoBES re-multiplies its
S-matrices, nuSQuIDS re-evolves, nuCraft re-integrates, and this library
re-solves. So a NuFast-Earth THROUGHPUT number remains flattering and
non-comparable, on the code whose author raised the objections.

`Set_Eigenvalue_Precision` **does** clear eigens, amplitudes and probabilities
together (verified at `NuFastEarth.cpp`), so a full invalidation is available.
But which of these is right is a **protocol decision, not a bug fix**, and it
is left open deliberately:

* have `bench.hpp::throughput` call `configure` before every repetition ---
  changes what the protocol measures for every code, including GLoBES, whose
  `configure` recomputes its mixing matrix;
* have the NuFast-Earth adapter invalidate its whole pipeline per
  `evaluate()` --- but the adapter is not told which protocol it is running
  under, so this would also change AMORTIZED, where refresh-only is correct
  and is the author's own definition;
* keep it and stamp the series as a cached refresh, never mixable with the
  other six codes' throughput.

AMORTIZED is unaffected under all three: refresh-only is exactly what the
author's `Atmospheric_Speed` times.

## Correction B: nuSQuIDS has an exact constant-density mode we do not enable

`Set_AllowConstantDensityOscillationOnlyEvolution` exists in the pinned wheel
(verified by probing it) and makes `EvolveState` skip the ODE entirely on a
`ConstantDensity` body, evolving each node algebraically. Measured here on
`CONST/60E`: it agrees with the ODE path to **7.2e-7**, which is the ODE's own
tolerance, and both paths are finite.

So on both constant-density grids nuSQuIDS is being run through a numerical
integrator when its authors ship a closed-form path for exactly that case.
This is objection LBL-3's shape --- an exact mode that exists and was never
tried --- aimed at a different code, and it is ours to fix, not theirs.

## Correction C: nuSQuIDS interpolates with an Akima spline, not linearly

Finding 10 above, and `_aligned_radii`'s own docstring, describe the table
handed to `nsq.EarthAtm` as piecewise-linear. The pinned source declares
`AkimaSpline inter_density` in `body.h`. This matters for the reference: an
as-handed reference for nuSQuIDS must reproduce that spline, including the
1e-9 km jump-holding node pairs, or Akima overshoot near the density jumps ---
small, structured, and shaped like solver error --- lands in its residual.

## Confirmed additions

* `conversions.mass_defect('nuSQuIDS')` raises `KeyError`: nothing extracts
  its `sqrt(2) G_F N_A`, so the builder has no potential to construct for it.
* `oscillation_parameters` has no sterile mixing at all, so the manifest's
  `energies_gev_3plus1` grid --- and the nuCraft exception, which only bites at
  3+1 --- cannot be run from the manifest as written.
* nuCraft's own `prem` model carries `rICore = 1121.5` (`NuCraft.py:122`)
  while its own profile table steps at 1221.5 (`:103`). A typo in the pinned
  source. The adapter's override to `earth.PREM_BOUNDARIES[0]` is right and
  must stay.
* `OBJECTIONS.md` says `test_bench_objections.py` reads its table. No such
  file exists; the reader is
  `test_bench_pipeline.py::test_every_objection_names_a_test_that_exists_or_is_declared_pending`.
* Objection coverage as it stands: **4 of 9** answerable with what is built
  (LBL-3, Earth-2, Earth-3, Other-2); 3 structurally blocked (LBL-1's looped
  control has no machinery, Earth-1, Earth-4); 2 waiting only on the
  orchestrator (LBL-2's stamp, Other-1's `manifest_sha256`).

---

# The accuracy axis, and four corrections that came with it

Everything above concerned whether the machinery could produce a fair number.
This section is what happened when it produced them. No timing run has
happened yet; the accuracy protocol reads no clock, so all of this was
measured on a machine that was not idle, which is exactly why it could be done
first.

## The reference basis, decided by measuring the alternative

The question was what to difference each code against. Two candidates:

* **as-handed** --- the discretised profile each code was actually given;
* **continuous** --- the real PREM, which is the physical Earth.

Measured, on one chord at three energies:

| code | vs as-handed | vs continuous |
|---|---|---|
| Prob3++ (1024 shells) | 1.3e-13 | **4.06e-7** |
| GLoBES (1024 shells) | 8.1e-13 | **4.06e-7** |
| NuOscProbExact (256 slabs/segment) | 6.1e-13 | 4.3e-7 |

The compiled codes land on the *same* 4.06e-7 to three figures, because that
number is the discretisation they were both handed and not a property of
either. Under as-handed they look six orders better than this library on
identical physics; under continuous all three sit in one class.

Both bases are internally symmetric --- the earlier unfairness came from
**mixing** them, not from either one. The choice was made on what the axis
means. As-handed asks "does a code evaluate its own approximation
self-consistently", which they all do, and it makes every discretisation dial
flat because the reference moves with the dial. Continuous asks how close the
answer is to the Earth, which is the question a user faces, and it keeps the
dials meaningful. **Continuous, for all seven.**

## Correction: the flag that does not mean what it says

`constant_shells` in NuFast-Earth's `Earth_Density` reads as though it selects
between constant and varying shells, and this file previously recorded that
the code "natively supports a continuous profile". It does not. `Mean_Density`
builds **one slab per interval between declared discontinuities** whatever the
flag says; the flag chooses only whether each slab carries the midpoint or the
mean density. The discontinuity list *is* the discretisation.

Measured the hard way: declaring the nine PREM boundaries with the
varying-shell flag --- which looks exactly like asking for continuous
treatment --- gives ten slabs for the whole chord and an error of **3.04e-2**,
five orders worse than the 1024-shell configuration it replaced. With the fine
cut restored it returns to 3.7e-7 to 3.0e-8, the same class as everyone else.

Neither density choice is uniformly better: mean wins at two of three energies
and loses threefold at the third. Three points settle nothing and the notebook
says so.

## Correction: a Romberg estimate reports the corner it discards

Two references appeared to floor eight orders short, at 3.3e-10 and 2.4e-10.
They were fine. The error estimate is the standard Romberg one --- how far the
corner moves when the last ladder rung is added --- which bounds the error of
the corner *one extrapolation level shallower* than the one returned. With a
four-rung ladder that shallower corner sat two levels back at 3e-10 while the
returned corner was already at 2e-13. **The quoted number was never the
reference's error**; it was the error of a corner computed and thrown away,
overstating the truth by about 1700x.

Worse, it hid a real repair. Cutting the chord at every spline knot --- which
a spline needs, since its density has a kink at each knot and a slab
straddling one has no convergent expansion --- appeared to buy a factor of
1.8. Measured against the value now established, the pre-repair corner was
9.6e-10 from the truth and the post-repair corner is 1.9e-13: three orders,
concealed by the estimator. A fifth rung fixes the estimate. nuCraft
5.79e-10 -> 1.91e-13, nuSQuIDS 2.78e-10 -> 4.99e-14.

The lesson is narrower than "check your work": a single summary number can be
wrong in a way the quantity it summarises is not, and the table it came from
is what says so.

## Correction: nuCraft does not floor at 3e-3

`manifest.json` asserted that floor and attributed it to the code's own
constants. Nobody had turned the dial. Swept at three flavors against its own
reference:

| numPrec | 1e-2 | 5e-4 (default) | 1e-6 | 1e-8 | 1e-10 |
|---|---|---|---|---|---|
| error | 1.6e-4 | 7.2e-6 | 1.3e-6 | 3.6e-8 | **7.9e-11** |

**This row is wrong and is superseded.** Re-measured by driving the
adapter in process at cos(thz) = -0.9 against a reference good to
2.4e-18, nuCraft's error falls 5.3e-4, 8.0e-6, 9.3e-7, 3.7e-7 and then
stops -- 4.0e-7 at numPrec 1e-10, no better than at 1e-8.  There is a
floor near 4e-7.  The correction above replaced one untested claim
(3e-3) with another (no floor above 7.9e-11); both were wrong, in
opposite directions, and both were about someone else's code.

Still falling threefold at the edge of the declared domain. That is zvode at
`atol = rtol = numPrec * 2e-3` and nothing else. The claim was wrong by seven
orders and it was aimed at another author's code.

Its constants exception survives but is now derived rather than asserted:
nuCraft's own comment defines one prefactor for `(2Ye, 0, 0, 1-Ye)`, so its two
constants must be exactly 2:1 by its own formula. `15.256e-5` is exactly twice
that prefactor at the atomic mass unit; `7.6525e-5` would need 928.5 MeV,
which is not the amu, the proton, the neutron, or the 0.939 its own comment
states. Inconsistent with its own derivation, most likely a digit typo for
`7.6325e-5` --- and irrelevant to every three-flavor figure, since the sterile
entries do not exist at three flavors.

## Correction: nuSQuIDS' exact mode is not uniformly better

This repository's adapter said the exact constant-density mode and the ODE
"agree to 7.2e-7, which is the ODE's own tolerance". Backwards. The algebraic
path is exact at 59 of 60 nodes, median 4e-16, but carries a
**tolerance-independent** excursion of 7.2e-7 between roughly 3.4 and 6.5 GeV;
the ODE at 1e-12 reaches 3.5e-12 everywhere. The 7.2e-7 is the exact mode's own
defect. Reproduced independently on a 12-energy grid at 2.3e-7 with median
2.9e-16. Any figure quoting nuSQuIDS on constant density must say which mode
produced it.

## What the accuracy machinery now does

* `runner.py` accepts `accuracy`, untimed exactly as `bench.hpp` is, and all
  three Python adapters implement `probabilities()`. The loader refuses an
  adapter without one, because a code that cannot be checked for accuracy
  should not be usable for speed.
* All **59 declared knob settings across seven codes run** --- probed before
  sweeping, because our own rtol domain once had ten of thirteen settings
  raising and a sweep that dies at the last rung wastes everything before it.
* `sweep_accuracy.py` measures every code over every knob, caching references
  on the key insight that a reference depends on the code and the grid point
  but **not** on the knob: the knob changes what the code does, not what the
  right answer is. A nine-point sweep costs one reference, not nine.
* Objection Earth-1 is answerable: `Problem` carries `n_layers` beside `knob`,
  and every artifact states `shells_total`.
* Objection LBL-1 is answerable: `--loop` drives NuFast-LBL's own entry point
  one energy at a time, and batched and looped return bit-identical
  probabilities --- the same code, the same physics, a different call pattern.
* Artifacts carry `conventions`, `profile_basis`, `n_layers`, `shells_total`
  and `looped`.
* ADVERSARIAL's sixth attack is closed: this library now **refuses to run** if
  its compiled kernels are not live, and records which kernel was entered and
  how many threads with every artifact.

First real accuracy data from the new pipeline, NuFast-LBL on constant
density, which answers objection LBL-3 outright:

| N_Newton | -1 (exact) | 0 | 1 | 2 | 3 |
|---|---|---|---|---|---|
| max deviation | 7.6e-16 | 1.3e-5 | 1.9e-9 | 1.0e-15 | 7.6e-16 |

The exact mode reaches round-off and three Newton steps match it, on the
setting the earlier comparison never tried.

## Deferred: BOTH nuCraft oscillogram cells

`amortized_nuCraft_OSC_100x100_knob6_scan-dmsq31.json` is deliberately not
measured.  Every other code carries a Dmsq31 cell on every grid it runs;
this is the one gap.

At 12800 profile nodes it costs about seventy-five minutes on its own, more
than the twenty-two other cells of its batch put together, and the machine
was wanted.  It feeds no figure -- nuCraft draws a single line, and the
oscillogram cells answer the hundred-by-hundred grid objection rather than
the parameter one -- and its chord pair WAS measured, giving the ratio of
about 1.01 that the table needs.

It is recorded here because the sensitivity table now covers six codes on
the chord and five on the oscillogram.  Publishing that without saying so
would read as nuCraft having been dropped for being slow, which is the exact
appearance this comparison exists to avoid.  Run the cell, or say in the
caption that the grid is not covered for it.

    python tests/bench/run_all.py --speed-only \
      --only "nuCraft_OSC_100x100_knob6_scan-dmsq31"

**Added later, and more serious than the above.** The delta_CP oscillogram
cell, `amortized_nuCraft_OSC_100x100_knob6.json`, is ALSO unmeasured at
12800 nodes.  A six-minute deadline killed it at 5515 s, roughly a minute
short, and it is recorded FAILED -- so the artifact on disk is the old
200-node one dated 08-31 05:46.

That one is not safe to defer quietly.  It feeds a figure (nuCraft's point
in the 100x100 comparison), and it leaves this code with its chord numbers
on the corrected profile and its oscillogram number on the profile whose
interpolation error was the plateau.  One code, two measurement sets.

Either re-run it, or drop nuCraft from the oscillogram comparison and say
so.  Do not publish it as it stands.

    python tests/bench/run_all.py --speed-only --force \
      --only "nuCraft_OSC_100x100"
