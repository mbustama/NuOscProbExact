# Every objection, and the measurement that answers it

Nine points were raised against the comparison by an author of NuFast, on
2026-08-07. Each is recorded here with the protocol, grid and knob that
answers it, the artifact field the answer lands in, and the test that fails if
that field is missing. `test_bench_pipeline.py::test_every_objection_names_a_test_that_exists_or_is_declared_pending`
reads this table, so an objection cannot be quietly dropped: removing a row is
a code change, and leaving a row unanswered fails a test.  (An earlier version
of this file named a `test_bench_objections.py` that has never existed.)

Verdicts were established by reading the pinned upstream sources and release
notes, not by taking either side's word. Two of our own earlier conclusions
were wrong and are recorded as such.

| id | objection | verdict | protocol / grid / knob | artifact field | test |
|---|---|---|---|---|---|
| LBL-1 | NuFast-LBL batches over E since v2.0.0; paper says it takes one energy per call | **holds against the paper, not the measurement** — our driver did batch; three text sites are wrong | THROUGHPUT, `CONST/N-sweep`, batched vs looped control | `fig08.series[NuFast-LBL].batched` and a `looped` series | `test_capabilities_claims_match_the_registry`, and PENDING `test_batched_codes_are_driven_batched` |
| LBL-2 | Precision limited by unit conversion, not by his refinement | **partly holds** — his G_F suspicion is wrong for us (we used his own `1.52588e-4`), but Fig. 9 has ħc deliberately native | accuracy, `CONST/60E`, `conventions: matched` vs `native` | `conventions` stamp on every accuracy artifact | `test_no_adapter_types_a_physical_constant`, and PENDING `test_conventions_stamp_present` |
| LBL-3 | `N_Newton < 0` gives exact eigenvalues and was never tested | **holds — his strongest point** | accuracy + AMORTIZED, `CONST/60E`, `N_Newton ∈ {-1,0,1,2}` | a point at every knob-domain endpoint | `test_precision_knobs_include_their_exact_modes`, and PENDING `test_every_registered_knob_extreme_was_swept` |
| Earth-1 | Is 256 shells total or per region? | **does not hold** — per major layer, 4×256 = 1024, as the paper says | AMORTIZED, `CHORD/12x1`, `n_layers` | `knob.n_layers` plus `shells_total` | `test_shell_count_semantics_recorded` |
| Earth-2 | Setters inside the timed loop defeat his caching | **holds, and worse than alleged** — we also constructed the engine and the Earth object per repetition, and our own `rhoYe` injected O(N²) work into his inner loop | AMORTIZED (all invariants hoisted by the harness contract) | `protocol.name == AMORTIZED`, plus the O(1) callback test | `test_no_adapter_owns_a_clock`, and PENDING `test_rhoye_is_constant_time`, `test_rhoye_costs_no_more_than_the_shipped_one` |
| Earth-3 | One zenith angle hides his zenith caching | **holds** — ~8× per-probability at 1024 layers on 100×100 | AMORTIZED, **`OSC/100x100`** alongside `CHORD/12x1` | a series under each grid stamp | `test_oscillogram_grid_present` |
| Earth-4 | NuFast-Earth does have a constant-density mode since v1.1.0 | **holds** — `single_trajectory_mode` is in the pinned source | accuracy + AMORTIZED, `CONST/60E` via single-trajectory mode | a NuFast-Earth series in the constant-density artifacts | `test_constant_density_capable_codes_appear_in_const_artifacts` |
| Other-1 | Pin the versions used | **holds** | n/a | `manifest_sha256` in every artifact; commit/sha256 per code | `test_every_code_is_pinned`, `test_every_code_records_how_latest_was_verified`, and PENDING `test_artifacts_share_one_manifest` |
| Other-2 | Repeat many times, report mean and a standard deviation | **holds** — every published figure was a minimum of a few reps | both protocols | `{mean, sd, min, n}`; no scalar-time field exists in the schema | `test_no_bare_scalar_times` |

## Things he did not raise that the audit found

These are ours, and are tracked here because they were worse than what he
alleged.

| id | problem | fix |
|---|---|---|
| OWN-1 | Six published numbers about external codes exist only in Python docstrings with no stored run — worst is the claim that nuCraft's floor is its constants, evidenced by an uncommitted "forced ratio" run the paper quotes as `about 3e-7` | six committed evidence artifacts, each with its generator; `test_every_benchmark_number_in_the_paper_has_an_artifact` |
| OWN-2 | The driver behind `speed_accuracy.json` is not in the repository, and `nufast_scan.json` cites a `tests/nufast_scan.md` that does not exist | `generated_by` must name a path that exists; `test_generated_by_exists` |
| OWN-3 | Three different values of one NuFast-Earth density factor coexist (applied `0.99209238`, correct `0.9920928368`, documented `0.9920938`) | one generator, `conversions.py`; literals forbidden in drivers |
| OWN-4 | Two datasets for the same code built with different flags, presented as one methodology | two declared build profiles; `test_artifacts_share_one_manifest` |
| OWN-5 | The paper states nuSQuIDS aligned floor `9e-7` while the drawn curve floors at `1.02e-6` | numbers rendered from artifacts, never typed |

## Corrections to our own earlier findings

- We did **not** drive NuFast-LBL one energy per call. The committed driver
  passes a vector of up to 30 000 energies. The paper's prose is wrong; the
  curve is fair. The unexplained factor-of-two fall in that figure **is**
  amortization, established by a looped control flat at ~56 ns.
- The `rhoYe` handed to NuFast-LBL was **not** derived from Prob3++'s
  constant. NuFast-LBL genuinely carries `1.52588e-4`, the same rounding
  Prob3++ uses; our factor was right.

## PENDING

These are named above but cannot exist yet, because each inspects an artifact
and no artifact exists until a measured run has happened.  They are listed here
so the gap between what is promised and what is built stays visible instead of
implied; `test_every_objection_names_a_test_that_exists_or_is_declared_pending`
fails if a name appears in the table without being either implemented or listed
here.

Needing artifacts from a run:
`test_batched_codes_are_driven_batched`,
`test_conventions_stamp_present`,
`test_every_registered_knob_extreme_was_swept`,
`test_shell_count_semantics_recorded`,
`test_oscillogram_grid_present`,
`test_constant_density_capable_codes_appear_in_const_artifacts`,
`test_artifacts_share_one_manifest`,
`test_no_bare_scalar_times`,
`test_generated_by_exists`,
`test_every_benchmark_number_in_the_paper_has_an_artifact`.

Needing the NuFast-Earth adapter to be timed against its shipped equivalent:
`test_rhoye_is_constant_time`,
`test_rhoye_costs_no_more_than_the_shipped_one`.

What exists today, and what each already forbids:
`test_every_code_is_pinned` (a month is not a version),
`test_every_code_records_how_latest_was_verified` (the paper's version is not
necessarily the latest),
`test_precision_knobs_include_their_exact_modes` (a sentinel must be in the
swept domain),
`test_no_adapter_owns_a_clock` (only the harness times),
`test_no_adapter_types_a_physical_constant` (one generator, no fourth value),
`test_capabilities_claims_match_the_registry` (no claiming an interface does
not exist while calling it),
`test_the_shared_parameter_set_has_exactly_one_home` (one set of mixing
parameters for every code),
`test_conversion_factors_are_derived_from_the_pinned_sources`.


## Our measured answers, moved here from outside the repository

These are **our** measurements, made on 2026-08-19 on the i5-1334U following
the author's own harness protocol.  They lived only in a private note until
now, which meant the sharpest numbers in this dispute existed nowhere that
could be re-checked and would have had to be measured again if that note were
lost.  His letter stays out of this repository, because it is private
correspondence; our own numbers are ours to publish and belong here.

Every figure below predates the rebuilt pipeline and was taken under the old
harness.  They are recorded as the answers that were given, to be replaced by
the sweep the new pipeline produces --- not treated as final.

### NuFast-LBL, driven fairly

* Batched, per probability: **37.1 +- 0.4 ns** on the plateau above ten
  energies, against our published 35 ns.  That endpoint stands.
* Looped, one energy per call: flat at **~56 ns**.  Batching therefore buys a
  clean **1.5x** by amortizing roughly 20 ns of per-call setup.  This is what
  establishes that the fall from 72 to 35 ns in the throughput figure **is**
  amortization --- which the paper explicitly says it has not established.
  The `--loop` control series now in `bench.hpp` exists to reproduce this
  under the new pipeline rather than cite it.
* Newton sweep at a thousand energies: `-1` 68.5 ns, `0` 37.7, `1` 44.2,
  `2` 50.3, `3` 56.0.  So its exact mode costs **1.82x** the cheapest setting,
  and the author's "only about 25% slower" does not hold for NuFast-LBL.  It
  does hold for NuFast-Earth, where the exact precision setting is +29%.

### NuFast-Earth: where our published 403.6 us went

Decomposed against the same driver run today:

* x1.34 machine state alone (the same driver now gives 301 us);
* **x2.7 our own harness** --- our density lookup scanned up to 1024
  discontinuities on every call, from inside his engine's inner loop.  His own
  O(1) model under the same protocol: 110.6 us;
* x1.10 setup left inside the timed region.

**The defensible one-call figure is about 110 us, not 403.6**, which weakens
the paper's "130 us ours against 400 us his" to rough parity at that corner.
Driven as designed, with the invariants hoisted and its caching working, it
reaches **0.057 us per probability** --- some 2000x below our published number
and independent of layer count.

Grid shape matters as much: at 1024 layers a 100x100 grid is about **8x
cheaper per probability** than the 12x1 we published, because trajectories and
eigensolves are shared across zenith angles.  That is objection Earth-3, in
numbers.

### Independent support for his round-off claim

At constant density, **GLoBES and NuFast-LBL at its exact setting agree to all
seventeen digits** --- both returning checksum 307.08408653736365.  Two
unrelated codes converging exactly is strong evidence that its exact mode
reaches double precision, which is what he claimed and what our sweep,
stopping at two Newton steps, could never have shown.

### Six published claims with no stored evidence

An audit of every external-code number in the paper against the stored data
found that every number in the throughput figure and every accuracy column in
the comparison table **is** supported.  These six are not, and each is a claim
about somebody else's code:

1. nuCraft's forced-ratio residual, published as "about 3e-7, which is what
   its solver is worth" --- exists only as `2.8e-7` in two source comments.
   This is the entire evidence for the sharpest claim we make against another
   author's code.  It is also now known to be wrong in a second way: nuCraft
   does not floor at 3e-3 either, and its error tracks its tolerance dial down
   to 7.9e-11.
2. nuCraft's constant-profile agreement, "about 3e-11" --- docstring only.
3. nuSQuIDS grid alignment, "3e-6 improving to 9e-7" --- unstored, and the
   drawn curve floors at 1.02e-6, which rounds to 1e-6.
4. nuSQuIDS: "the same value survives five different steppers" --- no evidence
   stored.
