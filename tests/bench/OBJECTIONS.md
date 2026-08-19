# Every objection, and the measurement that answers it

Nine points were raised against the comparison by an author of NuFast, on
2026-08-07. Each is recorded here with the protocol, grid and knob that
answers it, the artifact field the answer lands in, and the test that fails if
that field is missing. `test_bench_objections.py` reads this table, so an
objection cannot be quietly dropped: removing a row is a code change, and
leaving a row unanswered fails a test.

Verdicts were established by reading the pinned upstream sources and release
notes, not by taking either side's word. Two of our own earlier conclusions
were wrong and are recorded as such.

| id | objection | verdict | protocol / grid / knob | artifact field | test |
|---|---|---|---|---|---|
| LBL-1 | NuFast-LBL batches over E since v2.0.0; paper says it takes one energy per call | **holds against the paper, not the measurement** — our driver did batch; three text sites are wrong | THROUGHPUT, `CONST/N-sweep`, batched vs looped control | `fig08.series[NuFast-LBL].batched` and a `looped` series | `test_batched_codes_are_driven_batched`, `test_capabilities_match_registry` |
| LBL-2 | Precision limited by unit conversion, not by his refinement | **partly holds** — his G_F suspicion is wrong for us (we used his own `1.52588e-4`), but Fig. 9 has ħc deliberately native | accuracy, `CONST/60E`, `conventions: matched` vs `native` | `conventions` stamp on every accuracy artifact | `test_conventions_stamp_present`, `test_no_hand_typed_conversion_constants` |
| LBL-3 | `N_Newton < 0` gives exact eigenvalues and was never tested | **holds — his strongest point** | accuracy + AMORTIZED, `CONST/60E`, `N_Newton ∈ {-1,0,1,2}` | a point at every knob-domain endpoint | `test_every_registered_knob_extreme_was_swept` |
| Earth-1 | Is 256 shells total or per region? | **does not hold** — per major layer, 4×256 = 1024, as the paper says | AMORTIZED, `CHORD/12x1`, `n_layers` | `knob.n_layers` plus `shells_total` | `test_shell_count_semantics_recorded` |
| Earth-2 | Setters inside the timed loop defeat his caching | **holds, and worse than alleged** — we also constructed the engine and the Earth object per repetition, and our own `rhoYe` injected O(N²) work into his inner loop | AMORTIZED (all invariants hoisted by the harness contract) | `protocol.name == AMORTIZED`, plus the O(1) callback test | `test_driver_timer_ownership`, `test_rhoye_is_constant_time`, `test_rhoye_costs_no_more_than_the_shipped_one` |
| Earth-3 | One zenith angle hides his zenith caching | **holds** — ~8× per-probability at 1024 layers on 100×100 | AMORTIZED, **`OSC/100x100`** alongside `CHORD/12x1` | a series under each grid stamp | `test_oscillogram_grid_present` |
| Earth-4 | NuFast-Earth does have a constant-density mode since v1.1.0 | **holds** — `single_trajectory_mode` is in the pinned source | accuracy + AMORTIZED, `CONST/60E` via single-trajectory mode | a NuFast-Earth series in the constant-density artifacts | `test_constant_density_capable_codes_appear_in_const_artifacts` |
| Other-1 | Pin the versions used | **holds** | n/a | `manifest_sha256` in every artifact; commit/sha256 per code | `test_artifacts_share_one_manifest`, `test_every_code_is_pinned` |
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
