# The fair-comparison benchmark pipeline

Start here. This file is the entry point and deliberately holds no facts of its
own: everything authoritative lives in the files it points at, because a second
copy is a copy that drifts.

## Why this exists

An audit of the previous comparison found nine ways it was unfair or
unreproducible, and an author of one of the compared codes raised nine
objections of his own. Six of his held outright, one partly, and one against the paper's prose
rather than its measurements; `OBJECTIONS.md` is the authority on which. The worst defect was ours and he had
not noticed it: our driver's density lookup scanned up to 1024 discontinuities
on every call from inside his engine's inner loop, adding O(N^2) work to his
measurement and inflating his published point by about 2.7x.

## Where each thing is written down

| Question | File |
|---|---|
| What is measured, how, and with which pins | `manifest.json` |
| Which objection each measurement answers | `OBJECTIONS.md` |
| What the adversarial check must break | `ADVERSARIAL.md` |
| The invariants, executable | `../test_bench_pipeline.py` |
| The contract every compiled code implements | `bench.hpp` |
| The same contract for Python codes | `runner.py` |
| Where every physical constant comes from | `conversions.py` |
| Whether a run may be believed | `machine.py` |
| How the seven codes are fetched and built | `build.sh` |

Run `pytest tests/test_bench_pipeline.py` before reading anything. Twelve
invariants pass; each one forbids a specific mistake that was actually made.

## State of play

Built and proven: the harness owning every clock, all seven adapters, the
Python runner, the pinned build of all seven codes, the conversions generator,
the machine guard, and the invariant tests.

**Not yet built**, in dependency order:

1. `probabilities()` in the four C++ adapters **and an accuracy path on the
   Python side** --- `runner.py` accepts only `amortized|throughput` and no
   Python adapter has a `probabilities()` method, so that half is not done
   either. The scored channel is `numu->numu`; see `manifest.json`. The accuracy path was added to
   the contract after they were written, so they currently do not link --- which
   is the contract working, not a bug to route around.
2. Remove the input rescaling from the adapters --- **all of it, not only the
   density**. The manifest absorbs *both* hbar c and the matter normalisation
   into the reference, so `OUR_COSZ_HBARC_SCALE` (in `prob3.cpp`, `globes.cpp`,
   `nusquids.py`) and the Python-side scalings in `nucraft.py` and
   `nusquids.py` come out too. They still apply
   `OUR_PREM_MASS_DEFECT_*` to densities, which was right when all codes shared
   one reference. References are now per-code, so those factors belong to the
   reference builder and applying them in both places double-counts.
3. The per-code reference builder. Before designing it, read
   `memory/nuoscprob-prem-figure-referee.md`: the slab product converges at
   **second** order, the Richardson referee is `(4*P(256) - P(128))/3`, and the
   floor is ~1e-11. A first-order referee already cost one session. Commit
   b549a6d also decided a shared-problem comparison is kept for the Earth
   figure --- that decision currently survives only in that commit message. Each code is differenced against a 50-digit
   reference built with its own hbar c and matter normalisation and the shared
   mixing parameters. This is the piece where an error is invisible: a reference
   in the wrong convention makes a code look *better* than it is.
4. The orchestrator: the protocol x grid x knob matrix, artifact emission,
   canary bracketing, and the environment capture `machine.py` provides.
5. The adversarial check, per `ADVERSARIAL.md`, before any measurement.

## Two things a reader must not re-derive wrongly

Both were reported during this work and are **false**:

* *"The wrong rhoYe was handed to NuFast-LBL, derived from Prob3++'s
  constant."* No. NuFast-LBL genuinely carries `1.52588e-4`, the same rounding
  Prob3++ uses. The factor was right.
* *"We drove NuFast-LBL one energy per call."* No. The committed driver passes
  a vector of up to 30000 energies. The paper's prose was wrong, not the
  measurement --- and the unexplained factor-of-two in that figure **is**
  amortization, established by a looped control flat at ~56 ns.

## Before measuring

The CPU governor on this machine is `powersave`, which needs root to change and
is the largest single lever on timing stability. Either set it to `performance`
first or record it and lean on the canary. The session canary baseline is
0.0508 us per probability at ten thousand points.

Wall clock is the open question: 30 blocks x 2 timed protocols x 7 codes x knob
sweeps runs to hours, and the nuSQuIDS cells on the 100x100 grid are the
expensive ones. Decide tiering before starting, not during.

## Gaps a cold reader found, still open

A fresh agent was given only this directory and the manifest and asked to reach
the runs unaided.  It got most of the way and then listed what it could not
determine.  Four of those are fixed above and in `manifest.json`; these are not,
and each one is a place where a wrong guess is silent rather than loud.

- **`n_layers` cannot actually be swept.**  `Problem` carries a single `int
  knob`, and `nufast_earth.cpp` hard-codes 256 layers with a comment referring
  to a "knob's sibling" that does not exist.  The manifest's PREM layer sweep
  therefore cannot run as written.  Either widen `Problem` or drop the sweep,
  but do not leave the manifest promising it.
- **The orchestrator's cell matrix is undefined.**  Nothing says which code runs
  which grid under which protocol, nor how artifacts are named.  `run_all.py`
  has to invent it, and an invented matrix is not reproducible.
- **Our answers to the NuFast author's points are measurements held in one
  place only** --- `memory/nufast-author-claims-verbatim.md`, outside this
  repository and deliberately so, because his letter is private.  His prose
  must stay out; *our own numbers* are ours, and belong in `OBJECTIONS.md` in
  the repository where they can be re-checked.  Move them.  If that memory file
  is lost before they are moved, every one of them has to be measured again.
- **A loose thread nobody has pulled.**  `bench_nufast_earth`'s chord checksum
  sits about 0.05 per probability away from `prob3` and `globes` at *every*
  knob value.  A constant offset that survives refinement is not a convergence
  error; it is a convention or a channel mismatch, and it should be understood
  before it becomes a figure.
- **One number in `manifest.json` disagrees with what was measured.**  The
  manifest says nuCraft floors at 3e-3; the measurement recorded at 3+1 flavors
  was 3.7e-4.  Neither has been re-run to settle it, and the manifest should
  not be trusted here until one has been.
