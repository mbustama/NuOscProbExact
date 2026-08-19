# The fair-comparison benchmark pipeline

Start here. This file is the entry point and deliberately holds no facts of its
own: everything authoritative lives in the files it points at, because a second
copy is a copy that drifts.

## Why this exists

An audit of the previous comparison found nine ways it was unfair or
unreproducible, and an author of one of the compared codes raised nine
objections of his own. Five of his held. The worst defect was ours and he had
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

1. `probabilities()` in the four C++ adapters. The accuracy path was added to
   the contract after they were written, so they currently do not link --- which
   is the contract working, not a bug to route around.
2. Remove the input rescaling from the adapters. They still apply
   `OUR_PREM_MASS_DEFECT_*` to densities, which was right when all codes shared
   one reference. References are now per-code, so those factors belong to the
   reference builder and applying them in both places double-counts.
3. The per-code reference builder. Each code is differenced against a 50-digit
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
