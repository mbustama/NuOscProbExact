# The adversarial check, and what it must try to break

Written before the check runs, so its scope cannot shrink to whatever it
happened to find. The pipeline's tests were written by whoever wrote the code
they guard and therefore encode the same blind spots; this is the only pass that
addresses that.

The framing is deliberate. A reviewer asked to *review* returns "looks good" or
invents problems. A reviewer asked to **break a stated guarantee** either breaks
it, which is worth knowing, or fails to, which is evidence. Every item below is
phrased as something to falsify.

## Six things to attack

1. **Write an adapter that produces an unfair number and passes every test in
   `tests/test_bench_pipeline.py`.** If that is possible, the invariants are
   insufficient and the gap is the finding. Try: hoisting something that should
   be timed, claiming a batched entry point while looping, reaching a constant
   by a route the literal-grep misses.

2. **Check every reference is built with the same constants as the code it
   judges.** A reference in the wrong convention makes a code look *better* than
   it is, and nothing downstream would notice. For each code, verify against the
   pinned source that the reference used that code's own hbar c and its own
   matter normalisation, and that all codes got identical mixing parameters.
   nuCraft is the deliberate exception: its constants for units, but the correct
   isoscalar ratio of 1/2, because its own 0.50161 is an error rather than a
   convention and absorbing it would hide the defect.

3. **Find anything wrong left over from the previous timing setup.** That round
   had nine defects and some of its code survives. Look for: setters inside a
   timed region, an engine or profile object constructed per repetition, our own
   O(N^2) work injected into someone else's inner loop, a stale factor, a driver
   whose output file no longer matches what produced it, or a number in an
   artifact with no generator.

4. **Show the strategy cannot answer one of Peter's claims.** Read his nine
   points in `OBJECTIONS.md` and try to demonstrate that some claim remains
   unanswerable with what is built --- not merely unanswered yet. Partial
   coverage is the failure mode here: nine points answered eight ways.

5. **Find a code being run in a form slower or dumber than its best.** Every
   code must use its batched, optimised entry point and be swept to its best
   precision setting. Check each adapter against its own pinned headers for an
   entry point or knob we are not using.

6. **Find NuOscProbExact being under-represented.** The mirror-image failure,
   and the one nobody else would report. It must run batched, through the
   compiled kernels, with numba live. Already verified by hand, and worth
   re-checking because none of it is currently *asserted* at run time:

   | check | state |
   |---|---|
   | `fastkernels.HAVE_NUMBA` / `USE_NUMBA` / `available()` | all true |
   | numba threads | 12 |
   | `CONST/60E` enters | `probabilities_3nu_kernel` |
   | Earth grid enters | `earth_chords_3nu_kernel` |
   | batching | energies and zenith together |
   | invariants | hoisted into `setup()` |

   The gap: nothing asserts any of it during a run. If numba were absent, a
   threshold moved, or `USE_NUMBA` flipped, the pipeline would quietly measure
   the NumPy path and publish it as ours. Fix: assert at setup, and record the
   kernel actually entered plus the thread count in every artifact.

## Rules for the run

- **Before the measurements**, not after: a flaw found after an hours-long sweep
  costs the sweep.
- **No benchmarks.** Reading, throwaway adapters in the scratchpad, and the test
  suite only --- its own compiling would spoil the idle machine.
- **Its findings get verified before they are believed.** On this project's
  record roughly half of an agent's headline claims have needed correcting: a
  false alarm about a constant derived from the wrong code, and a claim that we
  drove NuFast-LBL one energy at a time when the driver batched. Each finding is
  checked against the source before it is repeated to anyone.
