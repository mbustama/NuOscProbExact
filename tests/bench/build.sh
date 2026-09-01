#!/usr/bin/env bash
# Clone, verify and build every external code at the commit or hash pinned in
# manifest.json, into a tree this script owns.  Idempotent: re-running it
# rebuilds nothing that is already present and correct.
#
# Two profiles are built for every compiled code, because one flag set cannot
# serve both axes: `speed` uses each project's own upstream flags, so we never
# publish a code slower than its authors do, and `accuracy` uses -O3 without
# fast-math, because value-changing flags make an accuracy measurement measure
# the flag instead of the code.
#
# Nothing here needs root.  GSL is taken from the system; everything else is
# built or unpacked under .bench-build/, and the Python side lives in
# .venv-bench/ from requirements.lock.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BUILD="$REPO/.bench-build"
SRC="$BUILD/src"; PREFIX="$BUILD/prefix"; BIN="$BUILD/bin"; LOG="$BUILD/log"
VENV="$REPO/.venv-bench"
MANIFEST="$REPO/tests/bench/manifest.json"
mkdir -p "$SRC" "$PREFIX" "$BIN" "$LOG"

pin() { python3 -c "
import json,sys
d=json.load(open('$MANIFEST'))
c=[c for c in d['codes'] if c['name']=='$1'][0]
print(c$2)
"; }

say() { printf '\n\033[1m== %s\033[0m\n' "$*"; }

# ---------------------------------------------------------------- git clones
clone_at() {                       # name url commit dir
  local name="$1" url="$2" commit="$3" dir="$SRC/$4"
  if [ -d "$dir/.git" ] && [ "$(git -C "$dir" rev-parse HEAD)" = "$commit" ]; then
    echo "  $name already at $commit"; return
  fi
  rm -rf "$dir"
  git clone -q "$url" "$dir"
  git -C "$dir" checkout -q "$commit"
  local got; got="$(git -C "$dir" rev-parse HEAD)"
  [ "$got" = "$commit" ] || { echo "  $name: WRONG COMMIT $got"; exit 1; }
  echo "  $name pinned at $commit"
}

# ------------------------------------------------------------------ tarballs
fetch_verify() {                   # name url sha256 file
  local name="$1" url="$2" want="$3" file="$SRC/$4"
  if [ -s "$file" ] && [ "$(sha256sum "$file" | cut -d' ' -f1)" = "$want" ]; then
    echo "  $name tarball already verified"; return
  fi
  curl -sSL -o "$file" "$url"
  local got; got="$(sha256sum "$file" | cut -d' ' -f1)"
  [ "$got" = "$want" ] || { echo "  $name: SHA256 MISMATCH"; echo "    want $want"; echo "    got  $got"; exit 1; }
  echo "  $name tarball sha256 verified"
}

say "1/7  NuFast-LBL"
clone_at NuFast-LBL "$(pin NuFast-LBL "['url']")" "$(pin NuFast-LBL "['pin']['commit']")" nufast-lbl
for prof in speed accuracy; do
  flags="-Ofast -ffast-math -std=c++17"; [ "$prof" = accuracy ] && flags="-O3 -std=c++17"
  # -Dmain=... because the released file carries its own demo main() at :307;
  # renaming it makes the translation unit linkable beside the harness's main
  # instead of merely compilable.
  g++ $flags -Dmain=nufast_lbl_demo_main -c "$SRC/nufast-lbl/NuFast_LBL.cpp" \
      -o "$BUILD/nufast_lbl.$prof.o"
  echo "  $prof profile ok  [$flags] -> nufast_lbl.$prof.o"
done

say "2/7  NuFast-Earth"
clone_at NuFast-Earth "$(pin NuFast-Earth "['url']")" "$(pin NuFast-Earth "['pin']['commit']")" nufast-earth
for prof in speed accuracy; do
  flags="-O3 -std=c++17"      # upstream Makefile CFlags; same for both profiles here
  out="$BUILD/nufast-earth-$prof"; mkdir -p "$out"
  # src/main.cpp is upstream's own demo entry point.  Skipped, because
  # bench.hpp owns main() -- linking both is a duplicate-symbol error, which is
  # precisely the mechanism that stops an adapter timing itself.
  for f in "$SRC/nufast-earth/src"/*.cpp; do
    [ "$(basename "$f")" = main.cpp ] && continue
    g++ $flags -c -I"$SRC/nufast-earth/include" "$f" -o "$out/$(basename "${f%.cpp}").o"
  done
  echo "  $prof profile ok  [$flags]  ($(ls "$out"/*.o | wc -l) objects)"
done

say "3/7  Prob3++"
clone_at Prob3++ "$(pin Prob3++ "['url']")" "$(pin Prob3++ "['pin']['commit']")" prob3
for prof in speed accuracy; do
  cf="-O2"; [ "$prof" = accuracy ] && cf="-O3"
  out="$BUILD/prob3-$prof"; mkdir -p "$out"; ( cd "$SRC/prob3"
    for c in mosc.c mosc3.c; do gcc $cf -c "$c" -o "$out/${c%.c}.o"; done      # MUST be C
    for c in EarthDensity.cc BargerPropagator.cc; do g++ $cf -c "$c" -o "$out/$(basename "${c%.cc}").o"; done )
  echo "  $prof profile ok  [$cf]  ($(ls "$out"/*.o | wc -l) objects)"
done

say "4/7  GLoBES"
fetch_verify GLoBES "$(pin GLoBES "['url']")" "$(pin GLoBES "['pin']['sha256']")" globes.tar.gz
if [ ! -f "$PREFIX/lib/libglobes.so" ] && [ ! -f "$PREFIX/lib/libglobes.a" ]; then
  tar xzf "$SRC/globes.tar.gz" -C "$SRC"
  ( cd "$SRC"/globes-3.2.18 \
    && ./configure --prefix="$PREFIX" --with-gsl-prefix=/usr >"$LOG/globes_configure.log" 2>&1 \
    && make -j"$(nproc)" >"$LOG/globes_make.log" 2>&1 \
    && make install >>"$LOG/globes_make.log" 2>&1 )
fi
ls "$PREFIX"/lib/libglobes* >/dev/null 2>&1 && echo "  installed into $PREFIX" || { echo "  GLoBES BUILD FAILED, see $LOG"; exit 1; }

say "5/7  nuCraft"
fetch_verify nuCraft "$(pin nuCraft "['url']")" "$(pin nuCraft "['pin']['sha256']")" nucraft.tar.gz
# The r22 tarball is FLAT -- NuCraft.py and friends sit at the top level with
# no wrapping directory -- so --strip-components=1 silently strips the
# filenames themselves and extracts nothing while still exiting 0.
rm -rf "$SRC/nucraft"; mkdir -p "$SRC/nucraft"
tar xzf "$SRC/nucraft.tar.gz" -C "$SRC/nucraft"
test -f "$SRC/nucraft/NuCraft.py" || { echo "  nuCraft: NuCraft.py not extracted"; exit 1; }
python3 -c "import os;print('  extracted:', sorted(f for f in os.listdir('$SRC/nucraft') if f.endswith(('.py','.txt'))))"

say "6/7  nuSQuIDS (pinned wheel; GSL/HDF5/SQuIDS bundled upstream)"
. "$VENV/bin/activate"
python -m pip install -q "nuSQuIDS==$(pin nuSQuIDS "['pin']['tag'].lstrip('v')")" 2>&1 | tail -2 || true
python - <<'PYV'
from importlib.metadata import version
import nuSQuIDS  # noqa: F401  -- the import itself is the check
got = version('nuSQuIDS')
print('  import ok, distribution version', got)
assert got == '1.13.3', 'nuSQuIDS is %s, manifest pins 1.13.3' % got
PYV

say "7/7  benchmark adapters (compile + link against the harness)"
# Each C++ adapter is compiled beside tests/bench/bench.hpp (which owns
# main() and every clock) and linked against the per-profile objects built
# above, into $BIN/bench_<code> (speed) and $BIN/bench_<code>_accuracy.
# conversions.h is regenerated first so no adapter can compile against a
# stale factor.  The exit status of every compile is checked for real: a
# failure names its log and stops the build.
python3 "$REPO/tests/bench/gen_conversions.py" > "$BUILD/conversions.h"
ADAPT="$REPO/tests/bench/adapters"
build_adapter() {                  # label command...
  local label="$1"; shift
  if "$@" > "$LOG/adapter_$label.log" 2>&1; then
    echo "  PASS $label -> $(basename "${!#}")"
  else
    local status=$?
    echo "  FAIL $label (exit $status, log: $LOG/adapter_$label.log)"
    exit "$status"
  fi
}
for prof in speed accuracy; do
  sfx=""; [ "$prof" = accuracy ] && sfx="_accuracy"

  # NuFast-LBL: upstream's own flags for speed, -O3 for accuracy, matching
  # the object it links against.
  lblflags="-Ofast -ffast-math -std=c++17"; [ "$prof" = accuracy ] && lblflags="-O3 -std=c++17"
  build_adapter "nufast_lbl.$prof" g++ $lblflags -I"$BUILD" \
    "$ADAPT/nufast_lbl.cpp" "$BUILD/nufast_lbl.$prof.o" \
    -o "$BIN/bench_nufast_lbl$sfx"

  # NuFast-Earth: -O3 -std=c++17 is both its upstream flag set and the
  # accuracy profile, so the two binaries differ only in provenance.
  build_adapter "nufast_earth.$prof" g++ -O3 -std=c++17 -I"$BUILD" \
    -I"$SRC/nufast-earth/include" "$ADAPT/nufast_earth.cpp" \
    "$BUILD/nufast-earth-$prof"/*.o -o "$BIN/bench_nufast_earth$sfx"

  # Prob3++: upstream -O2 for speed, -O3 for accuracy, as built above.
  p3flags="-O2 -std=c++17"; [ "$prof" = accuracy ] && p3flags="-O3 -std=c++17"
  build_adapter "prob3.$prof" g++ $p3flags -I"$BUILD" -I"$SRC/prob3" \
    "$ADAPT/prob3.cpp" "$BUILD/prob3-$prof"/*.o -o "$BIN/bench_prob3$sfx"
done
# GLoBES has ONE profile: the measured code lives in libglobes, built once
# by upstream's configure at its default -O2; relinking the thin adapter
# with different flags would not change the code being measured.  Recorded
# here rather than hidden, like the nuSQuIDS wheel exception.
build_adapter "globes" g++ -O2 -std=c++17 -I"$BUILD" -I"$PREFIX/include" \
  "$ADAPT/globes.cpp" -L"$PREFIX/lib" -Wl,-rpath,"$PREFIX/lib" \
  -lglobes -lgsl -lgslcblas -lm -o "$BIN/bench_globes"

say "our_prem.h (three adapters include it)"
python3 "$REPO/tests/external_drivers/gen_prem_header.py" > "$BUILD/our_prem.h"
test -s "$BUILD/our_prem.h" && echo "  generated $(wc -l < "$BUILD/our_prem.h") lines" || { echo "  FAILED"; exit 1; }

say "build record"
python3 - <<PY
import json, platform, subprocess, os
def sh(c):
    try: return subprocess.run(c, shell=True, capture_output=True, text=True).stdout.strip().splitlines()[0]
    except Exception: return "unknown"
def sh_all(c):
    try: return subprocess.run(c, shell=True, capture_output=True, text=True).stdout.splitlines()
    except Exception: return []
rec = {
  "schema": "bench-build/1",
  "built_on": sh("date -Is"),
  "repo_commit": sh("git -C '$REPO' rev-parse HEAD"),
  "manifest_sha256": sh("sha256sum '$MANIFEST' | cut -d' ' -f1"),
  "cxx": sh("g++ --version"), "cc": sh("gcc --version"),
  "gsl": sh("gsl-config --version"),
  "cpu": sh("grep -m1 'model name' /proc/cpuinfo | cut -d: -f2"),
  "cores": os.cpu_count(),
  "python": platform.python_version(),
  "profiles": {"speed": "each project's own upstream flags",
               "accuracy": "-O3 -std=c++17, no fast-math"},
  "exceptions": {"nuSQuIDS": "measured as distributed (upstream wheel); its flags are the wheel's, so the two-profile rule cannot apply to it -- recorded rather than hidden"},
  # Whether a built object can use more than one thread, asked of the object
  # rather than assumed.  The bench notes have long said the four compiled
  # codes link no threading library; this records the evidence next to the
  # build instead of asserting it in prose, and it says nothing about the
  # codes driven from Python, whose parallelism is NumPy's BLAS and is
  # captured per run by runner.thread_environment() instead.
  "threading_libraries": {
      obj: [ln.split()[0] for ln in sh_all("ldd '%s' 2>/dev/null" % obj)
            if any(t in ln for t in ("libgomp", "libomp", "pthread", "libopenblas", "libmkl"))]
      for obj in sorted(__import__("glob").glob("$BUILD/*.o")
                        + __import__("glob").glob("$BUILD/*.so")
                        + __import__("glob").glob("$BUILD/bench_*"))
  },
}
json.dump(rec, open("$BUILD/build_record.json", "w"), indent=2)
print(json.dumps(rec, indent=2))
PY
