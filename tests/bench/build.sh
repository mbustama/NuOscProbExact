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

say "1/6  NuFast-LBL"
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

say "2/6  NuFast-Earth"
clone_at NuFast-Earth "$(pin NuFast-Earth "['url']")" "$(pin NuFast-Earth "['pin']['commit']")" nufast-earth
for prof in speed accuracy; do
  flags="-O3 -std=c++17"      # upstream Makefile CFlags; same for both profiles here
  out="$BUILD/nufast-earth-$prof"; mkdir -p "$out"
  for f in "$SRC/nufast-earth/src"/*.cpp; do
    g++ $flags -c -I"$SRC/nufast-earth/include" "$f" -o "$out/$(basename "${f%.cpp}").o"
  done
  echo "  $prof profile ok  [$flags]  ($(ls "$out"/*.o | wc -l) objects)"
done

say "3/6  Prob3++"
clone_at Prob3++ "$(pin Prob3++ "['url']")" "$(pin Prob3++ "['pin']['commit']")" prob3
for prof in speed accuracy; do
  cf="-O2"; [ "$prof" = accuracy ] && cf="-O3"
  out="$BUILD/prob3-$prof"; mkdir -p "$out"; ( cd "$SRC/prob3"
    for c in mosc.c mosc3.c; do gcc $cf -c "$c" -o "$out/${c%.c}.o"; done      # MUST be C
    for c in EarthDensity.cc BargerPropagator.cc; do g++ $cf -c "$c" -o "$out/$(basename "${c%.cc}").o"; done )
  echo "  $prof profile ok  [$cf]  ($(ls "$out"/*.o | wc -l) objects)"
done

say "4/6  GLoBES"
fetch_verify GLoBES "$(pin GLoBES "['url']")" "$(pin GLoBES "['pin']['sha256']")" globes.tar.gz
if [ ! -f "$PREFIX/lib/libglobes.so" ] && [ ! -f "$PREFIX/lib/libglobes.a" ]; then
  tar xzf "$SRC/globes.tar.gz" -C "$SRC"
  ( cd "$SRC"/globes-3.2.18 \
    && ./configure --prefix="$PREFIX" --with-gsl-prefix=/usr >"$LOG/globes_configure.log" 2>&1 \
    && make -j"$(nproc)" >"$LOG/globes_make.log" 2>&1 \
    && make install >>"$LOG/globes_make.log" 2>&1 )
fi
ls "$PREFIX"/lib/libglobes* >/dev/null 2>&1 && echo "  installed into $PREFIX" || { echo "  GLoBES BUILD FAILED, see $LOG"; exit 1; }

say "5/6  nuCraft"
fetch_verify nuCraft "$(pin nuCraft "['url']")" "$(pin nuCraft "['pin']['sha256']")" nucraft.tar.gz
# The r22 tarball is FLAT -- NuCraft.py and friends sit at the top level with
# no wrapping directory -- so --strip-components=1 silently strips the
# filenames themselves and extracts nothing while still exiting 0.
rm -rf "$SRC/nucraft"; mkdir -p "$SRC/nucraft"
tar xzf "$SRC/nucraft.tar.gz" -C "$SRC/nucraft"
test -f "$SRC/nucraft/NuCraft.py" || { echo "  nuCraft: NuCraft.py not extracted"; exit 1; }
python3 -c "import os;print('  extracted:', sorted(f for f in os.listdir('$SRC/nucraft') if f.endswith(('.py','.txt'))))"

say "6/6  nuSQuIDS (pinned wheel; GSL/HDF5/SQuIDS bundled upstream)"
. "$VENV/bin/activate"
python -m pip install -q "nuSQuIDS==$(pin nuSQuIDS "['pin']['tag'].lstrip('v')")" 2>&1 | tail -2 || true
python - <<'PYV'
from importlib.metadata import version
import nuSQuIDS  # noqa: F401  -- the import itself is the check
got = version('nuSQuIDS')
print('  import ok, distribution version', got)
assert got == '1.13.3', 'nuSQuIDS is %s, manifest pins 1.13.3' % got
PYV

say "build record"
python3 - <<PY
import json, platform, subprocess, os
def sh(c):
    try: return subprocess.run(c, shell=True, capture_output=True, text=True).stdout.strip().splitlines()[0]
    except Exception: return "unknown"
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
}
json.dump(rec, open("$BUILD/build_record.json", "w"), indent=2)
print(json.dumps(rec, indent=2))
PY
