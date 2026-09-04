# Handover: build NuOscProbExact's Figure 11 in the Magnus paper, with Magnus on it

**To:** the Magnus session (`/home/mbustamante/Research/magnus`, package `magnuspy` v1.0.0)
**From:** the NuOscProbExact session (`/home/mbustamante/Research/NuOscProb/NuOscProbExact`, branch `dev-fair-benchmarks`, commit `26a372c`, 2026-09-01)

You are asked to reproduce **Figure 11 of the NuOscProbExact paper exactly** — same
colours, sizes, styles, layout, annotations, data — and to **add curves for Magnus**
to it. Nothing about the existing curves may move, be re-measured, re-scaled, or
re-styled. You are adding, not rebuilding.

This document is self-contained. You do not need the session that wrote it. Every
claim in it that you could get wrong silently is marked **VERIFY** with the command
that settles it. Run those first.

---

## 0. The one-paragraph version

Figure 11 is `speed_accuracy_combined.pdf`: three stacked speed–accuracy panels on
one shared time axis, drawn from three frozen JSON files by a single cell in
`notebooks/10_paper_figures.ipynb`. Copy the three JSON files, the whole
`tests/bench/` provenance tree, and the figure code. Measure Magnus on the same
grids under the same protocols on **the same laptop**, append its series to the
three JSON files under a new name, give it a colour/marker slot no existing code
uses, and redraw. The physics conventions need no work at all: Magnus and
NuOscProbExact carry byte-identical unit constants and PREM coefficients, which
means they share a reference. That is the single biggest simplification here, and
§3 tells you how to confirm it rather than trust it.

---

## 1. Six facts that will silently ruin the figure if you get them wrong

Not warnings — each of these has already cost someone a session in this project.

1. **The timings are machine-specific and the machine is this laptop.** Every
   number in Figure 11 was measured on an Intel Core i5-1334U (10 cores, 12
   threads, up to 4.6 GHz, 16 GB). If you measure Magnus anywhere else, the
   Magnus curve is not comparable to the others and the figure misrepresents
   performance. Magnus's own repo already records this hazard, in
   `notebooks/external_speed_accuracy.json` under `_caveat`: *"Timings are that
   machine, not this one, and no single rescaling maps between them."* There is
   no rescaling. Measure on this laptop or do not draw a joint figure.

2. **The accuracy axis is measured against a per-code reference, not a shared
   one.** Each code is differenced against its own 50-digit `mpmath` reference,
   built with *that code's* ħc and matter normalisation and the *shared* mixing
   parameters. A reference in the wrong convention makes a code look **better**
   than it is, and nothing in the output flags it. See §3 and §7.2.

3. **The scored channel for accuracy is `numu->numu`.** Not `numu->nue`. The
   `CONST/N-sweep` grid entry in `manifest.json` says `numu->nue`, but that is
   the *throughput* figure's problem statement, not the accuracy channel. Getting
   this wrong silently changes every accuracy number. `manifest.json ->
   scored_channel` is the authority.

4. **There are two constant-density measurement sets and they must never be
   mixed.** Figure 11's top panel reads `tests/speed_accuracy.json` — the
   benchmark pipeline's set, conventions matched, per-code references. Figure 10
   of the same paper reads `tests/const_density_scan.json`,
   `tests/nufast_scan.json` and `tests/nusquids_scan.json` — a 150-energy scan
   with ħc deliberately *unmatched*, which is what that figure is about. A number
   from one set must never appear in a sentence about the other. Take the top
   panel's numbers from `speed_accuracy.json` and nothing else.

5. **The 3+1 panel comes from a different pipeline than the other two.** Four
   flavours are measured by `tests/prem_scan.py` and by nothing else; the
   three-flavour panel comes from `tests/bench/`. They are kept in separate files
   on purpose — different measurement sets, against different references — and a
   single file holding both invites them being read as one. Keep them separate in
   Magnus too.

6. **Points are ordered along the curve by the dial, not by cost.** Sorting a
   flat-cost series by cost sorts it by measurement noise and draws a U that reads
   as two curves. `by_dial()` exists for exactly this. See §6.5.

---

## 2. What you are copying, exactly

All paths below are relative to the NuOscProbExact repo root,
`/home/mbustamante/Research/NuOscProb/NuOscProbExact`.

### 2.1 The three files the figure actually reads (mandatory)

| Source | Bytes | Feeds | Suggested destination in Magnus |
|---|---|---|---|
| `tests/speed_accuracy.json` | 22 199 | **top** panel (constant density) | `notebooks/external_speed_accuracy_const.json` |
| `tests/bench/earth_plane.json` | 24 659 | **middle** panel (PREM, 3 flavours) | `notebooks/external_earth_plane.json` |
| `tests/prem_speed_accuracy.json` | 26 548 | **bottom** panel (PREM, 3+1), key `sterile_3plus1` | `notebooks/external_prem_speed_accuracy.json` |

**Magnus already has older copies of two of these and they are stale — replace
them, do not merge.** `notebooks/external_prem_speed_accuracy.json` in Magnus has
`density_scale_nusquids = 0.99209238370534` and a `costhz_handed_to_nusquids`
field; the current file has `0.9920926641043969` and instead carries
`dm2_scale_nusquids`, `dm2_scale_nucraft` and `nucraft_nc_over_cc`. They are
different measurement sets. `notebooks/external_speed_accuracy.json` in Magnus has
6 series; the current `tests/speed_accuracy.json` has 8.

### 2.2 The provenance tree (mandatory — this is the "how it was generated" material)

Copy `tests/bench/` **entire**, preserving structure. It is 34 tracked files plus
517 artifacts:

| Group | Files | Bytes | What it is |
|---|---|---|---|
| `tests/bench/manifest.json` | 1 | 30 331 | **The authority.** Protocols, thread policy, reference conventions, grids, code pins, scored channel, oscillation parameters. |
| `tests/bench/README.md` | 1 | 6 914 | Entry point. **Partly stale — see §10.** |
| `tests/bench/OBJECTIONS.md` | 1 | 11 398 | Which measurement answers which objection from the NuFast author. |
| `tests/bench/FINDINGS.md` | 1 | 34 976 | What the audit found. |
| `tests/bench/ADVERSARIAL.md` | 1 | 4 140 | What the adversarial check must break. |
| `tests/bench/adapters/*` | 8 | 73 712 | One adapter per code (4 C++, 4 Python). How each code is actually driven. |
| `tests/bench/artifacts/*.json` | 517 | 808 582 | Every raw measurement cell, with its environment and statistics. |
| `tests/bench/*.py`, `bench.hpp`, `build.sh`, `requirements.lock` | ~14 | ~200 000 | The harness: runner, orchestrator, reference builder, conversions, machine guard, build script. |
| `tests/bench/accuracy_*.json`, `speed_accuracy_plane.json`, `scan_sensitivity.json`, `reference_audit.json` | 6 | 643 249 | Intermediate measurement sets that `emit_figures.py` reduces into the three files in §2.1. |

Also copy, because the middle and bottom panels depend on them:

| Source | Bytes | Why |
|---|---|---|
| `tests/external_drivers/README.md` | 5 517 | **How the three compiled codes were driven**, every trap hit, and the exact build commands. This is the single most useful prose document for a reader who wants to dispute the figure. |
| `tests/external_drivers/*` (16 files) | ~47 000 | The drivers themselves plus their raw `.txt` output. `tests/prem_scan.py` reads the `.txt` files and does not need the codes. |
| `tests/prem_scan.py` | 43 161 | Produces the 3+1 panel. |
| `tests/test_bench_pipeline.py` | 13 437 | Twelve invariants, each forbidding a mistake actually made. Run it after copying. |

**Total to copy: ~1.9 MB across ~570 files.** None of it is generated at build time;
all of it is committed and frozen.

### 2.3 The figure code

From `notebooks/make_notebooks.py` (this is a *generator* — the notebooks are built
from it, never edited directly):

| What | Anchor in the file | Lines at `26a372c` |
|---|---|---|
| `rcParams` block + `COLW` + `FIGSIZE_SQUARE` | search `COLW = 229.5/72.27` | 2327–2390 |
| `STYLE` (top panel styles) | `STYLE = {"NuOscProbExact": ("-o", "C3", 4.0),` | 3148–3153 |
| `const_points()` | `def const_points(name):` | 3162–3180 |
| `PREM_STYLE` (Earth panel styles) | `PREM_STYLE = {"NuOscProbExact (tolerance)":` | 3360–3381 |
| `by_dial()` | `def by_dial(points):` | 3383–3408 |
| `prem_axes()` | `def prem_axes(ax, panel, annotations, ...)` | 3410–3481 |
| **The figure cell itself** | markdown heading `## The three speed-accuracy planes together` | 3882–4063 |

Line numbers will drift; the anchors will not. Use the anchors.

---

## 3. What Magnus shares with NuOscProbExact — verified, but verify it again

This is the finding that makes the job tractable. The two libraries carry
**identical** unit constants and **identical** PREM coefficients:

| Constant | NuOscProbExact `src/globaldefs.py` | Magnus `src/magnus/globaldefs.py` |
|---|---|---|
| `CONV_KM_TO_INV_EV` | `5.06773e9` | `5.06773e9` |
| `GF` | `1.1663787e-23` | `1.1663787e-23` |
| `MASS_PROTON` | `938.272046e6` | `938.272046e6` |
| `MASS_NEUTRON` | `939.565379e6` | `939.565379e6` |
| `ELECTRON_FRACTION_EARTH_CRUST` | `0.5` | `0.5` |
| `VCC_EARTH_CRUST` | `sqrt(2)*GF*n_e` | `sqrt(2)*GF*n_e` |

and in `earth.py` both carry the same ten-row `_PREM_COEFFS` table and the same
`PREM_BOUNDARIES = [1221.5, 3480.0, 5701.0, 5771.0, 5971.0, 6151.0, 6346.6,
6356.0, 6368.0]`.

**Consequence:** Magnus needs **no** new reference, **no** density rescaling, and
**no** ħc correction. The NuOscProbExact reference *is* Magnus's reference. Magnus's
points drop onto the existing panels directly.

**VERIFY before relying on it** — run this from the Magnus repo and require it to
print `IDENTICAL` three times:

```bash
python3 - <<'PY'
import sys, numpy as np
sys.path.insert(0, '/home/mbustamante/Research/NuOscProb/NuOscProbExact/src')
import globaldefs as ngd, earth as nearth
sys.path.insert(0, '/home/mbustamante/Research/magnus/src')
from magnus import globaldefs as mgd, earth as mearth
for name in ('CONV_KM_TO_INV_EV','GF','MASS_PROTON','MASS_NEUTRON',
             'ELECTRON_FRACTION_EARTH_CRUST','VCC_EARTH_CRUST'):
    a, b = getattr(ngd, name), getattr(mgd, name)
    assert a == b, (name, a, b)
print('IDENTICAL constants')
assert np.array_equal(nearth._PREM_COEFFS, mearth._PREM_COEFFS)
print('IDENTICAL PREM coefficients')
assert np.array_equal(nearth.PREM_BOUNDARIES, mearth.PREM_BOUNDARIES)
print('IDENTICAL PREM boundaries')
PY
```

If any assertion fails, **stop** and treat Magnus as a code needing its own
reference, following `manifest.json -> reference_conventions` and
`tests/bench/reference.py`. Do not paper over a difference with a scale factor:
that is precisely the mistake the current pipeline was built to undo.

---

## 4. The machine, and how to know a run is admissible

- **Required CPU:** Intel Core i5-1334U, 10 cores / 12 threads, 4.6 GHz max, 16 GB.
  Confirm with `lscpu | grep 'Model name'`.
- **Governor:** this machine idles at `powersave`, which needs root to change and
  is the largest single lever on timing stability. Either set it to `performance`
  first, or record it and lean on the canary.
- **Session canary baseline:** 0.0508 µs per probability at ten thousand points.
  `tests/bench/machine.py` is the guard; it decides whether a run may be believed.
- **Per-cell admissibility:** the pipeline rejects a cell whose block coefficient
  of variation exceeds **0.10**. `block_cv` is recorded per point in
  `tests/speed_accuracy.json`. Magnus's cells must clear the same bar.
- Roughly **36 %** of the artifacts in the existing set are marked inadmissible and
  retained anyway, for audit. Do not silently drop yours; record them.

---

## 5. Thread policy — publish two Magnus series if Magnus is parallel

From `manifest.json -> thread_policy`, verbatim in substance:

- Each code runs with the parallelism it naturally uses. **No code is forced to one
  thread.** Forcing every code to one thread strips this library's parallelism
  while leaving the compiled codes untouched and publishes a number nobody would
  experience.
- The four compiled codes (GLoBES, NuFast-Earth, NuFast-LBL, Prob3++) are
  single-threaded **by construction** — verified with `ldd`: no `libgomp`, no
  threaded BLAS, no `pthread`. This was verified, not assumed.
- NuOscProbExact uses `NUMBA_NUM_THREADS`, 12 on this machine.
- **Every series must record the thread count the code actually used, not the
  policy.** A cross-code speed claim without both counts is not a claim.
- NuOscProbExact publishes both: `NuOscProbExact` and `NuOscProbExact (1 thread)`
  are both in `tests/speed_accuracy.json`. The 1-thread series is **measured but
  not drawn** on the top panel, because at constant density it lands at 1.1704 µs
  against 1.1703 — a second legend entry exactly under the first marker.

**For Magnus:** if Magnus is parallel, measure and freeze both series. Decide
separately whether to draw both — draw the second only if it is visually
distinguishable, and say so in the caption if you do.

Record for each Magnus run: `cpu_count`, `NUMBA_NUM_THREADS`, `OMP_NUM_THREADS`,
`MKL_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `NUMEXPR_NUM_THREADS`, and
`threadpoolctl.threadpool_info()`. `tests/bench/runner.py::thread_environment()`
returns exactly this; reuse it rather than rewriting it.

---

## 6. The figure specification — reproduce exactly

### 6.1 Style

Set these before anything is drawn. This is the group's standard `matplotlibrc`,
inlined so the figure does not depend on an external resource file. Sizes are
chosen for a figure included at the paper's `\columnwidth` beside 10 pt body text.

```python
COLW = 229.5/72.27          # \columnwidth of elsarticle 5p twocolumn, in inches
HAVE_LATEX = shutil.which('latex') is not None

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Palatino'],
    'text.usetex': HAVE_LATEX,
    'font.size': 8,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7.5,
    'axes.linewidth': 0.7,
    'lines.linewidth': 1.1,
    'xtick.top': True, 'xtick.bottom': True, 'xtick.direction': 'in',
    'ytick.left': True, 'ytick.right': True, 'ytick.direction': 'in',
    'xtick.major.size': 3.5, 'xtick.minor.size': 2.0,
    'ytick.major.size': 3.5, 'ytick.minor.size': 2.0,
    'xtick.major.width': 0.7, 'xtick.minor.width': 0.7,
    'ytick.major.width': 0.7, 'ytick.minor.width': 0.7,
    'xtick.major.pad': 3.0, 'ytick.major.pad': 3.0,
    'legend.frameon': True, 'legend.edgecolor': 'black',
    'legend.framealpha': 1.0, 'legend.fancybox': False,
    'legend.borderpad': 0.4, 'legend.handlelength': 1.9,
    'legend.handletextpad': 0.5, 'legend.columnspacing': 0.9,
    'legend.labelspacing': 0.35,
    'axes.grid': False,
    'axes.xmargin': 0.0,
    'figure.dpi': 150, 'savefig.dpi': 300,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})
if not HAVE_LATEX:
    plt.rcParams['font.serif'] = plt.rcParamsDefault['font.serif']
```

Without a LaTeX installation the labels fall back to matplotlib's mathtext: the
figure is the same figure, the typesetting is not. Palatino comes from the TeX
installation, not from matplotlib.

### 6.2 Canvas

```python
FIGW_WIDE = 469.75/72.27          # \textwidth: this is a figure* spanning both columns
fig, axes = plt.subplots(3, 1, sharex=True, figsize=(FIGW_WIDE, 6.46))
```

Shared x across all three panels; that sharing is the whole point of the figure, so
that a cost carries by eye from panel to panel.

```python
XLIM_ALL = (3.0e-2, 6.0e5)
LABEL_TR = (0.985, 0.95, "right", "top")
```

`XLIM_ALL` runs to `6e5`, not `1e5`. At `1e5` it silently dropped six **measured**
points off the right edge — Prob3++ at 49152 and 65536, GLoBES at 32768, 49152 and
65536, and the sterile panel's rtol 1e-8. Three of those are GLoBES's most accurate
settings, so the cut was not cosmetic. **Do not tighten it to make room for
Magnus.** If Magnus lands outside, widen it and say so.

Per-panel y limits: top `(1e-16, 1e-2)`, middle `(1e-10, 2e-1)`, bottom
`(2e-11, 1e-1)`. Y ticks: top `[10**k for k in range(-16,-1)]`, middle
`range(-10,0)`, bottom `range(-11,0)`, all at `labelsize=6.2`.

Every panel: `ax.set_axisbelow(True)` and
`ax.grid(True, axis="both", which="major", color="0.78", lw=1.0, zorder=0)` — one
vertical rule per decade of time, behind the data.

Bottom axis label `r"Time per probability [$\mu$s]"`. One shared y label:

```python
fig.supylabel(r"Error vs.\ the reference of each panel,  "
              r"$\max_\alpha |\Delta P_{\nu_\mu \to \nu_\alpha}|$",
              fontsize=9, x=0.030)
fig.subplots_adjust(hspace=0.06, left=0.10, right=0.995, top=0.995, bottom=0.055)
```

`x=0.030` is set by hand because the default sits on the tick labels at this width.
The y label says **"the reference of each panel"** and not "converged solution",
because the reference differs per panel and per code. Keep that wording.

### 6.3 Panel 1 — constant density

Source `tests/speed_accuracy.json`, selected through `const_points(name)`.
Problem: L = 1300 km, ρ = 3 g cm⁻³, three flavours, ν_μ row, largest error over
three channels and 60 energies from 0.6 to 20 GeV, against an exact 50-digit
reference.

```python
STYLE = {"NuOscProbExact":         ("-o", "C3", 4.0),
         "nuSQuIDS":               ("-v", "C2", 3.6),
         "NuFast-LBL":             ("-D", "C4", 3.2),
         "GLoBES":                 ("-*", "C6", 6.0),
         "Prob3++":                ("-P", "C5", 4.4),
         "Second-order expansion": ("-s", "C1", 3.6)}
```

`dict` order is legend order. Any series whose name starts with `NuOscProbExact`
gets `mfc="white", mew=1.0, zorder=5`; everything else `zorder=4`.

`best_at_this_accuracy` marks the fastest point at a given accuracy; only those
are drawn, because four of these codes have an inert dial at constant density and
plotting the rest stacks points at one height. Magnus must set this flag on its own
points, by the same rule.

`const_points()` keeps only points with `best_at_this_accuracy` true and a non-null
`max_abs_error`, then sorts by `(-max_abs_error, us_per_probability)` — **along the
curve by error, not by knob**. Sorting by knob put `N_Newton = -1` (the exact mode,
at the accurate end) before `0` (the least accurate), so the polyline ran bottom,
top, bottom and drew a closed loop reading as two NuFast-LBL curves.

Four of these codes have an **inert dial** at constant density — no profile to
subdivide, so every shell setting returns the same answer. Drawing them all says
nothing; drawing the slowest reports the code as slower than anyone would run it.
NuFast-LBL keeps its curve because its dial is the eigenvalue solve and does move
the accuracy.

`Second-order expansion` is not an external code: it is the α–s₁₃ series the paper
cites, evaluated by `tests/bench/adapters/second_order.py` in this library's
conventions against this library's reference. It used to be a hand-carried point
with no generator, no version and no formula recorded.

**NuFast-Earth is measured at constant density and the artifact is kept, but it is
not drawn here** — this is the constant-density comparison and an Earth propagator
does not belong on it.

Also on this panel: a double-precision rule at `2.2e-16`
(`color="0.5", ls=":", lw=0.7, zorder=1`) labelled `"Double precision"` at
`(5.0e1, 3.2e-16)`, `fontsize=6.0, color="0.4"`; seven annotations on NuFast-LBL and
nuSQuIDS; a boxed subtitle at `LABEL_TR`; legend `loc="lower right"`,
`bbox_to_anchor=(0.995, 0.03)`, `fontsize=5.6`, `framealpha=0.85`, `ncol=2`,
frame linewidth 0.7. Copy the annotation list verbatim from the source cell.

### 6.4 Panels 2 and 3 — the Earth planes

Both drawn through the **single** routine `prem_axes()`, deliberately: the
single-panel figures and this stacked one must not disagree about which
measurements they show. The caller owns the x range.

```python
PREM_STYLE = {"NuOscProbExact (tolerance)":      ("-o", "C3", 4.0),
              "NuOscProbExact":                  ("--s", "C3", 3.2),
              "NuOscProbExact (double-double)":  ("--s", "C3", 3.2),
              "NuOscProbExact (eigensolver)":    (":^", "C1", 3.8),
              "nuSQuIDS":                        ("-v", "C2", 3.6),
              "NuFast-Earth":                    ("-D", "C4", 3.2),
              "NuFast-Earth (dCP only)":         ("--h", "C4", 3.4),
              "GLoBES":                          ("-*", "C6", 6.0),
              "Prob3++":                         ("-P", "C5", 4.4),
              "nuCraft":                         ("-s", "C0", 3.4)}
```

The style is keyed off the **frozen series name**, not the legend label; `relabel`
renames for the legend only. The tolerance dial is the primary curve (the user
states the quantity the y axis measures); the explicit slab count is the same code
told the answer instead of asked to find it, drawn broken behind it.

**Panel 2** (`tests/bench/earth_plane.json`): PREM chord at cos θ_z = −0.9,
12 energies from 3 to 40 GeV, L = 11468 km. Legend `loc="lower left"`,
`bbox_to_anchor=(0.0716, 0.03)`, `ncol=2`. Eighteen annotations — copy verbatim.
`relabel` maps `NuOscProbExact -> "NuOscProbExact, $N_{\rm slabs}$"`,
`NuOscProbExact (tolerance) -> "NuOscProbExact, rtol"`,
`NuFast-Earth (dCP only) -> "NuFast-Earth ($\delta_{\rm CP}$ only)"`.

**NuFast-Earth appears twice on purpose.** It caches the work that depends on the
profile. A δ_CP scan invalidates none of it and costs 0.06 µs per probability
however finely the Earth is cut; moving Δm²₃₁ invalidates it and the same code
costs 4000 µs. The curve labelled plainly `NuFast-Earth` is the Δm²₃₁ scan; the
cheap one is labelled `NuFast-Earth ($\delta_{\rm CP}$ only)`. Same colour, different
marker — a hexagon, not a diamond, because sharing both colour and marker made them
hard to tell apart.

**Panel 3** (`tests/prem_speed_accuracy.json`, key `sterile_3plus1`): same chord
with a sterile state, 12 energies from 0.3 to 30 TeV, Δm²₄₁ = 1 eV², on
P(ν_μ→ν_μ) alone. Three codes hard-wire three flavours and cannot run it. Legend
`loc="lower left"`, `bbox_to_anchor=(0.012, 0.03)`, `ncol=1`. Restricted by
`only={"NuOscProbExact (double-double)", "NuOscProbExact (tolerance)", "nuSQuIDS",
"nuCraft"}` — the two root strategies are measured and frozen but only the default
is drawn, because the curves coincide to the last bit and plotting both put two
labels on one line.

### 6.5 `by_dial()` — do not skip this

`prem_axes` orders each series with `by_dial()`, which sorts by the **dial value**
parsed from `point["label"]`: descending if the smallest value is < 1 (a tolerance
runs loose to tight), ascending otherwise (a shell or slab count). Non-numeric
labels keep their input order.

Why it exists: NuFast-Earth's δ_CP curve spans 1.1× in cost (0.05745 to 0.06301 µs)
against 15× to 7880× for every other series, and consecutive points are separated
by less than their own standard deviations. Sorting that by cost sorts it by
measurement noise — it scrambled the layer counts into 16, 4, 64, 256, 8192, 4096,
1024, 1, so the line ran down to 8.0e-8 and back up to 1.8e-1, a U that reads as
two curves and is entirely an artifact of the ordering. By the dial it is what it
should be: a vertical line, cost pinned while the error falls six orders.

### 6.6 The leader lines

One `$N_{\rm layers} = 1$` label serves both NuFast-Earth curves, with a thin leader
to each. These are drawn **after** the axis limits are set, and are anchored to the
**label**, not to data coordinates: `savefig` re-lays the figure out with a tight
bounding box, which moves the axes under any data-space endpoint and leaves the
line detached from the text it is supposed to touch. `textcoords=_label` puts
`(0, 0.5)` at the middle of its left edge and `(1, 0.5)` at its right. Copy that
block verbatim; it is fiddly and the failure is silent.

### 6.7 Output

```python
fig.savefig(os.path.join(FIGDIR, "speed_accuracy_combined.pdf"))
```

`FIGDIR` comes from `os.environ.get('NUOSC_PAPER_FIGDIR', '.')`. In the paper the
figure is a `figure*` spanning both columns, placed `[t!]`.

---

## 7. Adding Magnus

### 7.1 Which panels Magnus belongs on

Decide per panel, on the same principle the existing figure uses — a code appears
where it can express the problem, and nowhere else:

- **Top (constant density):** yes, if Magnus has a constant-density path.
- **Middle (PREM, three flavours):** yes. This is Magnus's home ground.
- **Bottom (PREM, 3+1):** yes, if Magnus's sterile support reaches it. Magnus
  notebook 07 covers BSM sterile, so this is likely; confirm before promising it.

If Magnus cannot run a panel, say so in the caption the way the existing caption
does for nuCraft and NuFast-LBL. Do not fake a point.

### 7.2 What to measure

Follow `manifest.json` exactly. The three protocols:

- **ACCURACY** — untimed. Configure once, ask for the probabilities, difference
  against the 50-digit reference. No clock is read, so an accuracy number can never
  be mistaken for a speed. Normalisation: max absolute deviation over the grid.
- **AMORTIZED** — cost inside a fit or scan, with caching working as designed.
  Timed region is a scan of 25 steps over `dcp`, each step = configure(one varying
  parameter) + evaluate(whole grid). Everything invariant under the scan is hoisted
  out. Normalisation: µs per grid point = total / (steps × n_points). *This
  definition is the NuFast-Earth author's own harness
  (`src/examples/Speed.cpp::Atmospheric_Speed`), adopted deliberately so the
  definition is his rather than ours.*
- **THROUGHPUT** — what batching buys. Timed region is one request for N points.
  Codes that batch get one call; codes that cannot get a loop in their own language,
  labelled as a loop.

Grids, from `manifest.json -> grids`:

| Grid | Definition |
|---|---|
| `CONST/60E` | L = 1300 km, ρ = 3 g cm⁻³, 60 energies log-spaced 0.6–20 GeV |
| `CHORD/12x1` | cos θ_z = −0.9, 12 energies log 3–40 GeV; 3+1 variant 300–30000 GeV |
| `CONST/N-sweep` | L = 1300 km, ρ = 3 g cm⁻³, N swept 1…30000, channel `numu->nue` (throughput only) |
| `OSC/100x100` | 100 zenith × 100 energy, for the caching question |

Shared oscillation parameters, handed to every code and used to build every
reference — **normal ordering**:

```
s12sq = 0.31   s13sq = 0.0224   s23sq = 0.582
dcp   = 217°   dmsq21 = 7.39e-05 eV²   dmsq31 = 0.002525 eV²
3+1:  dmsq41 = 1.0 eV²
```

No adapter may contain a literal for any of these. `conversions.py` carries them
into C++ and `runner.py` into Python, so "the same values for all codes" is
structural rather than something to remember. Where a code wants a different
*parameterisation* of the same physics, derive it — never retype it. A hand-typed
`2.4511e-3` beside a `2.525e-3` is indistinguishable from a physics difference.

### 7.3 Magnus's dial

Every series on the Earth panels is a *dial swept over its range*, and the dial must
be named in the annotations. Pick Magnus's honestly: the Magnus expansion order,
the step count, a tolerance — whatever the user actually turns. Sweep it over the
range its interface offers, the way `manifest.json` sweeps shell counts from 1 to
65536 and tolerances from 1e-2 to 1e-12. Label both ends of the curve so no dial is
left unnamed.

### 7.4 The style slot for Magnus

Colours already taken: `C0` nuCraft, `C1` second-order expansion and the
eigensolver, `C2` nuSQuIDS, `C3` NuOscProbExact, `C4` both NuFast codes, `C5`
Prob3++, `C6` GLoBES. `"0.5"`/`"0.4"` are the double-precision rule and its label,
`"0.78"` the gridlines.

Markers already taken: `o v D * P s ^ h`.

**Recommendation: `("-X", "C9", 3.8)`** — cyan, unused, and `X` (filled x) is
unused and legible at 3.8 pt. Second choice `("-p", "C8", 4.0)` (olive pentagon).
**Check it visually before committing**: C9 against C2 green and C0 blue is the one
risk, and at these marker sizes it is a real one. If the Magnus paper wants Magnus
to read as the home code, plain black `"k"` is defensible and unambiguous — nothing
else on the plot is black except the axes and legend frames — but that is a
deliberate departure from "exactly as the one here" and should be the user's call,
not yours.

If Magnus needs two curves (e.g. two orders, or two thread counts), follow the
NuFast-Earth precedent: **same colour, different marker and linestyle**, and say in
the caption why one code appears twice.

### 7.5 Where Magnus's numbers go

Append to the copied JSON files, keeping the existing schema exactly:

- top panel → a new entry in `series` of the constant-density file, points carrying
  `label`, `knob`, `max_abs_error`, `us_per_probability`, `us_sd`, `block_cv`,
  `best_at_this_accuracy`.
- middle panel → a new entry in `series` of the Earth-plane file, points carrying
  `label`, `knob`, `max_abs_error`, `us_per_probability`, `us_sd`,
  `best_at_this_accuracy`.
- bottom panel → a new entry in `sterile_3plus1.series`, points carrying `label`,
  `max_abs_error`, `us_per_probability` and whatever dial field is Magnus's
  analogue of `n_slabs_per_segment` / `root_strategy`.

**Do not modify any existing series, value, or field.** The existing numbers are
frozen and are cited in the NuOscProbExact paper; a byte-level diff of the existing
series before and after your edit must be empty. Add a top-level provenance field
recording that Magnus's series were added, by what, and when.

**Register the new name in the style dict, and note that the two panels fail
differently if you forget.** The top panel loops `for name in STYLE`, so a series
absent from `STYLE` is **silently skipped** — the figure builds, looks right, and
is missing Magnus. The Earth panels do `order.index(x["name"])` and
`PREM_STYLE[series["name"]]`, so an unregistered name raises `ValueError` or
`KeyError` and you find out at once. After adding Magnus, count the legend entries
against the series you expect on every panel; do not trust that it drew.

---

## 8. Attribution — required

The existing curves are someone else's measurements redistributed into the Magnus
repo. Magnus already sets the precedent in
`notebooks/external_speed_accuracy.json`, which carries `_source`, `_attribution`,
`_note` and `_caveat`. Carry the same four fields on every file you copy:

- `_source`: `NuOscProbExact, <path>, commit 26a372c`
- `_attribution`: measured by the NuOscProbExact project on its own hardware
- `_note`: the problem statement, verbatim from the file's own `note`
- `_caveat`: the machine caveat from §1.1, unless you re-measured everything here

The Magnus paper's caption must say that the reference differs per code, that the
external timings come from the NuOscProbExact project, and which machine they were
taken on.

---

## 9. Acceptance checks

Run all of these. Each corresponds to a real failure.

1. `python3 -c "..."` from §3 prints `IDENTICAL` three times.
2. `pytest tests/test_bench_pipeline.py` — twelve invariants pass.
3. The three copied JSONs' existing `series` are byte-identical to the sources:
   `git diff --stat` on them shows only your added series.
4. The regenerated figure without Magnus is **pixel-identical** to
   NuOscProbExact's. Do this before adding Magnus; if it is not identical, the
   style or the data is wrong, and adding Magnus will hide it.

   ```bash
   # regenerate and re-execute — never edit a notebook by hand
   python3 notebooks/make_notebooks.py
   jupyter nbconvert --to notebook --execute --inplace notebooks/<the notebook>

   # compare against the original, page by page, at 150 dpi
   pdftoppm -r 150 -png A.pdf /tmp/a && pdftoppm -r 150 -png B.pdf /tmp/b
   compare -metric AE /tmp/a-1.png /tmp/b-1.png null: 2>&1   # expect 0
   ```
5. Every Magnus point has `block_cv <= 0.10`, or is recorded as inadmissible.
6. Every Magnus series records its actual thread count.
7. Magnus's accuracy is measured on `numu->numu`.
8. No point falls outside `XLIM_ALL`; if one does, widen the limit deliberately and
   note it, rather than letting it be clipped.
9. The legend has no two entries on one line, and no curve is hidden behind another.

---

## 10. Material that is stale or open — do not inherit it blindly

Handing you the tree without this would be handing you traps.

- **`tests/bench/README.md` is partly stale.** Its "State of play" section lists as
  "Not yet built" several things that have since been built — the accuracy path,
  the per-code reference builder, the orchestrator. The data files and 517
  artifacts exist and the figures are drawn from them. The README was last touched
  2026-08-19; `run_all.py` and `emit_figures.py` on 2026-09-01. Trust
  `manifest.json` and the code over that section.
- **`manifest.json` has one number that disagrees with what was measured.** It says
  nuCraft floors at 3e-3; the measurement recorded at 3+1 flavours was 3.7e-4.
  Neither has been re-run to settle it. Do not quote the manifest here.
- **`n_layers` cannot actually be swept** as the manifest promises: `Problem`
  carries a single `int knob` and `nufast_earth.cpp` hard-codes 256 layers.
- **A loose thread nobody has pulled:** `bench_nufast_earth`'s chord checksum sits
  about 0.05 per probability from `prob3` and `globes` at *every* knob value. A
  constant offset that survives refinement is a convention or channel mismatch, not
  a convergence error. It has not been understood.
- **The NuFast-LBL constant-density scan driver** was missing from the repository
  and has been reconstructed (`tests/external_drivers/nufast_scan.cpp`); it
  reproduces the frozen file to 4.4e-12 at `N_Newton=0` and 3.7e-12 at `N=2`. The
  metadata now points at the reconstruction.
- **Two claims that were reported during this work and are false.** Do not repeat
  either: (a) *"the wrong ρY_e was handed to NuFast-LBL"* — no, NuFast-LBL genuinely
  carries `1.52588e-4`, the same rounding Prob3++ uses, and the factor was right;
  (b) *"we drove NuFast-LBL one energy per call"* — no, the committed driver passes
  a vector of up to 30000 energies. The paper's prose was wrong, not the
  measurement.
- **nuCraft's constants are wrong, and this is deliberately not absorbed.** Its
  charged-current and sterile entries are `0.00015256` and `7.6525e-05`, ratio
  0.50161 where the isoscalar value is exactly 1/2. By nuCraft's own derivation
  (`NuCraft.py:205–212`) the two must be exactly 2:1, so the pair is inconsistent
  with its own documentation — most likely a digit typo for `7.6325e-5`. Its
  reference therefore uses nuCraft's constants for the units but the **correct**
  ratio: a reference that absorbed the rounding would forgive the defect and the
  floor it causes would silently vanish. **This bites only at 3+1**; at three
  flavours the sterile entry never enters.
- **Prob3++ does not reach round-off in vacuum** (5.7e-4). It is not a convention
  error: its ħc is hard-coded as `2.534` in `mosc.c`, four significant figures,
  implying 1.9731650e-7 against everyone else's 1.97327e-7. That is why it flattens
  near 3e-4 no matter how fine the shells.

---

## 11. What not to do

- Do not re-measure or re-scale any existing series.
- Do not force every code to one thread to "make it fair". §5 explains why.
- Do not build one shared reference for all codes. §1.2 and §10 explain why.
- Do not tighten `XLIM_ALL` to make room.
- Do not edit the notebooks directly — they are generated from
  `make_notebooks.py`, and a bare notebook fails CI. Edit the generator, then
  regenerate **and** re-execute.
- Do not quote a version as a date. Every external code is pinned by tag,
  commit or tarball hash; keep it that way.

---

## 12. Source-side provenance, for the record

- Repo: `/home/mbustamante/Research/NuOscProb/NuOscProbExact`
- Branch `dev-fair-benchmarks`, commit `26a372c`, 2026-09-01
- `notebooks/make_notebooks.py` is byte-identical between `b9c1bfd` and `26a372c`;
  the intervening commits touch prose and the bibliography only
- Figure 11 = `\label{fig:planes}` = `speed_accuracy_combined.pdf`
- Code pins: NuFast-LBL `v2.0.1` `a6eec95b3284`; NuFast-Earth `v1.2.0`
  `499392d49f21`; nuSQuIDS `v1.13.3` `104914da5a25`; nuCraft `r22` sha256
  `f406750025f8`; GLoBES `3.2.18` sha256 `7edd0fea6780`; Prob3++ `v3r20`
  `fa103b1e57fb`
- Two build profiles: **speed** uses each project's own upstream flags, so no code
  is published slower than its authors publish it; **accuracy** uses
  `-O3 -std=c++17` for everyone, no `-Ofast`, no `-ffast-math`, because
  value-changing flags make an accuracy measurement measure the flag. One exception,
  recorded rather than hidden: nuSQuIDS is measured as distributed, from its
  upstream wheel, so its flags are the wheel's.
