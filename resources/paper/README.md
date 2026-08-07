# The paper

The source of *{\tt NuOscProbExact}: a fast, general-purpose code to compute
exact two-, three-, and four-flavor neutrino oscillation probabilities*, the
Computer Physics Communications article that documents this library.

It lives here so that the paper and the code it describes travel together: the
figures in it are produced by `notebooks/10_paper_figures.ipynb` from data
frozen in `tests/`, and a claim in the text can be traced to the measurement
behind it without leaving the repository.

Nothing here is installed. The wheel is built from an explicit list of ten
modules under `src/`, and the source distribution from setuptools' default set
plus `MANIFEST.in`; neither reaches `resources/`. It is here for a clone, not
for `pip`.

## Building it

```bash
cd resources/paper
pdflatex main.tex && pdflatex main.tex
```

Twice, for the cross-references. There is no `bibtex` step: the bibliography
is an inlined `thebibliography` environment and `\bibliography{refs.bib}` is
commented out, so `refs.bib` is kept for reference but is *not* read at
compile time. A new citation goes in by hand at the end of the environment,
and the numbering is positional.

Everything the preamble loads is in a normal TeX Live install. `elsarticle.cls`
and `elsarticle-num.bst` are bundled anyway, so that the folder compiles on a
machine whose TeX Live was installed without the Elsevier bundle.

## The three versions

`main.tex` is the working copy, and it carries the revision in its own mark-up:

| | |
|---|---|
| `\add{...}` | added in this revision — green |
| `\del{...}` | removed in this revision — red, struck through |
| `\addstart` … `\addend` | the same, for spans too large to pass as a macro argument: a whole section, or anything containing a float, an equation or a listing |

`\markuptrue` in the preamble shows the colours; `\markupfalse` prints the same
file without them. Both states must compile before anything is handed on.

Floats do **not** inherit `\addstart`, so a table or figure inside such a span
needs its own `\addstart` after `\begin{table}`.

`python make_versions.py` derives two more files from it:

- **`main_clean.tex`** — the mark-up *resolved*, not merely uncoloured: every
  `\add{X}` becomes `X`, every `\del{X}` disappears, and the span switches are
  removed. This is the file for a journal, for arXiv, or for Overleaf. It needs
  nothing but Python.
- **`main_diff.tex`** — `latexdiff` of `baseline_cpc_v1.tex`, the published
  version, against `main_clean.tex`. The hand mark-up says what the author
  meant to change; the diff says what actually changed. Where the two disagree,
  one of them is wrong, which is the point of having both. Needs `latexdiff` on
  the `PATH`, or its location in `$LATEXDIFF`.

**`main_diff.tex` does not currently compile.** It is generated, but
`latexdiff` splits inline-math spans in the appendices across its own
`\DIFadd{}` boundaries and breaks the surrounding `$...$`. `--type=CFONT` takes
it from 147 LaTeX errors down to 8; `--math-markup=off`, `--math-markup=whole`,
and forcing the array and equation environments through `PICTUREENV` do not
clear the rest, which move around the appendices rather than disappearing. The
eight need hand-patching in the generated file, or a different tool. The clean
version is unaffected.

## The figures

All fifteen are in `figs/`, found through a single `\graphicspath{{figs/}}` so
that no `\includegraphics` names the folder. Twelve are produced by
`notebooks/10_paper_figures.ipynb`; `Notes_on_SU_n__probability_relations.pdf`
is external, and `prob_3nu_vs_energy_compare.jpg` is an old raster that nothing
includes — the `.pdf` beside it is what the paper uses.

To rebuild them:

```bash
cd notebooks
NUOSC_PAPER_FIGDIR=../resources/paper/figs jupyter nbconvert --to notebook \
    --execute --inplace --ExecutePreprocessor.timeout=2800 10_paper_figures.ipynb
```

Three to five minutes. Note that this rewrites all twelve PDFs whether or not
anything changed, and the bytes differ between runs, so `git status` will show
them as modified afterwards even when the figures are identical. That is
expected; commit them only when a figure actually changed.

The data behind the comparison figures is frozen in `tests/`, so none of the
external codes has to be installed to redraw them:

| file | what it holds |
|---|---|
| `tests/speed_accuracy.json` | the six-code constant-density plane |
| `tests/prem_speed_accuracy.json` | the two Earth planes, three-flavor and 3+1 |
| `tests/const_density_scan.json` | the 50-digit reference for the comparison figure |
| `tests/timing_other_codes.json` | the timings behind the performance figure |
| `tests/nusquids_scan.json`, `tests/nufast_scan.json` | those two codes against energy |

`tests/external_drivers/` holds the C and C++ drivers for the codes that cannot
be called from Python, with the raw output of each and a README recording every
convention that had to be matched first. The previous round of this comparison
built its drivers in a scratch directory and lost them, which is why they are
in the repository now.

## Timings

Every timing quoted in the paper was measured on one laptop: an Intel Core
i5-1334U, ten cores and twelve threads to 4.6 GHz, 16 GB of memory, Ubuntu
24.04 LTS, Python 3.12.7, numpy 1.26.4, numba 0.60.0, gcc 13.3. Absolute
numbers will differ elsewhere; the ratios between codes, which is what the
figures are about, are far more stable.
