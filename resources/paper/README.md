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

## The two versions

`main.tex` is the paper as it should read. **Write ordinary LaTeX in it.**
Nothing in it records what changed since the published version; that is worked
out mechanically:

```bash
cd resources/paper
python make_versions.py
```

which writes **`main_diff.tex`** — `latexdiff` of `baseline_cpc_v1.tex`, the
published version, against `main.tex` — and compiles both. Additions come out
blue, deletions red and struck through, both in the document's own font.
`main_diff.tex` is the file to send a journal alongside `main.tex` itself.
Needs `latexdiff` on the `PATH`, or its location in `$LATEXDIFF`.

Both compile: `main.pdf` at 28 pages, `main_diff.pdf` at 29, neither with
errors.

### The mark-up convention is retired

`main.tex` used to carry the revision by hand — `\add{...}` for new material,
`\del{...}` for removed, `\addstart`/`\addend` for spans too large to pass as a
macro argument — with `make_versions.py` resolving it into a `main_clean.tex`.
That is gone. Two files to keep in agreement became one, and the diff no longer
depends on the author having remembered to mark a change while typing it.

`make_versions.py` refuses to run if any of those four macros reappears in
`main.tex`, since they are no longer defined in the preamble and would compile
to an undefined-command error further downstream. If you paste text from an old
draft and hit that message, delete the macro and keep the words.

`ulem` is still loaded, and is the one leftover: `main_diff.tex` inherits this
preamble, and `latexdiff`'s deletions are struck through with `\sout`.

### Line numbers

On, for a referee. The preamble carries a switch:

```latex
\linenostrue         % \linenosfalse for the camera-ready
```

`elsarticle`'s own `review` option will *not* do this — it only sets preprint
mode and 1.5 spacing, and discards the `5p` two-column layout — so `lineno` is
loaded directly. `switch` puts the numbers on the outer edge of each column,
which two-column needs: without it the right column's numbers land in the
gutter. `mathlines` numbers display equations. Floats are never numbered, so
the figures, tables and listings are skipped, and `lineno` is known to drop
numbers on some `amsmath` display environments.

The diff marks changes the way `main.tex` does — additions blue, deletions
red and struck through, both in the document's own font and size.
`latexdiff`'s `CFONT` type instead sets additions in sans-serif and deletions
at `\scriptsize` without striking them, so `make_versions.py` redefines its
two macros; `ulem` is already loaded with `normalem`.

Five things had to be worked around to get `latexdiff` through this document,
and `make_versions.py` does all four. They are recorded here because each one
failed in a way that did not point at its own cause:

- **The bibliography is held out of the diff entirely.** `latexdiff` marks up
  the arguments of `\href` and `\path` inside `\bibitem`, and a marked-up URL
  sends `hyperref` into a recursion that exhausts TeX's input stack — reported
  after several minutes as a capacity error, at whatever line it happened to
  be reading. Note that this document has *two* `thebibliography`
  environments, a short one in the front matter and the real one; holding out
  only the first leaves the real one exposed, which is what took longest to
  see.
- **`\texorpdfstring` comes out of the two mark-up macros.** `latexdiff`
  routes them through it when `hyperref` is present, so bookmarks get plain
  text; inside a `\caption` inside a `minipage` inside a `figure*` — which is
  where the two speed-accuracy planes live — that recurses.
- **`--graphics-markup=none`**, since the default redefines `\includegraphics`
  through `\LetLtxMacro` to box changed figures, and that recurses here too.
- **`--disable-citation-markup`**, which is cosmetic rather than fatal: it
  otherwise rewrites every `\cite` into an `\hspace{0pt}%DIFAUXCMD`
  construction that is noise inside a float caption.

`repair()` also fixes two things `latexdiff` emits on a revision that moved
whole environments around — a `lstlisting` swallowed by a `\DIFadd{}`, and an
environment whose `\begin` and `\end` land in different added blocks. Both
report zero on the current source; they are kept because they cost nothing and
the next revision may reintroduce either.

A consequence worth knowing: because the bibliography is held out, **new
references are not marked as added in the diff**. They show where they are
cited, which is the useful place.

## The figures

All fifteen are in `figs/`, found through a single `\graphicspath{{figs/}}` so
that no `\includegraphics` names the folder. Fourteen are produced by
`notebooks/10_paper_figures.ipynb`; `Notes_on_SU_n__probability_relations.pdf`
is the one exception, and is external.

To rebuild them:

```bash
cd notebooks
NUOSC_PAPER_FIGDIR=../resources/paper/figs jupyter nbconvert --to notebook \
    --execute --inplace --ExecutePreprocessor.timeout=2800 10_paper_figures.ipynb
```

Three to five minutes. Note that this rewrites all fourteen PDFs whether or not
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
24.04.4 LTS, Python 3.12.7, numpy 1.26.4, numba 0.60.0, gcc 13.3. Absolute
numbers will differ elsewhere; the ratios between codes, which is what the
figures are about, are far more stable.
