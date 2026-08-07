#!/usr/bin/env python3
r"""Builds ``main_diff.tex``, the marked-up revision, from ``main.tex``.

``main.tex`` is the paper as it should read.  It carries no revision
mark-up of its own: it is written as ordinary LaTeX, and what changed is
worked out mechanically here rather than recorded by hand while typing.

``main_diff.tex``
    ``latexdiff`` of ``baseline_cpc_v1.tex`` --- the published version ---
    against ``main.tex``.  Additions come out blue, deletions red and
    struck through, both in the document's own font.  This is the file to
    send a journal alongside ``main.tex`` itself.

Run it from this directory::

    python make_versions.py

It needs ``latexdiff`` on the PATH, or its location in ``$LATEXDIFF``; it
is in TeX Live and in Debian's ``latexdiff`` package.  Both files are then
compiled twice, and their page and error counts printed.

An earlier version of this script also derived a ``main_clean.tex``, by
resolving ``\add{...}`` and ``\del{...}`` mark-up that ``main.tex`` used to
carry.  That convention is retired --- ``main.tex`` *is* the clean version
now --- and with it goes the job of keeping a hand-written account of the
revision in agreement with the real one.
"""

import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BASELINE = os.path.join(HERE, 'baseline_cpc_v1.tex')
# On the PATH by default; $LATEXDIFF overrides, for a copy unpacked
# somewhere of its own.
LATEXDIFF = os.environ.get('LATEXDIFF') or shutil.which('latexdiff')


BIB = re.compile(r'\\begin\{thebibliography\}.*?\\end\{thebibliography\}',
                 re.S)


def hide_bibliography(text, store):
    r"""Swaps the bibliography for a placeholder before diffing.

    ``latexdiff`` marks up the arguments of ``\href`` and ``\path`` inside
    ``\bibitem``, and a marked-up URL sends ``hyperref`` into a recursion
    that exhausts TeX's input stack rather than failing cleanly ---
    ``--append-safecmd`` does not reach it, and the compile simply hangs.

    A diff does not need a coloured bibliography anyway: new references
    show up where they are cited.  So the whole environment is held out of
    the comparison and put back afterwards, from the *new* version, which
    means the reference list in the diff is the current one and is not
    marked as changed.
    """
    def swap(m):
        store.append(m.group(0))
        return '\\BIBPLACEHOLDER{%d}' % (len(store)-1)
    # Every one of them: this document has two, a short list in the
    # front matter and the real one, and replacing only the first left
    # the real bibliography exposed to the diff.
    return BIB.sub(swap, text)


def restore_bibliography(text, bibliography):
    r"""Puts it back, stripping any mark-up latexdiff wrapped around it."""
    text = re.sub(r'\\DIFaddbegin\s*(\\BIBPLACEHOLDER\{\d+\})\s*\\DIFaddend',
                  lambda m: m.group(1), text)
    text = re.sub(r'\\DIF(?:add|del)\{\s*(\\BIBPLACEHOLDER\{\d+\})\s*\}',
                  lambda m: m.group(1), text)
    text = re.sub(r'%DIFDELCMD < *\\BIBPLACEHOLDER\{\d+\}[^\n]*\n', '', text)
    return re.sub(r'\\BIBPLACEHOLDER\{(\d+)\}',
                  lambda m: bibliography[int(m.group(1))], text)


ENVS = ('equation*', 'equation', 'eqnarray*', 'eqnarray')


def _group_end(text, i):
    r"""Index just past the brace group that opens at ``text[i] == '{'``."""
    depth = 0
    while i < len(text):
        if text[i] == '\\':
            i += 2
            continue
        if text[i] == '{':
            depth += 1
        elif text[i] == '}':
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1
    return -1


def _blocks(text):
    r"""Every ``\DIFaddbegin`` ... ``\DIFaddend`` span, as index pairs."""
    out = []
    for m in re.finditer(r'\\DIFaddbegin\b', text):
        stop = text.find('\\DIFaddend', m.end())
        if stop >= 0:
            out.append((m.end(), stop))
    return out


def restyle_markup(text):
    r"""Rewrites the four macros ``latexdiff`` marks changes with.

    Two changes, for two different reasons.

    *Appearance.*  ``--type=CFONT`` sets additions in blue sans-serif and
    deletions in red at ``\scriptsize``, so changed text is a different
    typeface and size from the text around it and deletions are not struck
    at all.  They are set in the document's own font instead: additions
    blue, deletions red and struck through with ``\sout``, which is what
    ``main.tex``'s own ``\add`` and ``\del`` do and what a reader of the
    revision expects.  ``ulem`` is already loaded, with ``normalem``.

    *Recursion.*  ``latexdiff`` notices ``hyperref`` and routes its two
    mark-up macros through ``\texorpdfstring`` so a bookmark gets the plain
    text.  Inside a ``\caption``, inside a ``minipage``, inside a
    ``figure*`` --- which is where the two speed-accuracy planes live ---
    that recurses until TeX's input stack is exhausted.  The bookmarks are
    not worth it.
    """
    text = text.replace(
        r'\providecommand{\DIFaddtex}[1]{{\protect\color{blue} \sf #1}}',
        r'\providecommand{\DIFaddtex}[1]{{\protect\color{blue}#1}}')
    text = text.replace(
        r'\providecommand{\DIFdeltex}[1]{{\protect\color{red} \scriptsize #1}}',
        r'\providecommand{\DIFdeltex}[1]{{\protect\color{red}\sout{#1}}}')
    text = text.replace(
        r'\providecommand{\DIFadd}[1]{\texorpdfstring{\DIFaddtex{#1}}{#1}}',
        r'\providecommand{\DIFadd}[1]{\DIFaddtex{#1}}')
    text = text.replace(
        r'\providecommand{\DIFdel}[1]{\texorpdfstring{\DIFdeltex{#1}}{}}',
        r'\providecommand{\DIFdel}[1]{\DIFdeltex{#1}}')
    return text


def repair(text):
    r"""Undoes the two ways ``latexdiff`` breaks this particular document.

    Neither is anything the author did; both are what the tool emits on a
    revision that moved whole environments around, and both are mechanical.

    *A listing inside a macro argument.*  ``latexdiff`` wraps added text in
    ``\DIFadd{...}``, and one span swallowed a ``lstlisting``.  That cannot
    work: the environment reads its body verbatim, and a macro argument is
    tokenised before it ever runs --- the same constraint that makes
    ``main.tex`` mark listing changes in the caption rather than the body.
    The wrapper comes off; the surrounding block still colours the listing.

    *An environment split across two added blocks.*  A ``\begin{equation*}``
    landed in one ``\DIFaddbegin`` block and its ``\end`` in another, which
    leaves the first block with an unclosed environment and a stray ``}``,
    and the second with an ``\end`` that opens nothing.  Both halves are
    repaired: the missing ``\end`` is put back, the stray brace dropped, and
    the orphaned ``\end`` removed.
    """
    out, i, unwrapped = [], 0, 0
    while True:
        m = re.compile(r'\\DIFadd\{').search(text, i)
        if not m:
            out.append(text[i:])
            break
        end = _group_end(text, m.end()-1)
        if end < 0 or '\\begin{lstlisting}' not in text[m.end():end-1]:
            out.append(text[i:m.end()])
            i = m.end()
            continue
        out.append(text[i:m.start()])
        out.append(text[m.end():end-1])
        unwrapped += 1
        i = end
    text = ''.join(out)

    closed = orphaned = 0
    for env in ENVS:
        begin, end_ = '\\begin{%s}' % env, '\\end{%s}' % env
        while True:                                   # missing \end
            for start, stop in _blocks(text):
                block = text[start:stop]
                if block.count(begin) <= block.count(end_):
                    continue
                head = text[:stop].rstrip()
                if head.endswith('}}'):
                    head = head[:-1]
                text = head + end_ + '\n' + text[stop:]
                closed += 1
                break
            else:
                break
        while True:                                   # orphaned \end
            for start, stop in _blocks(text):
                block = text[start:stop]
                if block.count(end_) <= block.count(begin):
                    continue
                k = text.index(end_, start)
                text = text[:k] + text[k+len(end_):]
                orphaned += 1
                break
            else:
                break

    print('  repaired: %d listing wrapper(s), %d unclosed and %d orphaned '
          'environment(s)' % (unwrapped, closed, orphaned))
    return text


def run(cmd, **kw):
    return subprocess.run(cmd, cwd=HERE, check=False,
                          capture_output=True, text=True, **kw)


def compile_twice(stem):
    for _ in range(2):
        r = run(['pdflatex', '-interaction=nonstopmode', stem + '.tex'])
    errors = [l for l in r.stdout.split('\n') if l.startswith('! ')]
    print('  %-16s %d pages, %d errors'
          % (stem + '.pdf', pages(stem + '.pdf'), len(errors)))
    for line in errors[:3]:
        print('     ', line)


def pages(path):
    r = run(['pdfinfo', path])
    for line in r.stdout.split('\n'):
        if line.startswith('Pages:'):
            return int(line.split()[1])
    return -1


def main():
    if not LATEXDIFF or not os.path.exists(LATEXDIFF):
        print('latexdiff is not on the PATH and $LATEXDIFF is unset; '
              'there is nothing to build')
        return 1
    source = open(os.path.join(HERE, 'main.tex')).read()
    # main.tex is the clean version.  If one of the retired macros creeps
    # back in -- pasted from the old source, most likely -- say so, rather
    # than emitting a diff with an undefined command inside it.
    for macro in (r'\add{', r'\del{', r'\addstart', r'\addend'):
        if macro in source:
            print('main.tex contains %s, but the mark-up convention was '
                  'retired; write plain LaTeX and let the diff find it'
                  % macro)
            return 1

    old_bib, new_bib = [], []
    tmp_old = os.path.join(HERE, '.diff_old.tex')
    tmp_new = os.path.join(HERE, '.diff_new.tex')
    open(tmp_old, 'w').write(
        hide_bibliography(open(BASELINE).read(), old_bib))
    open(tmp_new, 'w').write(hide_bibliography(source, new_bib))
    # --graphics-markup=none: latexdiff's default highlighting redefines
    # \includegraphics through \LetLtxMacro so it can box changed figures.
    # In this document that redefinition recurses, and TeX reports it as an
    # exhausted input stack -- after several minutes, at an
    # \includegraphics line, rather than as a syntax error, which is what
    # made it hard to place.
    #
    # --disable-citation-markup: latexdiff otherwise rewrites every \cite
    # into an \hspace{0pt}%DIFAUXCMD construction, which is noise inside a
    # float caption and buys nothing here.
    r = run(['perl', LATEXDIFF, '--type=CFONT', '--math-markup=off',
             '--disable-citation-markup', '--graphics-markup=none',
             '--flatten', '.diff_old.tex', '.diff_new.tex'])
    for tmp in (tmp_old, tmp_new):
        os.remove(tmp)
    if r.returncode != 0:
        print('latexdiff failed:', r.stderr[-400:])
        return 1
    diffed = repair(restyle_markup(
        restore_bibliography(r.stdout, new_bib)))
    open(os.path.join(HERE, 'main_diff.tex'), 'w').write(diffed)
    print('main_diff.tex written (%d bytes)' % len(diffed))

    compile_twice('main')
    compile_twice('main_diff')
    return 0


if __name__ == '__main__':
    sys.exit(main())
