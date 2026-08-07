#!/usr/bin/env python3
r"""Builds the two derived versions of the paper from ``main.tex``.

``main.tex`` is the working copy: it carries the revision in its own
mark-up, ``\add{...}`` for new material and ``\del{...}`` for removed, with
``\addstart``/``\addend`` for spans too large to pass as a macro argument.
Setting ``\markupfalse`` in the preamble prints it without colour, but the
macros are still there and the deleted text is still in the file.

This produces two things from it:

``main_clean.tex``
    the mark-up *resolved* rather than merely uncoloured --- every
    ``\add{X}`` replaced by ``X``, every ``\del{X}`` dropped, and the span
    switches removed.  This is the source to hand to a journal or to
    arXiv, and the one to paste into Overleaf if you want the paper
    without the revision showing.

``main_diff.tex``
    ``latexdiff`` of ``baseline_cpc_v1.tex`` --- the published version ---
    against ``main_clean.tex``, which marks the revision mechanically
    rather than by hand.  It is worth having both: the hand mark-up says
    what the author meant to change, and the diff says what actually
    changed.  Where they disagree, one of them is wrong.

Run it from this directory::

    python make_versions.py

``main_clean.tex`` needs nothing but Python.  ``main_diff.tex`` needs
``latexdiff`` on the PATH, or its location in ``$LATEXDIFF``; it is in
TeX Live and in Debian's ``latexdiff`` package.

Known limitation: ``main_diff.tex`` is generated but does not currently
compile.  ``latexdiff`` splits inline-math spans in the appendices across
its own ``\DIFadd{}`` boundaries, which breaks the surrounding ``$...$``.
``--type=CFONT`` takes it from 147 LaTeX errors to 8, and the remainder
move around the appendices rather than clearing; ``--math-markup=off``,
``--math-markup=whole`` and forcing the array and equation environments
through ``PICTUREENV`` do not fix it either.  The eight need hand-patching
in the generated file, or a different tool.
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


def _balanced(text, start):
    r"""Returns the index just past the brace group opening at `start`."""
    depth = 0
    i = start
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
    raise ValueError('unbalanced braces from position %d' % start)


def resolve(text, macro, keep):
    r"""Resolves every ``\macro{...}``, keeping the argument or dropping it.

    Brace counting rather than a regular expression, because the
    arguments contain braces of their own --- ``\add{$\mathbb{H}$}`` and
    worse.
    """
    out = []
    i = 0
    token = '\\' + macro + '{'
    while True:
        j = text.find(token, i)
        if j < 0:
            out.append(text[i:])
            return ''.join(out)
        out.append(text[i:j])
        end = _balanced(text, j + len(token) - 1)
        if keep:
            out.append(text[j+len(token):end-1])
        i = end


def clean(text):
    r"""The paper as it should read after the revision.

    Iterated to a fixed point, because the mark-up nests: a `\\add` span
    that contains a figure contains that figure's own `\\add` caption, and
    resolving the outer one only exposes the inner.
    """
    # \protect guards the macro, not the text: once the macro is gone the
    # \protect would run into the word after it, which is how the title
    # became \protectfast.
    text = text.replace('\\protect\\add{', '\\add{')
    text = text.replace('\\protect\\del{', '\\del{')
    for macro, keep in (('add', True), ('del', False)):
        while ('\\' + macro + '{') in text:
            before = text
            text = resolve(text, macro, keep=keep)
            if text == before:
                raise ValueError('resolving \\%s made no progress' % macro)
    # the span switches, where they are *used*
    text = re.sub(r'^[ \t]*\\add(start|end)[ \t]*$\n?', '', text, flags=re.M)
    text = re.sub(r'\\add(start|end)\b', '', text)
    return text


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


def detexorpdfstring(text):
    r"""Defines ``\DIFadd``/``\DIFdel`` without ``\texorpdfstring``.

    ``latexdiff`` notices ``hyperref`` and routes its two mark-up macros
    through ``\texorpdfstring`` so that a bookmark gets the plain text.
    Inside a ``\caption``, inside a ``minipage``, inside a ``figure*`` ---
    which is where the two speed-accuracy planes now live --- that recurses
    until TeX's input stack is exhausted.  The bookmarks are not worth it:
    the macros are defined to their own ``tex`` variants instead, which is
    what ``latexdiff`` uses when ``hyperref`` is absent.
    """
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
    source = open(os.path.join(HERE, 'main.tex')).read()

    # Only the body is resolved.  The preamble defines \\add and friends and
    # must keep them: cleaning the whole file turns
    # \\newcommand{\\addstart}{...} into \\newcommand{}{...}, which fails
    # in a way that takes a while to read back to its cause.
    split = source.index('\\begin{document}')
    body = clean(source[split:])
    resolved = source[:split] + body
    # the body only: the preamble keeps the definitions, and a comment
    # there explains the convention using the macros it documents
    for macro in (r'\add{', r'\del{', r'\addstart', r'\addend'):
        assert macro not in body, 'left behind: ' + macro
    open(os.path.join(HERE, 'main_clean.tex'), 'w').write(resolved)
    print('main_clean.tex written (%d bytes)' % len(resolved))

    if not LATEXDIFF or not os.path.exists(LATEXDIFF):
        print('latexdiff not on the PATH and $LATEXDIFF unset; '
              'main_clean.tex is written, the diff is skipped')
        compile_twice('main_clean')
        return 0
    old_bib, new_bib = [], []
    tmp_old = os.path.join(HERE, '.diff_old.tex')
    tmp_new = os.path.join(HERE, '.diff_new.tex')
    open(tmp_old, 'w').write(
        hide_bibliography(open(BASELINE).read(), old_bib))
    open(tmp_new, 'w').write(hide_bibliography(resolved, new_bib))
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
    diffed = repair(detexorpdfstring(
        restore_bibliography(r.stdout, new_bib)))
    open(os.path.join(HERE, 'main_diff.tex'), 'w').write(diffed)
    print('main_diff.tex written (%d bytes)' % len(diffed))

    compile_twice('main_clean')
    compile_twice('main_diff')
    return 0


if __name__ == '__main__':
    sys.exit(main())
