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
    r = run(['perl', LATEXDIFF, '--type=CFONT', '--math-markup=off',
             '--flatten', BASELINE, 'main_clean.tex'])
    if r.returncode != 0:
        print('latexdiff failed:', r.stderr[-400:])
        return 1
    open(os.path.join(HERE, 'main_diff.tex'), 'w').write(r.stdout)
    print('main_diff.tex written (%d bytes)' % len(r.stdout))

    compile_twice('main_clean')
    compile_twice('main_diff')
    return 0


if __name__ == '__main__':
    sys.exit(main())
