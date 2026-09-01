# -*- coding: utf-8 -*-
r"""Paragraph length and sentence count, measured from the LaTeX source.

Why this exists.  The audit detector reads `pdftotext -layout` output and finds
paragraph boundaries from indentation.  A displayed equation interrupts the
indentation, so one source paragraph containing three displays renders as many
indented blocks, and criteria 12 (paragraphs over 220 words) and 16a
(paragraphs over six sentences) never see it.  On 2026-08-25 that turned out to
have hidden a 523-word, 26-sentence paragraph in Sec. 3 and a 362-word,
13-sentence one in Sec. 5, among 31 paper-wide, through every audit pass run so
far -- all of which reported those two criteria clean.

In main.tex a paragraph is a blank line to a blank line.  Prose following a
display with no blank line between is still the same paragraph.  That is what
this measures.

    python3 code/check_paragraphs.py                 # the whole paper
    python3 code/check_paragraphs.py "Finding the"   # one section, all paragraphs
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TEX = os.path.join(HERE, os.pardir, os.pardir, "paper", "main.tex")

MAX_WORDS = 220        # criterion 12
MAX_SENTENCES = 6      # criterion 16a
BAND = (90, 150)       # criterion 16a, his own papers
MIN_REPORT = 30        # ignore fragments


def sections(text):
    out = []
    for m in re.finditer(r"\\section\*?\{", text):
        depth = 0
        for k in range(m.end() - 1, len(text)):
            if text[k] == "{":
                depth += 1
            elif text[k] == "}":
                depth -= 1
                if depth == 0:
                    break
        out.append((m.start(), " ".join(text[m.end():k].split())))
    return out


def section_at(bounds, pos):
    n, name = 0, "(front)"
    for i, (p, nm) in enumerate(bounds):
        if p <= pos:
            n, name = i + 1, nm
    return n, name


def strip_floats(text):
    """Drop float bodies, keeping their source offsets straight."""
    pieces, pos = [], 0
    for m in re.finditer(r"\\begin\{(figure|table|sidewaysfigure|sidewaystable|wrapfigure)\*?\}.*?\\end\{\1\*?\}", text, re.S):
        pieces.append((pos, text[pos:m.start()]))
        pos = m.end()
    pieces.append((pos, text[pos:]))
    return pieces


ABBREV = (r"(?<!\bet\sal)(?<!\be\.g)(?<!\bi\.e)(?<!\bcf)(?<!\bvs)"
          r"(?<!\bFig)(?<!\bEq)(?<!\bEqs)(?<!\bSec)(?<!\bSecs)(?<!\bRef)"
          r"(?<!\bRefs)(?<!\bApp)(?<!\bTab)(?<!\bNo)(?<!\bcirca)")


def sentences(text):
    """Split into sentences without breaking at abbreviations.

    "Barger et al. is the subject of Sec. 4" is one sentence, not three; before
    this guard the counts were inflated wherever a paragraph cited anyone.
    """
    parts = re.split(ABBREV + r"(?<=[.?!])\s+(?=[A-Z(])", text)
    return [x.strip() for x in parts if x.strip()]


def clean(par):
    t = re.sub(r"\\begin\{(equation|align|aligned|gather)\*?\}.*?\\end\{\1\*?\}",
               " EQ ", par, flags=re.S)
    t = re.sub(r"\\(?:eq)?ref\{[^}]*\}", "REF", t)
    t = re.sub(r"\\cite\w*(\[[^\]]*\])?\{[^}]*\}", "C", t)
    t = re.sub(r"\$[^$]*\$", "X", t)
    t = re.sub(r"\\(section|subsection|subsubsection|paragraph)\*?\{[^}]*\}", " ", t)
    # Keep the ARGUMENT of inline markup.  Deleting it wholesale swallowed the
    # "By sector:" / "By code:" labels of Sec. 11 and left sentences opening on a
    # bare colon, which made the sentence count unreadable.
    for _ in range(3):
        t = re.sub(r"\\(emph|textit|textbf|texttt|mbox|text)\{([^{}]*)\}", r"\2", t)
    t = re.sub(r"\\\w+\*?(\[[^\]]*\])?(\{[^}]*\})?", " ", t)
    return " ".join(t.split())


def split_lists(par):
    r"""Yield the pieces of a paragraph, treating a list as structure.

    An enumerate or itemize carries no blank lines, so a seven-item list reads
    as one 523-word paragraph unless the items are separated here.  A long
    ITEM is still worth flagging; a list of ordinary items is not a long
    paragraph.  (Found on 2026-08-25: the "seven evaluations" list of Sec. 3
    was the largest false positive this script produced.)
    """
    if not re.search(r"\\begin\{(enumerate|itemize|description)\}", par):
        yield par
        return
    body = re.sub(r"\\(begin|end)\{(enumerate|itemize|description)\}", " ", par)
    for piece in re.split(r"\\item\b", body):
        if piece.strip():
            yield piece


def paragraphs(text):
    # Section offsets must be measured on the same string the paragraphs are cut
    # from, or every paragraph is attributed to the wrong section.
    start = text.index(r"\section{")
    body = text[start:]
    if r"\appendix" in body:
        body = body[:body.index(r"\appendix")]
    bounds = sections(body)
    for base, chunk in strip_floats(body):
        off = 0
        for par in re.split(r"\n\s*\n", chunk):
            start = base + off
            off += len(par) + 2
            for piece in split_lists(par):
                c = clean(piece)
                if len(c.split()) < MIN_REPORT:
                    continue
                sents = sentences(c)
                yield section_at(bounds, start), len(c.split()), len(sents), c


def main():
    want = sys.argv[1] if len(sys.argv) > 1 else None
    text = open(TEX, encoding="utf-8").read()
    rows, over = [], []
    for (num, name), w, ns, c in paragraphs(text):
        if want and want.lower() not in name.lower():
            continue
        rows.append((num, name, w, ns, c))
        if w > MAX_WORDS or ns > MAX_SENTENCES:
            over.append((w, ns, num, name, c))
    if want:
        print("%-6s %-5s %s" % ("words", "sent", "opens"))
        for num, name, w, ns, c in rows:
            flag = "  <-- OVER" if (w > MAX_WORDS or ns > MAX_SENTENCES) else ""
            print("%-6d %-5d %s%s" % (w, ns, c[:66], flag))
        band = sum(1 for _, _, w, _, _ in rows if BAND[0] <= w <= BAND[1])
        print("\n%d paragraphs, %d in his 90-150 word band, %d over the limits"
              % (len(rows), band, len(over)))
        return 0
    over.sort(reverse=True)
    print("%-6s %-5s %-34s %s" % ("words", "sent", "section", "opens"))
    for w, ns, num, name, c in over:
        print("%-6d %-5d %-34s %s" % (w, ns, ("%d. %s" % (num, name))[:34], c[:64]))
    print("\n%d of %d paragraphs breach criterion 12 (>%d words) or 16a (>%d sentences)."
          % (len(over), len(rows), MAX_WORDS, MAX_SENTENCES))
    return 1 if over else 0


if __name__ == "__main__":
    sys.exit(main())
