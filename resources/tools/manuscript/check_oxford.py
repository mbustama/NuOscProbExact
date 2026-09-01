#!/usr/bin/env python3
r"""Find three-or-more-item lists written without the serial (Oxford) comma.

Written 2026-08-25.  This is deliberately a *reporting* tool, not a rewriter: an
earlier regex pass over the same problem inserted a comma into "what does and
does not repair" and had to be reverted from a backup.  The failure mode is that
" X and Y" preceded by a comma is only sometimes a list -- it is just as often a
compound predicate, a two-item list inside a longer sentence, a correlative
("both A and B"), or a pair of verbs sharing a subject.  So the script prints
candidates with context and leaves the judgement to a reader.

Math, citations, refs and \texttt{} bodies are blanked first, since the commas
inside them are not prose commas.

    python check_oxford.py [--context 90]
"""
import argparse
import re
from pathlib import Path

TEX = Path(__file__).resolve().parents[2] / "paper" / "main.tex"

# The canonical short list: three parallel chunks of one to four words each,
# closing on punctuation.  Deliberately low-recall and high-precision -- a looser
# pattern returned 417 candidates, almost all of them "clause, X and Y" rather
# than lists, which is how the earlier pass went wrong.
WORD = r"[A-Za-z][A-Za-z0-9'\-]*"
CHUNK = r"(?:%s(?:\s+%s){0,3})" % (WORD, WORD)
CAND = re.compile(r"\b(%s),\s+(%s)\s+(and|or)\s+(%s)\s*(?=[.,;:]|$)"
                  % (CHUNK, CHUNK, CHUNK))


def blank(s):
    """Replace the innards of math and markup with placeholders of the same shape."""
    s = re.sub(r'\$[^$]*\$', ' MATH ', s)
    s = re.sub(r'\\(cite|ref|eqref|citep)\{[^}]*\}', ' REF ', s)
    s = re.sub(r'\\texttt\{[^}]*\}', ' CODE ', s)
    s = re.sub(r'\\emph\{([^}]*)\}', r'\1', s)
    s = re.sub(r'\\[a-zA-Z]+', ' ', s)
    return s


def sentences(block):
    return re.split(r'(?<=[.!?])\s+', block)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--context", type=int, default=96)
    a = ap.parse_args()
    raw = TEX.read_text(encoding="utf-8")
    raw = re.sub(r'\\begin\{(table|figure|equation|align|tabular)\*?\}.*?'
                 r'\\end\{(table|figure|equation|align|tabular)\*?\}', ' ', raw, flags=re.S)
    lines = raw.split("\n")
    sec, out = "(front)", []
    # keep a line index so a hit can be located in the source
    for i, line in enumerate(lines, 1):
        m = re.match(r'\\section\{(.*?)\}', line)
        if m:
            sec = re.sub(r'[\\{}]', '', m.group(1))[:28]
        out.append((i, sec, line))
    # join into a flowing text but remember where each character came from
    flat, origin = [], []
    for i, sec, line in out:
        flat.append(line + " ")
        origin.extend([(i, sec)] * (len(line) + 1))
    # blank() changes length; re-derive positions on the unblanked text instead
    text_raw = "".join(flat)
    hits = []
    for m in CAND.finditer(blank(text_raw)):
        # locate the matched tail in the raw text to get a line number
        tail = m.group(4).strip()[:24]
        j = text_raw.find(tail, max(0, m.start() - 400))
        ln, sc = origin[min(j, len(origin) - 1)] if j >= 0 else (0, "?")
        s = max(0, m.start() - a.context)
        hits.append((ln, sc, " ".join(blank(text_raw)[s:m.end() + 6].split())))
    print("candidate three-item lists without a serial comma: %d\n" % len(hits))
    for ln, sc, ctx in hits:
        print("  line %-6d [%-28s] ...%s" % (ln, sc, ctx))
    return len(hits)


if __name__ == "__main__":
    raise SystemExit(0 if main() == 0 else 0)
