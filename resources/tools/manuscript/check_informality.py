#!/usr/bin/env python3
r"""Audit for colloquial register in a technical review.

Written 2026-08-26, after MB flagged "That silence has a cost", "each at a
price", "We can make this short", "has cost someone a result", and "Five
assumptions stand behind that calculation, made rather than argued for".

The fault is register, not sentence length: short plain declaratives are the
house style and are fine.  What does not belong is figurative, idiomatic, or
conversational writing.  Five families:

  A FIGURATIVE   personification and metaphor for abstractions -- silence that
                 has a cost, traps that punish, accuracy that drains.
  B COMMERCE     buying, paying, prices and bargains used for accuracy or cost.
  C INDEFINITE   colloquial indefinites and asides -- someone, nobody, anyone,
                 of course, after all, as it happens.
  D EVALUATIVE   casual value words -- cheap, wins, beats, worth having.
  E TAG          trailing aphoristic modifiers appended for effect -- "made
                 rather than argued for", "and no further", "each at a price".

Every hit is printed with context.  Some are legitimate: "lives in the adjoint
representation" is standard mathematics, and "costs" is literal when counting
operations.  The script reports; a reader judges.

    python check_informality.py [--section SUBSTRING]
"""
import argparse
import re
from pathlib import Path

TEX = Path(__file__).resolve().parents[2] / "paper" / "main.tex"

FAMILIES = [
 ("A FIGURATIVE", [
   r"\b(silence|accuracy|the field|the literature|memory|history|the formula|the method)\s+(has|have|pays?|paid|punish\w*|forgive\w*|remember\w*|forget\w*|want\w*|refus\w*|know\w*)\b",
   r"\bdrains?\b|\bdrained\b|\bdraining\b|\bbleeds?\b|\bstarv\w+\b",
   r"\b(traps?|the price|the cost)\s+(punish\w*|catch\w*|bite\w*|claim\w*)",
   r"\bsurvives? scrutiny\b|\bcomes? alive\b|\bstands? guard\b",
 ]),
 ("B COMMERCE", [
   r"\bat a price\b|\ba bad trade\b|\ba good trade\b|\bfor free\b|\bbuys? nothing\b",
   r"\b(accuracy|speed|headroom|conditioning|digits?)\s+(is|are|was|were)?\s*(bought|paid for|purchased)\b",
   r"\bbuys?\b(?!\s+(?:back|headroom))",
 ]),
 ("C INDEFINITE", [
   r"\b(someone|somebody|nobody|no one|anyone|everyone)\b",
   r"\bof course\b|\bafter all\b|\bas it happens\b|\bto be fair\b|\bneedless to say\b",
   r"\bone would (least )?expect\b|\ba reader would guess\b",
 ]),
 ("D EVALUATIVE", [
   r"\b(is|are|was|were|it is|seems?)\s+(cheap|costly|dear|pricey)\b",
   r"\bwins?\b|\bbeats?\b(?!\s+frequenc)|\bloses? out\b|\bcomes? out ahead\b",  # "beat frequencies" is a physics term
   r"\bworth (having|the effort|it)\b",
 ]),
 ("E TAG", [
   r",\s+made rather than argued for\b|,\s+each at a price\b|\band no further\b",
   r",\s+not \w+ (?:for|at|by) \w+\.",
   r"---\s*and nothing (?:more|else)\b|;\s*nothing more\.",
 ]),
]


def sections():
    raw = TEX.read_text(encoding="utf-8")
    for e in ("table", "figure", "tikzpicture", "tabular", "equation", "align", "description"):
        raw = re.sub(r'\\begin\{%s\*?\}.*?\\end\{%s\*?\}' % (e, e), ' ', raw, flags=re.S)
    raw = re.sub(r'(?m)^%.*$', '', raw)
    ch = re.split(r'\\section\{(.*?)\}', raw)
    out = []
    for i in range(1, len(ch) - 1, 2):
        body = " ".join(ch[i + 1].split())
        out.append((re.sub(r'[\\{}]', '', ch[i])[:36], body))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--section")
    a = ap.parse_args()
    grand = 0
    for title, body in sections():
        if a.section and a.section.lower() not in title.lower():
            continue
        rows = []
        for fam, pats in FAMILIES:
            for pat in pats:
                for m in re.finditer(pat, body, re.I):
                    ctx = " ".join(body[max(0, m.start() - 92):m.end() + 60].split())
                    rows.append((fam, m.group(0)[:26], ctx))
        if rows:
            print("\n=== %s : %d ===" % (title, len(rows)))
            seen = set()
            for fam, hit, ctx in rows:
                if ctx in seen:
                    continue
                seen.add(ctx)
                print("  [%s] %-26s ...%s" % (fam[0], hit, ctx))
            grand += len(rows)
    print("\ntotal hits: %d" % grand)


if __name__ == "__main__":
    main()
