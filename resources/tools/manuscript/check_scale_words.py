#!/usr/bin/env python3
r"""Scan for adjectives and adverbs of scale used in place of a number.

Written 2026-08-25 at MB's request: he would rather the review avoided "very"
and its family.  The distinction that matters in a physics review is not the
word but what it stands next to.

  QUANTIFIED   "roughly two orders of magnitude", "about $10^{-6}$", "a factor
               of three larger" -- the modifier attaches to a number and is
               doing honest work, hedging a stated value.

  BARE         "very good", "significantly better", "much smaller", "highly
               accurate" -- the modifier stands where a number could have, and
               the reader is told a size without being told the size.

Only BARE hits are a problem.  The script classifies each occurrence by looking
for a numeral, a power of ten, a percentage, or a spelled-out quantity within a
short window, and reports the two groups separately.

    python check_scale_words.py [--window 60]
"""
import argparse
import re
from pathlib import Path

TEX = Path(__file__).resolve().parents[2] / "paper" / "main.tex"

SCALE = r"""very quite rather somewhat fairly pretty extremely highly greatly
vastly enormously hugely immensely tremendously considerably substantially
significantly appreciably markedly exceedingly remarkably strikingly slightly
marginally mildly modestly moderately largely virtually practically relatively
comparatively""".split()

# a number, a power of ten, a percentage, an order-of-magnitude phrase
NUM = re.compile(r"\d|\\cdot|10\^|\bper ?cent|%|\border(s)? of magnitude\b|"
                 r"\b(one|two|three|four|five|six|seven|eight|nine|ten|dozen|"
                 r"half|twice|thrice)\b|\bfactor\b|\bdecade(s)?\b", re.I)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", type=int, default=60)
    a = ap.parse_args()
    raw = TEX.read_text(encoding="utf-8")
    raw = re.sub(r'\\begin\{(table|figure|equation|align|tabular)\*?\}.*?'
                 r'\\end\{(table|figure|equation|align|tabular)\*?\}', ' ', raw, flags=re.S)
    raw = re.sub(r'(?m)^%.*$', '', raw)
    # "rather than" often wraps across a line break, so test against the flattened
    # text rather than the raw line when deciding whether "rather" is a degree word.
    lines = raw.split("\n")
    sec = "(front)"
    quantified, bare = [], []
    for i, line in enumerate(lines, 1):
        m = re.match(r'\\section\{(.*?)\}', line)
        if m:
            sec = re.sub(r'[\\{}]', '', m.group(1))[:24]
        for w in SCALE:
            for mm in re.finditer(r"\b%s\b" % w, line, re.I):
                # "rather than" is a contrastive conjunction, not a degree adverb;
                # likewise "rather, ..." resuming a sentence.  Only degree uses count.
                if w == "rather":
                    nxt = " ".join((line[mm.end():] + " " +
                                    (lines[i] if i < len(lines) else "")).split())
                    if re.match(r"(than\b|,)", nxt):
                        continue
                lo = max(0, mm.start() - a.window)
                hi = min(len(line), mm.end() + a.window)
                ctx = " ".join(line[lo:hi].split())
                (quantified if NUM.search(line[lo:hi]) else bare).append((i, sec, w, ctx))
    words = len(re.findall(r"[A-Za-z][A-Za-z'\-]*", raw))
    for label, group in (("BARE -- a size given without the size", bare),
                         ("QUANTIFIED -- attached to a number, doing honest work", quantified)):
        print("\n=== %s: %d  (%.1f per 10k words) ===" % (label, len(group),
                                                          10000 * len(group) / words))
        for i, sc, w, ctx in group:
            print("  %-6d [%-24s] %-14s %s" % (i, sc, w, ctx))
    print("\n%d total over %d words" % (len(bare) + len(quantified), words))


if __name__ == "__main__":
    main()
