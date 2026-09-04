#!/usr/bin/env python3
r"""Triage every quoted number in the review by how badly it could be wrong.

Why this exists.  On 2026-08-26 four numbers in Secs. 3.3 and 7.1 failed
measurement in one afternoon -- a percentage that was configuration-dependent, a
magnitude taken from a regime its own law does not cover, a "within a factor of
two" that held for 2 of 18 configurations, and an O(1) that was 1e-8.  All four
had been written first and checked afterwards, which is the wrong order.  Before
re-measuring 485 numerals it is worth knowing which ones can actually hurt.

The classification, in decreasing order of risk:

  MEASURED   a mantissa of two or more significant figures in running prose.
             These are the run outputs: nothing in the physics fixes a second
             digit, so the number came from an execution and may not survive a
             different seed, spectrum, or library.  Every one of the four
             failures above was of this kind.
  SCALE      a bare power of ten, or a single leading digit.  An order-of-
             magnitude claim.  Wrong only if the exponent is wrong.
  ANCHORED   sits next to \cite, \ref or \eqref, or inside a table or display.
             Traceable to a source or to an equation, so a reader can check it.
  BENIGN     years, section and equation numbers, flavor counts, small integers.

Running prose only: displays, tabulars, captions of pure data tables and the
bibliography are excluded, since a number in a display is part of a derivation
and a number in a parameter table has a citation beside it.

    python3 check_numerals.py            # whole paper, ranked by section
    python3 check_numerals.py --list     # every MEASURED numeral, with context
"""
import argparse
import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TEX = os.path.join(HERE, os.pardir, os.pardir, "paper", "main.tex")

# a number with two or more significant figures, in any of the forms the review
# uses: 1.65, 2.3\cdot10^{-5}, 7.8e-5, 23\%
MANTISSA = re.compile(r"(?<![\d.])(\d+\.\d+|\d{2,})(?![\d.])")
POWER10 = re.compile(r"10\^\{?-?\d+\}?")
BENIGN_CTX = re.compile(r"\bCITE\b|\bREF\b")


def strip_noncprose(t):
    """Remove everything that is not running prose, keeping line structure."""
    t = re.sub(r"(?m)^\s*%.*$", "", t)
    for env in ("equation", "align", "aligned", "array", "gather", "tabular",
                "table", "figure", "thebibliography"):
        t = re.sub(r"\\begin\{%s\*?\}.*?\\end\{%s\*?\}" % (env, env), " ", t,
                   flags=re.S)
    # Inline math is NOT stripped: almost every quoted number in this review
    # sits inside $...$ ("$1.3$~eV$^2$"), so removing it hides exactly the
    # numerals this is for.  The first version of this script did strip it and
    # reported Sec. 7 as carrying none, which is how the bug was found.  Only
    # the symbols are dropped.  Citations and cross-references are turned into
    # markers FIRST -- stripping them as commands is what broke the ANCHORED
    # test on the second attempt, reporting zero anchored numerals paper-wide.
    t = re.sub(r"\\cite[a-z]*(\[[^\]]*\])?\{[^}]*\}", " CITE ", t)
    t = re.sub(r"\\(eq)?ref\{[^}]*\}", " REF ", t)
    t = re.sub(r"\\[a-zA-Z]+", " ", t)
    # Subscripts and superscripts are labels and exponents, not quoted values:
    # theta_{34}, 10^{-16}, Delta m^2_{31}.  A hand spot-check of the first
    # output found these were most of what it flagged.
    t = re.sub(r"[_^]\{[^{}]*\}", " ", t)
    t = re.sub(r"[_^]-?\w", " ", t)
    # LaTeX lengths and package version numbers are not quoted physics.
    t = re.sub(r"\\includegraphics(\[[^\]]*\])?\{[^}]*\}", " ", t)
    t = re.sub(r"(width|height|scale)\s*=\s*[\d.]+\s*\\?\w*", " ", t)
    t = re.sub(r"NuFIT~?[\d.]+", " NuFIT ", t)
    # A stated working precision is a choice, not a result.
    t = re.sub(r"\$?\d+\$?([~ ]|\\,)*-?digit", " ", t)
    return t


def sections(src):
    ms = list(re.finditer(r"(?m)^\\section\*?\{(.*?)\}", src))
    out = []
    for i, m in enumerate(ms):
        b = ms[i + 1].start() if i + 1 < len(ms) else len(src)
        out.append((m.group(1), m.start(), src[m.start():b]))
    return out


def classify(line, num):
    """Risk class for one numeral, given the prose line it sits in."""
    if BENIGN_CTX.search(line):
        return "ANCHORED"
    if re.match(r"^(19|20)\d\d$", num):
        return "BENIGN"
    if re.match(r"^\d+$", num) and int(num) <= 12:
        return "BENIGN"
    if "." in num or len(num.lstrip("0")) >= 2:
        return "MEASURED"
    return "SCALE"


def main(show_list):
    src = io.open(TEX, encoding="utf-8").read()
    rows, per = [], {}
    for title, off, body in sections(src):
        clean = strip_noncprose(body)
        n = {"MEASURED": 0, "SCALE": 0, "ANCHORED": 0, "BENIGN": 0}
        for line in clean.split("\n"):
            if not line.strip():
                continue
            for m in MANTISSA.finditer(line):
                k = classify(line, m.group(1))
                n[k] += 1
                if k == "MEASURED":
                    a = max(0, m.start() - 60)
                    rows.append((title, m.group(1),
                                 re.sub(r"\s+", " ", line[a:m.end() + 60])))
            n["SCALE"] += len(POWER10.findall(line))
        per[title] = n

    if show_list:
        for t, v, ctx in rows:
            print("  %-34s %-10s ...%s" % (t[:34], v, ctx))
        print()

    print("%-42s %8s %6s %8s %6s" % ("section", "MEASURED", "SCALE",
                                     "ANCHORED", "BENIGN"))
    tot = {"MEASURED": 0, "SCALE": 0, "ANCHORED": 0, "BENIGN": 0}
    for t, n in sorted(per.items(), key=lambda kv: -kv[1]["MEASURED"]):
        if sum(n.values()) == 0:
            continue
        print("%-42s %8d %6d %8d %6d"
              % (t[:42], n["MEASURED"], n["SCALE"], n["ANCHORED"], n["BENIGN"]))
        for k in tot:
            tot[k] += n[k]
    print("%-42s %8d %6d %8d %6d"
          % ("TOTAL", tot["MEASURED"], tot["SCALE"], tot["ANCHORED"],
             tot["BENIGN"]))
    print("\n%d numerals need re-measuring; the other %d do not."
          % (tot["MEASURED"], tot["SCALE"] + tot["ANCHORED"] + tot["BENIGN"]))
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true",
                    help="print every MEASURED numeral with its context")
    a = ap.parse_args()
    sys.exit(main(a.list))
