#!/usr/bin/env python3
r"""Scan for the marks of amateur or merely serviceable scientific prose.

Written 2026-08-25, after the apparatus-prose and machine-writing sweeps.  This
one looks for weakness rather than for a machine: hedges that quantify nothing,
nominalizations standing in for verbs, filler openings, redundant phrases,
unearned intensifiers, and vague pronoun reference.

Every category prints its hits with context.  Almost none of them is wrong on
its own -- "very" has legitimate uses, "note that" occasionally earns its place,
passive voice is right when the object is the subject of interest.  The signal
is density and repetition, so counts are reported per 10,000 words alongside the
raw hits.

    python check_weak_prose.py [--context 104]
"""
import argparse
import re
from pathlib import Path

TEX = Path(__file__).resolve().parents[2] / "paper" / "main.tex"

CATEGORIES = [
    ("HEDGE / FILLER WORD", [
        (r"\bvery\b", "very"), (r"\bquite\b", "quite"), (r"\bsomewhat\b", "somewhat"),
        (r"\bfairly\b", "fairly"), (r"\bbasically\b", "basically"),
        (r"\bessentially\b", "essentially"), (r"\bactually\b", "actually"),
        (r"\breally\b", "really"), (r"\bof course\b", "of course"),
        (r"\bclearly\b", "clearly"), (r"\bobviously\b", "obviously"),
        (r"\bit seems\b|\bit appears that\b", "it seems"),
        (r"\barguably\b", "arguably"),
    ]),
    ("UNQUANTIFIED INTENSIFIER", [
        (r"\bextremely\b", "extremely"), (r"\bhighly\b", "highly"),
        (r"\bgreatly\b", "greatly"), (r"\bvastly\b", "vastly"),
        (r"\bdramatically\b", "dramatically"), (r"\bsubstantially\b", "substantially"),
        (r"\bsignificantly\b", "significantly"), (r"\bconsiderably\b", "considerably"),
    ]),
    ("FILLER OPENING", [
        (r"\bit should be noted\b|\bit is worth (noting|mentioning)\b", "it should be noted"),
        (r"\bit is (clear|obvious|evident) that\b", "it is clear that"),
        (r"\bas (mentioned|noted|discussed) (above|earlier|previously)\b", "as mentioned above"),
        (r"\bas we (have )?(seen|saw|discussed)\b", "as we have seen"),
        (r"\bit is well known\b", "it is well known"),
        (r"(?m)^\s*(Also|Then|So|But|And)\s*,", "weak transition opener"),
        (r"\bnote that\b", "note that"),
    ]),
    ("REDUNDANCY", [
        (r"\bin order to\b", "in order to"),
        (r"\bdue to the fact that\b|\bowing to the fact\b", "due to the fact that"),
        (r"\ba total of\b", "a total of"),
        (r"\bcompletely (eliminat|remov)", "completely eliminate"),
        (r"\babsolutely (essential|necessary)\b", "absolutely essential"),
        (r"\bin spite of the fact\b", "in spite of the fact"),
        (r"\bat this point in time\b", "at this point in time"),
        (r"\bthe reason (is|being) because\b", "reason is because"),
    ]),
    ("NOMINALIZATION / DEAD VERB", [
        (r"\bmake use of\b", "make use of"),
        (r"\bperform(s|ed|ing)? (a|an|the) (calculation|comparison|analysis|evaluation)\b",
         "perform a calculation"),
        (r"\bcarr(y|ies|ied) out (a|an|the) (calculation|comparison|analysis)\b",
         "carry out a calculation"),
        (r"\bgives? rise to\b", "gives rise to"),
        (r"\bis in agreement with\b", "is in agreement with"),
        (r"\bhas the ability to\b", "has the ability to"),
        (r"\bserve(s)? (the purpose of|to provide)\b", "serves the purpose of"),
    ]),
    ("VAGUE WORD", [
        (r"\bthing(s)?\b", "thing(s)"), (r"\bstuff\b", "stuff"),
        (r"\bnice\b", "nice"), (r"\bgood\b", "good"), (r"\bbad\b", "bad"),
        (r"\ba number of\b", "a number of"), (r"\bvarious\b", "various"),
    ]),
    ("CLICHE", [
        (r"\btip of the iceberg\b", "tip of the iceberg"),
        (r"\bat the end of the day\b", "at the end of the day"),
        (r"\blast but not least\b", "last but not least"),
        (r"\bthe devil is in the detail", "devil in the details"),
        (r"\bhand in hand\b", "hand in hand"),
        (r"\ba double[- ]edged sword\b", "double-edged sword"),
    ]),
]


def clean(raw):
    raw = re.sub(r'\\begin\{(table|figure|equation|align|tabular)\*?\}.*?'
                 r'\\end\{(table|figure|equation|align|tabular)\*?\}', ' ', raw, flags=re.S)
    return raw


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--context", type=int, default=104)
    a = ap.parse_args()
    raw = clean(TEX.read_text(encoding="utf-8"))
    lines = raw.split("\n")
    words = len(re.findall(r"[A-Za-z][A-Za-z'\-]*", raw))
    total = 0
    for label, pats in CATEGORIES:
        hits, sec = [], "(front)"
        for i, line in enumerate(lines, 1):
            m = re.match(r'\\section\{(.*?)\}', line)
            if m:
                sec = re.sub(r'[\\{}]', '', m.group(1))[:24]
            if line.lstrip().startswith("%"):
                continue
            for pat, name in pats:
                if re.search(pat, line, re.I):
                    hits.append((i, sec, name, " ".join(line.split())[:a.context]))
        total += len(hits)
        print("\n=== %s: %d  (%.1f per 10k words) ===" % (label, len(hits),
                                                          10000 * len(hits) / words))
        for i, sc, name, t in hits:
            print("  %-6d [%-24s] %-24s %s" % (i, sc, name, t))
    print("\ntotal flags: %d over %d words" % (total, words))
    return total


if __name__ == "__main__":
    main()
    raise SystemExit(0)
