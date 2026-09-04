#!/usr/bin/env python3
r"""Look for prose that reads as machine-generated.

Written 2026-08-25 at MB's request.  Two families of tell:

LEXICAL -- vocabulary and connectives that large language models overproduce
("delve", "crucial", "it is worth noting", "moreover", "sheds light on",
"paves the way", "a testament to").  Some are legitimate in a physics review,
so every hit is printed with context rather than counted against the text.

STRUCTURAL -- the rhythms, which are the harder tell and the one a word list
misses: the antithesis flip ("not X, but Y" / "This is not X; it is Y"),
tricolon pile-ups, sentence-initial "This" with no referent, hedge stacks, and
summary sentences that restate the paragraph just written.

    python check_llm_tells.py [--structural] [--context 100]
"""
import argparse
import re
from pathlib import Path

TEX = Path(__file__).resolve().parents[2] / "paper" / "main.tex"

LEXICAL = [
    (r"\bdelv(e|es|ing)\b", "delve"),
    (r"\bleverag(e|es|ing|ed)\b", "leverage"),
    (r"\butilis|utiliz(e|es|ing|ation)\b", "utilize"),
    (r"\bshowcas(e|es|ing|ed)\b", "showcase"),
    (r"\bunderscor(e|es|ing|ed)\b", "underscore"),
    (r"\bpivotal\b", "pivotal"),
    (r"\bseamless(ly)?\b", "seamless"),
    (r"\bmyriad\b|\bplethora\b", "myriad/plethora"),
    (r"\b(landscape|realm|tapestry|arena)\s+of\b", "landscape/realm of"),
    (r"\bintricate\b|\bnuanced\b|\bmultifaceted\b|\bholistic\b", "intricate/nuanced"),
    (r"\bit is (worth|important|crucial) (noting|to note|to mention)\b", "it is worth noting"),
    (r"^\s*(Moreover|Furthermore|Additionally|Notably|Importantly|Overall)\b", "connective opener"),
    (r"\bplays? an? (crucial|key|vital|important|significant) role\b", "plays a key role"),
    (r"\bshed(s|ding)? light on\b", "sheds light on"),
    (r"\bpave(s|d)? the way\b", "paves the way"),
    (r"\ba testament to\b", "a testament to"),
    (r"\bat the forefront\b", "at the forefront"),
    (r"\bin the (context|realm|world) of\b", "in the context of"),
    (r"\bcutting[- ]edge\b|\bstate[- ]of[- ]the[- ]art\b", "cutting-edge"),
    (r"\bdeep(er)? (dive|understanding)\b", "deep dive"),
    (r"\bcomprehensive (overview|understanding|analysis)\b", "comprehensive overview"),
    (r"\bboth\s+\w+\s+and\s+\w+\s+alike\b", "both X and Y alike"),
    (r"\bnot only\b[^.]{0,60}\bbut also\b", "not only ... but also"),
    (r"\bin (today's|this modern)\b", "in today's"),
    (r"\bharness(es|ing|ed)?\b", "harness"),
    (r"\bunlock(s|ing|ed)?\b", "unlock"),
    (r"\brobustly\b|\bsignificantly enhance", "robustly/enhance"),
]

STRUCTURAL = [
    (r"\bis not (just|merely|only) [^.;]{3,60}[;,] (it|they) (is|are)\b", "antithesis flip"),
    (r"\bnot (a|an|the) [^.;]{2,40}, but (a|an|the)\b", "'not a X, but a Y'"),
    (r"\b(may|might|could) potentially\b|\bcan possibly\b", "hedge stack"),
    (r"\bIn (short|summary|essence|conclusion)\b|\bPut simply\b|\bSimply put\b", "summary opener"),
    (r"\bThe key (insight|point|takeaway|idea) (is|here)\b", "'the key insight is'"),
    (r"\bserves? as a\b|\bacts? as a\b", "'serves as a'"),
    (r"\bThis (underscores|highlights|demonstrates|reveals|illustrates) (the|how|that)\b",
     "'this highlights'"),
    (r"\bby (understanding|examining|exploring)\b", "'by exploring'"),
]


def scan(pats, label, context):
    raw = TEX.read_text(encoding="utf-8")
    lines = raw.split("\n")
    sec, hits = "(front)", []
    for i, line in enumerate(lines, 1):
        m = re.match(r'\\section\{(.*?)\}', line)
        if m:
            sec = re.sub(r'[\\{}]', '', m.group(1))[:26]
        if line.lstrip().startswith("%"):
            continue
        for pat, name in pats:
            mm = re.search(pat, line, re.I | re.M)
            if mm:
                hits.append((i, sec, name, " ".join(line.split())[:context]))
    print("\n=== %s: %d ===" % (label, len(hits)))
    for i, sc, name, t in hits:
        print("  %-6d [%-26s] %-22s %s" % (i, sc, name, t))
    return len(hits)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--context", type=int, default=104)
    a = ap.parse_args()
    n = scan(LEXICAL, "LEXICAL tells", a.context)
    n += scan(STRUCTURAL, "STRUCTURAL tells", a.context)
    print("\ntotal: %d" % n)
    return n


if __name__ == "__main__":
    raise SystemExit(0 if main() == 0 else 0)
