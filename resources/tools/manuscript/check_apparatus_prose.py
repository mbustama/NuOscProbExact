#!/usr/bin/env python3
r"""Find prose whose subject is the paper's own apparatus rather than its physics.

Written 2026-08-25, after Secs. 11 and 12 turned out to carry several paragraphs
about tables, columns and typography instead of about neutrinos.  The failure has
five forms -- inventory ("this section is four tables"), reading instructions
("the last column is the one to read"), bookkeeping (membership rules, citation
conventions, one-row exceptions), pre-emptive apology ("most readers will not
face any of this"), and counts that do not survive checking, which is the
diagnostic tell, since apparatus-prose is never audited the way physics-prose is.

Keyword lists catch only the phrasings one already knows, so this measures the
property instead: per paragraph, the ratio of apparatus words to physics words.
A high ratio does not convict a paragraph -- a sentence introducing a float is
supposed to name it -- it says where to read.

    python check_apparatus_prose.py [--top 25] [--min-words 40]
    python check_apparatus_prose.py --routing      # the reader-routing sweep

The --routing mode is the one to run before every commit.  Reader-routing --
"a reader who already knows the derivation can go straight to Sec. 3", "and can
skip the rest" -- is the form this review kept regenerating: six instances in
2026-08, two of them the same sentence sixty-six lines apart, and five of them
inviting the reader out of the section they open.  It is cheap to detect and
almost invisible on a read-through, so it gets its own pass.
"""
import argparse
import re
from pathlib import Path

TEX = Path(__file__).resolve().parents[2] / "paper" / "main.tex"

APPARATUS = r"""table tables column columns row rows entry entries cell cells list
lists caption captions panel panels figure figures section sections subsection
appendix index tabulated listed above below overleaf footnote item items line
lines page pages column's"""
PHYSICS = r"""neutrino neutrinos oscillation oscillations probability probabilities
eigenvalue eigenvalues eigenvector eigenvectors matter density flavor flavors
hamiltonian projector projectors mixing phase phases amplitude amplitudes energy
baseline matrix spectrum resonance vacuum mass masses potential coherence damping
error accuracy gap degeneracy exponential polynomial roots invariant unitary
propagation profile layer detector experiment"""


# Reader-routing: sentences that hand the reader an exit rather than content.
ROUTING = [
    (r"\b(a|the|any)\s+reader\s+who\b", "'a reader who ...'"),
    (r"\breaders?\s+(who|whose)\b", "'readers who/whose ...'"),
    (r"\b(can|may|could|might)\s+skip\b", "'can skip'"),
    (r"\bskip(s|ped|ping)?\s+(the|this|everything|ahead|past)\b", "'skip the ...'"),
    (r"\bgo\s+straight\s+to\b|\b(turn|proceed|jump)\s+(directly\s+)?to\s+Sec", "'go straight to'"),
    (r"\bis\s+not\s+obliged\b|\bneed\s+not\s+read\b|\bwithout\s+reading\b", "'is not obliged'"),
    (r"\btake\s+(them|it)\s+from\s+there\b", "'take it from there'"),
    (r"\bshould\s+(read|start\s+at|begin\s+with)\b", "'should read'"),
    (r"\bis\s+for\s+readers\b|\bfor\s+those\s+(who|looking)\b", "'is for readers ...'"),
    (r"\bwho\s+wants?\s+(the|only|just)\b", "'who wants the ...'"),
]


def routing_sweep():
    """Report every reader-routing construction, with its section."""
    s = TEX.read_text(encoding="utf-8")
    lines = s.split("\n")
    sec, hits = "(front matter)", []
    for i, line in enumerate(lines, 1):
        m = re.match(r'\\section\{(.*?)\}', line)
        if m:
            sec = re.sub(r'[\\{}]', '', m.group(1))
        for pat, name in ROUTING:
            if re.search(pat, line, re.I):
                hits.append((i, sec[:34], name, " ".join(line.split())[:88]))
                break
    print("reader-routing constructions: %d\n" % len(hits))
    for i, sc, name, t in hits:
        print("  %-6d %-34s %-22s %s" % (i, sc, name, t))
    return len(hits)


def tokens(block):
    return re.findall(r"[A-Za-z']+", block.lower())


def load_paragraphs():
    s = TEX.read_text(encoding="utf-8")
    s = re.sub(r'\\begin\{(table|figure)\*?\}.*?\\end\{(table|figure)\*?\}', ' ', s, flags=re.S)
    s = re.sub(r'\\begin\{description\}.*?\\end\{description\}', ' ', s, flags=re.S)
    out, sec, buf = [], "(front matter)", []
    for line in s.split("\n"):
        m = re.match(r'\\section\{(.*?)\}', line)
        if m:
            sec = re.sub(r'[\\{}]', '', m.group(1))
        if line.strip():
            buf.append(line)
        elif buf:
            out.append((sec, "\n".join(buf)))
            buf = []
    if buf:
        out.append((sec, "\n".join(buf)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--min-words", type=int, default=40)
    ap.add_argument("--routing", action="store_true",
                    help="sweep for reader-routing constructions and exit")
    a = ap.parse_args()
    if a.routing:
        raise SystemExit(1 if routing_sweep() else 0)
    app, phy = set(APPARATUS.split()), set(PHYSICS.split())
    rows = []
    for sec, block in load_paragraphs():
        if block.lstrip().startswith(("\\sub", "\\label", "%", "\\section",
                                      "\\paragraph", "\\begin", "\\end", "\\item")):
            continue
        # count structural references as apparatus too
        refs = len(re.findall(r'\\(?:ref|eqref)\{(?:tab|fig|sec|app):', block))
        w = tokens(block)
        if len(w) < a.min_words:
            continue
        na = sum(1 for x in w if x in app) + refs
        np_ = sum(1 for x in w if x in phy)
        rows.append((na / (np_ + 1.0), na, np_, len(w), sec,
                     " ".join(block.split())[:96]))
    rows.sort(reverse=True)
    print("%-6s %-4s %-4s %-5s %-30s %s" % ("ratio", "app", "phys", "words", "section", "opens"))
    for r, na, np_, nw, sec, head in rows[:a.top]:
        print("%-6.2f %-4d %-4d %-5d %-30s %s" % (r, na, np_, nw, sec[:30], head))


if __name__ == "__main__":
    main()
