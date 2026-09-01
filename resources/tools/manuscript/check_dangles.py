#!/usr/bin/env python3
r"""Criterion (13) over the whole document: dangling words on a line of their own.

The section-by-section detector lives in audit_section.py; this runs the same two
mechanical cases across all of main.pdf at once, which is how criterion (13) is
meant to be run -- last, and repeatedly, because every reflow moves the line
breaks and fixing one dangle routinely creates another.

  (i)  the final sentence of a paragraph left as a few words on their own,
       captions included;
  (ii) one or two words of a lead-in strandeded on the line before a display.

Case (iii) of the criterion -- a sentence resuming after a display whose trailing
inline formula strands a symbol -- is not detectable here, because the stranded
piece is math rather than prose.  It has to be read off the typeset page after
any equation a sentence runs on from.

    pdftotext -layout paper/main.pdf /tmp/main.txt
    python check_dangles.py /tmp/main.txt
"""
import io
import re
import sys

EQN = re.compile(r"\(\d+\.\d+\)\s*$")
MATH = re.compile(r"[=+∑∏∫≃≈≡⟨⟩±→"
                  r"∆√≤≥α-ωΔΨΣΠ∗|]")
HEAD = re.compile(r"^\s*(\d+(\.\d+)*\.?\s+[A-Z]|Table \d|Figure \d|Appendix)")
GAP = re.compile(r"\S {4,}\S")          # table column gap or equation layout


def load(path):
    """Return [(page_number, line), ...] with the original page numbers kept.

    Pages are skipped, never deleted, so a reported page number is the page of
    the PDF and can be looked at directly.
    """
    raw = io.open(path, encoding="utf-8", errors="replace").read().split("\f")
    out, in_refs = [], False
    for n, q in enumerate(raw, 1):
        if re.search(r"^\s*References\s*$", q, re.M):
            in_refs = True          # the reference list runs to the end
        if in_refs:
            continue
        if re.search(r"^\s*Key results\s*$", q, re.M):
            continue                # collected key-result list, not prose
        q = re.sub(r"^\s*Key result~?\s*\d+\.\d+\s*$.*?(?=\n\s*\n)", "", q,
                   flags=re.M | re.S)
        q = re.sub(r"^\s*\d+\.\d+(\.\d+)?\.?\s+[A-Z][^\n]{0,60}$", "", q, flags=re.M)
        for ln in q.split("\n"):
            out.append((n, ln.rstrip()))
    return out


def prose(x):
    """A justified body line of running prose, not a table row or an equation."""
    return bool(x.strip() and not EQN.search(x) and not HEAD.match(x)
                and not GAP.search(x)
                and len(x) - len(x.lstrip()) <= 20
                and len(re.findall(r"[A-Za-z]{3,}", x)) >= 1)


def page_of(pg, idx):
    return pg[:idx].count("\f") + 1


def main(path):
    rows = load(path)
    pages = [p for p, _ in rows]
    gl = [x for _, x in rows]

    ends = []
    for k, x in enumerate(gl):
        if not prose(x):
            continue
        txt = x.strip()
        if len(txt) >= 48:
            continue
        if k == 0 or not gl[k - 1].strip() or HEAD.match(gl[k - 1] or ""):
            continue
        prv = gl[k - 1]
        if not (prose(prv) and len(prv.rstrip()) >= 60):
            continue
        if not txt.endswith(('.', ';', ':', '!', '?')):
            continue
        ends.append((pages[k], txt, prv.strip()[-52:]))

    eqs = []
    for q, ln in enumerate(gl):
        if not EQN.search(ln):
            continue
        r = q - 1
        while r > 0:
            c = gl[r]
            if c.strip() and not EQN.search(c):
                ind = len(c) - len(c.lstrip())
                if (ind <= 10 and len(c.strip()) < 40
                        and re.search(r"[A-Za-z]{3}", c) and not MATH.search(c)):
                    eqs.append((pages[r], c.strip(), ln.strip()[:44]))
                if ind <= 10 and re.search(r"[A-Za-z]{4}", c):
                    break
            r -= 1

    print("(13) DANGLING WORDS -- whole document\n")
    print("  (i) paragraph-end dangles: %d" % len(ends))
    for p, t, prv in ends:
        print("      p%-4d %-46r after: ...%s" % (p, t, prv))
    print("\n  (ii) stranded lead-in before a display: %d" % len(eqs))
    for p, t, e in eqs:
        print("      p%-4d %-46r eq: %s" % (p, t, e))
    total = len(ends) + len(eqs)
    print("\n  %d total -> %s" % (total, "CLEAN" if total == 0 else "FIX AND RERUN"))
    return total


if __name__ == "__main__":
    sys.exit(1 if main(sys.argv[1]) else 0)
