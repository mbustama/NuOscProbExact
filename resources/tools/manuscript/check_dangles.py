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

    python check_dangles.py paper/main.pdf        # preferred
    python check_dangles.py /tmp/main.txt         # single-column only

GIVE IT THE PDF.  A two-column page extracted with ``pdftotext -layout`` puts
both columns on one line, so a short line in the left column is padded out by
whatever sits beside it and the tail becomes invisible.  Handed a .txt of this
manuscript the detector reported five hits, four of them bibliography
artifacts, while twenty-one prose tails stood on the page, one of them a single
character.  Handed the PDF it splits each page down the middle and reads the
columns separately, which is the only way the measure means anything.

The threshold scales to the measure of the column it is reading, as part D of
CHECKLIST.md requires, rather than to a fixed character count: appendix and
caption columns are set narrower than the body and a fixed cut misses tails
that are plainly short on the page.
"""
import io
import re
import subprocess
import sys

EQN = re.compile(r"\(\d+\.\d+\)\s*$")
MATH = re.compile(r"[=+∑∏∫≃≈≡⟨⟩±→"
                  r"∆√≤≥α-ωΔΨΣΠ∗|]")
HEAD = re.compile(r"^\s*(\d+(\.\d+)*\.?\s+[A-Z]|Table \d|Figure \d|Appendix)")
GAP = re.compile(r"\S {4,}\S")          # table column gap or equation layout


def columns_of(pdf):
    """Yield ``(page, column_letter, [lines])`` for a two-column PDF.

    The split is found rather than assumed.  Cutting the page at a geometric
    midpoint clips characters off the inner edge of a column -- it turned
    "no approximation covers." into a tail reading "n covers." -- so instead
    each page is extracted whole and the gutter located as the run of
    character columns that is blank on nearly every full line.  A page with no
    such run is treated as single-column and yielded as one block.
    """
    out = subprocess.check_output(["pdfinfo", pdf]).decode("utf-8", "replace")
    pages = int(re.search(r"Pages:\s+(\d+)", out).group(1))
    for n in range(1, pages + 1):
        txt = subprocess.check_output(
            ["pdftotext", "-layout", "-f", str(n), "-l", str(n), pdf, "-"],
            stderr=subprocess.DEVNULL).decode("utf-8", "replace")
        lines = [ln.rstrip("\n") for ln in txt.split("\n")]
        long_lines = [x for x in lines if len(x.strip()) >= 40]
        if len(long_lines) < 6:
            yield n, "-", lines
            continue
        width = max(len(x) for x in long_lines)
        blank = [i for i in range(width)
                 if sum(1 for x in long_lines if i >= len(x) or x[i] == " ")
                 >= 0.92 * len(long_lines)]
        runs, cur = [], []
        for i in blank:
            if cur and i == cur[-1] + 1:
                cur.append(i)
            else:
                if len(cur) >= 3:
                    runs.append(cur)
                cur = [i]
        if len(cur) >= 3:
            runs.append(cur)
        mid = [r for r in runs if 0.30 * width < r[0] < 0.70 * width]
        if not mid:
            yield n, "-", lines
            continue
        g0, g1 = mid[0][0], mid[0][-1] + 1
        yield n, "L", [x[:g0].rstrip() for x in lines]
        yield n, "R", [x[g1:].rstrip() for x in lines]


def measure_of(lines):
    """The column's own measure: its typical full line, not its longest."""
    lens = sorted(len(x.rstrip()) for x in lines if len(x.strip()) > 20)
    if len(lens) < 8:
        return None
    return lens[int(0.9 * (len(lens) - 1))]


BIB = re.compile(r"^\s*\[\d+\]|\bdoi:|\barXiv:|\bPhys\.\s|\bRev\.\s")


def is_bibliography(lines):
    """A column of reference entries, not of prose.

    The reference list and the appendices interleave: on the pages where
    the bibliography ends, it holds the left column while Appendix A has
    already started in the right.  A flag that latches at the ``References``
    heading and runs to the end of the document would swallow the
    appendices with it, so each column is judged on its own contents.
    """
    body = [x for x in lines if len(x.strip()) > 12]
    if len(body) < 6:
        return False
    return sum(1 for x in body if BIB.search(x)) >= 0.30 * len(body)


APPX = re.compile(r"^\s*Appendix\b")


def strip_leading_bibliography(lines):
    """Blank the reference entries a column may open with.

    On the pages where the appendices begin, the bibliography ends part
    way down a column and prose follows it, so the whole-column ratio
    test above does not fire.  Walk from the top instead, and stop at the
    first appendix heading.
    """
    body = [x for x in lines if x.strip()]
    if not body or not any(BIB.search(x) for x in body[:5]):
        return lines
    out, in_bib = [], True
    for ln in lines:
        if in_bib and APPX.match(ln):
            in_bib = False
        out.append("" if in_bib else ln)
    return out


def load_pdf(path):
    """``[(page, line), ...]`` read column by column, references dropped."""
    rows = []
    for n, _side, lines in columns_of(path):
        if is_bibliography(lines):
            continue
        lines = strip_leading_bibliography(lines)
        rows.append((n, None))          # column boundary, so tails do not span
        for ln in lines:
            # front-matter fields are not prose and cannot dangle
            if re.match(r"^\s*(Keywords?|PACS|MSC|Program summary)\s*:", ln):
                ln = ""
            rows.append((n, ln))
    return rows


def load(path):
    """Return [(page_number, line), ...] with the original page numbers kept.

    Pages are skipped, never deleted, so a reported page number is the page of
    the PDF and can be looked at directly.
    """
    if path.lower().endswith(".pdf"):
        return load_pdf(path)
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
    gl = ["" if x is None else x for x in (x for _, x in rows)]
    breaks = {k for k, (_, x) in enumerate(rows) if x is None}

    # Scale to the measure of the text actually being read, not to a constant.
    measure = measure_of([x for x in gl if x.strip()]) or 95
    # Part D: "Tails of a third to half a line are not what the criterion is
    # about.  Clear the one- and two-word tails; decide deliberately whether
    # the rest are worth the prose they would cost."  So flag below a third,
    # and list the third-to-half band separately rather than silently.
    tail_cut = int(0.34 * measure)
    soft_cut = int(0.50 * measure)
    full_cut = int(0.72 * measure)

    ends = []
    for k, x in enumerate(gl):
        if not prose(x):
            continue
        txt = x.strip()
        if len(txt) >= soft_cut:
            continue
        if k in breaks or (k - 1) in breaks:
            continue                 # a column edge is not a paragraph end
        if k == 0 or not gl[k - 1].strip() or HEAD.match(gl[k - 1] or ""):
            continue
        prv = gl[k - 1]
        if not (prose(prv) and len(prv.rstrip()) >= full_cut):
            continue
        if re.match(r"^(doi:|arXiv|https?:|10\.\d{4})", txt) or "doi.org" in txt:
            continue                 # a reference tail is not prose
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

    print("(13) DANGLING WORDS -- whole document")
    print("      measure %d chars, tail cut %d, full-line cut %d\n"
          % (measure, tail_cut, full_cut))
    hard = [e for e in ends if len(e[1]) < tail_cut]
    soft = [e for e in ends if len(e[1]) >= tail_cut]
    print("  (i) paragraph-end dangles to clear: %d" % len(hard))
    for p, t, prv in hard:
        print("      p%-4d %-46r after: ...%s" % (p, t, prv))
    print("\n  (i) third-to-half a line, decide deliberately: %d" % len(soft))
    for p, t, prv in soft:
        print("      p%-4d %-46r" % (p, t))
    print("\n  (ii) stranded lead-in before a display: %d" % len(eqs))
    for p, t, e in eqs:
        print("      p%-4d %-46r eq: %s" % (p, t, e))
    total = len(hard) + len(eqs)
    print("\n  %d total -> %s" % (total, "CLEAN" if total == 0 else "FIX AND RERUN"))
    return total


if __name__ == "__main__":
    sys.exit(1 if main(sys.argv[1]) else 0)
