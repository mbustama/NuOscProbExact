#!/usr/bin/env python3
r"""Measure active vs passive voice, and first-person density, section by section.

Written 2026-08-26.  Two questions, because they are not the same one:

  PASSIVE RATE -- how often a clause is "is/was/are/were/been/being + past
  participle".  Passive is not a fault; a physics review uses it correctly when
  the object is the subject of interest ("the eigenvalues are computed from").
  It becomes one when the agent matters and is hidden, which is why agentive
  passives (those carrying "by X") are counted separately from agentless ones.

  FIRST PERSON -- point 16b of the checklist asks for three to five "we" per
  1000 words through the body, not only the introduction.  That is a separate
  measurement from passive rate: prose can be active and still impersonal.

    python check_voice.py
"""
import re
from pathlib import Path

TEX = Path(__file__).resolve().parents[2] / "paper" / "main.tex"

BE = r"(?:is|are|was|were|be|been|being|becomes?|became)"
IRREG = ("given taken made shown written seen known found held built set put kept left lost "
         "drawn brought chosen done run cut read said told thought understood meant sent "
         "spent met felt drawn grown won begun broken driven forgotten hidden").split()
PART = r"(?:[a-z]+ed|" + "|".join(IRREG) + r")"
PASSIVE = re.compile(r"\b%s\s+(?:not\s+|also\s+|already\s+|then\s+|still\s+|never\s+|only\s+)?%s\b" % (BE, PART))
AGENTIVE = re.compile(r"\b%s\s+(?:\w+\s+){0,2}%s\b[^.;]{0,36}\bby\s+(?!the same|then)" % (BE, PART))


def sections():
    """Split on \\section before any macro stripping, then clean each body."""
    raw = TEX.read_text(encoding="utf-8")
    for e in ("table", "figure", "tikzpicture", "tabular", "equation", "align", "description"):
        raw = re.sub(r'\\begin\{%s\*?\}.*?\\end\{%s\*?\}' % (e, e), ' ', raw, flags=re.S)
    raw = re.sub(r'(?m)^%.*$', '', raw)
    chunks = re.split(r'\\section\{(.*?)\}', raw)
    out = []
    for i in range(1, len(chunks) - 1, 2):
        title, body = chunks[i], chunks[i + 1]
        body = re.sub(r'\$[^$]*\$', ' MATH ', body)
        body = re.sub(r'\\(cite|ref|eqref)\{[^}]*\}', ' REF ', body)
        body = re.sub(r'\\[a-zA-Z]+\{([^}]*)\}', r'\1', body)
        body = re.sub(r'\\[a-zA-Z]+', ' ', body)
        body = re.sub(r'[{}]', ' ', body)
        out.append((re.sub(r'[\\{}]', '', title)[:34], body))
    return out


def main():
    tot_w = tot_p = tot_a = tot_we = 0
    rows = []
    for title, body in sections():
        w = len(re.findall(r"[A-Za-z][A-Za-z'\-]*", body))
        if w < 300:
            continue
        p = len(PASSIVE.findall(body))
        ag = len(AGENTIVE.findall(body))
        we = len(re.findall(r"\b(?:we|our|us)\b", body, re.I))
        rows.append((title, w, p, ag, we))
        tot_w += w
        tot_p += p
        tot_a += ag
        tot_we += we
    print("%-34s %6s %7s %7s %8s %8s" % ("section", "words", "passive", "/1k", "we/1k", "agentive"))
    for t, w, p, ag, we in rows:
        flag = "  <-- below band" if 1000.0 * we / w < 3 else ""
        print("%-34s %6d %7d %7.1f %8.1f %8d%s" % (t, w, p, 1000.0 * p / w, 1000.0 * we / w, ag, flag))
    print("-" * 78)
    print("%-34s %6d %7d %7.1f %8.1f %8d" % ("TOTAL", tot_w, tot_p, 1000.0 * tot_p / tot_w,
                                             1000.0 * tot_we / tot_w, tot_a))
    print("\nagentive passives (carrying 'by X'): %d of %d (%.0f%%)"
          % (tot_a, tot_p, 100.0 * tot_a / max(tot_p, 1)))
    print("checklist 16b target for first person: 3-5 per 1000 words")


if __name__ == "__main__":
    main()
