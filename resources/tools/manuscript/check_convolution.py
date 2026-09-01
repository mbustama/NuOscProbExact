#!/usr/bin/env python3
r"""Convoluted sentences: the ones a reader has to unpick before reading.

The other checkers here catch words.  This one catches *shapes*.  Every
sentence the author flagged by hand shared one of four, and none of them
uses a single flagged word:

  1. abstract-noun predicate --- "The swing is not a matter of how the
     density inside each shell is chosen."  A negated abstract noun
     ("a matter of", "a question of") sitting over a passive nominal
     clause.  The sentence says what something is not, in terms of a
     noun that means nothing on its own.

  2. echoed head --- "A cost on the plane is the cost of one probability."
     The complement repeats the subject's head noun, so the copula
     defines a word with itself.

  3. buried verb --- "Each point on the plane is the time one probability
     takes at the setting drawn."  The real verb ("takes") is inside a
     noun phrase, and the sentence's own verb is only "is".

  4. clause chain --- "The costs on the plane are per probability at the
     settings drawn, which run far finer than a figure such as Fig. 7
     asks for, so they are not the cost of producing one."  Three
     clauses, a relative and a consecutive, and a pronoun whose referent
     is four clauses back.

Usage:  python3 check_convolution.py [main.tex]

Every hit is a candidate, not a verdict; a long sentence that reads
cleanly is not a fault.  Read each one aloud, which is the actual test.
"""

import re
import sys

ABSTRACT = (r"matter|question|function|consequence|feature|property|case|"
            r"issue|business|thing|aspect|point")
SUBORD = (r"\b(which|that|where|when|while|because|since|so|though|although|"
          r"whereas|unless|rather than|as though)\b")
COPULA = r"\b(is|are|was|were)\b"


def sentences(tex):
    """Yield (line_number, sentence) over prose, with markup neutralized."""
    src = tex.split("\n")
    for env in ("equation", "eqnarray", "align", "tabular", "lstlisting",
                "verbatim", "figure", "table"):
        tex = re.sub(r"\\begin\{%s\*?\}.*?\\end\{%s\*?\}" % (env, env),
                     " ", tex, flags=re.S)
    tex = re.sub(r"(?m)^\s*%.*$", " ", tex)
    tex = re.sub(r"\$[^$]*\$", " MATH ", tex)
    tex = re.sub(r"\\(cite|ref|equ|figu|Refs?)\*?\{[^}]*\}", " REF ", tex)
    tex = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?", " ", tex)
    tex = re.sub(r"[{}]", " ", tex)
    flat = re.sub(r"[ \t]+", " ", tex)
    for part in re.split(r"(?<=[.!?])\s+", flat):
        s = part.strip()
        if len(s.split()) < 8:
            continue
        key = " ".join(s.split()[:6])
        line = next((i + 1 for i, ln in enumerate(src) if key[:28] in ln), 0)
        yield line, s


def head_noun(phrase):
    """Crude head of a noun phrase: its last word before a preposition."""
    words = re.findall(r"[a-z]+", phrase.lower())
    stop = {"the", "a", "an", "this", "that", "these", "those", "each",
            "every", "one", "of", "on", "in", "at", "for", "to"}
    words = [w for w in words if w not in stop]
    return words[-1] if words else ""


def faults(s):
    out = []
    low = s.lower()

    # 1. abstract-noun predicate, especially negated
    m = re.search(r"%s\s+(not\s+)?(a|the)\s+(%s)\s+of\b" % (COPULA, ABSTRACT),
                  low)
    if m:
        out.append("abstract predicate: 'is %sa %s of'"
                   % ("not " if m.group(2) else "", m.group(4)))

    # 2. echoed head across the copula
    m = re.search(r"^(.{4,60}?)\s+%s\s+(the|a|an)\s+(.{4,48}?)\b" % COPULA,
                  low)
    if m:
        h1, h2 = head_noun(m.group(1)), head_noun(m.group(3))
        if h1 and h1 == h2:
            out.append("echoed head: '%s ... is the %s'" % (h1, h2))

    # 3. buried verb: copula, then a noun phrase carrying the real verb
    if re.search(r"%s\s+the\s+\w+(\s+\w+){0,3}\s+\w+\s+(takes|does|costs|"
                 r"needs|gives|makes|has|uses)\b" % COPULA, low):
        out.append("buried verb: the real verb sits inside a noun phrase")

    # 4. clause chain
    n_sub = len(re.findall(SUBORD, low))
    n_com = s.count(",")
    if n_sub >= 3 and n_com >= 2 and len(s.split()) >= 26:
        out.append("clause chain: %d subordinators, %d commas, %d words"
                   % (n_sub, n_com, len(s.split())))

    # 5. of-chains
    if len(re.findall(r"\bof\b", low)) >= 4 and len(s.split()) <= 45:
        out.append("of-chain: %d 'of' in %d words"
                   % (len(re.findall(r"\bof\b", low)), len(s.split())))

    return out


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "main.tex"
    tex = open(path, encoding="utf-8").read()
    total = 0
    for line, s in sentences(tex):
        f = faults(s)
        if not f:
            continue
        total += 1
        print("\n  line %-5d %s" % (line, "; ".join(f)))
        print("      %s" % (s if len(s) < 300 else s[:297] + "..."))
    print("\n  %d candidate%s" % (total, "" if total == 1 else "s"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
