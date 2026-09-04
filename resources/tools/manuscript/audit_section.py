#!/usr/bin/env python3
r"""Mechanical half of the audit criteria, for one section of main.tex.

The criteria are in CHECKLIST.md beside this file.  This covers the ones a
machine can see: (9) American English, (11) one-sentence paragraphs, (12)
over-long paragraphs, (13) dangling last lines in the typeset output, (14)
pointer-before-name, (15) the ", and" tic, and the measurable parts of (16).
Criteria 1-8, 10, 17, 19, 22 and most of 23 are read by hand; this only says
where to look.

Case (13) needs the typeset page, so it reads the extracted text rather than
the source.  Produce it first:

    pdftotext -layout ../../paper/main.pdf ../../paper/main.txt

Then, from this directory:

    python3 audit_section.py <substring of the \section{...} title>

NUOSC_TEX and NUOSC_TXT override either path, for auditing a manuscript that
does not live beside this script.
"""
import io
import os
import re
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_PAPER = os.path.join(_HERE, os.pardir, os.pardir, 'paper')

TEX = os.environ.get('NUOSC_TEX', os.path.join(_PAPER, 'main.tex'))
TXT = os.environ.get('NUOSC_TXT', os.path.join(_PAPER, 'main.txt'))

BRIT = [(r'\b\w+is(e|ed|es|ing|ation|ations)\b', 'BrE -ise'),
        (r'\bcolour|behaviour|favour|labour|neighbour(?!hood)', 'BrE -our'),
        (r'\bcentre|metre(?!r)|fibre|litre', 'BrE -re'),
        (r'\banalyse|catalyse|paralyse', 'BrE -yse'),
        (r'\b(model|label|travel|cancel|signal|total|fuel)l(ed|ing)\b', 'BrE -ll-'),
        (r'\bwhilst|amongst|learnt|amidst\b', 'BrE archaic'),
        (r'\bper cent\b', 'BrE "per cent" -> percent'),
        (r'\bdefence|offence|pretence|licence\b', 'BrE -ce'),
        (r'\btowards\b', 'prefer "toward"'),
        (r'\bgrey\b', 'BrE grey -> gray')]
# words that legitimately carry -ise/-yse etc. in AmE or are proper nouns
OK = re.compile(r'\b(precise|precised|concise|wise|otherwise|likewise|noise|'
                r'raise|raised|raises|rise|risen|rises|rising|promise|promised|'
                r'comprise|comprised|comprises|surprise|surprised|expertise|'
                r'exercise|exercised|revise|revised|devise|devised|advise|'
                r'supervise|arise|arises|arising|arisen|paradise|franchise|'
                r'guise|disguise|treatise|premise|premises|demise|excise|'
                r'incise|anise|chastise|despise|elise|louise|noisy|clockwise)\b',
                re.I)

def strip(t):
    t = re.sub(r'\\begin\{(tikzpicture|tabular\*?|verbatim|lstlisting)\}.*?'
               r'\\end\{\1\}', ' ', t, flags=re.S)
    # numbered/lettered lists are not paragraphs; 16f counts them separately
    t = re.sub(r'\\begin\{(enumerate|description)\}.*?\\end\{\1\}', ' LIST ', t, flags=re.S)
    t = re.sub(r'\\begin\{(equation|align|aligned|array|gather)\*?\}.*?'
               r'\\end\{\1\*?\}', ' EQ ', t, flags=re.S)
    t = re.sub(r'(?m)^\s*%.*$', '', t)
    return t

def sentences(p):
    p = re.sub(r'~?\\cite[a-z]*(\[[^\]]*\])?\{[^}]*\}', '', p)
    p = re.sub(r'\\(ref|eqref)\{[^}]*\}', 'X', p)
    p = re.sub(r'(Sec|Fig|Eq|Eqs|Secs|Figs|Tab|Refs|Ref|App|Apps|cf|i\.e|e\.g|vs|Prof|Dr)\.', r'\1<D>', p)
    p = re.sub(r'\$[^$]*\$', 'x', p)
    p = re.sub(r'\\[a-zA-Z]+\*?(\[[^\]]*\])?(\{[^{}]*\})?', ' ', p)
    p = re.sub(r'[{}]', ' ', p)
    p = re.sub(r'\s+', ' ', p).strip()
    return [s.replace('<D>', '.') for s in re.split(r'(?<=[.!?])\s+(?=[A-Z(])', p) if s.strip()]

def main(key):
    src = io.open(TEX, encoding='utf-8').read()
    ms = list(re.finditer(r'(?m)^\\section\{(.*?)\}', src))
    i = next(k for k, m in enumerate(ms) if key.lower() in m.group(1).lower())
    a = ms[i].start()
    b = ms[i+1].start() if i+1 < len(ms) else len(src)
    title = ms[i].group(1)
    body = src[a:b]
    line0 = src[:a].count('\n') + 1
    print("=" * 72)
    print("SECTION %d: %s   (main.tex lines %d-%d)"
                           % (i+1, title, line0, line0+body.count('\n')))
    print("=" * 72)

    clean = strip(body)
    paras = [p for p in re.split(r'\n\s*\n', clean) if p.strip()]
    # index of the last prose paragraph of each (sub)section, for item 11
    lastof = set()
    for k, p in enumerate(paras):
        if re.match(r'\s*\\(sub)*section', p) and k: lastof.add(k-1)
    if paras: lastof.add(len(paras)-1)

    print("\n(15) ', and' welds + sentence-initial 'And'")
    _f = re.sub(r'\\begin\{(tikzpicture|tabular\*?|verbatim|lstlisting)\}.*?\\end\{\1\}',
                ' ', body, flags=re.S)
    _f = re.sub(r'(?m)^\s*%.*$', '', _f)
    flat = re.sub(r'\s+', ' ', _f)
    hits = list(re.finditer(r', and |(?<=[.!?] )And ', flat))
    print("     %d instance(s), %.1f per 1k words"
          % (len(hits), 1000*len(hits)/max(1, len(flat.split()))))
    for m in hits:
        print("       ...%s [[%s]] %s..." % (flat[max(0, m.start()-58):m.start()],
              m.group(0).strip(), flat[m.end():m.end()+58]))

    print("\n(11) one-sentence paragraphs   (ok only if last of a (sub)section)")
    print("(12) paragraphs over 220 words")
    print("(14) paragraphs opening on a bare pointer word")
    for k, p in enumerate(paras):
        if re.match(r'\s*\\(sub)*section|\s*\\label|\s*\\caption', p.strip()): continue
        ss = sentences(p)
        w = sum(len(s.split()) for s in ss)
        if w < 12: continue
        if len(ss) == 1:
            print("     (11) %s 1 sentence, %dw: %s"
                  % ("[LAST-OK]" if k in lastof else "[FIX]    ", w, ss[0][:88]))
        if w > 220:
            print("     (12) [FIX] %d words, %d sentences: %s..." % (w, len(ss), ss[0][:70]))
        m = re.match(r'\s*(This|That|These|Those|It|They|Such|Both|Either)\b', ss[0])
        if m and not re.match(r'\s*(This|That)\s+(section|subsection|review|paper|row|table|figure|appendix)', ss[0]):
            print("     (14) [CHECK] opens on %-6s: %s" % (m.group(1), ss[0][:80]))

    print("\n(9) American English")
    found = 0
    for pat, lab in BRIT:
        for m in re.finditer(pat, re.sub(r'\s+',' ',clean), re.I):
            if OK.match(m.group(0)) or OK.search(m.group(0)): continue
            print("     %-22s %s" % (lab, re.sub(r'\s+',' ',clean)[max(0,m.start()-45):m.end()+35]))
            found += 1
    if not found: print("     clean")

    print("\n(16a) paragraphs over 6 sentences   (his band: 3-6)")
    over=[]
    for k,pp in enumerate(paras):
        if re.match(r'\s*\\(sub)*section|\s*\\(label|caption|begin|end|item)', pp.strip()): continue
        ss = sentences(pp)
        w = sum(len(x.split()) for x in ss)
        if w>=25 and len(ss)>6: over.append((len(ss),w,ss[0][:66]))
    tot=sum(1 for pp in paras if not re.match(r'\s*\\(sub)*section|\s*\\(label|caption|begin|end|item)',pp.strip())
            and sum(len(x.split()) for x in sentences(pp))>=25)
    print("     %d of %d paragraphs (%.0f%%)"%(len(over),tot,100*len(over)/max(1,tot)))
    for n_,w,t in over: print("       %2d sentences, %3dw: %s..."%(n_,w,t))

    print("\n(16b) first person   (his band: 3-5 'we' per 1000 words in the body)")
    full = re.sub(r'\\begin\{(tikzpicture|tabular\*?|verbatim|lstlisting)\}.*?\\end\{\1\}',
                  ' ', body, flags=re.S)
    full = re.sub(r'(?m)^\s*%.*$', '', full)
    we = len(re.findall(r'\b[Ww]e\b', full))
    wds = len(full.split())
    print("     %d instances, %.1f per 1k words  %s"
          %(we,1000*we/max(1,wds), "OK" if 1000*we/max(1,wds)>=3 else "<-- BELOW BAND"))

    print("\n(16c) orphaned displayed equations   (no lead-in running into them)")
    orph=0
    for m in re.finditer(r'\\begin\{(equation|align)\*?\}', body):
        before=re.sub(r'\s+',' ',body[max(0,m.start()-190):m.start()]).strip()
        end=re.search(r'\\end\{(equation|align)\*?\}',body[m.start():])
        post=re.sub(r'\s+',' ',body[m.start()+end.end():m.start()+end.end()+90]).strip() if end else ''
        if before.endswith('.') and not post.lower().startswith(('where','with','in which')):
            orph += 1
            print("       ...%s  [EQ]  %s..."%(before[-64:],post[:56]))
    if not orph: print("     clean")

    print("\n(16f) footnotes and bulleted lists in the body")
    fn = len(re.findall(r'\\footnote\{', body))
    it = len(re.findall(r'\\begin\{itemize\}', body))
    print("     %d footnote(s), %d itemize block(s)  %s"%(fn,it,"" if fn+it==0 else "<-- his six have none"))

    print("\n(17) statements that weaken the review")
    PAT=[r'[Nn]othing (in it|here|of this) is new', r'\bis not new\b', r'\bnothing new\b',
         r'\bnot original\b', r'\badds nothing\b', r'\blittle more than\b',
         r'\bno more than a\b', r'\bwe do not claim\b', r'\bmakes? no claim\b',
         r'\bis not a claim\b', r'\bmerely a\b', r'\bwe make no\b',
         r'\bnot the subject\b', r'\bcannot afford\b', r'\bwe are aware of\b',
         r'\bis standard\b', r'\bis routine\b', r'\bnot our\b', r'\bof course\b']
    n17=0
    flat17=re.sub(r'\s+',' ',clean)
    for pat in PAT:
        for m in re.finditer(pat, flat17):
            seg=flat17[max(0,m.start()-95):m.end()+95]
            print("       ...%s..." % seg)
            n17 += 1
    if not n17: print("     clean")



    print("\n(18) floats without a paragraph that opens on them")
    fl = re.findall(r"\\label\{((?:fig|tab):[^}]+)\}", body)
    nofloat = re.sub(r"\\begin\{(figure|table|sidewaysfigure|sidewaystable)\*?\}.*?"
                     r"\\end\{\1\*?\}", " ", src, flags=re.S)
    opens = set()
    OP = re.compile(r"^\s*(?:\\noindent\s*)?(?:Figures?|Tables?)~\\ref\{")
    for pp in re.split(r"\n\s*\n", nofloat):
        if OP.match(pp.strip()):
            for mm in re.finditer(r"\\ref\{((?:fig|tab):[^}]+)\}", pp.strip()[:170]):
                opens.add(mm.group(1))
    bad = [x for x in fl if x not in opens]
    print("     %d float(s) here, %d missing an opening paragraph" % (len(fl), len(bad)))
    for x in bad: print("       %s" % x)

    print("\n(13) DANGLING WORDS   [run last; rerun until clean]")
    seg = ""
    if os.path.exists(TXT):
        pg = io.open(TXT, encoding="utf-8", errors="replace").read()
        # drop the front-matter "Key results" list: its entries are wrapped
        # sentences whose short last lines are not paragraph ends in the body
        # A key-result block is set in a tinted panel: its last line is the end
        # of the panel, not a stranded prose tail.  Drop the collected list too.
        pg = "\f".join(q for q in pg.split("\f")
                        if not re.search(r"^\s*Key results\s*$", q, re.M))
        pg = re.sub(r"^\s*Key result~?\s*\d+\.\d+\s*$.*?(?=\n\s*\n)",
                    "", pg, flags=re.M | re.S)
        # a numbered heading directly above a display is a heading, not a
        # stranded lead-in
        pg = re.sub(r"^\s*\d+\.\d+(\.\d+)?\.?\s+[A-Z][^\n]{0,60}$",
                    "", pg, flags=re.M)
        hd = re.search(r"(?m)^\s*%d\.?\s+%s\s*$" % (i+1, re.escape(title)), pg) \
             or re.search(re.escape(title), pg)
        nxt = None
        if hd and i+1 < len(ms):
            nxt = re.search(r"(?m)^\s*%d\.?\s+%s\s*$"
                            % (i+2, re.escape(ms[i+1].group(1))), pg[hd.end():])
        if hd:
            seg = pg[hd.start(): hd.end()+nxt.start()] if nxt else pg[hd.start():hd.start()+40000]
    if not (os.path.exists(TXT) and seg):
        print("     (no typeset text)")
    else:
        gl = [x.rstrip() for x in seg.split("\n")]
        EQN  = re.compile(r"\(\d+\.\d+\)\s*$")
        MATH = re.compile(r"[=+\u2211\u220f\u222b\u2243\u2248\u2261\u27e8\u27e9\u00b1\u2192"
                          r"\u2206\u221a\u2264\u2265\u03b1-\u03c9\u0394\u03a8\u03a3\u03a0\u2217|]")
        HEAD = re.compile(r"^\s*(\d+(\.\d+)*\.?\s+[A-Z]|Table \d|Figure \d|Appendix)")
        GAP  = re.compile(r"\S {4,}\S")          # table column gap / equation layout
        # The body is justified: every interior line is stretched to the full
        # measure, so a short line simply IS a paragraph end.  No boundary
        # inference -- that is what missed list items and four-word tails.
        def prose(x):
            return (x.strip() and not EQN.search(x) and not HEAD.match(x)
                    and not GAP.search(x)
                    and len(x) - len(x.lstrip()) <= 20
                    and len(re.findall(r"[A-Za-z]{3,}", x)) >= 1)
        # absolute thresholds: the body measure is fixed by the class at ~108
        # columns, so no percentile is needed (and percentiles were being skewed
        # by the sideways figures).
        n_end = 0
        for k, x in enumerate(gl):
            if not prose(x): continue
            txt = x.strip()
            if len(txt) >= 48: continue
            if k == 0 or not gl[k-1].strip() or HEAD.match(gl[k-1] or ""): continue
            # a dangle sits under a FULL-measure prose line; that single test
            # excludes table bodies and figure labels without a window heuristic
            prv = gl[k-1]
            if not (prose(prv) and len(prv.rstrip()) >= 60): continue
            # a case-(i) dangle ENDS a sentence.  a lead-in to a display does not,
            # and is caught by case (ii) below instead.  this replaces a lookahead
            # that kept misfiring on the next list item's inline math.
            if not txt.endswith(('.', ';', ':', '!', '?')): continue
            print("       PARA END  %-30r after: ...%s"
                  % (txt, (gl[k-1] or "").strip()[-48:]))
            n_end += 1
        n_eq = 0
        for q, ln in enumerate(gl):
            if not EQN.search(ln): continue
            r = q - 1
            while r > 0:
                c = gl[r]
                if c.strip() and not EQN.search(c):
                    ind = len(c) - len(c.lstrip())
                    w = c.split()
                    if ind <= 10 and len(c.strip()) < 40 and re.search(r"[A-Za-z]{3}", c) and not MATH.search(c):
                        print("       BEFORE EQ %-28r  eq: %s" % (c.strip(), ln.strip()[:40]))
                        n_eq += 1
                    if ind <= 10 and re.search(r"[A-Za-z]{4}", c): break
                r -= 1
        print("     %d paragraph-end, %d before-equation  ->  %s"
              % (n_end, n_eq, "CLEAN" if n_end + n_eq == 0 else "RERUN AFTER FIXING"))

if __name__ == '__main__':
    main(' '.join(sys.argv[1:]) or 'Introduction')