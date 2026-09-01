# Manuscript audit criteria

A section-by-section standard for `resources/paper/main.tex`. **Apply every point,
every time. Do not choose among them.** Running a subset and reporting the section
done is the characteristic failure; when the author says "audit," all of it applies,
to captions exactly as to body text.

**Order:** everything else first, then criterion 13 last, iterated to a clean report.

Three sources are merged here. Criteria 1 through 19 and the standing rules come
from `NuOscProbReview/reports/audit-criteria.md`, the oldest and most general list.
Criteria 20 through 23 were added later in `magnus/resources/paper/audit-criteria.md`;
they are stricter, and 21 and 22 explicitly supersede earlier points. The project
rules in part F come from this paper's own audit handovers
(`~/Research/NuOscProb/HANDOVER-paper-audit-5.md` §2) and were each corrected at
least once by the author. **Where the sources disagree, the later one wins** — two
cases, both marked below.

---

## A. Correctness

**1. Mathematical correctness.** Recompute every number before touching the prose.
Not re-read — re-run. A restated claim triggers a fresh measurement; nothing is
grandfathered.

**2. Textual correctness.** Every quotation checked against the source, including
appendices, notes, and back matter. Absence from the pages you happened to open is
not absence from the source.

**3. Clarity.**

**4. No cryptic, confusing, obscure, or undefined words or passages.** Never make
the reader decode. Compressing an interpretation into a phrase is exactly how it
turns cryptic. If saying what you mean takes a full sentence with its terms
glossed, use one. Fix the section, not the sentence.

**5. Omissions.**

**6. Leftovers.** Stale numbers from regenerated data, superseded wording, dead
cross-references.

**7. Internal consistency within the section.**

**8. Consistency against every other section**, including tables, captions, and
appendices. Counts stated in prose must match what the figure draws.

## B. Language

**9. American English.**

**10. Condense** wherever it costs no clarity, correctness, or impact. The author's
standing gloss for this paper: edits are reader-first, and clarity and correctness
are never traded for length.

**11. No one-sentence paragraphs.** Merge with the previous or next. Exception: a
final paragraph of a section that summarizes it or prefigures the next.

**12. No very long paragraphs.** Break at an appropriate place, at the seam in the
argument. **Measure from the source file, where a paragraph is blank line to blank
line** — not from extracted PDF text, where a displayed equation splits one
paragraph into several blocks and the length check never sees it. Threshold 220
words.

**13. No dangling words on a line of their own.** Part D below; run last.

**14. Never put the pointer before the name.** No demonstrative whose noun has not
yet appeared.

**15. No ", and" tic.** Never weld two independent clauses with ", and"; never open
a sentence with "And". This does not touch the serial comma, which is mandatory —
in practice fewer than one in a hundred ", and" occurrences is the tic, so classify
by hand rather than by count.

For this paper the rule is stricter than the review's, and was escalated on
2026-08-13 from "use sparingly" to banned: the construction `<complete statement>,
and <clause that restates or draws out the consequence>` is not to be used, in the
paper or in what is written to the author. Ordinary coordination of two independent
actions is fine; the tic is the trailing clarification. Repair with a semicolon, a
full stop, an em dash, or by dropping the comma when it is really a compound
predicate — and vary the repair, because all-semicolons becomes its own tic.

**16. Style, read off the author's own published papers.** Measure, don't assume:

  a. *Paragraphs*: three to six sentences, 90–150 words, opening on a declarative
     topic sentence about four times in five.
  b. *Voice*: first person through the body, not only the introduction — roughly
     one "we" per two or three paragraphs, three to five per thousand words,
     densest in methods and results. Plain active verbs. Measured over 57,000 words
     of the author's arXiv papers, the genre matters: his physics-results letters run
     at 8.6 first-person per 1000 words, this code paper at 3.6. Never "note that".
  c. *Equations*: never orphaned. A lead-in clause runs grammatically into the
     display, and a "where …" clause defines every new symbol straight after.
  d. *Numbers*: units always explicit, hedging quantified rather than linguistic
     ("3 sigma", not "seems likely"). **Significant figures: this paper overrides
     the review.** The review asks for two or three; this paper wants **one, marked
     approximate** — "about", "some" — with exceptions only where the second digit
     carries the meaning (the 0.79%/0.75% mass-defect factors).
  e. *Captions*: two to five sentences, self-contained, panels and line styles
     named. **The body must interpret, never merely describe** — a caption that
     only lists what is plotted has failed. The headline must name what the plotted
     quantity *is*, in the reader's vocabulary, and claim only what the axes show.
     No fitted slopes or mantissas anywhere; scale statements instead. Captions
     must stand alone, because floats land before the text that defines their terms.
     Measured structure: a noun-phrase headline naming the plotted quantity and
     terminated by a period --- "Flavor content of the three active mass eigenstates."
     --- never an interpretation; then `{\it Top}:` / `{\it Bottom}:` / `\emph{Left}:`;
     then the configuration; then what the curves are, with the visual encoding in
     parentheses ("(black stars)", "(dashed)"); closing on a pointer.
  f. *Apparatus*: no footnotes; enumeration inline as "First … Second …" rather
     than bulleted.
  g. *Signposting*: forward references frequent, one or two per major section.
  h. *Section titles*: short, task-shaped, often gerund-led.
  i. *Habits*: em dashes for asides; italics for variables and defined terms,
     rarely for emphasis; parentheses for caveats. Conclusions restate the
     achievement and draw implications; they do not re-report results. The author
     reaches for parentheses where a drafter reaches for em dashes, so converting an
     aside from dashes to parentheses fixes both counts at once. This paper's own
     measured rates, which are the house baseline the Magnus paper is held to:
     first person 3.6, em dashes 4.6, parentheses 16.6, `\ie`/`\eg` 0.6, all per
     1000 words; sentences median 18 words; paragraphs median 114 words and 4
     sentences.

**17. No self-weakening statements about the work.** Cut negatives that cost the
manuscript standing for nothing in return. A negative earns its place when it sets
up a stronger positive, or when it is a real scope limit stated as precision rather
than apology. A scope statement is precision, never a concession. For this paper:
do not downplay `NuOscProbExact`, and watch for passages that concede without
bounding the concession.

**18. Every figure and table gets a full introduction in the text**, in a paragraph
that *opens* on it. This paper's form is literal: **"Figure N shows…"**, stating the
interpretation rather than the contents. The introducing paragraph says what the
float is *for*; it does not restate the caption.

**19. Merge adjacent short paragraphs** — under about 70 words, or two sentences or
fewer. Exception: one that opens or closes a section. Watch the trap: prose
following a display with no blank line still belongs to the paragraph before the
display, so merging into it can breach criterion 12 instead.

**20. No apparatus prose, no informality, no sentence that comments on the
manuscript.** A physics paper never turns to address the reader about itself. Cut
every sentence of the kind "Now the part that a paper should not soften", "the paper
should say so", "worth saying plainly", "one point worth making", "we do not hide
it", "read it as an envelope". Each is the writer talking about the writing. Say the
physics; let the reader judge the candor. The same applies to captions, where a
headline that editorializes must be replaced by one that names the plotted quantity.

**21. The "and" tic, in both of its forms. Supersedes 15.** Never weld a statement to
a clarification with ", and" — "the two live at different scales, and the lower one
does not move" — and never weld two independent statements with it — "the densities
are equal, and the probabilities differ". Split into two sentences, subordinate one
clause to the other, or use a semicolon. Classify every ", and" by hand; keep only
the serial comma and genuine compound predicates sharing one subject.

  **Detector trap.** The handover's regex opens with `[^.;:!?]{25,180},\s+and\s+` — the
  25-character floor hides every weld with a short first clause. "One limit is real, and
  it belongs to arithmetic rather than to either family" is 17 characters before the comma
  and survived a full sweep because of it. Run the pattern a second time with the floor
  dropped to about 5, filtering appositives and parallel gerund subjects by hand.

**22. Succinct first, but never at the expense of clarity or correctness. Governs 10,
12 and 19.** Criterion 10 condenses only where it costs nothing; this one says what
"costs nothing" means. Where the choice is between a shorter sentence and one a
reader can follow without decoding, write for the reader; length is secondary. A
passage that is correct and clear at 150 words beats one that is cryptic at 60, and
a term used is a term defined, however much space the gloss takes. Condense and merge
only when the result still reads to someone meeting the material for the first time,
not to someone who already knows what it refers to. Both halves bind: always strive
for succinct text, cut the word, the clause and the sentence that carry nothing, and
never pad to fill a line — but when brevity would cost clarity or correctness,
brevity loses. The test is two questions in order: "can a first-time reader follow
it?", then "can it be shorter without breaking that?"

**23. The sentence that does no work: announce nothing, and say nothing twice.** Two
faults with one diagnosis — a sentence or clause carries no information — and one
test: delete it and see whether anything changes. If nothing does, it was never
saying anything.

  a. *The announcement instead of the claim.* Do not tell the reader that what
     follows matters; say the thing and let it matter. "Two features of the panel are
     worth naming" is the announcement; "not all of that structure is MSW resonance"
     is the claim it was introducing. The same shape hides in "the distinction is not
     academic", "it is worth being plain about", "the practical consequence is worth
     stating", "the reason is worth stating", "one convention is worth naming", and
     "we name them rather than pursue them here". Every one can be deleted with its
     content promoted into its place. This is criterion 20 turned toward the
     commonest form the fault takes, and it is stubborn: it survives one pass and
     reappears in the sentences written to repair the first pass. Sweep for it
     section by section, not only over new text. **Exception:** "worth" used
     quantitatively — what a thing buys or costs, "worth a factor of three", "worth
     1.7 in slabs" — is not this fault and must not be swept away with it.

  b. *The empty second half.* A clause appended to a sentence must add something the
     first half did not have. Three ways it fails. It restates: "a measurement rather
     than a curiosity, resting on an argument that is old". It asserts a negation
     with no content: "named wrappers, none of which is privileged". Or it is true by
     definition: "only a term that grows with the energy leaves a signature that
     varies with it" — which is what "grows with the energy" already means, and which
     replaced a real statement about the term failing to saturate where the others
     do. A tautology is worse than a gap, because it reads as an explanation and
     leaves the reader believing the point was made.

## C. What the tools cover

Ten checkers live beside this file. They report; a reader judges. Nothing here
convicts a passage on its own.

| Tool | Criteria | What it measures |
|---|---|---|
| `audit_section.py` | 9, 11, 12, 13, 14, 15/21, 16a, 16b, 16c, 16f | the mechanical half, one section at a time |
| `check_paragraphs.py` | 12, 16a, 19 | length and sentence count **from the source**, where displays do not hide a paragraph |
| `check_dangles.py` | 13 | cases (i) and (ii) across the whole PDF at once |
| `check_voice.py` | 16b | passive rate, agentive vs agentless, and "we" per 1000 words |
| `check_oxford.py` | 15, 21 | serial-comma candidates, reporting only — an earlier rewriting pass corrupted the text and had to be restored from backup |
| `check_scale_words.py` | 16d | modifiers of scale, split into QUANTIFIED (attached to a number, fine) and BARE (standing where a number should be) |
| `check_informality.py` | 4, register | A figurative, B commerce, C indefinite, D evaluative, E trailing aphorism |
| `check_llm_tells.py` | 4, 23 | LEXICAL word list; STRUCTURAL — antithesis flip, tricolon pile-up, referent-less "This", hedge stacks, restating summaries. The antithesis flip is criterion 23b; the restating summary is 23a |
| `check_weak_prose.py` | 3, 10, 22 | hedges, nominalizations, filler openings, redundancy, intensifiers, vague pronouns; per 10,000 words |
| `check_apparatus_prose.py` | 20, 18 | ratio of apparatus words to physics words per paragraph; `--routing` finds "a reader who knows X can skip to Y" |

Criteria **1–8, 10, 16d, 16e, 16g, 16h, 16i, 17, 19, 22** are read by hand, and so
is **23** in practice — `check_llm_tells.py --structural` finds the antithesis flip
and the restating summary, but neither the announcement that could be deleted nor the
clause that is true by definition. The tools say where to look; they do not say
whether the passage is wrong.

Two invocation traps, both paid for:

- `check_dangles.py` takes the **extracted text**, not the PDF. Run
  `pdftotext -layout resources/paper/main.pdf <scratchpad>/main.txt` first.
- `audit_section.py` hardcodes the review repo's path. Run it with
  `NUOSC_TEX=resources/paper/main.tex` set.

**Five more checkers exist in `NuOscProbReview/tools/manuscript` and are not
imported here**, though they apply to this paper: `check_numerals.py` (triages
every quoted number by how badly it could be wrong — this paper has ~189),
`check_refs.py` and `check_authors.py` and `refcheck_report.py` (verify all 174
bib entries field by field against INSPIRE-HEP and Crossref), and
`hunt_author_refs.py` (papers by a given author the manuscript does not cite).
`check_key_results.py` does not apply: this paper has no `\keyresult` panels. The
`*_reference_diff.py` files are specific to the review's source texts.

## D. Criterion 13 — dangling words (run LAST, iterate to clean)

**Rerun after every fix until clean** — each reflow moves the line breaks, and
fixing one dangle routinely creates another.

Three cases:
  (i)   the final sentence of a paragraph left as a few words, captions included;
  (ii)  one or two words of a lead-in alone on the line before a displayed equation;
  (iii) a sentence resuming after a display whose trailing inline formula breaks
        and strands a symbol. The detector cannot see this one — the stranded piece
        is math, not prose. Read the typeset page after any equation a sentence
        runs on from.

**Mechanics, each learned by failing:**

- **Detect from justification, not paragraph structure.** The body is justified, so
  every interior line is stretched to the full measure and a short line simply *is*
  a paragraph end. Inferring paragraph boundaries from indentation misses every
  list item. Working test: a prose line ending in sentence punctuation, shorter
  than the local measure, sitting under a full-measure prose line.
- **Never filter on word count.** A single word alone on a line is the defect, not
  an exception. One-word tails have twice survived a "clean" report because the
  detector excluded them.
- **Scale the threshold to the local measure.** Appendix and glossary lines may run
  to 125 characters where the body runs to 95; a fixed cut misses tails that are
  plainly short on the page.
- **Fix by ADDING, never by cutting.** Cutting does not absorb a tail, it
  re-randomizes it. Sized cuts have turned a 29-character tail into "do not." (7)
  and a 27-character one into "sense." (6). Adding is predictable: a tail of n
  characters plus k of new text lands at n+k, so choose k ≈ 50−n and it clears in
  one pass.
- **Measure the tail, don't guess.** Print the length of each flagged tail and size
  the fix to it. Guessing shuffles which words land last; measuring ends it.
- **Add content, never filler.** If a paragraph has nothing true left to say, leave
  the tail. Padding to satisfy a character count is worse than the dangle.
- **The detector is necessary, not sufficient.** Treat CLEAN as "nothing else
  obvious" and still read the last page of each subsection.
- **Know when to stop.** Tails of a third to half a line are not what the criterion
  is about. Clear the one- and two-word tails; decide deliberately whether the rest
  are worth the prose they would cost.

## E. Standing rules that override local judgment

- **Never fabricate a reference, author, title, year, or identifier.** Where a
  source cannot be checked, leave a visible TODO or a non-printing verify field —
  never a plausible guess. Use the publisher's or indexer's canonical BibTeX, not
  a hand-composed entry.
- **Verify numerically before asserting.** Untested structural claims collapse under
  questioning. Sweep every free parameter — seed, configuration, regime. What is
  not stable becomes a scale statement or is cut.
- **Name the measurement that would falsify the sentence, and run that one.**
  Measure where the alternatives differ, never where they coincide.
- **Confirmatory checks are the danger.** Three modes, and the third is worst:
  one case measured and the class asserted; an inherited claim restated and never
  measured; and *a different quantity measured than the prose claims*. The third
  retires the doubt while proving nothing. Before reporting any defect, confirm the
  array or script actually backs the claim being tested, that units and orientation
  are right, and that the result is stable.
- **A finding is not reportable until it is internally validated.** If any of that
  is unconfirmed, say the check was inconclusive rather than reporting a defect.
- **Highlighted claims — headlines, panels, key results — get verified harder than
  body text, never softer.**
- **Provenance, not scholasticism.** Cite sources to establish where something came
  from; do not conduct a detailed scholarly discussion inside the manuscript.
- **No apparatus prose.** Never write about the manuscript's own furniture.
- **Concrete language.** Name the thing, not the abstraction. No coinages.
- **Don't patch — integrate.** An edit must read as though written with its
  surroundings, not bolted on.
- **Editorial calls belong to the author.** Report and recommend; do not decide.
- **Never leave an audit undone** unless the author says stop or pause.

## F. This paper's own rules

Each was corrected at least once by the author. They are not a softer version of
the above; they are what this manuscript has actually been failed on.

**Clarity**

1. **No compressed phrase standing in for a mechanism that has not been stated.**
   The one that recurs most. Six in a single session: "truncation, not an error",
   "performs no refinement at all", "a chord is a *stack*", "which three checks
   identify", "the division", "the matched density". Each parsed only for a reader
   who already knew the answer. Test every sentence with: *would the author have to
   ask me to unpack this?*
2. **Never invent a term mid-paragraph.** Say what a thing is instead of naming it.
3. **One word per object.** Settled vocabulary: **slabs** = what `NuOscProbExact`
   cuts along the chord; **shells** = concentric subdivisions of PREM's major
   layers, which is what the three compiled codes count; **density jumps** = PREM's
   own discontinuities; **pieces** = the neutral term when both are meant.
4. **Explain why a claim holds, don't assert it.** "Its error falls as 1/N²" needed
   the midpoint-rule argument spelled out.
5. **No metaphor that collides with physics vocabulary.** "The field changes" for
   "the set of codes changes" was rejected on sight.
6. **No personification of the apparatus.** "The bottom panel adds a sterile state,
   and the slab product does not care" is the form to catch — it fails this and
   criterion 15 at once.

**Structure**

7. **No structure announcements.** "Two things follow", "Three checks identify it",
   "Two obstructions, of different kinds", "Read plainly, Fig. 10 says two things",
   "Two divisions run through the table". Let the sentences do it.
8. **No "it is worth …" constructions.** ("worth more than" as a comparison of
   value is fine.) Criterion 23a is the general statement of this fault and carries
   the quantitative exception.
9. **No danglers at the end of paragraphs** — the author's phrasing for criterion
   13 case (i), applied while writing rather than only at the end.
10. **Sections must deliver on their titles.** The headline requirement. §11.1 once
    promised to account for each residual and covered three of six; §11.2 was titled
    "Speed against accuracy" and read as a ranking.
11. **No empty second halves.** "X, rather than Y" where Y carries no information —
    "what accuracy costs each code, rather than what each reaches at one setting" —
    is a tic. Keep the ones where the alternative is real and a reader might have
    assumed it. Criterion 23b is the general statement, and names two further forms:
    the contentless negation and the clause that is true by definition.

**Substance**

12. **Interpretation, not description** of what a figure plainly shows.
13. **Error-budget language, not removal language.** A section explains what a
    residual is made of; it does not "remove contributions".
14. **No claim about what an experiment can or cannot resolve.** Four have been cut
    as unverifiable. Use the measured conversion instead: at 1300 km through
    3 g/cm³, over 0.6–20 GeV, a δ_CP shift of 15° moves P(νμ→νe) by about 7e-3,
    0.1° by about 5e-5, and 0.001° by about 5e-7. §2.3 (`sec:why_exact`) carries it.
15. **Quote version numbers, not dates**, wherever a version exists.
16. **Two measurement sets, never mixed.** Fig. 9 is the 150-energy scan with ħc
    unmatched; Fig. 10 is `speed_accuracy.json` with conventions matched. A number
    from one never appears in a sentence about the other.
17. **Figures ship frozen data.** Every figure reads a committed JSON. A cosmetic
    change must leave all of them byte-identical; only a config or measurement
    change may move them.
18. **Detail is welcome when it earns its place.** The author likes detail and said
    so — but cut it when it detracts. Propose cuts, don't presume them.
19. **Rework, don't patch.** When something is wrong, rewrite the passage as one
    argument.
20. **Throw out rather than salvage.** A sentence rewritten three times that still
    fails should go. Recognize it before being told.

**Mechanics**

21. Single-line LaTeX paragraphs. Send long passages as files, never fenced blocks —
    the terminal soft-wraps them.
22. Oxford commas, American English.
23. The paper mixes `---` with and without surrounding spaces; match the local
    paragraph.

## G. Figures and tables

Measured off the author's published papers, and holding except where this paper has
been told otherwise.

- **Ticks inward on all four sides.** Axis labels as *quantity, symbol* [unit], with
  powers of ten folded into the bracket.
- **Curves annotated in place** where there is room; legends only for what cannot be.
  Legend boxes carry a thin black frame. Uncertainty is a filled band.
- **No blank space left or right of the data.** First word of every label, legend and
  title capitalized. Exponents in TeX form, `$2\times10^{-13}$`, never `2e-13`.
  Two-panel figures get 1:1 panels.
- **Floats are `[t!]`.** Never `[h]`, never `[H]`. Tables are all numbered and
  captioned, and every table caption opens on a bold phrase naming what it shows.
- **Gridlines: this paper diverges.** The measured convention is no gridlines
  anywhere. Figure 10 has them, thickened on the author's instruction of 2026-08-30
  because they could not be seen. Do not remove them to satisfy the general rule.

## H. Reviewing a manuscript cold

- Independence cannot be faked by someone who has already read it. A genuine cold
  read must come from a reader who has not seen the work, and must not be told what
  was previously audited or fixed — otherwise it confirms rather than reads.
- **Adversarially verify every finding before acting.** In practice roughly half do
  not survive a skeptic instructed to refute them and to default to refuted when
  uncertain. Acting on unverified findings damages correct passages.
- Distinguish *refuted*, *unverified*, and *confirmed*, and never let the second
  pass as either of the others.
