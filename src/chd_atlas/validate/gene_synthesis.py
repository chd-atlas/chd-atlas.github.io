# src/chd_atlas/validate/gene_synthesis.py
"""SYN001-SYN004: a gene's "In short" synthesis against its own file and the corpus.

Location format: the file's path with ``:synthesis`` appended, plus
``.sentences[<i>]`` where one sentence is at fault.

- SYN001 (error): a sentence cites an id that is neither an approved finding
  (``approved_by`` set) of the same file nor a lesion assertion whose ``gene``
  equals the file's gene. One issue per bad cite.
- SYN002 (error): a sentence's text holds a number, a standalone digit run or a
  word from `NUMBER_WORDS`. One issue per sentence, naming the first token found
  (a digit run is checked before the word list).
- SYN003 (error): more than `MAX_WORDS` words across all sentences; a word is a
  whitespace-separated token containing a letter or a digit.
- SYN004 (warning): the synthesis is approved and its ``approved_digest`` differs
  from `synthesis_digest` recomputed over the synthesis, the file's findings and
  the gene's own lesion assertions. An unapproved synthesis raises nothing: it is
  not published.

Pure functions over the loaded corpus. The model enforces sentence count and
approval-field arity; this module checks what a model cannot see.
"""

from __future__ import annotations

import re

from chd_atlas.corpus import Corpus
from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.models.gene_literature import synthesis_digest

MAX_WORDS = 120

NUMBER_WORDS = [
    "one",
    "two",
    "three",
    "four",
    "five",
    "six",
    "seven",
    "eight",
    "nine",
    "ten",
    "eleven",
    "twelve",
    "thirteen",
    "fourteen",
    "fifteen",
    "sixteen",
    "seventeen",
    "eighteen",
    "nineteen",
    "twenty",
    "dozen",
    "hundred",
    "half",
    "majority",
    "most",
    "dozens",
    "hundreds",
    "halves",
    "majorities",
]

# A digit run not glued to a letter, digit, dot, hyphen or colon on its left,
# nor to a letter or digit on its right: NKX2-5, 22q11.2, p.Gly80Arg, c.2T>C and
# HP:0001631 pass; "3-generation", "in 2003", ">50" and "p>0.05" do not.
_DIGITS = re.compile(r"(?<![A-Za-z0-9.\-:])\d+(?:[.,]\d+)*(?![A-Za-z0-9])")
_WORDS = re.compile(r"(?<![A-Za-z0-9])(" + "|".join(NUMBER_WORDS) + r")(?![A-Za-z0-9])", re.I)
_WORDLIKE = re.compile(r"[A-Za-z0-9]")


def _first_number(text: str) -> str | None:
    digits = _DIGITS.search(text)
    if digits:
        return digits.group(0)
    word = _WORDS.search(text)
    return word.group(0) if word else None


def validate_gene_synthesis(corpus: Corpus) -> list[ValidationIssue]:
    """Check every literature file's synthesis (SYN001-SYN004)."""
    issues: list[ValidationIssue] = []
    for path, literature in corpus.literature:
        synthesis = literature.synthesis
        if synthesis is None:
            continue
        own_assertions = [a for a in corpus.assertions if a.gene == literature.gene]
        citable: set[str] = {f.id for f in literature.findings if f.approved_by is not None}
        citable |= {a.id for a in own_assertions}

        words = 0
        for index, sentence in enumerate(synthesis.sentences):
            where = f"{path}:synthesis.sentences[{index}]"
            for cite in sentence.cites:
                if cite not in citable:
                    issues.append(
                        ValidationIssue(
                            "SYN001",
                            Severity.ERROR,
                            where,
                            f"cites {cite}, which is not an approved finding of this file "
                            f"or a lesion assertion of {literature.gene}",
                        )
                    )
            token = _first_number(sentence.text)
            if token is not None:
                issues.append(
                    ValidationIssue(
                        "SYN002",
                        Severity.ERROR,
                        where,
                        f"states a number: {token!r}; counts belong in the card's "
                        "literature line, not the narrative",
                    )
                )
            words += sum(1 for t in sentence.text.split() if _WORDLIKE.search(t))

        if words > MAX_WORDS:
            issues.append(
                ValidationIssue(
                    "SYN003",
                    Severity.ERROR,
                    f"{path}:synthesis",
                    f"{words} words; the limit is {MAX_WORDS}",
                )
            )
        if synthesis.approved_by is not None and synthesis.approved_digest != synthesis_digest(
            synthesis, literature.findings, own_assertions
        ):
            issues.append(
                ValidationIssue(
                    "SYN004",
                    Severity.WARNING,
                    f"{path}:synthesis",
                    "approval is stale: the text, its cites or drafted_by, the approved "
                    "findings (content, review mark or approval) or a lesion assertion of "
                    f"{literature.gene} changed since it was approved; re-approve with the "
                    f"command `chd-atlas literature review {path.stem}` prints",
                )
            )
    return issues
