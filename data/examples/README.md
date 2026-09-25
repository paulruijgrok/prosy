# data/examples

Synthetic sequences so the docs and tests run from a fresh clone. **None of
these are real binders, designs or lab data.**

| File | What it is |
|---|---|
| `nb01.fa` | A 121 aa VHH: germline-derived framework with randomised CDRs. CDR1 = 26-33, CDR2 = 51-57, CDR3 = 96-110 by construction. |
| `designed_enzyme.fa` | A 333 aa protein with the residue composition typical of a de novo design - Ala/Gly/Arg/Leu-rich, which is what makes naive codon optimization repetitive. |

Real project inputs live in `data/runs/`, which is gitignored.
