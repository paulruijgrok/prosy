# Anonymisation: what is synthetic, and how to purge the history

This repo is public. It carries **no real binder sequences, no real designs,
and no internal project or target names.**

## What was replaced

| Was | Now |
|---|---|
| 4 real VHH parents in the tests | Synthetic VHHs: germline-derived framework, randomised CDRs. CDR1 = 26-33, CDR2 = 51-57, CDR3 = 96-110 *by construction*, so the annotator is tested against known ground truth. |
| A real 1035 bp rejected design in `tests/test_synthesis.py` | A synthetic 333 aa protein with de-novo-like composition, codon-optimized the naive way. It reproduces the same failure mode at almost the same magnitude (70.2% repeat coverage in both cases). |
| Lab spreadsheet filename, internal antibody-set names, design-run codenames | `plasmid_database.xlsx`, `SetA`/`SetB`, `design_group_a…f`. |
| Dated run-folder paths naming internal projects | `data/runs/…`, which is gitignored. |

Destination vectors (`data/plasmids/FP01.fa`) were **kept** — a cloning vector
is not confidential design output.

The synthetic examples live in `data/examples/` and are what the docs and tests
use, so both run from a fresh clone.

## Staying clean

```bash
python tools/check_no_lab_data.py        # tracked files; also runs in the suite
```

`tests/test_no_lab_data.py` runs this on every `pytest`, so a real sequence
pasted into a test or doc fails the build. Keep real inputs in `data/runs/`.

## Purging what was already pushed

Sanitising the working tree does **nothing** about history: the real sequences
were present from the initial commit and are retrievable from any clone or from
GitHub's cached views. Check:

```bash
python tools/check_no_lab_data.py --history
```

If that reports hits, rewrite the history. This is destructive and rewrites
every commit SHA — coordinate with anyone holding a clone.

```bash
pip install git-filter-repo

# 1. Build the redaction list from the real material. This file contains the
#    confidential strings, so keep it outside the repo and delete it after.
python tools/build_redactions.py > /tmp/redactions.txt

# 2. Rewrite every commit. filter-repo refuses to run on a repo with a remote
#    unless --force, and removes the remote afterwards by design.
git filter-repo --replace-text /tmp/redactions.txt --force

# 3. Verify, then clean up.
python tools/check_no_lab_data.py --history
rm /tmp/redactions.txt

# 4. Re-add the remote and overwrite the published history.
git remote add origin git@github.com:<you>/prosy.git
git push --force --all
git push --force --tags
```

Afterwards:

- **Ask GitHub Support to purge cached views.** Force-pushing does not remove
  blobs already fetched by the GitHub API or shown in cached diff views; old
  commit SHAs can stay reachable by direct URL until GC runs on their side.
- Any existing fork or clone still holds the original objects.
- If the data is genuinely sensitive, treat it as disclosed and rotate
  accordingly rather than relying on the rewrite alone.

A cleaner alternative, if the history itself is not precious: start a fresh
repo from the sanitised tree and archive the original as private.
