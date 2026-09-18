# Verification pipeline (runs on Mordell)

```
submissions/inbox/*.json  ──verify.py──►  data/curves/<id>.json   (+ data/logs/<id>.log, <id>.conductor.log)
        ▲                                 data/rejected/<name>.json
fetch_issues.py (GitHub issues)           build.py ──► data/groups.json, data/curves.json, data/sources.json
import_paper.py (seed tables)
```

| file | role |
|---|---|
| `run_cycle.sh` | the whole cycle: `git pull` → `fetch_issues.py` → `verify.py` → `build.py` → commit → push. Uses a lock file; cron-able. |
| `fetch_issues.py` | turns open GitHub issues labelled `submission` into `submissions/inbox/issue-<n>.json`. |
| `verify.py` | validates a submission, runs `magma/verify_lib.m` and `magma/conductor_lib.m` under time limits, applies the acceptance rule (`acceptance()`: new (group, class), or smaller conductor than every census curve with that group and class; `"historical": true` in an inbox file bypasses it), writes the certificate, comments on and closes the issue. `--require-claim` rejects a curve whose computed group differs from the claimed one (used for validation runs). |
| `build.py` | rebuilds the JSON files the site reads from `data/curves/` and `knowledge.py`. `--check` fails if they are stale (used by CI). |
| `knowledge.py` | curated, cited facts: the bibliography, and the "infinitely many?" record for every group on the simple and the split side. |
| `import_paper.py` | converts the paper's `table1.txt` / `table2.txt` (in `data/knowledge/sources/`) into seed submissions; `transcribe_covers.py` extracted the 18 explicit covers of the paper's `verify_split_certificates.m` into `Genus2Torsion_covers.json`. |
| `magma/verify_lib.m` | the Magma checks (see its header): torsion subgroup, strict-prime certificate of geometric simplicity, irreducible-χ certificate of ℚ-simplicity, splitness certificates (involutions, Richelot isogenies over ℚ and quadratic fields, geometric automorphism group, submitted covers), the split signature, invariants, ℚ-isomorphism with census curves. |
| `magma/conductor_lib.m` | reduced minimal model, factored discriminant, conductor by Magma's `Conductor` (separate job: needs the discriminant factored; slow for split Jacobians whose curve has bad reduction where the Jacobian does not). Its exponent at 2 uses Ogg's formula (not guaranteed); when `verify_lib.m` finds an isogenous product E₁ × E₂ over ℚ, `verify.py` uses N(E₁)N(E₂) instead (rigorous) and records the disagreement. |
| `magma/family_2_certificate.m`, `magma/family_22212_certificate.m` | the exact-infinitude certificates computed for this site (logs in `data/knowledge/sources/`). |
| `magma/prepare_submission.m` | for submitters: run locally, prints the submission JSON (curve, torsion group + generators on the even model, class guess, conductor). |
| `magma/trivial_torsion_family.m` | the family certificate for the trivial group (log in `data/knowledge/sources/`). |
| `tests/` | `test_pipeline.py` (Python contract tests, run by CI), `run_magma_tests.sh` + `fixtures/` (known positives and negative controls for the Magma verifier; run by hand). |

## Running by hand

```sh
cd /home/fnajman/g2torsion_web
python3 -B pipeline/fetch_issues.py --dry-run          # what would be queued
pipeline/tests/run_magma_tests.sh                       # Magma smoke test in a scratch directory (~1 min)
python3 -B pipeline/verify.py                           # process the inbox, update issues
python3 -B pipeline/build.py                            # rebuild site data
python3 -B pipeline/tests/test_pipeline.py              # contract + certificate consistency tests
pipeline/run_cycle.sh                                   # everything, then push
```

Magma is always started with stdin from `/dev/null` (otherwise a script error leaves it waiting) and its
exit status is ignored (it is 0 even after errors): `verify.py` trusts only the JSON that the Magma
library writes and the `VERIFY_DONE` / `CONDUCTOR_DONE` markers in the logs.  One Magma job runs at a
time, under `--timeout` (default 3600 s) and `--mem-gb` (16); the conductor job under
`--conductor-timeout` (600 s) and 8 GB.

Every change to `magma/verify_lib.m` changes the SHA-256 recorded in new certificates; old
certificates keep the hash of the code that produced them.

## Adding knowledge

Facts about infinitely many realisations live in `knowledge.py` (`_SIMPLE_EXACT`, `_SIMPLE_FAMILY`,
`_SIMPLE_OPEN_NOTES`, `_SPLIT_FAMILY`, `_SPLIT_OPEN_NOTES`) with their sources in `SOURCES`; copies of
the source tables go to `data/knowledge/sources/` with a line in `PROVENANCE.md`.  After editing, run
`build.py` and the tests.

## Seed import

The 150 curves of the paper's tables were imported with `import_paper.py` and verified with
`verify.py --require-claim --no-github --no-lmfdb` (log: `data/knowledge/seed_import.log`): every row
reproduced the paper's torsion group, the Table 1 rows were certified geometrically simple, and the
Table 2 rows geometrically split (with the class over ℚ decided by the certificates found here).
LMFDB labels of the seed rows are those of the tables (production labels are permalinks, extended-
database labels are linked by equation); the LMFDB rate-limits automated lookups, so `verify.py`
attempts one lookup only for submissions without a label.
