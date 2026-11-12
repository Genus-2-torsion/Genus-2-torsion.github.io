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
| `fill_conductors.py` | reruns `magma/conductor_lib.m` for existing census curves under a longer time limit (default 6 h) and rewrites their minimal-model and conductor fields by the same rules as `verify.py` (`apply_conductor()`); used on 2026-10-20 for the ten curves that had no minimal model. |
| `build.py` | rebuilds the JSON files the site reads from `data/curves/` and `knowledge.py`. `--check` fails if they are stale (used by CI). |
| `knowledge.py` | curated, cited facts: the bibliography, class-specific impossibility results, and the "infinitely many?" record for every group on the simple and the split side. |
| `import_paper.py` | converts the paper's `table1.txt` / `table2.txt` (in `data/knowledge/sources/`) into seed submissions; `transcribe_covers.py` extracted the 18 explicit covers of the paper's `verify_split_certificates.m` into `Genus2Torsion_covers.json`. |
| `magma/verify_lib.m` | the Magma checks (see its header): torsion subgroup, strict-prime certificate of geometric simplicity, irreducible-χ certificate of ℚ-simplicity, splitness certificates (involutions, Richelot isogenies over ℚ and quadratic fields, geometric automorphism group, submitted covers), the split signature, invariants, ℚ-isomorphism with census curves. |
| `magma/conductor_lib.m` | reduced minimal model, factored discriminant, conductor by Magma's `Conductor` (separate job: needs the discriminant factored; slow for split Jacobians whose curve has bad reduction where the Jacobian does not). Its exponent at 2 uses Ogg's formula (not guaranteed); when `verify_lib.m` finds an isogenous product E₁ × E₂ over ℚ, `verify.py` uses N(E₁)N(E₂) instead (rigorous) and records the disagreement. | It writes the minimal model before calling `Conductor`, so a failed or timed-out conductor keeps it; `verify.py` accepts a conductor only if the job's JSON says `ok: true` (Magma carries on after an error, so the `CONDUCTOR_DONE` marker alone is not enough) and records why it is missing (time limit, or Magma's error message, e.g. `Fibre blowups over number fields not yet implemented`).
| `magma/conductor_by_prime.m` | the same job variables and JSON as `conductor_lib.m`, but the conductor is computed prime by prime, with the exponent at 2 on an isomorphic model `Model2` (default x -> 1/X): Magma's regular model at 2 can take hours on one model and seconds on another. Run by hand for a curve whose conductor job timed out (write the job in `<workdir>/<ID>/conductor.m` loading this file, run it, then `fill_conductors.py --reuse --workdir <workdir> <ID>`); used for 2.2.16.b and 2.24.b (Model2 = [0,1,1,1]) on 2026-11-10. |
| `magma/family_2_certificate.m`, `magma/family_22212_certificate.m` | the exact-infinitude certificates computed for this site (logs in `data/knowledge/sources/`). |
| `magma/prepare_submission.m` | for submitters: run locally, prints the submission JSON (curve, torsion group + generators on the even model, class guess, conductor). |
| `magma/trivial_torsion_family.m` | the family certificate for the trivial group (log in `data/knowledge/sources/`). |
| `magma/attribution_check.m` | compares disputed census curves with the curves printed in their cited sources (G2-invariants, isomorphism over ℚ); backs the credit errata of 2026-10-14 in `PROVENANCE.md` (log in `data/knowledge/sources/`). |
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

Class-specific impossibility results are returned by `exclusions_record` in `knowledge.py`.
Each exclusion has a reason, bibliography keys, and a proof URL; `build.py` publishes these
in the group's `exclusions` field and refuses to build if an exclusion contradicts a recorded
example in that class. These records are separate from the combined geometrically split
infinitude status. If both split classes are excluded, `infinite_record("split", key)` returns
`grade: "impossible"` with a proof URL, and the build also rejects a geometrically split example
whose class over ℚ is undecided. This currently applies to [13], [17], [22], [23], [26], [29], [31],
[33], [34], [39], [2,22], [2,26], [2,2,14], and [2,2,2,10]. For [13], [2,2,14], and [2,2,2,10] the
`gsplit` proofs need simplicity over ℚ (`QSIMPLE_ONLY`), so their combined notes also state the `qsplit` reason.
The exclusions for [30], [35], [60], [63], [70], [2,30], [5,5], [5,10], and [7,7] apply only to
`gsplit`; all have `qsplit` examples, so their combined infinitude records retain their existing
grades.

`QSimple228Obstruction` adds the principal-polarization bound that rational 7-torsion forces
rational 2-primary torsion of order at most four on a geometrically split surface simple over ℚ.
The predicate applies to any group with 7-torsion and order divisible by eight, including future
entries, while retaining earlier proof sources for already excluded groups. It newly marks
[2,28] impossible only in `gsplit`; its `qsplit` case and combined split grade stay open.
The proof is in `proofs/qsimple-228-obstruction.html`, with PDF/LaTeX and checks under
`data/knowledge/sources/qsimple_228_*`. Examples beyond the census, such as [56], are documented
on `exclusions.html` rather than added as curve records.

`QSimple2210Obstruction` (11 November 2026) proves that rational 5-torsion forces rational 2-rank at
most 2 on a principally polarized, geometrically split abelian surface simple over ℚ. The predicate
`order % 5 == 0 and rank(2) >= 3` applies to [2,2,10], [2,2,20] and [2,2,2,10] in the census and to
future entries such as [2,2,40]; groups that already had a `gsplit` exclusion (such as [2,2,30]) keep
their earlier proof. With the `qsplit` 2-rank obstruction, [2,2,2,10] becomes impossible on the
whole split side; [2,2,10] (with `qsplit` examples) and [2,2,20] (open in `qsplit`) keep their combined
grades. The proof is in `proofs/qsimple-2210-obstruction.html`, with PDF/LaTeX, the Magma certificate
and a supplementary validation under `data/knowledge/sources/qsimple_2210_*`.

`QSPLIT_ORDER256_EXCLUSIONS` records six exact groups excluded only in `qsplit`, even though
they have no census entries: [4,64], [8,32], [16,16], [2,2,64], [2,4,32], and [2,8,16].
They appear separately on `exclusions.html`; no artificial curve or census row is created.
The same `exclusions_record` check rejects any future contradictory `qsplit` entry.
Their combined split infinitude grades are unchanged. The proof and its downloadable
PDF/LaTeX are at `proofs/qsplit-order-256-obstruction.html`.

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
