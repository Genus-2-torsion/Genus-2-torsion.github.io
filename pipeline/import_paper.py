#!/usr/bin/env python3
"""Turn the paper's tables (data/knowledge/sources/table1.txt, table2.txt) into seed submissions.

Usage:  python3 pipeline/import_paper.py [--only-new]

Writes submissions/inbox/paper-t1-NNN.json (Table 1: geometrically simple) and paper-t2-NNN.json
(Table 2: geometrically split) with the curve, the group of the table, the class, the LMFDB label
of the table (a 2026-08 snapshot; re-resolved by the verifier), the citation keys of the Source(s)
column (which date the group) and the credit of the curve itself (see `credit`), and for the
COVER rows of Table 2 the explicit genus-1 cover transcribed from verify_split_certificates.m
(data/knowledge/sources/Genus2Torsion_covers.json).  The verifier then re-derives everything.
"""

from __future__ import annotations

import argparse
import ast
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import knowledge  # noqa: E402
from common import CURVES_DIR, KNOWLEDGE_DIR, SCHEMA_SUBMISSION, SUBMISSIONS_INBOX, SUBMISSIONS_PROCESSED, log, read_json, write_json  # noqa: E402

TABLE_URL = "https://github.com/AndrewVSutherland/Genus2Torsion/blob/508852a150e17c00b84b209fee23300df897871d/"
TABLE_DATE = "2026-08-30"

# Credit.  The Source(s) column of the paper's tables dates the *group*: it names the first known
# source of a realisation, which may give only a family whose torsion contains the group, or a
# different curve (see the "Source(s) column" paragraphs of the companion tables
# torsion_jac/paper/old/torsion_realizations.tex and split_torsion_table.tex).  It is kept, as citation
# keys, in "sources".  "discoverer" (discovery.by of the certificate) credits the displayed curve
# itself, which is usually the smallest-conductor curve of a database.  It defaults to the Source
# column; the rows below are those where the curve is known to come from elsewhere (errata of
# 2026-10-14 in data/knowledge/sources/PROVENANCE.md, check in data/knowledge/sources/attribution_check.log).
# The default is UNVERIFIED for the other rows: for most Table 1 rows the curve is a database curve
# while the column names the group's first source.
DISCOVERER_TEXT = {"BSSVY": "LMFDB genus 2 database (BSSVY 2016)",
                   "BookerSutherland": "Booker–Sutherland (extended database)"}
CURVE_CREDIT = {
    # split companion "HLP2000†, <witness>": HLP 2000 give a family with torsion containing the group
    # (containment only); the exact witness is a database curve ...
    "paper-t2-018": ["BookerSutherland"],   # [20]      363.a.3
    "paper-t2-024": ["BookerSutherland"],   # [30]      4356.j.1
    "paper-t2-025": ["BookerSutherland"],   # [35]      7436.a.2
    "paper-t2-027": ["BookerSutherland"],   # [40]      3168.i.1
    "paper-t2-030": ["BookerSutherland"],   # [60]      13500.r.4
    "paper-t2-038": ["BSSVY"],              # [2,12]    450.a.36450.1
    "paper-t2-043": ["BookerSutherland"],   # [2,24]    1440.c.2
    "paper-t2-049": ["BookerSutherland"],   # [3,12]    2700.c.1
    "paper-t2-056": ["BookerSutherland"],   # [5,10]    85500.cj.1
    "paper-t2-057": ["BSSVY"],              # [6,6]     196.a.21952.1 (= X0(28))
    "paper-t2-070": ["BookerSutherland"],   # [2,4,8]   573300.ir.1
    "paper-t2-071": ["BookerSutherland"],   # [2,6,6]   132300.eo.2
    # ... or a curve constructed in the paper ("HLP2000†, new")
    "paper-t2-028": ["BNSS2026"],           # [45]
    "paper-t2-058": ["BNSS2026"],           # [6,12]
    "paper-t2-060": ["BNSS2026"],           # [8,8]
    "paper-t2-068": ["BNSS2026"],           # [2,2,24]
    "paper-t2-077": ["BNSS2026"],           # [2,2,4,8]
    # [3,9]: companion "Leprevost1995†, BSSVY"; the curve is the production curve 2916.a.139968.1
    "paper-t2-048": ["BSSVY"],
    # [24]: 300.a.1 is neither of Leprévost's curves C24,1 and C24,2 (G2-invariants); extended database
    "paper-t2-020": ["BookerSutherland"],
}
# [27]: 7776.e.2 is not Leprévost's C27 (which is split over Q, census 27.c) but the order-27 curve
# y^2 + (6x^3+9x^2+6x-1)y = -3x^2 marked "New" in Table 2 of Howe 2015; Leprévost 1995 has no curve
# in this class (simple over Q, geometrically split), so the group-level source is replaced as well.
SOURCES_OVERRIDE = {"paper-t2-022": ["Howe2015"]}


def credit(name: str, keys: list[str]) -> tuple[list[str], str]:
    """(sources, discoverer) of the seed row `name` whose Source(s) column gave the keys `keys`."""
    sources = SOURCES_OVERRIDE.get(name, keys)
    who = CURVE_CREDIT.get(name, sources)
    return sources, "; ".join(DISCOVERER_TEXT.get(k, knowledge.SOURCES[k]["short"]) for k in who if k in knowledge.SOURCES)


def rows(fn: str):
    with open(KNOWLEDGE_DIR / "sources" / fn) as fh:
        for ln, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("|")
            assert len(parts) == 7, f"{fn}:{ln}"
            yield ln, parts


def already_done(name: str) -> bool:
    if (SUBMISSIONS_INBOX / f"{name}.json").exists() or list(SUBMISSIONS_PROCESSED.glob(f"{name}.*json")):
        return True
    for c in CURVES_DIR.glob("*.json"):
        src = read_json(c).get("source", {})
        if src.get("kind") == "import" and src.get("name") == name:
            return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only-new", action="store_true", help="skip rows already queued or verified")
    args = ap.parse_args()
    covers = read_json(KNOWLEDGE_DIR / "sources" / "Genus2Torsion_covers.json")["covers"]
    SUBMISSIONS_INBOX.mkdir(parents=True, exist_ok=True)
    n = 0
    for table, fn, cls in (("t1", "table1.txt", "simple"), ("t2", "table2.txt", "split")):
        for i, (ln, parts) in enumerate(rows(fn), 1):
            group = ast.literal_eval(parts[0])
            f, h = ast.literal_eval(parts[1])
            label = parts[2].strip()
            source_col = parts[4].strip()
            route = parts[5].strip()
            comment = parts[6].strip()
            keys, rm = knowledge.parse_source_column(source_col)
            name = f"paper-{table}-{i:03d}"
            if args.only_new and already_done(name):
                continue
            sources, discoverer = credit(name, keys)
            sub = {
                "schema": SCHEMA_SUBMISSION,
                "coeffs": [f, h],
                "group": group,
                "class": cls,
                "lmfdb_label": "" if label == "-" else label,
                "sources": sources,
                "rm": rm,
                "discoverer": discoverer,
                "reference": f"Balakrishnan–Najman–Shnidman–Sutherland 2026, Table {1 if table == 't1' else 2}, row {i}"
                             + (f" (source: {source_col})" if source_col else ""),
                "notes": "" if comment in ("-", "") else comment,
                "submitter": "seed import from the Genus2Torsion repository",
                "affiliation": "",
                "github": "",
                "date": TABLE_DATE,
                "route": route if route != "-" else "",
                "source": {"kind": "import", "name": name, "file": fn, "line": ln, "url": TABLE_URL + fn},
            }
            if route.startswith("COVER:"):
                c = covers[route.split(":")[1]]
                sub["cover"] = {"field": c["field"], "cremona": c["cremona"], "cubic": c["cubic"],
                                "p": c["p"], "q": c["q"], "hh": c["hh"], "note": c["note"] + " [" + c["source"] + "]"}
            write_json(SUBMISSIONS_INBOX / f"{name}.json", sub)
            n += 1
    log(f"queued {n} seed submissions in {SUBMISSIONS_INBOX}")


if __name__ == "__main__":
    main()
