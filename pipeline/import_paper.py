#!/usr/bin/env python3
"""Turn the paper's tables (data/knowledge/sources/table1.txt, table2.txt) into seed submissions.

Usage:  python3 pipeline/import_paper.py [--only-new]

Writes submissions/inbox/paper-t1-NNN.json (Table 1: geometrically simple) and paper-t2-NNN.json
(Table 2: geometrically split) with the curve, the group of the table, the class, the LMFDB label
of the table (a 2026-08 snapshot; re-resolved by the verifier), the citation keys, and for the
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
            sub = {
                "schema": SCHEMA_SUBMISSION,
                "coeffs": [f, h],
                "group": group,
                "class": cls,
                "lmfdb_label": "" if label == "-" else label,
                "sources": keys,
                "rm": rm,
                "discoverer": "; ".join(knowledge.SOURCES[k]["short"] for k in keys if k in knowledge.SOURCES),
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
