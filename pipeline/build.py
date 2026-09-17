#!/usr/bin/env python3
"""Build the static data files the web site reads.

  data/groups.json   one record per torsion group: its curves in each of the three classes (and the
                     undecided ones), the "infinitely many?" answers for the simple and the split
                     side (from pipeline/knowledge.py), and counts;
  data/curves.json   every accepted curve (the certificates in data/curves/, verbatim);
  data/sources.json  the bibliography.

Usage:  python3 pipeline/build.py [--check]     (--check: fail if the outputs would change)
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import knowledge  # noqa: E402
from common import (CLASSES, CURVES_DIR, DATA, REJECTED_DIR, UNDECIDED, group_bracket, group_label,  # noqa: E402
                    group_order, group_sort_key, key_to_group, poly_to_string, read_json, write_json)


def equation(f, h) -> dict:
    """Magma-style text of y^2 + h(x) y = f(x): the h-term is parenthesised only when h has several terms."""
    fs = poly_to_string(f)
    hs = poly_to_string(h)
    if hs == "0":
        text = f"y^2 = {fs}"
    elif hs == "1":
        text = f"y^2 + y = {fs}"
    elif " + " in hs or " - " in hs:
        text = f"y^2 + ({hs})*y = {fs}"
    else:
        text = f"y^2 + {hs}*y = {fs}"
    return {"f": fs, "h": hs, "text": text}


def cert_summary(c: dict) -> str:
    if c["simplicity"]["geometrically_simple"]:
        return f"geometrically simple: Frobenius at p = {c['simplicity']['strict_prime']} is strict (χ = {c['simplicity']['chi']})"
    parts = []
    if c["split"]["certificate"]:
        parts.append(c["split"]["certificate"]["detail"])
    if c["q_simple"]["certified"]:
        parts.append(f"simple over ℚ: χ_{c['q_simple']['prime']} = {c['q_simple']['chi']} is irreducible")
    if not parts:
        sig = c["split"]["signature"]
        parts.append("no certificate; " + ("the split signature holds (no good prime below 200 is strict)" if sig["all_good_primes_fail_strictness"] else "inconclusive Frobenius data"))
    return "; ".join(parts)


def curve_summary(c: dict) -> dict:
    cond = c["curve"].get("conductor") or ""
    return {
        "id": c["id"], "class": c["class"], "class_certified": c["class_certified"], "status": c["status"],
        "class_text": c["class_text"],
        "f": c["curve"]["f"], "h": c["curve"]["h"], "equation": equation(c["curve"]["f"], c["curve"]["h"]),
        "conductor": cond, "conductor_sort": int(cond) if cond.isdigit() else None,
        "conductor_note": c["curve"].get("conductor_note", ""),
        "discriminant_digits": c["curve"].get("discriminant_digits"),
        "lmfdb": c["lmfdb"], "sources": c.get("sources", []), "rm": c.get("rm", False),
        "discovered_by": c["discovery"]["by"], "year": c["discovery"]["year"], "reference": c.get("reference", ""),
        "new_group": c.get("new_group", False), "new_for_class": c.get("new_for_class", False),
        "certificate": cert_summary(c), "verified": c["dates"]["verified"], "submitted": c["dates"]["submitted"],
        "source_kind": (c.get("source") or {}).get("kind", ""),
    }


def build():
    curves = [read_json(p) for p in sorted(CURVES_DIR.glob("*.json"))]
    curves.sort(key=lambda c: (group_sort_key(c["group"]), c["id"]))
    by_key: dict[str, list] = {}
    for c in curves:
        by_key.setdefault(c["group_key"], []).append(c)
    keys = sorted(by_key, key=lambda k: group_sort_key(key_to_group(k)))
    groups = []
    for key in keys:
        inv = key_to_group(key)
        mine = by_key[key]
        classes = {cls: [] for cls in list(CLASSES) + ["undecided"]}
        for c in mine:
            cls = c["class"] if c["class"] in CLASSES else "undecided"
            classes[cls].append(curve_summary(c))
        for cls in classes:
            classes[cls].sort(key=lambda s: (s["conductor_sort"] is None, s["conductor_sort"] or 0, s["id"]))
        rec = {
            "key": key, "group": inv, "bracket": group_bracket(inv), "label": group_label(inv),
            "order": group_order(inv), "rank": len(inv),
            "classes": classes,
            "known": {cls: bool(classes[cls]) for cls in classes},
            "n_curves": len(mine),
            "infinite": {"simple": knowledge.infinite_record("simple", key), "split": knowledge.infinite_record("split", key)},
        }
        groups.append(rec)
    now = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    counts = {cls: sum(1 for g in groups if g["known"][cls]) for cls in list(CLASSES) + ["undecided"]}
    groups_out = {
        "generated": now,
        "conventions": {
            "group": "J(Q)_tors as a list of invariant factors [n1, ..., nr], n1 | ... | nr; [] is the trivial group; key = factors joined by '.', '1' for the trivial group",
            "classes": dict(CLASSES), "undecided": dict(UNDECIDED),
            "infinite": {"grades": knowledge.GRADE_TEXT, "labels": knowledge.GRADE_LABEL},
            "order": "groups are listed as in the paper's tables: by number of invariant factors, then lexicographically",
        },
        "n_groups": len(groups),
        "n_curves": len(curves),
        "n_certified": sum(1 for c in curves if c["status"] == "certified"),
        "n_verified": sum(1 for c in curves if c["status"] == "verified"),
        "n_groups_by_class": counts,
        "n_groups_any_split": sum(1 for g in groups if g["known"]["qsplit"] or g["known"]["gsplit"]),
        "n_infinite": {side: {grade: sum(1 for g in groups if g["infinite"][side]["grade"] == grade) for grade in ("exact", "family", "open")}
                       for side in ("simple", "split")},
        "n_rejected": len(list(REJECTED_DIR.glob("*.json"))),
        "groups": groups,
    }
    curves_out = {"generated": now, "curves": curves}
    sources_out = {"generated": now, "sources": knowledge.SOURCES}
    return groups_out, curves_out, sources_out


def stable(obj):
    o = json.loads(json.dumps(obj))
    o.pop("generated", None)
    return json.dumps(o, sort_keys=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    groups_out, curves_out, sources_out = build()
    targets = [(DATA / "groups.json", groups_out), (DATA / "curves.json", curves_out), (DATA / "sources.json", sources_out)]
    if args.check:
        changed = [str(p) for p, o in targets if not p.exists() or stable(read_json(p)) != stable(o)]
        if changed:
            sys.exit("out of date: " + ", ".join(changed))
        print("data files up to date")
        return
    for p, o in targets:
        write_json(p, o)
    print(f"wrote data/groups.json ({groups_out['n_groups']} groups), data/curves.json ({groups_out['n_curves']} curves: "
          f"{groups_out['n_certified']} certified, {groups_out['n_verified']} verified), data/sources.json")


if __name__ == "__main__":
    main()
