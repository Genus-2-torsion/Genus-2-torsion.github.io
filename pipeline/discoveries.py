"""Build the public timeline from the paper's tables and certified census curves."""

import ast
import csv
import io

from common import KNOWLEDGE_DIR, group_bracket, read_json
from import_paper import rows

SITE = "https://genus-2-torsion.github.io/"
PAPER = "https://arxiv.org/abs/2608.28543v1"


def build_discoveries(curves):
    paper_groups = {tuple(ast.literal_eval(fields[0]))
                    for fn in ("table1.txt", "table2.txt") for _, fields in rows(fn)}
    credits = read_json(KNOWLEDGE_DIR / "discovery_credits.json")
    recorded_groups = {tuple(c["group"]) for c in curves}
    new_groups = {}
    def report_order(curve):
        dates = curve["dates"]
        return (not curve["new_group"], dates.get("submitted") or dates["verified"][:10], curve["id"])

    for c in sorted(curves, key=report_order):
        group = tuple(c["group"])
        if group not in paper_groups:
            new_groups.setdefault(group, c["id"])
    first_ids = set(new_groups.values())

    entries = []
    for c in curves:
        credit = credits.get(c["id"], {})
        if c["id"] in first_ids:
            scope = "New group"
        elif (credit.get("class_first") and c["class_certified"] and c["new_for_class"]
              and tuple(c["group"]) in paper_groups):
            scope = "New in class"
        elif credit.get("smaller_conductor") and c["accepted_because"].startswith("conductor "):
            scope = "Smaller conductor"
        else:
            continue
        date = c["dates"].get("submitted")
        date_basis = "Reported"
        if not date:
            date = c["dates"]["verified"][:10]
            date_basis = "Verified"
        entries.append({
            "date": date,
            "date_basis": date_basis,
            "scope": scope,
            "group": group_bracket(c["group"]),
            "jacobian_class": c["class_text"],
            "person": credit.get("person", c["discovery"]["by"]),
            "model": credit.get("model"),
            "curve_id": c["id"],
            "curve_url": SITE + "curve.html?id=" + c["id"],
            "source_url": credit.get("source_url") or (c.get("source") or {}).get("url"),
        })
        if scope == "Smaller conductor":
            entries[-1]["conductor"] = c["curve"]["conductor"]
    entries.sort(key=lambda r: (r["date"], r["source_url"] or "", r["curve_id"]))
    return {"paper": PAPER, "paper_groups": len(paper_groups), "census_groups": len(recorded_groups),
            "new_groups": len(new_groups), "discoveries": entries}


def discoveries_csv(timeline):
    output = io.StringIO(newline="")
    fields = ("date", "date_basis", "scope", "group", "jacobian_class", "person", "model",
              "curve_id", "curve_url", "source_url", "conductor")
    writer = csv.DictWriter(output, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(timeline["discoveries"])
    return output.getvalue()
