#!/usr/bin/env python3
"""Transcribe the explicit genus-1 covers of verify_split_certificates.m (Genus2Torsion repository)
into data/knowledge/sources/Genus2Torsion_covers.json, the form in which import_paper.py attaches
them to the seed submissions of the COVER rows of table2.txt.

Usage:  python3 pipeline/transcribe_covers.py PATH/TO/verify_split_certificates.m

Each cover n is  < "Q", cremona_label, degree, p, q, hh, 0, 0 >  or  < "K", "", degree, p, q, hh, 0, 0 >
with coverK[n] = < K, g2, g3 >, K = NumberField(X^... ), for the map (x,y) -> (p/q, y*hh/q^2) from
y^2 = 4f + h^2 to E: v^2 = u^3 + a4 u + a6 (Weierstrass model of the Cremona curve) resp.
v^2 = 4u^3 - g2 u - g3.  Field elements are rewritten with w for K.1.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import KNOWLEDGE_DIR, RE_ELEMENT_W, RE_POLY_X_UPPER, write_json  # noqa: E402


def split_top(s: str) -> list[str]:
    """Split on commas at bracket depth 0."""
    out, depth, cur = [], 0, ""
    for c in s:
        if c in "([":
            depth += 1
        elif c in ")]":
            depth -= 1
        if c == "," and depth == 0:
            out.append(cur.strip()); cur = ""
        else:
            cur += c
    if cur.strip():
        out.append(cur.strip())
    return out


def elements(list_src: str) -> list[str]:
    """'PK![ K | a, b ]' or 'PxQ![a, b]' -> ['a', 'b'] with K.1 -> w."""
    inner = list_src[list_src.index("[") + 1: list_src.rindex("]")]
    if "|" in inner:
        inner = inner.split("|", 1)[1]
    items = [t.replace("K.1", "w").replace(" ", "") for t in split_top(inner)]
    for t in items:
        assert RE_ELEMENT_W.match(t), t
    return items


def main(path: str):
    src = Path(path).read_text()
    covers = {}
    current_field = None
    for line in src.splitlines():
        line = line.strip()
        m = re.match(r"K := NumberField\((.*)\);$", line)
        if m:
            current_field = m.group(1).replace(" ", "")
            assert RE_POLY_X_UPPER.match(current_field), current_field
            continue
        m = re.match(r"covers\[(\d+)\] := < \"(Q|K)\", \"([^\"]*)\", (\d+), (.*), PK!0, PK!0 >;$", line) or \
            re.match(r"covers\[(\d+)\] := < \"(Q|K)\", \"([^\"]*)\", (\d+), (.*), PxQ!0, PxQ!0 >;$", line)
        if m:
            n, over, label, deg, rest = m.groups()
            parts = split_top(rest)
            assert len(parts) == 3, (n, len(parts))
            p, q, hh = (elements(t) for t in parts)
            covers[int(n)] = {"over": over, "cremona": label, "degree": int(deg), "p": p, "q": q, "hh": hh,
                              "field": "" if over == "Q" else current_field}
            continue
        m = re.match(r"coverK\[(\d+)\] := < K, K!\((.*)\), K!\((.*)\) >;$", line)
        if m:
            n, g2, g3 = m.groups()
            g2 = g2.replace("K.1", "w").replace(" ", ""); g3 = g3.replace("K.1", "w").replace(" ", "")
            assert RE_ELEMENT_W.match(g2) and RE_ELEMENT_W.match(g3)
            covers[int(n)]["cubic"] = ["4", "0", "-(" + g2 + ")", "-(" + g3 + ")"]
            covers[int(n)]["note"] = "E: v^2 = 4u^3 - g2 u - g3 with g2 = " + g2 + ", g3 = " + g3
    for n, c in covers.items():
        if c["over"] == "Q":
            c["cubic"] = []   # filled from the Cremona label by the verifier: [1, 0, a4, a6] of WeierstrassModel
            c["note"] = f"degree-{c['degree']} map to the elliptic curve {c['cremona']} over Q"
        else:
            c["note"] = f"degree-{c['degree']} map to an elliptic curve over the field defined by {c['field']}; " + c["note"]
        c["source"] = "Genus2Torsion/verify_split_certificates.m, cover %d" % n
    assert len(covers) == 18 and all("cubic" in c for c in covers.values())
    out = KNOWLEDGE_DIR / "sources" / "Genus2Torsion_covers.json"
    write_json(out, {"description": "explicit genus-1 covers certifying geometric splitness, transcribed from verify_split_certificates.m of https://github.com/AndrewVSutherland/Genus2Torsion (commit 508852a); the map is (x,y) -> (p/q, y*hh/q^2) from y^2 = 4f + h^2 to v^2 = e3 u^3 + e2 u^2 + e1 u + e0; coefficients ascending; w = K.1",
                     "covers": {str(k): covers[k] for k in sorted(covers)}})
    print(f"wrote {out} ({len(covers)} covers)")


if __name__ == "__main__":
    main(sys.argv[1])
