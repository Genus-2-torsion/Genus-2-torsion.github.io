#!/usr/bin/env python3
"""Contract tests of the Python layer (no Magma needed; run by CI).

  python3 pipeline/tests/test_pipeline.py
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "pipeline"))

import knowledge  # noqa: E402
import verify  # noqa: E402
from common import (CLASSES, CURVES_DIR, KNOWLEDGE_DIR, coeffs_to_strings, group_key, group_sort_key, is_invariant_factors,  # noqa: E402
                    key_to_group, parse_poly, poly_to_string, read_json)


class PolyTests(unittest.TestCase):
    def test_parse(self):
        self.assertEqual(coeffs_to_strings(parse_poly("x^5 - 8*x^4 + 16*x^3 - x")), ["0", "-1", "0", "16", "-8", "1"])
        self.assertEqual(coeffs_to_strings(parse_poly("(3/2)*x^2 - x/2 + 1")), ["1", "-1/2", "3/2"])
        self.assertEqual(coeffs_to_strings(parse_poly("0")), ["0"])
        self.assertEqual(coeffs_to_strings(parse_poly("-x^6+x")), ["0", "1", "0", "0", "0", "0", "-1"])
        for bad in ("x^2*x", "x^2 x", "x +", "2**x", "y^2", "x^2; quit", "eval(1)"):
            with self.assertRaises(ValueError, msg=bad):
                parse_poly(bad)

    def test_roundtrip(self):
        s = "x^6 - 7*x^5 + 8*x^4 + 16*x^3 - 5*x^2 - 5*x"
        self.assertEqual(poly_to_string(coeffs_to_strings(parse_poly(s))), s)


class GroupTests(unittest.TestCase):
    def test_keys(self):
        self.assertEqual(group_key([2, 2, 12]), "2.2.12")
        self.assertEqual(group_key([]), "1")
        self.assertEqual(key_to_group("2.2.12"), [2, 2, 12])
        self.assertEqual(key_to_group("1"), [])
        self.assertTrue(is_invariant_factors([2, 2, 12]) and is_invariant_factors([]) and is_invariant_factors([31]))
        self.assertFalse(is_invariant_factors([2, 3]) or is_invariant_factors([1, 2]) or is_invariant_factors([2, 2, 2, 2, 2]))
        self.assertLess(group_sort_key([40]), group_sort_key([2, 2]))


class ValidateTests(unittest.TestCase):
    def test_accepts(self):
        v = verify.validate({"f": "x^5 - 8*x^4 + 16*x^3 - x", "h": "x", "group": [2, 4], "class": "simple"})
        self.assertEqual((v["f"], v["h"], v["group"], v["class"]), (["0", "-1", "0", "16", "-8", "1"], ["0", "1"], [2, 4], "simple"))
        v = verify.validate({"coeffs": [[0, -1, 0, 16, -8, 1], [0, 1]], "group": "[2,4]", "class": "geometrically split"})
        self.assertEqual(v["class"], "split")
        v = verify.validate({"f": "x^5+1", "cover": {"cremona": "19a1", "p": [1, 2], "q": [3], "hh": [1]}})
        self.assertEqual(v["cover"]["cremona"], "19a1")

    def test_rejects(self):
        for bad in ({"f": "x^7 + 1"}, {"f": "x^5 + 1; quit"}, {"f": "x^5+1", "group": [2, 3]}, {"f": "x^5+1", "class": "cm"},
                    {"f": "x^5+1", "cover": {"p": [1], "q": [1], "hh": [1]}}, {"f": "0"}):
            with self.assertRaises(verify.Reject, msg=str(bad)):
                verify.validate(bad)


class KnowledgeTests(unittest.TestCase):
    def test_sources_exist(self):
        for side in ("simple", "split"):
            for key in list(knowledge._SIMPLE_EXACT) + list(knowledge._SIMPLE_FAMILY) + list(knowledge._SPLIT_FAMILY) + list(knowledge._SPLIT_EXACT):
                rec = knowledge.infinite_record(side, key)
                self.assertIn(rec["grade"], ("exact", "family", "open"))
                for s in rec["sources"]:
                    self.assertIn(s, knowledge.SOURCES, f"{side} {key}: unknown source {s}")

    def test_tables_covered(self):
        """Every group of the paper's Table 1 has an explicit record on the simple side (exact, family, or an
        open note), and the record's 72 nontrivial groups plus the trivial group are exactly Table 1."""
        import ast
        t1 = []
        for line in (KNOWLEDGE_DIR / "sources" / "table1.txt").read_text().splitlines():
            if line and not line.startswith("#"):
                t1.append(group_key(ast.literal_eval(line.split("|")[0])))
        self.assertEqual(len(t1), 73)
        recorded = set(knowledge._SIMPLE_EXACT) | set(knowledge._SIMPLE_FAMILY) | set(knowledge._SIMPLE_OPEN_NOTES)
        # groups first realized for simple Jacobians after the paper, with a family certified here
        after_paper = {"3.12"}
        self.assertTrue(recorded <= set(t1) | after_paper, recorded - set(t1) - after_paper)
        # 39 of the record, minus [15] (containment only) and [4,8] (no certificate), plus the trivial group and
        # the groups certified here: [2,2,2,12], [3,12], [19], and the residue-class families of 2026-10-13
        # [13], [15], [16], [17], [24], [25], [26], [2,14], [2,16], [2,2,8], [2,2,14], and [2,4,4] (2026-10-14)
        self.assertEqual(len(knowledge._SIMPLE_EXACT), 53)
        self.assertEqual(len(knowledge._SIMPLE_FAMILY), 0)
        for key in ("13", "15", "16", "17", "19", "24", "25", "26", "2.14", "2.16", "2.2.8", "2.2.14", "2.2.2.12", "2.4.4"):
            self.assertEqual(knowledge.infinite_record("simple", key)["grade"], "exact", key)
        self.assertEqual(knowledge.infinite_record("simple", "4.8")["grade"], "open")
        t2 = []
        for line in (KNOWLEDGE_DIR / "sources" / "table2.txt").read_text().splitlines():
            if line and not line.startswith("#"):
                t2.append(group_key(ast.literal_eval(line.split("|")[0])))
        self.assertEqual(len(t2), 77)
        self.assertTrue(set(knowledge._SPLIT_FAMILY) <= set(t2))
        self.assertEqual(len(knowledge._SPLIT_FAMILY) + len(knowledge._SPLIT_OPEN_NOTES), 20)   # HLP Table 1
        self.assertEqual(list(knowledge._SPLIT_FAMILY_OTHER), ["48"])                             # Howe 2015
        self.assertEqual(knowledge.infinite_record("split", "2")["grade"], "family")             # inherited
        # split side: [8,8] certified here (family_88_split_certificate.log)
        self.assertEqual(list(knowledge._SPLIT_EXACT), ["8.8"])
        self.assertEqual(knowledge.infinite_record("split", "8.8")["grade"], "exact")
        for key in ("2.4.4", "8.8"):
            side = "simple" if key == "2.4.4" else "split"
            name = "family_244_certificate.log" if key == "2.4.4" else "family_88_split_certificate.log"
            self.assertIn(name, knowledge.infinite_record(side, key)["note"])
            self.assertTrue((KNOWLEDGE_DIR / "sources" / name).exists(), name)

    def test_subgroups(self):
        self.assertTrue(knowledge.is_subgroup([16], [32]) and knowledge.is_subgroup([2, 2, 6], [2, 2, 2, 6]) and knowledge.is_subgroup([], [2]))
        self.assertFalse(knowledge.is_subgroup([3, 3, 3], [3, 3]) or knowledge.is_subgroup([4], [2, 2]) or knowledge.is_subgroup([2, 2, 2, 2, 2], [2, 2, 2, 2]))

    def test_cite_parsing(self):
        self.assertEqual(knowledge.parse_source_column(r"\cite{BookerSutherland}, new"), (["BookerSutherland", "BNSS2026"], False))
        self.assertEqual(knowledge.parse_source_column(r"\cite{costa}$^{\text{RM}}$"), (["costa"], True))


class AcceptanceTests(unittest.TestCase):
    def curve(self, key, cls, cond):
        return {"group_key": key, "class": cls, "curve": {"conductor": cond}, "id": f"{key}.{cls[:1]}"}

    def test_rule(self):
        C = [self.curve("2.4", "simple", "997"), self.curve("2.4", "qsplit", "225"), self.curve("5", "simple", "")]
        acc = lambda cls, cert, inv, cond, hist=False: verify.acceptance(cls, cert, inv, cond, C, hist)[0]
        self.assertTrue(acc("simple", True, [2, 6], "5000"))          # new group
        self.assertTrue(acc("gsplit", True, [2, 4], "5000"))          # new for this class
        self.assertTrue(acc("simple", True, [2, 4], "996"))           # smaller conductor
        self.assertFalse(acc("simple", True, [2, 4], "997"))          # not smaller
        self.assertFalse(acc("simple", True, [2, 4], ""))             # conductor unknown, cannot compare
        self.assertTrue(acc("simple", True, [2, 4], "", True))        # historical flag
        self.assertTrue(acc("simple", True, [5], "12345"))            # known curve has no conductor
        self.assertFalse(acc("undecided", False, [2, 4], "10"))       # uncertified class, group known
        self.assertTrue(acc("undecided", False, [7], "10"))           # uncertified class, group new
        self.assertFalse(acc("split_undecided_over_Q", False, [2, 4], "10"))
        self.assertTrue(acc("split_undecided_over_Q", False, [5], "10"))

    def test_generators_validation(self):
        v = verify.validate({"f": "x^5+1", "generators": [["x", "0", 1], ["x^2 - 4*x + 1", "x", 2]]})
        self.assertEqual(v["generators"], [["x", "0", 1], ["x^2-4*x+1", "x", 2]])
        for bad in ([["x", "0", 3]], [["x", "0"]], [["x; quit", "0", 1]], [["x", "0", 1]] * 5):
            with self.assertRaises(verify.Reject):
                verify.validate({"f": "x^5+1", "generators": bad})


class CreditTests(unittest.TestCase):
    def test_seed_credit_matches_import(self):
        """The seed certificates carry the sources and the curve credit that import_paper.py assigns
        (the Source(s) column dates the group; discovery.by credits the curve: errata of 2026-10-14)."""
        import import_paper
        keys = {}
        for table, fn in (("t1", "table1.txt"), ("t2", "table2.txt")):
            for i, (_, parts) in enumerate(import_paper.rows(fn), 1):
                keys[f"paper-{table}-{i:03d}"] = knowledge.parse_source_column(parts[4].strip())[0]
        for name in list(import_paper.CURVE_CREDIT) + list(import_paper.SOURCES_OVERRIDE):
            self.assertIn(name, keys)
        for k in [k for v in list(import_paper.CURVE_CREDIT.values()) + list(import_paper.SOURCES_OVERRIDE.values()) for k in v]:
            self.assertIn(k, knowledge.SOURCES)
        n = 0
        for p in sorted(CURVES_DIR.glob("*.json")):
            c = read_json(p)
            src = c.get("source") or {}
            if src.get("kind") != "import":
                continue
            n += 1
            sources, who = import_paper.credit(src["name"], keys[src["name"]])
            self.assertEqual((c["sources"], c["discovery"]["by"]), (sources, who), c["id"])
        self.assertEqual(n, 150)
        # 27.b is Howe's curve, not Leprévost's C27 (which is the Q-split 27.c)
        self.assertEqual(read_json(CURVES_DIR / "27.b.json")["sources"], ["Howe2015"])
        self.assertNotIn("Leprévost", read_json(CURVES_DIR / "24.b.json")["discovery"]["by"])

    def test_credit_differs_flags_corrected_rows(self):
        """The pages show the curve's own credit next to the group's sources exactly for the seed rows
        whose curve credit was corrected (CURVE_CREDIT); everywhere else the sources credit the curve."""
        import build
        import import_paper
        flagged, corrected = set(), set()
        for p in sorted(CURVES_DIR.glob("*.json")):
            c = read_json(p)
            if build.credit_differs(c):
                flagged.add(c["id"])
            if (c.get("source") or {}).get("name") in import_paper.CURVE_CREDIT:
                corrected.add(c["id"])
        self.assertEqual(flagged, corrected)
        self.assertEqual(len(flagged), 19)
        self.assertIn("35.a", flagged)
        self.assertNotIn("27.b", flagged)   # its sources were replaced as well (SOURCES_OVERRIDE)


class CertificateTests(unittest.TestCase):
    def test_certificates_consistent(self):
        """Every certificate in data/curves/ is internally consistent."""
        for p in sorted(CURVES_DIR.glob("*.json")):
            c = read_json(p)
            self.assertEqual(c["id"], p.stem)
            self.assertEqual(c["group_key"], group_key(c["group"]))
            self.assertEqual(c["torsion"]["invariants"], c["group"])
            self.assertEqual(c["class_certified"], c["class"] in CLASSES)
            self.assertEqual(c["status"], "certified" if c["class_certified"] else "verified")
            if c["class"] == "simple":
                self.assertTrue(c["simplicity"]["geometrically_simple"] and not c["split"]["geometrically_split"])
            if c["class"] == "qsplit":
                self.assertTrue(c["split"]["split_over_Q"] and not c["q_simple"]["certified"])
                if c["split"].get("isogenous_product"):
                    self.assertEqual(c["curve"]["conductor"], c["split"]["isogenous_product"]["conductor"])
                    self.assertEqual(c["curve"]["conductor_source"], "isogenous product of elliptic curves")
            if c["class"] == "gsplit":
                self.assertTrue(c["split"]["geometrically_split"] and c["q_simple"]["certified"] and not c["split"]["split_over_Q"])
            self.assertTrue((ROOT / c["verification"]["log"]).exists(), c["verification"]["log"])
            self.assertIn("VERIFY_DONE", (ROOT / c["verification"]["log"]).read_text())


if __name__ == "__main__":
    unittest.main(verbosity=1)
