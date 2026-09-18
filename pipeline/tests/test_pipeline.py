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
            for key in list(knowledge._SIMPLE_EXACT) + list(knowledge._SIMPLE_FAMILY) + list(knowledge._SPLIT_FAMILY):
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
        self.assertTrue(recorded <= set(t1), recorded - set(t1))
        # 39 of the record, minus [15] (containment only), plus the trivial group and [2,2,2,12] (certified here)
        self.assertEqual(len(knowledge._SIMPLE_EXACT), 40)
        self.assertEqual(len(knowledge._SIMPLE_FAMILY), 4)   # [13], [17], [19], [15]
        self.assertEqual(knowledge.infinite_record("simple", "15")["grade"], "family")
        self.assertEqual(knowledge.infinite_record("simple", "2.2.2.12")["grade"], "exact")
        self.assertEqual(knowledge.infinite_record("simple", "16")["inherited_from"], "32")
        t2 = []
        for line in (KNOWLEDGE_DIR / "sources" / "table2.txt").read_text().splitlines():
            if line and not line.startswith("#"):
                t2.append(group_key(ast.literal_eval(line.split("|")[0])))
        self.assertEqual(len(t2), 77)
        self.assertTrue(set(knowledge._SPLIT_FAMILY) <= set(t2))
        self.assertEqual(len(knowledge._SPLIT_FAMILY) + len(knowledge._SPLIT_OPEN_NOTES), 20)   # HLP Table 1
        self.assertEqual(list(knowledge._SPLIT_FAMILY_OTHER), ["48"])                             # Howe 2015
        self.assertEqual(knowledge.infinite_record("split", "2")["grade"], "family")             # inherited

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
