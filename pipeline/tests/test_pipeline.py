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

    def test_lmfdb_conductor_override(self):
        ogg = {"conductor": "793600", "conductor_note": "the exponent of 2 in the conductor was computed with Ogg's formula, ..."}
        self.assertEqual(verify.lmfdb_conductor_override(ogg, {"label": "775.a.1"}), "775")
        self.assertEqual(verify.lmfdb_conductor_override(ogg, {"label": ""}), "")
        self.assertEqual(verify.lmfdb_conductor_override(ogg, {"label": "777.a.1"}), "")        # odd parts differ
        self.assertEqual(verify.lmfdb_conductor_override({"conductor": "793600", "conductor_note": ""}, {"label": "775.a.1"}), "")
        self.assertEqual(verify.lmfdb_conductor_override(None, {"label": "775.a.1"}), "")
        self.assertEqual(verify._factor_str(176400), "2^4 * 3^2 * 5^2 * 7^2")

    def test_apply_conductor_without_magma_conductor(self):
        """When Magma's conductor fails, the minimal model it wrote first is kept, the note says why, and an
        isogenous product's conductor is not described as agreeing with a Magma value that does not exist (2026-10-20)."""
        partial = {"ok": False, "minimal_model": {"f": ["-2"], "h": []}, "minimal_discriminant": "-28768731136",
                   "minimal_discriminant_factored": "-1 * 2^22 * 19^3"}
        why = verify.conductor_missing_reason(False, "Fibre blowups over number fields not yet implemented")
        curve = {}
        verify.apply_conductor(curve, None, {"label": "1521.d.2"}, None, why, partial)
        self.assertEqual((curve["minimal_discriminant"], curve["conductor"], curve["conductor_factored"], curve["conductor_source"]),
                         ("-28768731136", "1521", "3^2 * 13^2", "lmfdb label"))
        self.assertIn("failed on this curve: Fibre blowups", curve["conductor_note"])
        prod = {"ainvs": [["0", "1", "1", "-9", "-15"], ["0", "1", "1", "1", "0"]], "conductors": ["19", "19"], "conductor": "361"}
        curve = {}
        verify.apply_conductor(curve, None, {"label": "361.a.1"}, prod, why, partial)
        self.assertEqual((curve["conductor"], curve["conductor_source"]), ("361", "isogenous product of elliptic curves"))
        self.assertNotIn("Magma", curve["conductor_note"])
        curve = {}
        verify.apply_conductor(curve, None, {}, None, verify.conductor_missing_reason(True, ""), None)
        self.assertEqual((curve["minimal_model"], curve["conductor"]), (None, ""))
        self.assertIn("time limit", curve["conductor_note"])
        # a job that logged its marker after an error (Magma carries on) left only the minimal model: not a result
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            cfile = Path(d) / "conductor.json"
            cfile.write_text('{"ok":false,"minimal_model":{"f":["1"],"h":[]},"minimal_discriminant":"1","minimal_discriminant_factored":"1"}')
            self.assertIsNone(verify.conductor_result(cfile, "MINIMAL_MODEL_DONE\nCONDUCTOR_DONE"))
            cfile.write_text('{"ok":true,"minimal_model":{"f":["1"],"h":[]},"conductor":"361","conductor_factored":"19^2","conductor_note":"","cputime":1}')
            self.assertEqual(verify.conductor_result(cfile, "CONDUCTOR_DONE")["conductor"], "361")
            out = Path(d) / "conductor.stdout"
            out.write_text("System Error: User memory limit has been reached\n")
            self.assertEqual(verify.magma_error(out), "User memory limit has been reached")
        # every census curve has a minimal model (the ten missing ones were filled by fill_conductors.py on 2026-10-20)
        for p in sorted(CURVES_DIR.glob("*.json")):
            self.assertTrue(read_json(p)["curve"].get("minimal_model"), p.name)

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
                self.assertIn(rec["grade"], ("exact", "family", "open", "impossible"))
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
        # [13], [15], [16], [17], [24], [25], [26], [2,14], [2,16], [2,2,8], [2,2,14], [2,4,4] (2026-10-14),
        # and [3,6], [2,2,6], [2,2,12] (2026-10-18), [3,9], [4,8], [2,4,8], [28] (2026-10-19)
        self.assertEqual(len(knowledge._SIMPLE_EXACT), 60)
        self.assertEqual(knowledge.infinite_record("simple", "28")["grade"], "exact")   # Platonov–Fedorov 2018
        self.assertEqual(knowledge.infinite_record("simple", "2.4.8")["grade"], "exact")
        self.assertEqual(knowledge.infinite_record("simple", "3.9")["grade"], "exact")
        self.assertEqual(len(knowledge._SIMPLE_FAMILY), 0)
        for key in ("13", "15", "16", "17", "19", "24", "25", "26", "2.14", "2.16", "2.2.8", "2.2.14", "2.2.2.12", "2.4.4"):
            self.assertEqual(knowledge.infinite_record("simple", key)["grade"], "exact", key)
        self.assertEqual(knowledge.infinite_record("simple", "4.8")["grade"], "exact")   # open 2026-10-13 .. 2026-10-19
        # every simple-side exact grade names a residue-class certificate log present here (audit of 2026-10-18)
        import re
        for key, (note, _src) in knowledge._SIMPLE_EXACT.items():
            logs = re.findall(r"family_[A-Za-z0-9_]+_certificate\.log", note)
            self.assertTrue(logs, key)
            for name in logs:
                self.assertTrue((KNOWLEDGE_DIR / "sources" / name).exists(), (key, name))
        t2 = []
        for line in (KNOWLEDGE_DIR / "sources" / "table2.txt").read_text().splitlines():
            if line and not line.startswith("#"):
                t2.append(group_key(ast.literal_eval(line.split("|")[0])))
        self.assertEqual(len(t2), 77)
        self.assertTrue(set(knowledge._SPLIT_FAMILY) <= set(t2))
        self.assertEqual(len(knowledge._SPLIT_FAMILY) + len([k for k in knowledge._SPLIT_OPEN_NOTES if k != "11"]), 20)   # HLP Table 1
        # split [11] (2026-10-20): open, with a remark on why the construction of 11.b gives no family
        rec = knowledge.infinite_record("split", "11")
        self.assertEqual((rec["grade"], rec["sources"]), ("open", ["BNSS2026", "ThisSite"]))
        self.assertIn("split_11_remark.log", rec["note"])
        self.assertTrue((KNOWLEDGE_DIR / "sources" / "split_11_remark.log").exists())
        self.assertEqual(knowledge.infinite_record("split", "7.7")["sources"], ["HLP2000"])
        self.assertEqual(list(knowledge._SPLIT_FAMILY_OTHER), ["48"])                             # Howe 2015
        self.assertEqual(knowledge.infinite_record("split", "48")["grade"], "exact")             # Howe 2015 family, certified 2026-10-19
        self.assertEqual([k for k in knowledge._recorded_keys("split") if knowledge.infinite_record("split", k)["grade"] == "family"], [])
        # split side: [8,8], [2,2,2,16] and (2026-10-18) [14], [2,10], [2,14], [2,2,10], [2,2,2,6], [4,12], [18], [2,18] certified here (family_*_split_certificate.log)
        self.assertEqual(list(knowledge._SPLIT_EXACT)[:10], ["8.8", "2.2.2.16", "14", "2.10", "2.14", "2.2.10", "2.2.2.6", "4.12", "18", "2.18"])
        # the split-side sweep of 2026-10-19: 58 more, among them [3,24], [2,30], [2,48], [27], [36] (graded open before)
        self.assertEqual(len(knowledge._SPLIT_EXACT), 70)   # with [48] (Howe 2015) and [4,16] (2026-10-20)
        for key, (note, _src) in knowledge._SPLIT_EXACT.items():
            logs = re.findall(r"family_[A-Za-z0-9_]+_certificate\.log", note)
            self.assertTrue(logs, key)
            for name in logs:
                self.assertTrue((KNOWLEDGE_DIR / "sources" / name).exists(), (key, name))
        self.assertEqual(knowledge.infinite_record("split", "8.8")["grade"], "exact")
        # simple side: since 2026-10-18 no group is graded "family" (the last three got their own families)
        for key in ["3.6", "2.2.6", "2.2.12"]:
            self.assertEqual(knowledge.infinite_record("simple", key)["grade"], "exact")
        self.assertEqual(knowledge.infinite_record("split", "2.2.2.16")["grade"], "exact")
        names = {"2.4.4": "family_244_certificate.log", "8.8": "family_88_split_certificate.log", "2.2.2.16": "family_22216_split_certificate.log",
                 "14": "family_14_split_certificate.log", "2.10": "family_210_split_certificate.log", "2.14": "family_214_split_certificate.log",
                 "2.2.10": "family_2210_split_certificate.log", "2.2.2.6": "family_2226_split_certificate.log",
                 "4.12": "family_412_split_certificate.log", "18": "family_18_split_certificate.log",
                 "2.18": "family_218_split_certificate.log"}
        for key in names:
            side = "simple" if key == "2.4.4" else "split"
            name = names[key]
            self.assertIn(name, knowledge.infinite_record(side, key)["note"])
            self.assertTrue((KNOWLEDGE_DIR / "sources" / name).exists(), name)

    def test_exclusions(self):
        """The census groups excluded in each class, with the proof page of each exclusion (2026-10-21)."""
        keys = [e["key"] for e in read_json(ROOT / "data" / "groups.json")["groups"]]
        found = {cls: {} for cls in ("qsplit", "gsplit")}
        for key in keys:
            for cls, ex in knowledge.exclusions_record(key).items():
                found[cls][key] = ex["sources"][0]
                self.assertTrue((ROOT / ex["proof_url"]).exists(), ex["proof_url"])
                for s in ex["sources"]:
                    self.assertIn(s, knowledge.SOURCES, (key, cls, s))
        mazur = ["11", "13", "17", "19", "22", "23", "26", "29", "31", "33", "34", "39", "2.22", "2.26"]
        self.assertEqual(found["qsplit"], {**{k: "QSplitPrimeObstruction" for k in mazur},
                                           "2.2.14": "QSplitTwoRankObstruction", "2.2.2.10": "QSplitTwoRankObstruction"})
        self.assertEqual(found["gsplit"], {
            "22": "Geometric22Obstruction", "2.22": "Geometric22Obstruction",
            "23": "Geometric2329Obstruction", "29": "Geometric2329Obstruction",
            "26": "Geometric2639Obstruction", "2.26": "Geometric2639Obstruction", "39": "Geometric2639Obstruction",
            "17": "Geometric17Obstruction",
            "31": "Geometric31Obstruction", "33": "Geometric3334Obstruction", "34": "Geometric3334Obstruction",
            "13": "QSimple13Obstruction",
            "30": "QSimple30Obstruction", "60": "QSimple30Obstruction", "2.30": "QSimple30Obstruction",
            "35": "QSimple35Obstruction", "70": "QSimple35Obstruction",
            "63": "QSimple63Obstruction",
            "2.2.14": "QSimple2214Obstruction",
            "2.28": "QSimple228Obstruction",
            "2.2.10": "QSimple2210Obstruction", "2.2.20": "QSimple2210Obstruction",
            "2.2.2.10": "QSimple2210Obstruction",
            "5.5": "QSimpleRankTwoObstruction", "5.10": "QSimpleRankTwoObstruction", "7.7": "QSimpleRankTwoObstruction"})
        impossible = sorted(k for k in keys if knowledge.infinite_record("split", k)["grade"] == "impossible")
        self.assertEqual(impossible, sorted(["13", "17", "22", "23", "26", "29", "31", "33", "34", "39", "2.22", "2.26", "2.2.14",
                                             "2.2.2.10"]))
        # [17] has its own proof; [34] keeps the earlier proof for orders 33 and 34
        self.assertEqual(knowledge.infinite_record("split", "17")["proof_url"], "proofs/geometric-17-obstruction.html")
        rec = knowledge.infinite_record("split", "13")
        self.assertIn("Mazur1977", rec["sources"])
        self.assertTrue(rec["note"].startswith(knowledge.exclusions_record("13")["qsplit"]["reason"]))

    def test_qsimple_seven_two_primary_bound(self):
        """The principal-polarization bound changes only the Q-simple split class."""
        for key in ("2.28", "56", "112", "2.56"):
            with self.subTest(key=key):
                rec = knowledge.exclusions_record(key)
                self.assertEqual(set(rec), {"gsplit"})
                self.assertEqual(rec["gsplit"]["sources"][0], "QSimple228Obstruction")
                self.assertTrue((ROOT / rec["gsplit"]["proof_url"]).exists())
                self.assertIn("principal", rec["gsplit"]["reason"])
                self.assertIn("QSimple228Obstruction", knowledge.QSIMPLE_ONLY)
                self.assertNotEqual(knowledge.infinite_record("split", key)["grade"], "impossible")
        for key in ("7", "14", "28", "2.14"):
            self.assertNotIn("gsplit", knowledge.exclusions_record(key), key)
        # Keep the older proof for groups that already had an exclusion.
        self.assertEqual(knowledge.exclusions_record("2.2.14")["gsplit"]["sources"][0], "QSimple2214Obstruction")
        groups = read_json(ROOT / "data" / "groups.json")["groups"]
        group = next(g for g in groups if g["key"] == "2.28")
        self.assertTrue(group["known"]["simple"])
        self.assertFalse(group["known"]["qsplit"])
        self.assertFalse(group["known"]["gsplit"])
        self.assertEqual(group["infinite"]["split"]["grade"], "open")
        self.assertEqual(group["exclusions"], knowledge.exclusions_record("2.28"))
        log = (KNOWLEDGE_DIR / "sources" / "qsimple_228_checks.log").read_text()
        self.assertTrue(log.rstrip().endswith("VERIFIED_228_MODULAR_CERTIFICATE_DONE"))
        for bad in ("Runtime error", "Assertion failed", "User error", "Internal error"):
            self.assertNotIn(bad, log)
        for suffix in ("md", "tex", "pdf"):
            self.assertTrue((KNOWLEDGE_DIR / "sources" / f"qsimple_228_obstruction.{suffix}").exists())

    def test_qsimple_five_two_rank(self):
        """Rational 5-torsion forces 2-rank at most 2 on Q-simple, geometrically split Jacobians (2026-11-11)."""
        for key in ("2.2.10", "2.2.20", "2.2.2.10", "2.2.40", "2.10.10"):
            with self.subTest(key=key):
                rec = knowledge.exclusions_record(key)
                self.assertEqual(rec["gsplit"]["sources"][0], "QSimple2210Obstruction")
                self.assertTrue((ROOT / rec["gsplit"]["proof_url"]).exists())
                self.assertIn("principal", rec["gsplit"]["reason"])
        self.assertIn("QSimple2210Obstruction", knowledge.QSIMPLE_ONLY)
        # 2-rank at most 2, or no 5-torsion: no exclusion from this proof
        for key in ("5", "10", "20", "40", "2.10", "2.20", "2.2.2", "2.2.8"):
            with self.subTest(key=key):
                rec = knowledge.exclusions_record(key)
                self.assertNotEqual(rec.get("gsplit", {}).get("sources", [None])[0], "QSimple2210Obstruction")
        # groups that already had a gsplit exclusion keep their proof
        self.assertEqual(knowledge.exclusions_record("2.2.30")["gsplit"]["sources"][0], "QSimple30Obstruction")
        # combined split grades: [2,2,2,10] is impossible (qsplit by the 2-rank obstruction); the others keep theirs
        rec = knowledge.infinite_record("split", "2.2.2.10")
        self.assertEqual(rec["grade"], "impossible")
        self.assertEqual(rec["proof_url"], "proofs/qsimple-2210-obstruction.html")
        self.assertTrue(rec["note"].startswith(knowledge.exclusions_record("2.2.2.10")["qsplit"]["reason"]))
        self.assertEqual(knowledge.infinite_record("split", "2.2.10")["grade"], "exact")
        self.assertEqual(knowledge.infinite_record("split", "2.2.20")["grade"], "open")
        groups = {g["key"]: g for g in read_json(ROOT / "data" / "groups.json")["groups"]}
        for key in ("2.2.10", "2.2.20", "2.2.2.10"):
            self.assertFalse(groups[key]["known"]["gsplit"], key)
            self.assertEqual(groups[key]["exclusions"], knowledge.exclusions_record(key))
        self.assertTrue(groups["2.2.10"]["known"]["qsplit"])
        self.assertEqual(groups["2.2.2.10"]["infinite"]["split"]["grade"], "impossible")
        log = (KNOWLEDGE_DIR / "sources" / "qsimple_2210_checks.log").read_text()
        self.assertTrue(log.rstrip().endswith("ALL_CHECKS_PASSED"))
        for bad in ("Runtime error", "Assertion failed", "User error", "Internal error"):
            self.assertNotIn(bad, log)
        vlog = (KNOWLEDGE_DIR / "sources" / "qsimple_2210_validation.log").read_text()
        self.assertTrue(vlog.rstrip().endswith("VAL_DONE"))
        for name in ("qsimple_2210_obstruction.tex", "qsimple_2210_obstruction.pdf", "qsimple_2210_checks.m",
                     "qsimple_2210_validation.m"):
            self.assertTrue((KNOWLEDGE_DIR / "sources" / name).exists(), name)

    def test_qsimple_five_two_rank_future_entries_rejected(self):
        """The build rejects future Q-simple split examples with 5-torsion and 2-rank at least three."""
        import build
        import copy
        import json
        import tempfile
        from unittest.mock import patch
        template = next(read_json(p) for p in sorted(CURVES_DIR.glob("*.json"))
                        if read_json(p)["class"] == "gsplit")
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "future.json"
            for key in ("2.2.10", "2.2.20", "2.2.40"):
                c = copy.deepcopy(template)
                c.update(id=f"{key}.future", group_key=key, group=[int(n) for n in key.split(".")])
                c["torsion"]["invariants"] = c["group"]
                path.write_text(json.dumps(c))
                with self.subTest(key=key), patch.object(build, "CURVES_DIR", Path(d)):
                    with self.assertRaisesRegex(ValueError, "recorded gsplit example contradicts"):
                        build.build()

    def test_qsimple_seven_two_primary_future_entries_rejected(self):
        """The build rejects future contradictory Q-simple split examples."""
        import build
        import copy
        import json
        import tempfile
        from unittest.mock import patch
        template = next(read_json(p) for p in sorted(CURVES_DIR.glob("*.json"))
                        if read_json(p)["class"] == "gsplit")
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "future.json"
            for key in ("2.28", "56", "112", "2.56"):
                c = copy.deepcopy(template)
                c.update(id=f"{key}.future", group_key=key, group=[int(n) for n in key.split(".")])
                c["torsion"]["invariants"] = c["group"]
                path.write_text(json.dumps(c))
                with self.subTest(key=key), patch.object(build, "CURVES_DIR", Path(d)):
                    with self.assertRaisesRegex(ValueError, "recorded gsplit example contradicts"):
                        build.build()

    def test_order256_exclusions_beyond_census(self):
        """The new exact-group exclusions apply without census rows and only to qsplit."""
        excluded = {"4.64", "8.32", "16.16", "2.2.64", "2.4.32", "2.8.16"}
        self.assertEqual(knowledge.QSPLIT_ORDER256_EXCLUSIONS, excluded)
        for key in excluded:
            with self.subTest(key=key):
                rec = knowledge.exclusions_record(key)
                self.assertEqual(set(rec), {"qsplit"})
                self.assertEqual(rec["qsplit"]["sources"][0], "QSplitOrder256Obstruction")
                self.assertTrue((ROOT / rec["qsplit"]["proof_url"]).exists())
                self.assertTrue(set(rec["qsplit"]["sources"]) <= set(knowledge.SOURCES))
                self.assertNotEqual(knowledge.infinite_record("split", key)["grade"], "impossible")
        # The theorem does not settle these rank-four groups, or proper subgroups.
        for key in ("2.2.2.32", "2.2.4.16", "2.2.8.8", "4.32", "8.16", "2.4.16"):
            self.assertNotIn("qsplit", knowledge.exclusions_record(key), key)
        log = (KNOWLEDGE_DIR / "sources" / "qsplit_order256_matrix_checks.log").read_text()
        self.assertTrue(log.rstrip().endswith("ORDER256_TRACE_CERTIFICATE_DONE"))
        for suffix in ("md", "tex", "pdf"):
            self.assertTrue((KNOWLEDGE_DIR / "sources" / f"qsplit_order256_obstruction.{suffix}").exists())

    def test_order256_future_qsplit_entries_rejected(self):
        """The build catches each future contradiction even when its group is currently absent."""
        import build
        import copy
        import json
        import tempfile
        from unittest.mock import patch
        template = next(read_json(p) for p in sorted(CURVES_DIR.glob("*.json"))
                        if read_json(p)["class"] == "qsplit")
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "future.json"
            for key in ("4.64", "8.32", "16.16", "2.2.64", "2.4.32", "2.8.16"):
                c = copy.deepcopy(template)
                c.update(id=f"{key}.future", group_key=key, group=[int(n) for n in key.split(".")])
                c["torsion"]["invariants"] = c["group"]
                path.write_text(json.dumps(c))
                with self.subTest(key=key), patch.object(build, "CURVES_DIR", Path(d)):
                    with self.assertRaisesRegex(ValueError, "recorded qsplit example contradicts"):
                        build.build()

    def test_certificate_logs_complete(self):
        # every certificate log ends in its CERTIFIED marker and records no Magma error (a failed assert stops
        # `magma -b script.m` before the marker, but the exit status is 0, so the log is what must be checked)
        logs = sorted((KNOWLEDGE_DIR / "sources").glob("*certificate*.log"))
        self.assertGreater(len(logs), 100)
        for p in logs:
            text = p.read_text(errors="replace")
            self.assertIn("CERTIFIED", text, p.name)
            for bad in ("Runtime error", "Assertion failed", "User error", "Internal error"):
                self.assertNotIn(bad, text, p.name)

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
        self.assertEqual(len(flagged), 55)   # 19 corrected on 2026-10-14; 33 + 4 (Leprévost's notes) - 1 (3.9.a) on 2026-10-20
        self.assertIn("35.a", flagged)
        self.assertNotIn("27.b", flagged)   # its sources were replaced as well (SOURCES_OVERRIDE)
        self.assertNotIn("2.28.a", flagged)   # likewise: Platonov–Petrunin 2015's f28,2 (2026-10-20)
        self.assertNotIn("3.9.a", flagged)    # likewise: the first geometrically simple [3,9] in the census

    def test_seed_database_credits_match_production_check(self):
        """A seed curve credited to the production LMFDB database is isomorphic over Q to a production curve,
        and one credited to the extended database is not (seed_production_check.log, 2026-10-20)."""
        import re
        hits = {}
        for line in (KNOWLEDGE_DIR / "sources" / "seed_production_check.log").read_text().splitlines():
            m = re.match(r"ISO (\S+) (\S+) cands=\d+ \[(.*)\]$", line)
            if m:
                hits[m.group(1)] = [h.strip().strip('"') for h in m.group(3).split(",") if h.strip()]
        self.assertEqual(len(hits), 150)
        bssvy, bs = knowledge.credit_text(["BSSVY"]), knowledge.credit_text(["BookerSutherland"])
        n = 0
        for p in sorted(CURVES_DIR.glob("*.json")):
            c = read_json(p)
            if (c.get("source") or {}).get("kind") != "import":
                continue
            if c["discovery"]["by"] == bssvy:
                self.assertTrue(hits[c["id"]], c["id"]); n += 1
            if c["discovery"]["by"] == bs:
                self.assertEqual(hits[c["id"]], [], c["id"]); n += 1
        self.assertEqual(n, 107)   # 38 production, 69 extended
        self.assertTrue((KNOWLEDGE_DIR / "sources" / "seed_credit_check.log").read_text().rstrip().endswith("SEED_CREDIT_CHECK_DONE"))
        c = read_json(CURVES_DIR / "3.9.c.json")   # Leprévost 1995's C_3 (p. 303)
        self.assertEqual((c["sources"], c["discovery"]["by"], c["discovery"]["year"]), (["Leprevost1995"], "Leprévost 1995", 1995))

    def test_leprevost_3_9_family_is_geometrically_split(self):
        """Every member of Leprévost 1995's [3,9] family covers an elliptic curve (leprevost_3_9_split.log), so the
        geometrically simple [3,9] is dated by the database curve 3.9.a, not by Leprévost 1995 (2026-10-20)."""
        log = (KNOWLEDGE_DIR / "sources" / "leprevost_3_9_split.log").read_text()
        self.assertTrue(log.rstrip().endswith("is geometrically split"))
        self.assertIn("LEPREVOST_3_9_SPLIT_DONE", log)
        self.assertEqual(read_json(CURVES_DIR / "3.9.a.json")["sources"], ["BookerSutherland"])
        g = {e["key"]: e for e in read_json(ROOT / "data" / "groups.json")["groups"]}["3.9"]
        self.assertEqual(g["first_source"]["simple"], "BookerSutherland")
        self.assertEqual(g["first_source"]["gsplit"], "Leprevost1995")
        self.assertTrue((KNOWLEDGE_DIR / "sources" / "leprevost_notes_check.log").read_text().rstrip().endswith("LEPREVOST_NOTES_CHECK_DONE"))
        for i in ("15.a", "19.a", "22.a", "24.a", "26.a"):   # not curves of Leprévost's C. R. notes (leprevost_notes_check.log)
            self.assertNotIn("Leprévost", read_json(CURVES_DIR / f"{i}.json")["discovery"]["by"], i)


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

    def test_no_curve_left_undecided(self):
        """Every census curve has a certified class; the last undecided one, 2.2.2.2.b, was certified simple over Q
        by the elliptic product test on 2026-11-14."""
        for p in sorted(CURVES_DIR.glob("*.json")):
            c = read_json(p)
            self.assertIn(c["class"], CLASSES, c["id"])
        groups = read_json(ROOT / "data" / "groups.json")
        self.assertEqual(groups["n_verified"], 0)
        self.assertEqual([g["key"] for g in groups["groups"] if g["classes"]["undecided"]], [])

    def test_2222b_simple_over_Q(self):
        """2.2.2.2.b (784.c.1) is (2,2)-isogenous over Q to the Weil restriction of an elliptic curve over Q(sqrt(7)), and
        every chi_p is reducible; chi_5 = (T^2 - 5)^2 is not (T^2 - aT + 5)(T^2 - bT + 5) with integers a, b, so the
        Jacobian is simple over Q.  Re-verified with `verify.py --reverify 2.2.2.2.b` on 2026-11-14."""
        c = read_json(CURVES_DIR / "2.2.2.2.b.json")
        self.assertEqual((c["class"], c["class_certified"], c["status"]), ("gsplit", True, "certified"))
        self.assertTrue(c["split"]["geometrically_split"])
        self.assertFalse(c["split"]["split_over_Q"])
        self.assertIn("Weil restriction of an elliptic curve E over Q(sqrt(7))", c["split"]["certificate"]["detail"])
        self.assertEqual(c["simplicity"]["reducible_chi"], c["simplicity"]["admissible_primes"])
        qs = c["q_simple"]
        self.assertEqual((qs["certified"], qs["prime"], qs["chi"], qs["method"]),
                         (True, 5, "T^4 - 10*T^2 + 25", "chi_not_elliptic_product"))
        self.assertEqual(qs["detail"], "chi_5 = T^4 - 10*T^2 + 25 is not of the form (T^2 - a*T + 5)*(T^2 - b*T + 5) "
                                       "with integers a, b (s1 = 0, s2 = -10: s1^2 - 4*(s2 - 2*p) = 80 is not a square)")
        log = (ROOT / c["verification"]["log"]).read_text()
        self.assertIn("Q-simple: " + qs["detail"], log)
        self.assertIn("class: gsplit", log)
        # the curated fields and the acceptance record survive the re-verification
        self.assertEqual((c["sources"], c["discovery"]["by"], c["lmfdb"]["label"], c["accepted_because"], c["dates"]["submitted"]),
                         (["BookerSutherland"], "Booker–Sutherland (extended database)", "784.c.1", "new for this class", "2026-08-30"))
        g = {e["key"]: e for e in read_json(ROOT / "data" / "groups.json")["groups"]}["2.2.2.2"]
        row = next(s for s in g["classes"]["gsplit"] if s["id"] == "2.2.2.2.b")
        self.assertIn("simple over ℚ: χ_5 = T^4 - 10*T^2 + 25 is not of the form", row["certificate"])

    def test_elliptic_product_certificates(self):
        """Every certificate of simplicity over Q by the elliptic product test is rechecked: p is a good prime of the
        even model, chi_p = T^4 - s1 T^3 + s2 T^2 - p s1 T + p^2, and s1^2 - 4(s2 - 2p) is not a perfect square, so no
        integers a, b give chi_p = (T^2 - aT + p)(T^2 - bT + p), as they would at every good prime if J ~ E x E' over Q."""
        import math
        n = 0
        for p in sorted(CURVES_DIR.glob("*.json")):
            c = read_json(p)
            qs = c["q_simple"]
            if "method" not in qs:
                continue   # an irreducible chi_p (the certificates verified before 2026-11-14 have no "method")
            n += 1
            self.assertEqual((qs["method"], qs["certified"]), ("chi_not_elliptic_product", True), c["id"])
            q = qs["prime"]
            disc_g, g = int(c["curve"]["discriminant"]), [int(a) for a in c["curve"]["g"]]
            self.assertTrue(q % 2 and disc_g % q and g[-1] % q and len(g) in (6, 7), c["id"])
            cs = [int(a) for a in coeffs_to_strings(parse_poly(qs["chi"].replace("T", "x")))]
            self.assertEqual((len(cs), cs[4], cs[1], cs[0]), (5, 1, q * cs[3], q * q), c["id"])
            s1, s2 = -cs[3], cs[2]
            disc = s1 * s1 - 4 * (s2 - 2 * q)
            self.assertFalse(disc >= 0 and math.isqrt(disc) ** 2 == disc, c["id"])
            self.assertFalse(any(a * (s1 - a) == s2 - 2 * q for a in range(-4 * q, 4 * q + 1)), c["id"])
            self.assertIn(f"chi_{q} = {qs['chi']} is not of the form", qs["detail"])
            self.assertIn(f"(s1 = {s1}, s2 = {s2}: s1^2 - 4*(s2 - 2*p) = {disc} is not a square)", qs["detail"])
            self.assertNotEqual(c["class"], "qsplit", c["id"])
        self.assertGreaterEqual(n, 1)


class ReverifyTests(unittest.TestCase):
    def test_reverify_finds_every_submission(self):
        """`verify.py --reverify ID` finds the accepted submission of every census curve by its equation, keeps its Magma
        inputs and takes the curated fields from the certificate (61 certificates had their credits or notes corrected
        after their submission was processed)."""
        for p in sorted(CURVES_DIR.glob("*.json")):
            c = read_json(p)
            path, sub, cert = verify.reverify_submission(c["id"])
            self.assertEqual(cert["id"], c["id"])
            self.assertTrue(path.name.endswith(".accepted.json"), path.name)
            v = verify.validate(sub)
            self.assertEqual((v["f"], v["h"], v.get("group"), v.get("cover")),
                             (c["submitted"]["f"], c["submitted"]["h"], c["submitted"]["claimed_group"], c["submitted"]["cover"]), c["id"])
            self.assertEqual((v.get("sources", []), v.get("discoverer", ""), v.get("year"), v["historical"], v.get("notes", ""),
                              v.get("reference", ""), v.get("date", ""), v.get("lmfdb_label", "")),
                             (c["sources"], c["discovery"]["by"], c["discovery"]["year"], c["discovery"]["historical"], c["notes"],
                              c["reference"], c["dates"]["submitted"], c["submitted"]["lmfdb_label_given"]), c["id"])
        self.assertEqual(verify.reverify_submission("2.2.2.2.b")[0].name, "paper-t2-072.accepted.json")
        with self.assertRaises(SystemExit):
            verify.reverify_submission("2.2.2.2.zz")


if __name__ == "__main__":
    unittest.main(verbosity=1)
