#!/usr/bin/env python3
"""Verify submitted genus-2 curves with Magma and write certificates.

Usage
-----
  python3 pipeline/verify.py                 # process every file in submissions/inbox/
  python3 pipeline/verify.py FILE [FILE...]  # process the given submission files
  options: --timeout SECONDS (default 3600)  --mem-gb N (default 16)  --no-github
           --conductor-timeout SECONDS (default 600)  --no-conductor  --no-lmfdb
           --keep (do not move processed submissions)  --dry-run  --require-claim

For each submission the script
  1. validates the JSON (schema, whitelisted polynomial strings, invariant factors);
  2. writes a Magma job that loads pipeline/magma/verify_lib.m and runs it under `timeout`:
     torsion subgroup, geometric simplicity / Q-simplicity / splitness certificates, invariants,
     Q-isomorphism with the curves already in the census;
  3. runs a second Magma job (pipeline/magma/conductor_lib.m) for the minimal model and the
     conductor, under its own time limit (it needs the discriminant factored);
  4. records the LMFDB label (given with the submission, or by one lookup by equation on
     www.lmfdb.org for new submissions; a production label is a permalink, an extended-database
     label is linked by equation);
  5. writes data/curves/<id>.json (status "certified" when the torsion and the class are both
     certified, "verified" when only the torsion is) or data/rejected/<name>.json, keeps the Magma
     logs under data/logs/, and
  6. comments on / closes the GitHub issue the submission came from (unless --no-github).

Every verified curve is accepted (a duplicate over Q of a census curve is rejected); the
certificate records whether its (group, class) was new to the census at verification time.
Magma exits with status 0 even after an error, so success is judged from the JSON it wrote
(`ok: true`) and the VERIFY_DONE marker in the log, never from the exit code.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import knowledge  # noqa: E402
from common import (CLASSES, CURVES_DIR, LOGS, LOGS_DIR, MAX_STRING, RE_ELEMENT_W, RE_POLY_X_UPPER,  # noqa: E402
                    REJECTED_DIR, SCHEMA_CERTIFICATE, SCHEMA_REJECTED, SITE_URL, SUBMISSIONS_INBOX,
                    SUBMISSIONS_PROCESSED, UNDECIDED, WORK, coeff_lists_from_submission, github_request,
                    github_token, group_bracket, group_key, group_order, is_invariant_factors, lmfdb_jump_url,
                    lmfdb_lookup, log, read_json, repo_path, write_json)

VERIFY_LIB = Path(__file__).resolve().parent / "magma" / "verify_lib.m"
CONDUCTOR_LIB = Path(__file__).resolve().parent / "magma" / "conductor_lib.m"
MAGMA = shutil.which("magma") or "/usr/local/bin/magma"
HOST = "Mordell (Dept. of Mathematics, Univ. of Zagreb)"


class Reject(Exception):
    pass


# --------------------------------------------------------------------------- validation

def _coeff_strings(name, lst):
    if not isinstance(lst, list) or not lst:
        raise Reject(f"{name} must be a non-empty list of coefficients")
    out = []
    for c in lst:
        c = str(c).replace(" ", "")
        if not c or len(c) > MAX_STRING or not RE_ELEMENT_W.match(c):
            raise Reject(f"{name} contains characters outside the allowed set")
        out.append(c)
    return out


def validate(sub: dict) -> dict:
    """Return a normalised copy of the submission or raise Reject."""
    try:
        f, h = coeff_lists_from_submission(sub)
    except Exception as e:
        raise Reject(f"could not parse the curve: {e}")
    if len(f) > 7 or len(h) > 4:
        raise Reject("f must have degree at most 6 and h degree at most 3")
    if all(c == "0" for c in f) and all(c == "0" for c in h):
        raise Reject("the curve is empty")
    out = {"f": f, "h": h}
    grp = sub.get("group")
    if grp not in (None, "", []):
        try:
            inv = [int(n) for n in (grp if isinstance(grp, list) else str(grp).strip("[] ").split(","))]
        except Exception:
            raise Reject("the torsion group must be a list of invariant factors, e.g. [2, 2, 12] (or [] for the trivial group)")
        if not is_invariant_factors(inv):
            raise Reject(f"{inv} is not a list of invariant factors n1 | n2 | ... | n4 with n1 >= 2")
        out["group"] = inv
    cls = str(sub.get("class", "") or "").strip().lower()
    if cls in ("geometrically simple", "simple"):
        cls = "simple"
    elif cls in ("split", "geometrically split", "split over q", "qsplit", "gsplit"):
        cls = "split"
    elif cls in ("", "unknown", "?"):
        cls = ""
    else:
        raise Reject(f"unknown class {cls!r}: use 'simple', 'split' or leave it empty")
    out["class"] = cls
    cov = sub.get("cover")
    if cov:
        if not isinstance(cov, dict):
            raise Reject("cover must be an object {field, cubic | cremona, p, q, hh, note}")
        fld = str(cov.get("field", "") or "").replace(" ", "")
        if fld and (len(fld) > MAX_STRING or not RE_POLY_X_UPPER.match(fld)):
            raise Reject("cover.field must be a polynomial in X with rational coefficients")
        cremona = str(cov.get("cremona", "") or "").strip()
        if cremona and not cremona.replace(".", "").isalnum():
            raise Reject("cover.cremona must be a Cremona label such as 19a1")
        cubic = _coeff_strings("cover.cubic", cov["cubic"]) if cov.get("cubic") else []
        if cubic and len(cubic) != 4:
            raise Reject("cover.cubic must be [e3, e2, e1, e0] for v^2 = e3 u^3 + e2 u^2 + e1 u + e0")
        if not cubic and not cremona:
            raise Reject("cover needs either cubic = [e3, e2, e1, e0] or a Cremona label")
        if cremona and fld:
            raise Reject("a Cremona label describes a curve over Q: leave cover.field empty")
        out["cover"] = {"field": fld, "cremona": cremona, "cubic": cubic,
                        "p": _coeff_strings("cover.p", cov.get("p")), "q": _coeff_strings("cover.q", cov.get("q")),
                        "hh": _coeff_strings("cover.hh", cov.get("hh")), "note": str(cov.get("note", "") or "")[:500]}
    for k in ("submitter", "github", "reference", "notes", "date", "source", "affiliation", "discoverer", "year",
              "lmfdb_label", "sources", "rm", "route"):
        if k in sub:
            out[k] = sub[k]
    return out


# --------------------------------------------------------------------------- Magma

def magma_string(s: str) -> str:
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def magma_list(strs) -> str:
    return "[" + ", ".join(magma_string(s) for s in strs) + "]"


def write_job(v: dict, jobdir: Path, mem_gb: int, existing) -> Path:
    lines = ["SetColumns(0);",
             f"fcoeffs := {magma_list(v['f'])};",
             f"hcoeffs := {magma_list(v['h'])};",
             f"LogFile := {magma_string(str(jobdir / 'magma.log'))};",
             f"OutFile := {magma_string(str(jobdir / 'result.json'))};",
             f"MemGB := {mem_gb};"]
    if v.get("group") is not None:
        lines.append("ClaimedGroup := [" + ", ".join(str(n) for n in v["group"]) + "];")
    if v.get("route") == "Q2":
        lines.append("QuadDiscs := [2];")
    if v.get("cover"):
        c = v["cover"]
        lines.append("Cover := < " + ", ".join([magma_string(c["field"]), magma_list(c["cubic"]), magma_list(c["p"]),
                                                magma_list(c["q"]), magma_list(c["hh"]), magma_string(c["note"]),
                                                magma_string(c["cremona"])]) + " >;")
    if existing:
        lines.append("ExistingModels := [* " + ", ".join(
            "< " + magma_string(e["id"]) + ", " + magma_list(e["f"]) + ", " + magma_list(e["h"]) + " >" for e in existing) + " *];")
    lines.append(f'load "{VERIFY_LIB}";')
    lines.append("quit;")
    job = jobdir / "job.m"
    job.write_text("\n".join(lines) + "\n")
    return job


def write_conductor_job(v: dict, jobdir: Path, mem_gb: int) -> Path:
    lines = ["SetColumns(0);",
             f"fcoeffs := {magma_list(v['f'])};",
             f"hcoeffs := {magma_list(v['h'])};",
             f"LogFile := {magma_string(str(jobdir / 'conductor.log'))};",
             f"OutFile := {magma_string(str(jobdir / 'conductor.json'))};",
             f"MemGB := {mem_gb};",
             f'load "{CONDUCTOR_LIB}";',
             "quit;"]
    job = jobdir / "conductor.m"
    job.write_text("\n".join(lines) + "\n")
    return job


def run_magma(job: Path, timeout: int) -> tuple[int, bool]:
    """Run the job; return (exit code, timed_out).  stdin is /dev/null so an error never waits."""
    try:
        with open(job.with_suffix(".stdout"), "w") as out:
            proc = subprocess.run([MAGMA, "-b", str(job)], stdin=subprocess.DEVNULL, stdout=out,
                                  stderr=subprocess.STDOUT, timeout=timeout, cwd=job.parent)
        return proc.returncode, False
    except subprocess.TimeoutExpired:
        return -1, True


def denull(obj):
    """verify_lib.m writes {"null": true} for null."""
    if isinstance(obj, dict):
        if obj == {"null": True}:
            return None
        return {k: denull(x) for k, x in obj.items()}
    if isinstance(obj, list):
        return [denull(x) for x in obj]
    return obj


# --------------------------------------------------------------------------- post-processing

def existing_curves():
    return [read_json(p) for p in sorted(CURVES_DIR.glob("*.json"))]


def next_id(key: str, curves) -> str:
    used = {c["id"].split(".")[-1] for c in curves if c["group_key"] == key}
    k = 0
    while True:
        label, j = "", k
        while True:
            label = chr(ord("a") + j % 26) + label
            j = j // 26 - 1
            if j < 0:
                break
        if label not in used:
            return f"{key}.{label}"
        k += 1


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_certificate(v: dict, res: dict, cond: dict | None, curves, pid: str, now: str, lmfdb: dict) -> dict:
    inv = res["torsion"]["invariants"]
    key = group_key(inv)
    cls = res["class"]
    certified = cls in CLASSES
    known_classes = {c["class"] for c in curves if c["group_key"] == key and c["class"] in CLASSES}
    known_any = any(c["group_key"] == key for c in curves)
    curve = dict(res["curve"])
    if cond:
        curve.update({"minimal_model": cond["minimal_model"], "minimal_discriminant": cond["minimal_discriminant"],
                      "minimal_discriminant_factored": cond["minimal_discriminant_factored"],
                      "conductor": cond["conductor"], "conductor_factored": cond["conductor_factored"],
                      "conductor_note": cond["conductor_note"]})
    else:
        curve.update({"minimal_model": None, "minimal_discriminant": "", "minimal_discriminant_factored": "",
                      "conductor": "", "conductor_factored": "", "conductor_note": "not computed by Magma within the time limit"})
        if lmfdb.get("label"):
            # the first component of an LMFDB label is the conductor of the Jacobian
            curve["conductor"] = lmfdb["label"].split(".")[0]
            curve["conductor_note"] = "taken from the LMFDB label " + lmfdb["label"] + " (Magma did not finish the conductor computation within the time limit)"
    curve["conductor_source"] = "magma" if cond else ("lmfdb label" if curve["conductor"] else "")
    sources = v.get("sources") or []
    return {
        "schema": SCHEMA_CERTIFICATE,
        "id": pid,
        "group": inv, "group_key": key, "group_bracket": group_bracket(inv), "order": group_order(inv),
        "class": cls, "class_certified": certified,
        "class_text": CLASSES.get(cls) or UNDECIDED.get(cls, cls),
        "status": "certified" if certified else "verified",
        "new_group": not known_any,
        "new_for_class": certified and cls not in known_classes,
        "curve": curve,
        "lmfdb": lmfdb,
        "torsion": res["torsion"],
        "simplicity": res["simplicity"],
        "q_simple": res["q_simple"],
        "split": res["split"],
        "geometrically_isomorphic_to": res.get("same_g2", []),
        "discovery": {"by": v.get("discoverer", ""), "year": v.get("year", "")},
        "reference": v.get("reference", ""),
        "sources": sources,
        "rm": bool(v.get("rm", False)),
        "notes": v.get("notes", ""),
        "submitter": {"name": v.get("submitter", ""), "github": v.get("github", ""), "affiliation": v.get("affiliation", "")},
        "submitted": {"f": v["f"], "h": v["h"], "claimed_group": v.get("group"), "claimed_class": v.get("class", ""),
                      "cover": v.get("cover"), "lmfdb_label_given": v.get("lmfdb_label", "")},
        "source": v.get("source", {}),
        "dates": {"submitted": v.get("date", ""), "verified": now},
        "verification": {
            "host": HOST, "magma_version": res["magma_version"], "cputime_seconds": res["cputime"],
            "verify_lib_sha256": sha256_file(VERIFY_LIB),
            "log": f"data/logs/{pid}.log",
            "conductor_log": f"data/logs/{pid}.conductor.log" if cond else "",
            "conductor_cputime_seconds": cond["cputime"] if cond else None,
        },
    }


def lmfdb_record(v: dict, res: dict, do_lookup: bool) -> dict:
    """Link policy (as in the Genus2Torsion repository): a production label cond.class.disc.num is a
    permalink on www.lmfdb.org; an extended-database label cond.class.num (alpha.lmfdb.org, a dated
    snapshot) is not permanent, so the curve is linked by its equation.  Without a label, one lookup
    by equation is attempted on www.lmfdb.org (the site rate-limits automated requests, so this may
    fail); otherwise the curve is linked to the LMFDB search by equation."""
    f, h = res["curve"]["f"], res["curve"]["h"]
    given = str(v.get("lmfdb_label", "") or "").strip()
    if given and given != "-":
        parts = given.split(".")
        if len(parts) == 4:
            return {"label": given, "kind": "production", "url": f"https://www.lmfdb.org/Genus2Curve/Q/{given.replace('.', '/')}",
                    "snapshot_label": ""}
        if len(parts) == 3:
            return {"label": given, "kind": "alpha-snapshot", "url": lmfdb_jump_url(f, h, "alpha"), "snapshot_label": given}
    label = lmfdb_lookup(f, h) if do_lookup else None
    if label:
        return {"label": label, "kind": "production", "url": f"https://www.lmfdb.org/Genus2Curve/Q/{label.replace('.', '/')}",
                "snapshot_label": ""}
    return {"label": "", "kind": "jump", "url": lmfdb_jump_url(f, h, "www"), "snapshot_label": ""}


# --------------------------------------------------------------------------- GitHub feedback

def github_feedback(v: dict, outcome: str, body: str, dry: bool):
    src = v.get("source") or {}
    if src.get("kind") != "issue":
        return
    tok = github_token()
    if not tok:
        log("  no GitHub token: skipping issue update")
        return
    num = src["number"]
    if dry:
        log(f"  [dry-run] would comment on issue #{num} and label it {outcome}")
        return
    github_request("POST", repo_path(f"/issues/{num}/comments"), {"body": body}, token=tok)
    github_request("POST", repo_path(f"/issues/{num}/labels"), {"labels": [outcome]}, token=tok)
    github_request("PATCH", repo_path(f"/issues/{num}"), {"state": "closed"}, token=tok)
    log(f"  issue #{num}: commented, labelled '{outcome}', closed")


def site_url(pid: str) -> str:
    return f"{SITE_URL}/curve.html?id={pid}"


# --------------------------------------------------------------------------- main

def process(path: Path, args) -> str:
    name = path.stem
    log(f"== {path.name}")
    sub = read_json(path)
    jobdir = WORK / name
    if jobdir.exists():
        shutil.rmtree(jobdir)
    jobdir.mkdir(parents=True)
    v = None
    try:
        v = validate(sub)
        curves = existing_curves()
        existing = [{"id": c["id"], "f": c["curve"]["f"], "h": c["curve"]["h"]} for c in curves]
        job = write_job(v, jobdir, args.mem_gb, existing)
        if args.dry_run:
            log(f"  [dry-run] job written to {job}")
            return "dry"
        rc, timed_out = run_magma(job, args.timeout)
        result_file = jobdir / "result.json"
        if timed_out or not result_file.exists():
            raise Reject("Magma did not finish within the time limit "
                         f"({args.timeout} s); the torsion computation or the certificates took too long")
        res = denull(read_json(result_file))
        if not res.get("ok"):
            raise Reject("verification failed: " + res.get("error", "unknown error"))
        mlog = (jobdir / "magma.log").read_text() if (jobdir / "magma.log").exists() else ""
        if "VERIFY_DONE" not in mlog:
            raise Reject("Magma did not reach VERIFY_DONE (see log)")
        inv = res["torsion"]["invariants"]
        claim_note = ""
        if v.get("group") is not None and not res["torsion"]["claim_matches"]:
            claim_note = f"the submission claimed J(Q)_tors = {group_bracket(v['group'])}, but Magma computes {group_bracket(inv)}"
            if args.require_claim:
                raise Reject(claim_note)
            log(f"  note: {claim_note}; recorded with the computed group")
        if res.get("duplicates_over_Q"):
            raise Reject(f"the curve is isomorphic over Q to the census curve {res['duplicates_over_Q'][0]}")
        # geometrically isomorphic curves (same G2-invariants) already in the census
        g2 = res["curve"]["g2_invariants"]
        res["same_g2"] = [c["id"] for c in curves if c["curve"]["g2_invariants"] == g2]
        if v.get("cover") and not res["split"]["cover_accepted"]:
            log(f"  note: the submitted cover was not accepted: {res['split']['cover_error']}")
        # conductor job
        cond = None
        if not args.no_conductor:
            cjob = write_conductor_job(v, jobdir, min(args.mem_gb, 8))
            crc, ctimed = run_magma(cjob, args.conductor_timeout)
            cfile = jobdir / "conductor.json"
            clog = (jobdir / "conductor.log").read_text() if (jobdir / "conductor.log").exists() else ""
            if not ctimed and cfile.exists() and "CONDUCTOR_DONE" in clog:
                cond = read_json(cfile)
                log(f"  conductor {cond['conductor']}")
            else:
                log("  conductor not computed" + (" (timed out)" if ctimed else ""))
        lmfdb = lmfdb_record(v, res, not args.no_lmfdb)
        log(f"  LMFDB: {lmfdb['label'] or '-'} ({lmfdb['kind']})")
        now = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
        pid = next_id(group_key(inv), curves)
        cert = build_certificate(v, res, cond, curves, pid, now, lmfdb)
        write_json(CURVES_DIR / f"{pid}.json", cert)
        LOGS_DIR.mkdir(exist_ok=True)
        shutil.copy(jobdir / "magma.log", LOGS_DIR / f"{pid}.log")
        if cond and (jobdir / "conductor.log").exists():
            shutil.copy(jobdir / "conductor.log", LOGS_DIR / f"{pid}.conductor.log")
        log(f"  ACCEPTED as {pid} [{cert['status']}] {group_bracket(inv)} class={cert['class']}"
            + (" NEW GROUP" if cert["new_group"] else " new for this class" if cert["new_for_class"] else ""))
        cls_line = (f"- **class: {cert['class_text']}**" + (" (certified)" if cert["class_certified"] else " — not certified; the torsion is"))
        cert_txt = ""
        if cert["simplicity"]["geometrically_simple"]:
            cert_txt = f"geometrically simple: strict prime {cert['simplicity']['strict_prime']}, χ = {cert['simplicity']['chi']}"
        elif cert["split"]["certificate"]:
            cert_txt = cert["split"]["certificate"]["detail"]
        body = (f"Verified on Mordell and added to the census as **{pid}** (status {cert['status']}).\n\n"
                f"- **J(ℚ)tors = {group_bracket(inv)}** (order {cert['order']})\n{cls_line}\n"
                + (f"- {cert_txt}\n" if cert_txt else "")
                + (f"- LMFDB: {lmfdb['url']}\n" if lmfdb['url'] else "")
                + (f"- conductor {cond['conductor']}\n" if cond else "")
                + (f"\nNote: {claim_note}.\n" if claim_note else "")
                + ("\nThis is the first curve in the census with this torsion group.\n" if cert["new_group"]
                   else "\nThis is the first curve in the census with this torsion group in this class.\n" if cert["new_for_class"] else "")
                + f"\nPage: {site_url(pid)}")
        github_feedback(v, cert["status"], body, args.no_github)
        outcome = "accepted"
    except Reject as e:
        reason = str(e)
        log(f"  REJECTED: {reason}")
        rec = {"schema": SCHEMA_REJECTED, "source_file": path.name, "reason": reason, "submission": sub,
               "date": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")}
        if (jobdir / "magma.log").exists():
            rec["magma_log"] = (jobdir / "magma.log").read_text()[-4000:]
        write_json(REJECTED_DIR / f"{name}.json", rec)
        if v is not None:
            github_feedback(v, "rejected", f"The submission could not be certified: {reason}\n\n"
                            "You are welcome to correct the data and submit again.", args.no_github)
        outcome = "rejected"
    if not args.keep and not args.dry_run:
        SUBMISSIONS_PROCESSED.mkdir(parents=True, exist_ok=True)
        shutil.move(str(path), str(SUBMISSIONS_PROCESSED / f"{name}.{outcome}.json"))
    return outcome


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*")
    ap.add_argument("--timeout", type=int, default=3600)
    ap.add_argument("--mem-gb", type=int, default=16)
    ap.add_argument("--conductor-timeout", type=int, default=600)
    ap.add_argument("--no-conductor", action="store_true")
    ap.add_argument("--no-lmfdb", action="store_true")
    ap.add_argument("--no-github", action="store_true")
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--require-claim", action="store_true", help="reject when the computed group differs from the claimed one (validation runs)")
    args = ap.parse_args()
    files = [Path(f) for f in args.files] or sorted(SUBMISSIONS_INBOX.glob("*.json"))
    if not files:
        log("nothing to verify")
        return
    WORK.mkdir(exist_ok=True)
    LOGS.mkdir(exist_ok=True)
    CURVES_DIR.mkdir(parents=True, exist_ok=True)
    REJECTED_DIR.mkdir(parents=True, exist_ok=True)
    summary = {}
    for f in files:
        summary[f.name] = process(f, args)
    log("summary: " + ", ".join(f"{k}: {v}" for k, v in summary.items()))
    n_rej = sum(1 for v in summary.values() if v == "rejected")
    log(f"{len(summary)} processed, {n_rej} rejected")


if __name__ == "__main__":
    main()
