#!/usr/bin/env python3
"""Pull open submission issues from GitHub into submissions/inbox/ as JSON files.

Usage:  python3 pipeline/fetch_issues.py [--dry-run]

Reads the issues labelled `submission` that are still open and not yet queued or processed, parses
the issue-form body (### Label / value blocks) into the submission schema, and writes
submissions/inbox/issue-<number>.json.  Needs a GitHub token (GITHUB_TOKEN or ~/.git-credentials)
only for private repositories or higher rate limits; reading public issues works without one.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (SCHEMA_SUBMISSION, SUBMISSION_LABEL, SUBMISSIONS_INBOX, SUBMISSIONS_PROCESSED,  # noqa: E402
                    github_request, github_token, log, repo_path, write_json)

LABELS = {
    "f(x)": "f", "h(x)": "h", "torsion group": "group", "class of the jacobian": "class",
    "splitness certificate (optional)": "cover", "generators of the torsion subgroup (optional, recommended)": "generators",
    "discovered by": "discoverer", "year of discovery": "year",
    "reference": "reference", "your name and affiliation": "name", "notes": "notes",
}
DONE_LABELS = {"certified", "verified", "rejected"}


def parse_body(body: str) -> dict:
    fields = {}
    parts = re.split(r"^### (.+?)\s*$", body.replace("\r\n", "\n"), flags=re.M)
    for i in range(1, len(parts) - 1, 2):
        label = parts[i].strip().lower()
        value = parts[i + 1].strip()
        if value == "_No response_":
            value = ""
        if label in LABELS:
            fields[LABELS[label]] = value
    return fields


def coeff_or_poly(s: str):
    """'[0,-1,0,16]' -> list, else the polynomial string."""
    s = s.strip()
    if s.startswith("["):
        try:
            v = json.loads(s)
            if isinstance(v, list):
                return v
        except Exception:
            pass
    return s


def to_submission(issue: dict) -> dict:
    f = parse_body(issue.get("body") or "")
    fv, hv = coeff_or_poly(f.get("f", "")), coeff_or_poly(f.get("h", "") or "0")
    sub = {"schema": SCHEMA_SUBMISSION}
    if isinstance(fv, list) or isinstance(hv, list):
        sub["coeffs"] = [fv if isinstance(fv, list) else [], hv if isinstance(hv, list) else []]
        sub["f"], sub["h"] = (fv if isinstance(fv, str) else ""), (hv if isinstance(hv, str) else "")
    else:
        sub["f"], sub["h"] = fv, hv
    grp = f.get("group", "").strip()
    try:
        sub["group"] = [int(t) for t in grp.strip("[] ").split(",") if t.strip()] if grp else None
    except ValueError:
        sub["group"] = grp
    cls = f.get("class", "").strip().lower()
    sub["class"] = "simple" if cls.startswith("geometrically simple") else "split" if cls.startswith("geometrically split") else ""
    gens = f.get("generators", "").strip()
    if gens:
        gens = re.sub(r"^```(json)?\s*|\s*```$", "", gens)
        try:
            sub["generators"] = json.loads(gens)
        except Exception:
            sub["generators"] = [["invalid", gens, 0]]
    cov = f.get("cover", "").strip()
    if cov:
        cov = re.sub(r"^```(json)?\s*|\s*```$", "", cov)
        try:
            sub["cover"] = json.loads(cov)
        except Exception:
            sub["cover"] = {"invalid": cov}
    name = f.get("name", "").strip()
    sub["submitter"] = name.split(",")[0].strip() if name else (issue["user"]["login"] if issue.get("user") else "")
    sub["affiliation"] = ",".join(name.split(",")[1:]).strip() if "," in name else ""
    sub["github"] = issue["user"]["login"] if issue.get("user") else ""
    sub["reference"] = f.get("reference", "").strip()
    sub["discoverer"] = f.get("discoverer", "").strip()
    sub["year"] = int(f["year"].strip()) if f.get("year", "").strip().isdigit() else f.get("year", "").strip()
    sub["notes"] = f.get("notes", "").strip()
    sub["date"] = (issue.get("created_at") or "")[:10]
    sub["source"] = {"kind": "issue", "number": issue["number"], "url": issue.get("html_url", "")}
    return sub


def already_known(number: int) -> bool:
    for d in (SUBMISSIONS_INBOX, SUBMISSIONS_PROCESSED):
        if list(d.glob(f"issue-{number}.*json")) or (d / f"issue-{number}.json").exists():
            return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    tok = github_token()
    issues = []
    page = 1
    while True:
        batch = github_request("GET", repo_path("/issues"), token=tok,
                               params={"labels": SUBMISSION_LABEL, "state": "open", "per_page": 100, "page": page})
        if not batch:
            break
        issues.extend(i for i in batch if "pull_request" not in i)
        if len(batch) < 100:
            break
        page += 1
    SUBMISSIONS_INBOX.mkdir(parents=True, exist_ok=True)
    n_new = 0
    for issue in issues:
        labels = {l["name"] for l in issue.get("labels", [])}
        if labels & DONE_LABELS or already_known(issue["number"]):
            continue
        sub = to_submission(issue)
        target = SUBMISSIONS_INBOX / f"issue-{issue['number']}.json"
        if args.dry_run:
            log(f"[dry-run] would write {target}: group {sub.get('group')} by {sub['submitter']}")
        else:
            write_json(target, sub)
            log(f"queued issue #{issue['number']} -> {target.name}")
        n_new += 1
    log(f"{len(issues)} open submission issue(s), {n_new} newly queued")


if __name__ == "__main__":
    main()
