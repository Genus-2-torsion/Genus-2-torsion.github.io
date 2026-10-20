#!/usr/bin/env python3
"""Recompute the minimal model and conductor of census curves whose conductor job did not finish.

Usage
-----
  python3 pipeline/fill_conductors.py [--timeout SECONDS] [--mem-gb N] [--workdir DIR] [--reuse] ID [ID...]

For each curve this runs pipeline/magma/conductor_lib.m as verify.py does, under its own time limit
(default 6 hours instead of verify.py's 10 minutes), and rewrites the minimal-model and conductor fields
of data/curves/<ID>.json with verify.apply_conductor, by the same rules as for a new submission.  The
log is kept as data/logs/<ID>.conductor.log whenever the minimal model was computed.  With --reuse a
finished job already in WORKDIR/<ID> (conductor.m, conductor.log, conductor.json, conductor.stdout) is
read instead of rerun; a job that neither finished nor stopped with a Magma error is taken as timed out.
Rebuild afterwards with `python3 -B pipeline/build.py`.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import verify  # noqa: E402
from common import CURVES_DIR, LOGS_DIR, WORK, log, read_json, write_json  # noqa: E402


def fill(cid: str, workdir: Path, timeout: int, mem_gb: int, reuse: bool) -> None:
    path = CURVES_DIR / f"{cid}.json"
    cert = read_json(path)
    jobdir = workdir / cid
    jobdir.mkdir(parents=True, exist_ok=True)
    clog_p, cfile, cjob = jobdir / "conductor.log", jobdir / "conductor.json", jobdir / "conductor.m"
    timed_out = False
    if not (reuse and clog_p.exists()):
        for p in (clog_p, cfile):
            p.unlink(missing_ok=True)
        cjob = verify.write_conductor_job({"f": cert["curve"]["f"], "h": cert["curve"]["h"]}, jobdir, mem_gb)
        _, timed_out = verify.run_magma(cjob, timeout)
    clog = clog_p.read_text() if clog_p.exists() else ""
    partial, missing = None, ""
    cond = None if timed_out else verify.conductor_result(cfile, clog)
    if not cond:
        error = verify.magma_error(cjob.with_suffix(".stdout"))
        missing = verify.conductor_missing_reason(timed_out or not error, error)
        if cfile.exists() and "MINIMAL_MODEL_DONE" in clog:
            partial = read_json(cfile)
    verify.apply_conductor(cert["curve"], cond, cert.get("lmfdb") or {}, cert["split"].get("isogenous_product"), missing, partial)
    ver = cert["verification"]
    if (cond or partial) and clog_p.exists():
        LOGS_DIR.mkdir(exist_ok=True)
        shutil.copy(clog_p, LOGS_DIR / f"{cid}.conductor.log")
        ver["conductor_log"] = f"data/logs/{cid}.conductor.log"
    ver["conductor_cputime_seconds"] = cond["cputime"] if cond else None
    write_json(path, cert)
    log(f"{cid}: minimal model {'yes' if cert['curve']['minimal_model'] else 'no'}; conductor {cert['curve']['conductor'] or '-'} "
        f"({cert['curve']['conductor_source'] or 'none'})" + (f"; {missing}" if missing else ""))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ids", nargs="+")
    ap.add_argument("--timeout", type=int, default=6 * 3600)
    ap.add_argument("--mem-gb", type=int, default=16)
    ap.add_argument("--workdir", default=str(WORK / "conductors"))
    ap.add_argument("--reuse", action="store_true")
    args = ap.parse_args()
    for cid in args.ids:
        if not (CURVES_DIR / f"{cid}.json").exists():
            sys.exit(f"no census curve {cid}")
    for cid in args.ids:
        fill(cid, Path(args.workdir), args.timeout, args.mem_gb, args.reuse)


if __name__ == "__main__":
    main()
