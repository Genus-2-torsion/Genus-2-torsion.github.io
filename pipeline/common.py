"""Shared helpers for the genus-2 torsion census pipeline (paths, GitHub API, JSON I/O, groups)."""

from __future__ import annotations

import json
import os
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
CURVES_DIR = DATA / "curves"
LOGS_DIR = DATA / "logs"
REJECTED_DIR = DATA / "rejected"
KNOWLEDGE_DIR = DATA / "knowledge"
SUBMISSIONS_INBOX = ROOT / "submissions" / "inbox"
SUBMISSIONS_PROCESSED = ROOT / "submissions" / "processed"
PIPELINE = ROOT / "pipeline"
WORK = PIPELINE / "work"
LOGS = PIPELINE / "logs"

GITHUB_ORG = "Genus-2-torsion"
GITHUB_REPO = "Genus-2-torsion.github.io"
GITHUB_API = "https://api.github.com"
SUBMISSION_LABEL = "submission"
SITE_URL = "https://genus-2-torsion.github.io"

SCHEMA_SUBMISSION = "genus-2-torsion/submission/1"
SCHEMA_CERTIFICATE = "genus-2-torsion/curve/1"
SCHEMA_REJECTED = "genus-2-torsion/rejected/1"

# The three classes of the census (mutually exclusive), plus the two undecided states.
CLASSES = {
    "simple": "geometrically simple",
    "qsplit": "split over ℚ",
    "gsplit": "geometrically split, simple over ℚ",
}
UNDECIDED = {
    "split_undecided_over_Q": "geometrically split; whether it splits over ℚ is undecided",
    "undecided": "neither a simplicity nor a splitness certificate was found",
}

# Whitelists for strings that reach Magma's `eval`.  Nothing else is ever evaluated.
RE_RATIONAL = re.compile(r"^-?\d+(/\d+)?$")
RE_POLY_X = re.compile(r"^[0-9x\s+\-*/^()]+$")            # a polynomial in x with rational coefficients
RE_POLY_X_UPPER = re.compile(r"^[0-9X\s+\-*/^()]+$")      # field polynomial in X
RE_ELEMENT_W = re.compile(r"^[0-9w\s+\-*/^()]+$")         # element of the field, polynomial in w
MAX_STRING = 20000


def read_json(path: Path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, ensure_ascii=False, sort_keys=False)
        fh.write("\n")
    os.replace(tmp, path)


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


# --------------------------------------------------------------------------- groups

def group_key(inv) -> str:
    """[2, 2, 12] -> '2.2.12'; the trivial group [] -> '1'."""
    inv = [int(n) for n in inv]
    return ".".join(str(n) for n in inv) if inv else "1"


def key_to_group(key: str) -> list[int]:
    return [] if key == "1" else [int(t) for t in key.split(".")]


def group_label(inv) -> str:
    """Human label: [2, 2, 12] -> 'ℤ/2 ⊕ ℤ/2 ⊕ ℤ/12'; [] -> 'trivial'."""
    inv = [int(n) for n in inv]
    return " ⊕ ".join(f"ℤ/{n}" for n in inv) if inv else "trivial"


def group_bracket(inv) -> str:
    """The paper's notation: [2,2,12]; [] for the trivial group."""
    return "[" + ",".join(str(int(n)) for n in inv) + "]"


def group_order(inv) -> int:
    o = 1
    for n in inv:
        o *= int(n)
    return o


def group_sort_key(inv):
    """Table order of the paper: by number of invariant factors, then lexicographically."""
    inv = [int(n) for n in inv]
    return (len(inv), inv)


def is_invariant_factors(inv) -> bool:
    inv = [int(n) for n in inv]
    return all(n >= 2 for n in inv) and all(inv[i + 1] % inv[i] == 0 for i in range(len(inv) - 1)) and len(inv) <= 4


# --------------------------------------------------------------------------- polynomials

def parse_poly(s: str) -> list[Fraction]:
    """Parse a polynomial in x with rational coefficients (Magma/Sage style, e.g. '3*x^5 - x/2 + 1')
    into ascending coefficients.  Only + - * / ^ ( ) digits and x are accepted, and the grammar is
    restricted to sums of terms  [(]c[)] [*] [x[^k]] [/d]  so that nothing is ever evaluated."""
    s = s.replace(" ", "").replace("**", "^")
    if not s:
        raise ValueError("empty polynomial")
    if not RE_POLY_X.match(s):
        raise ValueError("polynomial contains characters outside the allowed set (digits, x, + - * / ^ ( ))")
    term_re = re.compile(r"(?:\(?(?P<coef>-?\d+(?:/\d+)?)\)?)?(?:\*?(?P<x>x(?:\^(?P<exp>\d+))?))?(?:/(?P<den>\d+))?")
    coeffs: dict[int, Fraction] = {}
    pos = 0
    sign = 1
    expect_term = True
    while pos < len(s):
        c = s[pos]
        if c == "+" and not expect_term:
            sign = 1; pos += 1; expect_term = True; continue
        if c == "-" and not expect_term:
            sign = -1; pos += 1; expect_term = True; continue
        if c == "-" and expect_term:
            sign = -sign; pos += 1; continue
        if not expect_term:
            raise ValueError(f"expected + or - at position {pos}")
        m = term_re.match(s, pos)
        if not m or not m.group(0) or (m.group("coef") is None and m.group("x") is None):
            raise ValueError(f"cannot parse polynomial near position {pos}")
        coef = Fraction(m.group("coef")) if m.group("coef") is not None else Fraction(1)
        if m.group("den") is not None:
            coef /= int(m.group("den"))
        deg = (int(m.group("exp")) if m.group("exp") else 1) if m.group("x") else 0
        coeffs[deg] = coeffs.get(deg, Fraction(0)) + sign * coef
        pos = m.end()
        sign = 1
        expect_term = False
    if expect_term:
        raise ValueError("polynomial ends with an operator")
    if not coeffs:
        return [Fraction(0)]
    n = max(coeffs)
    return [coeffs.get(i, Fraction(0)) for i in range(n + 1)]


def coeffs_to_strings(coeffs) -> list[str]:
    out = [str(Fraction(c)) for c in coeffs]
    while len(out) > 1 and out[-1] == "0":
        out.pop()
    return out


def poly_to_string(coeffs, var: str = "x") -> str:
    """Ascending coefficients (strings or numbers) -> 'a_n*x^n + ... + a_0' in Magma style."""
    cs = [Fraction(str(c)) for c in coeffs]
    terms = []
    for i in range(len(cs) - 1, -1, -1):
        c = cs[i]
        if c == 0:
            continue
        mag = abs(c)
        if i == 0:
            body = str(mag)
        else:
            xp = var if i == 1 else f"{var}^{i}"
            body = xp if mag == 1 else f"{mag}*{xp}"
        terms.append((c < 0, body))
    if not terms:
        return "0"
    s = ("-" if terms[0][0] else "") + terms[0][1]
    for neg, body in terms[1:]:
        s += (" - " if neg else " + ") + body
    return s


def coeff_lists_from_submission(sub: dict) -> tuple[list[str], list[str]]:
    """Accept either {"f": "...", "h": "..."} polynomial strings or {"coeffs": [[f],[h]]} ascending lists."""
    if "coeffs" in sub and sub["coeffs"]:
        fh = sub["coeffs"]
        if not (isinstance(fh, list) and len(fh) == 2 and all(isinstance(v, list) for v in fh)):
            raise ValueError("coeffs must be [[f],[h]] with ascending coefficient lists")
        f = [Fraction(str(c)) for c in fh[0]]
        h = [Fraction(str(c)) for c in fh[1]] if fh[1] else [Fraction(0)]
    else:
        f = parse_poly(str(sub.get("f", "")))
        hs = str(sub.get("h", "") or "0")
        h = parse_poly(hs)
    return coeffs_to_strings(f), coeffs_to_strings(h)


# --------------------------------------------------------------------------- GitHub

def github_token() -> str | None:
    """GITHUB_TOKEN from the environment, else the github.com entry of ~/.git-credentials."""
    tok = os.environ.get("GITHUB_TOKEN")
    if tok:
        return tok.strip()
    cred = Path.home() / ".git-credentials"
    if cred.exists():
        for line in cred.read_text().splitlines():
            if "github.com" in line:
                m = re.match(r"https://(?:[^:@]*:)?([^@]+)@github\.com", line.strip())
                if m:
                    return m.group(1)
    return None


def github_request(method: str, path: str, data=None, token: str | None = None, params=None):
    url = GITHUB_API + path
    if params:
        url += "?" + urllib.parse.urlencode(params)
    body = None
    headers = {"Accept": "application/vnd.github+json", "User-Agent": "genus-2-torsion-pipeline"}
    if token:
        headers["Authorization"] = f"token {token}"
    if data is not None:
        body = json.dumps(data).encode()
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=body, method=method, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            txt = resp.read().decode()
            return json.loads(txt) if txt else None
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"GitHub {method} {path} -> {e.code}: {e.read().decode()[:500]}") from e


def repo_path(suffix: str) -> str:
    return f"/repos/{GITHUB_ORG}/{GITHUB_REPO}{suffix}"


# --------------------------------------------------------------------------- LMFDB

def lmfdb_jump_url(f: list, h: list, host: str = "www") -> str:
    """The LMFDB 'jump' URL that looks a genus-2 curve up by its equation (any isomorphic model)."""
    eq = "[[" + ",".join(str(c) for c in f) + "],[" + ",".join(str(c) for c in h) + "]]"
    return f"https://{host}.lmfdb.org/Genus2Curve/Q/?jump=" + urllib.parse.quote(eq, safe="")


def lmfdb_lookup(f: list, h: list, timeout: int = 30) -> str | None:
    """Ask www.lmfdb.org for the curve; return its production label if the jump redirects to a
    curve page, else None (not in the production database, or offline)."""
    if any("/" in str(c) for c in f + h):
        return None
    url = lmfdb_jump_url(f, h)

    class NoRedirect(urllib.request.HTTPRedirectHandler):
        def redirect_request(self, req, fp, code, msg, headers, newurl):
            return None

    opener = urllib.request.build_opener(NoRedirect)
    req = urllib.request.Request(url, headers={"User-Agent": "genus-2-torsion-pipeline"})
    try:
        with opener.open(req, timeout=timeout) as resp:
            return None
    except urllib.error.HTTPError as e:
        if e.code in (301, 302, 303, 307, 308):
            loc = e.headers.get("Location", "")
            m = re.search(r"/Genus2Curve/Q/(\d+)/([a-z]+)/(\d+)/(\d+)/?$", loc)
            if m:
                return ".".join(m.groups())
        return None
    except Exception:
        return None
