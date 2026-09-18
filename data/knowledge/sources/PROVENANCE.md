# Provenance of the copied source tables

Every file in this directory is a verbatim copy (or a verbatim excerpt) of a
source outside this repository.  Nothing in it was edited; the facts the site
derives from them are in `pipeline/knowledge.py`, which cites them by key.

| file | copied from | version | what the site takes from it |
|---|---|---|---|
| `table1.txt` | https://github.com/AndrewVSutherland/Genus2Torsion (`table1.txt`), the certification repository of Balakrishnan–Najman–Shnidman–Sutherland, *Rational torsion on simple genus two Jacobians*, arXiv:2608.28543 | commit `508852a150e17c00b84b209fee23300df897871d` (2026-08-30) | the 73 seed curves with geometrically simple Jacobian (Table 1 of the paper): group, model `[[f],[h]]`, LMFDB label (2026-08 snapshot), source citation |
| `table2.txt` | same repository (`table2.txt`) | same commit | the 77 seed curves with geometrically split Jacobian (Table 2 of the paper), and the split-certificate route (`R1`, `Q2`, `COVER:n`) used by `verify_split_certificates.m`; the explicit covers for the 18 `COVER` rows were transcribed into `pipeline/import_paper.py` from that script |
| `Genus2Torsion_LICENSE` | same repository | same commit | licence of the copied data (MIT) |
| `HLP2000_Theorem1_Table1.txt` | E. W. Howe, F. Leprévost, B. Poonen, *Large torsion subgroups of split Jacobians of curves of genus two or three*, Forum Math. 12 (2000), 315–364; text extracted (pypdf) from arXiv:math/9809210v2, pages 1–3 | arXiv v2 (19 September 1998) | Theorem 1 and Table 1: the 20 groups G with a family of genus-2 curves over ℚ whose Jacobians *contain* G, and the type of the parameter space (the split-side "∞ ⊇" grade) |
| `census_infinity_record.md` | `/home/fnajman/torsion_jac/notes/census_infinity_record.md` (private repository `F-Najman/torsion_jac`) | commit `4d51d0d62a770bb83e9b37fc68ddb8d0a123ec20` (2026-08-21) | the group-by-group record of certified infinite families of geometrically simple Jacobians with exact torsion G (criteria L1–L3 stated in its §1); scripts `code/claude_census_infty_*.m` and logs live in that repository. Credited on the site as "Najman 2026 (unpublished)". |
| `trivial_torsion_family.log` | produced here by `pipeline/magma/trivial_torsion_family.m` (Magma) | — | the L1–L3 certificate for the trivial group: fibers of `y² = x⁵ + x + t` with trivial torsion, a strict prime each, and pairwise distinct G2-invariants |
| `family_2_certificate.log` | produced here by `pipeline/magma/family_2_certificate.m` | — | the certificate for [2] with the family `y² = x(x⁴ + x + t)` (replacing the record's `y² = x(x⁵ + x + t)`, which has no rational 2-torsion generically — see the errata below) |
| `family_22212_certificate.log` | produced here by `pipeline/magma/family_22212_certificate.m` | — | exact infinitude of geometrically simple [2,2,2,12] from the paper's family, by the congruence argument of the audit of 2026-09-18 |

## Errata to the copied record (`census_infinity_record.md`, copied verbatim; found by the audit of 2026-09-18, see `notes/`)

* §2.1: the [2] family `y² = x(x⁵ + x + t)` has generic 2-rank 0 (the quintic is irreducible over ℚ(t)); the fiber t = 2 is exact [2] only because x⁵ + x + 2 has the root −1. The site uses `y² = x(x⁴ + x + t)` instead.
* §3, [15]: "inherited from the cyclic [30] family" gives torsion ⊇ [30] ⊇ [15] only, not exact [15]; graded `family` on the site.
* §2.7: the marked divisor is described by `Q = 0, y = 0` on the model `y² = g`, on which the points over the roots of Q have y = ±x; the description mixes the coordinates of two models (not verified further here).
* §4, [2,2,2,12]: the discussion of infinitude is superseded by the paper's Theorem 3.3 and the exact-infinitude certificate above.
* §1: "excess torsion is a closed condition" is not the right justification for infinitely many exact fibers; see the about page for the residue-class argument.

The LMFDB labels in `table1.txt`/`table2.txt` are of two kinds: production labels
`cond.class.disc.num` (permalinks on www.lmfdb.org) and extended-database labels
`cond.class.num` (a 2026-08 snapshot of alpha.lmfdb.org, not permanent).  The
pipeline re-resolves production labels by asking www.lmfdb.org to look the curve
up by its equation (`/Genus2Curve/Q/?jump=[[f],[h]]`) and records the label it
redirects to; extended-database rows are linked by equation.
