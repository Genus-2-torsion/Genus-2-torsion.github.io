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

The LMFDB labels in `table1.txt`/`table2.txt` are of two kinds: production labels
`cond.class.disc.num` (permalinks on www.lmfdb.org) and extended-database labels
`cond.class.num` (a 2026-08 snapshot of alpha.lmfdb.org, not permanent).  The
pipeline re-resolves production labels by asking www.lmfdb.org to look the curve
up by its equation (`/Genus2Curve/Q/?jump=[[f],[h]]`) and records the label it
redirects to; extended-database rows are linked by equation.
