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
| `family_312_certificate.log` | produced here by `pipeline/magma/family_312_certificate.m` | — | exact infinitude of geometrically simple [3,12]: the line (A,B) = (−2/3, −10) of the A(12) chart, the rank-1 elliptic curve w² = (3z²+25)(3z²+121) and the congruence argument at 37 and 41 (found 2026-10-12, torsion_jac/notes/m312_a12_sieve_2026_10_12.md) |
| `family_19_certificate.log` | produced here by `pipeline/magma/family_19_certificate.m` | — | exact infinitude of geometrically simple [19]: the rational family F1 with a class of order 19 over ℚ(r) (census 19.a is r = −1), and the residue class r ≡ 2 (mod 3), where #J(𝔽₃) = 19 with a strict Frobenius (derived 2026-10-13, torsion_jac/notes/z38_2026_10_13.md) |
| `family_13_certificate.log` | produced here by `pipeline/magma/family_13_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [13]: relation-lattice family (3,2),(1,5), residue class t ≡ 1 (mod 3) with #J(𝔽₃) = 13, strict (found in the notes review of 2026-10-13, torsion_jac/notes/review_of_notes_2026_10_13.md) |
| `family_15_certificate.log` | produced here by `pipeline/magma/family_15_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [15]: relation-lattice family (5,0),(−3,3), t ≡ 1 (mod 3) with #J(𝔽₃) = 15, strict (found in the notes review of 2026-10-13, torsion_jac/notes/review_of_notes_2026_10_13.md) |
| `family_16_certificate.log` | produced here by `pipeline/magma/family_16_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [16]: relation-lattice family (3,2),(−2,4), t ≡ 2 (mod 3) with #J(𝔽₃) = 16, strict (found in the notes review of 2026-10-13, torsion_jac/notes/review_of_notes_2026_10_13.md) |
| `family_17_certificate.log` | produced here by `pipeline/magma/family_17_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [17]: relation-lattice family (3,2),(−1,5), t ≡ 2 (mod 3) with #J(𝔽₃) = 17, strict (found in the notes review of 2026-10-13, torsion_jac/notes/review_of_notes_2026_10_13.md) |
| `family_24_certificate.log` | produced here by `pipeline/magma/family_24_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [24]: relation-lattice family (4,1),(0,6) over a rational curve of quadratics, t ≡ 9 (mod 35) (found in the notes review of 2026-10-13, torsion_jac/notes/review_of_notes_2026_10_13.md) |
| `family_25_certificate.log` | produced here by `pipeline/magma/family_25_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [25]: relation-lattice family (5,0),(1,5) over the same curve, t ≡ 3 (mod 77) (found in the notes review of 2026-10-13, torsion_jac/notes/review_of_notes_2026_10_13.md) |
| `family_26_certificate.log` | produced here by `pipeline/magma/family_26_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [26]: the rational root-locus of the [13] family, t ≡ 17 (mod 35) (found in the notes review of 2026-10-13, torsion_jac/notes/review_of_notes_2026_10_13.md) |
| `family_214_certificate.log` | produced here by `pipeline/magma/family_214_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,14]: the contact-7 two-root surface, slice s = 1, t ≡ 2 (mod 65) (found in the notes review of 2026-10-13, torsion_jac/notes/review_of_notes_2026_10_13.md) |
| `family_216_certificate.log` | produced here by `pipeline/magma/family_216_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,16]: the [16] family at s = 3/4 + u², u ≡ 1 (mod 5) with #J(𝔽₅) = 32, strict (found in the notes review of 2026-10-13, torsion_jac/notes/review_of_notes_2026_10_13.md) |
| `family_228_certificate.log` | produced here by `pipeline/magma/family_228_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,2,8]: the odd M₁(8,2,2) family of the lab inventory, slice v = 3, u ≡ 11 (mod 371) (found in the notes review of 2026-10-13, torsion_jac/notes/review_of_notes_2026_10_13.md) |
| `family_2214_certificate.log` | produced here by `pipeline/magma/family_2214_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,2,14]: the contact-7 three-root locus over s₁ = −1/2, a genus-1 curve birational to y² = x³ − 12x (rank 1), congruence at 17 and 29 (found in the notes review of 2026-10-13, torsion_jac/notes/review_of_notes_2026_10_13.md) |
| `family_244_certificate.log` | produced here by `pipeline/magma/family_244_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,4,4]: Nicholls' family Λ₃₃₄ (member formula from the Genus2TorsionNotebook), its rational l₃-split locus, slice t = m = 2, [2,4,4] marked over ℚ(v) by two halving identities; v ≡ 2 (mod 589), #J(𝔽₁₉) = 288, #J(𝔽₃₁) = 1088, gcd 32, strict at 31 (found in the review of the Genus2TorsionNotebook, 2026-10-14) |
| `family_88_split_certificate.log` | produced here by `pipeline/magma/family_88_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [8,8]: the conjugation self-gluing of the X₁(8) curve from the Genus2TorsionNotebook (notes/claude_split_2266_1212_session_2026_08_13.md, §2.1), C_m: y² = F(x²) with explicit degree-2 maps onto E′ ≅ E_t; [8,8] marked over ℚ(m); m ≡ 2 (mod 4381), exactness from the group structures J(𝔽₁₃) = (ℤ/16)² and J(𝔽₃₃₇) = [8,8,40,40] (the gcd of #J(𝔽_p) is always divisible by 256 here) |
| `attribution_check.log` | produced here by `pipeline/magma/attribution_check.m` | — | the comparisons behind the credit errata of 2026-10-14 below: census curves 27.b, 24.b, 25.b, 3.9.b, 6.6.b, 5.5.a, 21.c, 27.c against the curves printed in Leprévost 1995 (C21, C24,1, C24,2, C25, C27, pp. 291–302 of the JTNB scan), Howe 2015 (Table 2, order 27, arXiv:1407.2654v3), HLP 2000 (X₀(22), Cor. 5) and the LMFDB, by G2-invariants and isomorphism over ℚ. It also shows that 21.c (LMFDB 676.a.5408.1, issue #7) is Leprévost's C21 and 27.c (issue #8) is his C27, which is why both are credited to Leprévost 1995 |

## Errata to the seed import: curve credits (review of 2026-10-14)

The Source(s) column of the paper's tables dates the *group*: it gives the first known source of a
realisation, which may only give a family whose torsion contains the group (companion tables
`torsion_jac/paper/old/torsion_realizations.tex` and `split_torsion_table.tex`, "Source(s) column"; in
the split companion a dagger marks such containment-only sources).  `import_paper.py` also copied it
into each curve's `discovery.by`.  The site keeps the column in `sources` (it dates the group, see the
about page), and `discovery.by` now credits the curve itself.  `pipeline/import_paper.py`
(`CURVE_CREDIT`, `SOURCES_OVERRIDE`) assigns these credits, and `test_seed_credit_matches_import` keeps the 150 seed
certificates consistent with it.  Corrected:

* **27.b** (7776.e.2), credited to Leprévost 1995: this is not Leprévost's C27 (G2-invariants
  [−173781261/512, 2737575/128, 106329/32] against [−1564031349/1024, 171810207/512, 10094907/256]).
  It is isomorphic over ℚ to y² + (6x³+9x²+6x−1)y = −3x², marked "New" in Table 2 of Howe 2015 (the
  fourth order-27 row).  Both `discovery.by` and `sources` are now Howe 2015.  C27 itself is split over ℚ:
  Leprévost's degree-3 map, p. 302, goes onto 54b3.  It has been in the census since 2026-10-14 as
  27.c, so Leprévost 1995 does not date the class "simple over ℚ, split over ℚ̄" of [27].  The paper's
  group-level credit of the split [27] to Leprévost 1995 is correct.
* **24.b** (300.a.1), credited to Leprévost 1995: this is neither C24,1 nor C24,2 (G2-invariants above).
  It is a curve of the extended database; no production LMFDB curve has conductor 300.  `discovery.by`
  is now Booker–Sutherland (extended database).  `sources` is unchanged: C24,2 is split over ℚ (its
  involution, p. 301), so Leprévost 1995 still dates the ℚ-split [24].
* **HLP 2000 as the discoverer of 17 database or new curves.** The split companion reads
  "HLP2000†, BookerSutherland" for [20], [30], [35], [40], [60], [2,24], [3,12], [5,10], [2,4,8] and
  [2,6,6], and "HLP2000†, BSSVY" for [2,12] and [6,6].  It reads "HLP2000†, new" for [45], [6,12],
  [8,8], [2,2,24] and [2,2,4,8].  `discovery.by` is now:
  * Booker–Sutherland (extended database) for 20.b, 30.b, 35.a, 40.b, 60.a, 2.24.a, 3.12.a, 5.10.a,
    2.4.8.b and 2.6.6.a;
  * LMFDB genus 2 database (BSSVY 2016) for 2.12.b and 6.6.b;
  * Balakrishnan–Najman–Shnidman–Sutherland 2026 for 45.a, 6.12.a, 8.8.a, 2.2.24.a and 2.2.4.8.a.

  `sources` (HLP 2000, group level) is unchanged for all of them.  6.6.b is X₀(28)
  (`SmallModularCurve(28)`); an earlier source for that curve's torsion was not looked for (UNVERIFIED).
* **3.9.b** (2916.e.2), credited to HLP 2000: the companion reads "Leprevost1995†, BSSVY", and the
  curve is isomorphic over ℚ to the production curve 2916.a.139968.1.  `discovery.by` is now LMFDB
  genus 2 database (BSSVY 2016).  `sources` is unchanged.  For [3,9], the paper's Table 2 cites HLP 2000 and its
  companion cites Leprévost 1995: this group-level difference is left as in `table2.txt`.

Checked and correct: 25.b is Leprévost's C25, and 5.5.a is X₀(22) (HLP 2000, Corollary 5 and the remark after it).
Not checked: 7.7.a and 63.a, which the companion credits to HLP 2000 without a dagger.  Left unchanged
(UNVERIFIED): the other 130 seed rows keep the Source column as `discovery.by`.  On the simple side, the
companion states that the displayed curves are smallest-conductor database curves while the column dates
the group.  So many Table 1 credits, for example 10.a "Flynn 1991" for 389.a.389.1, are probably
group-level too.  Deciding which database curves are published curves needs the sources' equations,
and not all of them are available here.  Note that the site shows `discovery.by` only for curves with
empty `sources`, so the corrected credits are in `data/curves.json` but not on the seed curves' pages.

## Errata to the copied record (`census_infinity_record.md`, copied verbatim; found by the audit of 2026-09-18, see `notes/`)

* §2.1: the [2] family `y² = x(x⁵ + x + t)` has generic 2-rank 0 (the quintic is irreducible over ℚ(t)); the fiber t = 2 is exact [2] only because x⁵ + x + 2 has the root −1. The site uses `y² = x(x⁴ + x + t)` instead.
* §3, [15]: "inherited from the cyclic [30] family" gives torsion ⊇ [30] ⊇ [15] only, not exact [15]; graded `family` on the site.
* §2.7: the marked divisor is described by `Q = 0, y = 0` on the model `y² = g`, on which the points over the roots of Q have y = ±x; the description mixes the coordinates of two models (not verified further here).
* §3, [4,8]: the record lists the tangent-cover family as certified, but there is no fiber certificate for it, and its only displayed simplicity witness (1849T⁴+38T²+1) is not strict; graded open on the site since 2026-10-13 (notes review).
* §4, [2,2,2,12]: the discussion of infinitude is superseded by the paper's Theorem 3.3 and the exact-infinitude certificate above.
* §1: "excess torsion is a closed condition" is not the right justification for infinitely many exact fibers; see the about page for the residue-class argument.

The LMFDB labels in `table1.txt`/`table2.txt` are of two kinds: production labels
`cond.class.disc.num` (permalinks on www.lmfdb.org) and extended-database labels
`cond.class.num` (a 2026-08 snapshot of alpha.lmfdb.org, not permanent).  The
pipeline re-resolves production labels by asking www.lmfdb.org to look the curve
up by its equation (`/Genus2Curve/Q/?jump=[[f],[h]]`) and records the label it
redirects to; extended-database rows are linked by equation.
