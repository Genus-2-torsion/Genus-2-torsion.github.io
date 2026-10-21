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
| `family_2_certificate.log` | produced here by `pipeline/magma/family_2_certificate.m` | — | the certificate for [2] with the family `y² = x(x⁴ + x + t)` (replacing the record's `y² = x(x⁵ + x + t)`, which has no rational 2-torsion generically — see the errata below); residue class t ≡ 2 (mod 15) since 2026-10-18 |
| `family_22212_certificate.log` | produced here by `pipeline/magma/family_22212_certificate.m` | — | exact infinitude of geometrically simple [2,2,2,12] from the paper's family, by the congruence argument of the audit of 2026-09-18; since 2026-10-18 it also checks the inputs from the paper (ψ(E) on the surface, [2,2,2,12] over the function field of E) |
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
| `family_22216_split_certificate.log` | produced here by `pipeline/magma/family_22216_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,2,2,16] (found 2026-10-16): the HLP (2,2)-gluing of E_{2,4}(t) and E_{2,8}(−20/31) on the full-gain fibre t = Aa², t + 1/4 = Bc², t − 1/4 = Ce², a 2-covering of y² = x³ − d²x, d = 554577910 (rank 1); [2,2,2,16] marked over the function field of the fibre; R = P + (d, 0) + 12kP, exactness from J(𝔽₁₁) ≅ J(𝔽₁₃) ≅ [2,2,2,16]; k = 0 is the census curve 2.2.2.16.a |
| `family_14_split_certificate.log` | produced here by `pipeline/magma/family_14_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [14] (found 2026-10-18, torsion_jac/notes/g33fam_2026_10_18.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 26b1 with the partners that have a rational 2-torsion point, a rational family over ℚ(v) with explicit degree-3 maps; v ≡ 0 (mod 497), #J(𝔽₇) = 70, #J(𝔽₇₁) = 4774, gcd 14; v = 0 is the census curve 14.b (26b1 × 52a1) |
| `family_210_split_certificate.log` | produced here by `pipeline/magma/family_210_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,10] (found 2026-10-18, torsion_jac/notes/g33fam_2026_10_18.md): the HLP (2,2)-gluing of the Kubert curve E₁₀(t) with its 5-isogenous curve along the induced isomorphism of 2-torsion (full gain); t ≡ 58 (mod 77), J(𝔽₇) ≅ J(𝔽₁₁) ≅ [10,10], and μ₅ ⊂ E₁[5] with 11 ≡ 1 (mod 5) excludes a second point of order 5 |
| `family_214_split_certificate.log` | produced here by `pipeline/magma/family_214_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,14] (found 2026-10-18, torsion_jac/notes/g33fam_2026_10_18.md): (3,3)-gluings of 294b2 with the partners that have full rational 2-torsion, parametrised by a genus-1 curve with Jacobian y² = x³ + 86436 (rank 1, P = (−392/9, 1666/27)); R₀ + 124kP, #J(𝔽₁₁) = 112, #J(𝔽₉₇) = 8820, gcd 28; v = 0 is the census curve 2.14.b (294b2 × 294c4) |
| `family_2210_split_certificate.log` | produced here by `pipeline/magma/family_2210_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,2,10] (found 2026-10-18, torsion_jac/notes/g33fam_2026_10_18.md): (3,3)-gluings of 66c1 with the partners that have full rational 2-torsion, parametrised by a genus-1 curve with Jacobian y² = x³ + 39204 (rank 1, P = (616, −15290)); R₀ + 48kP, #J(𝔽₂₃) = 720, #J(𝔽₄₇) = 2200, gcd 40; v = 0 is the census curve 2.2.10.b (66c1 × 462f2) |
| `family_2226_split_certificate.log` | produced here by `pipeline/magma/family_2226_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,2,2,6] (found 2026-10-18, torsion_jac/notes/f2226_2026_10_18.md): the HLP (2,2)-gluing of the Kubert curve E_{2,6}(s + 3) with its 3-isogenous curve along the induced isomorphism of 2-torsion (full gain), C_s: y² = (s² + 12)(s² − 4x²)((s − 2)²x² − (s + 6)²)((s + 2)²x² − (s − 6)²); s ≡ 4 (mod 13), J(𝔽₁₃) ≅ [2,2,6,6], and μ₃ ⊂ E₁[3] with 13 ≡ 1 (mod 3) excludes a second point of order 3; s = 4 is the census curve 2.2.2.6.b |
| `family_1_certificate.log` | produced here by `pipeline/magma/family_1_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [ ] (trivial torsion) (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the family y² = x⁵ + x + t on the class t ≡ 3 (mod 51); complements `trivial_torsion_family.log` |
| `family_3_certificate.log` | produced here by `pipeline/magma/family_3_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [3] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's contact form y² = (x³+x+1)² − t(x²+1)³ |
| `family_4_certificate.log` | produced here by `pipeline/magma/family_4_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [4] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's family y² = x·Q₄(x), α = 1 + s√2 |
| `family_5_certificate.log` | produced here by `pipeline/magma/family_5_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [5] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the contact-5 family y² = (1+x+x²)² − tx⁵ (replacing the record's reference, whose chart has torsion ⊇ [10]; errata below) |
| `family_6_certificate.log` | produced here by `pipeline/magma/family_6_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [6] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the contact-6 family h = t + x + x² + x³, y² = h² − (x−1)⁶ (the record's M(12) chart has no residue-class certificate; errata below) |
| `family_7_certificate.log` | produced here by `pipeline/magma/family_7_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [7] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's contact-7 chart, slice b = 1 |
| `family_8_certificate.log` | produced here by `pipeline/magma/family_8_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [8] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's A(8) chart, slice p = 2, r = 3 |
| `family_9_certificate.log` | produced here by `pipeline/magma/family_9_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [9] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's contact-9 chart |
| `family_10_certificate.log` | produced here by `pipeline/magma/family_10_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [10] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's contact-5 chart at b = 1 |
| `family_11_certificate.log` | produced here by `pipeline/magma/family_11_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [11] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): Flynn's family (as quoted by Daowsud–Schmidt) |
| `family_12_certificate.log` | produced here by `pipeline/magma/family_12_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [12] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's rational-Weierstrass subchart of the M(12) chart, slice r = 3 |
| `family_14_certificate.log` | produced here by `pipeline/magma/family_14_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [14] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's rational-root subfamily of the contact-7 chart |
| `family_18_certificate.log` | produced here by `pipeline/magma/family_18_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [18] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's rational-root subfamily of the contact-9 chart |
| `family_20_certificate.log` | produced here by `pipeline/magma/family_20_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [20] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's contact-5 + 4-torsion family (halving component) |
| `family_21_certificate.log` | produced here by `pipeline/magma/family_21_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [21] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): Leprévost's 1991 family f₂₁ |
| `family_22_certificate.log` | produced here by `pipeline/magma/family_22_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [22] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): Flynn's order-11 family on its rational-root curve |
| `family_23_certificate.log` | produced here by `pipeline/magma/family_23_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [23] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the Kuru–Sadek genus-2 family |
| `family_30_certificate.log` | produced here by `pipeline/magma/family_30_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [30] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's contact-5/contact-6 family, branch −1 |
| `family_32_certificate.log` | produced here by `pipeline/magma/family_32_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [32] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's reconstructed Elkies order-32 component, parametrized by w |
| `family_33_certificate.log` | produced here by `pipeline/magma/family_33_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [3,3] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): Bruin–Flynn–Testa 2014, Theorem 6, line r = 2, s = 3 |
| `family_44_certificate.log` | produced here by `pipeline/magma/family_44_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [4,4] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the Richelot quotient of the [2,2,4,4] family by the record's 'partner #5' kernel; group structures at 23, 59, 73 |
| `family_210_certificate.log` | produced here by `pipeline/magma/family_210_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,10] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): Elkies' 5-torsion chart with three rational Weierstrass points, slice (t, 2, 3) |
| `family_220_certificate.log` | produced here by `pipeline/magma/family_220_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,20] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's extra-2 locus of the [20] family |
| `family_224_certificate.log` | produced here by `pipeline/magma/family_224_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,2,4] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's family y² = x(x+4)(x+9)(x² + (2t+1)x + t²) |
| `family_2210_certificate.log` | produced here by `pipeline/magma/family_2210_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,2,10] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): Elkies' 5-torsion chart with four rational Weierstrass points |
| `family_2222_certificate.log` | produced here by `pipeline/magma/family_2222_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,2,2,2] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's family y² = x(x−1)(x−2)(x−3)(x−t) |
| `family_2224_certificate.log` | produced here by `pipeline/magma/family_2224_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,2,2,4] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the square-branch chart of BNSS2026 |
| `family_2226_certificate.log` | produced here by `pipeline/magma/family_2226_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,2,2,6] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's M(2,2,2,6) chart, line s = 7, m = 2 |
| `family_2228_certificate.log` | produced here by `pipeline/magma/family_2228_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,2,2,8] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the rational curve on the K3 surface M(2,2,2,8) of BNSS2026 |
| `family_2244_certificate.log` | produced here by `pipeline/magma/family_2244_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,2,4,4] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the rational curve in A(2,2,4,4) of BNSS2026 |
| `family_22210_certificate.log` | produced here by `pipeline/magma/family_22210_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,2,2,10] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): Elkies' Clebsch–Klein family, tangent section at (1, −8, −7, 5, 9) |
| `family_2_2_certificate.log` | produced here by `pipeline/magma/family_2_2_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,2] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's family y² = x(x−1)(x−2)(x³ + x + t) |
| `family_2_2_2_certificate.log` | produced here by `pipeline/magma/family_2_2_2_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,2,2] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's family y² = x(x−1)(x−2)(x² + x + t) |
| `family_2_4_certificate.log` | produced here by `pipeline/magma/family_2_4_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,4] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's family y² = x(x² + (2t+1)x + t²)(x² + 3x + 1) |
| `family_2_6_certificate.log` | produced here by `pipeline/magma/family_2_6_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,6] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's contact-6 chart with h = 1 + tx + x² + x³ |
| `family_2_8_certificate.log` | produced here by `pipeline/magma/family_2_8_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,8] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the Richelot quotient of the [2,2,2,8] family by the 'partner #4' kernel at t = 2; group structures at 17, 31, 157 (errata below) |
| `family_2_12_certificate.log` | produced here by `pipeline/magma/family_2_12_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,12] (audit of 2026-10-18, torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md): the record's line a = (1−r)/4 in the M(12) chart |
| `family_412_split_certificate.log` | produced here by `pipeline/magma/family_412_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [4,12] (found 2026-10-18, torsion_jac/notes/f412_2026_10_18.md): the HLP (2,2)-gluing of the Kubert curve E₁₂(t) with its 3-isogenous curve along the induced isomorphism of 2-torsion; the gain of order 4 descends from a quartic extension of ℚ(t); J(𝔽₇) ≅ [12,12] on t ≡ 5 (mod 7) and a μ₃ argument; census 4.12.a is t = 1/3 |
| `family_3_6_certificate.log` | produced here by `pipeline/magma/family_3_6_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [3,6] (found 2026-10-18, torsion_jac/notes/simplefam_2026_10_18.md): Bruin–Flynn–Testa's (ℤ/3)² family on a rational curve where F is a quintic with a rational root; class l ≡ 2 (mod 407), gcd 18, strict at 37 |
| `family_2_2_6_certificate.log` | produced here by `pipeline/magma/family_2_2_6_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,2,6] (found 2026-10-18, notes/simplefam_2026_10_18.md): y² = F₁F₂ with F₁ − F₂ = 2t and four rational roots; class t ≡ −2 (mod 145), gcd 24, strict at 29; t = −2 is census 2.2.6.a |
| `family_2_2_12_certificate.log` | produced here by `pipeline/magma/family_2_2_12_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of geometrically simple [2,2,12] (found 2026-10-18, notes/simplefam_2026_10_18.md): the BNSS mechanism F₁ − F₂ = (abc)² on a conic where only S + 2w is a square, with an explicit half of the irrational root pair; class m ≡ 3 (mod 1357), gcd 48, strict at 59 |
| `family_2_2_12_control_certificate.log` | produced here by `pipeline/magma/family_2_2_12_control_certificate.m` (shared argument in `pipeline/magma/family_cert_lib.m`) | — | a second family for simple [2,2,12], over a rank-1 elliptic curve (conductor 2660), through the census curve 2.2.12.a (found 2026-10-18, notes/simplefam_2026_10_18.md) |
| `family_18_split_certificate.log` | produced here by `pipeline/magma/family_18_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [18] (found 2026-10-18, torsion_jac/notes/s18_2026_10_18.md): (3,3)-gluings of 54b3 (torsion ℤ/9) with partners having a rational 2-torsion point, over a rational curve; group structures J(𝔽₅) ≅ ℤ/54, J(𝔽₁₃) ≅ (ℤ/18)² on v ≡ −2 (mod 65); the member v = 4 is census 36.b |
| `family_218_split_certificate.log` | produced here by `pipeline/magma/family_218_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,18] (found 2026-10-18, notes/s18_2026_10_18.md): (3,3)-gluings of 714i1 (torsion ℤ/9) with partners having full rational 2-torsion, over a genus-1 curve whose Jacobian y² = x³ + 289 has rank 1; group structures at 109, 251, 359 |
| `family_1_split_certificate.log` | produced here by `pipeline/magma/family_1_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [1] (found 2026-10-19, torsion_jac/notes/g33rat_2026_10_19.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 37a1 with all partners anti-symplectically 3-congruent to it, over the rational curve X_E^−(3); the residue class and the J(𝔽_p) bound are in the log |
| `family_2_split_certificate.log` | produced here by `pipeline/magma/family_2_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2] (found 2026-10-19, torsion_jac/notes/g33rat_2026_10_19.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 46a1 with all partners anti-symplectically 3-congruent to it, over the rational curve X_E^−(3); the residue class and the J(𝔽_p) bound are in the log |
| `family_3_split_certificate.log` | produced here by `pipeline/magma/family_3_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [3] (found 2026-10-19, torsion_jac/notes/g33rat_2026_10_19.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 19a3 with all partners anti-symplectically 3-congruent to it, over the rational curve X_E^−(3); the residue class and the J(𝔽_p) bound are in the log |
| `family_2_2_split_certificate.log` | produced here by `pipeline/magma/family_2_2_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,2] (found 2026-10-19, torsion_jac/notes/g33rat_2026_10_19.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 17a2 with all partners anti-symplectically 3-congruent to it, over the rational curve X_E^−(3); the residue class and the J(𝔽_p) bound are in the log |
| `family_4_split_certificate.log` | produced here by `pipeline/magma/family_4_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [4] (found 2026-10-19, torsion_jac/notes/g33rat_2026_10_19.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 17a4 with all partners anti-symplectically 3-congruent to it, over the rational curve X_E^−(3); the residue class and the J(𝔽_p) bound are in the log |
| `family_5_split_certificate.log` | produced here by `pipeline/magma/family_5_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [5] (found 2026-10-19, torsion_jac/notes/g33rat_2026_10_19.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 11a3 with all partners anti-symplectically 3-congruent to it, over the rational curve X_E^−(3); the residue class and the J(𝔽_p) bound are in the log |
| `family_6_split_certificate.log` | produced here by `pipeline/magma/family_6_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [6] (found 2026-10-19, torsion_jac/notes/g33rat_2026_10_19.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 14a4 with all partners anti-symplectically 3-congruent to it, over the rational curve X_E^−(3); the residue class and the J(𝔽_p) bound are in the log |
| `family_7_split_certificate.log` | produced here by `pipeline/magma/family_7_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [7] (found 2026-10-19, torsion_jac/notes/g33rat_2026_10_19.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 26b1 with all partners anti-symplectically 3-congruent to it, over the rational curve X_E^−(3); the residue class and the J(𝔽_p) bound are in the log |
| `family_2_2_2_split_certificate.log` | produced here by `pipeline/magma/family_2_2_2_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,2,2] (found 2026-10-19, torsion_jac/notes/g33rat_2026_10_19.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 17a2 with the partners that have a rational 2-torsion point, over a rational curve; the residue class and the J(𝔽_p) bound are in the log |
| `family_2_4_split_certificate.log` | produced here by `pipeline/magma/family_2_4_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,4] (found 2026-10-19, torsion_jac/notes/g33rat_2026_10_19.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 15a1 with all partners anti-symplectically 3-congruent to it, over the rational curve X_E^−(3); the residue class and the J(𝔽_p) bound are in the log |
| `family_8_split_certificate.log` | produced here by `pipeline/magma/family_8_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [8] (found 2026-10-19, torsion_jac/notes/g33rat_2026_10_19.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 15a4 with all partners anti-symplectically 3-congruent to it, over the rational curve X_E^−(3); the residue class and the J(𝔽_p) bound are in the log |
| `family_3_3_split_certificate.log` | produced here by `pipeline/magma/family_3_3_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [3,3] (found 2026-10-19, torsion_jac/notes/isopair_2026_10_19.md): the HLP (2,2)-gluing of the Hesse curve (ℤ/3 ⊕ μ₃ in E[3]) with E/μ₃ along the isogeny; the residue class and the J(𝔽_p) bound are in the log |
| `family_9_split_certificate.log` | produced here by `pipeline/magma/family_9_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [9] (found 2026-10-19, torsion_jac/notes/g33rat_2026_10_19.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 54b3 with all partners anti-symplectically 3-congruent to it, over the rational curve X_E^−(3); the residue class and the J(𝔽_p) bound are in the log |
| `family_10_split_certificate.log` | produced here by `pipeline/magma/family_10_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [10] (found 2026-10-19, torsion_jac/notes/g33rat_2026_10_19.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 66c1 with all partners anti-symplectically 3-congruent to it, over the rational curve X_E^−(3); the residue class and the J(𝔽_p) bound are in the log |
| `family_2_6_split_certificate.log` | produced here by `pipeline/magma/family_2_6_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,6] (found 2026-10-19, torsion_jac/notes/g33rat_2026_10_19.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 30a2 with all partners anti-symplectically 3-congruent to it, over the rational curve X_E^−(3); the residue class and the J(𝔽_p) bound are in the log |
| `family_12_split_certificate.log` | produced here by `pipeline/magma/family_12_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [12] (found 2026-10-19, torsion_jac/notes/g33rat_2026_10_19.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 90c3 with all partners anti-symplectically 3-congruent to it, over the rational curve X_E^−(3); the residue class and the J(𝔽_p) bound are in the log |
| `family_15_split_certificate.log` | produced here by `pipeline/magma/family_15_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [15] (found 2026-10-19, torsion_jac/notes/hlpP2_2026_10_19.md): the HLP (2,2)-gluing of E₅(t) with the curve of HLP Lemma 11 (the [21] construction with E₅ in place of E₇); the residue class and the J(𝔽_p) bound are in the log |
| `family_2_2_2_2_split_certificate.log` | produced here by `pipeline/magma/family_2_2_2_2_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,2,2,2] (found 2026-10-19, torsion_jac/notes/g33ell_2026_10_19.md): bielliptic curves with six rational Weierstrass points and the quotient 45a2 fixed, over a rational curve; the residue class and the J(𝔽_p) bound are in the log |
| `family_2_2_4_split_certificate.log` | produced here by `pipeline/magma/family_2_2_4_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,2,4] (found 2026-10-19, torsion_jac/notes/g33rat_2026_10_19.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 15a1 with the partners that have a rational 2-torsion point, over a rational curve; the residue class and the J(𝔽_p) bound are in the log |
| `family_2_8_split_certificate.log` | produced here by `pipeline/magma/family_2_8_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,8] (found 2026-10-19, torsion_jac/notes/g33rat_2026_10_19.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 210e2 with all partners anti-symplectically 3-congruent to it, over the rational curve X_E^−(3); the residue class and the J(𝔽_p) bound are in the log |
| `family_4_4_split_certificate.log` | produced here by `pipeline/magma/family_4_4_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [4,4] (found 2026-10-19, torsion_jac/notes/isopair_2026_10_19.md): the HLP (2,2)-gluing of an X₁(4) curve with a rational 3-isogeny with E/C₃, with a gain class of order 4; the residue class and the J(𝔽_p) bound are in the log |
| `family_16_split_certificate.log` | produced here by `pipeline/magma/family_16_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [16] (found 2026-10-19, torsion_jac/notes/pow2_2026_10_19.md): the HLP (2,2)-gluing of 15a4 with a one-parameter family of curves with a rational 2-torsion point, with a gain class of order 16 computed over a degree-8 field and descended to ℚ(s); the residue class and the J(𝔽_p) bound are in the log |
| `family_3_6_split_certificate.log` | produced here by `pipeline/magma/family_3_6_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [3,6] (found 2026-10-19, torsion_jac/notes/isopair_2026_10_19.md): the [6,6] gluing twisted by the swap of E[2] fixing the rational 2-torsion point; the 2-rank is fixed by Frobenius elements at 5, 11, 17 acting on the Weierstrass points; the residue class and the J(𝔽_p) bound are in the log |
| `family_20_split_certificate.log` | produced here by `pipeline/magma/family_20_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [20] (found 2026-10-19, torsion_jac/notes/hlpP2_2026_10_19.md): the HLP (2,2)-gluing of 66c1 with E₄(u) on a line of HLP's ℙ²; the residue class and the J(𝔽_p) bound are in the log |
| `family_21_split_certificate.log` | produced here by `pipeline/magma/family_21_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [21] (found 2026-10-19, torsion_jac/notes/hlpP2_2026_10_19.md): the HLP (2,2)-gluing of E₇(t) with the curve of HLP Lemma 11; the residue class and the J(𝔽_p) bound are in the log |
| `family_2_2_6_split_certificate.log` | produced here by `pipeline/magma/family_2_2_6_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,2,6] (found 2026-10-19, torsion_jac/notes/g33rat_2026_10_19.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 30a2 with the partners that have a rational 2-torsion point, over a rational curve; the residue class and the J(𝔽_p) bound are in the log |
| `family_2_12_split_certificate.log` | produced here by `pipeline/magma/family_2_12_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,12] (found 2026-10-19, torsion_jac/notes/g33rat_2026_10_19.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 90c3 with the partners that have a rational 2-torsion point, over a rational curve; the residue class and the J(𝔽_p) bound are in the log |
| `family_24_split_certificate.log` | produced here by `pipeline/magma/family_24_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [24] (found 2026-10-19, torsion_jac/notes/hlpEC_2026_10_19.md): the HLP (2,2)-gluing of 21a3 with E₆(s) on a rational conic (HLP §3.5, N = 8, N′ = 6); the residue class and the J(𝔽_p) bound are in the log |
| `family_5_5_split_certificate.log` | produced here by `pipeline/magma/family_5_5_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [5,5] (found 2026-10-19, torsion_jac/notes/isopair_2026_10_19.md): the HLP (2,2)-gluing of E₅(t⁵) with E/⟨P⟩ along the isogeny (HLP's Corollary 5 curve at t = 1); the residue class and the J(𝔽_p) bound are in the log |
| `family_3_9_split_certificate.log` | produced here by `pipeline/magma/family_3_9_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [3,9] (found 2026-10-19, torsion_jac/notes/isopair_2026_10_19.md): the HLP (2,2)-gluing of E₉(t) with E/⟨3P⟩ along the isogeny; the residue class and the J(𝔽_p) bound are in the log |
| `family_27_split_certificate.log` | produced here by `pipeline/magma/family_27_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [27] (graded open before) (found 2026-10-19, torsion_jac/notes/g33rat_27_2026_10_19.md): (3,3)-gluings of 54b3 on the genus-1 locus where a class of order 27 appears (a 3-power gain; an étale triple cover of z³ = 4g(w), Jacobian rank 1; coset argument); the residue class and the J(𝔽_p) bound are in the log |
| `family_30_split_certificate.log` | produced here by `pipeline/magma/family_30_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [30] (found 2026-10-19, torsion_jac/notes/hlpP2_2026_10_19.md): the HLP (2,2)-gluing of 66c1 with E₆(u) on a line of HLP's conic (§3.5); the residue class and the J(𝔽_p) bound are in the log |
| `family_2_2_2_4_split_certificate.log` | produced here by `pipeline/magma/family_2_2_2_4_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,2,2,4] (found 2026-10-19, torsion_jac/notes/g33ell_2026_10_19.md): bielliptic curves with six rational Weierstrass points and the quotient 21a1 fixed, over a rational curve; the residue class and the J(𝔽_p) bound are in the log |
| `family_2_2_8_split_certificate.log` | produced here by `pipeline/magma/family_2_2_8_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,2,8] (found 2026-10-19, torsion_jac/notes/g33rat_2026_10_19.md): (3,3)-gluings (Bröker–Howe–Lauter–Stevenhagen normal form) of 210e2 with the partners that have a rational 2-torsion point, over a rational curve; the residue class and the J(𝔽_p) bound are in the log |
| `family_2_4_4_split_certificate.log` | produced here by `pipeline/magma/family_2_4_4_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,4,4] (found 2026-10-19, torsion_jac/notes/isopair_2026_10_19.md): an HLP (2,2)-self-gluing of E_{2,4} on the conic where 1 − 16t² is a square; the residue class and the J(𝔽_p) bound are in the log |
| `family_2_16_split_certificate.log` | produced here by `pipeline/magma/family_2_16_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,16] (found 2026-10-19, torsion_jac/notes/pow2_2026_10_19.md): the HLP (2,2)-gluing of 210e2 with a one-parameter family of curves with full 2-torsion, with a gain class of order 16 descended to ℚ(s); the residue class and the J(𝔽_p) bound are in the log |
| `family_4_8_split_certificate.log` | produced here by `pipeline/magma/family_4_8_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [4,8] (found 2026-10-19, torsion_jac/notes/pow2_2026_10_19.md): the HLP (2,2)-gluing of E_{2,4}(u/4) with 210e2 sending 2P₁ to a non-special 2-torsion point; the pullbacks give [4,8]; the residue class and the J(𝔽_p) bound are in the log |
| `family_35_split_certificate.log` | produced here by `pipeline/magma/family_35_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [35] (found 2026-10-19, torsion_jac/notes/hlpEC_2026_10_19.md): the HLP (2,2)-gluing of 26b1 with E₅(t) over HLP's curve X′ (rank 2); the residue class and the J(𝔽_p) bound are in the log |
| `family_3_12_split_certificate.log` | produced here by `pipeline/magma/family_3_12_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [3,12] (found 2026-10-19, torsion_jac/notes/hlpP2_2026_10_19.md): the HLP (2,2)-gluing of 90c3 with E₆(u) on a line of HLP's conic (§3.5); the residue class and the J(𝔽_p) bound are in the log |
| `family_6_6_split_certificate.log` | produced here by `pipeline/magma/family_6_6_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [6,6] (found 2026-10-19, torsion_jac/notes/isopair_2026_10_19.md): the HLP (2,2)-gluing of a Hesse curve with ℤ/6 with E/μ₃ along the isogeny; the residue class and the J(𝔽_p) bound are in the log |
| `family_36_split_certificate.log` | produced here by `pipeline/magma/family_36_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [36] (graded open before) (found 2026-10-19, torsion_jac/notes/g33rat_27_2026_10_19.md): (3,3)-gluings of 90c3 on the genus-1 locus where the partner has a rational point of order 9 (an étale triple cover of z³ = 25g(w), Jacobian rank 1; coset argument); the residue class and the J(𝔽_p) bound are in the log |
| `family_40_split_certificate.log` | produced here by `pipeline/magma/family_40_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [40] (found 2026-10-19, torsion_jac/notes/hlpES_2026_10_19.md): the HLP (2,2)-gluing of E₈(t(s)) with 66c2 at HLP's point (§3.5); the residue class and the J(𝔽_p) bound are in the log |
| `family_45_split_certificate.log` | produced here by `pipeline/magma/family_45_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [45] (found 2026-10-19, torsion_jac/notes/hlpEC_2026_10_19.md): the HLP (2,2)-gluing of 15810u1 with E₅(t) over HLP's curve X′ (rank 2); the residue class and the J(𝔽_p) bound are in the log |
| `family_2_2_12_split_certificate.log` | produced here by `pipeline/magma/family_2_2_12_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,2,12] (found 2026-10-19, torsion_jac/notes/g33ell_2026_10_19.md): (3,3)-gluings of 90c3 with the partners that have full rational 2-torsion, over a genus-1 curve whose Jacobian y² = x³ + 9 has rank 1; the residue class and the J(𝔽_p) bound are in the log |
| `family_2_24_split_certificate.log` | produced here by `pipeline/magma/family_2_24_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,24] (found 2026-10-19, torsion_jac/notes/hlpP2_2026_10_19.md): the HLP (2,2)-gluing of 210e2 with E_{2,6}(t), ψ(4P₈) = 3P₆; the residue class and the J(𝔽_p) bound are in the log |
| `family_5_10_split_certificate.log` | produced here by `pipeline/magma/family_5_10_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [5,10] (found 2026-10-19, torsion_jac/notes/hlpES_2026_10_19.md): the HLP (2,2)-gluing of E₁₀(t) with E₁₀(u), u on the section 2σ of HLP's elliptic surface (§3.5); the residue class and the J(𝔽_p) bound are in the log |
| `family_2_30_split_certificate.log` | produced here by `pipeline/magma/family_2_30_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,30] (graded open before) (found 2026-10-19, torsion_jac/notes/hlpP2_2026_10_19.md): the HLP Prop. 12 gain locus of the [30] family, a genus-1 curve with Jacobian 69696gl1 (rank 2); the residue class and the J(𝔽_p) bound are in the log |
| `family_60_split_certificate.log` | produced here by `pipeline/magma/family_60_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [60] (found 2026-10-19, torsion_jac/notes/hlpEC_2026_10_19.md): the HLP (2,2)-gluing of 90c3 with E₁₀(t) over 900e1 (rank 1); the residue class and the J(𝔽_p) bound are in the log |
| `family_2_2_2_8_split_certificate.log` | produced here by `pipeline/magma/family_2_2_2_8_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,2,2,8] (found 2026-10-19, torsion_jac/notes/g33ell_2026_10_19.md): bielliptic curves with six rational Weierstrass points and the quotient 210e2 fixed, over a rational curve; the residue class and the J(𝔽_p) bound are in the log |
| `family_2_2_4_4_split_certificate.log` | produced here by `pipeline/magma/family_2_2_4_4_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,2,4,4] (found 2026-10-19, torsion_jac/notes/hlpES_2026_10_19.md): the HLP (2,2)-gluing of E_{2,4}(t) with E_{2,4}(u) on the full-gain locus (section 2R), with a gain class of order 4; the residue class and the J(𝔽_p) bound are in the log |
| `family_2_2_16_split_certificate.log` | produced here by `pipeline/magma/family_2_2_16_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,2,16] (found 2026-10-19, torsion_jac/notes/pow2_2026_10_19.md): the HLP (2,2)-gluing of E_{2,4}(t) with 210e2 on a rational conic (the [2,2,2,16] construction without its third condition), gain class of order 16; the residue class and the J(𝔽_p) bound are in the log |
| `family_2_4_8_split_certificate.log` | produced here by `pipeline/magma/family_2_4_8_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,4,8] (found 2026-10-19, torsion_jac/notes/hlpP2_2026_10_19.md): the HLP (2,2)-gluing of 210e2 with E_{2,8}(t), matching the points 4P₈; the residue class and the J(𝔽_p) bound are in the log |
| `family_2_6_6_split_certificate.log` | produced here by `pipeline/magma/family_2_6_6_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,6,6] (found 2026-10-19, torsion_jac/notes/hlpES_2026_10_19.md): the HLP (2,2)-gluing of E_{2,6}(t) with E_{2,6}(u) on the section 2σ of HLP's system (6) (§3.7); the residue class and the J(𝔽_p) bound are in the log |
| `family_3_24_split_certificate.log` | produced here by `pipeline/magma/family_3_24_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [3,24] (graded open before) (found 2026-10-19, torsion_jac/notes/hlpP2_2026_10_19.md): the HLP Prop. 12 gain locus of the [3,12] family, a genus-1 curve with Jacobian 4800cb1 (rank 1); the residue class and the J(𝔽_p) bound are in the log |
| `family_6_12_split_certificate.log` | produced here by `pipeline/magma/family_6_12_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [6,12] (found 2026-10-19, torsion_jac/notes/hlpES_2026_10_19.md): the HLP (2,2)-gluing of E₁₂(t) with E₁₂(u), u on the section 2σ of HLP's elliptic surface (§3.5); the residue class and the J(𝔽_p) bound are in the log |
| `family_2_2_24_split_certificate.log` | produced here by `pipeline/magma/family_2_2_24_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,2,24] (found 2026-10-19, torsion_jac/notes/hlpEC_2026_10_19.md): the HLP (2,2)-gluing of E_{2,6}(241/81) with E_{2,8}(u) over HLP's gain curve (Jacobian y² = x³ − 1681x, rank 2); the residue class and the J(𝔽_p) bound are in the log |
| `family_2_48_split_certificate.log` | produced here by `pipeline/magma/family_2_48_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,48] (graded open before) (found 2026-10-19, torsion_jac/notes/hlpP2_2026_10_19.md): the HLP Prop. 12 gain locus of the [2,24] family, a genus-1 curve with Jacobian 4800c1 (rank 1); the residue class and the J(𝔽_p) bound are in the log |
| `family_2_2_4_8_split_certificate.log` | produced here by `pipeline/magma/family_2_2_4_8_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [2,2,4,8] (found 2026-10-19, torsion_jac/notes/hlpES_2026_10_19.md): the HLP (2,2)-gluing of E_{2,8}(t) with E_{2,8}(u) on the section 2σ of HLP's surface (5) (§3.7); the residue class and the J(𝔽_p) bound are in the log |
| `family_48_split_certificate.log` | produced here by `pipeline/magma/family_48_split_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of split [48] (found 2026-10-19, torsion_jac/notes/hlpEC_48_2026_10_19.md): Howe 2015, Theorems 3.1 and 3.2 — the HLP (2,2)-gluing of 42a1 with E_s (torsion ℤ/6) over Howe's rank-2 curve D, with the gain of HLP Prop. 12; J(𝔽₅) = ℤ/48 on the class |
| `family_3_9_certificate.log` | produced here by `pipeline/magma/family_3_9_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of simple [3,9] (found 2026-10-19, torsion_jac/notes/simple39_2026_10_19.md): a pencil on the [9]-chart with a class of order 9, base-changed so that a second class of order 3 is rational; s ≡ 5 (mod 119), #J(𝔽₇) = 81, #J(𝔽₁₇) = 351, strict at 17 |
| `family_4_8_certificate.log` | produced here by `pipeline/magma/family_4_8_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of simple [4,8] (found 2026-10-19, torsion_jac/notes/simple248_2026_10_19.md): the [4,8] tangent cover restricted to a rational curve on which the half of [x, 0] is rational; τ ≡ 4 (mod 533), #J(𝔽₁₃) = 128, #J(𝔽₄₁) = 1888, strict at 41 |
| `family_2_4_8_certificate.log` | produced here by `pipeline/magma/family_2_4_8_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of simple [2,4,8] (found 2026-10-19, torsion_jac/notes/simple248_2026_10_19.md): the same rational curve of the tangent cover, over the genus-1 curve D: σ² = τ⁴ − 36τ² + 36 (Jacobian y² = x³ − 36x, rank 1) on which the quadratic factor B splits; J(𝔽₁₁) ≅ [2,4,16], J(𝔽₅₃) ≅ [2,4,328], strict at 53 |
| `family_28_certificate.log` | produced here by `pipeline/magma/family_28_certificate.m` (uses `pipeline/magma/family_cert_lib.m`) | — | exact infinitude of simple [28] (2026-10-19, torsion_jac/notes/pf28_2026_10_19.md): the family (8) of Platonov–Fedorov 2018 (Dokl. Math. 98, 468–471), transcribed and checked against the paper's relations (2), (3), its discriminant and its ten printed examples; D = [x − 1, −t⁵(t + 1)] has order 28 over ℚ(t); t ≡ 1 (mod 35), #J(𝔽₅) = 28, #J(𝔽₇) = 56, strict at 7 |
| `split_11_remark.log` | produced here by `pipeline/magma/split_11_remark.m` | — | remark for the open grade of split [11] (2026-10-20, torsion_jac/notes/split11_2026_10_19.md), not a certificate: the (2,2)-gluing of an elliptic curve over a quadratic field with a point of order 11 to its conjugate (the construction of 11.b, BNSS2026) needs a Galois-isomorphism E[2] → E^σ[2]; the parameters s on the X₁(11) model for which it exists lie under a curve Y with a degree-6 map to the s-line, whose reductions modulo 101, 103, 107 are absolutely irreducible of genus 8, so Y(ℚ) is finite; s = 4/5 (11.b) is the only such s to height 400, and a check without Y (does f^σ have a root in K[x]/(f)?) agrees to height 30 |
| `leprevost_3_9_split.log` | produced here by `pipeline/magma/leprevost_3_9_split.m` | — | every member of Leprévost's [3,9] family C_t (Leprévost 1995, pp. 302–304) is geometrically split (2026-10-20): over ℚ(√−3)(t), f_t = q·q′ and u = (x² − x)/q, v = y·hh/q² is a map of degree 3 onto c·v² = E(u), E(u) = Res_y(q′(y), u·q(y) − (y² − y)), of genus 1; it specialises at every rational t ≠ 0, 1. Checked at t = 3, 5, −1/2 against the L-polynomials modulo the primes p ≡ 1 (mod 3) below 200. Leprévost leaves the absolute simplicity of J(C_t) open (p. 303–304); Platonov–Zhgun–Petrunin 2013 (Theorem 5) settle t = 3 |
| `leprevost_notes_check.log` | produced here by `pipeline/magma/leprevost_notes_check.m` | — | the seed curves credited to Leprévost's C. R. notes (1991, 1993), compared with the full versions (equations read from the scans of Leprévost 1995, p. 286, and Manuscripta Math. 75 (1992), pp. 316–322): C₂₂, C₂₄,₁, C₂₄,₂, C₂₆ have points of the claimed orders and are not isomorphic over ℚ̄ to 22.a, 24.a, 26.a; the 1992 families of orders 13, 15, 17, 19 have points of the claimed orders, and a family search (validated on a moved member) finds 13.a at t = ±1, 17.a at t = −1, and no fiber for 15.a and 19.a, even over ℚ̄ |
| `attribution_check.log` | produced here by `pipeline/magma/attribution_check.m` | — | the comparisons behind the credit errata of 2026-10-14 below: census curves 27.b, 24.b, 25.b, 3.9.b, 6.6.b, 5.5.a, 21.c, 27.c against the curves printed in Leprévost 1995 (C21, C24,1, C24,2, C25, C27, pp. 291–302 of the JTNB scan), Howe 2015 (Table 2, order 27, arXiv:1407.2654v3), HLP 2000 (X₀(22), Cor. 5) and the LMFDB, by G2-invariants and isomorphism over ℚ. It also shows that 21.c (LMFDB 676.a.5408.1, issue #7) is Leprévost's C21 and 27.c (issue #8) is his C27, which is why both are credited to Leprévost 1995 |
| `seed_credit_check.log` | produced here by `pipeline/magma/seed_credit_check.m` | — | the comparisons behind the credit errata of 2026-10-20 below: census curves against the curves printed in Platonov–Petrunin 2012a (pp. 287–288), 2012b (p. 643), Platonov–Zhgun–Petrunin 2013 (p. 320), Platonov–Petrunin 2015 (p. 221), Elkies 2002 (N = 32), Bruin–Flynn–Testa 2014 (Thm 6, Examples 19–21) and Leprévost 1995 (C₃, p. 303), by isomorphism over ℚ, each transcription checked against the torsion its source claims; the one-parameter families of Flynn 1991 (Result 3), Leprévost 1997 (C₁,ₜ) and Leprévost 1995 (the [3,9] family Cₜ) searched for 9.a, 10.a, 20.a, 3.9.a, each search validated on a moved member; Flynn 1991's order-6 family (Result 1(b)) against 6.a; the rational Weierstrass points of 21.a |
| `seed_production_check.log` | produced here by `pipeline/magma/seed_production_check.m` | — | for each of the 150 seed curves, isomorphism over ℚ with every production LMFDB curve of the same conductor (read from devmirror.lmfdb.xyz on 2026-10-20): 48 seed curves are production curves, the other 102 are not |

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

## Errata to the seed import: curve credits (review of 2026-10-20)

The 130 seed rows left unchecked on 2026-10-14 were checked: the database credits against the production LMFDB
(`seed_production_check.log`), the paper credits against the equations printed in the sources (pages read from the
scans; `seed_credit_check.log`).  A paper is credited with a curve when it prints that curve (up to isomorphism over
ℚ); a member of a printed family that the paper does not single out is credited to the database it was taken from.
`sources` (the group's first realisation) is unchanged except for [2,28].

* **Database credits.** 17 curves credited to the production database (BSSVY 2016) are not in it (no production curve
  has their G2-invariants): 1.a, 3.b, 4.b, 5.b, 6.b, 7.b, 8.b, 9.b, 2.2.b, 2.4.b, 2.6.b, 2.8.b, 3.3.b, 3.6.b, 4.4.b,
  2.2.2.b, 2.2.4.b.  They are curves of the extended database, now credited to Booker–Sutherland.  The other 23
  production credits and all 38 extended-database credits were confirmed.
* **Not the cited paper's curve; now credited to the database.**
  * 14.a, 18.a: Platonov–Petrunin 2012a print one curve of each order (f₁₄, f₁₈, p. 287), reprinted by
    Platonov–Zhgun–Petrunin 2013; neither is these curves → production database.
  * 5.a, 7.a: prime conductors 277 and 461, so not modular curves, which Ogg 1973 treats (a genus-2 Jacobian that is a
    quotient of some J₁(N) has conductor M² or M₁M₂ with M, M₁, M₂ ≥ 11) → production database.
  * 6.a: not a member of Flynn 1991's order-6 family (Result 1(b), p. 435: every member has a rational Weierstrass point
    W and a rational point P with [P − W] of order 6; 6.a has no such pair) → production database.
  * 10.a: not a fiber of Flynn 1991's order-10 family (Result 3 with r = 1, p. 438) → production database.
  * 9.a, 20.a: the fibers t = 2 of Flynn 1991's order-9 family (Result 3 with r = 2) and t = 1 of Leprévost 1997's
    C₁,ₜ (p. 48), which the papers do not single out → production database.
  * 21.a: it has no rational Weierstrass point, while the order-21 curves of Leprévost 1991 all have one (Leprévost
    1995, p. 290, on his reference [2] = Leprévost 1991) → production database.
  * 32.a: not Elkies' printed order-32 curve → extended database.
  * 3.3.a: none of the curves printed by Bruin–Flynn–Testa 2014; it is the fiber (r, s, t) = (5/2, 16/9, 25/16) of
    their Theorem 6, which parametrises all curves with a rational (ℤ/3)² of trivial Weil pairing → extended database.
    (Their Example 19 is printed with (12x³ − 105)², which has torsion [3]; Theorem 6 gives (13x³ − 105)², with [3,3].)
  * 3.9.a: not a fiber of Leprévost 1995's [3,9] family Cₜ → extended database.
* **A different paper printed the curve.**
  * 11.a is f₁₁ of Platonov–Petrunin 2012a (p. 288); its conductor 353 is prime (see 5.a); Flynn 1990 and 1991 contain
    it only as an unmarked member of their order-11 families.
  * 36.b is f₃₆ of Platonov–Zhgun–Petrunin 2013 (p. 320), not the f₃₆ of Platonov–Petrunin 2012b.
  * 48.a is f₄₈,₁ of Platonov–Petrunin 2012b (p. 643, reprinted by Platonov–Zhgun–Petrunin 2013); it is none of
    Howe 2015's order-48 curves.
  * 19.a: Platonov–Zhgun–Petrunin 2013 print no curve of order 19, so the credit is Leprévost 1991 alone (not seen).
  * 2.28.a is f₂₈,₂ of Platonov–Petrunin 2015 (p. 221): they prove a point of order 28 and absolute simplicity, and its
    three rational linear factors give J(ℚ)[2] = (ℤ/2)², so J(ℚ)_tors = [2,28].  Both `discovery.by` and `sources`
    are now Platonov–Petrunin 2015, which predates the extended database (the paper's Table 1 credits Booker–Sutherland).
  * 3.9.c (submitted 2026-10-16 as an extended-database curve) is Leprévost 1995's C₃ (p. 303, the member t = 3 of his
    [3,9] family, with the printed discriminant −2¹⁸3¹²7), reprinted by Platonov–Zhgun–Petrunin 2013 (p. 321):
    `discovery` is now Leprévost 1995 (1995), and so are its `sources`.
* **Confirmed** (the credited paper prints the curve): 63.a and 7.7.a (HLP 2000, §3.6 and p. 16 of the preprint),
  13.a (Flynn 1990, App. 3.2, read in the author's version), 19.b (X₁(13), Mazur–Tate 1973), 23.a and 29.a (C₂₃ and C₂₉ of Leprévost 1995, the full
  version of the 1993 note), 27.a and 70.a (Howe 2015), 28.a (Platonov–Petrunin 2012a, Thm 4), 28.b
  (Platonov–Zhgun–Petrunin 2013, f₂₈,₂; Howe 2015's C₂₈,₁ reprints it), 30.a and 36.a (Platonov–Petrunin 2015,
  f₃₀,₃ and f₃₆), 33.a (Platonov–Petrunin 2012b), 34.a, 39.a, 40.a (Elkies 2002), 25.a (Nicholls 2018, Example 3.1.7),
  2.2.2.10.a (Elkies 2024, Example 2.12), 31.a (Costa et al., newform 1830.2.a.q).  The 8 curves credited to BNSS
  2026 are in neither database, as the paper's "new" says.
* **Not checked (UNVERIFIED), unchanged** (but see the second review below). The source could not be obtained: Ogg 1973 (21.b and 2.10.b are X₁(18) and
  X₁(16), whose torsion Leprévost 1995, p. 285, attributes to Ogg) and the C. R. notes Leprévost 1991 (15.a, 17.a,
  19.a, and the order-13 part of 13.a) and Leprévost 1993 (22.a, 24.a, 26.a).  Evidence from the same author's full
  versions (Manuscripta Math. 75 (1992), 303–326, for 1991; Leprévost 1995 for 1993): 15.a and 19.a are not fibers of
  the 1992 families of orders 15 and 19, 17.a is the fiber t = −1 of the order-17 family, and 22.a, 24.a, 26.a are
  none of Leprévost 1995's C₂₂, C₂₄,₁, C₂₄,₂, C₂₆.  17.a is also f₁₇ of Platonov–Petrunin 2012a.
* **Reprints** (the credit is unchanged): Platonov–Petrunin 2012a also print 17.a, 19.b (X₁(13)), 23.a and 29.a.
* Database credits say only that the curve is in that database; earlier appearances in papers not cited by the tables
  were not searched for.

## Second review of 2026-10-20: Leprévost's C. R. notes, Ogg 1973, and the simple [3,9]

The C. R. notes and Ogg 1973 still could not be read: Gallica and the AMS answer with a bot challenge, and MathSciNet
gives no review text from here.  zbMATH Open (api.zbmath.org, read 2026-10-20) has summaries of the notes:

* Leprévost 1991 (Zbl 0758.14017): "We construct here, for ℓ = 15, 17, 19 or 21 a family with one parameter t of curves
  of genus 2 … such that its jacobian has a point of order ℓ rational over ℚ(t)", with the same method in genus g as
  the full version, Manuscripta Math. 75 (1992), 303–326 (Zbl 0790.14021).  Leprévost 1991a (Zbl 0758.14016) gives a
  curve over ℚ(t) with a point of order 13; the 1992 paper says that its order-13 family is this one (p. 316).
* Leprévost 1993 (Zbl 0783.14016): a curve C_ℓ for ℓ = 21, 22, 23, 25, 26, 27 and 29, two curves for ℓ = 24, and the
  one-parameter [3,9] family.  These are the curves of Theorems 1.2.1–1.2.2 of the full version, Leprévost 1995
  ("Ce travail, annoncé pour l'essentiel dans [3]", p. 283).

Curve credits changed accordingly (the group-level `sources` are unchanged); the comparisons are in
`leprevost_notes_check.log`:

* 15.a and 19.a are not fibers of the 1992 families of orders 15 and 19, even over ℚ̄, so they are not curves of the
  1991 note: LMFDB genus 2 database (BSSVY 2016) for 15.a, Booker–Sutherland (extended database) for 19.a.
* 22.a, 24.a and 26.a are none of C₂₂, C₂₄,₁, C₂₄,₂, C₂₆ (Leprévost 1995, p. 286): LMFDB genus 2 database (BSSVY 2016)
  for 22.a and 24.a, Booker–Sutherland (extended database) for 26.a.
* Unchanged, UNVERIFIED: 13.a (Flynn 1990, confirmed, and Leprévost 1991a; 13.a is the fiber at t = 1 and at t = −1 of
  the 1992 order-13 family) and 17.a (Leprévost 1991; the fiber t = −1 of the 1992 order-17 family, and Platonov–Petrunin
  2012a's f₁₇).  Whether the notes single out these fibers can only be seen in the notes.
* Unchanged, UNVERIFIED: 21.b and 2.10.b (Ogg 1973).  They are X₁(18) and X₁(16); Leprévost 1995, pp. 284–285, cites
  Ogg's p. 226 for the torsion ℤ/21 and ℤ/2 × ℤ/10 of their Jacobians.  The paper is MR0337974 and Zbl 0273.14008 (no
  zbMATH review).

The simple [3,9].  The paper's Table 1 credits the geometrically simple [3,9] to Leprévost 1995, whose [3,9] family C_t
has Jacobians that are generically simple over ℚ (shown at t = 3).  Leprévost writes that he could not decide whether
they are absolutely simple (pp. 303–304), and Platonov–Zhgun–Petrunin 2013 (Theorem 5) showed that J(C₃) is not.  No
member is: `leprevost_3_9_split.log` gives a map of degree 3 from every C_t (t ≠ 0, 1) onto an elliptic curve over
ℚ(√−3).  So the `sources` of the seed curve 3.9.a (8436.a.2, an extended-database curve that is not a fiber of C_t) are
now Booker–Sutherland (extended database), and the census credits the first geometrically simple [3,9] to it (earlier
papers were not searched for).  Leprévost 1995 still dates the class "simple over ℚ, split over ℚ̄" of [3,9]
(3.9.c = C₃).

## Errata to the copied record (`census_infinity_record.md`, copied verbatim; found by the audit of 2026-09-18, see `notes/`)

* §2.1: the [2] family `y² = x(x⁵ + x + t)` has generic 2-rank 0 (the quintic is irreducible over ℚ(t)); the fiber t = 2 is exact [2] only because x⁵ + x + 2 has the root −1. The site uses `y² = x(x⁴ + x + t)` instead.
* §3, [15]: "inherited from the cyclic [30] family" gives torsion ⊇ [30] ⊇ [15] only, not exact [15]; graded `family` on the site.
* §2.7: the marked divisor is described by `Q = 0, y = 0` on the model `y² = g`, on which the points over the roots of Q have y = ±x; the description mixes the coordinates of two models (not verified further here).
* §3, [4,8]: the record lists the tangent-cover family as certified, but there is no fiber certificate for it, and its only displayed simplicity witness (1849T⁴+38T²+1) is not strict; graded open on the site since 2026-10-13 (notes review).
* §4, [2,2,2,12]: the discussion of infinitude is superseded by the paper's Theorem 3.3 and the exact-infinitude certificate above.
* §1: "excess torsion is a closed condition" is not the right justification for infinitely many exact fibers; see the about page for the residue-class argument.

Found by the audit of 2026-10-18 (torsion_jac/notes/oldcert_A, _B, _C_2026_10_18.md), which gave every family of the record an explicit residue-class certificate (rows above):

* §3, [5]: the record names no family ("five-torsion families", deferring to the lab inventory). The inventory's normalized chart f = h² − (1+a+b)²x⁵, h = 1 + ax + bx², has f(1) = 0 for all a, b, so [(1,0) − ∞] is a rational class of order 2 and the torsion contains [10]; exact [5] never occurs on it. The site certifies [5] with y² = (1+x+x²)² − tx⁵. (§2.2 says f(1) vanishes "when b = 1"; it vanishes for every b.)
* §3, [6]: the plain M(12) chart has generic torsion exactly [6], as stated; but its class of order 2 is divisible by 2 in J(𝔽_p) in every good residue class tested, so it has no explicit residue-class certificate. The site certifies [6] with the contact-6 family h = t + x + x² + x³, f = h² − (x−1)⁶.
* §3, [2,6]: "h cubic" needs h(0) = 1 (a root at x = 0) for the class of order 2; with a general cubic the chart gives [6].
* §2.6, [2,8]: "partner #4" is not one family. Magma's `RichelotIsogenousSurfaces` orders the partners by the size of the roots, so partner #4 is the kernel {[(−b²,0) − ∞], [(0,0) − (−a²,0)], [(−1,0) − (−t²,0)]} at t = 2, 3 and {[(−a²,0) − ∞], [(0,0) − (−b²,0)], [(−1,0) − (−t²,0)]} at t = 5, 7; both carry [2,8]. The site certifies the first.
* §2.6, [4,4]: "partner #5" is an index, not a definition; the certificate uses the kernel {[(−c²,0) − ∞], [(0,0) − (−b²,0)], [(−a²,0) − (−d²,0)]}, which is partner #5 at s = 5.

Found by the audit of 2026-10-18 of this site's own certificates (torsion_jac/notes/oldcert_D_2026_10_18.md):

* [2]: `family_2_certificate.m` of 2026-09-18 listed exact fibers and the class t ≡ 1 (mod 15), on which gcd(#J(𝔽₃), #J(𝔽₅)) = 4, so exactness on that class needed the Hilbert-irreducibility step of the about page. It now uses the class t ≡ 2 (mod 15), with gcd 2.
* [2,2,2,12]: the certificate cited BNSS Theorem 3.1 and the claim ψ(E) ⊂ S without checking them; step (0) now checks both over the function field of E.

## Split [4,16] infinitude (20 October 2026)

The complete proof is published at [`proofs/split-4-16.html`](../../../proofs/split-4-16.html),
with LaTeX source [`split_416_infinitude.tex`](split_416_infinitude.tex).
The standalone certificate is
[`pipeline/magma/family_416_split_certificate.m`](../../../pipeline/magma/family_416_split_certificate.m);
[`family_416_split_certificate.log`](family_416_split_certificate.log)
records its successful execution in Magma V2.30-1 on 20 October 2026.
The proof and certificate were developed in the local research directory
`random_research/g2_split_416_infinitude_2026_10_20`.

For every rational parameter, the explicit family is smooth, splits over ℚ,
and has exact torsion [4,16]. Exactness follows from Mazur's elliptic torsion
classification and an intersection of torsion Kummer images of order two.
The varying elliptic quotient has nonconstant j-invariant, giving infinitely
many geometric isomorphism classes. The construction uses
Howe–Leprévost–Poonen, Propositions 3–4, with the root correspondence specified
in the proof. The certificate checks the gluing formula, discriminant, Kummer
identities, parameterization, and the explicit rational isomorphism from the
member s = 1 to census example `4.16.a`. The general proof does not depend on
the two sample torsion computations in the log.

## Prime torsion obstruction for splitting over ℚ (20 October 2026)

[`proofs/qsplit-prime-obstruction.html`](../../../proofs/qsplit-prime-obstruction.html)
records the elementary consequence of Mazur's theorem that an abelian surface
split over ℚ has no rational torsion of prime order greater than 7. The argument
uses an elliptic subvariety and its elliptic quotient, both defined over ℚ.
The reference is B. Mazur, *Modular curves and the Eisenstein ideal*,
Publ. Math. IHÉS 47 (1977), 33–186.

This excludes [11], [13], [17], [19], [22], [23], [26], [29], [31], [33], [34],
[39], [2,22] and [2,26] in the **split over ℚ** class. The census cells and
group pages link to the proof. The exclusions are generated from the prime
divisors of the invariant factors by `pipeline/knowledge.py`; they are distinct
from the infinitude records for all geometrically split Jacobians.
The known geometrically split, ℚ-simple examples with torsion [11] and [19]
remain in their own class.

## Geometric [23] and [29] exclusions (20 October 2026)

The complete proof is published at
[`proofs/geometric-23-29-obstruction.html`](../../../proofs/geometric-23-29-obstruction.html),
with standalone LaTeX
[`geometric_23_29_obstruction.tex`](geometric_23_29_obstruction.tex).
It proves that every geometrically split abelian surface A/ℚ satisfies
A(ℚ)[23] = A(ℚ)[29] = 0, without a bound on its geometric splitting degree.
The proof was developed and independently reviewed in
`random_research/g2_gsplit_exclusions_2026_10_20/proof_23.md`.

Two nonisogenous elliptic factors descend to a field of degree at most two.
CM squares are excluded by Clark, *Bounds for torsion on abelian varieties with
integral moduli*, Theorem 5. In the non-CM square case the endomorphism-field
Galois group embeds in PGL₂(ℚ). The small-degree cases use
Derickx–Kamienny–Stein–Stoll, Theorem 1.2 (quartic torsion primes), and Khawaja,
Theorem 1.2 (degree-eight torsion primes). The remaining cases use a Galois-stable
maximal local matrix order, a change of Tate lattice preserving composition
factors, and Clifford theory. A hypothetical rational p-torsion point forces
the real cyclotomic field ℚ(ζ_p)⁺ into the endomorphism field, contradicting its
degree. This is a mathematical proof; it does not depend on a search or a
computed sample of curves.

The new exclusions apply to the geometrically split, ℚ-simple cells for [23]
and [29]. Their earlier split-over-ℚ exclusions retain the elementary Mazur
proof link. Since both split classes are excluded, the combined split
infinitude status is now `impossible`, with a link to the complete geometric proof.

## Curve data corrections (20 October 2026)

Ten curves had no minimal model: 1.c, 4.b, 7.c, 11.b, 2.4.c, 3.3.b, 3.6.c, 3.9.c, 2.2.2.c and 2.2.4.b. Their
conductor job had not finished, and the job wrote its results only at the end. `pipeline/fill_conductors.py`
reran `pipeline/magma/conductor_lib.m`, which now writes the reduced minimal model and factored discriminant
before the conductor (logs `data/logs/<id>.conductor.log`). All 187 curves now have a minimal model.

- For 1.c, 7.c, 11.b, 2.4.c and 3.3.b, Magma's genus-2 `Conductor` stops with "Fibre blowups over number fields
  not yet implemented". The census had said it "did not finish within the time limit".
- For 3.6.c, 3.9.c, 2.2.2.c, 4.b and 2.2.4.b, it did not finish. 3.6.c reached the 16 GB memory limit after
  638 s, and the other four were stopped after about 15 minutes.
- The conductors shown are unchanged. They come from the LMFDB label (1.c, 7.c, 2.4.c, 3.6.c, 3.9.c, 2.2.2.c,
  now with their factorisation) or from the elliptic factors N(E₁)N(E₂) (3.3.b, 4.b, 2.2.4.b). 11.b still has no
  conductor.
- The notes of 3.3.b, 4.b and 2.2.4.b said "(Magma's genus-2 Conductor agrees)", but Magma had not computed a
  conductor for them. The value compared was the LMFDB label's. `verify.py` now compares N(E₁)N(E₂) only with
  Magma's own result, and the three notes are corrected.
- The notes of 3.9.c said the curve "is due to the database cited under Reference". Since 20 October it is
  credited to Leprévost 1995 (C₃), and the note now says so.

## Further geometric exclusions (20 October 2026)

The following complete proofs, with standalone LaTeX sources, add five
exclusions to the geometrically split, ℚ-simple column:

* [31]: [`geometric-31-obstruction.html`](../../../proofs/geometric-31-obstruction.html)
  and [`geometric_31_obstruction.tex`](geometric_31_obstruction.tex).
* [33] and [34]:
  [`geometric-33-34-obstruction.html`](../../../proofs/geometric-33-34-obstruction.html)
  and [`geometric_33_34_obstruction.tex`](geometric_33_34_obstruction.tex).
* [35] and [70]: [`qsimple-35-obstruction.html`](../../../proofs/qsimple-35-obstruction.html)
  and [`qsimple_35_obstruction.tex`](qsimple_35_obstruction.tex).

The first two proofs exclude rational points of orders 31, 33, and 34
on every geometrically split abelian surface over ℚ. Their combined
split infinitude records are therefore `impossible`. The third proof
requires simplicity over ℚ and excludes a rational point of order 35,
so it applies to [35] and [70] only in `gsplit`. Both groups have
split-over-ℚ examples; neither receives a combined `impossible` grade.

The proofs were developed and independently reviewed in
`random_research/g2_gsplit_exclusions_2026_10_20/proof31.md`,
`proof_33.md`, `proof_34.md`, and `proof_35_qsimple.md`.
All allow arbitrary geometric splitting degrees.

For [31], a Galois-stable split quadratic order decomposes the residual
Tate module into two blocks exchanged over a quadratic field. A rational
31-point forces an impossible cyclotomic determinant. For [33] and [34],
the endomorphism-field argument reduces to elliptic isogeny classes over
fields of degree 1, 2, or 4. Factoring an isogeny by the primary part of
its kernel combines the two distinct torsion primes on one elliptic curve.
Derickx–Najman, *Classification of torsion of elliptic curves over quartic
fields*, Theorem 1.1, excludes orders 33 and 34 there.

For [35] and [70], polarization duality strengthens the cyclotomic
consequence of the residual tensor argument to ℚ(ζ_p) ⊆ L.
Applying this at 5 and 7 restricts the endomorphism field to degree
1, 2, or 4. The same coprime-isogeny argument then contradicts the
quartic exclusion of order 35. The CM square cases use Clark, Theorem 5;
[31] also uses the elliptic torsion-prime bounds of
Derickx–Kamienny–Stein–Stoll and Khawaja.

## Further exclusions (21 October 2026)

Seven complete proofs add exclusions in both split classes. Each proof page links its standalone
LaTeX source, and the Magma checks it relies on are listed below with their logs (Magma V2.30-1).

* `qsplit`, [2,2,14] and [2,2,2,10]:
  [`qsplit-two-rank-obstruction.html`](../../../proofs/qsplit-two-rank-obstruction.html). Mazur's
  theorem applied to an elliptic subvariety and its quotient: a rational point of order 7 forces
  rational 2-rank at most 2, and rational 2-rank 4 excludes points of orders 5 and 7.
* [22] and [2,22], every geometrically split surface:
  [`geometric-22-obstruction.html`](../../../proofs/geometric-22-obstruction.html) and
  [`geometric_22_obstruction.tex`](geometric_22_obstruction.tex). A rational 11-point forces
  simplicity over ℚ and excludes geometric squares (a residual tensor argument at 11 for non-CM
  squares, Le Fourn–Najman for quadratic ℚ-curves, Bourdon–Clark for CM squares). Two
  Galois-conjugate, non-isogenous factors over a quadratic field would then carry a point of order 22,
  contradicting Kenku–Momose, Theorem (0.1).
* [26], [2,26] and [39], every geometrically split surface:
  [`geometric-26-39-obstruction.html`](../../../proofs/geometric-26-39-obstruction.html) and
  [`geometric_26_39_obstruction.tex`](geometric_26_39_obstruction.tex). Endomorphism fields of degree
  at most 4 reduce to elliptic torsion over fields of degree at most 4 (Kenku–Momose,
  Derickx–Etropolski–van Hoeij–Morrow–Zureick-Brown, Derickx–Najman); CM squares fail a point count
  over 𝔽₂ or 𝔽₃; residual representations at 13 leave a cyclic sextic endomorphism field, where a
  rational 13-point would give a point of X₀(39)/w₃ ≅ 39a1 (rank 0) outside its four rational
  points. Checks: [`geometric_26_39_x039_w3.m`](geometric_26_39_x039_w3.m)
  ([log](geometric_26_39_x039_w3.log)); supplementary
  [`geometric_26_39_checks.m`](geometric_26_39_checks.m) ([log](geometric_26_39_checks.log)).
* [13], geometrically split and simple over ℚ:
  [`qsimple-13-obstruction.html`](../../../proofs/qsimple-13-obstruction.html) and
  [`qsimple_13_obstruction.tex`](qsimple_13_obstruction.tex). Rational 13-torsion is excluded on every
  ℚ-simple, geometrically split abelian surface except non-CM squares with endomorphism field of
  degree 2 or 3; in the cubic case no surface in the isogeny class is principally polarized, and in
  the quadratic case every principal polarization in the isogeny class is a product polarization (a
  lattice argument using Vukorepa's classification of quadratic isogenies). Checks:
  [`qsimple_13_cm_weil.m`](qsimple_13_cm_weil.m) ([log](qsimple_13_cm_weil.log)),
  [`qsimple_13_x0_39.m`](qsimple_13_x0_39.m) ([log](qsimple_13_x0_39.log)),
  [`qsimple_13_x0_26.m`](qsimple_13_x0_26.m) ([log](qsimple_13_x0_26.log)),
  [`qsimple_13_x1_13.m`](qsimple_13_x1_13.m) ([log](qsimple_13_x1_13.log), including a genus-2 curve
  of rank 0 proved unconditionally) and [`qsimple_13_x0_13.m`](qsimple_13_x0_13.m)
  ([log](qsimple_13_x0_13.log)).
* [30], [60] and [2,30], simple over ℚ only:
  [`qsimple-30-obstruction.html`](../../../proofs/qsimple-30-obstruction.html) and
  [`qsimple_30_obstruction.tex`](qsimple_30_obstruction.tex). Coprime torsion in an elliptic isogeny
  class and Derickx–Najman for non-isogenous factors and endomorphism fields of degree 1, 2 or 4; for
  CM squares, good reduction over a totally ramified extension of ℚ₂ and Tate's theorem exclude a point
  of order 15; for non-CM squares, polarization duality at 5 restricts the endomorphism field.
  Clark's Theorem 5 permits orders 30 and 60, so it does not settle the CM case.
* [63], simple over ℚ only:
  [`qsimple-63-obstruction.html`](../../../proofs/qsimple-63-obstruction.html) and
  [`qsimple_63_obstruction.tex`](qsimple_63_obstruction.tex). Residual representations at 7 force a
  cyclic endomorphism field of order 2, 3 or 6; the sextic case fails at 3, and in the cubic case
  Derickx–Etropolski–van Hoeij–Morrow–Zureick-Brown (Theorem A, §6.2) force the isogeny class of a
  twist of 162b1², on which #A(𝔽₅) = 21. Checks:
  [`qsimple_63_arithmetic.m`](qsimple_63_arithmetic.m) ([log](qsimple_63_arithmetic.log));
  consistency check with census 21.b: [`qsimple_63_census_21b.m`](qsimple_63_census_21b.m)
  ([log](qsimple_63_census_21b.log)). The point (2α−3, 2α−2) printed after Theorem A of
  arXiv:2007.13929v2 is not on 162b1; the script uses (1−2α, 2−2α), which has order 21.
* [5,5], [5,10] and [7,7], simple over ℚ only:
  [`qsimple-rank-two-obstruction.html`](../../../proofs/qsimple-rank-two-obstruction.html) and
  [`qsimple_rank_two_obstruction.tex`](qsimple_rank_two_obstruction.tex). For p ≥ 5, rational
  p-torsion of rank two and polarization duality restrict the residual composition factors to 1, χ
  and χ⁻¹; the Galois action on the geometric endomorphism algebra and complex conjugation then give a
  contradiction. The statement fails for p = 3 (census 3.3.c and 3.9.c).

The combined split infinitude records of [13], [22], [26], [39], [2,22] and [2,26] are now
`impossible`. The other groups above have examples in the class not covered by their proof, so their
combined records are unchanged. [26], [2,26] and [39] rest on the [26]/[39] proof, not on the [13]
proof.

The proofs of [22], [30], [5,5]/[5,10]/[7,7] were developed in
`random_research/g2_gsplit_exclusions_2026_10_20/` (`proof_22.md`, `proof_30_qsimple.md`,
`cm_15_obstruction.md`, `proof_odd_rank_two_qsimple.md`); those of [2,2,14], [2,2,2,10], [26], [39]
and [63] in `random_research/exclusions_2026_10_21/` (`prove_qsplit/`, `prove_gsplit/`), and [13] in
`random_research/gsplit13_2026_10_21/proof_13_gsplit.md`. Each was independently reviewed
(`exclusions_2026_10_21/review_*`), and each published page was checked against its reviewed source.

The same day, the proof pages of 20 October gained the published reference for Clark's Theorem 5
(Clark–Xarles, Canad. J. Math. 60 (2008), Theorem 1.4), Mazur's theorem number (8), one line on
semisimplicity in the [23]/[29] proof, and one line on why E and F are L-isogenous in the [35]/[70]
proof.

## LMFDB labels

The LMFDB labels in `table1.txt`/`table2.txt` are of two kinds: production labels
`cond.class.disc.num` (permalinks on www.lmfdb.org) and extended-database labels
`cond.class.num` (a 2026-08 snapshot of alpha.lmfdb.org, not permanent).  The
pipeline re-resolves production labels by asking www.lmfdb.org to look the curve
up by its equation (`/Genus2Curve/Q/?jump=[[f],[h]]`) and records the label it
redirects to; extended-database rows are linked by equation.
