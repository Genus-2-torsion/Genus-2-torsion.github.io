# 2026-09-18 — response to the audit of the site

An independent audit of the public snapshot of 17 September (report and Magma scripts in
`notes/audit_2026-09-18/`) raised seven points. Each was re-verified here before acting; all were
found correct. What changed:

| # | finding | verified by | change |
|---|---|---|---|
| 1 | nine conductors of ℚ-split curves too large by 2¹⁰ or 2²⁰ (Magma's Ogg-formula exponent at 2) | own Richelot products + elliptic conductors (`scratch/cond9.m`; the products match the LMFDB labels 121, 225, 315) | `verify_lib.m` records the product E₁ × E₂ over ℚ found by the Richelot search (also after an involution certificate); `verify.py` uses N(E₁)N(E₂) as the conductor (rigorous: isogeny invariant), keeps Magma's value in `conductor_magma`, and flags the disagreement. Seed re-run: 56 curves get a product conductor, exactly the audit's 9 change (5.b, 5.5.a → 121; 2.4.b, 2.8.b, 4.4.b, 4.8.b, 2.4.4.b → 225; 2.2.2.b → 315; 7.7.a → 53958996). |
| 2 | the record's [2] family y² = x(x⁵+x+t) has no rational 2-torsion generically (x⁵+x+t irreducible over ℚ(t)); its fiber t = 2 is exact only by an extra root | Magma: irreducibility over ℚ(t), fibers t = 3, 4 have trivial torsion | replaced by y² = x(x⁴+x+t) with the marked class [(0,0) − ∞] checked over ℚ(t), five exact strict fibers with distinct G2-invariants (`pipeline/magma/family_2_certificate.m`, log in `data/knowledge/sources/`) |
| 3 | exact infinitude for geometrically simple [2,2,2,12] follows from the paper's family by a congruence argument on the rank-1 elliptic base | every ingredient recomputed: good reduction of C_{P₁} at 29, 37, 71; orders 16, 10, 14 of Q mod p; #J(𝔽₂₉) = 864, #J(𝔽₃₇) = 1344 (gcd 96); strict Frobenius at 71; 2Q p-integral, ψ's denominators and ψ(2Q)'s coordinates p-adic units (this is what makes ψ((2+560k)Q) ≡ ψ(2Q) mod p for every k); G2 at ψ(2Q) ≠ G2 at ψ(4Q); explicit check at k = 1 | grade `family` → `exact` (`pipeline/magma/family_22212_certificate.m`, log in `data/knowledge/sources/`) |
| 4 | [2,2,2,2]: the geometrically split but ℚ-undecided curve 784.c.1 was dropped from the "geometrically split" counts and mis-described; no ℚ-split example | code review; the audit's curve y² = (x²−1)(x²−4)(x²−49) verified through the pipeline: torsion [2,2,2,2], involution x → −x, elliptic quotients of conductors 45 and 630, N = 28350 | group page and `n_groups_any_split` count the undecided-over-ℚ curves as geometrically split (77 groups), with correct wording; the curve is in the census as 2.2.2.2.c |
| 5 | `family` grades not propagated to subgroups; Howe 2015's order-48 family missing on the split side | abelian-group embedding criterion implemented (`knowledge.is_subgroup`); Howe's abstract confirms the rank-2 elliptic family | inherited `family` grades (lighter badge, "inherited from [H]"): simple side [16], [2,2,6], [2,2,8], [2,2,12], [2,4,4]; split side 34 groups; split [48] → `family` (Howe2015) |
| 6 | nine wrong DOIs, one broken thesis URL | Crossref metadata for every DOI in the bibliography; doi.org for the thesis | all corrected; the remaining DOIs were confirmed |
| 7 | wording: roots of unity "divide 12" (false: orders 5, 8, 10 occur; the test covers all n ≤ 12); "excess torsion is a closed condition" is not a valid justification; [15] rests on containment in [30]; "no group ruled out" needs the Weil-pairing restrictions; the record's §2.7 mixes coordinates of two models | reasoning; the lab inventory row for [15] indeed says only "subgroup of the cyclic [30] family" | about page rewritten: orders ≤ 12; the passage from generic exactness to infinitely many exact fibers now uses a residue class of a good pair of primes (uniform torsion bound, propagated strict prime) plus Hilbert irreducibility (thin sets miss infinitely many integers of a residue class); [15] downgraded to `family`; Question 3's restrictions stated. The §2.7 remark concerns the archived record (copied verbatim) and is noted here, not edited. |

Also from the audit's "what passed": all 150 torsion groups, all strict and ℚ-simplicity Frobenius
polynomials, all 77 geometric splitting certificates and the 55 involution certificates were
reproduced independently.

Still open after this response: the ℚ-status of 784.c.1 (a Weil restriction of a CM curve); all
other exact-infinitude certificates still rest on the record's private scripts; conductors of
ℚ-split curves without a Richelot product (13 cover-certified rows) and of non-split curves with
v₂(Δ) ≥ 12 keep Magma's Ogg-formula exponent at 2, flagged on their pages.
