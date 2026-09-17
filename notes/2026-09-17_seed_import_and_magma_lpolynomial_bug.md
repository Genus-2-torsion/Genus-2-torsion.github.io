# 2026-09-17 — seed import of the paper's tables, and a Magma L-polynomial bug

## Seed import

`pipeline/import_paper.py` queued the 150 rows of `table1.txt` / `table2.txt` (Genus2Torsion
repository, commit 508852a); `pipeline/verify.py --require-claim --no-github --no-lmfdb
--conductor-timeout 150` verified all of them (log: `data/knowledge/seed_import.log`):

| outcome | rows |
|---|---|
| torsion group equal to the table's group | 150 / 150 |
| geometrically simple (strict prime), Table 1 | 73 / 73 — same strict prime and χ as `verify_simple_certificates.m` in all 73 cases |
| split over ℚ (Table 2) | 68: involution over ℚ, Richelot (2,2)-isogeny over ℚ to a product, or one of the 13 covers to a Cremona curve |
| simple over ℚ, split over ℚ̄ (Table 2) | 8: [3], [11], [19], [21], [25], [27], [2,2], [2,10] — χ irreducible at some good prime, split by a Weil restriction, a geometric involution, or a cover over a number field |
| geometrically split, undecided over ℚ | 1: [2,2,2,2], 784.c.1 — (2,2)-isogenous over ℚ to Res_{ℚ(√7)/ℚ} E with E a twist by 3√7 of the CM curve j = 255³ (CM by ℤ[√−7]); χ_p is reducible at every good prime because of the CM, so ℚ-simplicity cannot be certified by an irreducible χ_p, and E is not a base change, so ℚ-splitness is not certified either |
| conductor not computed within 150 s (taken from the LMFDB label) | 4: paper-t2-004, -011, -046, -062 (split Jacobians whose curve has bad reduction where the Jacobian does not; Magma's regular-model computation is slow there) |

## Magma V2.29-4: the default L-polynomial algorithm is wrong for some genus-2 curves over 𝔽₃, 𝔽₅, 𝔽₇

Found because the first run rejected row t2-007 ([7], 676.b.3) with "#J(ℚ)tors = 7 does not divide
gcd #J(𝔽p) = 3".  For the working model

    y^2 = -11x^6 - 16x^5 + 8x^4 + 30x^3 + 8x^2 - 16x - 11    over F_5

`LPolynomial(C)` and `EulerFactor(Jacobian(C))` return `25T^4 + 4T^2 + 1` (i.e. #C(𝔽₅) = 6,
#J(𝔽₅) = 30), while brute-force counting and Magma's own `#Points(C)` give #C(𝔽₅) = 10,
#C(𝔽₂₅) = 36, hence `25T^4 + 20T^3 + 13T^2 + 4T + 1` and #J(𝔽₅) = 63 = 7·9 (consistent with the
rational 7-torsion).  `LPolynomial(C : Al := "Naive")` returns the correct polynomial.

Random test (`scratch`, 600 random squarefree quintics/sextics per prime): default ≠ naive for
1 curve at p = 3, 2 at p = 5, 1 at p = 7, none at p = 11, 13, ...; in every mismatch the naive
result agrees with `#Points`.  Examples: `2t^6 + t^5 + t^3 + 2t + 1` over 𝔽₃ (default
`9T^4 + 3T^3 + 4T^2 + T + 1`, correct `9T^4 + T^2 + 1`); `3t^6 + t^5 + 4t^4 + 3t^3 + 4t + 2` over 𝔽₅
(default `25T^4 + 5T^3 + 7T^2 + T + 1`, correct `25T^4 - 10T^3 + 3T^2 - 2T + 1`).

Consequences: `pipeline/magma/verify_lib.m` computes every χ_p with `Al := "Naive"` (which is also
~8× faster at p ≈ 500) and aborts if the result disagrees with `#Points(C)`, the functional
equation or the Weil bounds.  The paper's `verify_simple_certificates.m` uses the default
`EulerFactor`; its 73 strict-prime witnesses were re-derived here with the naive algorithm and
coincide exactly (same prime, same χ), so the paper's certificates are unaffected — but any
Frobenius-based certificate at p ∈ {3, 5, 7} computed with the default algorithm should be
re-checked, and the bug reported to Magma.
