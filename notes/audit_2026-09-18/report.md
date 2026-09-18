Audit of https://genus-2-torsion.github.io/ — 18 September 2026

I checked the public snapshot generated on 17 September 2026 at 15:52 UTC: 95 groups and 150 curves. The individual torsion groups check out. I found nine incorrect conductors, an incorrect family formula, a missing split example and display/counting errors, omissions in the infinitude information, and incorrect bibliography links. I also obtained an argument that resolves the site's stated open exact-infinitude question for [2,2,2,12].

All computations below used Magma V2.29-4. The scripts, logs, downloaded data, and source repositories are in this directory. The website itself was not modified.

**1. Nine conductor values are incorrect.** For each row below, an explicit isogeny over Q identifies the Jacobian with a product of two elliptic curves up to isogeny. Consequently its conductor is the product of the elliptic conductors. These are exact corrections, not estimates inferred from labels.

| Curve ID | Torsion | Displayed conductor | Correct conductor | Elliptic conductors |
|---|---|---:|---:|---|
| [5.b](https://genus-2-torsion.github.io/curve.html?id=5.b) | [5] | 126877696 | 121 | 11 × 11 |
| [2.4.b](https://genus-2-torsion.github.io/curve.html?id=2.4.b) | [2,4] | 230400 | 225 | 15 × 15 |
| [2.8.b](https://genus-2-torsion.github.io/curve.html?id=2.8.b) | [2,8] | 230400 | 225 | 15 × 15 |
| [4.4.b](https://genus-2-torsion.github.io/curve.html?id=4.4.b) | [4,4] | 230400 | 225 | 15 × 15 |
| [4.8.b](https://genus-2-torsion.github.io/curve.html?id=4.8.b) | [4,8] | 230400 | 225 | 15 × 15 |
| [5.5.a](https://genus-2-torsion.github.io/curve.html?id=5.5.a) | [5,5] | 123904 | 121 | 11 × 11 |
| [7.7.a](https://genus-2-torsion.github.io/curve.html?id=7.7.a) | [7,7] | 56580108189696 | 53958996 | 98826 × 546 |
| [2.2.2.b](https://genus-2-torsion.github.io/curve.html?id=2.2.2.b) | [2,2,2] | 322560 | 315 | 15 × 21 |
| [2.4.4.b](https://genus-2-torsion.github.io/curve.html?id=2.4.4.b) | [2,4,4] | 230400 | 225 | 15 × 15 |

The errors are entirely in the exponent of 2: the displayed conductor is too large by 2^20 for 5.b and 7.7.a, and by 2^10 for the other seven. These records already warn that the computation using Ogg's formula at 2 is not guaranteed, so the uncertainty is disclosed; the explicit elliptic factors now settle these nine cases. Eight also disagree with the conductor prefix of their stored alpha-LMFDB labels.

I compared the conductor products for 55 rational elliptic-product certificates; these nine disagree and the other 46 agree. This does not certify every conductor on the site. In particular, the additional conductor warnings at 2 remain relevant.

Reproduce the nine corrections with `magma -b audit_conductors.m`. The [script](audit_conductors.m) constructs the products afresh using `RichelotIsogenousSurfaces`; the [log](audit_conductors.log) ends with `NINE_CONDUCTOR_CORRECTIONS_VERIFIED`. The complete geometric splitting check is in [audit_split.log](audit_split.log), using the [paper's verifier](https://github.com/AndrewVSutherland/Genus2Torsion/blob/main/verify_split_certificates.m).

**2. The displayed [2] family does not generically have rational 2-torsion.** The [group page](https://genus-2-torsion.github.io/group.html?g=2), `pipeline/knowledge.py`, and §2.1 of the [archived infinitude record](https://genus-2-torsion.github.io/data/knowledge/sources/census_infinity_record.md) use

\[
y^2=x(x^5+x+t).
\]

Over Q(t), the quintic is irreducible. The branch-point orbit sizes are therefore 1 and 5. There is no Galois-stable unordered pair of branch points, so the generic Jacobian has no nonzero rational 2-torsion. Magma also gives trivial rational torsion at t=3,4,5.

The cited t=2 fiber really does have exact torsion [2] and a strict prime 17; it acquires an extra root at x=-1. The error is the claimed generic marked 2-torsion, not that sample. A single special fiber cannot supply the missing generic subgroup.

The likely intended elementary family is `y^2=x(x^4+x+t)`: its odd-degree model has the two rational Weierstrass points x=0 and infinity. I checked exact torsion [2] at t=1,2,3 and distinct G2 invariants, but have not substituted a complete replacement infinitude proof into the site. See [audit_families.log](audit_families.log), [audit_22212.m](audit_22212.m), and [check_22212_progression.log](check_22212_progression.log).

**3. Exact infinitude for geometrically simple [2,2,2,12] can be proved.** The [group page](https://genus-2-torsion.github.io/group.html?g=2.2.2.12) says it is open whether infinitely many curves in the known family have both exact torsion [2,2,2,12] and geometrically simple Jacobians. The following local argument establishes precisely this. This is a deduction from the published family plus new arithmetic checks, rather than a claim that the paper explicitly states the stronger result.

Use the family and elliptic curve in [BNSS, Theorems 3.1 and 3.3](https://arxiv.org/html/2608.28543v1#S3):

\[
E:y^2=x^3-21x-20,\qquad Q=(-3,4),
\]

and the rational map

\[
w=\frac{3x+y-6}{y+18},\quad
 t=-\frac{4x^2-7x+12y+16}{(x-8)(3x+2y+12)},
\]
\[
\psi(x,y)=[3t(1-t):3t-1:2w(3t^2+1):3t^2-3t+2:6t^2-3t+1].
\]

For P=[a:b:c:u:v] in its nondegenerate image, put
\[
C_P:Y^2=X(X-a^2)(X-b^2)(X-c^2)(X-u^2)(X-v^2).
\]
Its rational torsion contains G=[2,2,2,12], of order 96. The point Q has infinite order, and
\[
\psi(2Q)=[120:143:266:218:241]=P_1.
\]

Here are independent local computations for C at P1:

| p | Frobenius polynomial | #J(F_p) | Order of Q mod p |
|---:|---|---:|---:|
| 29 | T^4 + 22T^2 + 841 | 864 | 16 |
| 37 | T^4 − 26T^2 + 1369 | 1344 | 10 |
| 71 | T^4 − 8T^3 + 46T^2 − 568T + 5041 | 4512 | 14 |

The gcd of the first two point counts is 96. At 71 the polynomial is irreducible and a root pi satisfies `[Q(pi^n):Q]=4` for every n=2,...,12, certifying geometric simplicity of the reduction.

Set M=lcm(16,10,14)=560. Every point (2+560k)Q reduces to 2Q at all three primes. The map psi is regular there, and the corresponding genus-2 curves have the same good reductions, up to the coordinate scaling implicit in the projective parameters. Their rational torsion has order dividing 96 by reduction at 29 and 37, and already contains G of order 96. Thus it is exactly G. Their reduction at 71 proves geometric simplicity.

Finally, the moduli map is nonconstant: I checked that the G2 invariants at psi(2Q) and psi(4Q) differ. A nonconstant map from the elliptic parameter curve has finite fibers over its image, so the infinitely many points (2+560k)Q yield infinitely many geometrically nonisomorphic genus-2 curves. Finitely many excluded parameters do not affect this conclusion.

This justifies upgrading the simple-side badge from `family` to `exact`. See [audit_22212.m](audit_22212.m), [its log](audit_22212.log), and the [congruence check](check_22212_progression.m). The archived August 21 record's pessimistic discussion of this group is also superseded by the later published family and this argument.

**4. The [2,2,2,2] split information is incomplete, and its display is inconsistent with its own certificate.** The stored curve `2.2.2.2.b` has an explicit geometric splitting certificate: a Richelot isogeny to a Weil restriction from Q(sqrt(7)). Only whether it splits over Q is marked undecided. Nevertheless, the [group page](https://genus-2-torsion.github.io/group.html?g=2.2.2.2) displays no known geometrically split example and says there is no simplicity or splitness certificate. The JavaScript counts only the `qsplit` and `gsplit` buckets, dropping `split_undecided_over_Q` from that summary. The number of groups with a geometrically split example should be 77, not 76, using the current data alone.

Furthermore, an explicit example that splits over Q and has this exact torsion is
\[
C:y^2=(x^2-1)(x^2-4)(x^2-49).
\]
The six rational branch points supply [2,2,2,2]. Point counts give #J(F_11)=128 and #J(F_17)=400, whose gcd is 16, proving exactness. Its elliptic quotients are
\[
E_1:v^2=(u-1)(u-4)(u-49),\qquad
E_2:v^2=u(u-1)(u-4)(u-49),
\]
with maps (u,v)=(x^2,y) and (x^2,xy). The elliptic conductors are 45 and 630, hence N(J)=28350. This supplies the missing Q-split class. I do not claim that 28350 is minimal, or that the site's existing curve `2.2.2.2.b` itself splits over Q.

See [audit_missing_qsplit.m](audit_missing_qsplit.m), [its log](audit_missing_qsplit.log), and [check_small_examples.log](check_small_examples.log).

**5. Several `open` infinitude badges omit consequences of families already recorded.** The site's `family` grade means a family with torsion containing the stated group. Such a statement automatically descends to subgroups. In particular, these simple-side entries can at least receive that grade:

| Currently open group G | Recorded larger group H containing G |
|---|---|
| [16] | [32] |
| [2,2,6] | [2,2,2,6] |
| [2,2,8] | [2,2,2,8] |
| [2,2,12] | [2,2,2,12] |
| [2,4,4] | [2,2,4,4] |

This does not establish infinitely many curves with exact torsion G. It shows that the site's weaker family information is not propagated consistently. More generally, the JSON defines `open` as no family being recorded, which should not be read as a mathematical nonexistence claim.

On the split side, taking subgroups of the 18 recorded HLP family groups supplies 33 additional currently open entries:

`[]`, `[2]`, `[3]`, `[4]`, `[5]`, `[6]`, `[7]`, `[8]`, `[9]`, `[10]`, `[12]`, `[15]`, `[24]`, `[2,2]`, `[2,4]`, `[2,6]`, `[2,8]`, `[2,12]`, `[3,3]`, `[3,6]`, `[4,4]`, `[4,8]`, `[5,5]`, `[2,2,2]`, `[2,2,4]`, `[2,2,6]`, `[2,2,8]`, `[2,2,12]`, `[2,4,4]`, `[2,2,2,2]`, `[2,2,2,4]`, `[2,2,2,8]`, `[2,2,4,4]`.

There is also a direct literature omission for [48]: [Howe, Theorem 3.2 and Remark 3.3](https://arxiv.org/pdf/1407.2654) gives an infinite family parametrized by an elliptic curve of rank 2. These are Q-split Jacobians, by the 2-gluing construction of Theorem 3.1. Thus at least `family` is justified for split [48] and its subgroups, including [16]. Exact [48] infinitude is a separate statement.

**6. Nine bibliography DOI links are wrong, and a thesis link is broken.** These corrections were checked against publisher, author-repository, or author-supplied arXiv metadata. The reference IDs are those used by the site.

| Reference ID | Current incorrect DOI | Correct DOI and supporting source |
|---|---|---|
| Flynn1990 | 10.1016/0022-314X(90)90036-Q | **10.1016/0022-314X(90)90089-A** — [Oxford](https://ora.ox.ac.uk/objects/uuid%3Ad5489b8b-4744-4bc5-b535-f9daf97b0327) |
| Flynn1991 | 10.1007/BF01243920 | **10.1007/BF01243919** — [Oxford](https://ora.ox.ac.uk/objects/uuid%3A0231ccf2-3145-459c-ae47-7230eaf931d7) |
| Leprevost1995 | 10.5802/jtnb.145 | **10.5802/jtnb.144** — [journal](https://jtnb.centre-mersenne.org/item/JTNB_1995__7_1_283_0/) |
| Howe2015 | 10.1112/blms/bdu110 | **10.1112/blms/bdu107** — [arXiv metadata](https://arxiv.org/abs/1407.2654) |
| PP2012 | 10.1134/S1064562412020366 | **10.1134/S1064562412020330** — [publisher](https://www.pleiades.online/cgi-perl/search.pl?name=danmath&number=2&page=286&type=abstract&year=12) |
| PP2012b | 10.1134/S1064562412050146 | **10.1134/S1064562412050304** — [publisher](https://link.springer.com/article/10.1134/S1064562412050304) |
| PP2015 | 10.1134/S1064562415020222 | **10.1134/S1064562415020325** — [publisher](https://link.springer.com/article/10.1134/S1064562415020325) |
| DaowsudSchmidt | 10.1016/j.jnt.2017.11.001 | **10.1016/j.jnt.2017.11.014** — [publisher](https://www.sciencedirect.com/science/article/pii/S0022314X17304602) |
| Zywina2022 | 10.1007/s40993-022-00389-x | **10.1007/s40993-022-00391-0** — [publisher](https://link.springer.com/article/10.1007/s40993-022-00391-0) |

Some current links identify unrelated articles: `jtnb.145` is Cassou-Noguès–Taylor, the current Daowsud–Schmidt DOI is an article by Hou, and the current PP2012b DOI is an article by Denisov.

The `Nicholls2018` Oxford URL contains an incorrect UUID and returns 404. The [correct thesis record](https://ora.ox.ac.uk/objects/uuid%3A04cef70a-2ab9-44c2-8bbe-ca2ac33bfe41) has DOI **10.5287/ora-z58j05keq**.

**7. Some mathematical explanations need correction or additional justification.** These are distinct from errors in the 150 recorded torsion groups.

- The [About page](https://genus-2-torsion.github.io/about.html) says the relevant root-of-unity order divides 12. That divisibility assertion is false; orders 5, 8, or 10 need not divide 12. For example, the quartic CM field Q(zeta_5) contains fifth roots of unity. The implemented test checks every exponent from 2 through 12, so this wording mistake does not by itself invalidate that test.
- The About page and archived record describe excess rational torsion as a closed condition. Rational torsion can jump on a Zariski-dense thin set of rational parameters. Generic exactness can instead be justified using specialization, provided the generic marked subgroup has actually been constructed. Passing from generic exactness to infinitely many rational fibers with exact torsion needs an appropriate argument, for example local bounds together with Hilbert irreducibility for a rational base, or the explicit congruence argument above for the elliptic base. Sampling finitely many fibers alone is insufficient.
- The archived record's §2.7 writes `y^2=g`, with `g=4Q^3+(x^2-6x+1)Q^2+2x(x-1)Q+x^2`, but describes the marked divisor using `Q=0, y=0`. On this displayed model, reduction modulo Q gives `y^2=x^2`; the points have y=±x, not generically y=0. This appears to mix coordinates from the even and mixed models.
- The exact [15] entry is described as inherited from the cyclic [30] family. Containment proves only the weaker assertion. It cannot by itself produce a curve with exact torsion [15]. The record refers to an exactness discussion in a private inventory, so I classify this as a publicly unverified justification, not as a disproof of exact [15] infinitude.
- The About page's unrestricted statement that no group has been ruled out needs the admissibility restrictions from the cited Question 3. Literally, groups such as [3,3,3] are excluded: the rational 3-torsion must be isotropic for the Weil pairing, and has dimension at most 2. The question concerns groups remaining after such elementary restrictions.

**What passed, and what remains outside this audit.** I recomputed all 150 exact torsion groups with `TorsionSubgroup`, rational Weierstrass-point counts, discriminants of `4f+h^2`, G2 invariants, and all recorded strict or Q-simplicity Frobenius polynomials using naive finite-field point counting. These checks passed. I reran all 77 explicit geometric splitting certificates; all passed. The 55 rational involution cases were also checked independently, and the other 13 Q-split cases have verified rational elliptic covers. The eight geometrically split but Q-simple records have verified irreducible Frobenius certificates. I did not resolve the existing split-over-Q-undecided case.

I also recomputed every stored finite-field point-count gcd and checked that all 146 stored minimal-model equations are Q-isomorphic to the submitted curves, with the stated discriminants. This checks the equations and discriminants, not a global proof of their minimality. The site build consistency check and its 11 Python tests passed; those tests do not catch the mathematical and bibliographic issues above.

I have not certified all conductor exponents, minimum-conductor record claims, every alpha-LMFDB link, or all 40 exact-infinitude claims. Some supporting proofs/scripts are referenced only in a private repository. The results support the correctness of the individual torsion and geometric splitting/simplicity data checked here, rather than a blanket certification of every field on the site.

The principal reproducible outputs are [audit_curves.log](audit_curves.log), [audit_metadata.log](audit_metadata.log), [audit_split.log](audit_split.log), [audit_involutions.log](audit_involutions.log), and [audit_conductors.log](audit_conductors.log), with the corresponding `.m` files in this directory. The split verifier is in `paper/verify_split_certificates.m` and must be run from the `paper` directory so it can read `table2.txt`.
