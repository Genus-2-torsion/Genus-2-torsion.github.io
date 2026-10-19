/* family_2_30_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E1 x E2) with J(Q)_tors exactly [2,30].
   Found 2026-10-19 (torsion_jac notes/hlpP2_2026_10_19.md), as a by-product of the [30] certificate:
   the census curve 2.30.a (Booker-Sutherland; LMFDB 13068.f.2) is the member s = 3/11 of the [30] family.

   Source.  The [30] family of family_30_split_certificate.m: the Howe-Leprevost-Poonen (2,2)-gluing
   (HLP 2000, Prop. 4) of E1 = 66c1 = [1,0,0,-45,81] (torsion Z/10, P10 = (0,9), T1 = 5 P10) and
       E2 = E_6(u) : y^2 + (1 - u)xy - (u^2 + u)y = x^3 - (u^2 + u)x^2,   u = (33 s^2 - 1)/(9 - 33 s^2),
   (P6 = (0,0), T1' = 3 P6) along psi : E1[2] -> E2[2] (alpha_i -> beta_i; alpha_1, beta_1 the rational roots of
   the 2-division cubics, alpha_2 = (-b1 + w1 sqrt 33)/2, beta_2 = (-b2 + w2 sqrt 33)/2 as computed below),
   model C_s : y^2 = f_s(x).
   The gain.  By HLP 2000, Prop. 12, a point (P0, Q0) of E1 x E2 with 2P0 = T1 and 2Q0 = T1' maps to a rational
   point D of J = (E1 x E2)/graph(psi) iff the Kummer images of T1 and T1' agree under psi; with
   L = Q x Q(sqrt 33) this is the single condition
       N(s) = (beta_1 - beta_2)(alpha_1 - alpha_2)  is a square in K = Q(sqrt 33)
   (the Q-components are the norms of the K-components).  Over K(s), N = c (s - a)(s - b) g(s)^2 with a, b, and all
   zeros and poles of g, in K \ Q.  So for rational s the condition is (w1 + w2 sqrt 33)^2 = c (s - a)(s - b) with
   w1, w2 in Q: the intersection X of two quadrics in P^3, a curve of genus 1.  The map
       phi : E' -> X,  E' = 69696gl1 : y^2 = x^3 + 2904x - 74536  (rank 2, P = (385, 7623) of infinite order),
   given by four quartics (embedded below; checked to map E' into X), sends T = (22, 0) to a point over
   s = 3/11.  2D = image of (T1, T1') = 0, since (T1, T1') is in graph(psi); and D is not in the image of
   E1(Q) x E2(Q): otherwise 2 P0 = T1 would lie in 2 E1(Q), but E1(Q)[2^inf] = Z/2.  So J(Q)[2] has rank 2 and
   [2,30] = <D, D30> embeds in J_s(Q)_tors at every smooth member over X(Q).

   Argument.
   (1) Splitting and (2) the marked [30] over Q(s): as in family_30_split_certificate.m (repeated below).
   (3) The gain locus.  N, its odd part c (s - a)(s - b) and the quadrics of X are computed from the root labels
       of step (0); N/(c (s - a)(s - b)) is a square in K(s), and a, b, the zeros and poles of its square root
       are not rational.  phi maps E' into X (the quadrics pull back to multiples of the Weierstrass cubic).
       At every member below N(s) is checked to be a nonzero square in K.
   (4) Exactness.  For every 5-integral s = 3 (mod 5), C_s has the reduction of C_{3/11} at 5, with
       J(F_5) = [2,30].  Reduction is injective on J_s(Q)_tors, so J_s(Q)_tors embeds in [2,30].
   (5) Infinitely many members in the class.  E' has good reduction at 5 and P has order 8 mod 5, so
       R_k = T + 8kP reduces to T mod 5.  The quartics have integral coefficients and the fourth is a 5-unit at T,
       so s(R_k) is 5-integral and = s(T) = 3/11 = 3 (mod 5).  For every k, C_{s(R_k)} therefore has J(Q)_tors
       exactly [2,30].  The R_k are distinct and s has finite degree on E', so infinitely many s occur; the
       G2-invariants at s(T), s(T + 8P), s(T - 8P) differ, so s -> M_2 is nonconstant with finite fibres:
       infinitely many isomorphism classes.
   (6) Sanity: TorsionSubgroup at s = 3/11 (= 2.30.a, in the class) and at s(P) = 17/33 and s((25,117)) = 19/99
       (on the gain locus, outside the class) is [2,30].
   Run (from the repository root):
       magma -b pipeline/magma/family_2_30_split_certificate.m > data/knowledge/sources/family_2_30_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
SetQuitOnError(true);                   // fail closed: a failed assertion stops the run before the marker
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

function MonEC(G)                       // Y^2 = G(X), G a cubic, as a Weierstrass curve (X' = lc X, Y' = lc Y)
    cs := [ Coefficient(G, i) : i in [0..3] ];
    return EllipticCurve([0, cs[3], 0, cs[2]*cs[4], cs[1]*cs[4]^2]), cs[4];
end function;
function MonicCubic(E, x)               // y^2 = MonicCubic(x) after completing the square
    a1, a2, a3, a4, a6 := Explode(aInvariants(E));
    return x^3 + (a2 + a1^2/4)*x^2 + (a4 + a1*a3/2)*x + (a6 + a3^2/4);
end function;
function HLPSextic(al, be)              // HLP 2000, Prop. 4, with psi(al_i) = be_i
    L := Parent(al[1]); PL<X> := PolynomialRing(L);
    aa1 := (al[3]-al[2])^2/(be[3]-be[2]) + (al[2]-al[1])^2/(be[2]-be[1]) + (al[1]-al[3])^2/(be[1]-be[3]);
    bb1 := (be[3]-be[2])^2/(al[3]-al[2]) + (be[2]-be[1])^2/(al[2]-al[1]) + (be[1]-be[3])^2/(al[1]-al[3]);
    aa2 := al[1]*(be[3]-be[2]) + al[2]*(be[1]-be[3]) + al[3]*(be[2]-be[1]);
    bb2 := be[1]*(al[3]-al[2]) + be[2]*(al[1]-al[3]) + be[3]*(al[2]-al[1]);
    Df := ((al[1]-al[2])*(al[1]-al[3])*(al[2]-al[3]))^2;
    Dg := ((be[1]-be[2])*(be[1]-be[3])*(be[2]-be[3]))^2;
    AA := Dg*aa1/aa2; BB := Df*bb1/bb2;
    return -(AA*(al[2]-al[1])*(al[1]-al[3])*X^2 + BB*(be[2]-be[1])*(be[1]-be[3]))
           *(AA*(al[3]-al[2])*(al[2]-al[1])*X^2 + BB*(be[3]-be[2])*(be[2]-be[1]))
           *(AA*(al[1]-al[3])*(al[3]-al[2])*X^2 + BB*(be[1]-be[3])*(be[3]-be[2]));
end function;
// pullbacks of points (X, Y) of the Weierstrass models (scaled by c) to the Jacobian of y^2 = F(x^2)
function Pull1(JJ, Pt, c) PP := Parent((JJ!0)[1]); xx := PP.1; XX := Pt[1]/c; YY := Pt[2]/c; return elt< JJ | xx^2 - XX, PP!YY >; end function;
function Pull2(JJ, Pt, c) PP := Parent((JJ!0)[1]); xx := PP.1; XX := Pt[1]/c; YY := Pt[2]/c; return elt< JJ | xx^2 - 1/XX, (YY/XX)*xx >; end function;
// The largest finite abelian group that embeds in every group with the given invariants: for each prime l,
// the sorted l-primary exponents of a subgroup are bounded componentwise by those of the group.
function StructureBound(invlist)
    ls := &meet[ Seqset(&cat[ PrimeDivisors(n) : n in inv | n gt 1 ]) : inv in invlist ];
    out := [Integers()| ];
    for l in ls do
        es := [ Reverse(Sort([ Valuation(n, l) : n in inv | n mod l eq 0 ])) : inv in invlist ];
        r := Min([ #e : e in es ]);
        out cat:= [ l^Min([ e[i] : e in es ]) : i in [1..r] ];
    end for;
    return Invariants(AbelianGroup(out));
end function;
// f2 = h(lam x)/mu^2 for some lam, mu in the base field (h, f2 even sextics)
function SameUpToScaling(h, f2)
    okl, lam := IsSquare(Coefficient(f2, 2)*Coefficient(h, 0)/(Coefficient(f2, 0)*Coefficient(h, 2)));
    okm, mu := IsSquare(Coefficient(h, 0)/Coefficient(f2, 0));
    return okl and okm and &and[ Coefficient(f2, i)*mu^2 eq Coefficient(h, i)*lam^i : i in [0..6] ];
end function;

Fs<s> := FunctionField(Q); P<x> := PolynomialRing(Fs);
f := (-185971387536*s^13 + 101438938656*s^11 - 3794947200*s^10 - 16906489776*s^9 + 1379980800*s^8 + 890087616*s^7 - 104544000*s^6 - 13642992*s^5 + 2956800*s^4 - 142560*s^3 - 28800*s^2 + 3888*s)*x^6 + (104608905489*s^12 - 57059402994*s^10 - 2134657800*s^9 + 9509900499*s^8 + 776239200*s^7 - 960389100*s^6 - 58806000*s^5 + 91258695*s^4 - 13543200*s^3 + 1346622*s^2 - 16200*s - 2187)*x^4 + (2401490025*s^8 - 873269100*s^6 + 338879376*s^5 + 66156750*s^4 - 61614432*s^3 + 8320900*s^2 - 933552*s + 18225)*x^2 - 45163008*s^4 + 8211456*s^2 - 921600*s + 124416;
assert &and[ Coefficient(f, i) eq 0 : i in [1, 3, 5] ] and Degree(f) eq 6 and IsSquarefree(f);
printf "C_s : y^2 = f_s(x) = %o\n", f;

// ---------------------------------------------------------------- (0) source: the HLP gluing of 66c1 and E_6(u), u = (33 s^2 - 1)/(9 - 33 s^2)
E1 := EllipticCurve([Fs| 1, 0, 0, -45, 81]);
assert Conductor(EllipticCurve([Q| 1, 0, 0, -45, 81])) eq 66;
P10 := E1![0, 9];
assert 10*P10 eq E1!0 and &and[ (10 div l)*P10 ne E1!0 : l in PrimeDivisors(10) ];
u := (33*s^2 - 1)/(9 - 33*s^2);
E2 := EllipticCurve([Fs| 1 - u, -(u^2 + u), -(u^2 + u), 0, 0]);
P6 := E2![0, 0];
assert 6*P6 eq E2!0 and 2*P6 ne E2!0 and 3*P6 ne E2!0;
g1 := MonicCubic(E1, x); g2 := MonicCubic(E2, x);
lin1 := [ g[1] : g in Factorization(g1) | Degree(g[1]) eq 1 ]; q1 := [ g[1] : g in Factorization(g1) | Degree(g[1]) eq 2 ];
lin2 := [ g[1] : g in Factorization(g2) | Degree(g[1]) eq 1 ]; q2 := [ g[1] : g in Factorization(g2) | Degree(g[1]) eq 2 ];
assert #lin1 eq 1 and #q1 eq 1 and #lin2 eq 1 and #q2 eq 1;
ok1, w1 := IsSquare(Discriminant(q1[1])/(33)); ok2, w2 := IsSquare(Discriminant(q2[1])/(33));
assert ok1 and ok2 and not IsSquare(Discriminant(q1[1]));
K<z> := ext< Fs | x^2 - (33) >;
al := [ K!(-Coefficient(lin1[1], 0)), (-Coefficient(q1[1], 1) + w1*z)/2, (-Coefficient(q1[1], 1) - w1*z)/2 ];
be := [ K!(-Coefficient(lin2[1], 0)), (-Coefficient(q2[1], 1) + w2*z)/2, (-Coefficient(q2[1], 1) - w2*z)/2 ];
assert &and[ Evaluate(PolynomialRing(K)!g1, a) eq 0 : a in al ] and &and[ Evaluate(PolynomialRing(K)!g2, b) eq 0 : b in be ];
hK := HLPSextic(al, be);
assert &and[ IsCoercible(Fs, c) : c in Coefficients(hK) ];
h := P![ Fs!c : c in Coefficients(hK) ];
assert SameUpToScaling(h, f);
print "(0) the HLP (Prop. 4) gluing of E1 = 66c1 and E2 = E_6((33 s^2 - 1)/(9 - 33 s^2)) along psi has coefficients in Q(s),";
print "    and f_s = h(lambda x)/mu^2 with lambda, mu in Q(s)";

// ---------------------------------------------------------------- (1) the two elliptic quotients
F := P![ Coefficient(f, 2*i) : i in [0..3] ];
Fr := P![ Coefficient(F, 3 - i) : i in [0..3] ];
Ea, ca := MonEC(F); Eb, cb := MonEC(Fr);
oka, ma := IsIsomorphic(E2, Ea);
okb, mb := IsIsomorphic(E1, Eb);
assert oka and okb;
print "(1) phi1 = (x^2, y) : C_s -> Ea : Y^2 = F(X), phi2 = (1/x^2, y/x^3) : C_s -> Eb : Y^2 = X^3 F(1/X);";
print "    Ea is Q(s)-isomorphic to E2 = E_6(u) and Eb to E1 = 66c1: J is isogenous over Q to E1 x E2";

// ---------------------------------------------------------------- (2) the marked group over Q(s)
J := Jacobian(HyperellipticCurve(f));
D10 := Pull2(J, mb(P10), cb);
D6 := Pull1(J, ma(P6), ca);
assert 10*D10 eq J!0 and 5*D10 ne J!0 and 2*D10 ne J!0;
assert 6*D6 eq J!0 and 3*D6 ne J!0 and 2*D6 ne J!0 and 5*D10 eq 3*D6;
D30 := D10 + 2*D6;                         // (D10 + D6 has order 15: 15 (D10 + D6) = 5 D10 + 3 D6 = 2 (5 D10) = 0)
assert 30*D30 eq J!0 and &and[ (30 div l)*D30 ne J!0 : l in [2, 3, 5] ];
print "(2) over Q(s): D10 = phi2^*(P10) has order 10, D6 = phi1^*(P6) has order 6, 5 D10 = 3 D6, and D30 = D10 + 2 D6 has order 30";

// ---------------------------------------------------------------- (3) the gain locus
Kr<r> := QuadraticField(33);
FK<S> := FunctionField(Kr); PFK := PolynomialRing(Kr);
ev := func< a | Evaluate(Numerator(a), S)/Evaluate(Denominator(a), S) >;
function ToFK(c) cs := Eltseq(c); return ev(cs[1]) + ev(cs[2])*r; end function;     // K = Q(s)(z) -> Q(sqrt 33)(S), z -> r
assert r^2 eq 33;
N := ToFK((be[1] - be[2])*(al[1] - al[2]));
nu := PFK!Numerator(N); de := PFK!Denominator(N);
fnu := Factorization(nu); fde := Factorization(de);
assert &and[ Degree(g[1]) eq 1 and not IsCoercible(Q, Roots(g[1])[1][1]) : g in fnu cat fde ];   // no rational zeros or poles
oddp := (LeadingCoefficient(nu)/LeadingCoefficient(de))*&*[ PFK | g[1] : g in fnu cat fde | IsOdd(g[2]) ];
assert Degree(oddp) eq 2;
oksq, gg := IsSquare(N/Evaluate(oddp, S));
assert oksq;
printf "(3) N(s) = (beta_1 - beta_2)(alpha_1 - alpha_2) = c (s - a)(s - b) g(s)^2 over Q(sqrt 33)(s), with c (s - a)(s - b) = %o;\n", oddp;
print "    all zeros and poles of N are irrational: for rational s, N(s) is a square in Q(sqrt 33) iff c (s - a)(s - b) is";
c0, c1, c2 := Explode(Coefficients(oddp));
P3<SS, W1, W2, ZZ> := ProjectiveSpace(Q, 3);
q1X := W1^2 + 33*W2^2 - (Eltseq(c0)[1]*ZZ^2 + Eltseq(c1)[1]*SS*ZZ + Eltseq(c2)[1]*SS^2);
q2X := 2*W1*W2 - (Eltseq(c0)[2]*ZZ^2 + Eltseq(c1)[2]*SS*ZZ + Eltseq(c2)[2]*SS^2);
X := Curve(P3, [q1X, q2X]);
assert Genus(X) eq 1;
printf "    X : (w1 + w2 sqrt 33)^2 = c (s - a)(s - b), i.e. %o = 0 and %o = 0 in P^3(s : w1 : w2 : 1); genus 1\n", q1X, q2X;
Ep := EllipticCurve([Q| 0, 0, 0, 2904, -74536]);
assert Conductor(Ep) eq 69696;
PE<xe, ye, ze> := CoordinateRing(Ambient(Ep));
Wp := DefiningPolynomial(Ep);
dpz := [ 29792*xe^3*ye + 96864768*xe^3*ze + 1896048*xe^2*ye*ze + 180576*xe*ye^2*ze - 18506*ye^3*ze + 5485423680*xe^2*ze^2 + 643840032*xe*ye*ze^2 - 53714232*ye^2*ze^2 + 91851522048*xe*ze^3 - 14887075280*ye*ze^3 + 5045432901696*ze^4,
         150528*xe^3*ye + 116379648*xe^3*ze + 2076624*xe^2*ye*ze - 180576*xe*ye^2*ze - 182505*ye^3*ze + 10184258304*xe^2*ze^2 + 126611496*xe*ye*ze^2 - 214141356*ye^2*ze^2 + 606079517472*xe*ze^3 + 68465798808*ye*ze^3 + 10963353680928*ze^4,
         -632016*xe^2*ye*ze - 60192*xe*ye^2*ze + 1881*ye^3*ze - 2054916864*xe^2*ze^2 - 273670056*xe*ye*ze^2 + 1955052*ye^2*ze^2 - 67905405216*xe*ze^3 - 13506274584*ye*ze^3 + 696410760672*ze^4,
         41382*ye^3*ze + 192648456*ye^2*ze^2 + 930369255552*xe*ze^3 + 146693684016*ye*ze^3 + 18957810535488*ze^4 ];
for q in [q1X, q2X] do
    assert IsDivisibleBy(Evaluate(q, dpz), Wp);
end for;
print "    phi : E' = [0,0,0,2904,-74536] (69696gl1) -> X by four integral quartics: both quadrics pull back to multiples of the Weierstrass cubic";
Pg := Ep![385, 7623]; Tt := Ep![22, 0];
assert Order(Pg) eq 0 and 2*Tt eq Ep!0 and Tt ne Ep!0;
assert RankBounds(Ep) ge 1;
function SVal(R) v := Eltseq(R); return Evaluate(dpz[1], v)/Evaluate(dpz[4], v); end function;
function NAt(sv) return Evaluate(nu, Kr!sv)/Evaluate(de, Kr!sv); end function;
assert Evaluate(dpz[4], Eltseq(Tt)) eq 39425934157632 and SVal(Tt) eq 3/11;
printf "    P = %o has infinite order, T = %o has order 2; s(T) = %o\n", Pg, Tt, SVal(Tt);

// ---------------------------------------------------------------- (4) exactness at 5, (5) the members R_k = T + 8kP
p := 5; r0 := 3;
fbar := ReductionClass(f, p, r0);
inv5 := Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
n, chi := FrobData(fbar);
assert inv5 eq [2, 30];
printf "(4) p = 5, s = 3 mod 5: reduction %o (squarefree, degree %o); J(F_5) = %o, chi = %o\n", fbar, Degree(fbar), inv5, chi;
assert Conductor(Ep) mod p ne 0 and Order(ChangeRing(Ep, GF(p))!Pg) eq 8;
assert Integers()!Evaluate(dpz[4], Eltseq(Tt)) mod p ne 0;
assert &and[ IsCoercible(Integers(), cc) : cc in &cat[ Coefficients(d) : d in dpz ] ];
print "(5) E' has good reduction at 5 and P has order 8 mod 5: R_k = T + 8kP reduces to T mod 5, the quartics are integral and";
print "    the fourth is a 5-unit at T, so s(R_k) is 5-integral and = s(T) = 3/11 = 3 mod 5";
memb := [];
for k in [0, 1, -1] do
    sv := SVal(Tt + 8*k*Pg);
    assert Valuation(sv, p) ge 0 and GF(p)!sv eq GF(p)!r0;
    assert PolynomialRing(GF(p))!SpecPoly(f, sv) eq fbar;
    okk, _ := IsSquare(NAt(sv)); assert okk and NAt(sv) ne 0;
    Append(~memb, sv);
    printf "    k = %o: s(R_k) has %o digits, = 3 mod 5, N(s) is a nonzero square in Q(sqrt 33)\n", k, #Sprint(sv);
end for;
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, sv))) : sv in memb } eq #memb;
print "    the G2-invariants at s(T), s(T + 8P), s(T - 8P) are pairwise distinct: infinitely many isomorphism classes";

// ---------------------------------------------------------------- (6) sanity
C230a := HyperellipticCurve(PolynomialRing(Q)![12, -36, 18, 23, -13, -5, 1], PolynomialRing(Q)![0, 1, 1]);
assert IsIsomorphic(HyperellipticCurve(SpecPoly(f, 3/11)), C230a);
print "(6) the member s = 3/11 (in the class) is Q-isomorphic to the census curve 2.30.a";
for R in [Tt, Pg, Ep![25, 117]] do
    sv := SVal(R);
    okk, _ := IsSquare(NAt(sv)); assert okk;
    invm := TorsionInvariants(SpecPoly(f, sv));
    printf "(6) s = %o (%o mod 5): N(s) is a square in Q(sqrt 33); TorsionSubgroup %o\n", sv, GF(5)!sv, invm;
    assert invm eq [2, 30];
end for;
printf "time %o s\n", Cputime(T0);
print "FAMILY_2_30_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E1 x E2) and J(Q)_tors exactly [2,30]";
quit;
