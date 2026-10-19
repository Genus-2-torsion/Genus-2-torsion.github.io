/* family_3_24_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E1 x E2) with J(Q)_tors exactly [3,24].
   Found 2026-10-19 (torsion_jac notes/hlpP2_2026_10_19.md), as a by-product of the [3,12] certificate:
   the census curve 3.24.a (Booker-Sutherland; LMFDB 8100.e.4) is the member s = -1 of the [3,12] family.

   Source.  The [3,12] family of family_3_12_split_certificate.m: the Howe-Leprevost-Poonen (2,2)-gluing
   (HLP 2000, Prop. 4) of E1 = 90c3 = [1,-1,1,-122,1721] (torsion Z/12, P12 = (-9,-41), P4 = 3 P12) and
       E2 = E_6(u) : y^2 + (1 - u)xy - (u^2 + u)y = x^3 - (u^2 + u)x^2,   u = (-15 s^2 - 1)/(9 + 15 s^2),
   along psi : E1[2] -> E2[2] (alpha_i -> beta_i; alpha_1, beta_1 the rational roots of the 2-division cubics,
   alpha_2 = (-b1 + w1 sqrt -15)/2, beta_2 = (-b2 + w2 sqrt -15)/2 as computed below), model C_s : y^2 = f_s(x).
   The gain.  By HLP 2000, Prop. 12, a point (P0, Q0) of E1 x E2 with 2P0 = P4 and 2Q0 = T' = (beta_1, 0) maps
   to a rational point of J = (E1 x E2)/graph(psi) iff the Kummer images agree under psi; with
   L = Q x Q(sqrt -15) this is the single condition
       N(s) = (beta_1 - beta_2)(x(P4) - alpha_2)  is a square in K = Q(sqrt -15)
   (the Q-component follows by taking norms).  Over K(s), N = c (s - a)(s - b) g(s)^2 with a, b, and all zeros and
   poles of g, in K \ Q.  So for rational s the condition is (w1 + w2 sqrt -15)^2 = c (s - a)(s - b) with w1, w2
   in Q: the intersection X of two quadrics in P^3 (embedded below), a curve of genus 1.  The map
       phi : E' -> X,  E' = 4800cb1 : y^2 = x^3 + x^2 + 67x - 237  (rank 1, E'(Q) = Z/2 x Z, P = (273/16, 5025/64)),
   given by four quartics (embedded below; checked to map E' into X), sends T = (3, 0) to the point over
   s = -1.  So s(R) = phi(R)_1/phi(R)_4, R in E'(Q), are infinitely many members with the gain.
   The gained class D = image of (P0, Q0) has 2D = image of (P4, T') = image of (3 P4, 0) (the rational
   2-torsion points are glued), of order 4 and in <D12>; so D has order 8 and <D12, D6, D> = [3,24].

   Argument.
   (1) Splitting and (2) the marked [3,12] over Q(s): as in family_3_12_split_certificate.m (repeated below).
   (3) The gain locus.  N, its odd part c (s - a)(s - b) and the quadrics of X are computed from the root labels
       of step (0); N/(c (s - a)(s - b)) is a square in K(s), and a, b, the zeros and poles of its square root
       are not rational.  phi maps E' into X (the quadrics pull back to multiples of the Weierstrass cubic).
       At every member below N(s) is checked to be a nonzero square in K, and at s = -1 (census 3.24.a) the
       torsion subgroup (TorsionSubgroup) contains a class G8 of order 8 with 2 G8 in {3 D12, 9 D12}, as predicted.
       So [3,24] embeds in J_s(Q)_tors for every s = s(R), R in E'(Q), at which C_s is smooth.
   (4) Exactness.  For every 7-integral s = -1 (mod 7), C_s has the reduction of C_{-1} at 7, with
       J(F_7) = [3,24].  Reduction is injective on J_s(Q)_tors, so J_s(Q)_tors embeds in [3,24].
   (5) Infinitely many members in the class.  E' has good reduction at 7 and P has order 4 mod 7, so
       R_k = T + 4kP reduces to T mod 7.  The quartics have integral coefficients and the fourth is a 7-unit at T,
       so s(R_k) is 7-integral and = s(T) = -1 (mod 7).  For every k, C_{s(R_k)} therefore has J(Q)_tors exactly
       [3,24].  The R_k are distinct and s has finite degree on E', so infinitely many s occur; the G2-invariants at
       s(T) and s(T + 4P) differ, so s -> M_2 is nonconstant with finite fibres: infinitely many isomorphism classes.
   (6) Sanity: TorsionSubgroup at s = -1 (= 3.24.a, in the class) and at s(P) = 19/173 (on the gain locus,
       outside the class) is [3,24].
   Run (from the repository root):
       magma -b pipeline/magma/family_3_24_split_certificate.m > data/knowledge/sources/family_3_24_split_certificate.log */
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
f := (-303750000*s^13 - 364500000*s^11 + 131760000*s^10 - 133650000*s^9 + 105408000*s^8 - 15480000*s^7 + 17568000*s^6 - 522000*s^5 + 1093120*s^4 + 12000*s^3 + 23424*s^2 + 720*s)*x^6 + (284765625*s^12 + 341718750*s^10 + 123525000*s^9 + 125296875*s^8 + 98820000*s^7 - 81852300*s^6 + 16470000*s^5 - 38056545*s^4 - 8344800*s^3 + 1273614*s^2 + 21960*s - 675)*x^4 + (-231609375*s^8 - 185287500*s^6 + 35046000*s^5 - 30881250*s^4 + 14018400*s^3 + 30848676*s^2 - 467280*s - 41175)*x^2 + 51840000*s^4 + 20736000*s^2 - 22487040*s - 691200;
assert &and[ Coefficient(f, i) eq 0 : i in [1, 3, 5] ] and Degree(f) eq 6 and IsSquarefree(f);
printf "C_s : y^2 = f_s(x) = %o\n", f;

// ---------------------------------------------------------------- (0) source: the HLP gluing of 90c3 and E_6(u), u = (-15 s^2 - 1)/(9 + 15 s^2)
E1 := EllipticCurve([Fs| 1, -1, 1, -122, 1721]);
assert Conductor(EllipticCurve([Q| 1, -1, 1, -122, 1721])) eq 90;
P12 := E1![-9, -41];
assert 12*P12 eq E1!0 and &and[ (12 div l)*P12 ne E1!0 : l in PrimeDivisors(12) ];
u := (-15*s^2 - 1)/(9 + 15*s^2);
E2 := EllipticCurve([Fs| 1 - u, -(u^2 + u), -(u^2 + u), 0, 0]);
P6 := E2![0, 0];
assert 6*P6 eq E2!0 and 2*P6 ne E2!0 and 3*P6 ne E2!0;
g1 := MonicCubic(E1, x); g2 := MonicCubic(E2, x);
lin1 := [ g[1] : g in Factorization(g1) | Degree(g[1]) eq 1 ]; q1 := [ g[1] : g in Factorization(g1) | Degree(g[1]) eq 2 ];
lin2 := [ g[1] : g in Factorization(g2) | Degree(g[1]) eq 1 ]; q2 := [ g[1] : g in Factorization(g2) | Degree(g[1]) eq 2 ];
assert #lin1 eq 1 and #q1 eq 1 and #lin2 eq 1 and #q2 eq 1;
ok1, w1 := IsSquare(Discriminant(q1[1])/(-15)); ok2, w2 := IsSquare(Discriminant(q2[1])/(-15));
assert ok1 and ok2 and not IsSquare(Discriminant(q1[1]));
K<z> := ext< Fs | x^2 - (-15) >;
al := [ K!(-Coefficient(lin1[1], 0)), (-Coefficient(q1[1], 1) + w1*z)/2, (-Coefficient(q1[1], 1) - w1*z)/2 ];
be := [ K!(-Coefficient(lin2[1], 0)), (-Coefficient(q2[1], 1) + w2*z)/2, (-Coefficient(q2[1], 1) - w2*z)/2 ];
assert &and[ Evaluate(PolynomialRing(K)!g1, a) eq 0 : a in al ] and &and[ Evaluate(PolynomialRing(K)!g2, b) eq 0 : b in be ];
hK := HLPSextic(al, be);
assert &and[ IsCoercible(Fs, c) : c in Coefficients(hK) ];
h := P![ Fs!c : c in Coefficients(hK) ];
assert SameUpToScaling(h, f);
print "(0) the HLP (Prop. 4) gluing of E1 = 90c3 and E2 = E_6((-15 s^2 - 1)/(9 + 15 s^2)) along psi has coefficients in Q(s),";
print "    and f_s = h(lambda x)/mu^2 with lambda, mu in Q(s)";

// ---------------------------------------------------------------- (1) the two elliptic quotients
F := P![ Coefficient(f, 2*i) : i in [0..3] ];
Fr := P![ Coefficient(F, 3 - i) : i in [0..3] ];
Ea, ca := MonEC(F); Eb, cb := MonEC(Fr);
oka, ma := IsIsomorphic(E2, Ea);
okb, mb := IsIsomorphic(E1, Eb);
assert oka and okb;
print "(1) phi1 = (x^2, y) : C_s -> Ea : Y^2 = F(X), phi2 = (1/x^2, y/x^3) : C_s -> Eb : Y^2 = X^3 F(1/X);";
print "    Ea is Q(s)-isomorphic to E2 = E_6(u) and Eb to E1 = 90c3: J is isogenous over Q to E1 x E2";

// ---------------------------------------------------------------- (2) the marked group over Q(s)
J := Jacobian(HyperellipticCurve(f));
D12 := Pull2(J, mb(P12), cb);
D6 := Pull1(J, ma(P6), ca);
assert 12*D12 eq J!0 and 6*D12 ne J!0 and 4*D12 ne J!0;
assert 6*D6 eq J!0 and 3*D6 ne J!0 and 2*D6 ne J!0 and 6*D12 eq 3*D6;
assert 2*D6 notin { 4*D12, 8*D12 };          // the two elements of order 3 of <D12>
print "    <D12, D6> has order 12 * 6 / 2 = 36, 3-rank 2 and a cyclic 4-part: it is [3,12]";
print "(2) over Q(s): D12 = phi2^*(P12) has order 12, D6 = phi1^*(P6) has order 6, 6 D12 = 3 D6, 2 D6 not in <D12>";

// ---------------------------------------------------------------- (3) the gain locus
Kr<r> := QuadraticField(-15);
FK<S> := FunctionField(Kr); PFK := PolynomialRing(Kr);
ev := func< a | Evaluate(Numerator(a), S)/Evaluate(Denominator(a), S) >;
function ToFK(c) cs := Eltseq(c); return ev(cs[1]) + ev(cs[2])*r; end function;     // K = Q(s)(z) -> Q(sqrt -15)(S), z -> r
assert r^2 eq -15;
P4 := 3*P12;
assert 4*P4 eq E1!0 and 2*P4 ne E1!0;
N := ToFK((be[1] - be[2])*(K!P4[1] - al[2]));
nu := PFK!Numerator(N); de := PFK!Denominator(N);
fnu := Factorization(nu); fde := Factorization(de);
assert &and[ Degree(g[1]) eq 1 and not IsCoercible(Q, Roots(g[1])[1][1]) : g in fnu cat fde ];   // no rational zeros or poles
oddp := (LeadingCoefficient(nu)/LeadingCoefficient(de))*&*[ PFK | g[1] : g in fnu cat fde | IsOdd(g[2]) ];
assert Degree(oddp) eq 2;
oksq, gg := IsSquare(N/Evaluate(oddp, S));
assert oksq;
printf "(3) N(s) = (beta_1 - beta_2)(x(P4) - alpha_2) = c (s - a)(s - b) g(s)^2 over Q(sqrt -15)(s), with c (s - a)(s - b) = %o;\n", oddp;
print "    all zeros and poles of N are irrational: for rational s, N(s) is a square in Q(sqrt -15) iff c (s - a)(s - b) is";
c0, c1, c2 := Explode(Coefficients(oddp));
P3<SS, W1, W2, ZZ> := ProjectiveSpace(Q, 3);
q1X := W1^2 - 15*W2^2 - (Eltseq(c0)[1]*ZZ^2 + Eltseq(c1)[1]*SS*ZZ + Eltseq(c2)[1]*SS^2);
q2X := 2*W1*W2 - (Eltseq(c0)[2]*ZZ^2 + Eltseq(c1)[2]*SS*ZZ + Eltseq(c2)[2]*SS^2);
X := Curve(P3, [q1X, q2X]);
assert Genus(X) eq 1;
printf "    X : (w1 + w2 sqrt -15)^2 = c (s - a)(s - b), i.e. %o = 0 and %o = 0 in P^3(s : w1 : w2 : 1); genus 1\n", q1X, q2X;
Ep := EllipticCurve([Q| 0, 1, 0, 67, -237]);
assert Conductor(Ep) eq 4800;
PE<xe, ye, ze> := CoordinateRing(Ambient(Ep));
Wp := DefiningPolynomial(Ep);
dpz := [ 71968*xe^3*ye - 2611200*xe^3*ze - 558032*xe^2*ye*ze + 50400*xe*ye^2*ze - 75118*ye^3*ze - 7531200*xe^2*ze^2 + 13581856*xe*ye*ze^2 + 3393000*ye^2*ze^2 - 95030400*xe*ze^3 - 33616416*ye*ze^3 - 191625600*ze^4,
         249600*xe^3*ye - 15283200*xe^3*ze + 1383600*xe^2*ye*ze - 75600*xe*ye^2*ze - 244875*ye^3*ze + 49516800*xe^2*ze^2 + 33634200*xe*ye*ze^2 + 11842500*ye^2*ze^2 - 1199474400*xe*ze^3 - 114019200*ye*ze^3 + 2642918400*ze^4,
         327600*xe^2*ye*ze - 25200*xe*ye^2*ze + 2205*ye^3*ze - 14880000*xe^2*ze^2 - 738600*xe*ye*ze^2 + 884100*ye^2*ze^2 + 2580000*xe*ze^3 + 21782400*ye*ze^3 - 304320000*ze^4,
         3150*ye^3*ze + 7635000*ye^2*ze^2 - 498000000*xe*ze^3 - 102450000*ye*ze^3 + 2109000000*ze^4 ];
for q in [q1X, q2X] do
    assert IsDivisibleBy(Evaluate(q, dpz), Wp);
end for;
print "    phi : E' = [0,1,0,67,-237] (4800cb1) -> X by four integral quartics: both quadrics pull back to multiples of the Weierstrass cubic";
Pg := Ep![273/16, 5025/64]; Tt := Ep![3, 0];
assert Order(Pg) eq 0 and 2*Tt eq Ep!0 and Tt ne Ep!0;
assert RankBounds(Ep) ge 1;
function SVal(R) v := Eltseq(R); return Evaluate(dpz[1], v)/Evaluate(dpz[4], v); end function;
function NAt(sv) return Evaluate(nu, Kr!sv)/Evaluate(de, Kr!sv); end function;
assert Evaluate(dpz[4], Eltseq(Tt)) eq 615000000 and SVal(Tt) eq -1;
printf "    P = %o has infinite order, T = %o has order 2; s(T) = %o\n", Pg, Tt, SVal(Tt);

// ---------------------------------------------------------------- (4) exactness at 7, (5) the members R_k = T + 4kP
p := 7; r0 := 6;
fbar := ReductionClass(f, p, r0);
inv7 := Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
n, chi := FrobData(fbar);
assert inv7 eq [3, 24];
printf "(4) p = 7, s = 6 mod 7: reduction %o (squarefree, degree %o); J(F_7) = %o, chi = %o\n", fbar, Degree(fbar), inv7, chi;
assert Conductor(Ep) mod p ne 0 and Order(ChangeRing(Ep, GF(p))!Pg) eq 4;
assert Integers()!Evaluate(dpz[4], Eltseq(Tt)) mod p ne 0;
assert &and[ IsCoercible(Integers(), cc) : cc in &cat[ Coefficients(d) : d in dpz ] ];
print "(5) E' has good reduction at 7 and P has order 4 mod 7: R_k = T + 4kP reduces to T mod 7, the quartics are integral and";
print "    the fourth is a 7-unit at T, so s(R_k) is 7-integral and = s(T) = -1 mod 7";
memb := [];
for k in [0, 1, -1] do
    sv := SVal(Tt + 4*k*Pg);
    assert Valuation(sv, p) ge 0 and GF(p)!sv eq GF(p)!r0;
    assert PolynomialRing(GF(p))!SpecPoly(f, sv) eq fbar;
    okk, _ := IsSquare(NAt(sv)); assert okk and NAt(sv) ne 0;
    Append(~memb, sv);
    printf "    k = %o: s(R_k) has %o digits, = 6 mod 7, N(s) is a nonzero square in Q(sqrt -15)\n", k, #Sprint(sv);
end for;
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, sv))) : sv in memb } eq #memb;
print "    the G2-invariants at s(T), s(T + 4P), s(T - 4P) are pairwise distinct: infinitely many isomorphism classes";

// ---------------------------------------------------------------- (3), (6) the gain at members; sanity
f0 := SpecPoly(f, -1);
J0 := Jacobian(HyperellipticCurve(f0));
// a generator m D12 (m prime to 12) of <D12> whose Mumford representation over Q(s) is regular at s = -1
// (that of D12 itself has a pole there), specialized to the member
function RegularAt(D, s0) return &and[ Evaluate(Denominator(c), s0) ne 0 : c in Coefficients(D[1]) cat Coefficients(D[2]) ]; end function;
mreg := [ m : m in [1, 5, 7, 11] | RegularAt(m*D12, -1) ][1];
Dm := mreg*D12;
RX0 := PolynomialRing(Q);
D12_0 := elt< J0 | RX0![ Evaluate(c, -1) : c in Coefficients(Dm[1]) ], RX0![ Evaluate(c, -1) : c in Coefficients(Dm[2]) ] >;
printf "    (%o D12 is regular at s = -1)\n", mreg;
assert Order(D12_0) eq 12;
TG, mTG := TorsionSubgroup(J0);
assert Invariants(TG) eq [3, 24];
G8s := [ mTG(g) : g in TG | Order(g) eq 8 ];
assert #G8s gt 0 and &and[ 2*G in { 3*D12_0, 9*D12_0 } : G in G8s ];
print "(3) at s = -1: J(Q)_tors = [3,24]; every class G8 of order 8 has 2 G8 in {3 D12, 9 D12}, as Prop. 12 predicts";
C324a := HyperellipticCurve(PolynomialRing(Q)![4, 6, -4, -6, 5, -3, 1], PolynomialRing(Q)![0, 1, 1]);
assert IsIsomorphic(HyperellipticCurve(f0), C324a);
print "    the member s = -1 is Q-isomorphic to the census curve 3.24.a";
for R in [Pg] do
    sv := SVal(R);
    okk, _ := IsSquare(NAt(sv)); assert okk;
    invm := TorsionInvariants(SpecPoly(f, sv));
    printf "(6) s = %o (%o mod 7): N(s) is a square in Q(sqrt -15); TorsionSubgroup %o\n", sv, GF(7)!sv, invm;
    assert invm eq [3, 24];
end for;
printf "time %o s\n", Cputime(T0);
print "FAMILY_3_24_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E1 x E2) and J(Q)_tors exactly [3,24]";
quit;
