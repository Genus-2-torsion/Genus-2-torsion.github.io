/* family_3_6_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E x E1) with J(Q)_tors exactly [3,6].
   Found 2026-10-19 (torsion_jac notes/isopair_2026_10_19.md).

   Source.  A twisted Howe-Leprevost-Poonen (2,2)-gluing (HLP 2000, Prop. 4) of the curves of the [6,6] certificate:
       E : y^2 + 3 lam xy + (lam^3 - 1)y = x^3,  lam = (2t^3 + 1)/(3t^2)     (P = (0,0) of order 3, T of order 2),
   and E1 = E/C, C = mu_3, but glued along phi o s : E[2] -> E1[2], where phi : E -> E1 is the 3-isogeny and s is the
   automorphism of E[2] that fixes T and swaps the two other points of order 2 (it is Galois-equivariant; HLP use the same
   swap in the proof of Corollary 6).  The pullbacks still give Z/6 x Z/3, but the swap destroys the 2-torsion gain of the
   untwisted gluing, so J[2](Q) has rank 1.  The gluing has coefficients in Q(t); after removing square factors and
   scaling x it is the curve f_t below, f_t = F(x^2) = lc (x^2 - r) F2(x^2) over Q(t) with F2 quadratic, and
   phi1 : (x,y) -> (x^2, y) maps C_t onto Ea : Y^2 = F(X), phi2 : (x,y) -> (1/x^2, y/x^3) onto Eb : Y^2 = X^3 F(1/X);
   below, Ea is Q(t)-isomorphic to E and Eb to E1.
   Argument.
   (0) E, P, T, C = mu_3, E1, phi, as in the [6,6] certificate.
   (1) Splitting and gluing.  Ea = E and Eb = E1 over Q(t); phi1^* + phi2^* : Ea x Eb -> J is an isogeny over Q with
       kernel of order 4.  The gluing is along phi o s on E[2] (checked on the Weierstrass points).
   (2) The marked group over Q(t).  Da = phi1^*((P + T)') has order 6, Db = phi2^*(phi(P)') has order 3, and the 18
       classes i Da + j Db are distinct, so G = Z/6 x Z/3 = [3,6] is marked over Q(t) and embeds in J_t(Q)_tors.
   (3) Exactness on a residue class: the 935-adic class t = 234 (mod 5 * 11 * 17), t integral at 5, 11, 17.
       - At p = 5 (t = 4 mod 5) C_t has the reduction of C_4, with J(F_5) = [6,6].  So the odd part of J_t(Q)_tors is
         the marked [3,3], and its 2-part is elementary abelian, i.e. it is J_t(Q)[2].
       - J[2] is the group of even subsets of the six roots of f_t modulo complement, as a Galois module.  Over Q(t),
         f = lc (x^2 - r) F2(x^2), F2 = c2 X^2 + c1 X + c0, and c0/c2 = n^2 is a square in Q(t).  Label the roots
         sqrt r, -sqrt r, sqrt X2, -sqrt X2, sqrt X3, -sqrt X3 with sqrt X2 sqrt X3 = n.  The Galois action on them is then
         determined by the quadratic characters of r, d1 = c1^2 - 4 c0 c2 = c2^2 (X2 - X3)^2 and
         d2 = -c1/c2 + 2n = (sqrt X2 + sqrt X3)^2.  At p = 5, 11, 17 the values r, d1, d2 are p-adic units on the class,
         and their Legendre symbols give three Frobenius elements of (Z/2)^3 that generate it; only the trivial class
         and [x^2 - r, 0] are invariant under them (checked by enumerating the 16 classes).  So J_t(Q)[2] = Z/2.
       Hence J_t(Q)_tors is exactly [3,6].  (#J(F_p) = #E(F_p)^2 is a square and J(F_p)[2] has rank >= 2 at every p
       in the scan, so neither a gcd of point counts nor group structures alone can certify the 2-rank.)
   (4) Infinitely many.  The G2-invariants of the members differ, so t -> M_2 is nonconstant with finite fibres.
   (5) Sanity: TorsionSubgroup at t = 1/4, 234, -701 is [3,6]; the member t = 1/4 is the census curve 3.6.b.
   Run (from the repository root):
       magma -b pipeline/magma/family_3_6_split_certificate.m > data/knowledge/sources/family_3_6_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
f := (-6144*t^23 - 9216*t^22 - 35328*t^21 + 18432*t^20 - 32256*t^19 + 141696*t^18 - 62208*t^17 + 175104*t^16 - 225216*t^15 + 140736*t^14 - 240192*t^13 + 178800*t^12 - 141120*t^11 + 130752*t^10 - 72144*t^9 + 53568*t^8 - 24768*t^7 + 12240*t^6 - 3264*t^5 + 576*t^4 - 48*t^3)*x^6 + (-12288*t^24 - 30720*t^23 + 39936*t^22 + 425472*t^21 + 1717248*t^20 + 3337728*t^19 + 6301824*t^18 + 6992640*t^17 + 7895808*t^16 + 6741312*t^15 + 2701248*t^14 + 5370816*t^13 - 1990128*t^12 + 4271808*t^11 - 2171712*t^10 + 1966416*t^9 - 821088*t^8 + 421920*t^7 - 149880*t^6 + 49776*t^5 - 14208*t^4 + 3300*t^3 - 564*t^2 + 60*t - 3)*x^4 + (24576*t^24 + 374784*t^23 + 1864704*t^22 + 6259200*t^21 + 15753216*t^20 + 29253120*t^19 + 49529472*t^18 + 56095488*t^17 + 74174976*t^16 + 47838912*t^15 + 59512128*t^14 + 13627200*t^13 + 29967504*t^12 - 6181440*t^11 + 12720768*t^10 - 6334032*t^9 + 4365504*t^8 - 2026656*t^7 + 800664*t^6 - 248880*t^5 + 57504*t^4 - 9012*t^3 + 804*t^2 - 12*t - 3)*x^2 - 12288*t^24 + 55296*t^23 - 125952*t^22 + 330240*t^21 - 875520*t^20 + 1651200*t^19 - 2397312*t^18 + 3034368*t^17 - 3135744*t^16 + 2470848*t^15 - 1688256*t^14 + 964416*t^13 - 28176*t^12 - 597312*t^11 + 563712*t^10 - 291888*t^9 + 57888*t^8 + 77184*t^7 - 68736*t^6 + 12960*t^5 + 5232*t^4 - 2688*t^3 + 576*t^2 - 48*t;
assert &and[ Coefficient(f, i) eq 0 : i in [1, 3, 5] ] and Degree(f) eq 6 and IsSquarefree(f);
F := P![ Coefficient(f, 2*i) : i in [0..3] ];            // f = F(x^2)
printf "C_t : y^2 = f_t(x) = %o\n", f;
printf "    the irreducible factors of f_t over Q(t) have degrees %o\n", [ Degree(h[1]) : h in Factorization(f) ];

// ---------------------------------------------------------------- helpers
function MonicCubic(E)                       // the x-cubic of E after completing the square in y
    a1, a2, a3, a4, a6 := Explode(aInvariants(E));
    PP := PolynomialRing(BaseRing(E)); X := PP.1;
    return X^3 + (a2 + a1^2/4)*X^2 + (a4 + a1*a3/2)*X + (a6 + a3^2/4);
end function;
function MonEC(G)                            // a Weierstrass model of Y^2 = G(X), via X' = lc X, Y' = lc Y
    cs := [ Coefficient(G, i) : i in [0..3] ];
    return EllipticCurve([0, cs[3], 0, cs[2]*cs[4], cs[1]*cs[4]^2]), cs[4];
end function;
// pullbacks to J of points of Ea = MonEC(F) and Eb = MonEC(Fr), along phi1 = (x^2, y) and phi2 = (1/x^2, y/x^3)
function Pull1(JJ, Pt, c) PP := Parent((JJ!0)[1]); xx := PP.1; XX := Pt[1]/c; YY := Pt[2]/c; return elt< JJ | xx^2 - XX, PP!YY >; end function;
function Pull2(JJ, Pt, c) PP := Parent((JJ!0)[1]); xx := PP.1; XX := Pt[1]/c; YY := Pt[2]/c; return elt< JJ | xx^2 - 1/XX, (YY/XX)*xx >; end function;
// x(m(T)) for the point T of order 2 of the domain of m with x(T) = X, X in a residue ring Ft[x]/(g)
function XImage(m, X)
    a := aInvariants(Domain(m));
    vals := [ Evaluate(d, [X, -(a[1]*X + a[3])/2, Parent(X)!1]) : d in DefiningPolynomials(m) ];
    ok, inv := IsInvertible(vals[3]); assert ok;
    return vals[1]*inv;
end function;
// does the finite abelian group with invariants g embed in the one with invariants h (prime by prime)
function EmbedsIn(g, h)
    for l in &join[ Seqset(PrimeDivisors(n)) : n in g cat h | n gt 1 ] do
        a := Reverse(Sort([ Valuation(n, l) : n in g | n mod l eq 0 ]));
        b := Reverse(Sort([ Valuation(n, l) : n in h | n mod l eq 0 ]));
        if #a gt #b or &or[ a[i] gt b[i] : i in [1..#a] ] then return false; end if;
    end for;
    return true;
end function;

// ---------------------------------------------------------------- (0) source: E (Hesse form with a point of order 2), C = mu_3, E1 = E/C
lam := (2*t^3 + 1)/(3*t^2); mm := (t^3 - 1)/(3*t^2);        // lam^3 - 3 mm^2 lam - 2 mm^3 = 1
E := EllipticCurve([Ft| 3*lam, 0, lam^3 - 1, 0, 0]);
P3 := E![0, 0];
assert 3*P3 eq E!0 and P3 ne E!0;
T2 := E![-mm^2, -(3*lam*(-mm^2) + lam^3 - 1)/2];
assert T2 ne E!0 and 2*T2 eq E!0;
P6 := P3 + T2;
assert 6*P6 eq E!0 and 2*P6 ne E!0 and 3*P6 ne E!0;
x0 := -(lam^2 + lam + 1);
assert Evaluate(DivisionPolynomial(E, 3), x0) eq 0;
bb := 3*lam*x0 + lam^3 - 1; dsc := bb^2 + 4*x0^3;           // the points over x0: y^2 + bb y - x0^3 = 0
okmu, w0 := IsSquare(-3*dsc);
assert okmu and dsc ne 0;
E1, phi := IsogenyFromKernel(E, x - x0);
assert Degree(phi) eq 3;
dp := DefiningPolynomials(phi);
P6b := E1![ Evaluate(d, Eltseq(P6)) : d in dp ];
P3b := E1![ Evaluate(d, Eltseq(P3)) : d in dp ];
assert 6*P6b eq E1!0 and 2*P6b ne E1!0 and 3*P6b ne E1!0 and P3b eq 4*P6b;
print "(0) E : y^2 + 3 lam xy + (lam^3 - 1)y = x^3 with lam = (2t^3 + 1)/(3t^2): P = (0,0) of order 3 and T = (-mm^2, *) of order 2,";
print "    mm = (t^3 - 1)/(3t^2), so P6 = P + T has order 6; the 3-division polynomial vanishes at x0 = -(lam^2 + lam + 1),";
printf "    and -3 (disc of the y-quadratic at x0) = (%o)^2: C = <(x0, *)> is mu_3, E[3] = <P> + mu_3\n", w0;
print "    E1 = E/C, phi : E -> E1 of degree 3, and phi(P6) has order 6";

// ---------------------------------------------------------------- (1) the two elliptic quotients and the gluing
Fr := P![ Coefficient(F, 3 - i) : i in [0..3] ];
Ea, ca := MonEC(F); Eb, cb := MonEC(Fr);
oka, ma := IsIsomorphic(E, Ea);
okb, mb := IsIsomorphic(E1, Eb);
assert oka and okb;
print "(1) phi1 = (x^2, y) : C_t -> Ea : Y^2 = F(X), phi2 = (1/x^2, y/x^3) : C_t -> Eb : Y^2 = X^3 F(1/X);";
print "    Ea is Q(t)-isomorphic to E and Eb to E1, so J is isogenous over Q to E x E1 (phi1^* + phi2^* has kernel of order 4)";
g := MonicCubic(E);
A<th> := quo< P | g >;                        // th = x(T) for the generic point T of order 2 of E
fg := Factorization(g);
lin := [ h[1] : h in fg | Degree(h[1]) eq 1 ]; quad := [ h[1] : h in fg | Degree(h[1]) eq 2 ];
assert #lin eq 1 and #quad eq 1;               // E has exactly one rational point of order 2 over Q(t)
b1 := -Coefficient(lin[1], 0); u := Coefficient(quad[1], 1);
Sw := -u - x + ((2*b1 + u)/Evaluate(quad[1], b1))*quad[1];    // s: fixes x = b1, swaps the two roots of quad
assert Evaluate(Sw, b1) eq b1 and (Sw - (-u - x)) mod quad[1] eq 0 and Evaluate(Sw, Sw) mod g eq x mod g;
al := XImage(phi, Evaluate(Sw, th));          // x(phi(s(T)))
okA, iA := IsInvertible(A!ca); okB, iB := IsInvertible(A!cb); assert okA and okB;
Ra := XImage(ma, th)*iA;                      // the root R of F with (R, 0) = image of T on Y^2 = F(X)
Rb := XImage(mb, al)*iB;                      // the root of Fr with (Rb, 0) = image of phi(s(T)) on Y^2 = Fr(X)
assert Evaluate(F, Ra) eq 0 and Ra*Rb eq 1;
print "    for T in E[2] - 0 the two Weierstrass points x = +-sqrt(R) of C_t over the image (R, 0) of T in Ea lie over the";
print "    image (1/R, 0) of phi(s(T)) in Eb (checked in Q(t)[x]/(2-division polynomial of E)): C_t is the (2,2)-gluing";
print "    of E and E1 along phi o s on E[2], and J = (E x E1)/graph(phi o s on E[2]) (HLP 2000, Prop. 4)";

// ---------------------------------------------------------------- (2) the marked group over Q(t)
J := Jacobian(HyperellipticCurve(f));
Da := Pull1(J, ma(P6), ca); Db := Pull2(J, mb(P3b), cb);
assert 6*Da eq J!0 and 2*Da ne J!0 and 3*Da ne J!0 and Db ne J!0 and 3*Db eq J!0;
els := { i*Da + j*Db : i in [0..5], j in [0..2] };
assert #els eq 18;
print "(2) over Q(t): Da = phi1^*((P+T)') has order 6, Db = phi2^*(phi(P)') has order 3, and the 18 classes";
print "    i Da + j Db are distinct, so G = <Da, Db> = [3,6] is marked over Q(t)";

// ---------------------------------------------------------------- (3) exactness on the residue class t = 234 mod 935
members := [Q| 1/4, 234, -701, 1169, 234 + 935/2];
sanity := [Q| 1/4, 234, -701];
for m in members do for pp in [5, 11, 17] do assert Valuation(m, pp) ge 0; end for; assert GF(5)!m eq GF(5)!234 and GF(11)!m eq GF(11)!234 and GF(17)!m eq GF(17)!234; end for;
p := 5; r0 := 4;
fbar := ReductionClass(f, p, r0);
for m in members do assert PolynomialRing(GF(p))!SpecPoly(f, m) eq fbar; end for;
n, chi := FrobData(fbar);
inv := Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
assert inv eq [6, 6] and n eq 36;
printf "(3) p = %o, t = %o mod %o: reduction %o (squarefree, degree %o); J(F_%o) = %o, chi = %o\n", p, r0, p, fbar, Degree(fbar), p, inv, chi;
print "    so the odd part of J_t(Q)_tors is the marked [3,3] and its 2-part is J_t(Q)[2]";
// the 2-rank: Galois action on the six roots
fq := Factorization(f);
assert #fq eq 2;
s2 := [ h[1] : h in fq | Degree(h[1]) eq 2 ][1]; Q4 := [ h[1] : h in fq | Degree(h[1]) eq 4 ][1];
assert Coefficient(s2, 1) eq 0;
r := -Coefficient(s2, 0)/Coefficient(s2, 2);
F2 := P![ Coefficient(Q4, 2*i) : i in [0..2] ];
assert Evaluate(F2, x^2) eq Q4;
c0, c1, c2 := Explode([ Coefficient(F2, i) : i in [0..2] ]);
okn, nn := IsSquare(c0/c2);
assert okn;
d1 := c1^2 - 4*c0*c2; d2 := -c1/c2 + 2*nn;
assert (d2*(-c1/c2 - 2*nn))*c2^2 eq d1;
// all these specialize on the class: denominators are 5-adic units at t = 4 mod 5, and c0, c2 do not vanish there
for c in [r, Coefficient(s2, 2), c0, c1, c2, nn, d1, d2] do
    Nc, Dc := IntNumDen(c); assert GF(5)!Evaluate(Dc, 4) ne 0;
end for;
for c in [c0, c2] do Nc, Dc := IntNumDen(c); assert GF(5)!Evaluate(Nc, 4) ne 0; end for;
// labels 1..6 = sqrt r, -sqrt r, sqrt X2, -sqrt X2, sqrt X3, -sqrt X3 (sqrt X2 sqrt X3 = nn); v = (e_r, e_d1, e_d2)
function RootPerm(v)
    pr := v[1] eq 1 select [2, 1] else [1, 2];
    q := v[2] eq 0 select (v[3] eq 0 select [3, 4, 5, 6] else [4, 3, 6, 5])
                   else (v[3] eq 0 select [5, 6, 3, 4] else [6, 5, 4, 3]);
    return pr cat q;
end function;
cls := { { S, {1..6} diff S } : S in Subsets({1..6}) | #S in {0, 2, 4, 6} };
assert #cls eq 16;
frob := [];
for pr in [ <5, 4>, <11, 3>, <17, 13> ] do
    pp := pr[1]; rr := pr[2]; Fp := GF(pp);
    v := [];
    for c in [r, d1, d2] do
        Nc, Dc := IntNumDen(c);
        a := Fp!Evaluate(Nc, rr); b := Fp!Evaluate(Dc, rr);
        assert a ne 0 and b ne 0;
        Append(~v, IsSquare(a/b) select 0 else 1);
    end for;
    Append(~frob, v);
    printf "    p = %o, t = %o mod %o: r, d1, d2 are p-adic units with Legendre symbols giving Frobenius = %o\n", pp, rr, pp, v;
end for;
assert Rank(Matrix(GF(2), frob)) eq 3;
function Fixes(v, c)
    pv := RootPerm(v);
    return { Seqset([ pv[i] : i in S ]) : S in c } eq c;
end function;
inv2 := { c : c in cls | &and[ Fixes(v, c) : v in frob ] };
assert #inv2 eq 2 and { Seqset([Integers()|]), {1..6} } in inv2 and { {1, 2}, {3, 4, 5, 6} } in inv2;
print "    the three Frobenius elements fix only the trivial class and [x^2 - r, 0]: J_t(Q)[2] = Z/2";
print "    hence every t = 234 mod 935, integral at 5, 11, 17, has J_t(Q)_tors exactly [3,6]";

// ---------------------------------------------------------------- (4) nonconstant moduli map, (5) sanity
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } eq #members;
printf "(4) the G2-invariants of the members %o are pairwise distinct\n", members;
C0 := HyperellipticCurve(RX![27,-27,56,-13,21,8,6], RX![0,1,1]);
Cm := HyperellipticCurve(SpecPoly(f, 1/4));
assert IsIsomorphic(Cm, SimplifiedModel(C0));
print "    the member t = 1/4 is the census curve 3.6.b (isomorphic over Q)";
for m in sanity do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(5) TorsionSubgroup at t = %o: %o\n", m, invm;
    assert invm eq [3, 6];
end for;
printf "time %o s\n", Cputime(T0);
print "FAMILY_3_6_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E x E1) and J(Q)_tors exactly [3,6]";
quit;
