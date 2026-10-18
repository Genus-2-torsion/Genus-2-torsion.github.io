/* family_2226_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q with J(Q)_tors exactly [2,2,2,6].
   Found 2026-10-18 (torsion_jac notes/f2226_2026_10_18.md); the member s = 4 is the census curve
   2.2.2.6.b (LMFDB 900.b.2).

   Source.  The Howe-Leprevost-Poonen (2,2)-gluing (HLP 2000, Prop. 4) of the Kubert curve
       E2 = E_{2,6}(t) : y^2 + (1 - c)xy - by = x^3 - bx^2,   c = (10 - 2t)/(t^2 - 9),  b = c + c^2,
   (torsion Z/2 x Z/6, P = (0,0) of order 6) with its 3-isogenous curve E1 = E2/<2P>, along the isomorphism
   E1[2] -> E2[2] induced by the dual isogeny.  The Kummer maps commute with the isogeny, so the whole of
   E1[2] glues: the gain is full and J[2] is rational.  With s = t - 3 the gluing is Q(s)-isomorphic to
       C_s : y^2 = f_s(x) = (s^2 + 12)(s^2 - 4x^2)((s - 2)^2 x^2 - (s + 6)^2)((s + 2)^2 x^2 - (s - 6)^2),
   with the six rational Weierstrass points x = +-s/2, +-(s + 6)/(s - 2), +-(s - 6)/(s + 2).  f_s = F(x^2),
   so phi1 : (x,y) -> (x^2, y) maps C_s onto Ea : Y^2 = F(X) and phi2 : (x,y) -> (1/x^2, y/x^3) maps C_s
   onto Eb : Y^2 = X^3 F(1/X); below, Ea is Q(s)-isomorphic to E2 and Eb to E1.  (C_s = C_{-s}.)

   Argument.
   (1) Splitting.  Wherever f_s is squarefree of degree 6, Ea and Eb are elliptic curves and
       phi1^* + phi2^* : Ea x Eb -> J is an isogeny over Q with kernel of order 4.  So every member has J
       split over Q, and J[3] = phi1^*Ea[3] + phi2^*Eb[3] as Galois modules.
   (2) The marked group over Q(s).  The 15 classes [(x - r_i)(x - r_j), 0] of pairs of roots, with 0, are 16
       distinct classes killed by 2 (Mumford arithmetic over Q(s)): J[2] is rational.  D3 = phi1^*(Pa), where
       Pa is the point of order 3 of Ea (the image of 2P), has order 3.  So G = J[2] + <D3> = [2,2,2,6] is
       marked over Q(s), and it embeds in J_s(Q)_tors at every smooth member.
   (3) The Galois module Eb[3] contains mu_3.  The 3-division polynomial of Eb has a linear factor x - xK
       with Y^2 = -3 * (square) at xK: the subgroup {O, (xK, +-yK)} is defined over Q(sqrt -3), and complex
       conjugation inverts it.  (It is the kernel of the dual isogeny Eb = E1 -> E2.)
   (4) Exactness on a residue class.  For every 13-integral s = 4 (mod 13) the curve C_s has the reduction of
       C_4 at 13, with J(F_13) = [2,2,6,6].  Reduction is injective on J_s(Q)_tors (p odd, good reduction).
       - 2-part: J(F_13)[2^inf] = (Z/2)^4 = J[2], so J_s(Q)[2^inf] = (Z/2)^4.
       - No odd primes other than 3 divide #J(F_13) = 144, and J(F_13)[3^inf] = (Z/3)^2 has exponent 3.
       - 3-rank.  J_s(Q)[3] = phi1^*Ea(Q)[3] + phi2^*Eb(Q)[3], and each summand is cyclic (over Q an elliptic
         curve never has (Z/3)^2).  If Eb(Q)[3] were nonzero, then by (3) Eb[3] = Z/3 + mu_3 as a Galois
         module.  Frobenius at 13 = 1 (mod 3) would then fix Eb[3] pointwise and also fix the rational point
         Pa of Ea[3], so J(F_13)[3] would have rank >= 3.  It has rank 2.  So Eb(Q)[3] = 0 and
         J_s(Q)[3^inf] = <D3> = Z/3.
       Hence J_s(Q)_tors is exactly [2,2,2,6].
   (5) Infinitely many.  The G2-invariants of the members differ, so s -> M_2 is nonconstant with finite
       fibres, and the residue class gives infinitely many isomorphism classes.
   (6) Sanity: TorsionSubgroup at s = 4, 17, -9 is [2,2,2,6]; C_4 is Q-isomorphic to 2.2.2.6.b.
   Run (from the repository root):
       magma -b pipeline/magma/family_2226_split_certificate.m > data/knowledge/sources/family_2226_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

Fs<s> := FunctionField(Q); P<x> := PolynomialRing(Fs);
F := (s^2 + 12)*(s^2 - 4*x)*((s - 2)^2*x - (s + 6)^2)*((s + 2)^2*x - (s - 6)^2);
f := Evaluate(F, x^2);
assert Degree(f) eq 6 and IsSquarefree(f);
rts := [ s/2, -s/2, (s + 6)/(s - 2), -(s + 6)/(s - 2), (s - 6)/(s + 2), -(s - 6)/(s + 2) ];
assert f eq LeadingCoefficient(f)*&*[ x - r : r in rts ];
negs := func< c | Evaluate(Numerator(c), -s)/Evaluate(Denominator(c), -s) >;
assert f eq P![ negs(cf) : cf in Coefficients(f) ];
printf "C_s : y^2 = f_s(x) = %o\n", f;
print "f_s = F(x^2), squarefree of degree 6 over Q(s), with the six rational roots +-s/2, +-(s+6)/(s-2), +-(s-6)/(s+2); f_{-s} = f_s";

// ---------------------------------------------------------------- (0) source: the HLP gluing of E2 and E2/<2P>
tt := s + 3;
cK := (10 - 2*tt)/(tt^2 - 9); bK := cK + cK^2;
E2 := EllipticCurve([Fs| 1 - cK, -bK, -bK, 0, 0]);
P6 := E2![0, 0];
assert 6*P6 eq E2!0 and 2*P6 ne E2!0 and 3*P6 ne E2!0;
E1, phi := IsogenyFromKernel(E2, x - (2*P6)[1]);
function MonicCubic(E)
    a1, a2, a3, a4, a6 := Explode(aInvariants(E));
    return x^3 + (a2 + a1^2/4)*x^2 + (a4 + a1*a3/2)*x + (a6 + a3^2/4);
end function;
be := [ -Coefficient(h[1], 0) : h in Factorization(MonicCubic(E2)) ];
assert #be eq 3;
al := [ phi(Points(E2, z)[1])[1] : z in be ];          // x(phi(T_i)): the labels of E1[2] matching E2[2]
assert &and[ Evaluate(MonicCubic(E1), z) eq 0 : z in al ];
a1 := (al[3]-al[2])^2/(be[3]-be[2]) + (al[2]-al[1])^2/(be[2]-be[1]) + (al[1]-al[3])^2/(be[1]-be[3]);
bb1 := (be[3]-be[2])^2/(al[3]-al[2]) + (be[2]-be[1])^2/(al[2]-al[1]) + (be[1]-be[3])^2/(al[1]-al[3]);
a2 := al[1]*(be[3]-be[2]) + al[2]*(be[1]-be[3]) + al[3]*(be[2]-be[1]);
bb2 := be[1]*(al[3]-al[2]) + be[2]*(al[1]-al[3]) + be[3]*(al[2]-al[1]);
Df := ((al[1]-al[2])*(al[1]-al[3])*(al[2]-al[3]))^2;
Dg := ((be[1]-be[2])*(be[1]-be[3])*(be[2]-be[3]))^2;
AA := Dg*a1/a2; BB := Df*bb1/bb2;
hHLP := -(AA*(al[2]-al[1])*(al[1]-al[3])*x^2 + BB*(be[2]-be[1])*(be[1]-be[3]))
        *(AA*(al[3]-al[2])*(al[2]-al[1])*x^2 + BB*(be[3]-be[2])*(be[2]-be[1]))
        *(AA*(al[1]-al[3])*(al[3]-al[2])*x^2 + BB*(be[1]-be[3])*(be[3]-be[2]));
ratio := LeadingCoefficient(hHLP)/LeadingCoefficient(f);
assert hHLP eq ratio*f;
okr, sqr := IsSquare(ratio);
assert okr;
printf "(0) the HLP (Prop. 4) gluing of E2 = E_{2,6}(s+3) and E1 = E2/<2P> along the dual isogeny on 2-torsion is (%o)^2 f_s\n", sqr;

// ---------------------------------------------------------------- (1) the two elliptic quotients
function MonEC(G)      // Y^2 = G(X) cubic -> monic model V^2 = U^3 + c2 U^2 + c1 c3 U + c0 c3^2 (U = c3 X, V = c3 Y)
    cs := [ Coefficient(G, i) : i in [0..3] ];
    return EllipticCurve([0, cs[3], 0, cs[2]*cs[4], cs[1]*cs[4]^2]), cs[4];
end function;
Fr := P![ Coefficient(F, 3 - i) : i in [0..3] ];        // X^3 F(1/X)
Ea, ca := MonEC(F); Eb, cb := MonEC(Fr);
oka, ma := IsIsomorphic(E2, Ea);
okb, mb := IsIsomorphic(E1, Eb);
assert oka and okb;
print "(1) phi1 = (x^2, y) : C_s -> Ea : Y^2 = F(X) and phi2 = (1/x^2, y/x^3) : C_s -> Eb : Y^2 = X^3 F(1/X);";
print "    Ea is Q(s)-isomorphic to E2 = E_{2,6}(s+3) and Eb to E1 = E2/<2P>: J is isogenous over Q to E1 x E2";

// ---------------------------------------------------------------- (2) the marked group over Q(s)
J := Jacobian(HyperellipticCurve(f));
W := [ elt< J | (x - rts[1])*(x - rts[j]), 0 > : j in [2..5] ];
assert &and[ 2*w eq J!0 and w ne J!0 : w in W ];
J2 := { i1*W[1] + i2*W[2] + i3*W[3] + i4*W[4] : i1, i2, i3, i4 in [0..1] };
assert #J2 eq 16;
pts := [ <i, j> : i, j in [1..6] | i lt j ];
assert { elt< J | (x - rts[pr[1]])*(x - rts[pr[2]]), 0 > : pr in pts } join { J!0 } eq J2;
Pa := ma(2*P6);                                          // the point of order 3 on Ea (monic model)
Xa := Pa[1]/ca; Ya := Pa[2]/ca;                          // on Y^2 = F(X)
assert Evaluate(F, Xa) eq Ya^2;
D3 := elt< J | x^2 - Xa, Ya >;
assert D3 ne J!0 and 3*D3 eq J!0;
printf "(2) over Q(s): the 15 classes [(x - r_i)(x - r_j), 0] and 0 form J[2] = (Z/2)^4; D3 = phi1^*(Pa) = %o has order 3\n", D3;
print "    so G = J[2] + <D3> = [2,2,2,6] is marked over Q(s)";

// ---------------------------------------------------------------- (3) mu_3 inside Eb[3]
gb := MonicCubic(Eb);
lin := [ h[1] : h in Factorization(DivisionPolynomial(Eb, 3)) | Degree(h[1]) eq 1 ];
assert #lin ge 1;
xK := -Coefficient(lin[1], 0)/Coefficient(lin[1], 1);
yK2 := Evaluate(gb, xK);
okm, wK := IsSquare(-3*yK2);
assert okm and not IsSquare(yK2);
printf "(3) Eb: the 3-division polynomial has the root xK = %o, where yK^2 = -3 (%o)^2: {O, (xK, +-yK)} is defined over Q(sqrt -3) and inverted by conjugation, so it is mu_3\n", xK, wK/3;
assert #[ h : h in Factorization(DivisionPolynomial(Ea, 3)) | Degree(h[1]) eq 1 and IsSquare(Evaluate(MonicCubic(Ea), -Coefficient(h[1],0)/Coefficient(h[1],1))) ] ge 1;
print "    Ea has the rational point of order 3 (Pa)";

// ---------------------------------------------------------------- (4) exactness on the residue class s = 4 mod 13
p := 13; r0 := 4;
assert p mod 3 eq 1;
members := [Q| 4, 17, -9, 17/14];
fbar := ReductionClass(f, p, r0);
for m in members do
    assert GF(p)!m eq GF(p)!r0;
    assert PolynomialRing(GF(p))!SpecPoly(f, m) eq fbar;
end for;
n, chi := FrobData(fbar);
Ap := AbelianGroup(Jacobian(HyperellipticCurve(fbar)));
inv := Invariants(Ap);
assert #Ap eq n and inv eq [2, 2, 6, 6];
printf "(4) p = %o, s = %o mod %o: reduction %o (squarefree, degree %o); J(F_13) = %o (order %o), chi = %o\n",
    p, r0, p, fbar, Degree(fbar), inv, n, chi;
NwK, DwK := IntNumDen(wK);
assert GF(p)!Evaluate(NwK, r0) ne 0 and GF(p)!Evaluate(DwK, r0) ne 0;
print "    wK(4) is a unit mod 13, so at every member of the class the mu_3 of (3) specializes to a subgroup of order 3 of Eb";
print "    2-part (Z/2)^4 = J[2]; odd part (Z/3)^2 of exponent 3; 13 = 1 mod 3, so Eb(Q)[3] != 0 would force 3-rank >= 3 in J(F_13)";
print "    hence every 13-integral s = 4 mod 13 has J_s(Q)_tors exactly [2,2,2,6]";

// ---------------------------------------------------------------- (5) nonconstant moduli map, (6) sanity
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } ge 2;
printf "(5) the G2-invariants of the members %o are not all equal\n", members;
for m in [Q| 4, 17, -9] do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(6) TorsionSubgroup at s = %o: %o\n", m, invm;
    assert invm eq [2, 2, 2, 6];
end for;
Ccen := HyperellipticCurve(RX![0, 30, 39, -72, -111, -30], RX![0, 1, 1]);
assert IsIsomorphic(HyperellipticCurve(SpecPoly(f, 4)), Ccen);
print "(6) C_4 is Q-isomorphic to the census curve 2.2.2.6.b (LMFDB 900.b.2)";
printf "time %o s\n", Cputime(T0);
print "FAMILY_2226_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E1 x E2) and J(Q)_tors exactly [2,2,2,6]";
quit;
