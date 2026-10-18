/* family_210_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E1 x E2) with J(Q)_tors exactly [2,10].
   Found 2026-10-18 (torsion_jac notes/g33fam_2026_10_18.md).

   Source.  The Howe-Leprevost-Poonen (2,2)-gluing (HLP 2000, Prop. 4) of the Kubert X_1(10) curve
       E2 = E_10(t) : y^2 + (1 - c)xy - by = x^3 - bx^2,  b = t^3(t-1)(2t-1)/(t^2-3t+1)^2,  c = -t(t-1)(2t-1)/(t^2-3t+1),
   (torsion Z/10, P = (0,0) of order 10) with its 5-isogenous curve E1 = E2/<2P>, along the isomorphism
   E2[2] -> E1[2] induced by the 5-isogeny phi (the route of the [2,2,2,6] certificate, with 5 in place of 3).
   E2 has one rational 2-torsion point, so the HLP formula is evaluated over K = Q(t)(sqrt delta), delta the
   discriminant of the quadratic factor of the 2-division polynomial; the result has coefficients in Q(t).
   Up to a square factor it is
       C_t : y^2 = f_t(x) = F(x^2)     (F a cubic over Q(t), printed below),
   and f_t = lc (x - r)(x + r) q(x) q(-x) over Q(t), with r = (2t^2 - t - 1/2)/(t - 1/2) and q a quadratic.
   phi1 : (x,y) -> (x^2, y) maps C_t onto Ea : Y^2 = F(X), phi2 : (x,y) -> (1/x^2, y/x^3) onto Eb : Y^2 = X^3 F(1/X);
   below, Ea is Q(t)-isomorphic to E2 and Eb to E1.

   Argument.
   (1) Splitting.  Wherever f_t is squarefree of degree 6, phi1^* + phi2^* : Ea x Eb -> J is an isogeny over Q
       with kernel of order 4, so J is split over Q and J(Q)[l^inf] = Ea(Q)[l^inf] x Eb(Q)[l^inf] for odd l.
   (2) The marked group over Q(t).  c1 = [(x - r)(x + r), 0] and c2 = [q(x), 0] are distinct nonzero classes
       killed by 2 (Mumford arithmetic over Q(t)); c2 is the gain of the gluing (the Kummer maps commute with phi).
       D5 = phi1^*(Pa), Pa the image of 2P on Ea, has order 5.  So G = <c1, c2, D5> = [2,10] is marked over Q(t)
       and embeds in J_t(Q)_tors at every smooth member.
   (3) mu_5 inside E1[5].  The kernel of the dual isogeny E1 -> E2 is phi(E2[5]) = E2[5]/<2P>, which the Weil
       pairing identifies with Hom(<2P>, mu_5) = mu_5.  So if E1(Q)[5] were nonzero, E1[5] = Z/5 + mu_5 as a
       Galois module (mu_5 has no nonzero rational point), and Frobenius at a good prime p = 1 (mod 5) would act
       trivially on E1[5]; with the rational point 2P of E2 this gives J(F_p)[5] of rank >= 3.
   (4) Exactness on a residue class.  For every 7- and 11-integral t = 58 (mod 77) the curve C_t has the
       reduction of C_58 at 7 and at 11, with J(F_7) = J(F_11) = [10,10].  Reduction is injective on J_t(Q)_tors.
       - J(F_7) = [10,10]: the 2-part of J_t(Q)_tors has rank <= 2 and exponent 2, so it is <c1, c2>; the only
         odd prime is 5, with 5-part of exponent 5.
       - J(Q)[5^inf] = E1(Q)[5^inf] x E2(Q)[5^inf] by (1), and E2(Q)[5^inf] = Z/5 (Mazur).  At p = 11 = 1 (mod 5)
         J(F_11)[5] has rank 2, so by (3) E1(Q)[5] = 0, and J_t(Q)[5^inf] = <D5> = Z/5.
       Hence J_t(Q)_tors is exactly [2,10].  (#J(F_p) = #E2(F_p)^2 is a square, so a gcd of point counts cannot
       certify this; the group structures are used instead.)
   (5) Infinitely many.  The G2-invariants of the members differ, so t -> M_2 is nonconstant with finite
       fibres, and the residue class gives infinitely many isomorphism classes.
   (6) Sanity: TorsionSubgroup at t = 58, -19, 135 is [2,10].
   Run (from the repository root):
       magma -b pipeline/magma/family_210_split_certificate.m > data/knowledge/sources/family_210_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
f := (-4194304*t^16 + 25165824*t^15 - 66060288*t^14 + 103809024*t^13 - 122159104*t^12 + 137363456*t^11 - 155713536*t^10 + 149946368*t^9 - 106692608*t^8 + 52428800*t^7 - 16777216*t^6 + 3145728*t^5 - 262144*t^4)*x^6 + (16777216*t^18 - 100663296*t^17 + 264241152*t^16 - 318767104*t^15 + 59768832*t^14 + 355467264*t^13 - 669515776*t^12 + 959447040*t^11 - 1363673088*t^10 + 1691353088*t^9 - 1633419264*t^8 + 1163919360*t^7 - 566231040*t^6 + 141557760*t^5 + 25952256*t^4 - 38797312*t^3 + 15728640*t^2 - 3145728*t + 262144)*x^4 + (-33554432*t^18 - 201326592*t^17 + 1296039936*t^16 - 2860515328*t^15 + 3601858560*t^14 - 3136290816*t^13 + 2208825344*t^12 - 1223688192*t^11 - 379060224*t^10 + 2962227200*t^9 - 5397282816*t^8 + 5970591744*t^7 - 4479516672*t^6 + 2349858816*t^5 - 869007360*t^4 + 228589568*t^3 - 44040192*t^2 + 6291456*t - 524288)*x^2 + 16777216*t^18 - 234881024*t^17 + 1396703232*t^16 - 4630511616*t^15 + 9557770240*t^14 - 13429112832*t^13 + 14890303488*t^12 - 15979249664*t^11 + 17065050112*t^10 - 14936965120*t^9 + 8815116288*t^8 - 2623537152*t^7 - 339738624*t^6 + 625999872*t^5 - 213647360*t^4 + 11534336*t^3 + 11534336*t^2 - 3145728*t + 262144;
cz := GCD([Integers()!c : c in &cat[Coefficients(Numerator(cf)) : cf in Coefficients(f)]]);
sq := &*([Integers()| 1] cat [pp[1]^(pp[2] div 2) : pp in Factorization(cz)]);
f := f/sq^2;
assert &and[ Coefficient(f, i) eq 0 : i in [1, 3, 5] ] and Degree(f) eq 6 and IsSquarefree(f);
F := P![ Coefficient(f, 2*i) : i in [0..3] ];            // f = F(x^2)
r := (2*t^2 - t - 1/2)/(t - 1/2);
fa := Factorization(f);
assert #fa eq 4 and { Degree(g[1]) : g in fa } eq {1, 2};
qs := [ g[1] : g in fa | Degree(g[1]) eq 2 ];
assert (x - r)*(x + r)*qs[1]*qs[2]*LeadingCoefficient(f) eq f;
printf "C_t : y^2 = f_t(x) = F(x^2), F = %o\n", F;
print "f_t = lc (x - r)(x + r) q(x) q(-x) over Q(t), r = (2t^2 - t - 1/2)/(t - 1/2)";

// ---------------------------------------------------------------- (0) source: the HLP gluing of E2 = E_10(t) and E2/<2P>
bK := t^3*(t-1)*(2*t-1)/(t^2-3*t+1)^2; cK := -t*(t-1)*(2*t-1)/(t^2-3*t+1);
E2 := EllipticCurve([Ft| 1-cK, -bK, -bK, 0, 0]);
P10 := E2![0,0];
assert 10*P10 eq E2!0 and 5*P10 ne E2!0 and 2*P10 ne E2!0;
E1, phi := IsogenyFromKernel(E2, (x - (2*P10)[1])*(x - (4*P10)[1]));
assert Degree(phi) eq 5;
function MonicCubic(E)
    a1, a2, a3, a4, a6 := Explode(aInvariants(E));
    return x^3 + (a2 + a1^2/4)*x^2 + (a4 + a1*a3/2)*x + (a6 + a3^2/4);
end function;
g2 := MonicCubic(E2); g1 := MonicCubic(E1);
fa2 := Factorization(g2);
lin := [ g[1] : g in fa2 | Degree(g[1]) eq 1 ][1]; quad := [ g[1] : g in fa2 | Degree(g[1]) eq 2 ][1];
K<z> := ext< Ft | quad >; PK<X> := PolynomialRing(K);
be := [ K!(-Coefficient(lin, 0)), z, -Coefficient(quad, 1) - z ];      // roots: the x-coordinates of E2[2]
dp := DefiningPolynomials(phi);
a1, a2, a3, a4, a6 := Explode(aInvariants(E2));
function PhiX(xT) yT := -(a1*xT + a3)/2; vals := [ Evaluate(d, [xT, yT, 1]) : d in dp ]; return vals[1]/vals[3]; end function;
al := [ PhiX(b) : b in be ];                                           // x(phi(T_i)): the matching labels of E1[2]
assert &and[ Evaluate(PolynomialRing(K)!g1, a) eq 0 : a in al ] and #Seqset(al) eq 3;
aa1 := (al[3]-al[2])^2/(be[3]-be[2]) + (al[2]-al[1])^2/(be[2]-be[1]) + (al[1]-al[3])^2/(be[1]-be[3]);
bb1 := (be[3]-be[2])^2/(al[3]-al[2]) + (be[2]-be[1])^2/(al[2]-al[1]) + (be[1]-be[3])^2/(al[1]-al[3]);
aa2 := al[1]*(be[3]-be[2]) + al[2]*(be[1]-be[3]) + al[3]*(be[2]-be[1]);
bb2 := be[1]*(al[3]-al[2]) + be[2]*(al[1]-al[3]) + be[3]*(al[2]-al[1]);
Df := ((al[1]-al[2])*(al[1]-al[3])*(al[2]-al[3]))^2;
Dg := ((be[1]-be[2])*(be[1]-be[3])*(be[2]-be[3]))^2;
AA := Dg*aa1/aa2; BB := Df*bb1/bb2;
hHLP := -(AA*(al[2]-al[1])*(al[1]-al[3])*X^2 + BB*(be[2]-be[1])*(be[1]-be[3]))
        *(AA*(al[3]-al[2])*(al[2]-al[1])*X^2 + BB*(be[3]-be[2])*(be[2]-be[1]))
        *(AA*(al[1]-al[3])*(al[3]-al[2])*X^2 + BB*(be[1]-be[3])*(be[3]-be[2]));
hq := P![ Ft!c : c in Coefficients(hHLP) ];
ratio := LeadingCoefficient(hq)/LeadingCoefficient(f);
assert hq eq ratio*f;
okr, sqr := IsSquare(ratio);
assert okr;
print "(0) the HLP (Prop. 4) gluing of E2 = E_10(t) and E1 = E2/<2P> along phi on 2-torsion has coefficients in Q(t) and is (sqrt ratio)^2 f_t";

// ---------------------------------------------------------------- (1) the two elliptic quotients
function MonEC(G)
    cs := [ Coefficient(G, i) : i in [0..3] ];
    return EllipticCurve([0, cs[3], 0, cs[2]*cs[4], cs[1]*cs[4]^2]), cs[4];
end function;
Fr := P![ Coefficient(F, 3 - i) : i in [0..3] ];
Ea, ca := MonEC(F); Eb, cb := MonEC(Fr);
oka, ma := IsIsomorphic(E2, Ea);
okb, mb := IsIsomorphic(E1, Eb);
assert oka and okb;
print "(1) phi1 = (x^2, y) : C_t -> Ea : Y^2 = F(X), phi2 = (1/x^2, y/x^3) : C_t -> Eb : Y^2 = X^3 F(1/X);";
print "    Ea is Q(t)-isomorphic to E2 = E_10(t) and Eb to E1 = E2/<2P>: J is isogenous over Q to E1 x E2";

// ---------------------------------------------------------------- (2) the marked group over Q(t)
J := Jacobian(HyperellipticCurve(f));
c1 := elt< J | (x - r)*(x + r), 0 >; c2 := elt< J | qs[1], 0 >;
assert c1 ne J!0 and c2 ne J!0 and c1 ne c2 and 2*c1 eq J!0 and 2*c2 eq J!0;
Pa := ma(2*P10);
Xa := Pa[1]/ca; Ya := Pa[2]/ca;
assert Evaluate(F, Xa) eq Ya^2;
D5 := elt< J | x^2 - Xa, Ya >;
assert D5 ne J!0 and 5*D5 eq J!0;
printf "(2) over Q(t): c1 = [(x-r)(x+r), 0], c2 = [q, 0] are distinct nonzero classes killed by 2; D5 = phi1^*(Pa) has order 5\n";
print "    so G = <c1, c2, D5> = [2,10] is marked over Q(t)";

// ---------------------------------------------------------------- (4) exactness on the residue class t = 58 mod 77
primes := [7, 11]; residues := [58 mod 7, 58 mod 11];
members := [Q| 58, -19, 135, 58 + 77*13/2];
for i in [1..2] do
    p := primes[i]; r0 := residues[i];
    fbar := ReductionClass(f, p, r0);
    for m in members do
        assert GF(p)!m eq GF(p)!r0;
        assert PolynomialRing(GF(p))!SpecPoly(f, m) eq fbar;
    end for;
    n, chi := FrobData(fbar);
    inv := Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
    assert inv eq [10, 10] and n eq 100;
    printf "(4) p = %o, t = %o mod %o: reduction %o (squarefree, degree %o); J(F_%o) = %o, chi = %o\n", p, r0, p, fbar, Degree(fbar), p, inv, chi;
end for;
assert 11 mod 5 eq 1;
print "    p = 7: the 2-part of J_t(Q)_tors is <c1, c2> and the odd part is a 5-group of exponent 5";
print "    p = 11 = 1 mod 5: J(F_11)[5] has rank 2, so E1(Q)[5] = 0 by (3), and J_t(Q)[5^inf] = <D5>";
print "    hence every 7- and 11-integral t = 58 mod 77 has J_t(Q)_tors exactly [2,10]";

// ---------------------------------------------------------------- (5) nonconstant moduli map, (6) sanity
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } eq #members;
printf "(5) the G2-invariants of the members %o are pairwise distinct\n", members;
for m in [Q| 58, -19, 135] do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(6) TorsionSubgroup at t = %o: %o\n", m, invm;
    assert invm eq [2, 10];
end for;
printf "time %o s\n", Cputime(T0);
print "FAMILY_210_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E1 x E2) and J(Q)_tors exactly [2,10]";
quit;
