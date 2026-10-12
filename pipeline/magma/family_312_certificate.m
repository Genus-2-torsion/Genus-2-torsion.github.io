/* family_312_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [3,12].

   The family (found 2026-10-12, torsion_jac/notes/m312_a12_sieve_2026_10_12.md; census curve 3.12.b
   is its member z = 5).  On the line (A,B) = (-2/3,-10) of the order-12 chart A(12), written in the
   coordinates (A,z,B) = (p^2, z, p*r), the curve is C_z : y^2 = f_z with
       s = (z^2-4A+1)/(2z), t = (z^2+4A-1)^2/(8Az), M = ((s^2-1)(2B+1) - A(2st-4))/(4A),
       lambda = (4A-M^2)/(s^2-1), R = (A x^2 + (2B+1)x + B^2/A - 1)/lambda,
       f_z = R (R x^2 + 4(R+x-1)(R-1)),
   with the marked class P12 = [u4, R ell mod u4] + [u6, x R mod u6] of order 12 (u4 ~ 2(Ax+B)+M R,
   u6 ~ R+x-1, ell = s x + t).  On the square multiple lambda^4 f_z (an isomorphic model) the class
   Q3 = [u, G], u = x^2 + a x + b, satisfies lambda^4 f_z - l u^3 = G^2, where a, b, l are rational in
   lambda and G is defined over Q(E), E : w^2 = (3z^2+25)(3z^2+121).  Then div(y - G) = 3D - 3(oo+ + oo-),
   so Q3 has order 3 over Q(E).

   Argument.
   (1) Over Q(E): P12 has order 12 and Q3 has order 3. At z = 5 (w = 140) they are independent, so
       Z/12 x Z/3 is contained in J(Q(E)).
   (2) E has conductor 990 and rank 1; let P generate E(Q) mod torsion.  P has order 10 in E(F_37)
       and order 22 in E(F_41).  The quartic model w^2 = (3z^2+25)(3z^2+121) is smooth over F_37 and
       F_41, with the integral base point (5,140).  So for R = 110m P (m >= 1) the point of the quartic
       reduces to (5,140), i.e. z(R) = 5 mod 37 and mod 41.
   (3) The coefficients of f_z have denominators that are units at z = 5 mod 37 and mod 41, and
       f_5 is squarefree of degree 6 mod 37 and mod 41.  So C_{z(R)} has the same good reduction as
       C_5 = 3.12.b at 37 and 41: #J(F_37) = 1332, #J(F_41) = 1980, gcd 36 = #[3,12], and the
       Frobenius polynomial at 37 is strict.
   (4) Specialization of torsion sections is injective on smooth fibers, so J(Q) of C_{z(R)} contains
       Z/12 x Z/3. Its order divides 36 by (3), hence it equals [3,12]. J is geometrically simple
       (strict prime 37).
   (5) The G2-invariants at z = 5 and z = 13/21 differ, so the moduli map E -> M_2 is nonconstant
       with finite fibers: infinitely many of these curves are pairwise non-isomorphic.
   Run:  magma -b pipeline/magma/family_312_certificate.m > data/knowledge/sources/family_312_certificate.log
*/
SetColumns(0);
Q := Rationals(); Z := Integers();
Fz<z> := FunctionField(Q); P<x> := PolynomialRing(Fz);
A := Fz!(-2/3); B := Fz!(-10);
s := (z^2-4*A+1)/(2*z); t := (z^2+4*A-1)^2/(8*A*z);
M := ((s^2-1)*(2*B+1) - A*(2*s*t-4))/(4*A);
lam := (4*A - M^2)/(s^2-1);
Rr := (A*x^2 + (2*B+1)*x + B^2/A - 1)/lam;
f := Rr*(Rr*x^2 + 4*(Rr+x-1)*(Rr-1));
assert Degree(f) eq 6;
assert lam eq -(81*z^8 + 5292*z^6 + 98550*z^4 + 439500*z^2 + 1890625)/(9216*z^2);
assert 125/4 - lam eq (3*z^2+25)^3*(3*z^2+121)/(96*z)^2;

// (1a) the 3-torsion section Q3 over Q(E)
a := (18*lam + 855)/(lam + 30); b := (45*lam + 6795)/(lam + 30);
l := (lam^4/729 + 10*lam^3/81 + 100*lam^2/27 + 1000*lam/27)/(lam - 125/4);
u := x^2 + a*x + b;
fl := lam^4*f;
Pp := fl - l*u^3; c := LeadingCoefficient(Pp);
assert Degree(Pp) eq 6 and l ne 0 and Discriminant(u) ne 0;
okH, H := IsSquare(Pp/c);
assert okH and Degree(H) eq 3;
assert c eq lam^2*(lam-117)^2/(729*(125/4-lam));
assert IsSquare(c*(3*z^2+25)*(3*z^2+121));        // so G = sqrt(c) H lies in Q(E)[x]
print "over Q(E): lambda^4 f - l u^3 = G^2 with G a cubic in Q(E)[x], so Q3 = [u, G] has order 3";

// (1b) P12 of order 12 over Q(z)
ell := s*x + t;
u4 := 2*(A*x+B) + M*Rr; u4 := u4/LeadingCoefficient(u4);
u6 := Rr + x - 1; u6 := u6/LeadingCoefficient(u6);
J := Jacobian(HyperellipticCurve(f));
P12 := elt< J | u4, (Rr*ell) mod u4 > + elt< J | u6, (x*Rr) mod u6 >;
assert 12*P12 eq J!0 and 6*P12 ne J!0 and 4*P12 ne J!0;
print "over Q(z): P12 has order exactly 12";

// (1c) independence at z = 5, w = 140; the fiber is the census curve 3.12.b
R<X> := PolynomialRing(Q);
spec := func< g, z0 | R![ Evaluate(cc, z0) : cc in Coefficients(g) ] >;
c5 := Evaluate(c, 5); d5 := Evaluate((3*z^2+25)*(3*z^2+121), 5);
ok, w5 := IsSquare(d5); assert ok and w5 eq 140;
ok, sq := IsSquare(c5*d5); assert ok; sq := sq/w5;
f5 := spec(fl, 5); u5 := spec(u, 5); G5 := sq*spec(H, 5);
assert f5 - Evaluate(l, 5)*u5^3 eq G5^2;
J5 := Jacobian(HyperellipticCurve(f5));
Q3 := elt< J5 | u5, G5 mod u5 >;
L2 := Evaluate(lam, 5)^2;
P12_5 := elt< J5 | spec(u4,5), (L2*spec((Rr*ell) mod u4, 5)) mod spec(u4,5) >
       + elt< J5 | spec(u6,5), (L2*spec((x*Rr) mod u6, 5)) mod spec(u6,5) >;
assert Q3 ne J5!0 and 3*Q3 eq J5!0 and Order(P12_5) eq 12;
assert Q3 notin {J5!0, 4*P12_5, -4*P12_5};
C312b := HyperellipticCurve(229392*X^6 + 475068*X^5 - 1261612*X^4 - 7653746*X^3 + 231383*X^2 + 10277505*X + 19644345, X^2 + X);
assert IsIsomorphic(HyperellipticCurve(spec(f, 5)), C312b);
print "z = 5: <P12, Q3> = Z/12 x Z/3; the fiber is census curve 3.12.b";

// (2) the elliptic curve E and the points 110m P
P2<Zz,W,V> := ProjectiveSpace(Q, [1,2,1]);
Hq := Curve(P2, W^2 - (3*Zz^2+25*V^2)*(3*Zz^2+121*V^2));
E0, phi := EllipticCurve(Hq, Hq![5,140,1]);
E := MinimalModel(E0);
assert aInvariants(E) eq [1,-1,0,-675,-6075] and Conductor(E) eq 990;
r1, r2 := RankBounds(E); assert r1 eq 1 and r2 eq 1;
G, mG := MordellWeilGroup(E);
Pgen := [ mG(g) : g in OrderedGenerators(G) | Order(g) eq 0 ][1];
quartic := (3*X^2+25)*(3*X^2+121);
for p in [37, 41] do
    assert IsSquarefree(PolynomialRing(GF(p))!quartic) and IsSquare(GF(p)!9);
    assert 110 mod Order(ChangeRing(E, GF(p))!Pgen) eq 0;
end for;
printf "E = %o, conductor 990, rank 1, P = %o; ord(P mod 37) = %o, ord(P mod 41) = %o\n",
    aInvariants(E), Pgen, Order(ChangeRing(E, GF(37))!Pgen), Order(ChangeRing(E, GF(41))!Pgen);
ok, iso := IsIsomorphic(E, E0); assert ok;
zs := [ q[1]/q[3] : q in RationalPoints(Pullback(phi, iso(110*Pgen))) | q[3] ne 0 ];
zbig := [ zz : zz in zs | Max(Abs(Numerator(zz)), Denominator(zz)) gt 10^100 ];
assert #zbig ge 1 and &and[ GF(37)!zz eq 5 and GF(41)!zz eq 5 : zz in zbig ];
printf "sanity: z(110 P) has %o digits and is = 5 mod 37 and mod 41\n", #Sprint(Numerator(zbig[1]));

// (3) same good reduction as z = 5 at 37 and 41
Tt := PolynomialRing(Q).1;
nJ := [];
for p in [37, 41] do
    for cc in Coefficients(f) do
        dd := Denominator(cc); nn := Numerator(cc);
        k := LCM([Denominator(e) : e in Coefficients(nn)] cat [Denominator(e) : e in Coefficients(dd)]);
        assert Evaluate(PolynomialRing(Z)!(k*dd), 5) mod p ne 0;
    end for;
    fp := PolynomialRing(GF(p))![ GF(p)!Evaluate(cc, 5) : cc in Coefficients(f) ];
    assert Degree(fp) eq 6 and IsSquarefree(fp);
    Lp := LPolynomial(HyperellipticCurve(fp) : Al := "Naive");
    chi := Tt^4*Evaluate(PolynomialRing(Q)!Lp, 1/Tt); chi := PolynomialRing(Q)!chi;
    Append(~nJ, Evaluate(Lp, 1));
    strict := IsIrreducible(chi) and &and[ Degree(MinimalPolynomial(NumberField(chi).1^n)) eq 4 : n in [2..12] ];
    printf "p = %o: #J(F_p) = %o, chi = %o%o\n", p, Evaluate(Lp, 1), chi, strict select " (strict)" else "";
    if p eq 37 then assert strict; end if;
end for;
assert GCD(nJ) eq 36;

// (5) the moduli map is nonconstant
assert G2Invariants(HyperellipticCurve(spec(f, 5))) ne G2Invariants(HyperellipticCurve(spec(f, 13/21)));
print "gcd(#J(F_37), #J(F_41)) = 36 = #[3,12]; G2-invariants at z = 5 and z = 13/21 differ";
print "FAMILY_312_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [3,12]";
quit;
