/* family_28_certificate.m -- certificate that infinitely many pairwise non-isomorphic
   genus-2 curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [28].

   The family is that of V. P. Platonov and G. V. Fedorov, An infinite family of curves of genus 2
   over the field of rational numbers whose Jacobian varieties contain rational points of order 28,
   Dokl. Math. 98 (2018), no. 2, 468-471 (Dokl. Akad. Nauk 482 (2018), no. 4, 385-388),
   formula (8) on p. 470:
       C_t : y^2 = f_t(x) = a5 x^5 + a4 x^4 + a3 x^3 + a2 x^2 + a1 x + (t + 1)^4 (t^2 + t + 1)^2,
       a5 = -4t^5 (t + 1)(2t + 1)(t^2 + t + 1)^2,
       a4 = 4t^12 + 28t^11 + 80t^10 + 136t^9 + 137t^8 + 74t^7 + 11t^6 - 8t^5 + 9t^4 + 20t^3 + 15t^2 + 6t + 1,
       a3 = -4t^12 - 24t^11 - 66t^10 - 96t^9 - 56t^8 + 32t^7 + 70t^6 + 20t^5 - 62t^4 - 90t^3 - 62t^2 - 24t - 4,
       a2 = t^12 + 6t^11 + 13t^10 - 52t^8 - 102t^7 - 76t^6 + 30t^5 + 133t^4 + 150t^3 + 96t^2 + 36t + 6,
       a1 = 2(t + 1)^2 (t^2 + t + 1)(t^6 + 3t^5 + 2t^4 - 4t^3 - 7t^2 - 6t - 2),
   with disc = 256 t^33 (t + 1)^21 (t^2 + t + 1)^9 Q16(t) (p. 470).  It comes from the relations
   (2) w1^2 - f = (x - 1)^5 x and (3) w3^2 - w4^2 f = (x - 1)^2 x^6 (deg w1 = 3, deg w3 = 4, w4 = a(x - x0)),
   i.e. functions with divisors 5P + Q - 6oo and 2P + 6Q - 8oo, where P = (1, *), Q = (0, *);
   then 6(u1) - (u2) = 28(P - oo).  (The degree-4 polynomial printed as "w4" on p. 470 is w3.)

   Argument.
   (0) Over Q(t): (2) and (3) hold for the transcribed w1, w3, w4, and f_t = 4(t+1)^2(t^2+t+1)^4 f.
       D = [x - 1, -t^5 (t + 1)] = P - oo has 28D = 0, 14D != 0, 4D != 0, so D has order 28 and [28]
       is marked over Q(t).  14D = [x - (t+1)^2/(2t+1), 0] is the 2-torsion point of the linear factor
       of f_t; the quartic factor is irreducible over Q(t).
   (1) Residue certificate (family_cert_lib.m) at t = 1 mod 5 and t = 1 mod 7 (t = 1 mod 35):
       #J(F_5) = 28 and #J(F_7) = 56, gcd 28 = #G, strict at 7.  So every 5- and 7-integral
       t = 1 mod 35 gives a geometrically simple Jacobian with J(Q)_tors exactly [28].
       (The prime 5 alone gives exactness for every 5-integral t = 1 mod 5.)
   (2) TorsionSubgroup is [28] at t = 1, 36, -34, and the G2-invariants of the members differ, so
       infinitely many isomorphism classes.
   Run (from the repository root):
       magma -b pipeline/magma/family_28_certificate.m > data/knowledge/sources/family_28_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// (0) the family over Q(t), the relations (2), (3), and the class of order 28
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
a5 := -4*t^5*(t + 1)*(2*t + 1)*(t^2 + t + 1)^2;
a4 := 4*t^12 + 28*t^11 + 80*t^10 + 136*t^9 + 137*t^8 + 74*t^7 + 11*t^6 - 8*t^5 + 9*t^4 + 20*t^3 + 15*t^2 + 6*t + 1;
a3 := -4*t^12 - 24*t^11 - 66*t^10 - 96*t^9 - 56*t^8 + 32*t^7 + 70*t^6 + 20*t^5 - 62*t^4 - 90*t^3 - 62*t^2 - 24*t - 4;
a2 := t^12 + 6*t^11 + 13*t^10 - 52*t^8 - 102*t^7 - 76*t^6 + 30*t^5 + 133*t^4 + 150*t^3 + 96*t^2 + 36*t + 6;
a1 := 2*(t + 1)^2*(t^2 + t + 1)*(t^6 + 3*t^5 + 2*t^4 - 4*t^3 - 7*t^2 - 6*t - 2);
a0 := (t + 1)^4*(t^2 + t + 1)^2;
f := a5*x^5 + a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0;
Q16 := 4*t^16 + 32*t^15 + 179*t^14 + 693*t^13 + 1806*t^12 + 3283*t^11 + 4323*t^10 + 4224*t^9 + 3074*t^8
       + 1615*t^7 + 512*t^6 - 39*t^5 - 181*t^4 - 134*t^3 - 55*t^2 - 12*t - 1;
assert Degree(f) eq 5;
assert Discriminant(f) eq 256*t^33*(t + 1)^21*(t^2 + t + 1)^9*Q16;

x1 := 1/(1 - t^2);
x0 := (t^2 + t + 1)/(-t^4 - t^3 + t^2 + t + 1);
aa := (-t^4 - t^3 + t^2 + t + 1)/t^5;
w4 := aa*(x - x0);
w1 := x^3 - (2*t^6 + 6*t^5 + 15*t^4 + 25*t^3 + 25*t^2 + 15*t + 5)/(2*(t + 1)*(t^2 + t + 1)^2)*x^2
          + (t^6 + 3*t^5 + 10*t^4 + 18*t^3 + 19*t^2 + 12*t + 4)/(2*(t + 1)*(t^2 + t + 1)^2)*x
          - (t + 1)/(2*t^2 + 2*t + 2);
w3 := x^4 - (2*t^10 + 8*t^9 + 15*t^8 + 10*t^7 - 3*t^6 - 8*t^5 - 3*t^4 + 6*t^3 + 7*t^2 + 4*t + 1)/(2*t^5*(t + 1)*(t^2 + t + 1)^2)*x^3
          + (t^4 + t^3 - t^2 - t - 1)*(t^6 + 3*t^5 + 2*t^4 - 6*t^3 - 10*t^2 - 9*t - 3)/(2*t^5*(t + 1)*(t^2 + t + 1)^2)*x^2
          + (2*t^6 + 6*t^5 + 4*t^4 - 6*t^3 - 11*t^2 - 9*t - 3)/(2*t^5*(t + 1)*(t^2 + t + 1))*x
          + (t + 1)/(2*t^5);
f0 := w1^2 - (x - 1)^5*x;                                  // (2), with w2 = 1
assert w3^2 - w4^2*f0 eq (x - 1)^2*x^6;                    // (3)
assert f eq 4*(t + 1)^2*(t^2 + t + 1)^4*f0;               // (8)
print "over Q(t): relations (2), (3) of Platonov-Fedorov hold, f_t = 4(t+1)^2(t^2+t+1)^4 (w1^2 - (x-1)^5 x) is (8), disc = 256 t^33 (t+1)^21 (t^2+t+1)^9 Q16(t)";

J := Jacobian(HyperellipticCurve(f));
D := elt< J | x - 1, -t^5*(t + 1) >;                       // P - oo, P = (1, -t^5 (t+1))
DQ := elt< J | x, -(t + 1)^2*(t^2 + t + 1) >;              // Q - oo, Q = (0, -(t+1)^2 (t^2+t+1))
assert 28*D eq J!0 and 14*D ne J!0 and 4*D ne J!0;
assert DQ eq -5*D;
r := (t + 1)^2/(2*t + 1);
assert Evaluate(f, r) eq 0 and 14*D eq elt< J | x - r, 0 >;
fac := Factorization(f);
assert [ Degree(c[1]) : c in fac ] eq [1, 4];
print "over Q(t): D = [x - 1, -t^5 (t+1)] has order 28 (28D = 0, 14D != 0, 4D != 0); Q - oo = -5D; 14D = [x - (t+1)^2/(2t+1), 0]; f_t = (linear)(irreducible quartic)";

// (1), (2) residue certificate at t = 1 mod 5, t = 1 mod 7 (t = 1 mod 35)
g, st := ResidueCertificate(f, [5, 7], [1, 1], [Q| 1, 36, -34], 28, [28]);
for m in [Q| 36, -34] do
    inv := TorsionInvariants(SpecPoly(f, m));
    printf "sanity: t = %o has TorsionSubgroup %o\n", m, inv;
    assert inv eq [28];
end for;
print "FAMILY_28_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [28]";
quit;
