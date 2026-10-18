/* family_22_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [22].

   The family: the record's "rational-branch subfamilies of the order-11 families"
   (torsion_jac/notes/order22_from_order11.md, code/order22_from_order11_check.m), branch eps = +1 of
   Flynn's order-11 family
       F_t = x^6 + 2x^5 + (2t+3)x^4 + 2x^3 + (t^2+1)x^2 + 2t(1-t)x + t^2,
   on the curve t = (-s^2(s^2+1)(s^4-s^2+1) + 2s^5)/(s^2-1)^2, where F_t(s^2) = 0.

   Argument.
   (1) Over Q(s): the two points at infinity oo1, oo2 of y^2 = F are rational and Dinf = [oo1 - oo2]
       has order 11 (11 Dinf = 0, Dinf != 0).  W = (s^2, 0) is a rational Weierstrass point and
       DW = [W - oo1] satisfies 2 DW = -Dinf, so DW has order 22.
   (2) Residue certificate (family_cert_lib.m) at s = 2 mod 5: #J(F_5) = 22 = #G, strict.  So every
       5-integral s = 2 mod 5 gives J(Q)_tors exactly [22] and a geometrically simple Jacobian.
   (3) The G2-invariants of two members differ, so infinitely many isomorphism classes.
   Run (from the census root):
       magma -b pipeline/magma/family_22_certificate.m > data/knowledge/sources/family_22_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// (1) the family and the marked class over Q(s)
Fs<s> := FunctionField(Q); P<x> := PolynomialRing(Fs);
tF := (-s^2*(s^2 + 1)*(s^4 - s^2 + 1) + 2*s^5)/(s^2 - 1)^2;
f := x^6 + 2*x^5 + (2*tF + 3)*x^4 + 2*x^3 + (tF^2 + 1)*x^2 + 2*tF*(1 - tF)*x + tF^2;
assert Degree(f) eq 6 and Evaluate(f, s^2) eq 0;
C := HyperellipticCurve(f); J := Jacobian(C);
oo := PointsAtInfinity(C);
assert #oo eq 2;
Dinf := J![oo[1], oo[2]];
assert Dinf ne J!0 and 11*Dinf eq J!0;
DW := J![C![s^2, 0], oo[1]];
assert 2*DW eq -Dinf and 22*DW eq J!0 and 11*DW ne J!0;
print "over Q(s): [oo1 - oo2] has order 11, W = (s^2, 0) is a Weierstrass point, and [W - oo1] has order 22";

// (2), (3) residue certificate at s = 2 mod 5
g, st := ResidueCertificate(f, [5], [2], [2, 7, -3, 7/6], 22, [22]);
print "FAMILY_22_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [22]";
quit;
