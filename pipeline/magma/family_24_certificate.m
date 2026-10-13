/* family_24_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [24].

   Relation-lattice construction (method of torsion_jac/notes/z38_2026_10_13.md section 2; found in the
   notes review of 2026-10-13).  Quintic model y^2 = f, W = infinity, P_i = (i, A(i)), D_i = P_i - W, with
       B - A = (x-1) g1,  B + A = g2,  g1 g2 = lambda x^4 + (x-1)^5,  g1 quadratic,
       f = A^2 + lambda x^4 (x-1) = B^2 - (x-1)^6,  B(1) = A(1),
   so div(y - A) = 4P0 + P1 - 5W and div(y - B) = 6P1 - 6W: lattice (4,1), (0,6), determinant 24, cyclic.
   The quadratics g1 = x^2 - e1 x + e2 dividing some lambda x^4 + (x-1)^5 form a rational curve; a
   parametrization is e1 = (-33t^4 + 65t^3 - 42t^2 + 11t - 1)/t^4, e2 = (-9t^3 + 15t^2 - 7t + 1)/t^3
   (checked below: the remainder of lambda x^4 + (x-1)^5 mod g1 vanishes for one lambda in Q(t)).
   Certificate: t = 4 mod 5 gives #J(F_5) = 24; t = 2 mod 7 gives #J(F_7) = 48 with a strict Frobenius;
   so every t = 9 mod 35 (5- and 7-integral) has torsion exactly [24] and a geometrically simple
   Jacobian (family_cert_lib.m).
   Run:  magma -b pipeline/magma/family_24_certificate.m > data/knowledge/sources/family_24_certificate.log */
SetColumns(0);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
e1 := (-33*t^4 + 65*t^3 - 42*t^2 + 11*t - 1)/t^4; e2 := (-9*t^3 + 15*t^2 - 7*t + 1)/t^3;
g1 := x^2 - e1*x + e2;
r4 := x^4 mod g1; r5 := (x-1)^5 mod g1;
lam := -Coefficient(r5, 1)/Coefficient(r4, 1);
assert lam*r4 + r5 eq 0;
Rx := lam*x^4 + (x-1)^5; g2 := Rx div g1; assert g1*g2 eq Rx;
BmA := (x-1)*g1; BpA := g2;
A := (BpA - BmA)/2; B := (BpA + BmA)/2;
f := A^2 + lam*x^4*(x-1);
assert f eq B^2 - (x-1)^6 and Degree(A) le 2 and Degree(f) eq 5 and IsMonic(B) and Degree(B) eq 3;
assert Evaluate(B, 1) eq Evaluate(A, 1);
printf "lambda = %o\n", lam;
J := Jacobian(HyperellipticCurve(f));
D0 := elt< J | x, Evaluate(A, 0) >; D1 := elt< J | x - 1, Evaluate(A, 1) >;
assert 4*D0 + D1 eq J!0 and 6*D1 eq J!0;
assert 24*D0 eq J!0 and 12*D0 ne J!0 and 8*D0 ne J!0;
print "over Q(t): 4 D0 + D1 = 6 D1 = 0 and D0 has order 24";
g, st := ResidueCertificate(f, [5, 7], [4, 2], [9, 44, -26, 79], 24, [24]);
print "FAMILY_24_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [24]";
quit;
