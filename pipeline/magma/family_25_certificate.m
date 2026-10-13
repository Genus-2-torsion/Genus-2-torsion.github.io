/* family_25_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [25].

   Relation-lattice construction (method of torsion_jac/notes/z38_2026_10_13.md section 2; found in the
   notes review of 2026-10-13).  Quintic model y^2 = f, W = infinity, P0 = (0, A(0)), P1 = (1, B(1)), with
       B - A = x g1,  B + A = g2,  g1 g2 = lambda x^4 + (x-1)^5,  g1 quadratic,
       f = A^2 + lambda x^5 = B^2 - x (x-1)^5,  B(0) = A(0),
   so div(y - A) = 5P0 - 5W and div(y - B) = P0 + 5P1 - 6W: lattice (5,0), (1,5), determinant 25, cyclic.
   Same rational curve of quadratics g1 as in family_24_certificate.m.
   Certificate: t = 3 mod 7 gives #J(F_7) = 50; t = 3 mod 11 gives #J(F_11) = 75 with a strict
   Frobenius; so every t = 3 mod 77 (7- and 11-integral) has torsion exactly [25] and a geometrically
   simple Jacobian (family_cert_lib.m).
   Run:  magma -b pipeline/magma/family_25_certificate.m > data/knowledge/sources/family_25_certificate.log */
SetColumns(0);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
e1 := (-33*t^4 + 65*t^3 - 42*t^2 + 11*t - 1)/t^4; e2 := (-9*t^3 + 15*t^2 - 7*t + 1)/t^3;
g1 := x^2 - e1*x + e2;
r4 := x^4 mod g1; r5 := (x-1)^5 mod g1;
lam := -Coefficient(r5, 1)/Coefficient(r4, 1);
assert lam*r4 + r5 eq 0;
Rx := lam*x^4 + (x-1)^5; g2 := Rx div g1; assert g1*g2 eq Rx;
BmA := x*g1; BpA := g2;
A := (BpA - BmA)/2; B := (BpA + BmA)/2;
f := A^2 + lam*x^5;
assert f eq B^2 - x*(x-1)^5 and Degree(A) le 2 and Degree(f) eq 5 and IsMonic(B) and Degree(B) eq 3;
assert Evaluate(B, 0) eq Evaluate(A, 0);
printf "lambda = %o\n", lam;
J := Jacobian(HyperellipticCurve(f));
D0 := elt< J | x, Evaluate(A, 0) >; D1 := elt< J | x - 1, Evaluate(B, 1) >;
assert 5*D0 eq J!0 and D0 + 5*D1 eq J!0;
assert 25*D1 eq J!0 and 5*D1 ne J!0;
print "over Q(t): 5 D0 = D0 + 5 D1 = 0 and D1 has order 25";
g, st := ResidueCertificate(f, [7, 11], [3, 3], [3, 80, -74, 157], 25, [25]);
print "FAMILY_25_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [25]";
quit;
