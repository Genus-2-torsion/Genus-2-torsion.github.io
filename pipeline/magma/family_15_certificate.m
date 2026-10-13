/* family_15_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [15].

   Relation-lattice family (method of torsion_jac/notes/z38_2026_10_13.md section 2; found in the notes
   review of 2026-10-13).  Quintic model y^2 = f_t, W = infinity, P_i = (i, A(i)), D_i = P_i - W:
       B - A = t x^2 + (x-1)^3,  B + A = x^3,
       f_t = A^2 + t x^5 = B^2 - x^3 (x-1)^3,  with B(0) = -A(0).
   Then div(y - A) = 5P0 - 5W and div(y - B) = 3 iota(P0) + 3P1 - 6W, so 5D0 = -3D0 + 3D1 = 0: the lattice
   (5,0), (-3,3) has determinant 15 and entries with gcd 1, so <D0, D1> is cyclic of order 15.
   Certificate: every 3-integral t = 1 mod 3 has the reduction of t = 1 at 3, with #J(F_3) = 15 and a
   strict Frobenius (family_cert_lib.m).
   Run:  magma -b pipeline/magma/family_15_certificate.m > data/knowledge/sources/family_15_certificate.log */
SetColumns(0);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
BmA := t*x^2 + (x-1)^3; BpA := x^3;
A := (BpA - BmA)/2; B := (BpA + BmA)/2;
f := A^2 + t*x^5;
assert f eq B^2 - x^3*(x-1)^3 and Degree(A) le 2 and Degree(f) eq 5 and IsMonic(B) and Degree(B) eq 3;
assert Evaluate(B, 0) eq -Evaluate(A, 0);
printf "f_t = %o\n", f;
J := Jacobian(HyperellipticCurve(f));
D0 := elt< J | x, Evaluate(A, 0) >; D1 := elt< J | x - 1, Evaluate(B, 1) >;
assert 5*D0 eq J!0 and -3*D0 + 3*D1 eq J!0;
assert 15*D1 eq J!0 and 5*D1 ne J!0 and 3*D1 ne J!0;
print "over Q(t): 5 D0 = -3 D0 + 3 D1 = 0 and D1 has order 15";
g, st := ResidueCertificate(f, [3], [1], [4, 7, -2, 1/4], 15, [15]);
print "FAMILY_15_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [15]";
quit;
