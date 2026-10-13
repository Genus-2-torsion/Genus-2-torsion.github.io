/* family_13_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [13].

   Relation-lattice family (the method of torsion_jac/notes/z38_2026_10_13.md section 2; this family was
   found in the notes review of 2026-10-13, torsion_jac/notes/review_of_notes_2026_10_13.md).  Quintic
   model y^2 = f_t, W = infinity, P0 = (0, A(0)), P1 = (1, A(1)), D_i = P_i - W:
       B - A = x(x-1)^2,  B + A = (x-1)^3 + t x^2,
       f_t = A^2 + t x^3 (x-1)^2 = B^2 - x (x-1)^5,  with B(0) = A(0), B(1) = A(1).
   Then div(y - A) = 3P0 + 2P1 - 5W and div(y - B) = P0 + 5P1 - 6W, so 3D0 + 2D1 = D0 + 5D1 = 0, with
   determinant 13: <D0, D1> is cyclic of order 13 (checked over Q(t)).
   Certificate: every 3-integral t = 1 mod 3 has the reduction of t = 1 at 3, with #J(F_3) = 13 and a
   strict Frobenius.  The argument is in family_cert_lib.m.
   Run:  magma -b pipeline/magma/family_13_certificate.m > data/knowledge/sources/family_13_certificate.log */
SetColumns(0);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
BmA := x*(x-1)^2; BpA := (x-1)^3 + t*x^2;
A := (BpA - BmA)/2; B := (BpA + BmA)/2;
f := A^2 + t*x^3*(x-1)^2;
assert f eq B^2 - x*(x-1)^5 and Degree(A) le 2 and Degree(f) eq 5 and IsMonic(B) and Degree(B) eq 3;
assert Evaluate(B, 0) eq Evaluate(A, 0) and Evaluate(B, 1) eq Evaluate(A, 1);
printf "f_t = %o\n", f;
J := Jacobian(HyperellipticCurve(f));
D0 := elt< J | x, Evaluate(A, 0) >; D1 := elt< J | x - 1, Evaluate(A, 1) >;
assert 3*D0 + 2*D1 eq J!0 and D0 + 5*D1 eq J!0 and D1 ne J!0 and 13*D1 eq J!0;
print "over Q(t): 3 D0 + 2 D1 = D0 + 5 D1 = 0 and D1 has order 13";
g, st := ResidueCertificate(f, [3], [1], [4, 7, -2, 1/4], 13, [13]);
print "FAMILY_13_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [13]";
quit;
