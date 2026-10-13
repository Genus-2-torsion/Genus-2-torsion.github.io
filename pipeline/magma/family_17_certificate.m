/* family_17_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [17].

   Relation-lattice family (method of torsion_jac/notes/z38_2026_10_13.md section 2; found in the notes
   review of 2026-10-13).  Quintic model y^2 = f_t, W = infinity, P_i = (i, A(i)), D_i = P_i - W:
       lambda = -(t-1)^3/t^2,  (x-1)^3 + lambda x^2 = (x - t) q(x),
       B - A = (x-1)^2 (x-t),  B + A = x q(x),
       f_t = A^2 + lambda x^3 (x-1)^2 = B^2 - x (x-1)^5,  with B(0) = -A(0), B(1) = A(1).
   Then div(y - A) = 3P0 + 2P1 - 5W and div(y - B) = iota(P0) + 5P1 - 6W, so 3D0 + 2D1 = -D0 + 5D1 = 0,
   with determinant 17.
   Certificate: every 3-integral t = 2 mod 3 has the reduction of t = 2 at 3, with #J(F_3) = 17 and a
   strict Frobenius (family_cert_lib.m).
   Run:  magma -b pipeline/magma/family_17_certificate.m > data/knowledge/sources/family_17_certificate.log */
SetColumns(0);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
lam := -(t-1)^3/t^2;
cub := (x-1)^3 + lam*x^2;
q := cub div (x - t); assert q*(x - t) eq cub;
BmA := (x-1)^2*(x-t); BpA := x*q;
A := (BpA - BmA)/2; B := (BpA + BmA)/2;
f := A^2 + lam*x^3*(x-1)^2;
assert f eq B^2 - x*(x-1)^5 and Degree(A) le 2 and Degree(f) eq 5 and IsMonic(B) and Degree(B) eq 3;
assert Evaluate(B, 0) eq -Evaluate(A, 0) and Evaluate(B, 1) eq Evaluate(A, 1);
printf "f_t = %o\n", f;
J := Jacobian(HyperellipticCurve(f));
D0 := elt< J | x, Evaluate(A, 0) >; D1 := elt< J | x - 1, Evaluate(A, 1) >;
assert 3*D0 + 2*D1 eq J!0 and -D0 + 5*D1 eq J!0 and D1 ne J!0 and 17*D1 eq J!0;
print "over Q(t): 3 D0 + 2 D1 = -D0 + 5 D1 = 0 and D1 has order 17";
g, st := ResidueCertificate(f, [3], [2], [2, 5, -1, 1/2, 8], 17, [17]);
print "FAMILY_17_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [17]";
quit;
