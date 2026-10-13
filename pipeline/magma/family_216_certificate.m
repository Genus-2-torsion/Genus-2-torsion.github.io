/* family_216_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,16].

   On the [16] family of family_16_certificate.m (parameter s), the quadratic factor of f_s splits over
   Q exactly when s = 3/4 + u^2 (found in the notes review of 2026-10-13).  Over Q(u), f has two rational
   roots, so besides the cyclic class of order 16 it has a rational 2-torsion class T = [x - r, 0] outside
   it; together they give [2,16], marked over Q(u).
   Certificate: every 5-integral u = 1 mod 5 has the reduction of u = 1 at 5, with #J(F_5) = 32 = #[2,16]
   and a strict Frobenius (family_cert_lib.m).
   Run:  magma -b pipeline/magma/family_216_certificate.m > data/knowledge/sources/family_216_certificate.log */
SetColumns(0);
load "pipeline/magma/family_cert_lib.m";
Fu<u> := FunctionField(Q); P<x> := PolynomialRing(Fu);
s := 3/4 + u^2;
lam := -(s-1)^2/s;
R := lam*x + (x-1)^2;
q := R div (x - s); assert q*(x - s) eq R;
BmA := (x-1)^2*(x-s); BpA := x^2*q;
A := (BpA - BmA)/2; B := (BpA + BmA)/2;
f := A^2 + lam*x^3*(x-1)^2;
assert f eq B^2 - x^2*(x-1)^4 and Degree(A) le 2 and Degree(f) eq 5;
printf "f_u = %o\n", f;
J := Jacobian(HyperellipticCurve(f));
D0 := elt< J | x, Evaluate(A, 0) >; D1 := elt< J | x - 1, Evaluate(A, 1) >;
gen := (8*D1 ne J!0) select D1 else D0;
assert 16*gen eq J!0 and 8*gen ne J!0;
lin := [ g[1] : g in Factorization(f) | Degree(g[1]) eq 1 ];
assert #lin eq 2;
T := elt< J | lin[1], 0 >;
assert T ne J!0 and 2*T eq J!0 and T ne 8*gen;
printf "over Q(u): f has the rational roots %o; <gen> = Z/16 and T = [x - r, 0] is outside it: [2,16] is marked\n", [ -Coefficient(l, 0) : l in lin ];
g, st := ResidueCertificate(f, [5], [1], [1, 6, -4, 11, 1/6], 32, [2,16]);
print "FAMILY_216_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,16]";
quit;
