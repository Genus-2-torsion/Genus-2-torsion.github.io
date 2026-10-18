/* family_2_12_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,12].

   The record's line a = (1-r)/4 in the M(12) chart (census_infinity_record.md, section 3;
   torsion_jac/notes/m12_simple_route.md section 4a; torsion_jac/code/m12_z12x2_family.m):
   T = a x^2 - x + r, y^2 + (x-r)(T+1) y = a x^2 T(T+1), P = (0,0), with r = t.  With Y = 2y + (x-r)(T+1)
   the curve is Y^2 = W = ((x-r)(T+1))^2 + 4a x^2 T(T+1), which has the rational roots 2 and
   2(t+1)/(1-t) (the roots of T+1) and u = -2t/(t-2).
   D = [x(x - 2), b] with b(0) = -t(t+1), b(2) = 0 is P + (2,0) - K ~ P - (2,0), of order 12, and
   T2 = [(x - u)(x - 2), 0] is a 2-torsion class with T2 != 6D, so <D, T2> = Z/12 x Z/2 (checked over Q(t)).
   Certificate: t = 2 mod 5 gives #J(F_5) = 24; t = 2 mod 7 gives #J(F_7) = 48 with a strict Frobenius;
   so every 5- and 7-integral t = 2 mod 35 with t != 2 has torsion exactly [2,12] and a geometrically
   simple Jacobian (family_cert_lib.m).  (At t = 2 the root u is not defined; for t = 2 mod 5 the root u
   reduces to infinity and W reduces to a squarefree quintic, which is still good reduction.)
   Run:  magma -b pipeline/magma/family_2_12_certificate.m > data/knowledge/sources/family_2_12_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
r := t; a := (1 - r)/4;
T := a*x^2 - x + r; h := (x-r)*(T+1); F := a*x^2*T*(T+1);
W := h^2 + 4*F;
assert Degree(W) eq 6 and IsSquarefree(W);
w := Ft!2; u := -2*r/(r-2);
assert Evaluate(T+1, w) eq 0 and Evaluate(W, u) eq 0 and u ne w;
printf "W_t = %o\n", W;
J := Jacobian(HyperellipticCurve(W));
b := Evaluate(h, 0)*(w - x)/w;
D := elt< J | x*(x - w), b >;
T2 := elt< J | (x - u)*(x - w), 0 >;
assert 12*D eq J!0 and 6*D ne J!0 and 4*D ne J!0;
assert T2 ne J!0 and 2*T2 eq J!0 and T2 ne 6*D;
print "over Q(t): D = [x(x - 2), b] ~ P - (2,0) has order 12, and T2 = [(x - u)(x - 2), 0], u = -2t/(t-2), is a 2-torsion class outside <D>: [2,12] is marked";
g, st := ResidueCertificate(W, [5, 7], [2, 2], [37, 1/18, -33, 72], 24, [2,12]);
print "FAMILY_2_12_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,12]";
quit;
