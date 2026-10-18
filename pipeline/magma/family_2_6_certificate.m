/* family_2_6_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,6].

   The record's contact-6 chart (census_infinity_record.md, section 3; torsion_jac/notes/contact6_m36.md):
       h = 1 + a x + b x^2 + x^3,  f = h^2 - (x-1)^6,
   on the slice b = 1, a = t.  f has degree 5, and div(y - h) = 6P - 6 oo with P = (1, h(1)), so
   D6 = [x - 1, h(1)] has order 6.  Since h(0) = 1, f = (h - (x-1)^3) * x * ((h + (x-1)^3)/x) has the root
   x = 0, and T = [x, 0] has order 2.  3 D6 is the class [h - (x-1)^3, 0] of the quadratic factor, so
   T is not in <D6> and <D6, T> = Z/6 x Z/2 (checked over Q(t)).
   Certificate: t = 1 mod 3 gives #J(F_3) = 12; t = 0 mod 13 gives #J(F_13) = 120 with a strict
   Frobenius; so every 3- and 13-integral t = 13 mod 39 has torsion exactly [2,6] and a geometrically
   simple Jacobian (family_cert_lib.m).
   Run:  magma -b pipeline/magma/family_2_6_certificate.m > data/knowledge/sources/family_2_6_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
h := 1 + t*x + x^2 + x^3;
f := h^2 - (x-1)^6;
assert Degree(f) eq 5 and IsSquarefree(f) and Evaluate(f, 0) eq 0;
qA := h - (x-1)^3; assert Degree(qA) eq 2;
printf "f_t = %o\n", f;
J := Jacobian(HyperellipticCurve(f));
D6 := elt< J | x - 1, Evaluate(h, 1) >; T := elt< J | x, 0 >;
assert 6*D6 eq J!0 and 2*D6 ne J!0 and 3*D6 ne J!0;
assert 3*D6 eq elt< J | qA/LeadingCoefficient(qA), 0 >;
assert T ne J!0 and 2*T eq J!0 and T ne 3*D6;
print "over Q(t): D6 = [x - 1, h(1)] has order 6, 3 D6 = [h - (x-1)^3, 0], and T = [x, 0] is a 2-torsion class outside <D6>: [2,6] is marked";
g, st := ResidueCertificate(f, [3, 13], [1, 0], [13, -26, 52, -65], 12, [2,6]);
print "FAMILY_2_6_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,6]";
quit;
