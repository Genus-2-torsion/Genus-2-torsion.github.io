/* family_9_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [9].

   The record's contact-9 chart (census_infinity_record.md, section 3; torsion_jac/notes/contact9_family.md):
       h = 1 - (9/2)x + (63/8)x^2 - (105/16)x^3 + t x^4,  f_t = (h^2 + (x-1)^9)/x^4.
   Since x^4 f - h^2 = (x-1)^9, div(x^2 y - h) = 9P - 9 oo with P = (1, h(1)), so D = [P - oo] = [x - 1, h(1)]
   has order 9 (checked over Q(t)).
   Certificate: every 3-integral t = 1 mod 3 has the reduction of t = 1 at 3, with #J(F_3) = 9 and a
   strict Frobenius (family_cert_lib.m).
   Run:  magma -b pipeline/magma/family_9_certificate.m > data/knowledge/sources/family_9_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
h := 1 - (9/2)*x + (63/8)*x^2 - (105/16)*x^3 + t*x^4;
num := h^2 + (x-1)^9;
f := num div x^4;
assert f*x^4 eq num and Degree(f) eq 5 and IsSquarefree(f);
printf "f_t = %o\n", f;
J := Jacobian(HyperellipticCurve(f));
D := elt< J | x - 1, Evaluate(h, 1) >;
assert 9*D eq J!0 and 3*D ne J!0;
print "over Q(t): x^4 f - h^2 = (x-1)^9 and D = [x - 1, h(1)] has order 9";
g, st := ResidueCertificate(f, [3], [1], [1, 4, -2, 1/4], 9, [9]);
print "FAMILY_9_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [9]";
quit;
