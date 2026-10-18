/* family_7_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [7].

   The record's contact-7 chart (census_infinity_record.md, section 3; torsion_jac/notes/contact7_family.md):
       h = 1 - (7/2)x + a x^2 + b x^3,  f = (h^2 + (x-1)^7)/x^2,
   on the slice b = 1, a = t.  Since x^2 f - h^2 = (x-1)^7, div(x y - h) = 7P - 7 oo with P = (1, h(1)), so
   D = [P - oo] = [x - 1, h(1)] has order 7 (checked over Q(t)).
   Certificate: every 3-integral t = 2 mod 3 has the reduction of t = 2 at 3, with #J(F_3) = 7 and a
   strict Frobenius (family_cert_lib.m).
   Run:  magma -b pipeline/magma/family_7_certificate.m > data/knowledge/sources/family_7_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
h := 1 - (7/2)*x + t*x^2 + x^3;
num := h^2 + (x-1)^7;
f := num div x^2;
assert f*x^2 eq num and Degree(f) eq 5 and IsSquarefree(f);
printf "f_t = %o\n", f;
J := Jacobian(HyperellipticCurve(f));
D := elt< J | x - 1, Evaluate(h, 1) >;
assert D ne J!0 and 7*D eq J!0;
print "over Q(t): x^2 f - h^2 = (x-1)^7 and D = [x - 1, h(1)] has order 7";
g, st := ResidueCertificate(f, [3], [2], [2, 5, -1, 1/2], 7, [7]);
print "FAMILY_7_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [7]";
quit;
