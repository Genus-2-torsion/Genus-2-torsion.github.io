/* family_3_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [3].

   The record's family (census_infinity_record.md, section 2.2): the contact form f = h^2 - lambda Q^3 with
   h = x^3 + x + 1, Q = x^2 + 1, lambda = t:
       C_t : y^2 = f_t = (x^3 + x + 1)^2 - t (x^2 + 1)^3.
   With D = [x^2 + 1, h mod (x^2 + 1)] = [x^2 + 1, 1] (the two points over the roots of Q with y = h, minus
   the divisor at infinity), div(y - h) = 3D, so D has order 3 (checked over Q(t)).
   Certificate: every 3-integral t = 2 mod 3 has the reduction of t = 2 at 3, with #J(F_3) = 3 and a
   strict Frobenius (family_cert_lib.m).
   Run:  magma -b pipeline/magma/family_3_certificate.m > data/knowledge/sources/family_3_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
h := x^3 + x + 1; Qd := x^2 + 1;
f := h^2 - t*Qd^3;
assert Degree(f) eq 6 and IsSquarefree(f) and h mod Qd eq 1;
printf "f_t = %o\n", f;
J := Jacobian(HyperellipticCurve(f));
D := elt< J | Qd, h mod Qd >;
assert D ne J!0 and 3*D eq J!0;
print "over Q(t): D = [x^2 + 1, 1] has order 3";
g, st := ResidueCertificate(f, [3], [2], [2, 5, -1, 1/2], 3, [3]);
print "FAMILY_3_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [3]";
quit;
