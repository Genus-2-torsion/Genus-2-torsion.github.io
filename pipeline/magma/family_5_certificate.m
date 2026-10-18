/* family_5_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [5].

   The record gives no formula for [5] ("five-torsion families").  The inventory's [5] row
   (torsion_jac/notes/infinite_families_inventory.md, "Quintic contact") is the normalized chart
   f = h^2 - (1+a+b)^2 x^5, h = 1 + a x + b x^2, on which x = 1 is a root of f identically: it carries a
   rational 2-torsion class and has generic torsion containing [10], so it cannot give exact [5].
   This certificate uses the same contact-5 construction without that normalization, slice a = b = 1:
       C_t : y^2 = f_t = (1 + x + x^2)^2 - t x^5.
   With P = (0, 1), div(y - h) = 5P - 5 oo, so D = [P - oo] = [x, 1] has order 5 (checked over Q(t)).
   Certificate: every 3-integral t = 1 mod 3 has the reduction of t = 1 at 3, with #J(F_3) = 5 and a
   strict Frobenius (family_cert_lib.m).
   Run:  magma -b pipeline/magma/family_5_certificate.m > data/knowledge/sources/family_5_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
h := 1 + x + x^2;
f := h^2 - t*x^5;
assert Degree(f) eq 5 and IsSquarefree(f);
printf "f_t = %o\n", f;
J := Jacobian(HyperellipticCurve(f));
D := elt< J | x, 1 >;
assert D ne J!0 and 5*D eq J!0;
print "over Q(t): D = [x, 1] has order 5";
g, st := ResidueCertificate(f, [3], [1], [4, 7, -2, 1/4], 5, [5]);
print "FAMILY_5_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [5]";
quit;
