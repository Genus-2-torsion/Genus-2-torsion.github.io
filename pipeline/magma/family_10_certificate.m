/* family_10_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [10].

   The record's family (census_infinity_record.md, section 2.2): the contact-5 chart h = 1 + a x + b x^2,
   f = h^2 - (1+a+b)^2 x^5 at b = 1, a = t:
       C_t : y^2 = f_t = (1 + t x + x^2)^2 - (t+2)^2 x^5.
   With P = (0, 1), div(y - h) = 5P - 5 oo, so D5 = [x, 1] has order 5; f(1) = h(1)^2 - (t+2)^2 = 0
   identically, so T = [x - 1, 0] has order 2.  D5 + T has order 10 (checked over Q(t)).
   Certificate: t = 0 mod 5 gives #J(F_5) = 20, t = 1 mod 7 gives #J(F_7) = 30, both strict; so every
   5- and 7-integral t = 15 mod 35 has torsion exactly [10] and a geometrically simple Jacobian
   (family_cert_lib.m).
   Run:  magma -b pipeline/magma/family_10_certificate.m > data/knowledge/sources/family_10_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
h := 1 + t*x + x^2;
f := h^2 - (t+2)^2*x^5;
assert Degree(f) eq 5 and IsSquarefree(f) and Evaluate(f, 1) eq 0;
printf "f_t = %o\n", f;
J := Jacobian(HyperellipticCurve(f));
D5 := elt< J | x, 1 >; T := elt< J | x - 1, 0 >;
assert D5 ne J!0 and 5*D5 eq J!0 and T ne J!0 and 2*T eq J!0;
D10 := D5 + T;
assert 10*D10 eq J!0 and 5*D10 ne J!0 and 2*D10 ne J!0;
print "over Q(t): D5 = [x, 1] has order 5, T = [x - 1, 0] has order 2, so D5 + T has order 10";
g, st := ResidueCertificate(f, [5, 7], [0, 1], [15, -20, 50, 85], 10, [10]);
print "FAMILY_10_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [10]";
quit;
