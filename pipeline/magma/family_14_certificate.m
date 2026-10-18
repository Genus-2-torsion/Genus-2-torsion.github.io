/* family_14_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [14].

   The record's rational-root subfamily of the contact-7 chart (census_infinity_record.md, section 3;
   torsion_jac/notes/contact7_family.md, "Rational-root subfamily"): h = 1 - (7/2)x + a x^2 + b x^3,
   f = (h^2 + (x-1)^7)/x^2, with r = 1 - s^2 and h(r) = eps s^7, i.e.
       a = (eps s^7 - 1 + (7/2) r - b r^3)/r^2,
   so that f(r) = (s^14 + (-s^2)^7)/r^2 = 0.  Slice eps = 1, b = 1, parameter s = t.
   D7 = [x - 1, h(1)] has order 7 (div(x y - h) = 7P - 7 oo) and T = [x - r, 0] has order 2, so D7 + T has
   order 14 (checked over Q(t)).
   Certificate: t = 1 mod 5 gives #J(F_5) = 28; t = 3 mod 7 gives #J(F_7) = 70 with a strict Frobenius;
   so every 5- and 7-integral t = 31 mod 35 has torsion exactly [14] and a geometrically simple Jacobian
   (family_cert_lib.m).
   Run:  magma -b pipeline/magma/family_14_certificate.m > data/knowledge/sources/family_14_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
s := t; eps := 1; b := Ft!1; r := 1 - s^2;
a := (eps*s^7 - 1 + (7/2)*r - b*r^3)/r^2;
h := 1 - (7/2)*x + a*x^2 + b*x^3;
num := h^2 + (x-1)^7;
f := num div x^2;
assert f*x^2 eq num and Degree(f) eq 5 and IsSquarefree(f) and Evaluate(f, r) eq 0;
printf "f_t = %o\n", f;
J := Jacobian(HyperellipticCurve(f));
D7 := elt< J | x - 1, Evaluate(h, 1) >; T := elt< J | x - r, 0 >;
assert D7 ne J!0 and 7*D7 eq J!0 and T ne J!0 and 2*T eq J!0;
D14 := D7 + T;
assert 14*D14 eq J!0 and 7*D14 ne J!0 and 2*D14 ne J!0;
print "over Q(t): D7 = [x - 1, h(1)] has order 7, f(1 - t^2) = 0 and T = [x - (1 - t^2), 0] has order 2: D7 + T has order 14";
g, st := ResidueCertificate(f, [5, 7], [1, 3], [-4, 31, -39, 66], 14, [14]);
print "FAMILY_14_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [14]";
quit;
