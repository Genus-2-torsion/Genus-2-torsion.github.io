/* family_18_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [18].

   The record's rational-root subfamily of the contact-9 chart (census_infinity_record.md, section 3;
   torsion_jac/notes/contact9_family.md, "Rational-root subfamily"):
   h = 1 - (9/2)x + (63/8)x^2 - (105/16)x^3 + a x^4, f = (h^2 + (x-1)^9)/x^4, with r = 1 - s^2 and
   h(r) = eps s^9, i.e.
       a = (eps s^9 - (1 - (9/2) r + (63/8) r^2 - (105/16) r^3))/r^4,
   so that f(r) = (s^18 + (-s^2)^9)/r^4 = 0.  Here eps = 1 and s = t.
   D9 = [x - 1, h(1)] has order 9 (div(x^2 y - h) = 9P - 9 oo) and T = [x - r, 0] has order 2, so D9 + T
   has order 18 (checked over Q(t)).
   Certificate: t = 2 mod 7 gives #J(F_7) = 72; t = 10 mod 13 gives #J(F_13) = 198 with a strict
   Frobenius; so every 7- and 13-integral t = 23 mod 91 has torsion exactly [18] and a geometrically
   simple Jacobian (family_cert_lib.m).
   Run:  magma -b pipeline/magma/family_18_certificate.m > data/knowledge/sources/family_18_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
s := t; eps := 1; r := 1 - s^2;
a := (eps*s^9 - (1 - (9/2)*r + (63/8)*r^2 - (105/16)*r^3))/r^4;
h := 1 - (9/2)*x + (63/8)*x^2 - (105/16)*x^3 + a*x^4;
num := h^2 + (x-1)^9;
f := num div x^4;
assert f*x^4 eq num and Degree(f) eq 5 and IsSquarefree(f) and Evaluate(f, r) eq 0;
printf "f_t = %o\n", f;
J := Jacobian(HyperellipticCurve(f));
D9 := elt< J | x - 1, Evaluate(h, 1) >; T := elt< J | x - r, 0 >;
assert 9*D9 eq J!0 and 3*D9 ne J!0 and T ne J!0 and 2*T eq J!0;
D18 := D9 + T;
assert 18*D18 eq J!0 and 9*D18 ne J!0 and 6*D18 ne J!0;
print "over Q(t): D9 = [x - 1, h(1)] has order 9, f(1 - t^2) = 0 and T = [x - (1 - t^2), 0] has order 2: D9 + T has order 18";
g, st := ResidueCertificate(f, [7, 13], [2, 10], [1/4, 23, -68, 114], 18, [18]);
print "FAMILY_18_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [18]";
quit;
