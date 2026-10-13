/* family_214_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,14].

   The contact-7 chart h = 1 - (7/2)x + a x^2 + b x^3, f = (h^2 + (x-1)^7)/x^2, on which x^2 f - h^2 =
   (x-1)^7 gives div(xy - h) = 7P - 7W with P = (1, h(1)): the class [P - oo] has order 7.  The
   two-root surface of torsion_jac/notes/contact7_two_root_surface_2026_07_11.md: with
       V(z) = 2z^5 + 4z^4 + 6z^3 + 8z^2 + 10z + 5,  A(z) = V(z)/(2(z+1)^2),
       b = (A(t) - A(s))/(s^2 - t^2),  a = A(s) - b(1 - s^2),
   f has the two rational roots 1 - s^2 and 1 - t^2.  We use the slice s = 1, with parameter t.
   (Certificate found in the notes review of 2026-10-13.)
   Certificate: t = 2 mod 5 gives #J(F_5) = 28 = #[2,14]; t = 2 mod 13 gives a strict Frobenius; so every
   t = 2 mod 65 (5- and 13-integral) has torsion exactly [2,14] and a geometrically simple Jacobian
   (family_cert_lib.m).
   Run:  magma -b pipeline/magma/family_214_certificate.m > data/knowledge/sources/family_214_certificate.log */
SetColumns(0);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
V := func< z | 2*z^5 + 4*z^4 + 6*z^3 + 8*z^2 + 10*z + 5 >;
Az := func< z | V(z)/(2*(z+1)^2) >;
s := Ft!1;
b := (Az(t) - Az(s))/(s^2 - t^2);
a := Az(s) - b*(1 - s^2);
h := 1 - (7/2)*x + a*x^2 + b*x^3;
num := h^2 + (x-1)^7;
f := num div x^2; assert f*x^2 eq num;
assert Degree(f) eq 5;
assert Evaluate(f, 1 - s^2) eq 0 and Evaluate(f, 1 - t^2) eq 0;
J := Jacobian(HyperellipticCurve(f));
D7 := elt< J | x - 1, Evaluate(h, 1) >;
assert D7 ne J!0 and 7*D7 eq J!0;
T1 := elt< J | x - (1 - s^2), 0 >; T2 := elt< J | x - (1 - t^2), 0 >;
assert T1 ne J!0 and T2 ne J!0 and T1 ne T2 and 2*T1 eq J!0 and 2*T2 eq J!0;
print "over Q(t), s = 1: [P - oo] has order 7 and the roots 0, 1 - t^2 give (Z/2)^2: [2,14] is marked";
g, st := ResidueCertificate(f, [5, 13], [2, 2], [2, 67, -63, 132], 28, [2,14]);
print "FAMILY_214_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,14]";
quit;
