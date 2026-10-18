/* family_12_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [12].

   The record's rational-Weierstrass subchart of the M(12) chart (census_infinity_record.md, section 3;
   torsion_jac/notes/m12_simple_route.md): T = a x^2 - x + r, y^2 + (x-r)(T+1) y = a x^2 T(T+1),
   P = (0,0), with 1 - 4a(r+1) = z^2, i.e. a = (1 - z^2)/(4(r+1)).  Slice r = 3, parameter z = t.
   With Y = 2y + (x-r)(T+1) the curve is Y^2 = W = ((x-r)(T+1))^2 + 4a x^2 T(T+1), and W = (T+1) Q4.
   T+1 has the rational root w = 2(r+1)/(1-z), so W_ = (w, 0) is a rational Weierstrass point.
   P = (0, -r(r+1)) on Y^2 = W.  D = [x(x - w), b] with b(0) = -r(r+1), b(w) = 0 is P + W_ - K ~ P - W_,
   and has order 12 (checked over Q(t)).
   Certificate: t = 2 mod 7 gives #J(F_7) = 60; t = 4 mod 11 gives #J(F_11) = 168 with a strict
   Frobenius; so every 7- and 11-integral t = 37 mod 77 has torsion exactly [12] and a geometrically
   simple Jacobian (family_cert_lib.m).
   Run:  magma -b pipeline/magma/family_12_certificate.m > data/knowledge/sources/family_12_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
r := Ft!3; z := t;
a := (1 - z^2)/(4*(r+1));
T := a*x^2 - x + r; h := (x-r)*(T+1); F := a*x^2*T*(T+1);
W := h^2 + 4*F;
assert Degree(W) eq 6 and IsSquarefree(W) and W mod (T+1) eq 0;
w := 2*(r+1)/(1-z);
assert Evaluate(T+1, w) eq 0;
printf "W_t = %o\n", W;
J := Jacobian(HyperellipticCurve(W));
b := Evaluate(h, 0)*(w - x)/w;
D := elt< J | x*(x - w), b >;
assert 12*D eq J!0 and 6*D ne J!0 and 4*D ne J!0;
print "over Q(t): w = 8/(1-t) is a root of T+1, and D = [x(x - w), b] ~ P - (w, 0) has order 12";
g, st := ResidueCertificate(W, [7, 11], [2, 4], [-40, 37, 1/25, 114], 12, [12]);
print "FAMILY_12_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [12]";
quit;
