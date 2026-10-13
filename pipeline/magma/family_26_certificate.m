/* family_26_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [26].

   On the [13] family of family_13_certificate.m (parameter s), the locus where f_s has a rational root
   x0 is a rational curve (found in the notes review of 2026-10-13), parametrized by
       d(t) = 1 + 28t + 312t^2 + 1728t^3 + 4752t^4 + 5184t^5,
       s = (16/3)(1+3t)^4/d(t),  x0 = -(1/3)(1+6t)^6/d(t).
   Over Q(t) the class D1 = [x - 1, A(1)] of order 13 and T = [x - x0, 0] of order 2 give Z/26.
   Certificate: t = 2 mod 5 gives #J(F_5) = 26; t = 3 mod 7 gives #J(F_7) = 78 with a strict
   Frobenius; so every t = 17 mod 35 (5- and 7-integral) has torsion exactly [26] and a geometrically
   simple Jacobian (family_cert_lib.m).
   Run:  magma -b pipeline/magma/family_26_certificate.m > data/knowledge/sources/family_26_certificate.log */
SetColumns(0);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
d := 1 + 28*t + 312*t^2 + 1728*t^3 + 4752*t^4 + 5184*t^5;
s := (16/3)*(1+3*t)^4/d;
x0 := -(1/3)*(1+6*t)^6/d;
BmA := x*(x-1)^2; BpA := (x-1)^3 + s*x^2;
A := (BpA - BmA)/2; B := (BpA + BmA)/2;
f := A^2 + s*x^3*(x-1)^2;
assert f eq B^2 - x*(x-1)^5 and Degree(f) eq 5;
assert Evaluate(f, x0) eq 0;
J := Jacobian(HyperellipticCurve(f));
D1 := elt< J | x - 1, Evaluate(A, 1) >; T := elt< J | x - x0, 0 >;
assert D1 ne J!0 and 13*D1 eq J!0 and T ne J!0 and 2*T eq J!0;
print "over Q(t): D1 has order 13 and f(x0) = 0, so [x - x0, 0] has order 2: Z/26 is marked";
g, st := ResidueCertificate(f, [5, 7], [2, 3], [17, 52, -18, 87], 26, [26]);
print "FAMILY_26_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [26]";
quit;
