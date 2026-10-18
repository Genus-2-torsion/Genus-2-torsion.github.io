/* family_20_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [20].

   The family: the record's "contact-5 + 4-torsion family" (torsion_jac/notes/contact5_order40_family.md).
   It is the halving component b = (a^2-1)/2 of the contact-5 chart h = 1 + a x + b x^2,
   f = h^2 - (1+a+b)^2 x^5.  With t = a:
       h = 1 + t x + ((t^2-1)/2) x^2,   f_t = h^2 - ((t+1)^4/4) x^5.
   Since h(1) = (t+1)^2/2, x = 1 is a root of f_t.

   Argument.
   (1) Over Q(t): div(y - h) = 5 (0,1) - 5 oo, so D5 = [x, 1] has order 5 (checked: D5 != 0, 5 D5 = 0);
       H = [x^2 + 2x/(t+1), (t+2)x + 1] satisfies 2H = [x - 1, 0] != 0 and 4H = 0, so H has order 4.
       So <D5, H> = Z/20 is marked over Q(t).
   (2) Residue certificate (family_cert_lib.m) at t = 1 mod 3 and t = 1 mod 7, i.e. t = 1 mod 21:
       #J(F_3) = 20 and #J(F_7) = 40, gcd 20 = #G, and the Frobenius at 7 is strict.  So every 3- and
       7-integral t = 1 mod 21 gives J(Q)_tors exactly [20] and a geometrically simple Jacobian.
   (3) The G2-invariants of two members differ, so infinitely many isomorphism classes.
   Run (from the census root):
       magma -b pipeline/magma/family_20_certificate.m > data/knowledge/sources/family_20_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// (1) the family and the marked group over Q(t)
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
h := 1 + t*x + ((t^2-1)/2)*x^2;
f := h^2 - ((t+1)^4/4)*x^5;
assert Degree(f) eq 5 and Evaluate(f, 1) eq 0;
J := Jacobian(HyperellipticCurve(f));
D5 := elt< J | x, 1 >;
assert D5 ne J!0 and 5*D5 eq J!0;
H := elt< J | x^2 + 2*x/(t+1), (t+2)*x + 1 >;
assert 2*H eq elt< J | x - 1, 0 > and 2*H ne J!0 and 4*H eq J!0;
D20 := D5 + H;
assert 10*D20 ne J!0 and 4*D20 ne J!0 and 20*D20 eq J!0;
print "over Q(t): D5 = [x, 1] has order 5, H = [x^2 + 2x/(t+1), (t+2)x + 1] has order 4 (2H = [x-1, 0]), D5 + H has order 20";

// (2), (3) residue certificate at t = 1 mod 3, t = 1 mod 7
g, st := ResidueCertificate(f, [3, 7], [1, 1], [1, 22, -20, 22/43], 20, [20]);
print "FAMILY_20_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [20]";
quit;
