/* family_220_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,20].

   The family: the record's "extra-2 loci" of the contact-5 + 4-torsion family
   (torsion_jac/notes/contact5_order40_family.md, the 1+3 locus).  The [20] family is
       h = 1 + t x + ((t^2-1)/2) x^2,   f_t = h^2 - ((t+1)^4/4) x^5   (x = 1 is a root).
   With u = t + 1 and x = y/u the residual quartic f_t/(x-1) is equivalent to
   u (y^4 + 4y^3 + 8y^2 + 8y + 4) - 4y(y+1)^2, which has the root y = z on the curve
       t = -(z^4 + 4z + 4)/(z^4 + 4z^3 + 8z^2 + 8z + 4).
   So on this rational curve f has the second rational root x2 = z/(t+1).

   Argument.
   (1) Over Q(z): D5 = [x, 1] has order 5, H = [x^2 + 2x/(t+1), (t+2)x + 1] has order 4 with
       2H = [x - 1, 0], and T2 = [x - x2, 0] has order 2 with T2 != 2H.  So <H, T2> = Z/4 x Z/2 and
       <D5, H, T2> = Z/2 x Z/20 is marked over Q(z).
   (2) Residue certificate (family_cert_lib.m) at z = 1 mod 7: #J(F_7) = 40 = #G, strict.  So every
       7-integral z = 1 mod 7 gives J(Q)_tors exactly [2,20] and a geometrically simple Jacobian.
   (3) The G2-invariants of two members differ, so infinitely many isomorphism classes.
   Run (from the census root):
       magma -b pipeline/magma/family_220_certificate.m > data/knowledge/sources/family_220_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// (1) the family and the marked group over Q(z)
Fz<z> := FunctionField(Q); P<x> := PolynomialRing(Fz);
t := -(z^4 + 4*z + 4)/(z^4 + 4*z^3 + 8*z^2 + 8*z + 4);
h := 1 + t*x + ((t^2-1)/2)*x^2;
f := h^2 - ((t+1)^4/4)*x^5;
x2 := z/(t+1);
assert Degree(f) eq 5 and Evaluate(f, 1) eq 0 and Evaluate(f, x2) eq 0 and x2 ne 1;
J := Jacobian(HyperellipticCurve(f));
D5 := elt< J | x, 1 >;
assert D5 ne J!0 and 5*D5 eq J!0;
H := elt< J | x^2 + 2*x/(t+1), (t+2)*x + 1 >;
assert 2*H eq elt< J | x - 1, 0 > and 2*H ne J!0 and 4*H eq J!0;
T2 := elt< J | x - x2, 0 >;
assert T2 ne J!0 and 2*T2 eq J!0 and T2 ne 2*H;
print "over Q(z): D5 = [x, 1] has order 5, H has order 4 with 2H = [x-1, 0], T2 = [x - z/(t+1), 0] has order 2 and T2 != 2H: [2,20] is marked";

// (2), (3) residue certificate at z = 1 mod 7
g, st := ResidueCertificate(f, [7], [1], [1, 8, -6, 8/15], 40, [2,20]);
print "FAMILY_220_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,20]";
quit;
