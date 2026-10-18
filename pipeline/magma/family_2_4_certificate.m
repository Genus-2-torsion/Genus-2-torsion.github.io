/* family_2_4_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,4].

   The record's family (census_infinity_record.md, section 2.4; torsion_jac/code/claude_census_infty_round2c.m):
       C_t : y^2 = f_t = x (x^2 + (2t+1)x + t^2)(x^2 + 3x + 1)   (degree 5),
   both quadratics of the form q_al = (x + al)^2 + x, with al = t and al = 1.
   Argument.
   (1) Halving identity: q_t q_1 = u^2 + x l^2 with u = (x + t)(x + 1) - x = x^2 + t x + t, l = 2x + t + 1,
       so f - (x l)^2 = x u^2 and D = [u, x l mod u] satisfies 2D = [x, 0].  With T = [x^2 + 3x + 1, 0]
       (order 2, not equal to 2D), the 8 sums i D + j T (i in Z/4, j in Z/2) are distinct, so [2,4] is
       marked over Q(t).
   (2) Residue certificate (family_cert_lib.m) at t = 3 mod 7 and t = 3 mod 17, i.e. every 7- and
       17-integral t = 3 mod 119: #J(F_7) = 56, #J(F_17) = 368, gcd 8 = #[2,4]; the Frobenius at 17 is
       strict.  The record's fiber t = 3 lies in the class and has TorsionSubgroup [2,4].
   (3) The G2-invariants of the members differ, so t -> M_2 is nonconstant.
   (Named family_2_4 because family_24_certificate.m is the census certificate for [24].)
   Run (from the repository root):
       magma -b pipeline/magma/family_2_4_certificate.m > data/knowledge/sources/family_2_4_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
qt := x^2 + (2*t+1)*x + t^2; q1 := x^2 + 3*x + 1;
f := x*qt*q1;
u := x^2 + t*x + t; l := 2*x + t + 1;
assert qt*q1 eq u^2 + x*l^2;
J := Jacobian(HyperellipticCurve(f));
D := elt< J | u, (x*l) mod u >;
T0 := elt< J | x, 0 >; T := elt< J | q1, 0 >;
assert T0 ne J!0 and 2*D eq T0 and 4*D eq J!0;
assert T ne J!0 and 2*T eq J!0 and T ne T0;
assert #{ i*D + j*T : i in [0..3], j in [0..1] } eq 8;
print "over Q(t): D = [x^2 + t x + t, x(2x + t + 1) mod u] has 2D = [x, 0] != 0, T = [x^2 + 3x + 1, 0] != 2D, and the 8 sums i D + j T are distinct: [2,4] is marked";
g, st := ResidueCertificate(f, [7, 17], [3, 3], [Q| 3, 122, -116, 241], 8, [2,4]);
print "FAMILY_2_4_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,4]";
quit;
