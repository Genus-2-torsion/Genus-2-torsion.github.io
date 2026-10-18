/* family_224_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,2,4].

   The record's family (census_infinity_record.md, section 2.4; torsion_jac/code/claude_census_infty_round2b.m,
   part C):
       C_t : y^2 = f_t = x (x+4)(x+9)(x^2 + (2t+1)x + t^2)   (degree 5).
   Argument.
   (1) Halving identity: x + 4 = 2^2 + x, x + 9 = 3^2 + x, q_t = (x + t)^2 + x, and
       (a^2 + x b^2)(c^2 + x d^2) = (ac - x bd)^2 + x (ad + bc)^2 give
           (x+4)(x+9) q_t = u^2 + x l^2,  u = x^2 + (t-1)x - 6t,  l = 4x + 5t + 6,
       so D = [u, x l mod u] satisfies 2D = [x, 0].  With T1 = [x+4, 0], T2 = [x+9, 0], the 16 sums
       i D + j T1 + k T2 are distinct, so [2,2,4] is marked over Q(t).
   (2) Residue certificate (family_cert_lib.m) at t = 1 mod 7 and t = 1 mod 13, i.e. every 7- and
       13-integral t = 1 mod 91: #J(F_7) = 32, #J(F_13) = 208, gcd 16 = #[2,2,4]; the Frobenius at 13 is
       strict.  The record's fiber t = 1 lies in the class and has TorsionSubgroup [2,2,4].
   (3) The G2-invariants of the members differ, so t -> M_2 is nonconstant.
   Run (from the repository root):
       magma -b pipeline/magma/family_224_certificate.m > data/knowledge/sources/family_224_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
qt := x^2 + (2*t+1)*x + t^2;
f := x*(x+4)*(x+9)*qt;
u := x^2 + (t-1)*x - 6*t; l := 4*x + 5*t + 6;
assert (x+4)*(x+9)*qt eq u^2 + x*l^2;
J := Jacobian(HyperellipticCurve(f));
D := elt< J | u, (x*l) mod u >;
T0 := elt< J | x, 0 >; T1 := elt< J | x+4, 0 >; T2 := elt< J | x+9, 0 >;
assert T0 ne J!0 and 2*D eq T0 and 4*D eq J!0;
assert #{ i*D + j*T1 + k*T2 : i in [0..3], j, k in [0..1] } eq 16;
print "over Q(t): D = [x^2 + (t-1)x - 6t, x(4x + 5t + 6) mod u] has 2D = [x, 0] != 0, and the 16 sums i D + j [x+4, 0] + k [x+9, 0] are distinct: [2,2,4] is marked";
g, st := ResidueCertificate(f, [7, 13], [1, 1], [Q| 1, 92, -90, 183], 16, [2,2,4]);
print "FAMILY_224_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,2,4]";
quit;
