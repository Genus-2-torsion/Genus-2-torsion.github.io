/* family_2_2_2_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,2,2].

   The record's family (census_infinity_record.md, section 2.1; torsion_jac/code/claude_census_infty_certify_b.m):
       C_t : y^2 = f_t = x(x-1)(x-2)(x^2 + x + t)   (degree 5).
   Argument.
   (1) Over Q(t): the classes [x, 0], [x-1, 0], [x-2, 0] of the rational Weierstrass points are
       independent of order 2 (the 8 sums are distinct), so [2,2,2] is marked.
   (2) Residue certificate (family_cert_lib.m) at t = 3 mod 7 and t = 3 mod 83, i.e. every 7- and
       83-integral t = 3 mod 581: #J(F_7) = 32, #J(F_83) = 6824, gcd 8 = #[2,2,2]; the Frobenius at 83 is
       strict.  The record's fiber t = 3 lies in the class and has TorsionSubgroup [2,2,2].
   (3) The G2-invariants of the members differ, so t -> M_2 is nonconstant.
   (Named family_2_2_2 because family_222 would be the name for the group [2,22].)
   Run (from the repository root):
       magma -b pipeline/magma/family_2_2_2_certificate.m > data/knowledge/sources/family_2_2_2_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
f := x*(x-1)*(x-2)*(x^2 + x + t);
assert Degree(f) eq 5;
J := Jacobian(HyperellipticCurve(f));
Ts := [ elt< J | x - r, 0 > : r in [0, 1, 2] ];
assert &and[ 2*T eq J!0 : T in Ts ];
assert #{ i*Ts[1] + j*Ts[2] + k*Ts[3] : i, j, k in [0..1] } eq 8;
print "over Q(t): [x, 0], [x-1, 0], [x-2, 0] are independent classes of order 2: [2,2,2] is marked";
g, st := ResidueCertificate(f, [7, 83], [3, 3], [Q| 3, 584, -578, 1165], 8, [2,2,2]);
print "FAMILY_2_2_2_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,2,2]";
quit;
