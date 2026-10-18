/* family_2222_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,2,2,2].

   The record's family (census_infinity_record.md, section 2.1; torsion_jac/code/claude_census_infty_certify.m):
       C_t : y^2 = f_t = x(x-1)(x-2)(x-3)(x-t)   (degree 5, five rational Weierstrass points).
   Argument.
   (1) Over Q(t): the classes [x - r, 0], r = 0, 1, 2, 3, are independent of order 2 (the 16 sums are
       distinct), so [2,2,2,2] = J[2] is marked.
   (2) Residue certificate (family_cert_lib.m) at t = 5 mod 7 and t = 5 mod 11, i.e. every 7- and
       11-integral t = 5 mod 77: #J(F_7) = 48, #J(F_11) = 80, gcd 16 = #[2,2,2,2]; the Frobenius at 11
       is strict.  The record's fiber t = 5 lies in the class and has TorsionSubgroup [2,2,2,2].
   (3) The G2-invariants of the members differ, so t -> M_2 is nonconstant.
   Run (from the repository root):
       magma -b pipeline/magma/family_2222_certificate.m > data/knowledge/sources/family_2222_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
f := x*(x-1)*(x-2)*(x-3)*(x-t);
assert Degree(f) eq 5;
J := Jacobian(HyperellipticCurve(f));
Ts := [ elt< J | x - r, 0 > : r in [0, 1, 2, 3] ];
assert &and[ 2*T eq J!0 : T in Ts ];
assert #{ i*Ts[1] + j*Ts[2] + k*Ts[3] + m*Ts[4] : i, j, k, m in [0..1] } eq 16;
print "over Q(t): [x - r, 0], r = 0, 1, 2, 3, are independent classes of order 2: [2,2,2,2] is marked";
g, st := ResidueCertificate(f, [7, 11], [5, 5], [Q| 5, 82, -72, 159], 16, [2,2,2,2]);
print "FAMILY_2222_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,2,2,2]";
quit;
