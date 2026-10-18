/* family_2_2_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,2].

   The record's family (census_infinity_record.md, section 2.1; torsion_jac/code/claude_census_infty_certify_b.m):
       C_t : y^2 = f_t = x(x-1)(x-2)(x^3 + x + t)   (degree 6).
   Argument.
   (1) Over Q(t): the root-pair classes T1 = [x(x-1), 0] and T2 = [x(x-2), 0] are distinct and nonzero,
       of order 2, so [2,2] is marked.
   (2) Residue certificate (family_cert_lib.m) at t = 3 mod 7 and t = 3 mod 17, i.e. every 7- and
       17-integral t = 3 mod 119: #J(F_7) = 72, #J(F_17) = 316, gcd 4 = #[2,2]; the Frobenius at 17 is
       strict.  The record's fiber t = 3 lies in the class and has TorsionSubgroup [2,2].
   (3) The G2-invariants of the members differ, so t -> M_2 is nonconstant.
   (Named family_2_2 because family_22 would be the name for the cyclic group [22].)
   Run (from the repository root):
       magma -b pipeline/magma/family_2_2_certificate.m > data/knowledge/sources/family_2_2_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
f := x*(x-1)*(x-2)*(x^3 + x + t);
assert Degree(f) eq 6;
J := Jacobian(HyperellipticCurve(f));
T1 := elt< J | x*(x-1), 0 >; T2 := elt< J | x*(x-2), 0 >;
assert T1 ne J!0 and T2 ne J!0 and T1 ne T2 and 2*T1 eq J!0 and 2*T2 eq J!0;
print "over Q(t): T1 = [x(x-1), 0] and T2 = [x(x-2), 0] are distinct classes of order 2: [2,2] is marked";
g, st := ResidueCertificate(f, [7, 17], [3, 3], [Q| 3, 122, -116, 241], 4, [2,2]);
print "FAMILY_2_2_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,2]";
quit;
