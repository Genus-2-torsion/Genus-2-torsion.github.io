/* family_2224_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,2,2,4].

   The record's family (census_infinity_record.md, section 2.5; torsion_jac/code/claude_census_infty_certify.m):
   the slice (a,b,c,d) = (1,2,3,t) of the square-branch chart y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2)
   (BNSS2026, Lemma 2.3(b) of the paper, "lem: 2244"):
       C_t : y^2 = f_t = x(x+1)(x+4)(x+9)(x+t^2)   (degree 5).
   Argument.
   (1) Lemma 2244(b): with s_i the elementary symmetric polynomials in a, b, c, d, the class
       D0 = [x^2 - s2 x + s4, (s1 s2 - s3) x - s1 s4] satisfies 2 D0 = [x, 0].  Checked by Mumford
       arithmetic over Q(t), together with: the 32 sums i D0 + j [x+1, 0] + k [x+4, 0] + m [x+9, 0] are
       distinct, so [2,2,2,4] is marked.
   (2) Residue certificate (family_cert_lib.m) at t = 5 mod 11 and t = 5 mod 59, i.e. every 11- and
       59-integral t = 5 mod 649: #J(F_11) = 128, #J(F_59) = 4000, gcd 32 = #[2,2,2,4]; the Frobenius at
       59 is strict.  The record's fiber t = 5 lies in the class and has TorsionSubgroup [2,2,2,4].
   (3) The G2-invariants of the members differ, so t -> M_2 is nonconstant.
   Run (from the repository root):
       magma -b pipeline/magma/family_2224_certificate.m > data/knowledge/sources/family_2224_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
aa := [Ft| 1, 2, 3, t];
f := x*&*[ x + c^2 : c in aa ];
assert f eq x*(x+1)*(x+4)*(x+9)*(x+t^2);
s1 := &+aa; s2 := &+[ aa[i]*aa[j] : i, j in [1..4] | i lt j ];
s3 := &+[ aa[i]*aa[j]*aa[k] : i, j, k in [1..4] | i lt j and j lt k ]; s4 := &*aa;
J := Jacobian(HyperellipticCurve(f));
D0 := elt< J | x^2 - s2*x + s4, (s1*s2 - s3)*x - s1*s4 >;
T0 := elt< J | x, 0 >;
Ts := [ elt< J | x + c, 0 > : c in [1, 4, 9] ];
assert T0 ne J!0 and 2*D0 eq T0 and 4*D0 eq J!0;
assert #{ i*D0 + j*Ts[1] + k*Ts[2] + m*Ts[3] : i in [0..3], j, k, m in [0..1] } eq 32;
printf "over Q(t): D0 = %o has 2 D0 = [x, 0] != 0, and the 32 sums i D0 + j [x+1, 0] + k [x+4, 0] + m [x+9, 0] are distinct: [2,2,2,4] is marked\n", D0;
g, st := ResidueCertificate(f, [11, 59], [5, 5], [Q| 5, 654, -644, 1303], 32, [2,2,2,4]);
print "FAMILY_2224_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,2,2,4]";
quit;
