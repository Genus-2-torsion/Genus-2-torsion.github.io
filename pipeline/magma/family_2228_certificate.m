/* family_2228_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,2,2,8].

   The record's family (census_infinity_record.md, section 2.5; torsion_jac/code/claude_census_infty_certify.m,
   family 8): the rational curve on the K3 surface M(2,2,2,8) of BNSS2026, on the square-branch chart
       (a, b, c, d) = (-4t^2(t+1)/(t^2+t+1)^2, -t/(t+1), 1, t),
       C_t : y^2 = f_t = x(x+a^2)(x+b^2)(x+c^2)(x+d^2)   (degree 5).
   (The paper's signs (4t^2(t+1)/(t^2+t+1)^2, t/(t+1), -1, -t) give the same curve.)
   Argument.
   (1) Lemma 2244(b) of BNSS2026: D0 = [x^2 - s2 x + s4, (s1 s2 - s3) x - s1 s4] satisfies 2 D0 = [x, 0].
       On this curve s2^2 = 4 s4 (the K3 equation (ab+ac+ad+bc+bd+cd)^2 = 4abcd), so
       x^2 - s2 x + s4 = (x - x0)^2 with x0 = s2/2, and D0 = 2[P - oo] for the rational point
       P = (x0, (s1 s2 - s3) x0 - s1 s4).  Mumford checks over Q(t): 2[P - oo] = D0, 4[P - oo] = [x, 0],
       8[P - oo] = 0, and the 64 sums i[P - oo] + j[x+a^2, 0] + k[x+b^2, 0] + m[x+1, 0] (i in Z/8;
       j, k, m in Z/2) are distinct, so [2,2,2,8] is marked.
   (2) Residue certificate (family_cert_lib.m) at t = 2 mod 13 and t = 2 mod 31, i.e. every 13- and
       31-integral t = 2 mod 403: #J(F_13) = 256, #J(F_31) = 1088, gcd 64 = #[2,2,2,8]; the Frobenius at
       31 is strict.  The record's fiber t = 2 lies in the class and has TorsionSubgroup [2,2,2,8].
   (3) The G2-invariants of the members differ, so t -> M_2 is nonconstant.
   Run (from the repository root):
       magma -b pipeline/magma/family_2228_certificate.m > data/knowledge/sources/family_2228_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
aa := [ -4*t^2*(t+1)/(t^2+t+1)^2, -t/(t+1), Ft!1, t ];
f := x*&*[ x + c^2 : c in aa ];
assert Degree(f) eq 5;
J := Jacobian(HyperellipticCurve(f));
s1 := &+aa; s2 := &+[ aa[i]*aa[j] : i, j in [1..4] | i lt j ];
s3 := &+[ aa[i]*aa[j]*aa[k] : i, j, k in [1..4] | i lt j and j lt k ]; s4 := &*aa;
assert s2^2 eq 4*s4;
x0 := s2/2; y0 := (s1*s2 - s3)*x0 - s1*s4;
assert Evaluate(f, x0) eq y0^2;
P8 := elt< J | x - x0, y0 >;
D0 := elt< J | x^2 - s2*x + s4, (s1*s2 - s3)*x - s1*s4 >;
T0 := elt< J | x, 0 >;
assert 2*P8 eq D0 and 4*P8 eq T0 and T0 ne J!0 and 8*P8 eq J!0;
printf "over Q(t): P = (%o, %o) is a rational point; 2[P - oo] = D0, 4[P - oo] = [x, 0] != 0, 8[P - oo] = 0\n", x0, y0;
Ts := [ elt< J | x + aa[i]^2, 0 > : i in [1..3] ];
assert #{ i*P8 + j*Ts[1] + k*Ts[2] + m*Ts[3] : i in [0..7], j, k, m in [0..1] } eq 64;
print "over Q(t): the 64 sums i[P - oo] + j[x+a^2, 0] + k[x+b^2, 0] + m[x+1, 0] are distinct: [2,2,2,8] is marked";
g, st := ResidueCertificate(f, [13, 31], [2, 2], [Q| 2, 405, -401], 64, [2,2,2,8]);
print "FAMILY_2228_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,2,2,8]";
quit;
