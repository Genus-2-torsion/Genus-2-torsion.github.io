/* family_4_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [4].

   The record's family (census_infinity_record.md, section 2.4; torsion_jac/code/claude_census_infty_round2b.m,
   part A): with q_al(x) = x^2 + (2 al + 1) x + al^2 = (x + al)^2 + x and al = 1 + s sqrt2,
       Q4 = q_al * q_alc (the norm from Q(sqrt2)(s) to Q(s); alc = 1 - s sqrt2),
       C_s : y^2 = f_s = x Q4(x)   (degree 5).
   Argument.
   (1) Halving identity.  From (a^2 + x b^2)(c^2 + x d^2) = (ac - x bd)^2 + x (ad + bc)^2:
           Q4 = u^2 + x l^2,  u = (x + al)(x + alc) - x = x^2 + x + 1 - 2s^2,  l = 2x + 2.
       So f - (x l)^2 = x u^2, div(y - x l) = (0,0) + 2(P1 + P2) - 5 oo, and D = [u, x l mod u] satisfies
       2D = [x, 0], a class of order 2.  Checked by Mumford arithmetic over Q(s): D has order 4, so [4]
       is marked.
   (2) Residue certificate (family_cert_lib.m) at s = 1 mod 3 and s = 1 mod 19, i.e. every 3- and
       19-integral s = 1 mod 57: #J(F_3) = 4, #J(F_19) = 412, gcd 4 = #[4]; the Frobenius at 19 is strict.
       The record's fiber s = 1 lies in the class and has TorsionSubgroup [4].
   (3) The G2-invariants of the members differ, so s -> M_2 is nonconstant.
   (The library's log lines call the family parameter t; here it is s.)
   Run (from the repository root):
       magma -b pipeline/magma/family_4_certificate.m > data/knowledge/sources/family_4_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Fs<s> := FunctionField(Q); P<x> := PolynomialRing(Fs);
// the norm construction of the record, over L = Q(s)(sqrt2)
L<w> := ext< Fs | x^2 - 2 >;
PL<X> := PolynomialRing(L);
al := 1 + s*w; alc := 1 - s*w;
Q4L := (X^2 + (2*al + 1)*X + al^2)*(X^2 + (2*alc + 1)*X + alc^2);
assert &and[ Eltseq(c)[2] eq 0 : c in Coefficients(Q4L) ];
Q4 := P![ Eltseq(c)[1] : c in Coefficients(Q4L) ];
assert Q4 eq x^4 + 6*x^3 + (11 - 4*s^2)*x^2 + (6 - 4*s^2)*x + (1 - 2*s^2)^2;
f := x*Q4;
printf "f_s = %o\n", f;
u := x^2 + x + 1 - 2*s^2; l := 2*x + 2;
assert Q4 eq u^2 + x*l^2;
J := Jacobian(HyperellipticCurve(f));
D := elt< J | u, (x*l) mod u >;
T := elt< J | x, 0 >;
assert T ne J!0 and 2*D eq T and 4*D eq J!0;
print "over Q(s): Q4 = u^2 + x l^2 with u = x^2 + x + 1 - 2s^2, l = 2x + 2; D = [u, x l mod u] has 2D = [x, 0] != 0: [4] is marked";
g, st := ResidueCertificate(f, [3, 19], [1, 1], [Q| 1, 58, -56, 115], 4, [4]);
print "FAMILY_4_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [4]";
quit;
