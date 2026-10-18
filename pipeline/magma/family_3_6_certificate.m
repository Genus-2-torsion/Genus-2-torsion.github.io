/* family_3_6_certificate.m -- certificate that infinitely many pairwise non-isomorphic
   genus-2 curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [3,6].

   The family (torsion_jac/notes/simplefam_2026_10_18.md).  Bruin-Flynn-Testa, Acta Arith. 165 (2014),
   Theorem 6 (as in family_33_certificate.m): with parameters r, s, t,
       H1 = x^2 + r x + t,  lambda1 = 4s,     H2 = x^2 + x + r,  lambda2 = 4st,
       G1 = (s - st - 1)x^3 + 3s(r-t)x^2 + 3sr(r-t)x - st^2 + sr^3 + t,
       G2 = (s - st + 1)x^3 + 3s(r-t)x^2 + 3sr(r-t)x - st^2 + sr^3 - t,
       F = G1^2 + lambda1 H1^3 = G2^2 + lambda2 H2^3,
   and T1 = [H1, G1], T2 = [H2, G2] span (Z/3)^2.
   - The coefficient of x^6 in F is s^2(1-t)^2 + 2s(1+t) + 1; it vanishes on t = tau^2,
     s = -1/(1+tau)^2, where F is a quintic (rational Weierstrass point at infinity).
   - On that surface, F(x0) = 0 for a finite x0 forces H1(x0), H2(x0) to be squares; writing
     r = alpha2^2 - x0^2 - x0 and tau = (C/lambda - lambda)/2 with C = x0(alpha2^2 - x0^2), one
     component of the condition is lambda(x0-2alpha2)(x0-alpha2)(x0+alpha2) - x0(x0-alpha2)^2
     = lambda^2(lambda-1), which contains the line x0 = 2l - 1, alpha2 = l - 1, lambda = l.
   On this line
       (r, s, t) = (1 - 3l^2, -1/(9 l^2 (l-1)^2), (3l^2 - 3l + 1)^2),
   F is a quintic with the rational root x0 = 2l - 1 and an irreducible quartic cofactor: factor type
   [1,4] of a quintic, rational 2-rank 1.  So J(Q) contains (Z/3)^2 x Z/2 = [3,6].

   Argument.
   (1) Over Q(l): T1 = [H1, G1 mod H1], T2 = [H2, G2 mod H2] have order 3, T2 is not in <T1>, and
       T0 = [x - x0, 0] has order 2.  The 18 sums are distinct: [3,6] is marked.
   (2) Residue certificate (family_cert_lib.m) at l = 2 mod 11 and l = 2 mod 37 (l = 2 mod 407):
       #J(F_11) = 126, #J(F_37) = 1746, gcd 18 = #G, strict at 37.  So every 11- and 37-integral
       l = 2 mod 407 gives J(Q)_tors exactly [3,6] and a geometrically simple Jacobian.
   (3) The G2-invariants of two members differ, so infinitely many isomorphism classes.
   No census control: the census curve 3.6.a has no rational Weierstrass point, so it is not a member.
   Run (from the repository root):
       magma -b pipeline/magma/family_3_6_certificate.m > data/knowledge/sources/family_3_6_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// (0) the BFT identity over Q(r,s,t), and the quintic locus
F3<r3,s3,t3> := FunctionField(Q, 3); P3<x3> := PolynomialRing(F3);
BFT := func< r, s, tt, xx | [ xx^2 + r*xx + tt, 4*s, xx^2 + xx + r, 4*s*tt,
    (s - s*tt - 1)*xx^3 + 3*s*(r-tt)*xx^2 + 3*s*r*(r-tt)*xx - s*tt^2 + s*r^3 + tt,
    (s - s*tt + 1)*xx^3 + 3*s*(r-tt)*xx^2 + 3*s*r*(r-tt)*xx - s*tt^2 + s*r^3 - tt ] >;
H1, l1, H2, l2, G1, G2 := Explode(BFT(r3, s3, t3, x3));
assert G1^2 + l1*H1^3 eq G2^2 + l2*H2^3;
assert Coefficient(G1^2 + l1*H1^3, 6) eq s3^2*(1-t3)^2 + 2*s3*(1+t3) + 1;
print "over Q(r,s,t): G1^2 + 4s H1^3 = G2^2 + 4st H2^3 (BFT Theorem 6); x^6-coefficient s^2(1-t)^2 + 2s(1+t) + 1";

// (1) the line and the marked group over Q(l)
Fl<l> := FunctionField(Q); P<x> := PolynomialRing(Fl);
rl := 1 - 3*l^2; sl := -1/(9*l^2*(l-1)^2); tl := (3*l^2 - 3*l + 1)^2;
// derivation trace: tau = -(3l^2 - 3l + 1), x0 = 2l - 1, alpha2 = l - 1, lambda = l
tau := -(3*l^2 - 3*l + 1); al2 := l - 1; xx0 := 2*l - 1; Cl := xx0*(al2^2 - xx0^2);
assert tl eq tau^2 and sl eq -1/(1 + tau)^2 and rl eq al2^2 - xx0^2 - xx0 and tau eq (Cl/l - l)/2;
H1, l1, H2, l2, G1, G2 := Explode(BFT(rl, sl, tl, x));
f := G1^2 + l1*H1^3;
assert f eq G2^2 + l2*H2^3 and Degree(f) eq 5;
x0 := 2*l - 1;
assert Evaluate(f, x0) eq 0;
fac := Factorization(f);
assert [ <Degree(e[1]), e[2]> : e in fac ] eq [ <1,1>, <4,1> ] and fac[1][1] eq x - x0;
printf "(r, s, t) = (%o, %o, %o); f_l has the root x0 = %o and an irreducible quartic cofactor\n", rl, sl, tl, x0;
J := Jacobian(HyperellipticCurve(f));
T1 := elt< J | H1, G1 mod H1 >;
T2 := elt< J | H2, G2 mod H2 >;
T0 := elt< J | x - x0, 0 >;
assert T1 ne J!0 and 3*T1 eq J!0 and T2 ne J!0 and 3*T2 eq J!0;
assert T2 ne T1 and T2 ne -T1;
assert T0 ne J!0 and 2*T0 eq J!0;
grp := { i*T1 + j*T2 + k*T0 : i, j in [0..2], k in [0..1] };
assert #grp eq 18;
print "over Q(l): T1 = [H1, G1], T2 = [H2, G2] (order 3, independent) and T0 = [x - x0, 0] (order 2) generate a group of order 18 = [3,6]";

// (2), (3) residue certificate at l = 2 mod 11, l = 2 mod 37 (l = 2 mod 407)
g, st := ResidueCertificate(f, [11, 37], [2, 2], [2, 409, -405, 816], 18, [3,6]);
for m in [409, -405] do
    printf "sanity: l = %o has TorsionSubgroup %o\n", m, TorsionInvariants(SpecPoly(f, m));
end for;
print "FAMILY_3_6_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [3,6]";
quit;
