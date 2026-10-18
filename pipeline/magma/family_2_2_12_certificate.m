/* family_2_2_12_certificate.m -- certificate that infinitely many pairwise non-isomorphic
   genus-2 curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,2,12].

   The family (torsion_jac/notes/simplefam_2026_10_18.md).  As in BNSS2026, put
       F2 = (x - a^2)(x - b^2)(x - c^2),   F1 = x (x^2 - S x + w^2),   S = a^2 + b^2 + c^2,
       w^2 = a^2 b^2 + a^2 c^2 + b^2 c^2,
   so that F1 - F2 = (abc)^2 is constant and D3 = oo+ - oo- has order 3 on y^2 = f = F1 F2.  BNSS2026 take
   both roots s, t of g = x^2 - S x + w^2 to be squares, which gives [2,2,2,12].  A 2-descent
   computation suggests that the class T_st = [g, 0] is already divisible by 2 in J(Q) when only S + 2w is
   a square (one sign): the descent conditions at the roots 0, a^2, b^2, c^2 of F1 F2 / g hold since
   g(0) = w^2, g(a^2) = (bc)^2, g(b^2) = (ac)^2, g(c^2) = (ab)^2, and at the roots of g they say that s is a
   square in Q(s).  On the family below this is checked directly: D4 is an explicit half of T_st over Q(m).
   With g irreducible the rational 2-rank is 3 (factor type [1,1,1,1,2]), and J(Q) contains
   [2,2,4] x Z/3 = [2,2,12].
   The conditions w^2 = a^2b^2 + a^2c^2 + b^2c^2, G^2 = S + 2w describe the quartic surface
   a^4 + b^4 + c^4 + G^4 = 2 (sum of the six products of two squares) (Descartes' relation for the
   squares a^2, b^2, c^2, G^2).  It contains the conic {G = a + b, c^2 = 4(a^2 + ab + b^2)}, parametrized by
       a = m^2 - 1,  b = 2m + 1,  c = 2(m^2 + m + 1),  w = -(2a^2 + ab + 2b^2),  G = a + b.
   This conic is the family.

   Argument.
   (1) Over Q(m): D3 = oo+ - oo- has order 3; [x(x-a^2), 0], [x(x-b^2), 0] and T_st = [g, 0] have order 2;
       D4 = [u, v] below (from the identity ell^2 g - u^2 = gamma^2 x(x-a^2)(x-b^2)(x-c^2), u monic
       quadratic, ell linear, v = -ell g / gamma mod u) has 2 D4 = T_st, so order 4.  The 48 sums
       i D3 + j D4 + k [x(x-a^2),0] + l [x(x-b^2),0] are distinct: [2,2,12] is marked.
   (2) Residue certificate (family_cert_lib.m) at m = 3 mod 23 and m = 3 mod 59 (m = 3 mod 1357):
       #J(F_23) = 528, #J(F_59) = 3984, gcd 48 = #G, strict at 59.  So every 23- and 59-integral
       m = 3 mod 1357 gives J(Q)_tors exactly [2,2,12] and a geometrically simple Jacobian.
   (3) The G2-invariants of two members differ, so infinitely many isomorphism classes.
   Run (from the repository root):
       magma -b pipeline/magma/family_2_2_12_certificate.m > data/knowledge/sources/family_2_2_12_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// (0) the conic lies on the Descartes-type quartic surface
Fm<m> := FunctionField(Q); P<x> := PolynomialRing(Fm);
a := m^2 - 1; b := 2*m + 1; c := 2*(m^2 + m + 1);
S := a^2 + b^2 + c^2; w := -(2*a^2 + a*b + 2*b^2); G := a + b;
assert w^2 eq a^2*b^2 + a^2*c^2 + b^2*c^2 and G^2 eq S + 2*w;
v4 := [a, b, c, G];
assert &+[ e^4 : e in v4 ] eq 2*&+[ v4[i]^2*v4[j]^2 : i, j in [1..4] | i lt j ];
print "over Q(m): w^2 = a^2b^2 + a^2c^2 + b^2c^2 and S + 2w = (a+b)^2";

// (1) the curve and the marked group over Q(m)
g := x^2 - S*x + w^2;
hq := x*(x - a^2)*(x - b^2)*(x - c^2);
F1 := x*g; F2 := (x - a^2)*(x - b^2)*(x - c^2);
assert F1 - F2 eq (a*b*c)^2;
f := F1*F2;
assert Degree(f) eq 6 and IsSquarefree(f) and LeadingCoefficient(f) eq 1;
assert [ <Degree(e[1]), e[2]> : e in Factorization(f) ] eq [ <1,1>, <1,1>, <1,1>, <1,1>, <2,1> ];
printf "f_m = x (x^2 - S x + w^2)(x - a^2)(x - b^2)(x - c^2) with a = %o, b = %o, c = %o, S = %o, w = %o\n", a, b, c, S, w;
C := HyperellipticCurve(f); J := Jacobian(C);
D3 := (C![1, 1, 0]) - (C![1, -1, 0]);
assert D3 ne J!0 and 3*D3 eq J!0;
// the halving identity
u := x^2 - (3*m^4 + 5*m^3 + 9*m^2 + 6*m + 4)*x - (2*m^7 + 2*m^6 + 5*m^5 + 4*m^4 + m^3 - 5*m^2 - 6*m - 3);
ell := -(m^4 + 2*m^3 + 4*m^2 + 2)/(2*(m^3 - 1))*x + m^3 - 1;
gam := (m^4 + 2*m^3 + 2*m^2 + 4*m)/(2*(m^3 - 1));
assert ell^2*g - u^2 eq gam^2*hq;
// (how u, ell were found: g(e) = w^2, (bc)^2, (ac)^2, (ab)^2 at e = 0, a^2, b^2, c^2, and the identity
//  forces u(e) = +-sqrt(g(e)) ell(e) there; the sign pattern (+,+,+,-) gives lc(ell)^2 - 1 = gamma^2)
es := [Fm!0, a^2, b^2, c^2]; sq := [w, b*c, a*c, a*b]; sg := [1, 1, 1, -1];
assert &and[ Evaluate(g, es[i]) eq sq[i]^2 and Evaluate(u, es[i]) eq sg[i]*sq[i]*Evaluate(ell, es[i]) : i in [1..4] ];
v := (-ell*g/gam) mod u;
D4 := elt< J | u, v >;
Tst := elt< J | g, 0 >;
assert Tst ne J!0 and 2*Tst eq J!0 and 2*D4 eq Tst;
printf "ell = %o, gamma = %o\nD4 = [%o, %o] has 2 D4 = [x^2 - S x + w^2, 0] != 0, so order 4\n", ell, gam, u, v;
Ta := elt< J | x*(x - a^2), 0 >; Tb := elt< J | x*(x - b^2), 0 >;
assert 2*Ta eq J!0 and 2*Tb eq J!0;
grp := { i*D3 + j*D4 + k*Ta + l*Tb : i in [0..2], j in [0..3], k, l in [0..1] };
assert #grp eq 48;
print "over Q(m): D3 (order 3), D4 (order 4), [x(x-a^2),0], [x(x-b^2),0] generate a group of order 48 = Z/3 x Z/4 x Z/2 x Z/2 = [2,2,12]";

// (2), (3) residue certificate at m = 3 mod 23, m = 3 mod 59 (m = 3 mod 1357)
gg, st := ResidueCertificate(f, [23, 59], [3, 3], [3, -1354, 1360, 2717], 48, [2,2,12]);
printf "sanity: m = %o has TorsionSubgroup %o\n", -1354, TorsionInvariants(SpecPoly(f, -1354));
print "FAMILY_2_2_12_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,2,12]";
quit;
