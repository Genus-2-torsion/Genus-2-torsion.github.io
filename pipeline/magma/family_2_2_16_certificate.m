/* family_2_2_16_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,2,16].

   Construction (record: torsion_jac/notes/simple2216_2026_11_14.md; search scripts fam16h*.m, sections*.m, ratfam*.m).
   Quintic model y^2 = f(x), W = oo, P = (0, A(0)), T0 = (1, 0), Q = [P - W], H = [R1 + R2 - 2W].
   (I)  f = A^2  - x^4 (x - 1),  A(1) = 0:        div(y - A)  = 4P + T0 - 5W,          so 4Q = [T0 - W];
   (II) f = A'^2 - x q(x)^2,     A'(0) = -A(0):   div(y - A') = 2R1 + 2R2 + iota(P) - 5W, so 2H = Q.
   Hence H has order 16 (8H = [T0 - W]).  This is the configuration of the census curve 2.2.16.b
   (Sutherland): there the order-8 class [P - W] is halvable and 4[P - W] is a rational Weierstrass class.
   (It is not on the M_1(8,2,2) chart of Genus2TorsionNotebook/notes/claude_ov_2216_2026_07_25.md, whose
   4Q is the class of the quadratic pair.)  Writing A + A' = x L1, A - A' = M gives M L1 = -K,
   K = q^2 - x^4 + x^3; K has a rational root rho iff rho(rho - 1) is a square; A(1) = 0 is then automatic
   up to scaling y.  This gives a rational surface of [16]-curves, parameters (tau, q1):
       rho = 1/(1 - tau^2), sigma = tau rho, q = x^2 + q1 x + rho(sigma - rho - q1).
   A further rational root x1 of f/(x - 1) (with q(x1) != 0) must be x1 = xi^2, x1 - 1 = eta^2 (by (I), (II)),
   xi = (s^2 + 1)/(2s), eta = (s^2 - 1)/(2s), and the root condition is LINEAR in q1: q1 = Phi(tau, s)
   (below).  So [2,16] is a rational surface (tau, s), and [2,2,16] is the surface
   S : Phi(tau, s1) = Phi(tau, s2), quadratic in tau, a double cover of P^1 x P^1 (s1, s2) branched along a
   (6,6) curve with two conjugate quadruple points (i, -i), (-i, i).  The pencil of (1,1)-curves
   lambda = (1 - s1 s2)/(s1 + s2) through them makes S an elliptic surface over the lambda-line (genus-1
   fibres w^2 = quartic in s1); its generic fibre has Mordell-Weil rank >= 2 over Q(lambda) (sections_rank.m:
   the sections s1 = 0 and s1 = 1 specialize at lambda = 1/2 to points independent in E(Q)/2E(Q) together
   with E[2]).  The fibre lambda = 1/2 through 2.2.16.b is 25200bf2, rank 2 (2-descent).  The family below
   is the section 2P, P = (s1, w) = (1, 2(lambda-1)^2(lambda+1)), origin (0, (lambda+1)^2(lambda-1)):
       s1 = (lambda - 1)^2 (lambda + 1)/(2 lambda),  s2 = (1 - lambda s1)/(s1 + lambda),
       tau = -4 lambda^2 (lambda^3 - lambda^2 + lambda + 1)/(lambda^8 - 2 lambda^7 + 2 lambda^5 + 2 lambda^4
             + 2 lambda^3 + 4 lambda^2 - 2 lambda + 1),   q1 = Phi(tau, s1) = Phi(tau, s2).
   The census curve 2.2.16.b is the member lambda = 2.

   Argument.
   (1) Over Q(lambda): the identities (I), (II) hold; f = -(x - 1)(x - x1)(x - x2) q2 with q2 irreducible,
       x_i = xi(s_i)^2.  H = [q, A' mod q] has 2H = Q, 4Q = T0 = [x - 1, 0] != 0, so H has order 16;
       T1 = [x - x1, 0], T2 = [x - x2, 0] have order 2, and the 64 sums i H + j T1 + k T2 are distinct:
       [2,2,16] is marked over Q(lambda), and it specializes injectively to every smooth fibre.
   (2) Residue certificate (family_cert_lib.m) at lambda = 12 mod 31 and lambda = 8 mod 37
       (lambda = 415 mod 1147): J(F_31) = [2,2,272] (order 1088) and J(F_37) = [2,2,304] (order 1216),
       gcd 64 = #[2,2,16], strict at 31 and at 37.  So every 31- and 37-integral lambda = 415 mod 1147
       gives a geometrically simple Jacobian with J(Q)_tors exactly [2,2,16].
   (3) The G2-invariants of members of the class differ, so infinitely many isomorphism classes.
   Run (from the repository root):
       magma -b pipeline/magma/family_2_2_16_certificate.m < /dev/null > data/knowledge/sources/family_2_2_16_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
SetQuitOnError(true);                   // fail closed: a failed assertion stops the run before the marker
T0time := Cputime();
load "pipeline/magma/family_cert_lib.m";

// ---------------------------------------------------------------- (0) the family over Q(lambda)
Fl<lam> := FunctionField(Q); P<x> := PolynomialRing(Fl);
Phi := func< tt, ss | (tt^3*ss^4 + 2*tt^3*ss^2 + tt^3 + 2*tt^2*ss^2 + 2*tt^2 - tt*ss^4 + tt*ss^3 - 3*tt*ss^2 - tt*ss - ss^3 - ss^2 - ss - 1)
                      / (tt*(tt + 1)*(ss - 1)*(tt*ss^2 + tt - ss^2 + 1)) >;
PhiExp := func< tt, ss | (tt^3*ss^4 + 2*tt^3*ss^2 + tt^3 + 2*tt^2*ss^2 + 2*tt^2 - tt*ss^4 + tt*ss^3 - 3*tt*ss^2 - tt*ss - ss^3 - ss^2 - ss - 1)
                      / (tt^3*ss^3 - tt^3*ss^2 + tt^3*ss - tt^3 + 2*tt^2*ss - 2*tt^2 - tt*ss^3 + tt*ss^2 + tt*ss - tt) >;
F2<TT,SS> := FunctionField(Q, 2);
assert Phi(TT, SS) eq PhiExp(TT, SS);          // the factored denominator (no unit lost)
s1 := (lam - 1)^2*(lam + 1)/(2*lam);
s2 := (1 - lam*s1)/(s1 + lam);
assert lam eq (1 - s1*s2)/(s1 + s2);
tau := -4*lam^2*(lam^3 - lam^2 + lam + 1)/(lam^8 - 2*lam^7 + 2*lam^5 + 2*lam^4 + 2*lam^3 + 4*lam^2 - 2*lam + 1);
rho := 1/(1 - tau^2); sigma := tau*rho;
assert sigma^2 eq rho*(rho - 1);
q1 := Phi(tau, s1);
assert q1 eq Phi(tau, s2);
q0 := rho*(sigma - rho - q1);
q := x^2 + q1*x + q0;
K := q^2 - x^4 + x^3;
assert Degree(K) eq 3 and Evaluate(K, rho) eq 0;
N := K div (x - rho); assert N*(x - rho) eq K;
c := (rho - 1)/Evaluate(N, 1);
M := c*N; L1 := x - rho;
A0 := (M + x*L1)/2; Ap0 := (x*L1 - M)/2;
ok, gam := IsSquare(-c); assert ok;
A := A0/gam; Ap := Ap0/gam;
f := A^2 - x^4*(x - 1);
assert f eq Ap^2 - x*q^2;                                         // (I) and (II)
assert Degree(A) le 2 and Degree(Ap) le 2 and Evaluate(A, 1) eq 0 and Evaluate(Ap, 0) eq -Evaluate(A, 0) and Evaluate(A, 0) ne 0;
xi1 := (s1^2 + 1)/(2*s1); xi2 := (s2^2 + 1)/(2*s2);
x1 := xi1^2; x2 := xi2^2;
assert Evaluate(f, 1) eq 0 and Evaluate(f, x1) eq 0 and Evaluate(f, x2) eq 0;
assert IsSquare(x1 - 1) and IsSquare(x2 - 1);
assert Degree(f) eq 5 and IsSquarefree(f) and LeadingCoefficient(f) eq -1;
fac := Factorization(f);
assert Sort([ Degree(t[1]) : t in fac ]) eq [1, 1, 1, 2] and &and[ t[2] eq 1 : t in fac ];
q2 := [ t[1] : t in fac | Degree(t[1]) eq 2 ][1];
assert f eq -(x - 1)*(x - x1)*(x - x2)*q2;
printf "(0) f_lambda = A^2 - x^4 (x - 1) = A'^2 - x q^2 with A(1) = 0, A'(0) = -A(0); f = -(x - 1)(x - x1)(x - x2) q2,\n";
printf "    x1 = %o,\n    x2 = %o\n", x1, x2;
printf "    max degree in lambda of the coefficients of f: %o\n", Max([ Max(Degree(Numerator(cc)), Degree(Denominator(cc))) : cc in Coefficients(f) ]);

// ---------------------------------------------------------------- (1) the marked group over Q(lambda)
J := Jacobian(HyperellipticCurve(f));
Qc := elt< J | x, Evaluate(A, 0) >;          // [P - W], P = (0, A(0))
H := elt< J | q, Ap mod q >;                 // [R1 + R2 - 2W], R1, R2 over the roots of q
Tw0 := elt< J | x - 1, 0 >; Tw1 := elt< J | x - x1, 0 >; Tw2 := elt< J | x - x2, 0 >;
assert 2*H eq Qc;
assert 4*Qc eq Tw0 and Tw0 ne J!0 and 2*Tw0 eq J!0;
assert 16*H eq J!0 and 8*H ne J!0 and 8*H eq Tw0;
assert Tw1 ne J!0 and 2*Tw1 eq J!0 and Tw2 ne J!0 and 2*Tw2 eq J!0;
grp := { i*H + j*Tw1 + k*Tw2 : i in [0..15], j in [0..1], k in [0..1] };
assert #grp eq 64;
print "(1) over Q(lambda): 2H = Q, 4Q = [x - 1, 0] != 0, so H has order 16; [x - x1, 0], [x - x2, 0] have order 2;";
print "    the 64 sums i H + j T1 + k T2 are distinct: [2,2,16] is marked over Q(lambda)";

// ---------------------------------------------------------------- (2) residue certificate
members := [Q| -23/11, 52/25, 29/36];
g, st := ResidueCertificate(f, [31, 37], [12, 8], members, 64, [2,2,16]);
for i in [1..2] do
    p := [31, 37][i]; r0 := [12, 8][i];
    fbar := ReductionClass(f, p, r0);
    invp := Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
    printf "    p = %o: J(F_%o) = %o; reduction factors as %o\n", p, p, invp, [ Degree(t[1]) : t in Factorization(fbar) ];
end for;
printf "(2) every 31- and 37-integral lambda = 415 mod 1147 has J(Q)_tors exactly [2,2,16] and a geometrically simple Jacobian\n";

// ---------------------------------------------------------------- (3) sanity
for m in members[2..#members] do
    inv := TorsionInvariants(SpecPoly(f, m));
    printf "sanity: lambda = %o has TorsionSubgroup %o\n", m, inv;
    assert inv eq [2,2,16];
end for;
RX2<X2> := PolynomialRing(Q);
fS := 360000*X2^6 - 11736600*X2^5 + 128899980*X2^4 - 496170642*X2^3 + 166739244*X2^2 + 519681459*X2 - 284225007;
CS := HyperellipticCurve(fS, X2^2 + X2);
assert IsIsomorphic(HyperellipticCurve(SpecPoly(f, 2)), CS);
print "(3) the member lambda = 2 (outside the class) is isomorphic to the census curve 2.2.16.b";
printf "time %o s\n", Cputime(T0time);
print "FAMILY_2_2_16_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,2,16]";
quit;
