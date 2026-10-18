/* family_2_2_12_control_certificate.m -- a second family for [2,2,12], over a rank-1 elliptic
   curve, that contains the census curve 2.2.12.a (control for family_2_2_12_certificate.m).
   It is a full certificate in its own right: infinitely many pairwise non-isomorphic genus-2 curves over Q
   with geometrically simple Jacobian and J(Q)_tors exactly [2,2,12].

   Construction (torsion_jac/notes/simplefam_2026_10_18.md, section 4).  2.2.12.a is y^2 = F1 F2 with
   F1 - F2 constant (so D3 = oo+ - oo- has order 3), F1 with three rational roots, F2 = (linear)(irreducible
   quadratic), and its class of order 4 halves the root pair {root of F1, root of F2}.  Normalize that pair
   to {0, 1}:
       F1 = x (x - e1)(x - e2),   F2 = F1 - kappa = (x - 1) q(x),   kappa = F1(1),   f = F1 F2,
   and halve T01 = [x(x-1), 0].  The 2-descent conditions for T01 are e1(1-e1) = square, e2(1-e2) = square
   (e_i = 1/(1 + k_i^2)), N := 3 k1^2 k2^2 + k1^2 + k2^2 = square (the norm of beta(1-beta), beta a root of q),
   and beta(1-beta) a square in Q(beta) (one more square condition).  2.2.12.a is (k1, k2) = (1/2, 1/3).
   On the line k1 = 1/2 the norm condition is the conic (2m)^2 = 7 k2^2 + 1, parametrized by
   k2 = 2 tau/(7 - tau^2), and the last condition is the genus-1 curve
       C1 : A^2 = tau^4 + tau^3 - 17 tau^2 + 7 tau + 49,
   whose Jacobian E: y^2 = x^3 + x^2 - 285 x + 1900 (conductor 2660) has rank 1.  The point R0 = (-7, 35) of
   C1 gives 2.2.12.a.

   Argument.
   (1) Over the function field Q(C1): D3 = oo+ - oo- has order 3; Ta = [x(x-e1), 0], Tb = [x(x-e2), 0],
       T01 = [x(x-1), 0] have order 2; D4 = [u, v] (computed below from the identity
       gamma^2 h - ell^2 g = u^2, g = x(x-1), h = f/g, v = -ell g/gamma mod u) has 2 D4 = T01.  The 48 sums
       i D3 + j D4 + k Ta + l Tb are distinct: [2,2,12] is marked, and it specializes injectively to every
       smooth fiber.
   (2) E has rank 1 (RankBounds); with origin R0, the point P coming from (tau, A) = (3, -5) has infinite
       order, and order 12 in E(F_17) and 36 in E(F_233).  The quartic is squarefree mod 17 and 233, so for
       R = 36k P the point of C1 reduces to R0, i.e. tau(R) = -7 mod 17 and mod 233.
   (3) The coefficients of f have denominators that are units at tau = -7 mod 17 and mod 233, and f_{-7}
       is squarefree of degree 6 there (ReductionClass).  So C_{tau(R)} has the reduction of 2.2.12.a at 17
       and 233: #J(F_17) = 240, #J(F_233) = 56496, gcd 48 = #[2,2,12], strict at 233.  Hence
       J(Q)_tors = [2,2,12] and J is geometrically simple for every R = 36k P (k >= 1).
   (4) The G2-invariants at tau = -7 and tau = 3 differ, so infinitely many isomorphism classes.
   Run (from the repository root):
       magma -b pipeline/magma/family_2_2_12_control_certificate.m > data/knowledge/sources/family_2_2_12_control_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// the family over Q(tau)
Ftau<tau> := FunctionField(Q); Px<x> := PolynomialRing(Ftau);
k1 := Ftau!(1/2); k2 := 2*tau/(7 - tau^2);
e1 := 1/(1 + k1^2); e2 := 1/(1 + k2^2);
F1 := x*(x - e1)*(x - e2); kap := Evaluate(F1, 1); F2 := F1 - kap;
fT := F1*F2;
assert Degree(fT) eq 6 and IsSquarefree(fT) and LeadingCoefficient(fT) eq 1;
printf "f_tau = x (x - e1)(x - e2) (x (x - e1)(x - e2) - kappa), e1 = 4/5, e2 = %o, kappa = %o\n", e2, kap;

// (1) the marked group over the function field of C1
PA<AA> := PolynomialRing(Ftau);
quartic := tau^4 + tau^3 - 17*tau^2 + 7*tau + 49;
K<A> := FunctionField(AA^2 - quartic);
P<X> := PolynomialRing(K);
f := P![ K!c : c in Coefficients(fT) ];
eK1 := K!e1; eK2 := K!e2; k2K := K!k2;
m := (1 + tau*k2K)/2;
assert m^2 eq 3*(K!k1)^2*k2K^2 + (K!k1)^2 + k2K^2;      // the norm condition holds on the conic
q := (f div (X*(X - eK1)*(X - eK2))) div (X - 1);
g := X*(X - 1); h := (X - eK1)*(X - eK2)*q;
assert f eq g*h and Degree(q) eq 2;
s1 := K!(2/5); s2 := k2K/(1 + k2K^2);
assert s1^2 eq eK1*(1 - eK1) and s2^2 eq eK2*(1 - eK2);
L<be> := ext< K | q >;
ok, sb := IsSquare(be*(1 - be)); assert ok;              // beta(1-beta) is a square in Q(C1)(beta)
sb0 := Eltseq(sb)[1]; sb1 := Eltseq(sb)[2];
q1 := Coefficient(q, 1); q0 := Coefficient(q, 0);
// gamma^2 h - ell^2 g = u^2 forces u(e) = +-s(e) ell(e) at the roots e of h; solve for u = X^2 + u1 X + u0,
// ell = al + bl X (signs (+,+,+); the assertions below check the result)
rows := [ [eK1, K!1, -s1, -s1*eK1], [eK2, K!1, -s2, -s2*eK2],
          [K!0, K!1, -sb0, sb1*q0], [K!1, K!0, -sb1, -(sb0 - sb1*q1)] ];
rhs := [ -eK1^2, -eK2^2, q0, q1 ];
sol := Solution(Transpose(Matrix(K, rows)), Vector(K, rhs));
u1, u0, al, bl := Explode(Eltseq(sol));
u := X^2 + u1*X + u0; ell := al + bl*X;
ok, gam := IsSquare(1 + bl^2); assert ok;
assert gam^2*h - ell^2*g eq u^2;
C := HyperellipticCurve(f); J := Jacobian(C);
D4 := elt< J | u, (-ell*g/gam) mod u >;
T01 := elt< J | g, 0 >;
assert T01 ne J!0 and 2*T01 eq J!0 and 2*D4 eq T01;
D3 := (C![1, 1, 0]) - (C![1, -1, 0]);
assert D3 ne J!0 and 3*D3 eq J!0;
Ta := elt< J | X*(X - eK1), 0 >; Tb := elt< J | X*(X - eK2), 0 >;
assert 2*Ta eq J!0 and 2*Tb eq J!0;
grp := { i*D3 + j*D4 + k*Ta + l*Tb : i in [0..2], j in [0..3], k, l in [0..1] };
assert #grp eq 48;
print "(1) over Q(C1): D3 (order 3), D4 (2 D4 = [x(x-1), 0], order 4), [x(x-e1),0], [x(x-e2),0] generate a group of order 48 = [2,2,12]";

// (2) the elliptic curve and the points 36k P
P2<Zz,W,V> := ProjectiveSpace(Q, [1,2,1]);
Hq := Curve(P2, W^2 - (Zz^4 + Zz^3*V - 17*Zz^2*V^2 + 7*Zz*V^3 + 49*V^4));
R0 := Hq![-7, 35, 1];
E0, phi := EllipticCurve(Hq, R0);
E := MinimalModel(E0);
assert aInvariants(E) eq [0, 1, 0, -285, 1900] and Conductor(E) eq 2660;
r1, r2 := RankBounds(E); assert r1 eq 1 and r2 eq 1;
ok, iso := IsIsomorphic(E0, E); assert ok;
Pt := iso(phi(Hq![3, -5, 1]));
assert Order(Pt) eq 0;
RT<T> := PolynomialRing(Q);
quarticQ := T^4 + T^3 - 17*T^2 + 7*T + 49;
for p in [17, 233] do
    assert IsSquarefree(PolynomialRing(GF(p))!quarticQ);
    assert 36 mod Order(ChangeRing(E, GF(p))!Pt) eq 0;
end for;
printf "(2) E = %o, conductor 2660, rank 1; P = %o has infinite order, order %o mod 17 and %o mod 233\n",
    aInvariants(E), Pt, Order(ChangeRing(E, GF(17))!Pt), Order(ChangeRing(E, GF(233))!Pt);
ok, inv := IsIsomorphic(E, E0); assert ok;
taus := [ pt[1]/pt[3] : pt in RationalPoints(Pullback(phi, inv(36*Pt))) | pt[3] ne 0 ];
assert #taus ge 1 and &and[ GF(17)!tt eq GF(17)!(-7) and GF(233)!tt eq GF(233)!(-7) : tt in taus ];
printf "sanity: tau(36 P) has %o digits and is = -7 mod 17 and mod 233\n", #Sprint(Numerator(taus[1]));

// (3) same good reduction as tau = -7 (census curve 2.2.12.a) at 17 and 233
nJ := []; anystrict := false;
for p in [17, 233] do
    fbar := ReductionClass(fT, p, -7);
    assert PolynomialRing(GF(p))!SpecPoly(fT, taus[1]) eq fbar;
    n, chi := FrobData(fbar); st := StrictChi(chi);
    anystrict := anystrict or st; Append(~nJ, n);
    printf "p = %o, tau = -7 mod %o: reduction squarefree of degree %o, #J(F_%o) = %o, chi = %o%o\n",
        p, p, Degree(fbar), p, n, chi, st select " (strict)" else "";
end for;
assert GCD(nJ) eq 48 and anystrict;
printf "(3) gcd of #J(F_p) = %o = #[2,2,12]\n", GCD(nJ);

// control and sanity
f7 := SpecPoly(fT, -7);
RX2<X2> := PolynomialRing(Q);
Ccensus := HyperellipticCurve(X2^6 - 7*X2^5 + 8*X2^4 + 16*X2^3 - 5*X2^2 - 5*X2, X2^2 + X2);   // 2.2.12.a
assert IsIsomorphic(HyperellipticCurve(f7), Ccensus);
print "control: the fiber at R0 (tau = -7) is Q-isomorphic to census curve 2.2.12.a";
inv := TorsionInvariants(f7); assert inv eq [2,2,12];
inv3 := TorsionInvariants(SpecPoly(fT, 3)); assert inv3 eq [2,2,12];
printf "sanity: tau = -7 and tau = 3 have TorsionSubgroup %o, %o\n", inv, inv3;
// (4) the moduli map is nonconstant
assert G2Invariants(HyperellipticCurve(f7)) ne G2Invariants(HyperellipticCurve(SpecPoly(fT, 3)));
print "(4) G2-invariants at tau = -7 and tau = 3 differ";
print "FAMILY_2_2_12_CONTROL_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,2,12], in a family through census curve 2.2.12.a";
quit;
