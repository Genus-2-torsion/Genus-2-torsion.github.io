/* family_2_4_8_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,4,8].
   (torsion_jac copy: code/simple248_family_2_4_8_certificate.m; notes/simple248_2026_10_19.md.)

   Construction.  The M_1(8,2^w) odd model of torsion_jac/notes/m18_m14_halving.md, on the chart n = 1:
       C : y^2 = f(x) = x A(x) B(x),
       A = x^2 + (R^3 + 4R^2 t + R - 8Rt + 4t) x + R^4,
       B = (R + 2 + 4t) x^2 + (R^2 + 4R + 1 + 8t) x + (2R^2 + R + 4t),
       t = (2R^2 + (1 - w^2) R - 2w^2)/(4(w^2 - 1)),
   where P8 = (-R, y0) - oo has order 8 and T0 = [x, 0] has order 2.  The 2-descent (Cassels map) of T0
   reduces the [4,8] tangent cover to two square conditions,
       (R^2 - 1)(R^2 - w^2) = square,   (R - 1)(w - 1)(Rw + 3R + 3w + 1) = square,
   and B splits over Q iff in addition (R - 1)(w + 1)(Rw - 3R + 3w - 1) = square.  On the line
       w = (13R - 12)/(12R - 13),   R = 3(tau^2 - 4)/(2(tau^2 - 9)),
   the first two conditions hold identically in Q(tau), and the third becomes
       D : sigma^2 = tau^4 - 36 tau^2 + 36,
   a genus-1 curve whose Jacobian E : y^2 = x^3 - 36x (conductor 576) has rank 1.  The point
   R0 = (tau, sigma) = (12/35, 6906/1225) of D gives a geometrically simple member with
   (R, w) = (2378/3627, 2522/3723).  (The census curve 2.4.8.a is the point (R, w) = (-16/11, 14/11),
   where A splits instead of B; it is not on D.)

   Argument.
   (1) Over the function field Q(D) = Q(tau)(sigma): P8 has order 8; the half H = [u, x(Mx + N) mod u]
       given by h - x(Mx + N)^2 = c4 u^2 (h = AB, c4 = lc(h), u = x^2 + Ux + V, V = -R^2 w) satisfies
       2H = T0; TB = [x - beta, 0] (beta a root of B in Q(D)) has order 2.  The 64 sums
       i P8 + j H + k TB are distinct: [2,4,8] is marked over Q(D), and it specializes injectively to
       every smooth fiber.  (Already over Q(tau), P8 and H give [4,8].)
   (2) E has rank 1 (RankBounds); with origin R0, the point P coming from (tau, sigma) = (-1, 1) has
       infinite order, order 6 in E(F_11) and 20 in E(F_53).  The quartic is squarefree mod 11 and 53, so
       for R = 60k P (k >= 1) the point of D reduces to R0, i.e. tau(R) = 12/35 mod 11 and mod 53.
   (3) The coefficients of f have denominators that are units at tau = 6 mod 11 and tau = 17 mod 53, and
       the reductions are squarefree of degree 5 (ReductionClass).  So C_{tau(R)} has the reductions of the
       member at R0: J(F_11) = [2,4,16] (order 128), J(F_53) = [2,4,328] (order 2624), gcd 64 = #[2,4,8],
       strict at 53.  Hence J(Q)_tors = [2,4,8] and J is geometrically simple for every R = 60k P.
   (4) The G2-invariants at tau(R0) and tau(60P) differ, so the moduli map D -> M_2 is nonconstant, and the
       points 60kP give infinitely many isomorphism classes.
   Run (from the census repository root):
       magma -b pipeline/magma/family_2_4_8_certificate.m < /dev/null > data/knowledge/sources/family_2_4_8_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
SetQuitOnError(true);                   // fail closed: a failed assertion stops the run before the marker
T0time := Cputime();
load "pipeline/magma/family_cert_lib.m";

// ---------------------------------------------------------------- (0) the family over Q(tau)
Ft<tau> := FunctionField(Q); Pt<xt> := PolynomialRing(Ft);
Rt := 3*(tau^2 - 4)/(2*(tau^2 - 9));
wt := (13*Rt - 12)/(12*Rt - 13);
tt := (2*Rt^2 + (1 - wt^2)*Rt - 2*wt^2)/(4*(wt^2 - 1));
At := xt^2 + (Rt^3 + 4*Rt^2*tt + Rt - 8*Rt*tt + 4*tt)*xt + Rt^4;
Bt := (Rt + 2 + 4*tt)*xt^2 + (Rt^2 + 4*Rt + 1 + 8*tt)*xt + (2*Rt^2 + Rt + 4*tt);
fT := xt*At*Bt;
assert Degree(fT) eq 5 and IsSquarefree(fT);
// the square conditions of the construction (not used by the proof, which is the Mumford arithmetic below)
assert IsSquare((Rt^2 - 1)*(Rt^2 - wt^2)) and IsSquare((Rt - 1)*(wt - 1)*(Rt*wt + 3*Rt + 3*wt + 1));
assert IsSquare((Rt - 1)*(wt + 1)*(Rt*wt - 3*Rt + 3*wt - 1)/(tau^4 - 36*tau^2 + 36));
printf "(0) C_tau : y^2 = x A B with R = %o, w = (13R - 12)/(12R - 13), t = %o\n", Rt, tt;

// ---------------------------------------------------------------- (1) the marked group
// over Q(tau): P8 and the half H of T0
Jt := Jacobian(HyperellipticCurve(fT));
ok, y0 := IsSquare(Evaluate(fT, -Rt)); assert ok;
P8t := elt< Jt | xt + Rt, Pt!y0 >;
assert 8*P8t eq Jt!0 and 4*P8t ne Jt!0;
T0t := elt< Jt | xt, 0 >;
hT := At*Bt; c := [ Coefficient(hT, i) : i in [0..4] ];
V := -Rt^2*wt;
// the tangent condition as a polynomial in U (written in the variable xt of the same ring Q(tau)[.])
FU := 4*(c[4] - 2*c[5]*xt)*(c[2] - 2*c[5]*xt*V) - (c[3] - c[5]*(xt^2 + 2*V))^2;
found := false;
for rt in Roots(FU) do
    U0 := rt[1];
    okM, M := IsSquare(c[4] - 2*c[5]*U0); okN, N := IsSquare(c[2] - 2*c[5]*U0*V);
    if not (okM and okN) then continue; end if;
    for e in [1, -1] do
        if hT - xt*(M*xt + e*N)^2 eq c[5]*(xt^2 + U0*xt + V)^2 then
            uH := xt^2 + U0*xt + V; vH := (xt*(M*xt + e*N)) mod uH; found := true; break rt;
        end if;
    end for;
end for;
assert found;
Ht := elt< Jt | uH, vH >;
assert T0t ne Jt!0 and 2*T0t eq Jt!0 and 2*Ht eq T0t;
assert #{ i*P8t + j*Ht : i in [0..7], j in [0..3] } eq 32;
printf "(1) over Q(tau): P8 = [x + R, y0] has order 8, H = [u, v] with u = %o has 2H = [x, 0]; <P8, H> = [4,8]\n", uH;
// over Q(D): add the root of B
PA<AA> := PolynomialRing(Ft);
K<sg> := FunctionField(AA^2 - (tau^4 - 36*tau^2 + 36));
PK<x> := PolynomialRing(K);
toK := func< g | PK![ K!a : a in Coefficients(g) ] >;      // coefficient-wise coercion
f := toK(fT); BK := toK(Bt);
J := Jacobian(HyperellipticCurve(f));
P8 := elt< J | toK(xt + Rt), toK(Pt!y0) >;
H := elt< J | toK(uH), toK(vH) >;
T0 := elt< J | x, 0 >;
rB := Roots(BK); assert #rB eq 2;
TB := elt< J | x - rB[1][1], 0 >;
assert 8*P8 eq J!0 and 4*P8 ne J!0 and 2*H eq T0 and T0 ne J!0 and 2*T0 eq J!0;
assert TB ne J!0 and 2*TB eq J!0;
grp := { i*P8 + j*H + k*TB : i in [0..7], j in [0..3], k in [0..1] };
assert #grp eq 64;
print "    over Q(D), D : sigma^2 = tau^4 - 36 tau^2 + 36: B splits, TB = [x - beta, 0] has order 2, and the 64 sums i P8 + j H + k TB are distinct: [2,4,8] is marked";

// ---------------------------------------------------------------- (2) the elliptic base and the points 60k P
P2<Zz,W,Vv> := ProjectiveSpace(Q, [1,2,1]);
Hq := Curve(P2, W^2 - (Zz^4 - 36*Zz^2*Vv^2 + 36*Vv^4));
R0 := Hq![12, 6906, 35];
E0, phi := EllipticCurve(Hq, R0);
E := MinimalModel(E0);
assert aInvariants(E) eq [0, 0, 0, -36, 0] and Conductor(E) eq 576;
r1, r2 := RankBounds(E); assert r1 eq 1 and r2 eq 1;
ok, iso := IsIsomorphic(E0, E); assert ok;
Pe := iso(phi(Hq![-1, 1, 1]));
assert Order(Pe) eq 0;
RT<T> := PolynomialRing(Q);
quarticQ := T^4 - 36*T^2 + 36;
primes := [11, 53]; tau0 := Q!12/35;
for p in primes do
    assert IsSquarefree(PolynomialRing(GF(p))!quarticQ) and Valuation(tau0, p) eq 0;
    assert 60 mod Order(ChangeRing(E, GF(p))!Pe) eq 0;
end for;
printf "(2) E = %o, conductor 576, rank 1; P = %o has infinite order, order %o mod 11 and %o mod 53\n",
    aInvariants(E), Pe, Order(ChangeRing(E, GF(11))!Pe), Order(ChangeRing(E, GF(53))!Pe);
ok, inv := IsIsomorphic(E, E0); assert ok;
// Pullback also returns a base point of the defining equations of phi (tau = -12/35, phi(q) != target):
// keep only the points that phi really maps to 60 P
tgt := inv(60*Pe);
pts60 := [];
for q in RationalPoints(Pullback(phi, tgt)) do
    good := false;
    try good := q[3] ne 0 and phi(q) eq tgt; catch err good := false; end try;
    if good then Append(~pts60, q); end if;
end for;
assert #pts60 eq 1;
tau1 := pts60[1][1]/pts60[1][3];
assert &and[ Valuation(tau1, p) ge 0 and GF(p)!tau1 eq GF(p)!tau0 : p in primes ];
assert IsSquare(Evaluate(quarticQ, tau1));
printf "    sanity: tau(60 P) has %o digits, lies on D, and is = 12/35 mod 11 and mod 53\n", #Sprint(Numerator(tau1));

// ---------------------------------------------------------------- (3) same good reduction as tau0 at 11 and 53
nJ := []; anystrict := false;
for p in primes do
    r0 := Integers()!(GF(p)!tau0);
    fbar := ReductionClass(fT, p, r0);
    assert PolynomialRing(GF(p))!SpecPoly(fT, tau0) eq fbar;
    assert PolynomialRing(GF(p))!SpecPoly(fT, tau1) eq fbar;
    n, chi := FrobData(fbar); st := StrictChi(chi);
    invp := Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
    assert &*invp eq n;
    anystrict := anystrict or st; Append(~nJ, n);
    printf "(3) p = %o, tau = %o mod %o: reduction %o (squarefree, degree %o), J(F_%o) = %o, #J(F_%o) = %o, chi = %o%o\n",
        p, r0, p, fbar, Degree(fbar), p, invp, p, n, chi, st select " (strict)" else "";
end for;
assert GCD(nJ) eq 64 and anystrict;
printf "    gcd of #J(F_p) = %o = #[2,4,8]: every R = 60k P (k >= 1) has J(Q)_tors exactly [2,4,8], geometrically simple\n", GCD(nJ);

// ---------------------------------------------------------------- (4) nonconstant moduli map, (5) sanity
f0 := SpecPoly(fT, tau0); f1 := SpecPoly(fT, tau1);
assert G2Invariants(HyperellipticCurve(f0)) ne G2Invariants(HyperellipticCurve(f1));
print "(4) the G2-invariants at tau(R0) = 12/35 and tau(60 P) differ";
printf "(5) the member at R0 (tau = 12/35, R = %o, w = %o): y^2 = %o\n", Evaluate(Rt, tau0), Evaluate(wt, tau0), f0;
inv0 := TorsionInvariants(f0); assert inv0 eq [2, 4, 8];
printf "    TorsionSubgroup at tau = 12/35: %o\n", inv0;
inv6 := TorsionInvariants(SpecPoly(fT, 6)); assert inv6 eq [2, 4, 8];
printf "    TorsionSubgroup at tau = 6 ((R, w) = (16/9, 4/3), a point of D off the class): %o\n", inv6;
printf "time %o s\n", Cputime(T0time);
print "FAMILY_2_4_8_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,4,8], over the rank-1 elliptic curve D : sigma^2 = tau^4 - 36 tau^2 + 36";
quit;
