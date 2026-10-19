/* family_35_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E x F_R) with J(Q)_tors exactly [35].
   Found 2026-10-19 (torsion_jac notes/hlpEC_2026_10_19.md).

   Source.  HLP 2000, section 3.6 (the [35] entry of Table 1): the (2,2)-gluing (HLP Prop. 4) of the fixed curve
       E = E_7(-1) = 26b1 : y^2 + (1 - c)xy - by = x^3 - bx^2, b = -2, c = 2   (torsion Z/7, P_E = (0,0))
   with Kubert X_1(5) curves F_t = E_5(t) (b = c = t, (0,0) of order 5) along a Galois-module isomorphism
   F_t[2] -> E[2].  Neither curve has a rational 2-torsion point.  The pairs (t, psi) form a genus-1 twist X' of
   X_1(2,10); HLP show that it is 3-isogenous to X'' : y^2 = -26 t(t^2 - 11t - 1) (conductor 54080, rank 2).
   Here X' is rebuilt explicitly: it is the curve of Tschirnhausen data (c1, c2) with charpoly(phi(a)) equal to the
   depressed 2-division cubic of F_t, and it is isomorphic to
       EM : y^2 = x^3 + x^2 - 901x - 7461   (conductor 54080, rank 2, MW = Z/2 x Z^2, e.g. (-9,0), (91,820), (-18,-57)),
   the 3-isogenous curve of X'' (X' -> X'' is the dual isogeny).  The functions t, c1, c2 on EM (below, as elements
   A(x) + B(x) y of Q(EM)) were found by factoring the sextic resolvent over Q(EM).  HLP's point t = 1/26 is the image of O.

   The gluing in closed form.  Let E : y^2 = f(z), f = z^3 + pE z + qE with roots a1, a2, a3, and let
   phi(z) = c1 z + c2 (z^2 + 2 pE/3) be a Tschirnhausen transformation with c2 != 0.  Then
   F : y^2 = g(z), g = charpoly of phi(a) = prod (z - phi(a_i)), has F[2] = E[2] as Galois modules via
   psi : (a_i, 0) -> (phi(a_i), 0), and the HLP (Prop. 4) gluing of E and F along psi is, up to a square,
       C : y^2 = h(x) = c2 * Fh(c1 x^2 - 1, c2 x^2),   Fh(X, Z) = X^3 + pE X Z^2 + qE Z^3
   (with w_i = c1 - c2 a_i: a1/a2 = -1/(c2^2 f(c1/c2)), b1/b2 = c2 in HLP's notation, and HLP's sextic
   becomes Delta_f^4 c2^3 W^2 h(x), W = w1 w2 w3; checked numerically against Prop. 4 in the lab note).
   The certificate does not rely on that derivation: it checks the two quotient maps directly.
       phi1 : (x, y) -> (x^2, y)       maps C onto Ea : Y^2 = Ha(U) = c2 Fh(c1 U - 1, c2 U),
       phi2 : (x, y) -> (1/x^2, y/x^3) maps C onto Eb : Y^2 = Hb(U) = c2 Fh(c1 - U, c2) = c2^4 f((c1 - U)/c2),
   so Eb = E (z = (c1 - U)/c2, y = Y/c2^2), and Ea is checked to be isomorphic to F.

   Argument.
   (1) Splitting.  At every rational point R of EM where h_R is squarefree of degree 6, phi1^* + phi2^* is an isogeny
       Ea x Eb -> J over Q, with Eb = E and Ea = F_t(R) (the identities below hold in Q(EM); where C_R is smooth, Ea_R
       is smooth, so t(R) is not a cusp of X_1(5) and the isomorphism specializes).
   (2) The marked group.  phi1, phi2 have prime degree 2 and C is not of genus 1, so they do not factor through
       isogenies and phi1^*, phi2^* are injective: phi2^*(P_E) has order 7 and phi1^*((0,0) of F_t) has order 5, so
       Z/35 embeds in J_R(Q).  Mumford arithmetic at the members below confirms it.
   (3) Exactness.  For R = R0 (mod 3), R0 = (91, -820) (t = 32/13), the curve C_R has good reduction at 3, with
       #J(F_3) = 35 (necessarily: 7 | #E(F_3) <= 7 and 5 | #F(F_3) <= 7).  J_R(Q)_tors injects into J(F_3), so it is
       exactly Z/35.
   (4) Infinitely many members.  P = R0 has order 3 in EM(F_3), so R_k = R0 + 3k P = (1 + 3k) P reduces to R0 mod 3;
       x(R_k), y(R_k) are 3-integral and congruent to x(R0), y(R0), and the coefficients c1, c2 = A(x) + B(x) y have
       denominators that are units at x(R0) mod 3, so C_{R_k} has the reduction of (3).  The R_k are distinct, and the
       G2-invariants of C_R are not constant on EM, so R -> M_2 has finite fibres and the R_k give infinitely many
       isomorphism classes.
   (5) Sanity: TorsionSubgroup at R0 (t = 32/13) and at (-9, 0) (t = -26) is [35].  The census curve 35.a is the
       member at (-25, 8) (t = -2/13), checked below by G2-invariants; it is not = R0 mod 3, but the bound in (3)
       holds for every member with good reduction at 3.
   Run (from the repository root):
       magma -b pipeline/magma/family_35_split_certificate.m > data/knowledge/sources/family_35_split_certificate.log */
SetColumns(0);
SetQuitOnError(true);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

function Kubert(b, c) return EllipticCurve([1 - c, -b, -b, 0, 0]); end function;
// y^2 + a1 xy + a3 y = x^3 + a2 x^2 + a4 x + a6  ->  Y^2 = z^3 + A z + B, z = x + sh, Y = y + (a1 x + a3)/2
function DepCubic(E)
    b2, b4, b6 := Explode(bInvariants(E)[1..3]);
    A2 := b2/4; A1 := b4/2; A0 := b6/4;
    return A1 - A2^2/3, A0 - A1*A2/3 + 2*A2^3/27, A2/3;
end function;
function MonEC(G)                                 // Y^2 = G(U), G a cubic, as a Weierstrass curve
    cs := [ Coefficient(G, i) : i in [0..3] ];
    return EllipticCurve([0, cs[3], 0, cs[2]*cs[4], cs[1]*cs[4]^2]), cs[4];
end function;
function OrderUpTo(P, n)
    for k in [1..n] do if k*P eq Parent(P)!0 then return k; end if; end for; return -1;
end function;
// The largest finite abelian group that embeds in every group with the given invariants.
function StructureBound(invlist)
    ls := &meet[ Seqset(&cat[ PrimeDivisors(n) : n in inv | n gt 1 ]) : inv in invlist ];
    out := [Integers()| ];
    for l in ls do
        es := [ Reverse(Sort([ Valuation(n, l) : n in inv | n mod l eq 0 ])) : inv in invlist ];
        r := Min([ #e : e in es ]);
        out cat:= [ l^Min([ e[i] : e in es ]) : i in [1..r] ];
    end for;
    return Invariants(AbelianGroup(out));
end function;

// ---------------------------------------------------------------- (0) the fixed curve E
E := Kubert(-2, 2);
assert Conductor(E) eq 26 and Invariants(TorsionSubgroup(E)) eq [7];
PE := E![0, 0]; assert Order(PE) eq 7;
pE, qE, shE := DepCubic(E);
aE := aInvariants(E);
zP := shE; wP := aE[3]/2;                          // P_E in the model y^2 = f(z)
fE := X^3 + pE*X + qE;
assert Evaluate(fE, zP) eq wP^2;
printf "(0) E = %o (conductor %o, torsion Z/%o, P_E = (0,0)); f(z) = %o\n", aE, 26, 7, fE;

// ---------------------------------------------------------------- (1) the base curve X' and the family over its function field
EM := EllipticCurve([0, 1, 0, -901, -7461]);
aM := aInvariants(EM); assert aM[1] eq 0 and aM[3] eq 0;
Fx<xx> := FunctionField(Q); Py<YY> := PolynomialRing(Fx);
K3<yy> := FunctionField(YY^2 - (xx^3 + aM[2]*xx^2 + aM[4]*xx + aM[5]));
tK := (-31*xx^4 - 4340*xx^3 - 60202*xx^2 + 2060012*xx + 6289745)/(13*xx^6 - 624*xx^5 + 8866*xx^4 + 66404*xx^3 - 2350907*xx^2 + 5272228*xx + 190297588)*yy + (xx^6 + 1197*xx^5 + 56030*xx^4 - 764462*xx^3 - 41528219*xx^2 + 85050081*xx + 2965565276)/(26*xx^6 - 1248*xx^5 + 17732*xx^4 + 132808*xx^3 - 4701814*xx^2 + 10544456*xx + 380595176);
c1K := (-987*xx^10 - 101194*xx^9 - 7839349*xx^8 + 528230572*xx^7 + 6322541939*xx^6 - 261586607530*xx^5 - 2041240435051*xx^4 + 18539162411544*xx^3 + 338538174658872*xx^2 + 1170195981745072*xx - 3861977298803824)/(1352*xx^12 - 129792*xx^11 + 4959136*xx^10 - 74706112*xx^9 - 523118544*xx^8 + 33987884736*xx^7 - 311270359712*xx^6 - 3649787998272*xx^5 + 76810516874312*xx^4 + 3872053164096*xx^3 - 6935579802725184*xx^2 + 16052676348577024*xx + 289705375988941952)*yy + (89*xx^12 + 10820*xx^11 + 1126304*xx^10 - 25895524*xx^9 - 3376807042*xx^8 + 21849685436*xx^7 + 2223285044532*xx^6 + 1962301668724*xx^5 - 413673972824063*xx^4 - 3528715576685664*xx^3 + 10716422264455604*xx^2 + 278530313321107280*xx + 963536501917589408)/(1352*xx^12 - 129792*xx^11 + 4959136*xx^10 - 74706112*xx^9 - 523118544*xx^8 + 33987884736*xx^7 - 311270359712*xx^6 - 3649787998272*xx^5 + 76810516874312*xx^4 + 3872053164096*xx^3 - 6935579802725184*xx^2 + 16052676348577024*xx + 289705375988941952);
c2K := (103*xx^10 - 22774*xx^9 - 2645151*xx^8 + 89060388*xx^7 + 2092194729*xx^6 - 25648973518*xx^5 - 805123071585*xx^4 + 289273213416*xx^3 + 70436929259808*xx^2 + 389573936613944*xx + 934967479083344)/(338*xx^12 - 32448*xx^11 + 1239784*xx^10 - 18676528*xx^9 - 130779636*xx^8 + 8496971184*xx^7 - 77817589928*xx^6 - 912446999568*xx^5 + 19202629218578*xx^4 + 968013291024*xx^3 - 1733894950681296*xx^2 + 4013169087144256*xx + 72426343997235488)*yy + (-11*xx^12 - 160*xx^11 + 406474*xx^10 - 90228*xx^9 - 738480602*xx^8 - 4503571508*xx^7 + 390909212064*xx^6 + 5086847406140*xx^5 - 64588011093291*xx^4 - 1140525529248364*xx^3 - 841095407952210*xx^2 + 47068461926634872*xx + 183933515112366328)/(338*xx^12 - 32448*xx^11 + 1239784*xx^10 - 18676528*xx^9 - 130779636*xx^8 + 8496971184*xx^7 - 77817589928*xx^6 - 912446999568*xx^5 + 19202629218578*xx^4 + 968013291024*xx^3 - 1733894950681296*xx^2 + 4013169087144256*xx + 72426343997235488);
printf "(1) base X' = EM : %o, conductor %o\n", aM, Conductor(EM);
// (1a) charpoly(c1 a + c2 (a^2 + 2pE/3)) = z^3 + a1(t) z + a0(t), the depressed 2-division cubic of the Kubert curve F_t
R2<u1, u2> := PolynomialRing(Q, 2);
Ma := Matrix(Q, 3, 3, [0, 1, 0, 0, 0, 1, -qE, -pE, 0]);
MM := u1*ChangeRing(Ma, R2) + u2*(ChangeRing(Ma*Ma, R2) + (2*pE/3)*IdentityMatrix(R2, 3));
cp := CharacteristicPolynomial(MM);
assert Coefficient(cp, 2) eq 0;
s2 := Coefficient(cp, 1); NN := -Coefficient(cp, 0);
Ft<t> := FunctionField(Q);
bF := t; cF := t;                                // F_t = E_5(t), Kubert, (0,0) of order 5
FKub := Kubert(bF, cF);
a1t, a0t, sht := DepCubic(FKub);
assert Evaluate(s2, [c1K, c2K]) eq Evaluate(a1t, tK) and Evaluate(NN, [c1K, c2K]) eq -Evaluate(a0t, tK);
print "(1a) over Q(X'): charpoly(c1 a + c2 (a^2 + 2pE/3)) = z^3 + a1(t) z + a0(t), the depressed 2-division cubic of the";
print "     Kubert curve F_t of the family below, at t = t(R): F_R[2] = E[2] as Galois modules via (a_i,0) -> (phi(a_i),0)";
// (1b) the curve and its two quotients
PK<x> := PolynomialRing(K3);
FhK := func< XX, ZZ | XX^3 + pE*XX*ZZ^2 + qE*ZZ^3 >;
h := c2K*FhK(c1K*x^2 - 1, c2K*x^2);
assert Degree(h) eq 6 and &and[ Coefficient(h, i) eq 0 : i in [1, 3, 5] ];
Ha := PK![ Coefficient(h, 2*i) : i in [0..3] ];                  // h = Ha(x^2)
Hb := PK![ Coefficient(h, 6 - 2*i) : i in [0..3] ];              // x^6 h(1/x) = Hb(x^2)
assert Hb eq c2K*FhK(c1K - x, c2K);
assert Hb eq c2K^4*(((c1K - x)/c2K)^3 + pE*(c1K - x)/c2K + qE);
Ea, la := MonEC(Ha);
FtK := Kubert(Evaluate(bF, tK), Evaluate(cF, tK));
c4a, c6a := Explode(cInvariants(Ea)); c4f, c6f := Explode(cInvariants(FtK));
assert c4a^3*c6f^2 eq c4f^3*c6a^2;                                 // equal j-invariants
lam := (c6a/c6f)/(c4a/c4f);
assert c4a eq lam^2*c4f and c6a eq lam^3*c6f and IsSquare(lam);    // c4a = u^4 c4f, c6a = u^6 c6f with u^2 = lam
print "(1b) C_R : y^2 = h(x) = c2 Fh(c1 x^2 - 1, c2 x^2); phi2 maps C onto Eb : Y^2 = c2^4 f((c1 - U)/c2) = E,";
print "     phi1 maps C onto Ea : Y^2 = Ha(U), and Ea is isomorphic to F_t over Q(X') (c4(Ea) = lam^2 c4(F_t), c6(Ea) = lam^3 c6(F_t)";
print "     with lam a square in Q(X')): J(C_R) is isogenous over Q to E x F_t(R)";

// ---------------------------------------------------------------- reduction of the family at a point of X'
function RedAt(c, p, R)                           // c = A + B yy in K3, R = (x0, y0) p-integral, reduced mod p
    A, B := Explode(Eltseq(c));
    val := GF(p)!0; x0 := GF(p)!R[1]; y0 := GF(p)!R[2];
    for pr in [ <A, GF(p)!1>, <B, y0> ] do
        N, D := IntNumDen(pr[1]);
        assert Evaluate(PolynomialRing(GF(p))!D, x0) ne 0;
        val +:= pr[2]*Evaluate(PolynomialRing(GF(p))!N, x0)/Evaluate(PolynomialRing(GF(p))!D, x0);
    end for;
    return val;
end function;
function SpecAt(R)                                 // h, t, c1, c2 at a rational point R of X'
    ev := func< c | Evaluate(Eltseq(c)[1], R[1]) + Evaluate(Eltseq(c)[2], R[1])*R[2] >;
    c1 := ev(c1K); c2 := ev(c2K); tv := ev(tK);
    hR := c2*((c1*X^2 - 1)^3 + pE*(c1*X^2 - 1)*(c2*X^2)^2 + qE*(c2*X^2)^3);
    return hR, tv, c1, c2;
end function;

// ---------------------------------------------------------------- (3) exactness on the residue class R = R0 mod primes
P := EM![91, -820]; R0 := EM![91, -820];
assert Order(P) eq 0;
primes := [3];
invs := []; nJ := [];
for p in primes do
    assert Valuation(Discriminant(EM), p) eq 0 and Valuation(R0[1], p) ge 0 and Valuation(R0[2], p) ge 0;
    assert &and[ Valuation(c, p) ge 0 : c in [pE, qE] ];
    c1b := RedAt(c1K, p, R0); c2b := RedAt(c2K, p, R0);
    Pp<xp> := PolynomialRing(GF(p));
    hb := c2b*((c1b*xp^2 - 1)^3 + (GF(p)!pE)*(c1b*xp^2 - 1)*(c2b*xp^2)^2 + (GF(p)!qE)*(c2b*xp^2)^3);
    assert Degree(hb) eq 6 and IsSquarefree(hb);
    n, chi := FrobData(hb);
    A := AbelianGroup(Jacobian(HyperellipticCurve(hb)));
    assert #A eq n;
    Append(~invs, Invariants(A)); Append(~nJ, n);
    printf "(3) p = %o, R = R0 mod %o: reduction %o (squarefree, degree 6); #J(F_p) = %o, J(F_p) = %o\n", p, p, hb, n, Invariants(A);
end for;
g := GCD(nJ);
printf "(3) gcd of #J(F_p) = %o = #G\n", g;
assert g eq 35;
print "    J(Q)_tors injects into J(F_3) (good reduction at the odd prime 3), and contains Z/35 by (2): it is exactly Z/35";

// ---------------------------------------------------------------- (4) infinitely many members
ords := [ Order(ChangeRing(EM, GF(p))!Eltseq(P)) : p in primes ];
mstep := LCM(ords);
printf "(4) P = %o has orders %o mod %o; for R_k = R0 + %o k P the point R_k reduces to R0 mod every p in %o,\n", P, ords, primes, mstep, primes;
print "    so x(R_k), y(R_k) are p-integral and congruent to x(R0), y(R0), and C_{R_k} has the reductions of (3)";
R1 := R0 + mstep*P;
h1 := SpecAt(R1);
for i in [1..#primes] do
    p := primes[i];
    c1b := RedAt(c1K, p, R0); c2b := RedAt(c2K, p, R0);
    Pp<xp> := PolynomialRing(GF(p));
    hb := c2b*((c1b*xp^2 - 1)^3 + (GF(p)!pE)*(c1b*xp^2 - 1)*(c2b*xp^2)^2 + (GF(p)!qE)*(c2b*xp^2)^3);
    assert &and[ Valuation(c, p) ge 0 : c in Coefficients(h1) ] and Pp!h1 eq hb;
end for;
printf "(4) sanity: h at R_1 = R0 + %o P has coefficients of up to %o digits and the same reductions\n", mstep, Max([ #Sprint(c) : c in Coefficients(h1) ]);
h0 := SpecAt(R0);
G20 := G2Invariants(HyperellipticCurve(h0));
k2 := 0;
for k in [2..6] do
    hk := SpecAt(k*P);
    if Degree(hk) eq 6 and Discriminant(hk) ne 0 and G2Invariants(HyperellipticCurve(hk)) ne G20 then k2 := k; break; end if;
end for;
assert k2 gt 0;
printf "(4) C at %o P has G2-invariants different from C at R0: R -> M_2 is nonconstant on X', so the R_k give infinitely many isomorphism classes\n", k2;

// ---------------------------------------------------------------- (2), (5) the marked group and TorsionSubgroup at members
function MemberCheck(R)
    hR, tv, c1, c2 := SpecAt(R);
    C := HyperellipticCurve(hR); J := Jacobian(C); PJ := Parent((J!0)[1]); xj := PJ.1;
    // P_E on Eb : Y^2 = c2^4 f((c1 - U)/c2): U = c1 - c2 zP, Y = c2^2 wP; pullback by phi2
    Ub := c1 - c2*zP; Yb := c2^2*wP;
    D0 := elt< J | xj^2 - 1/Ub, (Yb/Ub)*xj >;
    // the point (0,0) of F_t on Ea : Y^2 = Ha(U) (Ha = h(sqrt U)); pullback by phi1
    HaR := RX![ Coefficient(hR, 2*i) : i in [0..3] ];
    EaR, laR := MonEC(HaR);
    FR := Kubert(Evaluate(bF, tv), Evaluate(cF, tv));
    ok, iso := IsIsomorphic(FR, EaR); assert ok;
    Q1 := iso(FR![0, 0]);
    D1 := elt< J | xj^2 - Q1[1]/laR, PJ!(Q1[2]/laR) >;
    o0 := OrderUpTo(D0, 100); o1 := OrderUpTo(D1, 100); o01 := OrderUpTo(D0 + D1, 200);
    return hR, tv, o0, o1, o01;
end function;
for R in [EM![91, -820], EM![-9, 0]] do
    hR, tv, o0, o1, o01 := MemberCheck(R);
    inv := TorsionInvariants(hR);
    printf "(2),(5) R = %o (t = %o): phi2^*(P_E) has order %o, phi1^*((0,0) of F_t) has order %o, their sum has order %o; TorsionSubgroup %o\n", R, tv, o0, o1, o01, inv;
    assert o0 eq 7 and o1 eq 5 and o01 eq 35;
    assert inv eq [35];
end for;
// the census curve 35.a is the member at (-25, 8)
hc := HyperellipticCurve(RX![160, 80, 166, 182, 116, 30, 10], RX![0, 1, 1]);
hm, tm := SpecAt(EM![-25, 8]);
assert tm eq -2/13 and G2Invariants(HyperellipticCurve(hm)) eq G2Invariants(hc);
print "(5) the member at (-25, 8) (t = -2/13) has the G2-invariants of the census curve 35.a";
printf "time %o s\n", Cputime(T0);
print "FAMILY_35_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E x F_R) and J(Q)_tors exactly [35]";
quit;
