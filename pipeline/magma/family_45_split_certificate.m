/* family_45_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E x F_R) with J(Q)_tors exactly [45].
   Found 2026-10-19 (torsion_jac notes/hlpEC_2026_10_19.md).

   Source.  HLP 2000, section 3.6 (the [45] entry of Table 1): the (2,2)-gluing (HLP Prop. 4) of the fixed curve
       E = E_9(-5) : y^2 + (1 - c)xy - by = x^3 - bx^2, b = -4650, c = -150   (conductor 15810 = 15810u1, torsion Z/9)
   with Kubert X_1(5) curves F_t = E_5(t) along a Galois-module isomorphism F_t[2] -> E[2] (no rational 2-torsion on
   either side).  The pairs (t, psi) form a genus-1 twist X' of X_1(2,10), 3-isogenous to HLP's
   X'' : y^2 = Delta_9(-5) t(t^2 - 11t - 1) (conductor 13838400, rank 2).  Rebuilt explicitly, X' is
       EM : y^2 = x^3 - 1153200x - 327701000   (conductor 13838400, rank 2, MW = Z/2 x Z^2, e.g. (-310,0), (1209,-6727)),
   with the functions t, c1, c2 of Q(EM) below (found by factoring the sextic resolvent over Q(EM)); HLP's point
   t = 93/10 is the image of O.

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
   (2) The marked group.  phi1^*, phi2^* are injective (prime degree 2, C not of genus 1): phi2^*(P_E) has order 9 and
       phi1^*((0,0) of F_t) has order 5, so Z/45 embeds in J_R(Q).  Mumford arithmetic at the members confirms it.
   (3) Exactness.  For R = R0 (mod 7), R0 = (1209, -6727) (t = 27/640), the curve C_R has good reduction at 7, with
       #J(F_7) = 45 (#E(F_7) = 9 and #F(F_7) = 5: E[2](F_7) = F[2](F_7) = 0).  So J_R(Q)_tors is exactly Z/45.
   (4) Infinitely many members.  P = R0 has order 2 in EM(F_7), so R_k = R0 + 2k P reduces to R0 mod 7, and C_{R_k}
       has the reduction of (3) (x, y 7-integral and congruent; the denominators of c1, c2 = A(x) + B(x) y are units at
       x(R0) mod 7).  The G2-invariants are not constant on EM: infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup at R0 (t = 27/640) and at (-310, 0) (t = -10/93) is [45].  The census curve 45.a is the
       member at (-310, 0) (= R0 mod 7, so in the class) and at O (HLP's t = 93/10); checked below by G2-invariants.
       In fact every member with good reduction at 7 has #J(F_7) = 45.
   Run (from the repository root):
       magma -b pipeline/magma/family_45_split_certificate.m > data/knowledge/sources/family_45_split_certificate.log */
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
E := Kubert(-4650, -150);
assert Conductor(E) eq 15810 and Invariants(TorsionSubgroup(E)) eq [9];
PE := E![0, 0]; assert Order(PE) eq 9;
pE, qE, shE := DepCubic(E);
aE := aInvariants(E);
zP := shE; wP := aE[3]/2;                          // P_E in the model y^2 = f(z)
fE := X^3 + pE*X + qE;
assert Evaluate(fE, zP) eq wP^2;
printf "(0) E = %o (conductor %o, torsion Z/%o, P_E = (0,0)); f(z) = %o\n", aE, 15810, 9, fE;

// ---------------------------------------------------------------- (1) the base curve X' and the family over its function field
EM := EllipticCurve([0, 0, 0, -1153200, -327701000]);
aM := aInvariants(EM); assert aM[1] eq 0 and aM[3] eq 0;
Fx<xx> := FunctionField(Q); Py<YY> := PolynomialRing(Fx);
K3<yy> := FunctionField(YY^2 - (xx^3 + aM[2]*xx^2 + aM[4]*xx + aM[5]));
tK := (-3813*xx^4 - 18912480*xx^3 - 8794303200*xx^2 + 11813680632000*xx + 1126843383360000)/(5*xx^6 + 77190*xx^5 + 512410005*xx^4 + 1789867112800*xx^3 + 3336134648568000*xx^2 + 2878067937041280000*xx + 900195693641244800000)*yy + (93*xx^6 + 934743*xx^5 + 1716970650*xx^4 - 1228414010400*xx^3 - 2214688321932000*xx^2 - 286080799448640000*xx + 14117698554403200000)/(10*xx^6 + 154380*xx^5 + 1024820010*xx^4 + 3579734225600*xx^3 + 6672269297136000*xx^2 + 5756135874082560000*xx + 1800391387282489600000);
c1K := (-39060961*xx^10 - 1770755142096*xx^9 - 17459949717035066*xx^8 - 64162147420445094916*xx^7 - 78685506057381623783081*xx^6 + 54554856157467123772765820*xx^5 + 200310874690051110239891259100*xx^4 + 110808346113270366123548054576000*xx^3 - 68106036278183609279709213350240000*xx^2 - 86875128950999146179764172651084800000*xx - 23480597845210298916535982023363136000000)/(4590000*xx^12 + 141720840000*xx^11 + 2034727933140000*xx^10 + 17810031285701640000*xx^9 + 105064190185593762990000*xx^8 + 436620037397724346204800000*xx^7 + 1299131551558006303031712000000*xx^6 + 2759682274201873454400207360000000*xx^5 + 4104387409676984215035925824000000000*xx^4 + 4117359560214694068977961056256000000000*xx^3 + 2623574808195840561203796547768320000000000*xx^2 + 951350706068447357366082461658316800000000000*xx + 148780679865704406165862340531564544000000000000)*yy + (-454511*xx^12 + 4640162584*xx^11 + 162528071610514*xx^10 + 1063314458886081364*xx^9 + 2442894021928427962449*xx^8 + 642096278174929598291100*xx^7 - 5172753475460646854404398300*xx^6 - 6602421611547691674457632540000*xx^5 - 260565793124607685783840668000000*xx^4 + 4591969928189217943825747978396000000*xx^3 + 3325521004199908127879374889513248000000*xx^2 + 849628791145252386048447281056892800000000*xx + 52236581590190020104850940835982873600000000)/(4590000*xx^12 + 141720840000*xx^11 + 2034727933140000*xx^10 + 17810031285701640000*xx^9 + 105064190185593762990000*xx^8 + 436620037397724346204800000*xx^7 + 1299131551558006303031712000000*xx^6 + 2759682274201873454400207360000000*xx^5 + 4104387409676984215035925824000000000*xx^4 + 4117359560214694068977961056256000000000*xx^3 + 2623574808195840561203796547768320000000000*xx^2 + 951350706068447357366082461658316800000000000*xx + 148780679865704406165862340531564544000000000000);
c2K := (-961*xx^10 - 43140096*xx^9 - 423994405466*xx^8 - 1553367426405316*xx^7 - 1890165509333779481*xx^6 + 1362505637007191173820*xx^5 + 4900787202066719035299100*xx^4 + 2714006425292232057628976000*xx^3 - 1640665194041929247180006240000*xx^2 - 2104336991323944058058736204800000*xx - 569639327781429508107449001536000000)/(382500*xx^12 + 11810070000*xx^11 + 169560661095000*xx^10 + 1484169273808470000*xx^9 + 8755349182132813582500*xx^8 + 36385003116477028850400000*xx^7 + 108260962629833858585976000000*xx^6 + 229973522850156121200017280000000*xx^5 + 342032284139748684586327152000000000*xx^4 + 343113296684557839081496754688000000000*xx^3 + 218631234016320046766983045647360000000000*xx^2 + 79279225505703946447173538471526400000000000*xx + 12398389988808700513821861710963712000000000000)*yy + (-11*xx^12 + 113584*xx^11 + 3944663014*xx^10 + 25718558694964*xx^9 + 58701684239018049*xx^8 + 13929865741480331100*xx^7 - 127983489598386936918300*xx^6 - 162263292105492401928540000*xx^5 - 7204834268522111135868000000*xx^4 + 111512424857095827447000796000000*xx^3 + 81131521364439224504343324448000000*xx^2 + 20877069391348516102296441212800000000*xx + 1314022382955617099451241010713600000000)/(382500*xx^12 + 11810070000*xx^11 + 169560661095000*xx^10 + 1484169273808470000*xx^9 + 8755349182132813582500*xx^8 + 36385003116477028850400000*xx^7 + 108260962629833858585976000000*xx^6 + 229973522850156121200017280000000*xx^5 + 342032284139748684586327152000000000*xx^4 + 343113296684557839081496754688000000000*xx^3 + 218631234016320046766983045647360000000000*xx^2 + 79279225505703946447173538471526400000000000*xx + 12398389988808700513821861710963712000000000000);
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
P := EM![1209, -6727]; R0 := EM![1209, -6727];
assert Order(P) eq 0;
primes := [7];
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
assert g eq 45;
print "    J(Q)_tors injects into J(F_7) (good reduction at the odd prime 7), and contains Z/45 by (2): it is exactly Z/45";

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
for R in [EM![1209, -6727], EM![-310, 0]] do
    hR, tv, o0, o1, o01 := MemberCheck(R);
    inv := TorsionInvariants(hR);
    printf "(2),(5) R = %o (t = %o): phi2^*(P_E) has order %o, phi1^*((0,0) of F_t) has order %o, their sum has order %o; TorsionSubgroup %o\n", R, tv, o0, o1, o01, inv;
    assert o0 eq 9 and o1 eq 5 and o01 eq 45;
    assert inv eq [45];
end for;
// the census curve 45.a is the member at (-310, 0), which is = R0 mod 7
hc := HyperellipticCurve(RX![168300000000, 0, 49996210000, 0, 29240200, 0, 13981]);
hm, tm := SpecAt(EM![-310, 0]);
assert G2Invariants(HyperellipticCurve(hm)) eq G2Invariants(hc);
assert GF(7)!(-310) eq GF(7)!R0[1] and GF(7)!0 eq GF(7)!R0[2];
print "(5) the member at (-310, 0) (= R0 mod 7) has the G2-invariants of the census curve 45.a";
printf "time %o s\n", Cputime(T0);
print "FAMILY_45_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E x F_R) and J(Q)_tors exactly [45]";
quit;
