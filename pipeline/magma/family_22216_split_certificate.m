/* family_22216_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E1 x E2) with J(Q)_tors exactly [2,2,2,16].
   Found 2026-10-16 (torsion_jac notes gain2_loci_2026_10_15.md, gain2_22216_fibres_2026_10_16.md,
   gain2_22216_certificate_2026_10_16.md); the member C_{t0} below is the census curve 2.2.2.16.a.

   Source.  The (2,2)-gluing of E1 = E_{2,4}(t) and E2 = E_{2,8}(u) along alpha = identity on the labelled roots of
   2-torsion, with full gain g = 4 (torsion_jac notes gain2_loci_2026_10_15.md, gain2_22216_fibres_2026_10_16.md):
   the locus a012_W1 is  X(u) t = sq,  V(u)(t^2 - 1/16) = sq,  2u(u+1/2)(u^2-1/8)(u^2+u/2+1/8)(t + 1/4) = sq,
   with X(u) = (u^2-1/8)(u^2+u+1/8), V(u) = u(u+1/2)(u+1/4).  Fix u = -20/31.  Then the conditions read
       t = A a^2,   t + 1/4 = B c^2,   t - 1/4 = C e^2,     A = -1788961, B = 11195, C = -138818,
   a genus-1 curve D_u (a 2-covering of E_d : y^2 = x^3 - d^2 x, d = 554577910 = sqf(ABC), via x = 4 d t).
   E1 = E_{2,4}(t) : y^2 + xy - b y = x^3 - b x^2, b = t^2 - 1/16, P1 = (0,0) of order 4, roots
   e = (t^2 - 1/16, -t/2 - 1/8, t/2 - 1/8) (special root first);  E2 = E_{2,8}(-20/31) (HLP Table 3 Tate form),
   P2 = (0,0) of order 8, roots e' (special root first).  HLP 2000 Prop. 4 glues them to y^2 = hHLP(x).

   The family (derived here).  With A(t), B(t) the HLP scalars, B(t)/A(t) = kappa t (t^2-1/16)^2 for a constant kappa,
   and the factors of hHLP factor further, so that hHLP = (A(t) t (t^2-1/16)^2)^2 hc with
       C_t : y^2 = hc_t(x) = -A(t) (k1 t - x^2) (k2 (t-1/4)^2 - x^2) (x^2 + k3 (t+1/4)^2),
       k1 = -34405369600/2257783039, k2 = 32083974400/2217373921, k3 = -12036503521/41990400,
   A(t) a quadratic polynomial in t (printed below); the coefficients of hc have degree <= 7 in t.
   (x,y) -> (x^2, y) maps C_t onto a curve Q(t)-isomorphic to E2, and (x,y) -> (1/x^2, y/x^3) onto one
   Q(t)-isomorphic to E1(t) (checked below with IsIsomorphic over Q(t)), so every smooth member has J split over Q,
   isogenous to E1(t) x E2.

   Argument.
   (1) The marked group.  Write pull1, pull2 for the pullbacks of points of E1, E2 to J (through the two maps and the
       Q(t)-isomorphisms; pull2(X,Y) = [x^2 - X, Y], pull1(X,Y) = [x^2 - 1/X, (Y/X) x] on the quotient models).
       (a) Over Q(t): G0 = <pull1(P1), pull1(T1), pull2(P2), pull2(T2)> (T1, T2 non-special 2-torsion points), and
           pull1(T) + pull2(alpha(T)) = 0 for the three T in E1[2]: the kernel on 2-torsion is graph(alpha).
       (b) Over Q(n), n the parameter of the conic B c^2 - C e^2 = 1/2 (t = B c^2 - 1/4): the 16-torsion class.
           Over F(n), F = Q(sqrt(2B), sqrt(-2C)), the half R1 of P1 on E1 has x(R1) = r1 r2 + r1 r3 + r2 r3 with
           r1 = 2 r2 r3, r2 = c sqrt(2B)/2, r3 = e sqrt(-2C)/2 (r_i^2 = -e_i), and a half R2 of the matching point of
           order 8 on E2 is defined over F.  D16 = pull1(R1) + pull2(R2) is computed over F(n) by Mumford
           arithmetic; its coefficients lie in Q(n) (Galois descent made explicit).  Over Q(n): 16 D16 = 0,
           8 D16 != 0 and 2 D16 lies in G0.
       (c) Over Q(a), t = A a^2: the class D2 from the non-special 2-torsion.  Over F'(a), F' = Q(i, sqrt(-A)), a half
           of T1 = (e_2, 0) is defined (r1 = i (t + 1/4), r3 = a sqrt(-A)); with a half of the matching 2-torsion
           point of E2 (defined over F'), D2 = pull1(.) + pull2(.) has coefficients in Q(a).  It turns out to be a
           class of order 2 (its v is 0: a quadratic factor of hc over Q(a)).
       All three fields lie in K = Q(D_u), so G = <G0, D16, D2> is a finite subgroup of J(K)[16].  At every point
       R of D_u(Q) over which C_t is smooth, J[16] is finite etale over the local ring, so specialization is
       injective on J(K)[16] and G embeds in J_{t(R)}(Q)_tors.  The specialization at the point over
       t0 = -4333968559/38784504316 (the census curve 2.2.2.16.a) has 128 elements, 16 of them killed by 2, exponent 16: G is
       [2,2,2,16] (the only abelian group of order 2^7 with 2-rank 4 and exponent 16).
   (2) Exactness on a residue class.  On E_d (d = 554577910) let P = (-47963453395683600/226291849, ...) and
       T3 = (d, 0); R0 = P + T3 has x(R0) = 4 d t0.  P has order 6 mod 11 and order 4 mod 13, so for every k the
       point R_k = R0 + 12 k P reduces to R0 mod 11 and mod 13.  R_k - R0 lies in 2 E_d(Q), so delta(R_k) = delta(R0)
       and t_k = x(R_k)/(4d) lies on D_u: the three locus conditions hold.  t_k is 11- and 13-integral and
       congruent to t0 mod 11 and 13 (11, 13 do not divide 2d and R0 is not O mod p).  The coefficients of hc are
       polynomials in t whose denominators are prime to 11 and 13, and hc_{t0} mod 11 and mod 13 are squarefree of
       degree 6 (ReductionClass of family_cert_lib.m), so every C_{t_k} has good reduction at 11 and 13 with the
       reduction of C_{t0}: J(F_11) = J(F_13) = [2,2,2,16], of order 128.  J_{t_k}(Q)_tors injects into J(F_11)
       (p odd, good reduction), so it has order <= 128, and it contains G: it is exactly [2,2,2,16].
   (3) Infinitely many.  P has infinite order, so the R_k are pairwise distinct and t_k takes infinitely many values
       (x takes each value at most twice).  The G2-invariants of C_{t0} and C_{t(3P+T3)} differ, so t -> M_2 is
       nonconstant on D_u, with finite fibres: infinitely many Q-isomorphism classes.
   Remark (not used above): E_d has rank exactly 1 (Cassels-Tate pairing on Sel^2) and P is saturated (Sage/eclib)
   (torsion_jac data/gain2c_2026_10_16/rank_fibre_554577910.log, saturation_sage.log), so the points of D_u(Q) lie over t((2k+1)P + T3), k >= 0, and t0 is
   the one of smallest height.
   (4) Sanity: TorsionSubgroup at t0 and at t(3P+T3) (numerator of 90 digits; skip with Tors2 := 0) is [2,2,2,16]; the
       quotients of C_{t0}, C_{t1} are Q-isomorphic to E1(t), E2; C_{t0} is Q-isomorphic to the census curve 2.2.2.16.a.
   No step relies on the gain formula of HLP Prop. 12 or on the gain2 Phase A tables: the marked group is built and
   checked by Mumford arithmetic over Q(t), Q(n), Q(a), and its structure is read off a specialization.
   Run:  magma -b pipeline/magma/family_22216_split_certificate.m > data/knowledge/sources/family_22216_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
if not assigned Tors2 then Tors2 := 1; elif Type(Tors2) eq MonStgElt then Tors2 := StringToInteger(Tors2); end if;
T0 := Cputime();

function Sqf(q)
    q := Rationals()!q; n := Numerator(q)*Denominator(q);
    return Sign(n)*&*([1] cat [p[1] : p in Factorization(Abs(n)) | p[2] mod 2 eq 1]);
end function;
function Ord2(g)      // order of a class killed by 64, else -1 (Order is not available over function fields)
    z := Parent(g)!0; h := g; k := 1;
    while k le 64 do if h eq z then return k; end if; h := 2*h; k := 2*k; end while;
    return -1;
end function;

// ---------------------------------------------------------------- (0) the fixed curve E2 and the fibre constants
uv := -20/31;
Xu := (uv^2 - 1/8)*(uv^2 + uv + 1/8); Vu := uv*(uv + 1/2)*(uv + 1/4);
Wu := uv*(uv + 1/2)*(uv^2 - 1/8)*(uv^2 + uv/2 + 1/8);
Acl := Sqf(Xu); Bcl := Sqf(2*Wu); Ccl := Sqf(Vu*Bcl);
assert [Acl, Bcl, Ccl] eq [-1788961, 11195, -138818];
dd := Abs(Sqf(Acl*Bcl*Ccl)); assert dd eq 554577910;
bu := (16*uv^3+16*uv^2+6*uv+1)/(8*uv^2-1)^2; cu := (16*uv^3+16*uv^2+6*uv+1)/(2*uv*(4*uv+1)*(8*uv^2-1));
E2 := EllipticCurve([Q| 1-cu, -bu, -bu, 0, 0]);
e2 := [Q| (uv^3/4 + 3*uv^2/16 + uv/32)/(uv^4 - uv^2/4 + 1/64),
          (-uv^4/4 - uv^3/4 - 3*uv^2/32 - uv/64)/(uv^4 + uv^3/2 - uv^2/16 - uv/16 - 1/128),
          (uv^3/4 + 3*uv^2/16 + uv/16 + 1/128)/(uv^4 - uv^2/8)];
P2 := E2![0,0];
assert Order(P2) eq 8 and (4*P2)[1] eq e2[1];
T2s := [E2![z, -((1-cu)*z - bu)/2] : z in e2];          // 2-torsion points (e'_1 special)
assert &and[2*T eq E2!0 : T in T2s];
assert Invariants(TorsionSubgroup(E2)) eq [2,8];
printf "E2 = E_{2,8}(%o) = %o, torsion [2,8], P2 = (0,0) of order 8\n", uv, aInvariants(E2);
printf "fibre u = %o: t = A a^2, t + 1/4 = B c^2, t - 1/4 = C e^2 with (A, B, C) = (%o, %o, %o); d = %o\n", uv, Acl, Bcl, Ccl, dd;

// ---------------------------------------------------------------- (1) the family over Q(t)
Ft<t> := FunctionField(Q); Pt<x> := PolynomialRing(Ft);
al := [t^2 - 1/16, -t/2 - 1/8, t/2 - 1/8]; be := [Ft!z : z in e2];
b1 := t^2 - 1/16;
E1 := EllipticCurve([Ft| 1, -b1, -b1, 0, 0]);
P1 := E1![0,0];
assert 4*P1 eq E1!0 and 2*P1 ne E1!0 and (2*P1)[1] eq al[1];
T1s := [E1![z, -(z - b1)/2] : z in al];
assert &and[2*T eq E1!0 and T ne E1!0 : T in T1s];
// HLP Prop. 4 with alpha = identity on the labelled roots
a1 := (al[3]-al[2])^2/(be[3]-be[2]) + (al[2]-al[1])^2/(be[2]-be[1]) + (al[1]-al[3])^2/(be[1]-be[3]);
bb1 := (be[3]-be[2])^2/(al[3]-al[2]) + (be[2]-be[1])^2/(al[2]-al[1]) + (be[1]-be[3])^2/(al[1]-al[3]);
a2 := al[1]*(be[3]-be[2]) + al[2]*(be[1]-be[3]) + al[3]*(be[2]-be[1]);
bb2 := be[1]*(al[3]-al[2]) + be[2]*(al[1]-al[3]) + be[3]*(al[2]-al[1]);
Df := ((al[1]-al[2])*(al[1]-al[3])*(al[2]-al[3]))^2;
Dg := ((be[1]-be[2])*(be[1]-be[3])*(be[2]-be[3]))^2;
AA := Dg*a1/a2; BB := Df*bb1/bb2;
hHLP := -(AA*(al[2]-al[1])*(al[1]-al[3])*x^2 + BB*(be[2]-be[1])*(be[1]-be[3]))
        *(AA*(al[3]-al[2])*(al[2]-al[1])*x^2 + BB*(be[3]-be[2])*(be[2]-be[1]))
        *(AA*(al[1]-al[3])*(al[3]-al[2])*x^2 + BB*(be[1]-be[3])*(be[3]-be[2]));
kap := BB/(AA*t*(t^2 - 1/16)^2);
assert kap in Q;
kq := [Q| kap*(be[2]-be[1])*(be[1]-be[3]), kap*(be[3]-be[2])*(be[2]-be[1]), kap*(be[1]-be[3])*(be[3]-be[2])];
assert kq eq [-34405369600/2257783039, 32083974400/2217373921, -12036503521/41990400];
assert Denominator(AA) eq 1 and Degree(Numerator(AA)) eq 2;
hc := -AA*(kq[1]*t - x^2)*(kq[2]*(t - 1/4)^2 - x^2)*(x^2 + kq[3]*(t + 1/4)^2);
assert hHLP eq (AA*t*(t^2 - 1/16)^2)^2*hc;
assert &and[Denominator(c) eq 1 and Degree(Numerator(c)) le 7 : c in Coefficients(hc)];
printf "A(t) = %o\n", AA;
printf "C_t : y^2 = hc(x) = %o\n", hc;
print "hHLP = (A(t) t (t^2-1/16)^2)^2 hc: C_t is Q(t)-isomorphic to the HLP gluing of E1(t) and E2 along alpha = id";
// the two quotients
PXX<XX> := PolynomialRing(Ft);
Hq := PXX![Coefficient(hc, 2*i) : i in [0..3]];          // hc(x) = Hq(x^2)
Hr := PXX![Coefficient(hc, 6 - 2*i) : i in [0..3]];      // x^6 hc(1/x) = Hr(x^2)
function MonEC(Hp) c := Coefficients(Hp); return EllipticCurve([0, c[3], 0, c[2]*c[4], c[1]*c[4]^2]); end function;
EA := MonEC(Hq); EB := MonEC(Hr);      // Y^2 = Hq(X), Y^2 = Hr(X) in monic form (x = c4 X, y = c4 Y)
E2t := ChangeRing(E2, Ft);
okA, mA := IsIsomorphic(E2t, EA);
okB, mB := IsIsomorphic(E1, EB);
assert okA and okB;
print "(x,y) -> (x^2, y) maps C_t onto Y^2 = Hq(X) ~ E2 and (x,y) -> (1/x^2, y/x^3) onto Y^2 = Hr(X) ~ E1(t) (Q(t)-isomorphisms): J split over Q";
// apply an isomorphism given by IsomorphismData, transported to another field by sub (convention fixed on P1)
dA := IsomorphismData(mA); dB := IsomorphismData(mB);
im := mB(P1); r,s,tt,uu := Explode(dB);
conv := (uu^2*P1[1] + r eq im[1] and uu^3*P1[2] + s*uu^2*P1[1] + tt eq im[2]) select 2 else 1;
if conv eq 1 then assert (P1[1]-r)/uu^2 eq im[1] and (P1[2] - s*(P1[1]-r) - tt)/uu^3 eq im[2]; end if;
function ApplyIso(dat, pt, sub)
    r,s,tt,uu := Explode([sub(z) : z in dat]);
    if conv eq 1 then return [(pt[1]-r)/uu^2, (pt[2] - s*(pt[1]-r) - tt)/uu^3]; end if;
    return [uu^2*pt[1] + r, uu^3*pt[2] + s*uu^2*pt[1] + tt];
end function;
// the family and the two pullbacks over a field FF containing tF (the image of t)
function Build(FF, tF)
    sub := func< f | Evaluate(Numerator(f), tF)/Evaluate(Denominator(f), tF) >;
    PF<xF> := PolynomialRing(FF);
    hF := PF![sub(c) : c in Coefficients(hc)];
    J := Jacobian(HyperellipticCurve(hF));
    cA := sub(Coefficients(Hq)[4]); cB := sub(Coefficients(Hr)[4]);
    pull1 := function(pt)        // E1 -> J through (1/x^2, y/x^3)
        if pt eq Parent(pt)!0 then return J!0; end if;
        im := ApplyIso(dB, [pt[1], pt[2]], sub); Xv := im[1]/cB; Yv := im[2]/cB;
        return elt< J | xF^2 - 1/Xv, (Yv/Xv)*xF >;
    end function;
    pull2 := function(pt)        // E2 -> J through (x^2, y)
        if pt eq Parent(pt)!0 then return J!0; end if;
        im := ApplyIso(dA, [FF!pt[1], FF!pt[2]], sub); Xv := im[1]/cA; Yv := im[2]/cA;
        return elt< J | xF^2 - Xv, PF!Yv >;
    end function;
    return J, pull1, pull2;
end function;

// ---------------------------------------------------------------- (1a) G0 over Q(t)
Jt, pl1, pl2 := Build(Ft, t);
gA := pl1(P1); gT := pl1(T1s[2]); gB := pl2(P2); gU := pl2(T2s[2]);
assert [Ord2(gA), Ord2(gT), Ord2(gB), Ord2(gU)] eq [4, 2, 8, 2];
for i in [1..3] do assert pl1(T1s[i]) + pl2(T2s[i]) eq Jt!0; end for;
assert &and[pl1(T1s[i]) + pl2(T2s[j]) ne Jt!0 : i, j in [1..3] | i ne j];
print "(1a) over Q(t): pull1(P1), pull1(T1), pull2(P2), pull2(T2) have orders 4, 2, 8, 2; pull1(T) + pull2(alpha T) = 0 exactly for alpha = id";
G0 := [i*gA + j*gT + k*gB + l*gU : i in [0..3], j in [0..1], k in [0..7], l in [0..1]];
function InG0(D, sub)        // index of D among the specializations of G0, or 0
    for i in [1..#G0] do
        g := G0[i];
        gs := elt< Parent(D) | Parent(D[1])![sub(c) : c in Coefficients(g[1])], Parent(D[1])![sub(c) : c in Coefficients(g[2])], g[3] >;
        if gs eq D then return i; end if;
    end for;
    return 0;
end function;

function InQ(z)              // z in FF(m) (FF a number field) has rational coefficients
    nu := Numerator(z); de := Denominator(z); lc := LeadingCoefficient(de);
    return &and[ IsCoercible(Q, c) : c in Coefficients(nu/lc) cat Coefficients(de/lc) ];
end function;
function Descend(D, Kq)      // a class over FF(m) with coefficients in Q(m) -> the class over Q(m) = Kq
    toQ := func< z | Evaluate(PolynomialRing(Q)![Q!c : c in Coefficients(Numerator(z)/LeadingCoefficient(Denominator(z)))], Kq.1)
                    / Evaluate(PolynomialRing(Q)![Q!c : c in Coefficients(Denominator(z)/LeadingCoefficient(Denominator(z)))], Kq.1) >;
    return [PolynomialRing(Kq)![toQ(c) : c in Coefficients(D[1])], PolynomialRing(Kq)![toQ(c) : c in Coefficients(D[2])]], D[3];
end function;

// ---------------------------------------------------------------- (1b) the class of order 16 over Q(n)
Pn<nv> := PolynomialRing(Q);
Pc2<cc,ee,ww> := ProjectiveSpace(Q, 2);
Cn := Conic(Pc2, Bcl*cc^2 - Ccl*ee^2 - (1/2)*ww^2);
okc := HasRationalPoint(Cn); assert okc;
eqs := DefiningEquations(Parametrization(Cn));
Qn<n> := FunctionField(Q);
cn := Evaluate(eqs[1], [n, 1])/Evaluate(eqs[3], [n, 1]); en := Evaluate(eqs[2], [n, 1])/Evaluate(eqs[3], [n, 1]);
assert Bcl*cn^2 - Ccl*en^2 eq 1/2;
tn := Bcl*cn^2 - 1/4; assert tn - 1/4 eq Ccl*en^2;
printf "(1b) conic B c^2 - C e^2 = 1/2: c = %o, e = %o, t = %o\n", cn, en, tn;
Py<yv> := PolynomialRing(Q);
F := AbsoluteField(NumberField([yv^2 - 2*Bcl, yv^2 + 2*Ccl]));
rho := Sqrt(F!(2*Bcl)); sig := Sqrt(F!(-2*Ccl));
FFn<m> := FunctionField(F);
embn := func< f | Evaluate(Numerator(f), m)/Evaluate(Denominator(f), m) >;
tF := embn(tn); cF := embn(cn); eF := embn(en);
JF, pF1, pF2 := Build(FFn, tF);
bF := tF^2 - 1/16;
r2 := cF*rho/2; r3 := eF*sig/2; r1 := 2*r2*r3;
assert r1^2 eq -(tF^2 - 1/16) and r2^2 eq tF/2 + 1/8 and r3^2 eq 1/8 - tF/2;
E1F := EllipticCurve([FFn| 1, -bF, -bF, 0, 0]);
xR := r1*r2 + r1*r3 + r2*r3; YR := (r1 + r2)*(r1 + r3)*(r2 + r3);
R1 := E1F![xR, YR - (xR - bF)/2];
if 2*R1 ne E1F![0,0] then R1 := -R1; end if;
assert 2*R1 eq E1F![0,0];
E2F := ChangeRing(E2, F);
D16 := 0; Q2used := 0;
for Q2 in [k*P2 + T : k in [1,3,5,7], T in [E2!0] cat T2s] do
    hs := DivisionPoints(E2F!Q2, 2);
    if #hs eq 0 then continue; end if;
    D := pF1(R1) + pF2(hs[1]);
    if &and[InQ(c) : c in Coefficients(D[1]) cat Coefficients(D[2])] then D16 := D; Q2used := Q2; break; end if;
end for;
assert D16 cmpne 0;
Jn, pn1, pn2 := Build(Qn, tn);
uvn, dg := Descend(D16, Qn);
D16n := elt< Jn | uvn[1], uvn[2], dg >;
assert Ord2(D16n) eq 16;
subn := func< f | Evaluate(Numerator(f), tn)/Evaluate(Denominator(f), tn) >;
i16 := InG0(2*D16n, subn); assert i16 gt 0;
printf "(1b) over F(n), F = Q(sqrt(%o), sqrt(%o)): R1 = half of P1 on E1, R2 = half of %o on E2; D16 = pull1(R1) + pull2(R2) has coefficients in Q(n):\n     D16 = %o\n", 2*Bcl, -2*Ccl, Q2used, D16n;
print "(1b) over Q(n): D16 has order 16 and 2 D16 lies in G0";

// ---------------------------------------------------------------- (1c) the class from the non-special 2-torsion over Q(a)
Qa<a> := FunctionField(Q);
ta := Acl*a^2;
F2 := AbsoluteField(NumberField([yv^2 + 1, yv^2 + Acl]));
ii := Sqrt(F2!(-1)); sA := Sqrt(F2!(-Acl));
FFa<w> := FunctionField(F2);
tW := Acl*w^2;
JW, pW1, pW2 := Build(FFa, tW);
eW := [tW^2 - 1/16, -tW/2 - 1/8, tW/2 - 1/8]; bW := tW^2 - 1/16;
E1W := EllipticCurve([FFa| 1, -bW, -bW, 0, 0]);
// half of T = (e_2, 0): r1 = sqrt(e2 - e1) = i (t + 1/4), r2 = 0, r3 = sqrt(e2 - e3) = sqrt(-t) = w sqrt(-A)
s1 := ii*(tW + 1/4); s3 := w*sA;
assert s1^2 eq eW[2] - eW[1] and s3^2 eq eW[2] - eW[3];
xS := eW[2] + s1*s3; YS := s1*s3*(s1 + s3);
S1 := E1W![xS, YS - (xS - bW)/2];
TW := E1W![eW[2], -(eW[2] - bW)/2];
if 2*S1 ne TW then S1 := -S1; end if;
assert 2*S1 eq TW;
E2F2 := ChangeRing(E2, F2);
D2 := 0; Q2b := 0;
for Q2 in [k*(2*P2) + T : k in [0..3], T in T2s] do
    hs := DivisionPoints(E2F2!Q2, 2);
    if #hs eq 0 then continue; end if;
    D := pW1(S1) + pW2(hs[1]);
    if &and[InQ(c) : c in Coefficients(D[1]) cat Coefficients(D[2])] then D2 := D; Q2b := Q2; break; end if;
end for;
assert D2 cmpne 0;
Ja, pa1, pa2 := Build(Qa, ta);
uva, dga := Descend(D2, Qa);
D2a := elt< Ja | uva[1], uva[2], dga >;
suba := func< f | Evaluate(Numerator(f), ta)/Evaluate(Denominator(f), ta) >;
o2 := Ord2(D2a); assert o2 in {2, 4, 8, 16};
i2 := InG0(2*D2a, suba); assert i2 gt 0;
printf "(1c) over F'(a), F' = Q(i, sqrt(%o)), t = %o a^2: half of (e_2,0) on E1 plus half of %o on E2 gives D2 with coefficients in Q(a):\n     D2 = %o (order %o; 2 D2 in G0)\n", -Acl, Acl, Q2b, D2a, o2;

// ---------------------------------------------------------------- (1d) specialization at the point of 2.2.2.16.a
t0 := -4333968559/38784504316;
okA0, a0 := IsSquare(t0/Acl); okB0, c0 := IsSquare((t0 + 1/4)/Bcl); okC0, e0 := IsSquare((t0 - 1/4)/Ccl);
assert okA0 and okB0 and okC0;
// the parameter n0 of (c0, e0) (or a sign change of it) on the conic; a0 for (1c)
function ConicParam(cv, ev)
    g := GCD(Evaluate(eqs[1], [nv, 1]) - cv*Evaluate(eqs[3], [nv, 1]), Evaluate(eqs[2], [nv, 1]) - ev*Evaluate(eqs[3], [nv, 1]));
    if Degree(g) eq 1 then return true, -Coefficient(g, 0)/Coefficient(g, 1); end if;
    return false, _;
end function;
P0<x0> := PolynomialRing(Q);
f0 := P0![Evaluate(c, t0) : c in Coefficients(hc)];
J0 := Jacobian(HyperellipticCurve(f0));
function Spec(D, val)        // specialize a class over Q(n) / Q(a) / Q(t) at a rational value
    ev := func< z | Evaluate(Numerator(z), val)/Evaluate(Denominator(z), val) >;
    return elt< J0 | P0![ev(c) : c in Coefficients(D[1])], P0![ev(c) : c in Coefficients(D[2])], D[3] >;
end function;
gens0 := [Spec(g, t0) : g in [gA, gT, gB, gU]];
done16 := false;
for sc in [1, -1] do
    for se in [1, -1] do
        if done16 then continue; end if;
        okn, n0 := ConicParam(sc*c0, se*e0);
        if not okn then continue; end if;
        try Append(~gens0, Spec(D16n, n0)); done16 := true; catch e; end try;
    end for;
end for;
assert done16;
done2 := false;
for sa in [1, -1] do
    try Append(~gens0, Spec(D2a, sa*a0)); done2 := true; catch e; end try;
    if done2 then break; end if;
end for;
assert done2;
G := {J0!0};
repeat
    old := #G;
    G := G join {g + h : g in G, h in gens0};
until #G eq old;
n2 := #{g : g in G | 2*g eq J0!0};
ex := Max([Ord2(g) : g in G]);
printf "(1d) at t0 = %o (a0 = %o, c0 = %o, e0 = %o): the specialized generators generate %o classes, %o killed by 2, exponent %o\n", t0, a0, c0, e0, #G, n2, ex;
assert #G eq 128 and n2 eq 16 and ex eq 16;
print "(1d) so G = <G0, D16, D2> is Z/2 x Z/2 x Z/2 x Z/16, marked over K = Q(D_u), and it embeds in J_{t(R)}(Q)_tors for every R in D_u(Q) with C_{t(R)} smooth";

// ---------------------------------------------------------------- (2) the residue class
E := EllipticCurve([0,0,0,-dd^2,0]);
P := E![-47963453395683600/226291849, -25397917058104776424224600/3404108284507];
T3 := E![dd, 0];
assert Order(P) eq 0;                         // P has infinite order (E_d(Q)_tors = E_d[2])
R0 := P + T3;
assert R0[1]/(4*dd) eq t0;
function Del(R) return [Sqf(R[1]), Sqf(R[1] + dd), Sqf(R[1] - dd)]; end function;
assert Del(R0) eq [Sqf(dd*Acl), Sqf(dd*Bcl), Sqf(dd*Ccl)];
primes := [11, 13];
ords := [Order(ChangeRing(E, GF(p))!P) : p in primes];
assert ords eq [6, 4];
nstep := 12;
assert &and[nstep mod o eq 0 : o in ords] and IsEven(nstep);
for p in primes do
    assert (2*dd) mod p ne 0;
    assert Valuation(R0[1], p) ge 0;           // R0 is not O mod p
end for;
printf "(2) E_d, d = %o; P has order %o mod 11, 13; R_k = R0 + %o k P, R0 = P + (d,0), x(R0) = 4 d t0\n", dd, ords, nstep;
R1k := R0 + nstep*P;
t1k := R1k[1]/(4*dd);
for p in primes do assert GF(p)!t1k eq GF(p)!t0; end for;
assert IsSquare(t1k/Acl) and IsSquare((t1k + 1/4)/Bcl) and IsSquare((t1k - 1/4)/Ccl);
printf "(2) sanity k = 1: t(R_1) has %o digits, = t0 mod 11 and mod 13, and satisfies the three locus conditions\n", #Sprint(Numerator(t1k));
for p in primes do
    r0 := Z!(GF(p)!t0);
    fbar := ReductionClass(hc, p, r0);
    assert PolynomialRing(GF(p))!SpecPoly(hc, t0) eq fbar;
    assert PolynomialRing(GF(p))!SpecPoly(hc, t1k) eq fbar;
    nJ, chi := FrobData(fbar);
    Ap := AbelianGroup(Jacobian(HyperellipticCurve(fbar)));
    printf "(2) p = %o, t = %o mod %o: reduction %o (squarefree, degree %o); J(F_%o) = %o (order %o), chi = %o\n",
        p, r0, p, fbar, Degree(fbar), p, Invariants(Ap), nJ, chi;
    assert nJ eq 128 and #Ap eq 128 and Invariants(Ap) eq [2,2,2,16];
end for;
print "(2) every C_{t(R_k)} has the reductions of C_{t0} at 11 and 13, #J(F_11) = 128 = #G: J(Q)_tors is exactly [2,2,2,16]";

// ---------------------------------------------------------------- (3) infinitely many isomorphism classes
R3 := 3*P + T3;
assert R3 - R0 eq 2*P;
t3 := R3[1]/(4*dd);
assert IsSquare(t3/Acl) and IsSquare((t3 + 1/4)/Bcl) and IsSquare((t3 - 1/4)/Ccl);
f3 := SpecPoly(hc, t3);
assert G2Invariants(HyperellipticCurve(f0)) ne G2Invariants(HyperellipticCurve(f3));
printf "(3) t(3P + (d,0)) has %o digits; the G2-invariants of C_{t0} and C_{t(3P+T3)} differ\n", #Sprint(Numerator(t3));

// ---------------------------------------------------------------- (4) sanity
inv0 := TorsionInvariants(f0);
printf "(4) TorsionSubgroup at t0: %o\n", inv0;
assert inv0 eq [2,2,2,16];
if Tors2 eq 1 then
    inv3 := TorsionInvariants(f3);
    printf "(4) TorsionSubgroup at t(3P+T3): %o\n", inv3;
    assert inv3 eq [2,2,2,16];
end if;
for tv in [t0, t3] do
    fv := SpecPoly(hc, tv);
    Hq0 := PolynomialRing(Q)![Coefficient(fv, 2*i) : i in [0..3]];
    Hr0 := PolynomialRing(Q)![Coefficient(fv, 6 - 2*i) : i in [0..3]];
    E1v := EllipticCurve([Q| 1, -(tv^2-1/16), -(tv^2-1/16), 0, 0]);
    assert IsIsomorphic(MonEC(Hq0), E2) and IsIsomorphic(MonEC(Hr0), E1v);
end for;
print "(4) at t0 and t(3P+T3) the two quotients are Q-isomorphic to E2 and E1(t)";
// the census curve 2.2.2.16.a (its even model y^2 = 4f + h^2)
fB := P0![921049431076909147251530730180268944, -321215744373328039928102646572568, -1625435349401383665152370863687,
          997294032762556482905418146, 38866477002321559938721, -169349659544271600];
assert IsIsomorphic(HyperellipticCurve(f0), HyperellipticCurve(fB));
print "(4) C_{t0} is Q-isomorphic to the census curve 2.2.2.16.a";
printf "time %o s\n", Cputime(T0);
print "FAMILY_22216_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E1 x E2) and J(Q)_tors exactly [2,2,2,16]";
quit;
