/* family_2_2_16_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E1 x E2) with J(Q)_tors exactly [2,2,16].
   Found 2026-10-19 (torsion_jac notes/pow2_2026_10_19.md).

   Source.  The Howe-Leprevost-Poonen (2,2)-gluing (HLP 2000, Prop. 4) of E1 = E_{2,4}(t) with E2 = 210e2
   (torsion Z/2 x Z/8) along alpha = identity on the labelled roots (special roots matched: "sp -> sp"), with gain
   g = 2 (torsion_jac notes gain2_loci_2026_10_15.md, section 2: [2,4] x [2,8], g = 2, sp -> sp gives [2,2,16]).
   This is the construction of family_22216_split_certificate.m with one of its three locus conditions dropped
   (there E2 = E_{2,8}(-20/31) and t = A a^2 is imposed as well, which adds a fourth independent 2-torsion class).
       E1 = E_{2,4}(t) : y^2 + xy - b y = x^3 - b x^2, b = t^2 - 1/16, P1 = (0,0) of order 4, monic roots
            (t^2 - 1/16, -t/2 - 1/8, t/2 - 1/8) (special root first);
       E2 = 210e2 : y^2 + xy = x^3 - 1070x + 7812, P8 = (244, 3658) of order 8, monic roots (28, -36, 31/4)
            (special root x(4 P8) = 28 first).
   For every t the gluing is (derived here; HLP's formula, checked below)
       C_t : y^2 = f_t(x) = -(28t - 25)(100t - 7)(324x^2 - (4t - 1)^2)(1024x^2 - (4t + 1)^2)(175x^2 - 4t),
   and over Q(t) it carries [2,2,8].  The gain: the Kummer class of P1 is (1/16 - t^2, (t + 1/4)/2, (1/4 - t)/2)
   and that of Q2 = P8 + (-36, 18) is (-6, 10, -15) mod squares, so they agree iff
       t + 1/4 = 5 c^2  and  t - 1/4 = 30 e^2,
   a conic with the point (c, e) = (2/5, 1/10) (t = 11/20).  Lines through it give
       c = (12n^2 - 6n + 2)/(30n^2 - 5),  t = t(n) = 5c^2 - 1/4,  and t(n) - 1/4 = 30 e(n)^2 (checked below).
   The member n = 2/3 (t = 11/20, E1(11/20) = 330b2) is the census curve 2.2.16.a.

   Argument.
   (1) Splitting.  phi1 : (x, y) -> (x^2, y) and phi2 : (x, y) -> (1/x^2, y/x^3) map C_t onto Ea : Y^2 = F(X)
       (f_t = F(x^2)) and Eb : Y^2 = X^3 F(1/X), which are Q(t)-isomorphic to E2 and E1(t) (checked below).
       Wherever f_t is squarefree of degree 6, phi1^* + phi2^* : Ea x Eb -> J is an isogeny over Q with kernel of
       order 4, so J is isogenous over Q to E1(t) x E2.
   (2) The marked group.  Over Q(n): gA = pull(P1) (order 4), gB = pull(P8) (order 8), Z1 = pull((-36, 18)),
       Z2 = gA - 2 gB (2 gA = pull(2 P1) = pull(4 P8) = 4 gB because the gluing matches the special points).
       Over F(n), F = Q(sqrt 10, sqrt -15): the half R1 of P1 on E1 has x(R1) = r1 r2 + r1 r3 + r2 r3 with
       r2 = c sqrt(10)/2, r3 = e sqrt(-15), r1 = +-2 r2 r3 (r_i^2 = -e_i), and the halves of Q2 on E2 are defined
       over F.  D16 = pull(R1) + pull(R2) is computed by Mumford arithmetic over F(n); its coefficients lie in Q(n)
       (Galois descent made explicit).  Over Q(n): D16 has order 16, 2 D16 = gA + pull(Q2), and the 64 sums
       i D16 + j Z1 + k Z2 (0 <= i < 16, 0 <= j, k < 2) are distinct, so G = <D16, Z1, Z2> = Z/16 x Z/2 x Z/2 is
       marked over Q(n) and embeds in J(Q)_tors at every smooth member.
   (3) Exactness on a residue class.  #J(F_p) = #E1(F_p) #E2(F_p) is divisible by 8 * 16 = 128 at every good p,
       so the group structure is used.  For n = 5 mod 13 and n = 7 mod 19 (n 13- and 19-integral) the curve has
       fixed reductions with J(F_13) = [2,8,16] (2-part partition (4,3,1), odd part 1) and J(F_19) = [2,2,2,48]
       ((4,1,1,1), odd part 3).  The termwise minimum of the partitions is (4,1,1), that of [2,2,16], and the odd
       parts are coprime; J(Q)_tors embeds in both and contains G, so it is exactly [2,2,16].  The census curve
       (n = 2/3) lies in this class.
   (4) Infinitely many.  The G2-invariants of the members differ, so n -> M_2 is nonconstant with finite fibres,
       and the residue class gives infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup at n = 83 and n = -164 is [2,2,16]; n = 2/3 is Q-isomorphic to 2.2.16.a.
   Run (from the repository root):
       magma -b pipeline/magma/family_2_2_16_split_certificate.m > family_2_2_16_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

function MonEC(G)       // Y^2 = c3 X^3 + c2 X^2 + c1 X + c0  ->  monic model (x, y) = (c3 X, c3 Y)
    c := [ Coefficient(G, i) : i in [0..3] ];
    return EllipticCurve([0, c[3], 0, c[2]*c[4], c[1]*c[4]^2]), c[4];
end function;
function Ord2(g)      // order of a class killed by 64, else -1
    z := Parent(g)!0; h := g; k := 1;
    while k le 64 do if h eq z then return k; end if; h := 2*h; k := 2*k; end while;
    return -1;
end function;
function Sqf(q)
    q := Rationals()!q; n := Numerator(q)*Denominator(q);
    return Sign(n)*&*([1] cat [p[1] : p in Factorization(Abs(n)) | p[2] mod 2 eq 1]);
end function;
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

// ---------------------------------------------------------------- (0) E2, the gluing over Q(t), the conic
E2 := EllipticCurve([Q| 1, 0, 0, -1070, 7812]);
assert Conductor(E2) eq 210 and Invariants(TorsionSubgroup(E2)) eq [2, 8];
P8 := E2![244, 3658]; T2a := E2![-36, 18]; T2b := E2![31/4, -31/8]; T2s := E2![28, -14];
assert Order(P8) eq 8 and 4*P8 eq T2s and 2*T2a eq E2!0 and 2*T2b eq E2!0;
be := [Q| 28, -36, 31/4];
Q2 := P8 + T2a;
assert [Sqf(Q2[1] - r) : r in be] eq [-6, 10, -15];
Ft<t> := FunctionField(Q); Px<x> := PolynomialRing(Ft);
al := [t^2 - 1/16, -t/2 - 1/8, t/2 - 1/8];
assert [-a : a in al] eq [1/16 - t^2, (t + 1/4)/2, (1/4 - t)/2];
print "E2 = 210e2, monic roots (28, -36, 31/4) (special first); Q2 = P8 + (-36,18) has Kummer class (-6, 10, -15)";
print "E1 = E_{2,4}(t), monic roots (t^2-1/16, -t/2-1/8, t/2-1/8); P1 = (0,0) has Kummer class (1/16 - t^2, (t+1/4)/2, (1/4-t)/2)";
aa1 := (al[3]-al[2])^2/(be[3]-be[2]) + (al[2]-al[1])^2/(be[2]-be[1]) + (al[1]-al[3])^2/(be[1]-be[3]);
bb1 := (be[3]-be[2])^2/(al[3]-al[2]) + (be[2]-be[1])^2/(al[2]-al[1]) + (be[1]-be[3])^2/(al[1]-al[3]);
aa2 := al[1]*(be[3]-be[2]) + al[2]*(be[1]-be[3]) + al[3]*(be[2]-be[1]);
bb2 := be[1]*(al[3]-al[2]) + be[2]*(al[1]-al[3]) + be[3]*(al[2]-al[1]);
Df := ((al[1]-al[2])*(al[1]-al[3])*(al[2]-al[3]))^2;
Dg := ((be[1]-be[2])*(be[1]-be[3])*(be[2]-be[3]))^2;
AA := Dg*aa1/aa2; BB := Df*bb1/bb2;
hHLP := -(AA*(al[2]-al[1])*(al[1]-al[3])*x^2 + BB*(be[2]-be[1])*(be[1]-be[3]))
        *(AA*(al[3]-al[2])*(al[2]-al[1])*x^2 + BB*(be[3]-be[2])*(be[2]-be[1]))
        *(AA*(al[1]-al[3])*(al[3]-al[2])*x^2 + BB*(be[1]-be[3])*(be[3]-be[2]));
ft := -(28*t - 25)*(100*t - 7)*(324*x^2 - (4*t - 1)^2)*(1024*x^2 - (4*t + 1)^2)*(175*x^2 - 4*t);
oksq, rt := IsSquare(Ft!(hHLP/ft));
assert oksq;
print "C_t : y^2 = f_t(x) = -(28t - 25)(100t - 7)(324x^2 - (4t-1)^2)(1024x^2 - (4t+1)^2)(175x^2 - 4t)";
printf "the HLP (Prop. 4) curve over Q(t) is (%o)^2 f_t\n", rt;
// the two quotients, over Q(t)
F3 := Px![ Coefficient(ft, 2*k) : k in [0..3] ];
Fr := Px![ Coefficient(F3, 3 - k) : k in [0..3] ];
assert Evaluate(F3, x^2) eq ft;
Ea, ca := MonEC(F3); Eb, cb := MonEC(Fr);
bt := t^2 - 1/16;
E1t := EllipticCurve([Ft| 1, -bt, -bt, 0, 0]);
okA, mA := IsIsomorphic(ChangeRing(E2, Ft), Ea); okB, mB := IsIsomorphic(E1t, Eb);
assert okA and okB;
print "(1) phi1 = (x^2, y): C_t -> Ea ~ E2 and phi2 = (1/x^2, y/x^3): C_t -> Eb ~ E1(t) over Q(t): J is isogenous over Q to E1(t) x E2";
// the conic t + 1/4 = 5 c^2, t - 1/4 = 30 e^2
Qn<n> := FunctionField(Q); Pn<xn> := PolynomialRing(Qn);
cn := (12*n^2 - 6*n + 2)/(30*n^2 - 5);
tn := 5*cn^2 - 1/4;
oke, en := IsSquare((tn - 1/4)/30);
assert oke;
printf "conic: c(n) = %o, e(n) = %o, t(n) = %o\n", cn, en, tn;
assert Evaluate(tn, 2/3) eq 11/20;

// ---------------------------------------------------------------- (2) the marked group over Q(n)
subn := func< c | Evaluate(Numerator(c), tn)/Evaluate(Denominator(c), tn) >;
f := Pn![ subn(c) : c in Coefficients(ft) ];
dA := IsomorphismData(mA); dB := IsomorphismData(mB);
P1t := E1t![0, 0];
im := mB(P1t); r, ss, tt, uu := Explode(dB);
conv := (uu^2*P1t[1] + r eq im[1] and uu^3*P1t[2] + ss*uu^2*P1t[1] + tt eq im[2]) select 2 else 1;
if conv eq 1 then assert (P1t[1]-r)/uu^2 eq im[1] and (P1t[2] - ss*(P1t[1]-r) - tt)/uu^3 eq im[2]; end if;
function ApplyIso(dat, pt, sub)
    r, ss, tt, uu := Explode([sub(c) : c in dat]);
    if conv eq 1 then return [(pt[1]-r)/uu^2, (pt[2] - ss*(pt[1]-r) - tt)/uu^3]; end if;
    return [uu^2*pt[1] + r, uu^3*pt[2] + ss*uu^2*pt[1] + tt];
end function;
// the family and the two pullbacks over a field FF containing tF (the image of t)
function Build(FF, tF)
    sub := func< c | Evaluate(Numerator(c), tF)/Evaluate(Denominator(c), tF) >;
    PF<xF> := PolynomialRing(FF);
    JJ := Jacobian(HyperellipticCurve(PF![ sub(c) : c in Coefficients(ft) ]));
    cA := sub(ca); cB := sub(cb);
    pull1 := function(pt)        // E1 -> Eb -> J through (1/x^2, y/x^3)
        if pt eq Parent(pt)!0 then return JJ!0; end if;
        im := ApplyIso(dB, [pt[1], pt[2]], sub); Xv := im[1]/cB; Yv := im[2]/cB;
        return elt< JJ | xF^2 - 1/Xv, (Yv/Xv)*xF >;
    end function;
    pull2 := function(pt)        // E2 -> Ea -> J through (x^2, y)
        if pt eq Parent(pt)!0 then return JJ!0; end if;
        im := ApplyIso(dA, [FF!pt[1], FF!pt[2]], sub); Xv := im[1]/cA; Yv := im[2]/cA;
        return elt< JJ | xF^2 - Xv, PF!Yv >;
    end function;
    return JJ, pull1, pull2;
end function;
J, pl1, pl2 := Build(Qn, tn);
assert J eq Jacobian(HyperellipticCurve(f));
E1n := EllipticCurve([Qn| 1, -(tn^2 - 1/16), -(tn^2 - 1/16), 0, 0]);
gA := pl1(E1n![0, 0]); gB := pl2(P8); Z1 := pl2(T2a);
Z2 := gA - 2*gB;
assert [Ord2(gA), Ord2(gB), Ord2(Z1), Ord2(Z2)] eq [4, 8, 2, 2];
alt := [tn^2 - 1/16, -tn/2 - 1/8, tn/2 - 1/8];
T1n := [ E1n![a, -(a - (tn^2 - 1/16))/2] : a in alt ];
T2l := [ T2s, T2a, T2b ];
for k in [1..3] do assert pl1(T1n[k]) + pl2(T2l[k]) eq J!0; end for;
print "(2) over Q(n): gA = pull(P1), gB = pull(P8) of orders 4, 8; pull(T) + pull(alpha T) = 0 on 2-torsion (alpha = id, sp -> sp), so 2 gA = 4 gB";
Py := RX;
F := AbsoluteField(NumberField([X^2 - 10, X^2 + 15]));
r10 := Sqrt(F!10); r15 := Sqrt(F!(-15));
FFn<m> := FunctionField(F);
embn := func< c | Evaluate(Numerator(c), m)/Evaluate(Denominator(c), m) >;
tF := embn(tn); cF := embn(cn); eF := embn(en);
JF, pF1, pF2 := Build(FFn, tF);
bF := tF^2 - 1/16;
E1F := EllipticCurve([FFn| 1, -bF, -bF, 0, 0]);
r2 := cF*r10/2; r3 := eF*r15;
assert r2^2 eq tF/2 + 1/8 and r3^2 eq 1/8 - tF/2;
R1 := 0;
for sg in [1, -1] do
    r1 := sg*2*r2*r3;
    assert r1^2 eq -(tF^2 - 1/16);
    xR := r1*r2 + r1*r3 + r2*r3; YR := (r1 + r2)*(r1 + r3)*(r2 + r3);
    Rc := E1F![xR, YR - (xR - bF)/2];
    if 2*Rc eq E1F![0, 0] then R1 := Rc; break; end if;
    if 2*(-Rc) eq E1F![0, 0] then R1 := -Rc; break; end if;
end for;
assert R1 cmpne 0;
function InQ(c)
    nu := Numerator(c); de := Denominator(c); lc := LeadingCoefficient(de);
    return &and[ IsCoercible(Q, a) : a in Coefficients(nu/lc) cat Coefficients(de/lc) ];
end function;
PQn := PolynomialRing(Q);
function ToQn(c)
    nu := Numerator(c); de := Denominator(c); lc := LeadingCoefficient(de);
    return Evaluate(PQn![Q!a : a in Coefficients(nu/lc)], n)/Evaluate(PQn![Q!a : a in Coefficients(de/lc)], n);
end function;
E2F := ChangeRing(E2, F);
D16F := 0; Qused := 0;
for k in [1, 3, 5, 7] do
    for T in [E2!0, T2a, T2b, T2s] do
        hs := DivisionPoints(E2F!Eltseq(k*P8 + T), 2);
        if #hs eq 0 then continue; end if;
        D := pF1(R1) + pF2(hs[1]);
        if &and[ InQ(c) : c in Coefficients(D[1]) cat Coefficients(D[2]) ] then D16F := D; Qused := k*P8 + T; break k; end if;
    end for;
end for;
assert D16F cmpne 0;
D16 := elt< J | Pn![ ToQn(c) : c in Coefficients(D16F[1]) ], Pn![ ToQn(c) : c in Coefficients(D16F[2]) ], D16F[3] >;
assert Ord2(D16) eq 16 and 2*D16 eq gA + pl2(Qused);
GG := { i*D16 + j*Z1 + l*Z2 : i in [0..15], j, l in [0..1] };
assert #GG eq 64;
printf "(2) over F(n), F = Q(sqrt 10, sqrt -15): R1 = half of P1 on E1, R2 = half of %o on E2; D16 = pull(R1) + pull(R2) has coefficients in Q(n):\n     D16 = %o\n", Qused, D16;
print "    over Q(n): D16 has order 16, 2 D16 = gA + pull(Q2), and the 64 sums i D16 + j Z1 + k Z2 are distinct: G = [2,2,16] is marked";

// ---------------------------------------------------------------- (3) exactness on a residue class
primes := [13, 19]; residues := [5, 7];
nres := CRT(residues, primes); M := &*primes;
members := [Q| 2/3, nres, nres - M, 2/3 + M];
invs := [];
for k in [1..#primes] do
    p := primes[k]; r0 := residues[k];
    fbar := ReductionClass(f, p, r0);
    for mm in members do
        assert Valuation(mm, p) ge 0 and GF(p)!mm eq GF(p)!r0;
        assert PolynomialRing(GF(p))!SpecPoly(f, mm) eq fbar;
    end for;
    nJ, chi := FrobData(fbar);
    Ap := AbelianGroup(Jacobian(HyperellipticCurve(fbar)));
    assert #Ap eq nJ;
    Append(~invs, Invariants(Ap));
    printf "(3) p = %o, n = %o mod %o: reduction %o (squarefree, degree %o), J(F_%o) = %o (order %o), chi = %o\n",
        p, r0, p, fbar, Degree(fbar), p, Invariants(Ap), nJ, chi;
end for;
assert invs eq [[2, 8, 16], [2, 2, 2, 48]];
bnd := StructureBound(invs);
printf "(3) the largest group embedding in all J(F_p) has invariants %o\n", bnd;
assert bnd eq [2, 2, 16];
printf "(3) every n = %o mod %o, 13- and 19-integral, has J(Q)_tors = G = [2,2,16]\n", nres, M;

// ---------------------------------------------------------------- (4) infinitely many, (5) sanity
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, mm))) : mm in members } eq #members;
printf "(4) the G2-invariants of the members %o are pairwise distinct\n", members;
for mm in [members[2], members[3]] do
    invm := TorsionInvariants(SpecPoly(f, mm));
    printf "(5) TorsionSubgroup at n = %o: %o\n", mm, invm;
    assert invm eq [2, 2, 16];
end for;
C2216a := HyperellipticCurve(RX![9, -33, 3, 61, 3, -33, -11], RX![0, 1, 1]);     // census 2.2.16.a
assert IsIsomorphic(HyperellipticCurve(SpecPoly(f, 2/3)), C2216a);
E1c := EllipticCurve([Q| 1, -(11^2/20^2 - 1/16), -(11^2/20^2 - 1/16), 0, 0]);
assert CremonaReference(MinimalModel(E1c)) eq "330b2";
print "(5) the member n = 2/3 (t = 11/20, E1 = 330b2) is Q-isomorphic to the census curve 2.2.16.a";
printf "time %o s\n", Cputime(T0);
print "FAMILY_2_2_16_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E1 x E2) and J(Q)_tors exactly [2,2,16]";
quit;
