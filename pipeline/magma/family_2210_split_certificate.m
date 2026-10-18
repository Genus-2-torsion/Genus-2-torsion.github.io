/* family_2210_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E x F_v) with J(Q)_tors exactly [2,2,10].
   Found 2026-10-18 (torsion_jac notes/g33fam_2026_10_18.md); the member v = 0 is the census curve 2.2.10.b
   (LMFDB 30492.q.1), which is the (3,3)-gluing of 66c1 and 462f2.

   Source.  (3,3)-gluings of the fixed curve E = 66c1 (torsion Z/10) with partners F that have full rational
   2-torsion, in the normal form of Broker-Howe-Lauter-Stevenhagen (arXiv:1403.6911, Alg. 5.4 and App. A,
   normalization b = a):
       C : y^2 = tw h1 h2,   h1 = x^3 + 3ax + 2a,   h2 = 2dx^3 + 3cx^2 + 1,   d = (1 - 12ac)/(16a),
   with the degree-3 maps phi1 : (x,y) -> (nn1/h1, y mm1/h1^2) onto Y^2 = tw f1(X) and phi2 : (x,y) ->
   (nn2/h2, y mm2/h2^2) onto Y^2 = tw f2(X), certified by the identities h1^4 f1(nn1/h1) = h1 h2 mm1^2 and
   h2^4 f2(nn2/h2) = h1 h2 mm2^2.  w parametrizes the rational curve j(E0(a,c)) = j(E), v the rational cover
   "f2 has a rational root X(v)", with v = 0 at the BHLS point (a, c) = (-2197/1225, -775/10816) of 66c1 x 462f2.
   tw is the twist with E0^tw = E.  F_v = F0^tw has full rational 2-torsion exactly when v lies on the genus-1 curve
       D : s^2 = -2v^4 + 19v^3 - 12v^2 - 80v + 196,
   whose Jacobian is EM : y^2 = x^3 + 39204 (a twist of 36a1), with the point P = (616, -15290) of infinite order.

   Argument.
   (1) Splitting.  As in family_214_split_certificate.m: at every member of the class below, phi1, phi2 are
       morphisms over Q onto E and F_v, and J is isogenous over Q to E x F_v.
   (2) The marked group.  phi_* phi^* = [3] makes phi1^* injective on E(Q)_tors = Z/10 and phi2^* injective on
       F_v(Q)[2] = (Z/2)^2.  #E(F_23) = 30 != 24 = #F_v(F_23), so E and F_v are not isogenous and
       phi1_* phi2^* in Hom(F_v, E) is 0; if phi1^*(A) = phi2^*(B) with A, B of order prime to 3, applying phi1_*
       gives 3A = 0, so A = 0 and then B = 0.  Hence Z/10 x (Z/2)^2 = [2,2,10] embeds in J_v(Q).  At the members
       below this is also checked by Mumford arithmetic (40 classes).
   (3) Exactness on a residue class.  For every 23- and 47-integral v = 0 (mod 1081), C_v has the reductions of
       C_0 at 23 and 47, with #J(F_23) = 720 and #J(F_47) = 2200 (each equal to #E(F_p) #F(F_p)).  J_v(Q)_tors
       injects into both, so its order divides gcd(720, 2200) = 40: it is exactly [2,2,10].
   (4) Infinitely many members.  D -> EM sends B = (0, -14) to O and R0 = (0, 14) to (792/49, -71478/343).
       P has order 12 mod 23 and 48 mod 47, so R_k = R0 + 48 k P reduces to R0 mod 23 and 47; v = NX/NZ with
       NZ(R0) a unit, so v(R_k) is 23- and 47-integral and = 0: each R_k gives a member.  The R_k are distinct,
       v has degree 2 on D, and the member v2 = 5560/3033 (over R0 + 2P) has G2-invariants different from those
       of C_0, so v -> M_2 is nonconstant on D: infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup at v = 0 and v = 5560/3033 is [2,2,10]; C_0 is Q-isomorphic to the census curve
       2.2.10.b.
   Run (from the repository root):
       magma -b pipeline/magma/family_2210_split_certificate.m > data/knowledge/sources/family_2210_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

// ---------------------------------------------------------------- BHLS (3,3) normal form, b = a (arXiv:1403.6911, Alg. 5.4, App. A)
function BHLSParts(aa, cc)
    K := Parent(aa); P<x> := PolynomialRing(K);
    bb := aa; dd := (1 - 12*aa*cc)/(16*aa);
    D1 := aa^3 + bb^2; D2 := cc^3 + dd^2;
    f1 := x^3 + 12*(2*aa^2*dd - bb*cc)*x^2 + 12*(16*aa*dd^2 + 3*cc^2)*D1*x + 512*D1^2*dd^3;
    f2 := x^3 + 12*(2*bb*cc^2 - aa*dd)*x^2 + 12*(16*bb^2*cc + 3*aa^2)*D2*x + 512*D2^2*bb^3;
    h1 := x^3 + 3*aa*x + 2*bb; h2 := 2*dd*x^3 + 3*cc*x^2 + 1;
    nn1 := 12*D1*(-2*dd*x + cc); mm1 := D1*(16*dd*x^3 - 12*cc*x^2 - 1);
    nn2 := 12*D2*x^2*(aa*x - 2*bb); mm2 := D2*(x^3 + 12*aa*x - 16*bb);
    return h1, h2, f1, f2, nn1, mm1, nn2, mm2, dd, D1, D2;
end function;
function CoverIdentity(f, nn, hh, mm, h1, h2)      // hh^4 f(nn/hh) = h1 h2 mm^2, f a monic cubic
    c := Coefficients(f);
    return nn^3*hh + c[3]*nn^2*hh^2 + c[2]*nn*hh^3 + c[1]*hh^4 eq h1*h2*mm^2;
end function;
function TwistedCurve(f, tw)                       // Y^2 = tw f(X) in Weierstrass form (X' = tw X, Y' = tw Y)
    c := Coefficients(f);
    return EllipticCurve([0, c[3]*tw, 0, c[2]*tw^2, c[1]*tw^3]);
end function;
function TwistTo(f1, E)                            // tw with c4, c6 of E0^tw equal to those of E
    c40, c60 := Explode(cInvariants(EllipticCurve(f1))); c4, c6 := Explode(cInvariants(E));
    return (c6/c60)/(c4/c40);
end function;
function SqfRep(z)                                 // squarefree polynomial representative of z mod squares, and the root of z/rep
    num := Numerator(z)*Denominator(z); K := Parent(z); R := Parent(num);
    lc := LeadingCoefficient(num); n := Numerator(lc)*Denominator(lc);
    sq := &*([Integers()| 1] cat [p[1]^(p[2] div 2) : p in Factorization(Abs(n))]);
    rep := (n div sq^2)*&*([R| 1] cat [p[1] : p in Factorization(num) | IsOdd(p[2])]);
    ok, rt := IsSquare(K!(z/rep)); assert ok;
    return rep, rt;
end function;
function Hom3(f, CR, wt)                           // homogenize a polynomial in x for the weighted ring (1, 3, 1)
    X := CR.1; Z := CR.3;
    return &+([CR| 0] cat [Coefficient(f, i)*X^i*Z^(wt - i) : i in [0..Degree(f)]]);
end function;
function OrderUpTo(P, n)
    for k in [1..n] do if k*P eq Parent(P)!0 then return k; end if; end for; return -1;
end function;


// ---------------------------------------------------------------- (0) the fixed curve E
E := EllipticCurve([1,0,0,-45,81]);
assert Conductor(E) eq 66 and Invariants(TorsionSubgroup(E)) eq [10];
TE, mTE := TorsionSubgroup(E);
printf "E = 66c1 = %o, torsion %o\n", aInvariants(E), Invariants(TE);

// ---------------------------------------------------------------- (1) the family over Q(v)
Fw<w> := FunctionField(Q); Fv := Fw; v := w;
aW := (-5799431716/6731694363*w^6 + 97934044/11010295*w^5 - 402644/10805*w^4 + 499319908756/6171329775*w^3 - 2889836880657408/30126374851625*w^2 + 130017249758196006912/2206003798510240625*w - 474729982221303869865984/32306925629182473953125)/(w^4 - 9241706372/1234265955*w^3 + 1575442112633856/78328574614225*w^2 - 25679923358150688768/1147121975225325125*w + 23261769128843889623433216/2839132624292555811000625);
cW := (-2243898121/23197726864*w^2 + 116102084048352/276031201812665*w - 24950882719015197696/52552200357104226025)/(w^4 - 1310434/190385*w^3 + 615915249121/36246448225*w^2 - 9405182997554688/530829234255125*w + 51511499806999117824/7773994135666305625);
wV := Evaluate((391296/190385*w^3 - 3521664/190385*w^2 + 3130368/190385*w)/(w^3 - 49/6*w^2 + 52/9*w + 40/9), w);
XV := Evaluate((2673/32*w^16 - 180873/64*w^15 + 20708325/512*w^14 - 164557305/512*w^13 + 3190535601/2048*w^12 - 597596329/128*w^11 + 1469816315/192*w^10 - 14929475/24*w^9 - 235195125/8*w^8 + 74953868*w^7 - 96934376*w^6 + 67881792*w^5 - 17282848*w^4 - 8600064*w^3 + 7422976*w^2 - 5496832/3*w + 450560/3)/(w^18 - 97/3*w^17 + 5611/12*w^16 - 113680/27*w^15 + 2149486/81*w^14 - 9251281/81*w^13 + 22007836/81*w^12 + 19680544/81*w^11 - 441926320/81*w^10 + 2071636864/81*w^9 - 1687062784/27*w^8 + 1297641472/27*w^7 + 5303050240/27*w^6 - 61330817024/81*w^5 + 104307113984/81*w^4 - 102021988352/81*w^3 + 59633041408/81*w^2 - 19381878784/81*w + 2671771648/81), w);
aV := Evaluate(aW, wV); cV := Evaluate(cW, wV);
h1, h2, f1, f2, nn1, mm1, nn2, mm2, dV, D1, D2 := BHLSParts(aV, cV);
assert jInvariant(EllipticCurve(f1)) eq jInvariant(E);
print "(1) j(E0(a(v), c(v))) = j(E) identically in v";
assert CoverIdentity(f1, nn1, h1, mm1, h1, h2) and CoverIdentity(f2, nn2, h2, mm2, h1, h2);
print "(1) cover identities h1^4 f1(nn1/h1) = h1 h2 mm1^2 and h2^4 f2(nn2/h2) = h1 h2 mm2^2 hold over Q(v)";
tw := TwistTo(f1, E);
E0t := TwistedCurve(f1, tw); Ft := TwistedCurve(f2, tw);
assert cInvariants(E0t) eq cInvariants(E);
assert Evaluate(f2, XV) eq 0;
PX<x> := Parent(h1);
qq := f2 div (x - XV); assert f2 eq (x - XV)*qq;
drep, drt := SqfRep(Discriminant(qq));
Pv<V> := PolynomialRing(Q);
dpol := Pv!Evaluate(drep, V);
assert dpol eq Pv!Evaluate(-2*V^4 + 19*V^3 - 12*V^2 - 80*V + 196, V);
printf "(1) F = F0^tw has the 2-torsion point (tw X(v), 0); f2 = (X - X(v)) q(X), and disc(q) = dpol(v) drt(v)^2 with dpol = %o\n", dpol;
print "    so F_v has full rational 2-torsion exactly when dpol(v) is a square: the base is the genus-1 curve D : s^2 = dpol(v)";
rep, rt := SqfRep(tw);
ff := rep*h1*h2; den := LCM([Denominator(t) : t in Coefficients(ff)]); ff := ff*den^2;
gc := GCD([Pv!Evaluate(Numerator(t), V) : t in Coefficients(ff)]);
sq := &*([Pv| 1] cat [fa[1]^(fa[2] div 2) : fa in Factorization(gc)]);
sqv := Evaluate(sq, v); ff := ff/sqv^2;
assert &and[Denominator(t) eq 1 : t in Coefficients(ff)];
kk := rt*sqv/den;
assert ff*kk^2 eq tw*h1*h2;
printf "(1) C_v : y^2 = f_v(x), coefficients of degrees %o in v; y' = kk y maps C_v to y'^2 = tw h1 h2\n", [Degree(Numerator(t)) : t in Coefficients(ff)];
print "(1) phi1 : (x,y) -> (tw nn1/h1, tw kk y mm1/h1^2) onto E0^tw ~ E, phi2 : (x,y) -> (tw nn2/h2, tw kk y mm2/h2^2) onto F_v";

// ---------------------------------------------------------------- (2) exactness on the residue class v = 0 mod [23, 47]
primes := [23, 47]; r0 := 0;
finite := [aV, cV, dV, D1, D2, tw, kk, XV, drt];
nonzero := [aV, D1, D2, tw, kk, drt, Discriminant(f1), Discriminant(f2)];
nJ := []; noniso := false;
for p in primes do
    for z in finite do Nz, Dz := IntNumDen(z); assert (GF(p)!Evaluate(Dz, r0)) ne 0; end for;
    for z in nonzero do Nz, Dz := IntNumDen(z); assert (GF(p)!Evaluate(Nz, r0)) ne 0; end for;
    fbar := ReductionClass(ff, p, r0);
    n, chi := FrobData(fbar);
    Ebar := ChangeRing(E, GF(p));
    Fbar := EllipticCurve([GF(p)!Evaluate(Numerator(t), r0)/Evaluate(Denominator(t), r0) : t in aInvariants(Ft)]);
    assert n eq #Ebar*#Fbar;
    if #Ebar ne #Fbar then noniso := true; end if;
    Append(~nJ, n);
    printf "(2) p = %o, v = %o mod %o: reduction %o (squarefree, degree %o); #J(F_p) = %o = #E(F_p) #F(F_p) = %o * %o; J(F_p) = %o\n",
        p, r0, p, fbar, Degree(fbar), n, #Ebar, #Fbar, Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
end for;
assert GCD(nJ) eq 40;
assert noniso;
printf "(2) gcd of #J(F_p) = %o = #G; #E(F_p) != #F(F_p) at some p, so E and F_v are not isogenous at any member of the class\n", GCD(nJ);

// ---------------------------------------------------------------- (3) infinitely many points of D in the class
H := HyperellipticCurve(dpol);
s0 := 14;
assert Evaluate(dpol, 0) eq s0^2;
B := H![0, -s0, 1]; R0 := H![0, s0, 1];
EH, mH := EllipticCurve(H, B);
EM, mM := MinimalModel(EH);
iso := mH*mM;                       // H -> EM, B -> O
inv := Inverse(iso);
dp := DefiningPolynomials(inv);     // EM -> H (weights 1, 2, 1 on H)
R0e := iso(R0);
printf "(3) D : s^2 = %o is a genus-1 curve; its Jacobian EM = %o (minimal), with B = (0, -%o) -> O and R0 = (0, %o) -> %o\n", dpol, aInvariants(EM), s0, s0, R0e;
Pinf := EM![616, -15290];                  // fixed point of infinite order (Generators(EM) is not deterministic between runs)
assert Order(Pinf) eq 0;
printf "(3) EM(Q) has the point of infinite order P = %o (rank bounds %o)\n", Pinf, [RankBounds(EM)];
ords := [];
for p in primes do
    assert Integers()!Discriminant(EM) mod p ne 0;
    EMp := ChangeRing(EM, GF(p));
    Append(~ords, Order(EMp!Pinf));
    assert Valuation(R0e[1], p) ge 0 and Valuation(R0e[2], p) ge 0;
end for;
mstep := LCM(ords);
// v on EM near R0e: v = dp[1]/dp[3]; clear denominators and check the denominator is a unit at R0e
PZ := Parent(dp[1]);
cden := LCM([Denominator(cc) : cc in Coefficients(dp[1]) cat Coefficients(dp[3])]);
NX := cden*dp[1]; NZ := cden*dp[3];
vx := Evaluate(NX, [R0e[1], R0e[2], 1]); vz := Evaluate(NZ, [R0e[1], R0e[2], 1]);
assert vx eq 0 and vz ne 0;
for p in primes do assert Valuation(vz, p) eq 0; assert &and[Valuation(cc, p) ge 0 : cc in Coefficients(NX) cat Coefficients(NZ)]; end for;
printf "(3) P has orders %o mod %o; for R_k = R0 + %o k P the point R_k reduces to R0 mod every p, where v = NX/NZ has NZ a p-unit,
", ords, primes, mstep;
print "    so v(R_k) is p-integral and = 0 mod p: every R_k gives a member of the class with F_v of full rational 2-torsion";
R1 := R0e + mstep*Pinf;
v1 := Evaluate(NX, [R1[1], R1[2], 1])/Evaluate(NZ, [R1[1], R1[2], 1]);
assert IsSquare(Evaluate(dpol, v1)) and &and[Valuation(v1, p) ge 1 : p in primes];
printf "(3) sanity: v(R_1) has %o digits, lies on D and is 0 mod %o
", #Sprint(v1), primes;
// a second member over a different point of X_E (the points of D over w(0) all give the same curve C_0)
f0 := PolynomialRing(Q)![Evaluate(t, 0) : t in Coefficients(ff)];
g20 := G2Invariants(HyperellipticCurve(f0));
v2 := 0; k2 := 0;
for k in [1..8] do
    Rk := R0e + k*Pinf;
    if Rk eq EM!0 then continue; end if;
    im := inv(Rk);
    if im[3] eq 0 then continue; end if;
    vk := im[1]/im[3];
    assert IsSquare(Evaluate(dpol, vk));
    fk := PolynomialRing(Q)![Evaluate(t, vk) : t in Coefficients(ff)];
    if Degree(fk) ge 5 and Discriminant(fk) ne 0 and G2Invariants(HyperellipticCurve(fk)) ne g20 then v2 := vk; k2 := k; break; end if;
end for;
assert k2 gt 0;
printf "(3) the point R0 + %o P of D lies over v2 = %o, where C_v2 has G2-invariants different from those of C_0
", k2, v2;

// ---------------------------------------------------------------- (4) the marked group and members
function Member(v0)
    s := func< z | Evaluate(z, v0) >;
    R<xx> := PolynomialRing(Q);
    fv := R![s(t) : t in Coefficients(ff)];
    return fv, s;
end function;
function PullClass(fv, s, nn, hh, mm, X0, Y0)       // the class phi^*(X0,Y0) - phi^*(O) on y^2 = fv
    R := Parent(fv); C := HyperellipticCurve(fv); J := Jacobian(C);
    CRl := CoordinateRing(Ambient(C)); Y := CRl.2;
    nnv := R![s(t) : t in Coefficients(nn)]; hhv := R![s(t) : t in Coefficients(hh)]; mmv := R![s(t) : t in Coefficients(mm)];
    k := s(kk);
    N := nnv - X0*hhv; assert Degree(N) eq 3;
    RN := quo<R | N>; okinv, inv := IsInvertible(RN!(k*mmv)); assert okinv;
    vN := R!((RN!(Y0*hhv^2))*inv);
    assert (vN^2 - fv) mod N eq 0;
    D := Divisor(C, ideal<CRl | Hom3(N, CRl, 3), Y - Hom3(vN, CRl, 3)>);
    W := Divisor(C, ideal<CRl | Hom3(hhv, CRl, 3), Y>);
    assert Degree(D) eq 3 and Degree(W) eq 3;
    return J!(D - W);
end function;
function TwoTorsionClasses(fv)                     // the rational 2-torsion classes [q, 0], q | fv of degree 2
    R := Parent(fv); J := Jacobian(HyperellipticCurve(fv));
    fa := [g[1] : g in Factorization(fv)];
    qs := [g : g in fa | Degree(g) eq 2] cat [fa[i]*fa[j] : i, j in [1..#fa] | i lt j and Degree(fa[i]) eq 1 and Degree(fa[j]) eq 1];
    return {J| elt< J | q, 0 > : q in qs } join {J!0};
end function;
members := [Q| 0] cat [v2];
G2s := [];
for v0 in members do
    fv, s := Member(v0);
    assert Degree(fv) ge 5 and Discriminant(fv) ne 0;
    Append(~G2s, G2Invariants(HyperellipticCurve(fv)));
    Ev := EllipticCurve([s(t) : t in aInvariants(E0t)]);
    ok, isoE := IsIsomorphic(E, Ev); assert ok;
    Q0 := isoE(mTE(TE.1)); tv := s(tw);
    cE := PullClass(fv, s, nn1, h1, mm1, Q0[1]/tv, Q0[2]/tv);
    T2 := TwoTorsionClasses(fv);
    G := { k*cE + t : k in [0..Order(TE.1)-1], t in T2 };
    if #Sprint(v0) le 12 then
        invm := TorsionInvariants(fv);
        printf "(4) v = %o: phi1^*(generator of E(Q)_tors) has order %o, J(Q)[2] has %o elements, they generate %o classes; TorsionSubgroup = %o\n",
            v0, OrderUpTo(cE, 20), #T2, #G, invm;
        assert invm eq [2, 2, 10];
    else
        dv := Integers()!(Discriminant(fv)*LeadingCoefficient(fv)*Denominator(v0));
        gp := 0; np := 0;
        for p in PrimesInInterval(5, 400) do
            if dv mod p eq 0 then continue; end if;
            fp := PolynomialRing(GF(p))!fv;
            if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
            n, chi := FrobData(fp); gp := GCD(gp, n); np +:= 1;
        end for;
        printf "(4) v = v2: phi1^*(generator of E(Q)_tors) has order %o, J(Q)[2] has %o elements, they generate %o classes; gcd of #J(F_p) over %o good p < 400 = %o\n",
            OrderUpTo(cE, 20), #T2, #G, np, gp;
        assert gp eq 40;
    end if;
    assert OrderUpTo(cE, 20) eq Order(TE.1) and #G eq 40;
end for;
assert #Seqset(G2s) ge 2;
printf "(4) the members %o (points of D) give %o distinct G2-invariants: v -> M_2 is nonconstant on D, so the R_k give infinitely many isomorphism classes\n", members, #Seqset(G2s);
fc, s := Member(0);
R<xx> := PolynomialRing(Q);
Ccensus := HyperellipticCurve(R![-20,-12,66,-2,-45,3,9], R![0,1,1]);
assert IsIsomorphic(HyperellipticCurve(fc), Ccensus);
print "(4) the member v = 0 is Q-isomorphic to the census curve 2.2.10.b";
printf "time %o s\n", Cputime(T0);
print "FAMILY_2210_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E x F_v) and J(Q)_tors exactly [2,2,10]";
quit;
