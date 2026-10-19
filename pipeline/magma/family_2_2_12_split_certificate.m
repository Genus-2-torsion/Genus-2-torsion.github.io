/* family_2_2_12_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E x F_v) with J(Q)_tors exactly [2,2,12].
   Found 2026-10-19 (torsion_jac notes/g33ell_2026_10_19.md).

   Source.  (3,3)-gluings of the fixed curve E = 90c3 (torsion Z/12) with partners F that have full rational
   2-torsion, in the normal form of Broker-Howe-Lauter-Stevenhagen (arXiv:1403.6911, Alg. 5.4 and App. A,
   normalization b = a), exactly as in family_218_split_certificate.m:
       C : y^2 = tw h1 h2,   h1 = x^3 + 3ax + 2a,   h2 = 2dx^3 + 3cx^2 + 1,   d = (1 - 12ac)/(16a),
   with the degree-3 maps phi1 : (x,y) -> (nn1/h1, y mm1/h1^2) onto Y^2 = tw f1(X) and phi2 : (x,y) ->
   (nn2/h2, y mm2/h2^2) onto Y^2 = tw f2(X), certified by the identities h1^4 f1(nn1/h1) = h1 h2 mm1^2 and
   h2^4 f2(nn2/h2) = h1 h2 mm2^2.  w parametrizes the rational curve j(E0(a,c)) = j(E) (the partners
   anti-symplectically 3-congruent to E), v the rational cover "f2 has a rational root X(v)".  tw is the twist
   with E0^tw = E.  Then f2 = (X - X(v)) q(X), and F_v = F0^tw has full rational 2-torsion exactly when disc(q)
   is a square, i.e. when v lies on the genus-1 curve
       D : s^2 = v^4 - 4v^3 + 6v^2 - 652v + 649,
   whose Jacobian is EM : y^2 = x^3 + 9 (a twist of 36a1), with the point P = (3, 6) of infinite order and the
   3-torsion point T = (0, 3).  90c3 is the Z/12 curve of smallest conductor; for the next one, 210b5, this
   Jacobian is y^2 + y = x^3 + 12, also of rank 1.  The rational 3-torsion of E gives a rational 3-torsion point
   on F_v, so F_v(Q) contains Z/2 x Z/6, and by Mazur's theorem F_v(Q)_tors = Z/2 x Z/6 at every member.

   Argument.
   (1) Splitting.  The identities hold over Q(v); at every member of the residue class below, a, c, d, D1, D2,
       tw, kk, X(v), drt are finite and a, D1, D2, tw, kk, drt, disc f1, disc f2 are nonzero, so phi1, phi2 are
       morphisms over Q onto E and F_v, and J is isogenous over Q to E x F_v.
   (2) The marked group.  phi1 has prime degree 3 and C is not of genus 1, so phi1 does not factor through an
       isogeny and phi1^* : E -> J is injective; so phi1^*(E(Q)_tors) = Z/12 in J_v(Q).  phi2_* phi2^* = [3] on
       F_v, so phi2^* is injective on F_v(Q)[2] = (Z/2)^2 (v on D).  The kernel of phi1^* + phi2^* : E x F_v -> J
       lies in E[3] x F_v[3] (it is the graph of the gluing), so the two images meet trivially in 2-power order,
       and Z/12 + (Z/2)^2 = [2,2,12] embeds in J_v(Q).  At the members below this is also checked by Mumford
       arithmetic: the pullback of the point (-9, 49) of order 12, computed from divisors, has order 12, and with
       the 8 classes of J(Q)[2] (from the factors of f_v) it generates 48 classes.
   (3) Exactness on a residue class.  #J(F_p) = #E(F_p) #F(F_p) is divisible by 12 * 12, so the group structures
       are used.  For every 13- and 83-integral v = 101/5 (mod 13*83), C_v has the reductions of C_{101/5} at
       these primes, with J(F_13) = [2,6,12] (2-primary part [2,2,4], 3-primary part (Z/3)^2) and
       J(F_83) = [2,2,2,1008] (3-primary part Z/9, cyclic).  J_v(Q)_tors embeds in both, and the largest group
       that does is [2,2,12]: it is exactly [2,2,12].
   (4) Infinitely many members.  D -> EM sends B = (1, 0) to O; the member v = 101/5 lies over R0 = 2P.
       P has orders 21 and 7 mod 13 and 83, so R_k = R0 + 21 k P reduces to R0 mod both primes.  On EM,
       v = NX/NZ with integral NX, NZ and NZ(R0) a unit at the two primes, so v(R_k) is p-integral and
       = v(R0) = 101/5 mod p: each R_k gives a member of the class.  The R_k are distinct and v has degree 2
       on D, so there are infinitely many v.  The member over 3P (v = 8567/629) has other G2-invariants than
       C_{101/5}, so v -> M_2 is nonconstant on D with finite fibres: infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup at v = 101/5 and at v = 8567/629 is [2,2,12].  The census curve 2.2.12.b
       (J isogenous to 15a x 90c) is not a member: its partner would lie in the isogeny class 15a, and
       j(F_v) = j(E') has no rational solution v for the eight curves E' of 15a (see the note).
   Run (from the repository root):
       magma -b pipeline/magma/family_2_2_12_split_certificate.m > data/knowledge/sources/family_2_2_12_split_certificate.log */
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
// The largest finite abelian group that embeds in every group with the given invariants: for each prime l,
// the sorted l-primary exponents of a subgroup are bounded componentwise by those of the group.
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
E := EllipticCurve([1,-1,1,-122,1721]);
assert Conductor(E) eq 90 and Invariants(TorsionSubgroup(E)) eq [12];
P12 := E![-9, 49];
assert Order(P12) eq 12;
printf "E = 90c3 = %o, torsion %o, P12 = %o of order 12\n", aInvariants(E), Invariants(TorsionSubgroup(E)), P12;

// ---------------------------------------------------------------- (1) the family over Q(v)
Fv<v> := FunctionField(Q); w := v;
aW := (-50176/17523*w^6 - 15680/1947*w^5 - 5476/649*w^4 - 10445/2596*w^3 - 110889/127204*w^2 - 10935/127204*w - 19683/6232996)/(w^4 + 7019/18172*w^3 - 128223/127204*w^2 - 238383/890428*w - 98415/890428);
cW := (-5841/200704*w^2 + 26487/702464*w + 203391/9834496)/(w^4 + 15/8*w^3 + 13617/12544*w^2 + 1215/6272*w + 6561/614656);
wV := (v^3 - 3*v^2 + 480/7*v - 1600/7)/(v^3 - 3*v^2 - 160*v);
XV := (243/204800*v^18 - 2187/102400*v^17 - 243/8192*v^16 + 4131/2048*v^15 - 372519/40960*v^14 + 1358613/102400*v^13 - 249607413/204800*v^12 + 1034451/80*v^11 + 14714865/512*v^10 - 182563227/256*v^9 + 6181518321/2560*v^8 + 28499769/2*v^7 - 32705370*v^6 + 32659200*v^5 + 248832000*v^4)/(v^16 - 16*v^15 + 123*v^14 - 1412*v^13 + 11165*v^12 - 53088*v^11 + 253120*v^10 - 527098*v^9 - 1934757*v^8 + 29728268*v^7 - 2978687*v^6 - 443073246*v^5 + 175839947*v^4 + 9074102080*v^3 - 22295942400*v^2 - 92867584000*v + 106332160000);
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
dpol := V^4 - 4*V^3 + 6*V^2 - 652*V + 649;
ok, cc := IsSquare(dpol/(Pv!Evaluate(drep, V))); assert ok;
drt := drt/cc;
assert Discriminant(qq) eq Evaluate(dpol, v)*drt^2;
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

// ---------------------------------------------------------------- (3) exactness on the residue class v = 101/5 mod [13, 83]
primes := [13, 83]; vres := 101/5;
finite := [aV, cV, dV, D1, D2, tw, kk, XV, drt];
nonzero := [aV, D1, D2, tw, kk, drt, Discriminant(f1), Discriminant(f2)];
invs := [];
for p in primes do
    r0 := Integers()!(GF(p)!vres);
    for z in finite do Nz, Dz := IntNumDen(z); assert (GF(p)!Evaluate(Dz, r0)) ne 0; end for;
    for z in nonzero do Nz, Dz := IntNumDen(z); assert (GF(p)!Evaluate(Nz, r0)) ne 0; end for;
    fbar := ReductionClass(ff, p, r0);
    n, chi := FrobData(fbar);
    Ebar := ChangeRing(E, GF(p));
    Fbar := EllipticCurve([GF(p)!Evaluate(Numerator(t), r0)/Evaluate(Denominator(t), r0) : t in aInvariants(Ft)]);
    assert n eq #Ebar*#Fbar;
    A := AbelianGroup(Jacobian(HyperellipticCurve(fbar)));
    assert #A eq n;
    Append(~invs, Invariants(A));
    printf "(3) p = %o, v = %o mod %o: reduction %o (squarefree, degree %o); #J(F_p) = %o = #E(F_p) #F(F_p) = %o * %o; J(F_p) = %o\n",
        p, r0, p, fbar, Degree(fbar), n, #Ebar, #Fbar, Invariants(A);
end for;
bnd := StructureBound(invs);
printf "(3) the largest group embedding in all J(F_p) has invariants %o\n", bnd;
assert bnd eq [2, 2, 12];
printf "(3) every p-integral v = %o mod %o (p in %o) has J(Q)_tors = G = [2,2,12]\n", vres, &*primes, primes;

// ---------------------------------------------------------------- (4) infinitely many points of D in the class
H := HyperellipticCurve(dpol);
B := H![1, 0, 1];
EH, mH := EllipticCurve(H, B);
EM, mM := MinimalModel(EH);
iso := mH*mM;                       // H -> EM, B -> O
inv := Inverse(iso);
dp := DefiningPolynomials(inv);     // EM -> H (weights 1, 2, 1 on H)
assert aInvariants(EM) eq [0, 0, 0, 0, 9];
Pinf := EM![3, 6]; T3 := EM![0, 3];        // fixed points (Generators(EM) is not deterministic between runs)
assert Order(Pinf) eq 0 and Order(T3) eq 3;
R0e := 2*Pinf;
im0 := inv(R0e); assert im0[1]/im0[3] eq vres;
printf "(4) D : s^2 = dpol(v) is a genus-1 curve; its Jacobian EM = %o (minimal), with B = (1, 0) -> O; P = %o has infinite order (rank bounds %o), T = %o has order 3; the member v = %o lies over R0 = 2P = %o\n",
    aInvariants(EM), Pinf, [RankBounds(EM)], T3, vres, R0e;
ords := [];
for p in primes do
    assert Integers()!Discriminant(EM) mod p ne 0;
    EMp := ChangeRing(EM, GF(p));
    Append(~ords, Order(EMp!Pinf));
    assert Valuation(R0e[1], p) ge 0 and Valuation(R0e[2], p) ge 0;
end for;
mstep := LCM(ords);
PZ := Parent(dp[1]);
cden := LCM([Denominator(c) : c in Coefficients(dp[1]) cat Coefficients(dp[3])]);
NX := cden*dp[1]; NZ := cden*dp[3];
vx := Evaluate(NX, [R0e[1], R0e[2], 1]); vz := Evaluate(NZ, [R0e[1], R0e[2], 1]);
assert vz ne 0 and vx/vz eq vres;
for p in primes do assert Valuation(vz, p) eq 0; assert &and[Valuation(c, p) ge 0 : c in Coefficients(NX) cat Coefficients(NZ)]; end for;
printf "(4) P has orders %o mod %o; for R_k = R0 + %o k P the point R_k reduces to R0 mod every p, where v = NX/NZ has NZ a p-unit,\n", ords, primes, mstep;
print "    so v(R_k) is p-integral and = v(R0) mod p: every R_k gives a member of the class with F_v of full rational 2-torsion";
R1 := R0e + mstep*Pinf;
v1 := Evaluate(NX, [R1[1], R1[2], 1])/Evaluate(NZ, [R1[1], R1[2], 1]);
assert IsSquare(Evaluate(dpol, v1)) and &and[Valuation(v1 - vres, p) ge 1 : p in primes];
printf "(4) sanity: v(R_1) has %o digits, lies on D and is = %o mod %o\n", #Sprint(v1), vres, primes;
f0 := PolynomialRing(Q)![Evaluate(t, vres) : t in Coefficients(ff)];
g20 := G2Invariants(HyperellipticCurve(f0));
v2 := 0; k2 := 0;
for k in [2..8] do
    Rk := k*Pinf;
    im := inv(Rk);
    if im[3] eq 0 then continue; end if;
    vk := im[1]/im[3];
    assert IsSquare(Evaluate(dpol, vk));
    fk := PolynomialRing(Q)![Evaluate(t, vk) : t in Coefficients(ff)];
    if Degree(fk) ge 5 and Discriminant(fk) ne 0 and G2Invariants(HyperellipticCurve(fk)) ne g20 then v2 := vk; k2 := k; break; end if;
end for;
assert k2 gt 0;
printf "(4) the point %o P of D lies over v2 = %o, where C_v2 has G2-invariants different from those of C_{%o}\n", k2, v2, vres;

// ---------------------------------------------------------------- (2), (4) the marked group and members
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
members := [Q| vres, v2];
G2s := [];
for v0 in members do
    fv, s := Member(v0);
    assert Degree(fv) ge 5 and Discriminant(fv) ne 0;
    Append(~G2s, G2Invariants(HyperellipticCurve(fv)));
    Ev := EllipticCurve([s(t) : t in aInvariants(E0t)]);
    ok, isoE := IsIsomorphic(E, Ev); assert ok;
    Q0 := isoE(P12); tv := s(tw);
    cE := PullClass(fv, s, nn1, h1, mm1, Q0[1]/tv, Q0[2]/tv);
    T2 := TwoTorsionClasses(fv);
    G := { k*cE + t : k in [0..11], t in T2 };
    printf "(2) v = %o: phi1^*(P12) has order %o, J(Q)[2] has %o elements, they generate %o classes\n", v0, OrderUpTo(cE, 30), #T2, #G;
    assert OrderUpTo(cE, 30) eq 12 and #T2 eq 8 and #G eq 48;
end for;
assert #Seqset(G2s) eq 2;
printf "(4) the members %o (points of D) have distinct G2-invariants: v -> M_2 is nonconstant on D, so the R_k give infinitely many isomorphism classes\n", members;
for v0 in [Q| vres, v2] do
    fv, s := Member(v0);
    invm := TorsionInvariants(fv);
    printf "(5) TorsionSubgroup at v = %o: %o\n", v0, invm;
    assert invm eq [2, 2, 12];
end for;
printf "time %o s\n", Cputime(T0);
print "FAMILY_2_2_12_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E x F_v) and J(Q)_tors exactly [2,2,12]";
quit;
