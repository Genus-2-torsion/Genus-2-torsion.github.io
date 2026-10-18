/* family_218_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E x F_v) with J(Q)_tors exactly [2,18].
   Found 2026-10-18 (torsion_jac notes/s18_2026_10_18.md).

   Source.  (3,3)-gluings of the fixed curve E = 714i1 (torsion Z/9) with partners F that have full rational
   2-torsion, in the normal form of Broker-Howe-Lauter-Stevenhagen (arXiv:1403.6911, Alg. 5.4 and App. A,
   normalization b = a):
       C : y^2 = tw h1 h2,   h1 = x^3 + 3ax + 2a,   h2 = 2dx^3 + 3cx^2 + 1,   d = (1 - 12ac)/(16a),
   with the degree-3 maps phi1 : (x,y) -> (nn1/h1, y mm1/h1^2) onto Y^2 = tw f1(X) and phi2 : (x,y) ->
   (nn2/h2, y mm2/h2^2) onto Y^2 = tw f2(X), certified by the identities h1^4 f1(nn1/h1) = h1 h2 mm1^2 and
   h2^4 f2(nn2/h2) = h1 h2 mm2^2.  w parametrizes the rational curve j(E0(a,c)) = j(E) (the partners
   anti-symplectically 3-congruent to E), v the rational cover "f2 has a rational root X(v)".  tw is the twist
   with E0^tw = E.  Then f2 = (X - X(v)) q(X), and F_v = F0^tw has full rational 2-torsion exactly when disc(q)
   is a square, i.e. when v lies on the genus-1 curve
       D : s^2 = -5183v^4 - 60435076v^3 - 234870845178v^2 - 304157704212868v + 304392635498305,
   whose Jacobian is EM : y^2 = x^3 + 289 (a twist of 36a1, as Gamma(3) cap Gamma(2) = Gamma(6) predicts), with
   the point P = (-4, -15) of infinite order and the 3-torsion point T = (0, 17).  For the Z/9 curves 54b3,
   1482l1, 1554 (Kubert parameter -2), 1890, 17670, 92910 this Jacobian has analytic rank 0; 714i1 is the Z/9 curve of
   smallest conductor for which it has positive rank.  The rational 3-torsion of E gives a rational 3-torsion
   point on F_v, so F_v(Q) contains Z/2 x Z/6.

   Argument.
   (1) Splitting.  The identities hold over Q(v); at every member of the residue class below, a, c, d, D1, D2,
       tw, kk, X(v), drt are finite and a, D1, D2, tw, kk, drt, disc f1, disc f2 are nonzero, so phi1, phi2 are
       morphisms over Q onto E and F_v, and J is isogenous over Q to E x F_v.
   (2) The marked group.  phi1 has prime degree 3 and C is not of genus 1, so phi1 does not factor through an
       isogeny and phi1^* : E -> J is injective; so phi1^*(E(Q)_tors) = Z/9 in J_v(Q).  phi2_* phi2^* = [3] on
       F_v, so phi2^* is injective on F_v(Q)[2] = (Z/2)^2 (v on D); the orders are coprime, so [2,18] embeds in
       J_v(Q).  At the members below this is also checked by Mumford arithmetic (the pullback of the point
       (24, -180) of order 9, computed from divisors, and the classes of the quadratic factors of f_v generate
       36 classes).
   (3) Exactness on a residue class.  #J(F_p) = #E(F_p) #F(F_p) is divisible by 9 * 12, so the group
       structures are used.  For every 109-, 251- and 359-integral v = -3023 (mod 109*251*359), C_v has the
       reductions of C_{-3023} at these primes, with J(F_109) = [2,2,54,54] (2-primary exponent 2),
       J(F_251) = [2,31320] (3-primary part Z/27, cyclic) and J(F_359) = [2,2,180,180] (3-primary exponent 9,
       and no 29-torsion).  J_v(Q)_tors embeds in all three, and the largest group that does is [2,18]: it is
       exactly [2,18].
   (4) Infinitely many members.  D -> EM sends B = (1, 0) to O; the member v = -3023 lies over R0 = P + T.
       P has orders 6, 18, 9 mod 109, 251, 359, so R_k = R0 + 18 k P reduces to R0 mod these primes.  On EM,
       v = NX/NZ with integral NX, NZ and NZ(R0) a unit at the three primes, so v(R_k) is p-integral and
       = v(R0) = -3023 mod p: each R_k gives a member of the class.  The R_k are distinct and v has degree 2
       on D, so there are infinitely many v.  The member over 2P has other G2-invariants than C_{-3023}, so
       v -> M_2 is nonconstant on D with finite fibres: infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup at v = -3023 and at v = -51395/13 (over P) is [2,18].  No census curve is a
       member: every member has the isogeny factor 714i1, so 714 divides its conductor, while the census curve
       2.18.b has conductor 64980 (its Jacobian is isogenous to 114a x 570k; see the note).
   Run (from the repository root):
       magma -b pipeline/magma/family_218_split_certificate.m > data/knowledge/sources/family_218_split_certificate.log */
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
E := EllipticCurve([1,0,0,108,11664]);
assert Conductor(E) eq 714 and Invariants(TorsionSubgroup(E)) eq [9];
P9 := E![24, -180];
assert Order(P9) eq 9;
printf "E = 714i1 = %o, torsion %o, P9 = %o of order 9\n", aInvariants(E), Invariants(TorsionSubgroup(E)), P9;

// ---------------------------------------------------------------- (1) the family over Q(v)
Fv<v> := FunctionField(Q); w := v;
aW := (-5637866810112/13032150911*w^6 + 2591*w^5 - 109377500284479167/16913600430336*w^4 + 7388828400737427230506600618847/858209638551186373167218688*w^3 - 18576341852732664313035077123843855807/2880701087216333744690137325961216*w^2 + 74736338099248096670494312668433557213030431/29008432372882590718663852351580694183936*w - 375907077053230448259127180108932930751356677895551/876337866986310692436834656208245446672184967168)/(w^5 - 50655970591057/16913600430336*w^4 + 2566027356522031669928377249/858209638551186373167218688*w^3 - 11471745670501305941733921521124673/11522804348865334978760549303844864*w^2);
cW := (-1/10368*w^4 - 13032150911/135308803442688*w^3 + 880264950033565947380543/4577118072272993990225166336*w^2)/(w^4 - 33766303010401/8456800215168*w^3 + 1710684960960340235644961473/286069879517062124389072896*w^2 - 5734766164821316553160831692308321/1440350543608166872345068662980608*w + 28844592087706714268812934260298555466241/29008432372882590718663852351580694183936);
wV := (-39096452733/10069921*v^2 - 303896727093609/10069921*v - 590597974153284474/10069921)/(v^3 + 7773*v^2 + 15100995*v - 15108769);
XV := (-38078453403648/5183*v^14 - 2071886728145891328/5183*v^13 - 52341605981535494292480/5183*v^12 - 813734247481220344573808640/5183*v^11 - 8697563375516691929925772391424/5183*v^10 - 67610716675260301871823717630535680/5183*v^9 - 394184227156361580006140457925492549632/5183*v^8 - 1751003480559514242220725540988122513997824/5183*v^7 - 5955155856535719520173402663995115340402710528/5183*v^6 - 15431248563554686804251972638428456956702585741312/5183*v^5 - 29989571657479318012997813037863716476617921976038400/5183*v^4 - 42387288151605545194117992087796088505024293846028349440/5183*v^3 - 41188070706943271240642571290167771116401590295992418509824/5183*v^2 - 24629212200274983553488279958274572984494083147040889202706432/5183*v - 6837705696887869480885887532488417765197103478002029505156194304/5183)/(v^18 + 342447289/5183*v^17 + 10647567531357/5183*v^16 + 206914017268120112/5183*v^15 + 2814686150995328406808/5183*v^14 + 28443136493109347685154332/5183*v^13 + 221097239256814800487409011784/5183*v^12 + 1350376535187155479390961935586800/5183*v^11 + 6560577853360700754838120312021811454/5183*v^10 + 25498724633632604368408556313637579454086/5183*v^9 + 79283416834062913276383696964526724113350682/5183*v^8 + 196090557354738593690887452138995605297275210576/5183*v^7 + 5366958896577615387307394356459670865335601162560/73*v^6 + 569582818874776925172756224616692869794005995780945340/5183*v^5 + 632414243035168249080034663633710086512425358356645804288/5183*v^4 + 491450702705764056969167280722665067791290014028019046119696/5183*v^3 + 238601344989720042119126297331730186564689799600550401738668011/5183*v^2 + 54429744566730015899320113295130656895839571402139304591399775337/5183*v - 54668837995406648809790758052502247118587917873639628620819944895/5183);
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
dpol := -5183*V^4 - 60435076*V^3 - 234870845178*V^2 - 304157704212868*V + 304392635498305;
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

// ---------------------------------------------------------------- (3) exactness on the residue class v = -3023 mod [109, 251, 359]
primes := [109, 251, 359]; vres := -3023;
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
assert bnd eq [2, 18];
printf "(3) every p-integral v = %o mod %o (p in %o) has J(Q)_tors = G = [2,18]\n", vres, &*primes, primes;

// ---------------------------------------------------------------- (4) infinitely many points of D in the class
H := HyperellipticCurve(dpol);
B := H![1, 0, 1];
EH, mH := EllipticCurve(H, B);
EM, mM := MinimalModel(EH);
iso := mH*mM;                       // H -> EM, B -> O
inv := Inverse(iso);
dp := DefiningPolynomials(inv);     // EM -> H (weights 1, 2, 1 on H)
assert aInvariants(EM) eq [0, 0, 0, 0, 289];
Pinf := EM![-4, -15]; T3 := EM![0, 17];        // fixed points (Generators(EM) is not deterministic between runs)
assert Order(Pinf) eq 0 and Order(T3) eq 3;
R0e := Pinf + T3;
im0 := inv(R0e); assert im0[1]/im0[3] eq vres;
printf "(4) D : s^2 = dpol(v) is a genus-1 curve; its Jacobian EM = %o (minimal), with B = (1, 0) -> O; P = %o has infinite order (rank bounds %o), T = %o has order 3; the member v = %o lies over R0 = P + T = %o\n",
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
    Q0 := isoE(P9); tv := s(tw);
    cE := PullClass(fv, s, nn1, h1, mm1, Q0[1]/tv, Q0[2]/tv);
    T2 := TwoTorsionClasses(fv);
    G := { k*cE + t : k in [0..8], t in T2 };
    printf "(2) v = %o: phi1^*(P9) has order %o, J(Q)[2] has %o elements, they generate %o classes\n", v0, OrderUpTo(cE, 20), #T2, #G;
    assert OrderUpTo(cE, 20) eq 9 and #T2 eq 4 and #G eq 36;
end for;
assert #Seqset(G2s) eq 2;
printf "(4) the members %o (points of D) have distinct G2-invariants: v -> M_2 is nonconstant on D, so the R_k give infinitely many isomorphism classes\n", members;
for v0 in [Q| vres, -51395/13] do
    fv, s := Member(v0);
    invm := TorsionInvariants(fv);
    printf "(5) TorsionSubgroup at v = %o: %o\n", v0, invm;
    assert invm eq [2, 18];
end for;
printf "time %o s\n", Cputime(T0);
print "FAMILY_2_18_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E x F_v) and J(Q)_tors exactly [2,18]";
quit;
