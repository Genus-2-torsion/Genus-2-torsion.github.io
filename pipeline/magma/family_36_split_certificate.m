/* family_36_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E x F_w) with J(Q)_tors exactly [36].
   Found 2026-10-19 (torsion_jac notes/g33rat_27_2026_10_19.md).

   Source.  The family of family_12_split_certificate.m: (3,3)-gluings, in the normal form of Broker-Howe-Lauter-
   Stevenhagen (arXiv:1403.6911, Alg. 5.4, App. A, b = a), of E = 90c3 (torsion Z/12, P12 = (-9,49), P3 = 4 P12)
   with all partners F_w anti-symplectically 3-congruent to E, over the rational curve X_E^-(3) with parameter w;
   C_w : y^2 = f_w(x) with explicit degree-3 maps phi1 onto E and phi2 onto F_w, J ~ E x F_w over Q.  F_w has
   the rational 3-torsion point T = psi(P3).  Where F_w also has a rational point Q of order 9 with 3Q = T,
   phi1^*(P12) (order 12) and phi2^*(Q) (order 9) generate Z/36 in J(Q).  Both pullbacks are injective (prime
   degree), the 3-parts overlap in phi2^*(T) = -phi1^*(P3), and the 2-part comes from E.  This locus has genus 1
   and positive rank.

   The locus.  Z := {(w, Q) : Q in F_w, 3Q = T} is the untwisted form of the modular curve of Gamma(3) cap
   Gamma1(9) over the w-line.  It is explicit: x(3Q) = x(T) is irreducible of degree 9 over Q(w), and the script
   computes genus 1.  Z(Q) = {w : F_w(Q) has a point of order 9 over T}.  Translation by the rational point T acts
   freely on Z, and Z/<T> is the mu_3-torsor {phi(Q)}, phi : F_w -> F_w/<T>.  Its Kummer class is the descent
   class of T for the dual isogeny, 1/A3 mod cubes, from the model y^2 + A1 xy + A3 y = x^3 with T = (0,0).  The
   script checks A3 = 25 g(w) mod cubes, so Z/<T> = C0 : z^3 = 25 g(w), g = 890428w^3 + 1234359w^2 + 336798w + 98415,
   a smooth plane cubic of genus 1 with Jacobian y^2 + y = x^3 - 1 of rank 1.  The base points R0 = (0, 135)
   (w = 0) and S = (-81/119, 2700/17) (w = -81/119, the census curve 36.b) lift to Z(Q): F_w has a point of order 9
   over T.  q : Z -> C0 is an étale cyclic triple cover of genus-1 curves (quotient by a free translation).  With
   origins over R0 it is an isogeny, so q(Z(Q)) = R0 + H for a subgroup H of C0(Q); S lifting gives S - R0 in H.

   Argument.
   (1) Splitting.  As in family_12_split_certificate.m, on the class below.
   (2) The marked group.  Every R_k = R0 + 4k(S - R0) lies in q(Z(Q)) by the coset argument above, so F_{w(R_k)}
       has a rational point Q of order 9 with 3Q = T, and phi1^*(E(Q)_tors) + <phi2^*(Q)> = Z/36 in J(Q).
       Checked at the members k = 0, 1: x(Q) is a rational root of Z at w(R_k), and 3Q = T.  At w = 0 and
       w = -81/119 Mumford arithmetic gives the subgroup [36], and TorsionSubgroup is [36].
   (3) Exactness.  C0 has good reduction at 11 and 37, and S - R0 has order 4 in C0(F_11) and C0(F_37).  So every
       R_k reduces to R0 and w(R_k) = 0 (mod 11 * 37), p-integral.  Every 11- and 37-integral w = 0 (mod 407) has
       the reductions of C_0 at 11 and 37, where J(F_11) = Z/108 (cyclic) and J(F_37) = [36,36] (exponent 36).
       The largest group embedding in both is [36], so J(Q)_tors = [36] on the class.
   (4) Infinitely many.  The R_k are distinct (S - R0 has infinite order).  The G2-invariants at k = 0, 1 differ,
       so R -> M_2 is nonconstant on C0 with finite fibres.
   (5) Sanity and control: TorsionSubgroup [36] at w = 0 and at w = -81/119, which is Q-isomorphic to the census
       curve 36.b.
   Run (from the repository root):
       magma -b pipeline/magma/family_36_split_certificate.m > data/knowledge/sources/family_36_split_certificate.log */
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
    ls := &meet[ Seqset(&cat([[Integers()| ]] cat [ PrimeDivisors(n) : n in inv | n gt 1 ])) : inv in invlist ];
    out := [Integers()| ];
    for l in ls do
        es := [ Reverse(Sort([ Valuation(n, l) : n in inv | n mod l eq 0 ])) : inv in invlist ];
        r := Min([ #e : e in es ]);
        out cat:= [ l^Min([ e[i] : e in es ]) : i in [1..r] ];
    end for;
    return Invariants(AbelianGroup(out));
end function;
// The subgroup of J(Q) generated by classes of finite order (as a list of its elements; small groups only).
function GeneratedSet(J, gens, bnd)
    S := [J!0];
    for g in gens do
        o := OrderUpTo(g, bnd); assert o gt 0;
        new := S;
        for s in S do for k in [1..o-1] do
            y := s + k*g;
            if not exists{z : z in new | z eq y} then Append(~new, y); end if;
        end for; end for;
        S := new;
    end for;
    return S;
end function;
// A finite abelian group is determined by the numbers #A[m]; compare a list of elements with invariants G.
function SameGroup(S, G)
    A := AbelianGroup(G);
    if #S ne #A then return false; end if;
    return &and[ #[z : z in S | m*z eq Parent(z)!0] eq #[a : a in A | m*a eq A!0] : m in Divisors(#A) ];
end function;

// ---------------------------------------------------------------- (0) the fixed curve E
G := [36];
E := EllipticCurve([1,-1,1,-122,1721]);
TE := [12];
assert Conductor(E) eq 90 and Invariants(TorsionSubgroup(E)) eq TE;
gensE := [E![-9, 49]];
assert [Order(P) : P in gensE] eq [12];
if #gensE gt 0 then
    assert #{ &+[k[i]*gensE[i] : i in [1..#gensE]] : k in CartesianProduct([[0..Order(P)-1] : P in gensE]) } eq &*TE;   // the generators generate E(Q)_tors
end if;
printf "E = 90c3 = %o, conductor %o, torsion %o; generators used for the pullbacks: %o (2-torsion classes come from f)\n", aInvariants(E), Conductor(E), TE, gensE;

// ---------------------------------------------------------------- (1) the family over Q(w)
Fw<w> := FunctionField(Q);
aW := (-50176/17523*w^6 - 15680/1947*w^5 - 5476/649*w^4 - 10445/2596*w^3 - 110889/127204*w^2 - 10935/127204*w - 19683/6232996)/(w^4 + 7019/18172*w^3 - 128223/127204*w^2 - 238383/890428*w - 98415/890428);
cW := (-5841/200704*w^2 + 26487/702464*w + 203391/9834496)/(w^4 + 15/8*w^3 + 13617/12544*w^2 + 1215/6272*w + 6561/614656);
aP := aW; cP := cW;
h1, h2, f1, f2, nn1, mm1, nn2, mm2, dP, D1, D2 := BHLSParts(aP, cP);
assert jInvariant(EllipticCurve(f1)) eq jInvariant(E);
print "(1) j(E0(a(w), c(w))) = j(E) identically in w";
assert CoverIdentity(f1, nn1, h1, mm1, h1, h2) and CoverIdentity(f2, nn2, h2, mm2, h1, h2);
print "(1) cover identities h1^4 f1(nn1/h1) = h1 h2 mm1^2 and h2^4 f2(nn2/h2) = h1 h2 mm2^2 hold over Q(w)";
tw := TwistTo(f1, E);
E0t := TwistedCurve(f1, tw); Ft := TwistedCurve(f2, tw);
assert cInvariants(E0t) eq cInvariants(E);
print "(1) E0^tw has the c-invariants of E; F_w = F0^tw";
rep, rt := SqfRep(tw);
PX<x> := Parent(h1);
ff := rep*h1*h2; den := LCM([Denominator(t) : t in Coefficients(ff)]); ff := ff*den^2;
Pv<V> := PolynomialRing(Q);
gc := GCD([Pv!Evaluate(Numerator(t), V) : t in Coefficients(ff)]);
sq := &*([Pv| 1] cat [fa[1]^(fa[2] div 2) : fa in Factorization(gc)]);
sqv := Evaluate(sq, w); ff := ff/sqv^2;
assert &and[Denominator(t) eq 1 : t in Coefficients(ff)];
kk := rt*sqv/den;
assert ff*kk^2 eq tw*h1*h2;
printf "(1) C_w : y^2 = f_w(x), coefficients of degrees %o in w; y' = kk y maps C_w to y'^2 = tw h1 h2\n", [Degree(Numerator(t)) : t in Coefficients(ff)];
printf "(1) phi1 : (x,y) -> (tw nn1/h1, tw kk y mm1/h1^2) onto E0^tw ~ E, phi2 : (x,y) -> (tw nn2/h2, tw kk y mm2/h2^2) onto F\n";
fac := Factorization(ff);
assert Degree(ff) eq 6 and &and[g[2] eq 1 : g in fac];
lins := [g[1] : g in fac | Degree(g[1]) eq 1];
quads := [g[1] : g in fac | Degree(g[1]) eq 2] cat [lins[i]*lins[j] : j in [i+1..#lins], i in [1..#lins]];
printf "(1) factor degrees of f_w over Q(w): %o; %o quadratic divisors, so #J(Q(w))[2] = %o\n", [Degree(g[1]) : g in fac], #quads, 1 + #quads;
assert 1 + #quads eq 2^1;


// ---------------------------------------------------------------- (G) the locus where F_w has a rational point of order 9
PXw<X> := PolynomialRing(Fw);
a1F, a2F, a3F, a4F, a6F := Explode(aInvariants(Ft));
fa3 := Factorization(PXw!DivisionPolynomial(Ft, 3));
lin := [f[1] : f in fa3 | Degree(f[1]) eq 1]; assert #lin eq 1;
xT := -Coefficient(lin[1], 0)/Coefficient(lin[1], 1);
okT, yT := IsSquare(xT^3 + a2F*xT^2 + a4F*xT + a6F); assert okT;
T := Ft![xT, yT]; assert T ne Ft!0 and 3*T eq Ft!0;
function FlexForm(K, ai, x0, y0)       // coordinates x = X + x0, y = Y + s X + y0 putting the flex (x0,y0) at (0,0)
    b1, b2, b3, b4, b6 := Explode(ai);
    s := (3*x0^2 + 2*b2*x0 + b4 - b1*y0)/(2*y0 + b1*x0 + b3);
    PR<XX, YY> := PolynomialRing(K, 2);
    xx := XX + x0; yy := YY + s*XX + y0;
    Gp := yy^2 + b1*xx*yy + b3*yy - (xx^3 + b2*xx^2 + b4*xx + b6);
    A1 := MonomialCoefficient(Gp, XX*YY); A3 := MonomialCoefficient(Gp, YY);
    assert Gp eq YY^2 + A1*XX*YY + A3*YY - XX^3;
    return A1, A3, s;
end function;
function IsCubeFn(z)        // z in Q(w)^* a cube?
    ok1, c1 := IsPower(LeadingCoefficient(Numerator(z))/LeadingCoefficient(Denominator(z)), 3);
    return ok1 and &and[f[2] mod 3 eq 0 : f in Factorization(Numerator(z)) cat Factorization(Denominator(z))];
end function;
A1F, A3F := FlexForm(Fw, aInvariants(Ft), xT, yT);
Pw<W> := PolynomialRing(Q);
g := 890428*W^3 + 1234359*W^2 + 336798*W + 98415;
assert IsIrreducible(g);
assert IsCubeFn(A3F/(25*Evaluate(g, w)));
print "(G1) F_w has the rational 3-torsion point T; in the model y^2 + A1 xy + A3 y = x^3 with T = (0,0), A3(w) = 25 g(w) * (cube in Q(w)), g = 890428w^3 + 1234359w^2 + 336798w + 98415 (irreducible)";
// the untwisted locus Z = {(w, x(Q)) : Q in F_w, 3Q = +-T}
b2F := 4*a2F; b4F := 2*a4F; b6F := 4*a6F; b8F := 4*a2F*a6F - a4F^2;
psi3 := 3*X^4 + b2F*X^3 + 3*b4F*X^2 + 3*b6F*X + b8F;
psi4f := 2*X^6 + b2F*X^5 + 5*b4F*X^4 + 10*b6F*X^3 + 10*b8F*X^2 + (b2F*b8F - b4F*b6F)*X + (b4F*b8F - b6F^2);
phi3 := X*psi3^2 - psi4f*(4*X^3 + b2F*X^2 + 2*b4F*X + b6F);
gZ := phi3 - xT*psi3^2;
fZ := Factorization(gZ);
assert Degree(gZ) eq 9 and #fZ eq 1 and fZ[1][2] eq 1;
assert Genus(FunctionField(gZ/LeadingCoefficient(gZ))) eq 1;
print "(G2) Z : x(3Q) = x(T) is irreducible of degree 9 over Q(w) and has genus 1 (the modular curve of Gamma(3) cap Gamma1(9) over the w-line)";
P2<WW, ZZ, UU> := ProjectivePlane(Q);
gh := 890428*WW^3 + 1234359*WW^2*UU + 336798*WW*UU^2 + 98415*UU^3;
C0 := Curve(P2, ZZ^3 - 25*gh);
assert Genus(C0) eq 1;
R0 := [Q| 0, 135]; S := [Q| -81/119, 2700/17];
for Pt in [R0, S] do assert Pt[2]^3 eq 25*Evaluate(g, Pt[1]); end for;
A0, mA0 := EllipticCurve(C0, C0![R0[1], R0[2], 1]);
A0m, isoA := MinimalModel(A0);
PS := isoA(mA0(C0![S[1], S[2], 1]));
assert Order(PS) eq 0;
printf "(G3) C0 : z^3 = 25 g(w) has genus 1 and contains R0 = (0, 135) and S = (-81/119, 2700/17); Jac(C0) = %o, S - R0 = %o of infinite order\n", aInvariants(A0m), PS;
function Point9(v0)    // a rational Q on F_v0 with 3Q = T, from a rational root of gZ at v0
    Fv := EllipticCurve([Evaluate(t, v0) : t in aInvariants(Ft)]);
    Tv := Fv![Evaluate(xT, v0), Evaluate(yT, v0)];
    gv := PolynomialRing(Q)![Evaluate(cf, v0) : cf in Coefficients(gZ)];
    for r in Roots(gv) do
        av := aInvariants(Fv); assert av[1] eq 0 and av[3] eq 0;
        ok, yr := IsSquare(r[1]^3 + av[2]*r[1]^2 + av[4]*r[1] + av[5]);
        if ok then
            for Qv in [Fv![r[1], yr], Fv![r[1], -yr]] do if 3*Qv eq Tv then return true, Qv, Fv, Tv; end if; end for;
        end if;
    end for;
    return false, Fv!0, Fv, Tv;
end function;
for Pt in [R0, S] do
    ok, Qv, Fv := Point9(Pt[1]); assert ok and Order(Qv) eq 9;
    printf "(G3) base point w = %o: F_w has the rational point %o of order 9 with 3Q = T\n", Pt[1], Qv;
end for;
function Third(K, c, g, P1, P2)
    w1, z1 := Explode(P1); w2, z2 := Explode(P2);
    gK := PolynomialRing(K)!g; gp := Derivative(gK);
    if w1 ne w2 then lam := (z2 - z1)/(w2 - w1);
    else assert z1 eq z2 and z1 ne 0; lam := c*Evaluate(gp, w1)/(3*z1^2); end if;
    XK := PolynomialRing(K).1;
    cub := (lam*(XK - w1) + z1)^3 - c*gK;
    assert Degree(cub) eq 3;
    w3 := -Coefficient(cub, 2)/Coefficient(cub, 3) - w1 - w2;
    return [K| w3, lam*(w3 - w1) + z1];
end function;
function Add(K, c, g, O, P1, P2) return Third(K, c, g, O, Third(K, c, g, P1, P2)); end function;
function Neg(K, c, g, O, P1) return Third(K, c, g, Third(K, c, g, O, O), P1); end function;
function Mul(K, c, g, O, k, P1)
    if k lt 0 then return Mul(K, c, g, O, -k, Neg(K, c, g, O, P1)); end if;
    Rk := O; B := P1;
    while k gt 0 do
        if IsOdd(k) then Rk := (Rk eq O) select B else Add(K, c, g, O, Rk, B); end if;
        k div:= 2; if k gt 0 then B := Add(K, c, g, O, B, B); end if;
    end while;
    return Rk;
end function;
// the members: R_k = R0 + 4k (S - R0) in the group law of C0 with origin R0
for p in [11, 37] do
    Fp := GF(p);
    P2p<Wp, Zp, Up> := ProjectivePlane(Fp);
    C0p := Curve(P2p, Zp^3 - 25*(890428*Wp^3 + 1234359*Wp^2*Up + 336798*Wp*Up^2 + 98415*Up^3));
    assert not IsSingular(C0p);
    Ap, mAp := EllipticCurve(C0p, C0p![Fp!R0[1], Fp!R0[2], 1]);
    assert Order(mAp(C0p![Fp!S[1], Fp!S[2], 1])) eq 4;
end for;
print "(G4) C0 has good reduction at 11 and 37, and S - R0 has order 4 in C0(F_p) for both: every R_k = R0 + 4k(S - R0) reduces to R0, so w(R_k) is p-integral and w(R_k) = 0 mod 11 * 37";
memK := [0, 1];
Rs := [Mul(Q, 25, g, R0, 4*k, S) : k in memK];
for Rk in Rs do assert Rk[2]^3 eq 25*Evaluate(g, Rk[1]); end for;
ws := [Rk[1] : Rk in Rs];
printf "(G4) members R_k for k = %o: w = %o\n", memK, ws;

// ---------------------------------------------------------------- (3) exactness on the class w = 0 mod 11 * 37
primes := [11, 37]; r0 := 0;
finite := [aP, cP, dP, D1, D2, tw, kk];
nonzero := [aP, D1, D2, tw, kk, Discriminant(f1), Discriminant(f2)];
invs := []; fbars := [* *];
for p in primes do
    rp := Integers()!(GF(p)!r0);
    for z in finite do Nz, Dz := IntNumDen(z); assert (GF(p)!Evaluate(Dz, rp)) ne 0; end for;
    for z in nonzero do Nz, Dz := IntNumDen(z); assert (GF(p)!Evaluate(Nz, rp)) ne 0; end for;
    fbar := ReductionClass(ff, p, rp); Append(~fbars, fbar);
    n, chi := FrobData(fbar);
    Ebar := ChangeRing(E, GF(p));
    Fbar := EllipticCurve([GF(p)!Evaluate(Numerator(t), rp)/Evaluate(Denominator(t), rp) : t in aInvariants(Ft)]);
    assert n eq #Ebar*#Fbar;
    Ab := AbelianGroup(Jacobian(HyperellipticCurve(fbar)));
    assert #Ab eq n;
    Append(~invs, Invariants(Ab));
    printf "(3) p = %o, w = %o mod %o: reduction %o (squarefree, degree %o); #J(F_p) = %o = #E(F_p) #F(F_p) = %o * %o; J(F_p) = %o\n",
        p, rp, p, fbar, Degree(fbar), n, #Ebar, #Fbar, Invariants(Ab);
end for;
assert invs eq [[108], [36, 36]];
bnd := StructureBound(invs);
printf "(3) the largest group embedding in all J(F_p) has invariants %o\n", bnd;
assert bnd eq G;
print "(3) every p-integral w = 0 mod 407 (p in [11, 37]) has J(Q)_tors inside [36]";

// ---------------------------------------------------------------- (2), (4), (5) members
function Member(v0)
    s := func< z | Evaluate(z, v0) >;
    R<xx> := PolynomialRing(Q);
    fv := R![s(t) : t in Coefficients(ff)];
    return fv, s;
end function;
function PullClass(fv, s, kv, nn, hh, mm, X0, Y0)   // the class phi^*(X0,Y0) - phi^*(O) on y^2 = fv, where fv kv^2 = tw h1 h2 at the member
    R := Parent(fv); C := HyperellipticCurve(fv); J := Jacobian(C);
    CRl := CoordinateRing(Ambient(C)); Y := CRl.2;
    nnv := R![s(t) : t in Coefficients(nn)]; hhv := R![s(t) : t in Coefficients(hh)]; mmv := R![s(t) : t in Coefficients(mm)];
    k := kv;
    N := nnv - X0*hhv; assert Degree(N) eq 3;
    RN := quo<R | N>; okinv, inv := IsInvertible(RN!(k*mmv)); assert okinv;
    vN := R!((RN!(Y0*hhv^2))*inv);                     // y = Y0 hh^2/(k mm) on the fibre
    assert (vN^2 - fv) mod N eq 0;
    D := Divisor(C, ideal<CRl | Hom3(N, CRl, 3), Y - Hom3(vN, CRl, 3)>);
    W := Divisor(C, ideal<CRl | Hom3(hhv, CRl, 3), Y>);
    assert Degree(D) eq 3 and Degree(W) eq 3;
    return J!(D - W);
end function;
G2s := [];
for i in [1..#ws] do
    v0 := ws[i];
    for j in [1..#primes] do p := primes[j]; assert Valuation(v0, p) ge 0 and (GF(p)!v0) eq GF(p)!r0; end for;
    fv, s := Member(v0);
    assert Degree(fv) eq 6 and Discriminant(fv) ne 0;
    for j in [1..#primes] do assert PolynomialRing(GF(primes[j]))!fv eq fbars[j]; end for;
    Append(~G2s, G2Invariants(HyperellipticCurve(fv)));
    ok, Qv, Fv := Point9(v0); assert ok and Order(Qv) eq 9;
    printf "(2) member k = %o, w = %o: F_w has a rational point Q of order 9 with 3Q = T (x(Q) a rational root of Z at w)\n", memK[i], v0;
end for;
assert #Seqset(G2s) eq #ws;
printf "(4) the G2-invariants of the members k = %o are distinct: R -> M_2 is nonconstant on C0, so the R_k give infinitely many isomorphism classes\n", memK;
// Mumford check of [36] at the member w = 0 and at the base point S; TorsionSubgroup there
for v0 in [Q| 0, S[1]] do
    fv, s := Member(v0);
    d := LCM([Denominator(c) : c in Coefficients(fv)]); fvI := d^2*fv; kv := s(kk)/d;
    Ev := EllipticCurve([s(t) : t in aInvariants(E0t)]); ok, iso := IsIsomorphic(E, Ev); assert ok;
    Fv := EllipticCurve([s(t) : t in aInvariants(Ft)]);
    tv := s(tw);
    ok, Qv := Point9(v0); assert ok;
    P12 := gensE[1]; k := 1; while GCD(k, 12) ne 1 or iso(k*P12)[1] eq 0 do k +:= 1; end while;
    c12 := PullClass(fvI, s, kv, nn1, h1, mm1, iso(k*P12)[1]/tv, iso(k*P12)[2]/tv);
    c9 := PullClass(fvI, s, kv, nn2, h2, mm2, Qv[1]/tv, Qv[2]/tv);
    Jm := Parent(c12);
    assert OrderUpTo(c12, 36) eq 12 and OrderUpTo(c9, 36) eq 9;
    Sg := GeneratedSet(Jm, [c12, c9], 36);
    assert SameGroup(Sg, G);
    invm := TorsionInvariants(fv); assert invm eq G;
    printf "(5) w = %o: phi1^*(P12) (order 12) and phi2^*(Q) (order 9) generate a subgroup of J(Q) isomorphic to [36] (Mumford arithmetic); TorsionSubgroup = %o\n", v0, invm;
end for;
R<xx> := PolynomialRing(Q);
fc, s := Member(S[1]);
assert IsIsomorphic(HyperellipticCurve(fc), HyperellipticCurve(R![9,-9,-1,-5,5,-1,1], R![0,1,1]));
print "(5) control: the base point w = -81/119 (not in the class) is Q-isomorphic to the census curve 36.b";
printf "time %o s\n", Cputime(T0);
print "FAMILY_36_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E x F_w) and J(Q)_tors exactly [36]";
quit;
