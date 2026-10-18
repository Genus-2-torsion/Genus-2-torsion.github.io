/* family_14_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E x F_v) with J(Q)_tors exactly [14].
   Found 2026-10-18 (torsion_jac notes/g33fam_2026_10_18.md); the member v = 0 is the census curve 14.b
   (LMFDB 1352.a.1), which is the (3,3)-gluing of 26b1 and 52a1.

   Source.  (3,3)-gluings of the fixed curve E = 26b1 (torsion Z/7) with partners F that have a rational
   2-torsion point, in the normal form of Broker-Howe-Lauter-Stevenhagen (arXiv:1403.6911, Alg. 5.4 and App. A,
   normalization b = a):
       C : y^2 = tw h1 h2,   h1 = x^3 + 3ax + 2a,   h2 = 2dx^3 + 3cx^2 + 1,   d = (1 - 12ac)/(16a),
   with the degree-3 maps phi1 : (x,y) -> (nn1/h1, y mm1/h1^2) onto Y^2 = tw f1(X) and phi2 : (x,y) ->
   (nn2/h2, y mm2/h2^2) onto Y^2 = tw f2(X), certified by the identities h1^4 f1(nn1/h1) = h1 h2 mm1^2 and
   h2^4 f2(nn2/h2) = h1 h2 mm2^2.  The curve j(E0(a,c)) = j(E) in the (a,c)-plane is rational (it is the
   curve of partners anti-symplectically 3-congruent to E); w parametrizes it (a(w), c(w) of degrees 6, 4).
   The cover "f2 has a rational root X" is rational too (Gamma(3) cap Gamma_0(2) has genus 0); v parametrizes it
   (w(v) of degree 3, X(v) of degree 18), with v = 0 at the BHLS point (a, c) = (2/25, 175/128) of 26b1 x 52a1.
   tw is the twist with E0^tw = E (equal c-invariants), and F_v = F0^tw has the 2-torsion point (tw X(v), 0).

   Argument.
   (1) Splitting.  The identities hold over Q(v); at every member of the residue class below, a, c, d, D1, D2,
       tw, kk, X(v) are finite and a, D1, D2, tw, kk, disc f1, disc f2 are nonzero, so phi1, phi2 are morphisms
       over Q onto E and F_v, and J is isogenous over Q to E x F_v.
   (2) The marked group.  phi1_* phi1^* = [3] on E and phi2_* phi2^* = [3] on F_v, so phi1^* is injective on
       E(Q)[7] = Z/7 and phi2^* on F_v(Q)[2]; the orders are coprime, so Z/14 embeds in J_v(Q) at every member.
       Over Q(v), f_v has a quadratic factor q and [q, 0] has order 2 (Mumford arithmetic over Q(v)).  At the
       members v = 0 and v = 497 the class phi1^*(P7) - phi1^*(O) is computed from divisors: it has order 7,
       and its sum with [q, 0] has order 14.
   (3) Exactness on a residue class.  For every 7- and 71-integral v = 0 (mod 497), C_v has the reductions of
       C_0 at 7 and 71, with #J(F_7) = 70 and #J(F_71) = 4774 (each equal to #E(F_p) #F(F_p)).  J_v(Q)_tors
       injects into both (p odd, good reduction), so its order divides gcd(70, 4774) = 14: it is exactly [14].
   (4) Infinitely many.  The G2-invariants of the members v = 0 and v = 497 differ, so v -> M_2 is nonconstant
       with finite fibres, and the residue class gives infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup at v = 0 and v = -1 is [14]; C_0 is Q-isomorphic to the census curve 14.b.
   Run (from the repository root):
       magma -b pipeline/magma/family_14_split_certificate.m > data/knowledge/sources/family_14_split_certificate.log */
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
E := EllipticCurve([1,-1,1,-3,3]);
assert Conductor(E) eq 26 and Invariants(TorsionSubgroup(E)) eq [7];
TE, mTE := TorsionSubgroup(E);
printf "E = 26b1 = %o, torsion %o\n", aInvariants(E), Invariants(TE);

// ---------------------------------------------------------------- (1) the family over Q(v)
Fv<v> := FunctionField(Q); w := v;
aW := (-1/9*w^6 + 1849/2400*w^5 - 2273969/1920000*w^4 - 6379403159/4608000000*w^3 + 1271148671/600000000*w^2 + 577777369/234375000*w + 1397415032/2197265625)/(w^6 - 7539/800*w^5 + 22230027/640000*w^4 - 31536346753/512000000*w^3 + 2682153733857/51200000000*w^2 - 455866040907/25600000000*w + 36106863257/51200000000);
cW := (-3/16*w^4 + 13017/6400*w^3 - 58959729/10240000*w^2 + 31558929/5120000*w - 23065329/10240000)/(w^4 - 1849/400*w^3 + 1129137/640000*w^2 + 1033591/125000*w + 1249924/390625);
wV := (1534043/1656800*v^3 - 7253577/1656800*v^2 + 5502849/1656800*v + 2009417/331360)/(v^3 - 765/2071*v^2 + 6069/2071*v + 4913/2071);
XV := (-3851809607954257291719/815705480820877915888640*v^18 - 17487261347439959680293/81570548082087791588864*v^17 + 1780241438769185158194753/815705480820877915888640*v^16 - 183151422270803908940499/50981592551304869743040*v^15 - 184985141046736477178511/40785274041043895794432*v^14 - 2890275766575198440931999/101963185102609739486080*v^13 + 2589869195718107443309845/40785274041043895794432*v^12 + 6927794658822271713903099/50981592551304869743040*v^11 - 58740869626080084702069057/407852740410438957944320*v^10 - 11688381431539738424132751/40785274041043895794432*v^9 + 46317775197128553782620047/407852740410438957944320*v^8 + 15474167072195640918529563/50981592551304869743040*v^7 - 2483972829767286484280199/203926370205219478972160*v^6 - 2889550848852812357438355/20392637020521947897216*v^5 - 2583071272155124577917467/203926370205219478972160*v^4 + 1004710491079945343343/79287080173102441280*v^3 - 11403798977254406295390831/815705480820877915888640*v^2 + 2613251987298937838035143/407852740410438957944320*v + 1645380880891923823948053/163141096164175583177728)/(v^18 + 6976757730415926/23374728824617465*v^17 - 6940618748342548779897/796587383614138589735*v^16 - 13574057890422116143632/796587383614138589735*v^15 + 10722494484471946419948/796587383614138589735*v^14 + 92783732977489488629784/796587383614138589735*v^13 - 1056015965868330519156/159317476722827717947*v^12 - 103245966552696798503088/159317476722827717947*v^11 - 660534164032577175257982/796587383614138589735*v^10 + 1244849862482286620292668/796587383614138589735*v^9 + 5866390260401711629571346/796587383614138589735*v^8 + 9875760709591888569448848/796587383614138589735*v^7 + 9643026379574418949400412/796587383614138589735*v^6 + 5877606945593116276270872/796587383614138589735*v^5 + 1924293867630601240209612/796587383614138589735*v^4 - 1511359371197182738704/796587383614138589735*v^3 - 298244445005962236920049/796587383614138589735*v^2 - 24817207856590102925310/159317476722827717947*v - 14063084452067724991009/796587383614138589735);
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
TFt := Ft![tw*XV, 0];
assert 2*TFt eq Ft!0 and TFt ne Ft!0;
print "(1) E0^tw has the c-invariants of E; F = F0^tw has the rational 2-torsion point (tw X(v), 0)";
rep, rt := SqfRep(tw);
PX<x> := Parent(h1);
ff := rep*h1*h2; den := LCM([Denominator(t) : t in Coefficients(ff)]); ff := ff*den^2;
Pv<V> := PolynomialRing(Q);
gc := GCD([Pv!Evaluate(Numerator(t), V) : t in Coefficients(ff)]);
sq := &*([Pv| 1] cat [fa[1]^(fa[2] div 2) : fa in Factorization(gc)]);
sqv := Evaluate(sq, v); ff := ff/sqv^2;
assert &and[Denominator(t) eq 1 : t in Coefficients(ff)];
kk := rt*sqv/den;
assert ff*kk^2 eq tw*h1*h2;
printf "(1) C_v : y^2 = f_v(x), coefficients of degrees %o in v; y' = kk y maps C_v to y'^2 = tw h1 h2\n", [Degree(Numerator(t)) : t in Coefficients(ff)];
printf "(1) phi1 : (x,y) -> (tw nn1/h1, tw kk y mm1/h1^2) onto E0^tw ~ E, phi2 : (x,y) -> (tw nn2/h2, tw kk y mm2/h2^2) onto F\n";
fac := Factorization(ff);
printf "(1) factor degrees of f_v over Q(v): %o\n", [Degree(g[1]) : g in fac];
quads := [g[1] : g in fac | Degree(g[1]) eq 2] cat [g[1]*h[1] : g, h in fac | Degree(g[1]) eq 1 and Degree(h[1]) eq 1 and g[1] ne h[1]];
assert #quads ge 1;
qv := quads[1];
Jv := Jacobian(HyperellipticCurve(ff));
c2v := elt< Jv | qv, 0 >;
assert c2v ne Jv!0 and 2*c2v eq Jv!0;
print "(1) marked over Q(v): the class c2 = [q, 0] of a quadratic factor q of f_v has order 2 (Mumford arithmetic over Q(v))";

// ---------------------------------------------------------------- (2) exactness on the residue class v = 0 mod [7, 71]
primes := [7, 71]; r0 := 0;
finite := [aV, cV, dV, D1, D2, tw, kk, XV];                       // regular at the class: denominators are units
nonzero := [aV, D1, D2, tw, kk, Discriminant(f1), Discriminant(f2)];  // and these do not vanish
nJ := [];
for p in primes do
    for z in finite do Nz, Dz := IntNumDen(z); assert (GF(p)!Evaluate(Dz, r0)) ne 0; end for;
    for z in nonzero do Nz, Dz := IntNumDen(z); assert (GF(p)!Evaluate(Nz, r0)) ne 0; end for;
    fbar := ReductionClass(ff, p, r0);
    n, chi := FrobData(fbar);
    Ebar := ChangeRing(E, GF(p));
    Fbar := EllipticCurve([GF(p)!Evaluate(Numerator(t), r0)/Evaluate(Denominator(t), r0) : t in aInvariants(Ft)]);
    assert n eq #Ebar*#Fbar;
    Append(~nJ, n);
    printf "(2) p = %o, v = %o mod %o: reduction %o (squarefree, degree %o); #J(F_p) = %o = #E(F_p) #F(F_p) = %o * %o; J(F_p) = %o\n",
        p, r0, p, fbar, Degree(fbar), n, #Ebar, #Fbar, Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
end for;
assert GCD(nJ) eq 14;
printf "(2) gcd of #J(F_p) = %o = #G: every p-integral v = 0 mod %o (p in %o) has J(Q)_tors = G = [14]\n", GCD(nJ), &*primes, primes;

// ---------------------------------------------------------------- (3), (4) members
CR := 0;
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
    vN := R!((RN!(Y0*hhv^2))*inv);                     // y = Y0 hh^2/(k mm) on the fibre
    assert (vN^2 - fv) mod N eq 0;
    D := Divisor(C, ideal<CRl | Hom3(N, CRl, 3), Y - Hom3(vN, CRl, 3)>);
    W := Divisor(C, ideal<CRl | Hom3(hhv, CRl, 3), Y>);
    assert Degree(D) eq 3 and Degree(W) eq 3;
    return J!(D - W);
end function;
members := [0, 497];
G2s := [];
for v0 in members do
    fv, s := Member(v0);
    assert Degree(fv) ge 5 and Discriminant(fv) ne 0;
    Append(~G2s, G2Invariants(HyperellipticCurve(fv)));
    Ev := EllipticCurve([s(t) : t in aInvariants(E0t)]); Fvv := EllipticCurve([s(t) : t in aInvariants(Ft)]);
    ok, iso := IsIsomorphic(E, Ev); assert ok;
    classes := [];
    for P in [mTE(g) : g in Generators(TE)] do
        Q0 := iso(P); tv := s(tw);
        Append(~classes, PullClass(fv, s, nn1, h1, mm1, Q0[1]/tv, Q0[2]/tv));
    end for;
    R2 := Parent(fv); Jm := Jacobian(HyperellipticCurve(fv));
    Append(~classes, elt< Jm | R2![s(c) : c in Coefficients(qv)], 0 >);
    o := [OrderUpTo(c, 14) : c in classes];
    assert o eq [7, 2] and OrderUpTo(classes[1] + classes[2], 14) eq 14;
    printf "(4) v = %o: phi1^*(E(Q)_tors) and [q, 0] have orders %o; their sum has order 14 (Mumford arithmetic on the member)\n", v0, o;
end for;
assert #Seqset(G2s) eq #members;
printf "(3) the G2-invariants of the members %o are pairwise distinct: v -> M_2 is nonconstant, so the residue class gives infinitely many isomorphism classes\n", members;
G2t := [];
for v0 in [0, -1] do
    fv, s := Member(v0);
    invm := TorsionInvariants(fv);
    Append(~G2t, G2Invariants(HyperellipticCurve(fv)));
    printf "(5) TorsionSubgroup at v = %o: %o
", v0, invm;
    assert invm eq [14];
end for;
assert #Seqset(G2t) eq #G2t;
print "(5) these members are pairwise non-isomorphic (distinct G2-invariants)";
fc, s := Member(0);
R<xx> := PolynomialRing(Q);
Ccensus := HyperellipticCurve(R![0,-4,1,-3,1,-1], R![0,0,1,1]);
assert IsIsomorphic(HyperellipticCurve(fc), Ccensus);
print "(4) the member v = 0 is Q-isomorphic to the census curve 14.b";
printf "time %o s\n", Cputime(T0);
print "FAMILY_14_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E x F_v) and J(Q)_tors exactly [14]";
quit;
