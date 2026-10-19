/* family_2_2_8_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E x F_v) with J(Q)_tors exactly [2,2,8].
   Found 2026-10-19 (torsion_jac notes/g33rat_2026_10_19.md; builder code/g33rat_build.m).

   Source.  (3,3)-gluings of the fixed curve E = 210e2 (torsion Z/2 x Z/8) with partners F that have a rational
   2-torsion point, in the normal form of Broker-Howe-Lauter-Stevenhagen (arXiv:1403.6911, Alg. 5.4 and App. A,
   normalization b = a):
       C : y^2 = tw h1 h2,   h1 = x^3 + 3ax + 2a,   h2 = 2dx^3 + 3cx^2 + 1,   d = (1 - 12ac)/(16a),
   with the degree-3 maps phi1 : (x,y) -> (nn1/h1, y mm1/h1^2) onto Y^2 = tw f1(X) and phi2 : (x,y) ->
   (nn2/h2, y mm2/h2^2) onto Y^2 = tw f2(X), certified by the identities h1^4 f1(nn1/h1) = h1 h2 mm1^2 and
   h2^4 f2(nn2/h2) = h1 h2 mm2^2.  The curve j(E0(a,c)) = j(E) in the (a,c)-plane is rational (the partners
   anti-symplectically 3-congruent to E); w parametrizes it (a(w), c(w) of degrees 6, 4).  The cover "f2 has a
   rational root X" is rational too (Gamma(3) cap Gamma_0(2) has genus 0, and the cover has odd degree 3 over the
   w-line); v parametrizes it.  tw is the twist with E0^tw = E (equal c-invariants), and F_v = F0^tw has the
   2-torsion point (tw X(v), 0).

   Argument.
   (1) Splitting.  The identities hold over Q(v); at every member of the residue class below, a, c, d, D1, D2,
       tw, kk, X(v) are finite and a, D1, D2, tw, kk, disc f1, disc f2 are nonzero, so phi1, phi2 are
       morphisms over Q onto E and F_v, and J is isogenous over Q to E x F_v.
   (2) The marked group.  phi1 has prime degree 3 and C is not of genus 1, so phi1 does not factor through an
       isogeny and phi1^* : E -> J is injective; so phi1^*(E(Q)_tors) = [2,8] embeds in J_v(Q) at every member.
       Over Q(v), f_v has the factor degrees printed below, so J(Q(v))[2] has 2-rank 3 (its nonzero rational points
       are the classes [u, 0] of the quadratic divisors u of f_v; each is checked to have order 2 by Mumford arithmetic
       over Q(v)).  These classes stay distinct at every member (f_v squarefree), so J_v(Q)[2] has 2-rank >= 3,
       one more than phi1^*(E(Q)_tors); a class c of order 2 outside phi1^*(E(Q)_tors) gives
       phi1^*(E(Q)_tors) + <c> = [2,2,8] in J_v(Q).  At the members the pullbacks of the generators of E(Q)_tors
       ((34, 88) of order 8, (4, 58) of order 8) and the classes [u, 0] are computed by Mumford arithmetic and generate a group isomorphic to [2,2,8].
   (3) Exactness on a residue class.  For every 11- and 13-integral v = -2 (mod 143), C_v has the
       reductions of C_{-2} at 11 and 13, with J(F_11) = [2, 2, 40] and J(F_13) = [2, 8, 16]
       (each of order #E(F_p) #F(F_p)).  J_v(Q)_tors embeds in both (p odd, good reduction).  If A embeds in G_1
       and G_2, then for every prime l the sorted l-primary exponents of A are bounded componentwise by those of each
       G_i; the largest such group (StructureBound) is [2,2,8], which is the marked group, so J_v(Q)_tors = [2,2,8].
   (4) Infinitely many.  The G2-invariants of the members v = -2 and v = 141 differ, so v -> M_2 is
       nonconstant with finite fibres, and the residue class gives infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup at v = -2 and v = 141 is [2,2,8].
       The census Q-split curve 2.2.8.b is not a member: 210e2 is not an isogeny factor of its Jacobian (compared
       L-polynomials; g33rat note, census_ctrl.log).
   Run (from the repository root):
       magma -b pipeline/magma/family_2_2_8_split_certificate.m > data/knowledge/sources/family_2_2_8_split_certificate.log */
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
G := [2, 2, 8];
E := EllipticCurve([1,0,0,-1070,7812]);
TE := [2, 8];
assert Conductor(E) eq 210 and Invariants(TorsionSubgroup(E)) eq TE;
gensE := [E![34, 88], E![4, 58]];
assert [Order(P) : P in gensE] eq [8, 8];
if #gensE gt 0 then
    assert #{ &+[k[i]*gensE[i] : i in [1..#gensE]] : k in CartesianProduct([[0..Order(P)-1] : P in gensE]) } eq &*TE;   // the generators generate E(Q)_tors
end if;
printf "E = 210e2 = %o, conductor %o, torsion %o; generators used for the pullbacks: %o (2-torsion classes come from f)\n", aInvariants(E), Conductor(E), TE, gensE;

// ---------------------------------------------------------------- (1) the family over Q(v)
Fv<v> := FunctionField(Q); w := v;
aW := (-186410361755524/18129006550947*w^6 - 9202268932/557112767*w^5 - 516484/51361*w^4 - 8976152828/3121259331*w^3 - 945066567864275712/2393555647531367041*w^2 - 1234699882918915719168/48486256752042902149537*w - 69428221442031346129346494464/111545893534574906495450088321121)/(w^4 + 2903221244/3121259331*w^3 + 603580475049972480/2393555647531367041*w^2 + 1234699882918915719168/48486256752042902149537*w + 69428221442031346129346494464/111545893534574906495450088321121);
cW := (-6043002183649/745641447022096*w^2 - 2480990955582312/944028674520412417*w - 837049351876342425363456/2171801435614082796196532161)/(w^4 + 21694/20257*w^3 + 17021850332257/46602590438881*w^2 + 39695855289316992/944028674520412417*w + 3348197407505369701453824/2171801435614082796196532161);
wV := (2342952/6826609*v^2 - 2342952/6826609)/(v^3 - 337/72*v^2 + 3*v);
XV := (30625/4608*v^16 - 10320625/110592*v^15 + 51565058125/95551488*v^14 - 1921132740625/1146617856*v^13 + 1493576018790625/495338913792*v^12 - 1707124260625/573308928*v^11 + 3901410625/5308416*v^10 + 7358574663125/3439853568*v^9 - 1631233639590625/495338913792*v^8 + 2054888040625/859963392*v^7 - 93966044375/95551488*v^6 + 72244375/331776*v^5 - 30625/1536*v^4)/(v^18 - 337/36*v^17 + 129121/5184*v^16 + 337/9*v^15 - 8236951/20736*v^14 + 1678846043/1492992*v^13 - 9906029/6912*v^12 - 80523791/62208*v^11 + 429799793279/35831808*v^10 - 52029105469/1492992*v^9 + 441633698495/6718464*v^8 - 43355961819889/483729408*v^7 + 4925949272993/53747712*v^6 - 68744187747823/967458816*v^5 + 5755236359614561/139314069504*v^4 - 5607812276881/322486272*v^3 + 178067880211/35831808*v^2 - 142124695/165888*v + 51361/768);
aP := Evaluate(aW, wV); cP := Evaluate(cW, wV);
h1, h2, f1, f2, nn1, mm1, nn2, mm2, dP, D1, D2 := BHLSParts(aP, cP);
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
assert Degree(ff) eq 6 and &and[g[2] eq 1 : g in fac];
lins := [g[1] : g in fac | Degree(g[1]) eq 1];
quads := [g[1] : g in fac | Degree(g[1]) eq 2] cat [lins[i]*lins[j] : j in [i+1..#lins], i in [1..#lins]];
printf "(1) factor degrees of f_v over Q(v): %o; %o quadratic divisors, so #J(Q(v))[2] = %o\n", [Degree(g[1]) : g in fac], #quads, 1 + #quads;
assert 1 + #quads eq 2^3;
Jv := Jacobian(HyperellipticCurve(ff));
for u in quads do cu := elt< Jv | u, 0 >; assert cu ne Jv!0 and 2*cu eq Jv!0; end for;
print "(1) marked over Q(v): each class [u, 0] has order 2 (Mumford arithmetic over Q(v)); J(Q(v))[2] has 2-rank 3, one more than E(Q)_tors";

// ---------------------------------------------------------------- (3) exactness on the residue class v = -2 mod 143
primes := [11, 13]; r0 := -2;
finite := [aP, cP, dP, D1, D2, tw, kk, XV];                       // regular at the class: denominators are units
nonzero := [aP, D1, D2, tw, kk, Discriminant(f1), Discriminant(f2)];  // and these do not vanish
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
    A := AbelianGroup(Jacobian(HyperellipticCurve(fbar)));
    assert #A eq n;
    Append(~invs, Invariants(A));
    printf "(3) p = %o, v = %o mod %o: reduction %o (squarefree, degree %o); #J(F_p) = %o = #E(F_p) #F(F_p) = %o * %o; J(F_p) = %o\n",
        p, rp, p, fbar, Degree(fbar), n, #Ebar, #Fbar, Invariants(A);
end for;
assert invs eq [[2, 2, 40], [2, 8, 16]];
bnd := StructureBound(invs);
printf "(3) the largest group embedding in all J(F_p) has invariants %o\n", bnd;
assert bnd eq G;
printf "(3) every p-integral v = %o mod %o (p in %o) has J(Q)_tors = G = [2,2,8]\n", r0, &*primes, primes;

// ---------------------------------------------------------------- (2), (4) members
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
members := [Q| -2, 141];
for v0 in members do for i in [1..#primes] do
    p := primes[i]; assert Denominator(v0) eq 1 and (GF(p)!v0) eq (GF(p)!r0);
end for; end for;
G2s := [];
for v0 in members do
    fv, s := Member(v0);
    assert Degree(fv) eq 6 and Discriminant(fv) ne 0;
    for i in [1..#primes] do assert PolynomialRing(GF(primes[i]))!fv eq fbars[i]; end for;
    Append(~G2s, G2Invariants(HyperellipticCurve(fv)));
    Ev := EllipticCurve([s(t) : t in aInvariants(E0t)]);
    ok, iso := IsIsomorphic(E, Ev); assert ok;
    tv := s(tw);
    R2 := Parent(fv); Jm := Jacobian(HyperellipticCurve(fv));
    cls := [Jm| ];
    for P in gensE do       // a generator k P (k prime to the order) whose image has X != 0, so that the fibre has degree 3
        k := 1; while GCD(k, Order(P)) ne 1 or iso(k*P)[1] eq 0 do k +:= 1; end while;
        Q0 := iso(k*P);
        Append(~cls, PullClass(fv, s, nn1, h1, mm1, Q0[1]/tv, Q0[2]/tv));
    end for;
    assert [OrderUpTo(c, 24) : c in cls] eq [Order(P) : P in gensE];
    cl2 := [elt< Jm | R2![s(c) : c in Coefficients(u)], 0 > : u in quads];
    S := GeneratedSet(Jm, cls cat cl2, 24);
    assert SameGroup(S, G);
    Fm := EllipticCurve([s(t) : t in aInvariants(Ft)]);
    printf "(2) v = %o: the pullbacks phi1^*(P) of the generators have orders %o; with the %o classes [u, 0] they generate a subgroup of J(Q) of order %o isomorphic to %o (Mumford arithmetic on the member); F_v has torsion %o\n", v0, [OrderUpTo(c, 24) : c in cls], #cl2, #S, G, Invariants(TorsionSubgroup(Fm));
end for;
assert #Seqset(G2s) eq #members;
printf "(4) the G2-invariants of the members %o are pairwise distinct: v -> M_2 is nonconstant, so the residue class gives infinitely many isomorphism classes\n", members;
G2t := [];
for v0 in [Q| -2, 141] do
    fv, s := Member(v0);
    invm := TorsionInvariants(fv);
    Append(~G2t, G2Invariants(HyperellipticCurve(fv)));
    printf "(5) TorsionSubgroup at v = %o: %o\n", v0, invm;
    assert invm eq G;
end for;
assert #Seqset(G2t) eq #G2t;
print "(5) these members are pairwise non-isomorphic (distinct G2-invariants)";
printf "time %o s\n", Cputime(T0);
print "FAMILY_2_2_8_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E x F_v) and J(Q)_tors exactly [2,2,8]";
quit;
