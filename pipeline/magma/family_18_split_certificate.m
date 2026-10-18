/* family_18_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E x F_v) with J(Q)_tors exactly [18].
   Found 2026-10-18 (torsion_jac notes/s18_2026_10_18.md).

   Source.  (3,3)-gluings of the fixed curve E = 54b3 (torsion Z/9) with partners F that have a rational
   2-torsion point, in the normal form of Broker-Howe-Lauter-Stevenhagen (arXiv:1403.6911, Alg. 5.4 and App. A,
   normalization b = a):
       C : y^2 = tw h1 h2,   h1 = x^3 + 3ax + 2a,   h2 = 2dx^3 + 3cx^2 + 1,   d = (1 - 12ac)/(16a),
   with the degree-3 maps phi1 : (x,y) -> (nn1/h1, y mm1/h1^2) onto Y^2 = tw f1(X) and phi2 : (x,y) ->
   (nn2/h2, y mm2/h2^2) onto Y^2 = tw f2(X), certified by the identities h1^4 f1(nn1/h1) = h1 h2 mm1^2 and
   h2^4 f2(nn2/h2) = h1 h2 mm2^2.  The curve j(E0(a,c)) = j(E) in the (a,c)-plane is rational (the partners
   anti-symplectically 3-congruent to E); w parametrizes it (a(w), c(w) of degrees 6, 4).  The cover "f2 has a
   rational root X" is rational too (Gamma(3) cap Gamma_0(2) has genus 0, and the cover has odd degree 3 over
   the w-line); v parametrizes it (w(v) of degree 3, X(v) of degree 18).  tw is the twist with E0^tw = E (equal
   c-invariants), and F_v = F0^tw has the 2-torsion point (tw X(v), 0).  The rational 3-torsion point 3P of E
   gives a rational 3-torsion point on F_v (the Galois isomorphism E[3] -> F_v[3]), so F_v(Q) contains Z/6.

   Argument.
   (1) Splitting.  The identities hold over Q(v); at every member of the residue class below, a, c, d, D1, D2,
       tw, kk, X(v) are finite and a, D1, D2, tw, kk, disc f1, disc f2 are nonzero, so phi1, phi2 are morphisms
       over Q onto E and F_v, and J is isogenous over Q to E x F_v.
   (2) The marked group.  phi1 has prime degree 3 and C is not of genus 1, so phi1 does not factor through an
       isogeny and phi1^* : E -> J is injective (equivalently J = (E x F_v)/Gamma with Gamma the graph of
       E[3] -> F_v[3], which meets E x 0 trivially); so phi1^*(E(Q)_tors) = Z/9 in J_v(Q).  Over Q(v), f_v
       has a quadratic factor q and c2 = [q, 0] has order 2 (Mumford arithmetic over Q(v)), so Z/18 embeds in
       J_v(Q) at every member.  At the members v = -2 and v = 63 the class phi1^*(P) - phi1^*(O), P = (9, 19)
       of order 9, is computed from divisors: it has order 9, and its sum with [q, 0] has order 18.
   (3) Exactness on a residue class.  #J(F_p) = #E(F_p) #F(F_p) is divisible by 9 * 6, so a gcd of orders
       cannot cap the torsion at 18; the group structures are used.  For every 5- and 13-integral v = -2
       (mod 65), C_v has the reductions of C_{-2} at 5 and 13, with J(F_5) = Z/54 (cyclic) and
       J(F_13) = Z/18 x Z/18 (exponent 18).  J_v(Q)_tors embeds in both (p odd, good reduction), so it is
       cyclic of order dividing 18: it is exactly [18].  (In general: if A embeds in G_1 and G_2, the sorted
       l-primary exponents of A are bounded componentwise by those of each G_i; the script computes this bound.)
   (4) Infinitely many.  The G2-invariants of the members v = -2 and v = 63 differ, so v -> M_2 is
       nonconstant with finite fibres, and the residue class gives infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup at v = -2 and v = -1 is [18] (the generic torsion of the family is [18]; the
       rational 3-torsion of F_v brings no extra 3-torsion into J, see the note).  No census curve is a member:
       every member has the isogeny factor 54b3, so 54 divides its conductor, while the census curve 18.b has
       conductor 490 (its Jacobian is isogenous to 14a x 35a; see the note).  Control: the member v = 4 (outside
       the class) is the census curve 36.b, the gluing of 54b3 and 90c3; there F_4 has torsion Z/12 and
       J(Q)_tors = [36].
   Run (from the repository root):
       magma -b pipeline/magma/family_18_split_certificate.m > data/knowledge/sources/family_18_split_certificate.log */
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
E := EllipticCurve([1,-1,1,-14,29]);
assert Conductor(E) eq 54 and Invariants(TorsionSubgroup(E)) eq [9];
P9 := E![9, 19];
assert Order(P9) eq 9;
printf "E = 54b3 = %o, torsion %o, P = %o of order 9\n", aInvariants(E), Invariants(TorsionSubgroup(E)), P9;

// ---------------------------------------------------------------- (1) the family over Q(v)
Fv<v> := FunctionField(Q); w := v;
aW := (-1024/1971*w^6 - 64/73*w^5 + 28/73*w^4 + 261/292*w^3 - 63/292*w^2 - 81/292*w + 27/292)/(w^4 - 85/292*w^3 - 585/292*w^2 + 513/292*w - 135/292);
cW := (-657/4096*w^2 + 513/2048*w - 81/4096)/(w^4 + 9/8*w^3 - 207/256*w^2 - 81/128*w + 81/256);
wV := (v^3 - 3*v^2 - 16)/(v^3 - 3*v^2 - 16*v);
XV := (27/2048*v^18 - 243/1024*v^17 + 3213/2048*v^16 - 2133/512*v^15 + 8829/2048*v^14 - 33723/1024*v^13 + 230283/2048*v^12 + 5589/8*v^11 - 969651/256*v^10 + 95553/128*v^9 + 4911273/256*v^8 - 19575/2*v^7 - 21249*v^6 + 36288*v^5 + 27648*v^4)/(v^16 - 16*v^15 + 123*v^14 - 692*v^13 + 3101*v^12 - 10752*v^11 + 29488*v^10 - 59818*v^9 + 61083*v^8 + 95660*v^7 - 358607*v^6 + 75810*v^5 + 941867*v^4 - 439520*v^3 - 1858368*v^2 + 324608*v + 1196032);
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
assert #quads eq 1;
qv := quads[1];
Jv := Jacobian(HyperellipticCurve(ff));
c2v := elt< Jv | qv, 0 >;
assert c2v ne Jv!0 and 2*c2v eq Jv!0;
print "(1) marked over Q(v): the class c2 = [q, 0] of the quadratic factor q of f_v has order 2 (Mumford arithmetic over Q(v)); it is the only rational 2-torsion class";

// ---------------------------------------------------------------- (3) exactness on the residue class v = -2 mod [5, 13]
primes := [5, 13]; vres := -2;
finite := [aV, cV, dV, D1, D2, tw, kk, XV];                       // regular at the class: denominators are units
nonzero := [aV, D1, D2, tw, kk, Discriminant(f1), Discriminant(f2)];  // and these do not vanish
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
assert bnd eq [18];
printf "(3) every p-integral v = %o mod %o (p in %o) has J(Q)_tors = G = [18]\n", vres, &*primes, primes;

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
members := [Q| -2, 63];
for v0 in members do for p in primes do assert Valuation(v0, p) ge 0 and (GF(p)!v0) eq (GF(p)!vres); end for; end for;
G2s := [];
for v0 in members do
    fv, s := Member(v0);
    assert Degree(fv) ge 5 and Discriminant(fv) ne 0;
    Append(~G2s, G2Invariants(HyperellipticCurve(fv)));
    Ev := EllipticCurve([s(t) : t in aInvariants(E0t)]);
    ok, iso := IsIsomorphic(E, Ev); assert ok;
    Q0 := iso(P9); tv := s(tw);
    c9 := PullClass(fv, s, nn1, h1, mm1, Q0[1]/tv, Q0[2]/tv);
    R2 := Parent(fv); Jm := Jacobian(HyperellipticCurve(fv));
    c2 := elt< Jm | R2![s(c) : c in Coefficients(qv)], 0 >;
    o := [OrderUpTo(c9, 18), OrderUpTo(c2, 18)];
    assert o eq [9, 2] and OrderUpTo(c9 + c2, 18) eq 18;
    Fm := EllipticCurve([s(t) : t in aInvariants(Ft)]);
    printf "(2) v = %o: phi1^*(P) and [q, 0] have orders %o; their sum has order 18 (Mumford arithmetic on the member); F_v = %o has torsion %o\n", v0, o, aInvariants(MinimalModel(Fm)), Invariants(TorsionSubgroup(Fm));
end for;
assert #Seqset(G2s) eq #members;
printf "(4) the G2-invariants of the members %o are pairwise distinct: v -> M_2 is nonconstant, so the residue class gives infinitely many isomorphism classes\n", members;
G2t := [];
for v0 in [-2, -1] do
    fv, s := Member(v0);
    invm := TorsionInvariants(fv);
    Append(~G2t, G2Invariants(HyperellipticCurve(fv)));
    printf "(5) TorsionSubgroup at v = %o: %o\n", v0, invm;
    assert invm eq [18];
end for;
assert #Seqset(G2t) eq #G2t;
print "(5) these members are pairwise non-isomorphic (distinct G2-invariants)";
fc, s := Member(4);
R<xx> := PolynomialRing(Q);
Ccensus := HyperellipticCurve(R![9,-9,-1,-5,5,-1,1], R![0,1,1]);
assert IsIsomorphic(HyperellipticCurve(fc), Ccensus);
F4 := EllipticCurve([s(t) : t in aInvariants(Ft)]);
inv4 := TorsionInvariants(fc);
assert inv4 eq [36] and Invariants(TorsionSubgroup(F4)) eq [12];
printf "(5) control: the member v = 4 (not in the class) is Q-isomorphic to the census curve 36.b; there F_4 = %o has torsion [12] and J(Q)_tors = %o\n", aInvariants(MinimalModel(F4)), inv4;
printf "time %o s\n", Cputime(T0);
print "FAMILY_18_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E x F_v) and J(Q)_tors exactly [18]";
quit;
