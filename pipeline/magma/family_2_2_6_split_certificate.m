/* family_2_2_6_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E x F_v) with J(Q)_tors exactly [2,2,6].
   Found 2026-10-19 (torsion_jac notes/g33rat_2026_10_19.md; builder code/g33rat_build.m).

   Source.  (3,3)-gluings of the fixed curve E = 30a2 (torsion Z/2 x Z/6) with partners F that have a rational
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
       isogeny and phi1^* : E -> J is injective; so phi1^*(E(Q)_tors) = [2,6] embeds in J_v(Q) at every member.
       Over Q(v), f_v has the factor degrees printed below, so J(Q(v))[2] has 2-rank 3 (its nonzero rational points
       are the classes [u, 0] of the quadratic divisors u of f_v; each is checked to have order 2 by Mumford arithmetic
       over Q(v)).  These classes stay distinct at every member (f_v squarefree), so J_v(Q)[2] has 2-rank >= 3,
       one more than phi1^*(E(Q)_tors); a class c of order 2 outside phi1^*(E(Q)_tors) gives
       phi1^*(E(Q)_tors) + <c> = [2,2,6] in J_v(Q).  At the members the pullbacks of the generators of E(Q)_tors
       ((1, 2) of order 6, (13, 38) of order 6) and the classes [u, 0] are computed by Mumford arithmetic and generate a group isomorphic to [2,2,6].
       F_v also inherits a rational point of order 3 from E (the Galois isomorphism E[3] -> F_v[3]); no extra
       3-torsion appears in J on the residue class (step (3) bounds the 3-part).
   (3) Exactness on a residue class.  For every 13- and 17-integral v = 3 (mod 221), C_v has the
       reductions of C_{3} at 13 and 17, with J(F_13) = [2, 6, 12] and J(F_17) = [2, 2, 54]
       (each of order #E(F_p) #F(F_p)).  J_v(Q)_tors embeds in both (p odd, good reduction).  If A embeds in G_1
       and G_2, then for every prime l the sorted l-primary exponents of A are bounded componentwise by those of each
       G_i; the largest such group (StructureBound) is [2,2,6], which is the marked group, so J_v(Q)_tors = [2,2,6].
   (4) Infinitely many.  The G2-invariants of the members v = 3 and v = -218 differ, so v -> M_2 is
       nonconstant with finite fibres, and the residue class gives infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup at v = 3 and v = -218 is [2,2,6].
       Control: the member v = -1 (outside the class; w(-1) = -70/11) is Q-isomorphic to the census curve
       2.2.6.b, with J(Q)_tors = [2,2,6]; its reductions do not give a structure bound [2,2,6] at small primes.
   Run (from the repository root):
       magma -b pipeline/magma/family_2_2_6_split_certificate.m > data/knowledge/sources/family_2_2_6_split_certificate.log */
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
G := [2, 2, 6];
E := EllipticCurve([1,0,1,-19,26]);
TE := [2, 6];
assert Conductor(E) eq 30 and Invariants(TorsionSubgroup(E)) eq TE;
gensE := [E![1, 2], E![13, 38]];
assert [Order(P) : P in gensE] eq [6, 6];
if #gensE gt 0 then
    assert #{ &+[k[i]*gensE[i] : i in [1..#gensE]] : k in CartesianProduct([[0..Order(P)-1] : P in gensE]) } eq &*TE;   // the generators generate E(Q)_tors
end if;
printf "E = 30a2 = %o, conductor %o, torsion %o; generators used for the pullbacks: %o (2-torsion classes come from f)\n", aInvariants(E), Conductor(E), TE, gensE;

// ---------------------------------------------------------------- (1) the family over Q(v)
Fv<v> := FunctionField(Q); w := v;
aW := (-192/43561*w^6 + 1776/6223*w^5 - 5956/889*w^4 + 77293/1143*w^3 - 104230/381*w^2 + 181300/381*w - 343000/1143)/(w^4 - 37333/1143*w^3 + 79030/381*w^2 - 181300/381*w + 343000/1143);
cW := (-43561/2304*w^2 + 63455/576*w - 60025/144)/(w^4 - 259/6*w^3 + 78841/144*w^2 - 63455/36*w + 60025/36);
wV := (v^3 - 7*v^2 - 27*v + 81)/(v^3 - 7*v^2 + 54/7*v);
XV := (175/17496*v^18 - 1225/2916*v^17 + 43325/5832*v^16 - 155575/2187*v^15 + 1118425/2916*v^14 - 1352825/1458*v^13 - 16186975/8748*v^12 + 1956325/81*v^11 - 21380825/216*v^10 + 236606825/972*v^9 - 763901425/1944*v^8 + 11243300/27*v^7 - 1657775/6*v^6 + 102900*v^5 - 16200*v^4)/(v^16 - 112/3*v^15 + 2009/3*v^14 - 23054/3*v^13 + 567077/9*v^12 - 3505460/9*v^11 + 16814419/9*v^10 - 572293414/81*v^9 + 573025453/27*v^8 - 4099626664/81*v^7 + 23068744531/243*v^6 - 33377126258/243*v^5 + 36128177029/243*v^4 - 1036889812/9*v^3 + 59863153*v^2 - 18577242*v + 2592324);
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

// ---------------------------------------------------------------- (3) exactness on the residue class v = 3 mod 221
primes := [13, 17]; r0 := 3;
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
assert invs eq [[2, 6, 12], [2, 2, 54]];
bnd := StructureBound(invs);
printf "(3) the largest group embedding in all J(F_p) has invariants %o\n", bnd;
assert bnd eq G;
printf "(3) every p-integral v = %o mod %o (p in %o) has J(Q)_tors = G = [2,2,6]\n", r0, &*primes, primes;

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
members := [Q| 3, -218];
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
for v0 in [Q| 3, -218] do
    fv, s := Member(v0);
    invm := TorsionInvariants(fv);
    Append(~G2t, G2Invariants(HyperellipticCurve(fv)));
    printf "(5) TorsionSubgroup at v = %o: %o\n", v0, invm;
    assert invm eq G;
end for;
assert #Seqset(G2t) eq #G2t;
print "(5) these members are pairwise non-isomorphic (distinct G2-invariants)";
fc, s := Member(-1);
R<xx> := PolynomialRing(Q);
Ccensus := HyperellipticCurve(R![0,1,-4,-3,3,2], R![0,1,0,1]);
assert IsIsomorphic(HyperellipticCurve(fc), Ccensus);
invc := TorsionInvariants(fc);
assert invc eq [2,2,6];
printf "(5) control: the member v = -1 (not in the class) is Q-isomorphic to the census curve 2.2.6.b; J(Q)_tors = %o\n", invc;
printf "time %o s\n", Cputime(T0);
print "FAMILY_2_2_6_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E x F_v) and J(Q)_tors exactly [2,2,6]";
quit;
