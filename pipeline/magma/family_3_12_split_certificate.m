/* family_3_12_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E1 x E2) with J(Q)_tors exactly [3,12].
   Found 2026-10-19 (torsion_jac notes/hlpP2_2026_10_19.md).

   Source.  The Z/3 x Z/12 entry of Howe-Leprevost-Poonen 2000, Table 1 (section 3.5): the (2,2)-gluing
   (Prop. 4) of a curve with torsion Z/12 and a curve with torsion Z/6 whose 2-torsion fields agree
   (Delta_12(t) w^2 = (u + 1)(9u + 1), a conic over Q(t) with the point u = -1), restricted to the line
   t = const: E1 is the fixed curve 90c3 = [1,-1,1,-122,1721] (torsion Z/12, P12 = (-9,-41)), whose
   non-rational 2-torsion is defined over Q(sqrt -15), and E2 = E_6(u) is the Kubert curve
       E_6(u) : y^2 + (1 - u)xy - (u^2 + u)y = x^3 - (u^2 + u)x^2,   P6 = (0,0) of order 6,
   with u = (-15 s^2 - 1)/(9 + 15 s^2) (the line w = s(u + 1) through u = -1), so that (u + 1)(9u + 1) is
   -15 times a square and the non-rational 2-torsion of E2 is defined over Q(s)(sqrt -15).  The gluing is
   along the Galois isomorphism psi : E1[2] -> E2[2] sending the rational point to the rational point and
   (-b1 + sqrt D1)/2 to (-b2 + sqrt D2)/2, with sqrt D_i = sqrt(D_i/-15) sqrt -15 (the other psi is s -> -s).
   The HLP sextic is computed over Q(s)(sqrt -15); its coefficients lie in Q(s), and x -> lambda x, y -> mu y
   (lambda, mu in Q(s)) takes it to the model
       C_s : y^2 = f_s(x)   (embedded below; f_s = F(x^2), coefficients of degree <= 13 in s).
   phi1 : (x,y) -> (x^2, y) maps C_s onto Ea : Y^2 = F(X), phi2 : (x,y) -> (1/x^2, y/x^3) onto
   Eb : Y^2 = X^3 F(1/X); below, Ea is Q(s)-isomorphic to E2 and Eb to E1.

   Argument.
   (1) Splitting.  Wherever f_s is squarefree of degree 6, phi1^* + phi2^* : Ea x Eb -> J is an isogeny over Q
       (the involutions (x,y) -> (-x,y) and (x,y) -> (-x,-y) compose to the hyperelliptic involution), so J is
       isogenous over Q to E1 x E2.
   (2) The marked group over Q(s).  D12 = phi2^*(image of P12) has order 12, D6 = phi1^*(image of P6) has
       order 6, 6 D12 = 3 D6 (the rational 2-torsion points are glued), and 2 D6 is not in <D12> (Mumford
       arithmetic over Q(s)).  So <D12, D6> = [3,12] embeds in J_s(Q)_tors at every smooth member.
   (3) Exactness on a residue class.  #J(F_p) = #E1(F_p) #E2(F_p) is always divisible by 72, so the group
       structure is used.  For every 7- and 11-integral s = -9 (mod 77), C_s has the reductions of C_{-9} at 7
       and 11, with J(F_7) = [12,12] (2-part of exponent 4, 3-part of exponent 3) and J(F_11) = [3,24] (2-part
       cyclic).  J_s(Q)_tors embeds in both (p odd, good reduction), so it embeds in [3,12]: it is exactly
       [3,12].  (StructureBound: for every prime l the sorted l-primary exponents of a subgroup are bounded
       componentwise by those of the group.)
   (4) Infinitely many.  The G2-invariants of the members differ, so s -> M_2 is nonconstant with finite
       fibres, and the residue class gives infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup at s = -9, 68, -86 is [3,12].  The member s = 1/3 (outside the class) is the census
       curve 3.12.a.  The member s = -1, outside the class, is the census curve 3.24.a, with J(Q)_tors = [3,24]: a
       2-power gain of the gluing (see family_3_24_split_certificate.m).
   Run (from the repository root):
       magma -b pipeline/magma/family_3_12_split_certificate.m > data/knowledge/sources/family_3_12_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
SetQuitOnError(true);                   // fail closed: a failed assertion stops the run before the marker
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

function MonEC(G)                       // Y^2 = G(X), G a cubic, as a Weierstrass curve (X' = lc X, Y' = lc Y)
    cs := [ Coefficient(G, i) : i in [0..3] ];
    return EllipticCurve([0, cs[3], 0, cs[2]*cs[4], cs[1]*cs[4]^2]), cs[4];
end function;
function MonicCubic(E, x)               // y^2 = MonicCubic(x) after completing the square
    a1, a2, a3, a4, a6 := Explode(aInvariants(E));
    return x^3 + (a2 + a1^2/4)*x^2 + (a4 + a1*a3/2)*x + (a6 + a3^2/4);
end function;
function HLPSextic(al, be)              // HLP 2000, Prop. 4, with psi(al_i) = be_i
    L := Parent(al[1]); PL<X> := PolynomialRing(L);
    aa1 := (al[3]-al[2])^2/(be[3]-be[2]) + (al[2]-al[1])^2/(be[2]-be[1]) + (al[1]-al[3])^2/(be[1]-be[3]);
    bb1 := (be[3]-be[2])^2/(al[3]-al[2]) + (be[2]-be[1])^2/(al[2]-al[1]) + (be[1]-be[3])^2/(al[1]-al[3]);
    aa2 := al[1]*(be[3]-be[2]) + al[2]*(be[1]-be[3]) + al[3]*(be[2]-be[1]);
    bb2 := be[1]*(al[3]-al[2]) + be[2]*(al[1]-al[3]) + be[3]*(al[2]-al[1]);
    Df := ((al[1]-al[2])*(al[1]-al[3])*(al[2]-al[3]))^2;
    Dg := ((be[1]-be[2])*(be[1]-be[3])*(be[2]-be[3]))^2;
    AA := Dg*aa1/aa2; BB := Df*bb1/bb2;
    return -(AA*(al[2]-al[1])*(al[1]-al[3])*X^2 + BB*(be[2]-be[1])*(be[1]-be[3]))
           *(AA*(al[3]-al[2])*(al[2]-al[1])*X^2 + BB*(be[3]-be[2])*(be[2]-be[1]))
           *(AA*(al[1]-al[3])*(al[3]-al[2])*X^2 + BB*(be[1]-be[3])*(be[3]-be[2]));
end function;
// pullbacks of points (X, Y) of the Weierstrass models (scaled by c) to the Jacobian of y^2 = F(x^2)
function Pull1(JJ, Pt, c) PP := Parent((JJ!0)[1]); xx := PP.1; XX := Pt[1]/c; YY := Pt[2]/c; return elt< JJ | xx^2 - XX, PP!YY >; end function;
function Pull2(JJ, Pt, c) PP := Parent((JJ!0)[1]); xx := PP.1; XX := Pt[1]/c; YY := Pt[2]/c; return elt< JJ | xx^2 - 1/XX, (YY/XX)*xx >; end function;
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
// f2 = h(lam x)/mu^2 for some lam, mu in the base field (h, f2 even sextics)
function SameUpToScaling(h, f2)
    okl, lam := IsSquare(Coefficient(f2, 2)*Coefficient(h, 0)/(Coefficient(f2, 0)*Coefficient(h, 2)));
    okm, mu := IsSquare(Coefficient(h, 0)/Coefficient(f2, 0));
    return okl and okm and &and[ Coefficient(f2, i)*mu^2 eq Coefficient(h, i)*lam^i : i in [0..6] ];
end function;

Fs<s> := FunctionField(Q); P<x> := PolynomialRing(Fs);
f := (-303750000*s^13 - 364500000*s^11 + 131760000*s^10 - 133650000*s^9 + 105408000*s^8 - 15480000*s^7 + 17568000*s^6 - 522000*s^5 + 1093120*s^4 + 12000*s^3 + 23424*s^2 + 720*s)*x^6 + (284765625*s^12 + 341718750*s^10 + 123525000*s^9 + 125296875*s^8 + 98820000*s^7 - 81852300*s^6 + 16470000*s^5 - 38056545*s^4 - 8344800*s^3 + 1273614*s^2 + 21960*s - 675)*x^4 + (-231609375*s^8 - 185287500*s^6 + 35046000*s^5 - 30881250*s^4 + 14018400*s^3 + 30848676*s^2 - 467280*s - 41175)*x^2 + 51840000*s^4 + 20736000*s^2 - 22487040*s - 691200;
assert &and[ Coefficient(f, i) eq 0 : i in [1, 3, 5] ] and Degree(f) eq 6 and IsSquarefree(f);
printf "C_s : y^2 = f_s(x) = %o\n", f;

// ---------------------------------------------------------------- (0) source: the HLP gluing of 90c3 and E_6(u), u = (-15 s^2 - 1)/(9 + 15 s^2)
E1 := EllipticCurve([Fs| 1, -1, 1, -122, 1721]);
assert Conductor(EllipticCurve([Q| 1, -1, 1, -122, 1721])) eq 90;
P12 := E1![-9, -41];
assert 12*P12 eq E1!0 and &and[ (12 div l)*P12 ne E1!0 : l in PrimeDivisors(12) ];
u := (-15*s^2 - 1)/(9 + 15*s^2);
E2 := EllipticCurve([Fs| 1 - u, -(u^2 + u), -(u^2 + u), 0, 0]);
P6 := E2![0, 0];
assert 6*P6 eq E2!0 and 2*P6 ne E2!0 and 3*P6 ne E2!0;
g1 := MonicCubic(E1, x); g2 := MonicCubic(E2, x);
lin1 := [ g[1] : g in Factorization(g1) | Degree(g[1]) eq 1 ]; q1 := [ g[1] : g in Factorization(g1) | Degree(g[1]) eq 2 ];
lin2 := [ g[1] : g in Factorization(g2) | Degree(g[1]) eq 1 ]; q2 := [ g[1] : g in Factorization(g2) | Degree(g[1]) eq 2 ];
assert #lin1 eq 1 and #q1 eq 1 and #lin2 eq 1 and #q2 eq 1;
ok1, w1 := IsSquare(Discriminant(q1[1])/(-15)); ok2, w2 := IsSquare(Discriminant(q2[1])/(-15));
assert ok1 and ok2 and not IsSquare(Discriminant(q1[1]));
K<z> := ext< Fs | x^2 - (-15) >;
al := [ K!(-Coefficient(lin1[1], 0)), (-Coefficient(q1[1], 1) + w1*z)/2, (-Coefficient(q1[1], 1) - w1*z)/2 ];
be := [ K!(-Coefficient(lin2[1], 0)), (-Coefficient(q2[1], 1) + w2*z)/2, (-Coefficient(q2[1], 1) - w2*z)/2 ];
assert &and[ Evaluate(PolynomialRing(K)!g1, a) eq 0 : a in al ] and &and[ Evaluate(PolynomialRing(K)!g2, b) eq 0 : b in be ];
hK := HLPSextic(al, be);
assert &and[ IsCoercible(Fs, c) : c in Coefficients(hK) ];
h := P![ Fs!c : c in Coefficients(hK) ];
assert SameUpToScaling(h, f);
print "(0) the HLP (Prop. 4) gluing of E1 = 90c3 and E2 = E_6((-15 s^2 - 1)/(9 + 15 s^2)) along psi has coefficients in Q(s),";
print "    and f_s = h(lambda x)/mu^2 with lambda, mu in Q(s)";

// ---------------------------------------------------------------- (1) the two elliptic quotients
F := P![ Coefficient(f, 2*i) : i in [0..3] ];
Fr := P![ Coefficient(F, 3 - i) : i in [0..3] ];
Ea, ca := MonEC(F); Eb, cb := MonEC(Fr);
oka, ma := IsIsomorphic(E2, Ea);
okb, mb := IsIsomorphic(E1, Eb);
assert oka and okb;
print "(1) phi1 = (x^2, y) : C_s -> Ea : Y^2 = F(X), phi2 = (1/x^2, y/x^3) : C_s -> Eb : Y^2 = X^3 F(1/X);";
print "    Ea is Q(s)-isomorphic to E2 = E_6(u) and Eb to E1 = 90c3: J is isogenous over Q to E1 x E2";

// ---------------------------------------------------------------- (2) the marked group over Q(s)
J := Jacobian(HyperellipticCurve(f));
D12 := Pull2(J, mb(P12), cb);
D6 := Pull1(J, ma(P6), ca);
assert 12*D12 eq J!0 and 6*D12 ne J!0 and 4*D12 ne J!0;
assert 6*D6 eq J!0 and 3*D6 ne J!0 and 2*D6 ne J!0 and 6*D12 eq 3*D6;
assert 2*D6 notin { 4*D12, 8*D12 };          // the two elements of order 3 of <D12>
print "    <D12, D6> has order 12 * 6 / 2 = 36, 3-rank 2 and a cyclic 4-part: it is [3,12]";
print "(2) over Q(s): D12 = phi2^*(P12) has order 12, D6 = phi1^*(P6) has order 6, 6 D12 = 3 D6, 2 D6 not in <D12>";

// ---------------------------------------------------------------- (3) exactness on the residue class s = -9 mod 77
primes := [7, 11]; s0 := -9; modulus := 77;
members := [Q| -9, 68, -86, 145, 222];
invs := [];
for p in primes do
    r0 := Integers()!(GF(p)!s0);
    fbar := ReductionClass(f, p, r0);
    for m in members do
        assert GF(p)!m eq GF(p)!r0;
        assert PolynomialRing(GF(p))!SpecPoly(f, m) eq fbar;
    end for;
    n, chi := FrobData(fbar);
    inv := Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
    assert &*inv eq n;
    Append(~invs, inv);
    printf "(3) p = %o, s = %o mod %o: reduction %o (squarefree, degree %o); J(F_%o) = %o, chi = %o\n", p, r0, p, fbar, Degree(fbar), p, inv, chi;
end for;
assert invs eq [[12, 12], [3, 24]];
bnd := StructureBound(invs);
printf "    the largest group embedding in all J(F_p) has invariants %o\n", bnd;
assert bnd eq [3, 12];
printf "    hence every s = %o mod %o, integral at %o, has J_s(Q)_tors exactly [3,12]\n", s0, modulus, primes;

// ---------------------------------------------------------------- (4) nonconstant moduli map, (5) sanity
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } eq #members;
printf "(4) the G2-invariants of the members %o are pairwise distinct\n", members;
for m in [Q| -9, 68, -86] do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(5) TorsionSubgroup at s = %o: %o\n", m, invm;
    assert invm eq [3, 12];
end for;
inv1 := TorsionInvariants(SpecPoly(f, -1));
printf "(5) outside the class: TorsionSubgroup at s = -1: %o\n", inv1;
assert inv1 eq [3, 24];
Ccen := HyperellipticCurve(PolynomialRing(Q)![-8, 12, 6, -11, -3, 3, 1], PolynomialRing(Q)![1, 1, 1]);
assert IsIsomorphic(HyperellipticCurve(SpecPoly(f, 1/3)), Ccen);
print "(5) the member s = 1/3 (outside the class) is Q-isomorphic to the census curve 3.12.a";
printf "time %o s\n", Cputime(T0);
print "FAMILY_3_12_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E1 x E2) and J(Q)_tors exactly [3,12]";
quit;
