/* family_20_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E1 x E2) with J(Q)_tors exactly [20].
   Found 2026-10-19 (torsion_jac notes/hlpP2_2026_10_19.md).

   Source.  The Z/20 entry of Howe-Leprevost-Poonen 2000, Table 1 (section 3.5): the (2,2)-gluing (Prop. 4) of
   a curve with torsion Z/10 and a curve with a point of order 4 whose 2-torsion fields agree (Delta_10(t) y^2 =
   16u + 1), restricted to the line t = const: E1 is the fixed curve 66c1 = [1,0,0,-45,81] (torsion Z/10,
   P10 = (0,9)), whose non-rational 2-torsion is defined over Q(sqrt 33), and
       E2 = E_4(u) : y^2 + xy - uy = x^3 - ux^2,   u = (33 s^2 - 1)/16,   P4 = (0,0) of order 4,
   whose non-rational 2-torsion is defined over Q(s)(sqrt(16u + 1)) = Q(s)(sqrt 33).  The gluing is along the
   Galois isomorphism psi : E1[2] -> E2[2] sending the rational point to the rational point and the root
   (-b1 + sqrt D1)/2 of the quadratic factor of the 2-division cubic of E1 to (-b2 + sqrt D2)/2 for E2, where
   sqrt D1 = sqrt(D1/33) sqrt 33, sqrt D2 = sqrt(D2/33) sqrt 33 with sqrt(D2/33) in Q(s) (the other psi is
   s -> -s).  The HLP sextic is computed over Q(s)(sqrt 33); its coefficients lie in Q(s), and
   x -> lambda x, y -> mu y (lambda, mu in Q(s)) takes it to the model
       C_s : y^2 = f_s(x)   (embedded below; f_s = F(x^2), coefficients of degree <= 7 in s).
   phi1 : (x,y) -> (x^2, y) maps C_s onto Ea : Y^2 = F(X), phi2 : (x,y) -> (1/x^2, y/x^3) onto
   Eb : Y^2 = X^3 F(1/X); below, Ea is Q(s)-isomorphic to E2 and Eb to E1.

   Argument.
   (1) Splitting.  Wherever f_s is squarefree of degree 6, phi1^* + phi2^* : Ea x Eb -> J is an isogeny over Q
       (the involutions (x,y) -> (-x,y) and (x,y) -> (-x,-y) compose to the hyperelliptic involution), so J is
       isogenous over Q to E1 x E2.
   (2) The marked group over Q(s).  D10 = phi2^*(image of P10) has order 10, D4 = phi1^*(image of P4) has
       order 4, and 5 D10 = 2 D4 (the rational 2-torsion points are glued); D20 = D10 + D4 has order 20
       (Mumford arithmetic over Q(s)).  So Z/20 embeds in J_s(Q)_tors at every smooth member.
   (3) Exactness on a residue class.  #J(F_p) = #E1(F_p) #E2(F_p) is always divisible by 40, so the group
       structure is used.  For every 5- and 19-integral s = 2 (mod 95), C_s has the reductions of C_2 at 5 and
       19, with J(F_5) = Z/80 (cyclic) and J(F_19) = [20,20] (exponent 20).  J_s(Q)_tors embeds in both (p odd,
       good reduction), so it is cyclic of order dividing 20: it is exactly [20].  (StructureBound: for every
       prime l the sorted l-primary exponents of a subgroup are bounded componentwise by those of the group.)
   (4) Infinitely many.  The G2-invariants of the members differ, so s -> M_2 is nonconstant with finite
       fibres, and the residue class gives infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup at s = 2, 97, -93 is [20].  (The member s = 1, outside the class, has
       J(Q)_tors = [40]: a 2-power gain of the gluing; see the note.)
   Run (from the repository root):
       magma -b pipeline/magma/family_20_split_certificate.m > data/knowledge/sources/family_20_split_certificate.log */
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
f := (1293732*s^7 - 217800*s^6 - 39204*s^5 + 13200*s^4 - 1188*s^3 - 200*s^2 + 36*s)*x^6
   + (-2910897*s^6 - 490050*s^5 + 958881*s^4 - 207900*s^3 + 29057*s^2 - 450*s - 81)*x^4
   + (2205225*s^4 - 2567268*s^3 + 503350*s^2 - 77796*s + 2025)*x^2 + 1368576*s^2 - 230400*s + 41472;
assert &and[ Coefficient(f, i) eq 0 : i in [1, 3, 5] ] and Degree(f) eq 6 and IsSquarefree(f);
printf "C_s : y^2 = f_s(x) = %o\n", f;

// ---------------------------------------------------------------- (0) source: the HLP gluing of 66c1 and E_4(u), u = (33 s^2 - 1)/16
E1 := EllipticCurve([Fs| 1, 0, 0, -45, 81]);
assert Conductor(EllipticCurve([Q| 1, 0, 0, -45, 81])) eq 66;
P10 := E1![0, 9];
assert 10*P10 eq E1!0 and 5*P10 ne E1!0 and 2*P10 ne E1!0;
u := (33*s^2 - 1)/16;
E2 := EllipticCurve([Fs| 1, -u, -u, 0, 0]);
P4 := E2![0, 0];
assert 4*P4 eq E2!0 and 2*P4 ne E2!0;
g1 := MonicCubic(E1, x); g2 := MonicCubic(E2, x);
lin1 := [ g[1] : g in Factorization(g1) | Degree(g[1]) eq 1 ]; q1 := [ g[1] : g in Factorization(g1) | Degree(g[1]) eq 2 ];
lin2 := [ g[1] : g in Factorization(g2) | Degree(g[1]) eq 1 ]; q2 := [ g[1] : g in Factorization(g2) | Degree(g[1]) eq 2 ];
assert #lin1 eq 1 and #q1 eq 1 and #lin2 eq 1 and #q2 eq 1;
ok1, w1 := IsSquare(Discriminant(q1[1])/33); ok2, w2 := IsSquare(Discriminant(q2[1])/33);
assert ok1 and ok2 and not IsSquare(Discriminant(q1[1]));
K<z> := ext< Fs | x^2 - 33 >;
al := [ K!(-Coefficient(lin1[1], 0)), (-Coefficient(q1[1], 1) + w1*z)/2, (-Coefficient(q1[1], 1) - w1*z)/2 ];
be := [ K!(-Coefficient(lin2[1], 0)), (-Coefficient(q2[1], 1) + w2*z)/2, (-Coefficient(q2[1], 1) - w2*z)/2 ];
assert &and[ Evaluate(PolynomialRing(K)!g1, a) eq 0 : a in al ] and &and[ Evaluate(PolynomialRing(K)!g2, b) eq 0 : b in be ];
hK := HLPSextic(al, be);
assert &and[ IsCoercible(Fs, c) : c in Coefficients(hK) ];
h := P![ Fs!c : c in Coefficients(hK) ];
assert SameUpToScaling(h, f);
print "(0) the HLP (Prop. 4) gluing of E1 = 66c1 and E2 = E_4((33s^2 - 1)/16) along psi has coefficients in Q(s),";
print "    and f_s = h(lambda x)/mu^2 with lambda, mu in Q(s)";

// ---------------------------------------------------------------- (1) the two elliptic quotients
F := P![ Coefficient(f, 2*i) : i in [0..3] ];
Fr := P![ Coefficient(F, 3 - i) : i in [0..3] ];
Ea, ca := MonEC(F); Eb, cb := MonEC(Fr);
oka, ma := IsIsomorphic(E2, Ea);
okb, mb := IsIsomorphic(E1, Eb);
assert oka and okb;
print "(1) phi1 = (x^2, y) : C_s -> Ea : Y^2 = F(X), phi2 = (1/x^2, y/x^3) : C_s -> Eb : Y^2 = X^3 F(1/X);";
print "    Ea is Q(s)-isomorphic to E2 = E_4(u) and Eb to E1 = 66c1: J is isogenous over Q to E1 x E2";

// ---------------------------------------------------------------- (2) the marked group over Q(s)
J := Jacobian(HyperellipticCurve(f));
D10 := Pull2(J, mb(P10), cb);
D4 := Pull1(J, ma(P4), ca);
assert 10*D10 eq J!0 and 5*D10 ne J!0 and 2*D10 ne J!0;
assert 4*D4 eq J!0 and 2*D4 ne J!0 and 5*D10 eq 2*D4;
D20 := D10 + D4;
assert 20*D20 eq J!0 and 10*D20 ne J!0 and 4*D20 ne J!0;
print "(2) over Q(s): D10 = phi2^*(P10) has order 10, D4 = phi1^*(P4) has order 4, 5 D10 = 2 D4, and D20 = D10 + D4 has order 20";

// ---------------------------------------------------------------- (3) exactness on the residue class s = 2 mod 95
primes := [5, 19]; s0 := 2; modulus := 95;
members := [Q| 2, 97, -93, 192, 287];
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
assert invs eq [[80], [20, 20]];
bnd := StructureBound(invs);
printf "    the largest group embedding in all J(F_p) has invariants %o\n", bnd;
assert bnd eq [20];
printf "    hence every 5- and 19-integral s = %o mod %o has J_s(Q)_tors exactly [20]\n", s0, modulus;

// ---------------------------------------------------------------- (4) nonconstant moduli map, (5) sanity
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } eq #members;
printf "(4) the G2-invariants of the members %o are pairwise distinct\n", members;
for m in [Q| 2, 97, -93] do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(5) TorsionSubgroup at s = %o: %o\n", m, invm;
    assert invm eq [20];
end for;
inv1 := TorsionInvariants(SpecPoly(f, 1));
printf "(5) outside the class: TorsionSubgroup at s = 1: %o\n", inv1;
assert inv1 eq [40];
printf "time %o s\n", Cputime(T0);
print "FAMILY_20_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E1 x E2) and J(Q)_tors exactly [20]";
quit;
