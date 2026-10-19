/* family_6_6_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E x E1) with J(Q)_tors exactly [6,6].
   Found 2026-10-19 (torsion_jac notes/isopair_2026_10_19.md).

   Source.  The Howe-Leprevost-Poonen (2,2)-gluing (HLP 2000, Prop. 4) of the Hesse-form curve
       E : y^2 + 3 lam xy + (lam^3 - 1)y = x^3,  lam = (2t^3 + 1)/(3t^2)     (P = (0,0) of order 3),
   on the cover of the Hesse line where E has a rational point T of order 2 (x(T) = -mm^2, mm = (t^3 - 1)/(3t^2): the
   cover (lam + mm)^2 (lam - 2 mm) = 1 is rational, with parameter t = lam + mm), with E1 = E/C, C = mu_3 the subgroup
   over x = -(lam^2 + lam + 1), along the isomorphism E[2] -> E1[2] induced by the 3-isogeny phi : E -> E1.
   phi(P + T) has order 6, so E(Q) and E1(Q) both contain Z/6.  The gluing has coefficients in Q(t); after removing
   square factors and scaling x it is the curve f_t below, f_t = F(x^2) = lc (x^2 - r) q(x) q(-x) over Q(t), and
   phi1 : (x,y) -> (x^2, y) maps C_t onto Ea : Y^2 = F(X), phi2 : (x,y) -> (1/x^2, y/x^3) onto Eb : Y^2 = X^3 F(1/X);
   below, Ea is Q(t)-isomorphic to E and Eb to E1.
   Argument.
   (0) E, P, T, C = mu_3 (the two points of C - 0 are defined over Q(t)(sqrt -3) and swapped by its conjugation), E1, phi.
   (1) Splitting and gluing.  Ea = E and Eb = E1 over Q(t); phi1^* + phi2^* : Ea x Eb -> J is an isogeny over Q with
       kernel of order 4.  The gluing is along phi on E[2] (checked on the Weierstrass points).
   (2) The marked group over Q(t).  Da = phi1^*((P + T)') and Db = phi2^*(phi(P + T)') have order 6; the class
       c2 = [q, 0] of a quadratic factor q of f_t is the gain of the gluing (the Kummer maps commute with phi).  The 36
       classes i Da + j (2 Db) + k c2 (0 <= i < 6, 0 <= j < 3, 0 <= k < 2) are distinct, so G = [6,6] is marked over
       Q(t) and embeds in J_t(Q)_tors at every smooth member.
   (3) Exactness on a residue class.  For every 5-integral t = 3 (mod 5) the curve C_t has the reduction of C_3 at 5,
       with J(F_5) = [6,6] of order 36.  Reduction is injective on J_t(Q)_tors, so #J_t(Q)_tors divides 36 = #G and
       J_t(Q)_tors is exactly [6,6].
   (4) Infinitely many.  The G2-invariants of the members differ, so t -> M_2 is nonconstant with finite fibres.
   (5) Sanity: TorsionSubgroup at t = 3, -2, 8 is [6,6]; the member t = -2 is the census curve 6.6.b.
   Run (from the repository root):
       magma -b pipeline/magma/family_6_6_split_certificate.m > data/knowledge/sources/family_6_6_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
f := (24*t^11 + 108*t^10 + 162*t^9 + 198*t^8 + 108*t^7 + 81*t^6 + 18*t^5 + 27*t^4 + 3*t^2)*x^6 + (-48*t^12 - 360*t^11 - 1620*t^10 - 3822*t^9 - 2970*t^8 - 1620*t^7 - 1935*t^6 - 270*t^5 - 405*t^4 - 24*t^3 - 45*t^2 - 3)*x^4 + (96*t^12 + 1800*t^11 + 8100*t^10 + 14934*t^9 + 14850*t^8 + 8100*t^7 + 7515*t^6 + 1350*t^5 + 2025*t^4 + 48*t^3 + 225*t^2 + 6)*x^2 - 48*t^12 + 72*t^11 + 324*t^10 - 906*t^9 + 594*t^8 + 324*t^7 - 477*t^6 + 54*t^5 + 81*t^4 - 24*t^3 + 9*t^2 - 3;
assert &and[ Coefficient(f, i) eq 0 : i in [1, 3, 5] ] and Degree(f) eq 6 and IsSquarefree(f);
F := P![ Coefficient(f, 2*i) : i in [0..3] ];            // f = F(x^2)
printf "C_t : y^2 = f_t(x) = %o\n", f;
printf "    the irreducible factors of f_t over Q(t) have degrees %o\n", [ Degree(h[1]) : h in Factorization(f) ];

// ---------------------------------------------------------------- helpers
function MonicCubic(E)                       // the x-cubic of E after completing the square in y
    a1, a2, a3, a4, a6 := Explode(aInvariants(E));
    PP := PolynomialRing(BaseRing(E)); X := PP.1;
    return X^3 + (a2 + a1^2/4)*X^2 + (a4 + a1*a3/2)*X + (a6 + a3^2/4);
end function;
function MonEC(G)                            // a Weierstrass model of Y^2 = G(X), via X' = lc X, Y' = lc Y
    cs := [ Coefficient(G, i) : i in [0..3] ];
    return EllipticCurve([0, cs[3], 0, cs[2]*cs[4], cs[1]*cs[4]^2]), cs[4];
end function;
// pullbacks to J of points of Ea = MonEC(F) and Eb = MonEC(Fr), along phi1 = (x^2, y) and phi2 = (1/x^2, y/x^3)
function Pull1(JJ, Pt, c) PP := Parent((JJ!0)[1]); xx := PP.1; XX := Pt[1]/c; YY := Pt[2]/c; return elt< JJ | xx^2 - XX, PP!YY >; end function;
function Pull2(JJ, Pt, c) PP := Parent((JJ!0)[1]); xx := PP.1; XX := Pt[1]/c; YY := Pt[2]/c; return elt< JJ | xx^2 - 1/XX, (YY/XX)*xx >; end function;
// x(m(T)) for the point T of order 2 of the domain of m with x(T) = X, X in a residue ring Ft[x]/(g)
function XImage(m, X)
    a := aInvariants(Domain(m));
    vals := [ Evaluate(d, [X, -(a[1]*X + a[3])/2, Parent(X)!1]) : d in DefiningPolynomials(m) ];
    ok, inv := IsInvertible(vals[3]); assert ok;
    return vals[1]*inv;
end function;
// does the finite abelian group with invariants g embed in the one with invariants h (prime by prime)
function EmbedsIn(g, h)
    for l in &join[ Seqset(PrimeDivisors(n)) : n in g cat h | n gt 1 ] do
        a := Reverse(Sort([ Valuation(n, l) : n in g | n mod l eq 0 ]));
        b := Reverse(Sort([ Valuation(n, l) : n in h | n mod l eq 0 ]));
        if #a gt #b or &or[ a[i] gt b[i] : i in [1..#a] ] then return false; end if;
    end for;
    return true;
end function;

// ---------------------------------------------------------------- (0) source: E (Hesse form with a point of order 2), C = mu_3, E1 = E/C
lam := (2*t^3 + 1)/(3*t^2); mm := (t^3 - 1)/(3*t^2);        // lam^3 - 3 mm^2 lam - 2 mm^3 = 1
E := EllipticCurve([Ft| 3*lam, 0, lam^3 - 1, 0, 0]);
P3 := E![0, 0];
assert 3*P3 eq E!0 and P3 ne E!0;
T2 := E![-mm^2, -(3*lam*(-mm^2) + lam^3 - 1)/2];
assert T2 ne E!0 and 2*T2 eq E!0;
P6 := P3 + T2;
assert 6*P6 eq E!0 and 2*P6 ne E!0 and 3*P6 ne E!0;
x0 := -(lam^2 + lam + 1);
assert Evaluate(DivisionPolynomial(E, 3), x0) eq 0;
bb := 3*lam*x0 + lam^3 - 1; dsc := bb^2 + 4*x0^3;           // the points over x0: y^2 + bb y - x0^3 = 0
okmu, w0 := IsSquare(-3*dsc);
assert okmu and dsc ne 0;
E1, phi := IsogenyFromKernel(E, x - x0);
assert Degree(phi) eq 3;
dp := DefiningPolynomials(phi);
P6b := E1![ Evaluate(d, Eltseq(P6)) : d in dp ];
P3b := E1![ Evaluate(d, Eltseq(P3)) : d in dp ];
assert 6*P6b eq E1!0 and 2*P6b ne E1!0 and 3*P6b ne E1!0 and P3b eq 4*P6b;
print "(0) E : y^2 + 3 lam xy + (lam^3 - 1)y = x^3 with lam = (2t^3 + 1)/(3t^2): P = (0,0) of order 3 and T = (-mm^2, *) of order 2,";
print "    mm = (t^3 - 1)/(3t^2), so P6 = P + T has order 6; the 3-division polynomial vanishes at x0 = -(lam^2 + lam + 1),";
printf "    and -3 (disc of the y-quadratic at x0) = (%o)^2: C = <(x0, *)> is mu_3, E[3] = <P> + mu_3\n", w0;
print "    E1 = E/C, phi : E -> E1 of degree 3, and phi(P6) has order 6";

// ---------------------------------------------------------------- (1) the two elliptic quotients and the gluing
Fr := P![ Coefficient(F, 3 - i) : i in [0..3] ];
Ea, ca := MonEC(F); Eb, cb := MonEC(Fr);
oka, ma := IsIsomorphic(E, Ea);
okb, mb := IsIsomorphic(E1, Eb);
assert oka and okb;
print "(1) phi1 = (x^2, y) : C_t -> Ea : Y^2 = F(X), phi2 = (1/x^2, y/x^3) : C_t -> Eb : Y^2 = X^3 F(1/X);";
print "    Ea is Q(t)-isomorphic to E and Eb to E1, so J is isogenous over Q to E x E1 (phi1^* + phi2^* has kernel of order 4)";
g := MonicCubic(E);
A<th> := quo< P | g >;                        // th = x(T) for the generic point T of order 2 of E
al := XImage(phi, th);                        // x(phi(T))
okA, iA := IsInvertible(A!ca); okB, iB := IsInvertible(A!cb); assert okA and okB;
Ra := XImage(ma, th)*iA;                      // the root R of F with (R, 0) = image of T on Y^2 = F(X)
Rb := XImage(mb, al)*iB;                      // the root of Fr with (Rb, 0) = image of phi(T) on Y^2 = Fr(X)
assert Evaluate(F, Ra) eq 0 and Ra*Rb eq 1;
print "    for T in E[2] - 0 the two Weierstrass points x = +-sqrt(R) of C_t over the image (R, 0) of T in Ea lie over the";
print "    image (1/R, 0) of phi(T) in Eb (checked in Q(t)[x]/(2-division polynomial of E)): C_t is the (2,2)-gluing";
print "    of E and E1 along phi|E[2], and J = (E x E1)/graph(phi|E[2]) (HLP 2000, Prop. 4)";

// ---------------------------------------------------------------- (2) the marked group over Q(t)
J := Jacobian(HyperellipticCurve(f));
Da := Pull1(J, ma(P6), ca); Db := Pull2(J, mb(P6b), cb);
assert 6*Da eq J!0 and 2*Da ne J!0 and 3*Da ne J!0 and 6*Db eq J!0 and 2*Db ne J!0 and 3*Db ne J!0;
fq := Factorization(f);
assert #fq eq 3 and &and[ Degree(h[1]) eq 2 : h in fq ];
ev := [ h[1] : h in fq | Coefficient(h[1], 1) eq 0 ]; od := [ h[1] : h in fq | Coefficient(h[1], 1) ne 0 ];
assert #ev eq 1 and #od eq 2 and Evaluate(od[1], -x) eq LeadingCoefficient(od[1])/LeadingCoefficient(od[2])*od[2];
c2 := elt< J | od[1], 0 >;
assert c2 ne J!0 and 2*c2 eq J!0 and 3*Da eq elt< J | ev[1], 0 > and 3*Db eq 3*Da;
els := { i*Da + j*(2*Db) + k*c2 : i in [0..5], j in [0..2], k in [0..1] };
assert #els eq 36;
print "(2) over Q(t): Da = phi1^*((P+T)') and Db = phi2^*(phi(P+T)') have order 6, with 3 Da = 3 Db = [x^2 - r, 0];";
print "    c2 = [q, 0], q an odd quadratic factor of f_t, is the gain of the gluing (order 2); the 36 classes";
print "    i Da + j (2 Db) + k c2 are distinct, so G = <Da, 2 Db, c2> = Z/6 x Z/3 x Z/2 = [6,6] is marked over Q(t)";

// ---------------------------------------------------------------- (3) exactness on the residue class t = 3 mod 5
p := 5; r0 := 3;
members := [Q| 3, -2, 8, 13, -7, 4/3];
sanity := [Q| 3, -2, 8];
fbar := ReductionClass(f, p, r0);
for m in members do
    assert Valuation(m, p) ge 0 and GF(p)!m eq GF(p)!r0;
    assert PolynomialRing(GF(p))!SpecPoly(f, m) eq fbar;
end for;
n, chi := FrobData(fbar);
inv := Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
assert inv eq [6, 6] and n eq 36;
printf "(3) p = %o, t = %o mod %o: reduction %o (squarefree, degree %o); J(F_%o) = %o, chi = %o\n", p, r0, p, fbar, Degree(fbar), p, inv, chi;
print "    #J_t(Q)_tors divides #J(F_5) = 36 = #G, so every 5-integral t = 3 mod 5 has J_t(Q)_tors exactly [6,6]";

// ---------------------------------------------------------------- (4) nonconstant moduli map, (5) sanity
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } eq #members;
printf "(4) the G2-invariants of the members %o are pairwise distinct\n", members;
C0 := HyperellipticCurve(RX![1,3,6,7,6,3,1], RX![0,1,1]);
Cm := HyperellipticCurve(SpecPoly(f, -2));
assert IsIsomorphic(Cm, SimplifiedModel(C0));
print "    the member t = -2 is the census curve 6.6.b (isomorphic over Q)";
for m in sanity do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(5) TorsionSubgroup at t = %o: %o\n", m, invm;
    assert invm eq [6, 6];
end for;
printf "time %o s\n", Cputime(T0);
print "FAMILY_6_6_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E x E1) and J(Q)_tors exactly [6,6]";
quit;
