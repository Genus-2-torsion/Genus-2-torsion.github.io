/* family_3_3_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E x E1) with J(Q)_tors exactly [3,3].
   Found 2026-10-19 (torsion_jac notes/isopair_2026_10_19.md).

   Source.  The Howe-Leprevost-Poonen (2,2)-gluing (HLP 2000, Prop. 4) of the Hesse-form curve
       E = E_t : y^2 + 3t xy + (t^3 - 1)y = x^3     (P = (0,0) of order 3, a1^3 - 27 a3 = 27),
   whose 3-torsion is E[3] = <P> + C with C = mu_3 the subgroup over x = -(t^2 + t + 1), with E1 = E/C, along the
   isomorphism E[2] -> E1[2] induced by the 3-isogeny phi : E -> E1.  (The Hesse line parametrizes curves with a
   rational point of order 3 and an independent rational 3-isogeny; Gamma1(3) cap Gamma0(9), conjugate to
   Gamma1(3) cap Gamma^0(3), has genus 0.)  phi(P) has order 3, so E(Q) and E1(Q) both have a point of order 3.  The gluing has coefficients in Q(t); after removing square factors and
   scaling x it is the curve f_t below, f_t = F(x^2) with F a cubic over Q(t), and
   phi1 : (x,y) -> (x^2, y) maps C_t onto Ea : Y^2 = F(X), phi2 : (x,y) -> (1/x^2, y/x^3) onto Eb : Y^2 = X^3 F(1/X);
   below, Ea is Q(t)-isomorphic to E and Eb to E1.
   Argument.
   (0) E, P, C = mu_3 (the two points of C - 0 are defined over Q(t)(sqrt -3) and swapped by its conjugation), E1, phi.
   (1) Splitting and gluing.  Ea = E and Eb = E1 over Q(t); phi1^* + phi2^* : Ea x Eb -> J is an isogeny over Q with
       kernel of order 4.  The Weierstrass points over the image of T in E[2] lie over the image of phi(T): the gluing is
       along phi on E[2].
   (2) The marked group over Q(t).  Da = phi1^*(P') and Db = phi2^*(phi(P)') (primes: images on Ea, Eb) have order 3
       and Db is not in <Da>, so G = <Da, Db> = [3,3] is marked over Q(t) and embeds in J_t(Q)_tors at every smooth member.
   (3) Exactness on a residue class.  For every 5-integral t = 2 (mod 5) the curve C_t has the reduction of C_2 at 5,
       with J(F_5) = [3,3] of order 9.  Reduction is injective on J_t(Q)_tors (p = 5 odd, good reduction), so
       #J_t(Q)_tors divides 9 = #G and J_t(Q)_tors is exactly [3,3].
   (4) Infinitely many.  The G2-invariants of the members differ, so t -> M_2 is nonconstant with finite fibres.
   (5) Sanity: TorsionSubgroup at t = 2, 7, -3 is [3,3]; the member t = -8 is the census curve 3.3.b.
   Run (from the repository root):
       magma -b pipeline/magma/family_3_3_split_certificate.m > data/knowledge/sources/family_3_3_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
f := (t^3 + 3*t^2 + 3*t + 2)*x^6 + (-3*t^4 - 15*t^3 - 45*t^2 - 69*t - 30)*x^4 + (6*t^4 + 75*t^3 + 225*t^2 + 273*t + 150)*x^2 - 3*t^4 + 3*t^3 + 9*t^2 - 15*t + 6;
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

// ---------------------------------------------------------------- (0) source: E (Hesse form), C = mu_3, E1 = E/C
E := EllipticCurve([Ft| 3*t, 0, t^3 - 1, 0, 0]);
P3 := E![0, 0];
assert 3*P3 eq E!0 and P3 ne E!0;
x0 := -(t^2 + t + 1);
assert Evaluate(DivisionPolynomial(E, 3), x0) eq 0;
bb := 3*t*x0 + t^3 - 1; dsc := bb^2 + 4*x0^3;               // the points over x0: y^2 + bb y - x0^3 = 0
okmu, w0 := IsSquare(-3*dsc);
assert okmu and dsc ne 0;
E1, phi := IsogenyFromKernel(E, x - x0);
assert Degree(phi) eq 3;
dp := DefiningPolynomials(phi);
P3b := E1![ Evaluate(d, Eltseq(P3)) : d in dp ];
assert P3b ne E1!0 and 3*P3b eq E1!0;
print "(0) E : y^2 + 3t xy + (t^3 - 1)y = x^3, P = (0,0) of order 3; the 3-division polynomial vanishes at x0 = -(t^2+t+1)";
printf "    and -3 (disc of the y-quadratic at x0) = (%o)^2: C = <(x0, *)> is mu_3, E[3] = <P> + mu_3\n", w0;
print "    E1 = E/C, phi : E -> E1 of degree 3, and phi(P) has order 3";

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
Da := Pull1(J, ma(P3), ca); Db := Pull2(J, mb(P3b), cb);
assert Da ne J!0 and 3*Da eq J!0 and Db ne J!0 and 3*Db eq J!0;
assert Db notin [J!0, Da, 2*Da];
print "(2) over Q(t): Da = phi1^*(P') and Db = phi2^*(phi(P)') have order 3 and Db is not in <Da>,";
print "    so G = <Da, Db> = [3,3] is marked over Q(t)";

// ---------------------------------------------------------------- (3) exactness on the residue class t = 2 mod 5
p := 5; r0 := 2;
members := [Q| 2, 7, -3, 12, 7/6, -8];
sanity := [Q| 2, 7, -3];
fbar := ReductionClass(f, p, r0);
for m in members do
    assert Valuation(m, p) ge 0 and GF(p)!m eq GF(p)!r0;
    assert PolynomialRing(GF(p))!SpecPoly(f, m) eq fbar;
end for;
n, chi := FrobData(fbar);
inv := Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
assert inv eq [3, 3] and n eq 9;
printf "(3) p = %o, t = %o mod %o: reduction %o (squarefree, degree %o); J(F_%o) = %o, chi = %o\n", p, r0, p, fbar, Degree(fbar), p, inv, chi;
print "    #J_t(Q)_tors divides #J(F_5) = 9 = #G, so every 5-integral t = 2 mod 5 has J_t(Q)_tors exactly [3,3]";

// ---------------------------------------------------------------- (4) nonconstant moduli map, (5) sanity
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } eq #members;
printf "(4) the G2-invariants of the members %o are pairwise distinct\n", members;
C0 := HyperellipticCurve(RX![-2,6,-8,6,-8,6,-2], RX![0]);
Cm := HyperellipticCurve(SpecPoly(f, -8));
assert IsIsomorphic(Cm, SimplifiedModel(C0));
print "    the member t = -8 is the census curve 3.3.b (isomorphic over Q)";
for m in sanity do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(5) TorsionSubgroup at t = %o: %o\n", m, invm;
    assert invm eq [3, 3];
end for;
printf "time %o s\n", Cputime(T0);
print "FAMILY_3_3_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E x E1) and J(Q)_tors exactly [3,3]";
quit;
