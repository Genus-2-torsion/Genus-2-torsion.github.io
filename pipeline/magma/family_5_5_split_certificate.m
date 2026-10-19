/* family_5_5_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E x E1) with J(Q)_tors exactly [5,5].
   Found 2026-10-19 (torsion_jac notes/isopair_2026_10_19.md).

   Source.  The Howe-Leprevost-Poonen (2,2)-gluing (HLP 2000, Prop. 4) of the X_1(5) curve in Tate normal form at b = c = t^5,
       E : y^2 + (1 - t^5)xy - t^5 y = x^3 - t^5 x^2     (P = (0,0) of order 5),
   with E1 = E/<P>, along the isomorphism E[2] -> E1[2] induced by the 5-isogeny phi : E -> E1.  The kernel of the dual
   isogeny is phi(E[5]) = mu_5 (Weil pairing), and on the cover b = t^5 of X_1(5) the curve E1 has a rational point Q1
   of order 5 (given explicitly below), so E1[5] = <Q1> + mu_5 and E = E1/mu_5.  (Gamma1(5) cap Gamma0(25), conjugate to
   Gamma1(5) cap Gamma^0(5), has genus 0.)  The member t = 1 is the census curve 5.5.a; its factors are X_1(11) and
   X_0(11), as in HLP's Corollary 5.  The gluing has coefficients in Q(t); after removing square factors and scaling x it is the curve f_t below,
   f_t = F(x^2) with F a cubic over Q(t), and phi1 : (x,y) -> (x^2, y) maps C_t onto Ea : Y^2 = F(X),
   phi2 : (x,y) -> (1/x^2, y/x^3) onto Eb : Y^2 = X^3 F(1/X); below, Ea is Q(t)-isomorphic to E and Eb to E1.
   Argument.
   (0) E, P, E1 = E/<P>, phi, and the point Q1 of order 5 of E1.
   (1) Splitting and gluing.  Ea = E and Eb = E1 over Q(t); phi1^* + phi2^* : Ea x Eb -> J is an isogeny over Q with
       kernel of order 4.  The gluing is along phi on E[2] (checked on the Weierstrass points).
   (2) The marked group over Q(t).  Da = phi1^*(P') and Db = phi2^*(Q1') have order 5 and Db is not in <Da>, so
       G = <Da, Db> = [5,5] is marked over Q(t) and embeds in J_t(Q)_tors at every smooth member.
   (3) Exactness on a residue class.  For every 3-integral t = 1 (mod 3) the curve C_t has the reduction of C_1 at 3,
       with J(F_3) = [5,5] of order 25.  Reduction is injective on J_t(Q)_tors (p = 3 odd, good reduction), so
       #J_t(Q)_tors divides 25 = #G and J_t(Q)_tors is exactly [5,5].
   (4) Infinitely many.  The G2-invariants of the members differ, so t -> M_2 is nonconstant with finite fibres.
   (5) Sanity: TorsionSubgroup at t = 1, 4, -2 is [5,5]; the member t = 1 is the census curve 5.5.a.
   Run (from the repository root):
       magma -b pipeline/magma/family_5_5_split_certificate.m > data/knowledge/sources/family_5_5_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
f := (-t^20 - t^10)*x^6 + (t^30 - 6*t^25 - 6*t^20 - 6*t^10 + 6*t^5 + 1)*x^4 + (-2*t^30 - 36*t^25 + 15*t^20 + 15*t^10 + 36*t^5 - 2)*x^2 + t^30 - 22*t^25 + 120*t^20 + 120*t^10 + 22*t^5 + 1;
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

// ---------------------------------------------------------------- (0) source: E = E_5(t^5), E1 = E/<P>, Q1 in E1(Q(t))
s := t^5;
E := EllipticCurve([Ft| 1 - s, -s, -s, 0, 0]);
P5 := E![0, 0];
assert 5*P5 eq E!0 and P5 ne E!0;
E1, phi := IsogenyFromKernel(E, (x - P5[1])*(x - (2*P5)[1]));
assert Degree(phi) eq 5;
dp := DefiningPolynomials(phi);
assert E1![ Evaluate(d, Eltseq(P5)) : d in dp ] eq E1!0;
Q1 := E1![ t^8 + t^7 + 2*t^6 - 2*t^5 + 5*t^4 - 3*t^3 + 2*t^2 - t,
           t^13 + 3*t^11 - t^10 + 5*t^9 - 5*t^8 + 11*t^7 - 16*t^6 + 14*t^5 - 10*t^4 + 6*t^3 - 3*t^2 + t ];
assert 5*Q1 eq E1!0 and Q1 ne E1!0;
printf "(0) E : y^2 + (1 - t^5)xy - t^5 y = x^3 - t^5 x^2, P = (0,0) of order 5; E1 = E/<P> = %o (Velu)\n", aInvariants(E1);
printf "    E1 has the point Q1 = %o of order 5 over Q(t)\n", Q1;

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
Da := Pull1(J, ma(P5), ca); Db := Pull2(J, mb(Q1), cb);
assert Da ne J!0 and 5*Da eq J!0 and Db ne J!0 and 5*Db eq J!0;
assert Db notin [ i*Da : i in [0..4] ];
print "(2) over Q(t): Da = phi1^*(P') and Db = phi2^*(Q1') have order 5 and Db is not in <Da>,";
print "    so G = <Da, Db> = [5,5] is marked over Q(t)";

// ---------------------------------------------------------------- (3) exactness on the residue class t = 1 mod 3
p := 3; r0 := 1;
members := [Q| 1, 4, -2, 7, 1/4];
sanity := [Q| 1, 4, -2];
fbar := ReductionClass(f, p, r0);
for m in members do
    assert Valuation(m, p) ge 0 and GF(p)!m eq GF(p)!r0;
    assert PolynomialRing(GF(p))!SpecPoly(f, m) eq fbar;
end for;
n, chi := FrobData(fbar);
inv := Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
assert inv eq [5, 5] and n eq 25;
printf "(3) p = %o, t = %o mod %o: reduction %o (squarefree, degree %o); J(F_%o) = %o, chi = %o\n", p, r0, p, fbar, Degree(fbar), p, inv, chi;
print "    #J_t(Q)_tors divides #J(F_3) = 25 = #G, so every 3-integral t = 1 mod 3 has J_t(Q)_tors exactly [5,5]";

// ---------------------------------------------------------------- (4) nonconstant moduli map, (5) sanity
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } eq #members;
printf "(4) the G2-invariants of the members %o are pairwise distinct\n", members;
C0 := HyperellipticCurve(RX![-2,4,2,5,2,1], RX![1,1,1,1]);
Cm := HyperellipticCurve(SpecPoly(f, 1));
assert IsIsomorphic(Cm, SimplifiedModel(C0));
print "    the member t = 1 is the census curve 5.5.a (isomorphic over Q)";
for m in sanity do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(5) TorsionSubgroup at t = %o: %o\n", m, invm;
    assert invm eq [5, 5];
end for;
printf "time %o s\n", Cputime(T0);
print "FAMILY_5_5_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E x E1) and J(Q)_tors exactly [5,5]";
quit;
