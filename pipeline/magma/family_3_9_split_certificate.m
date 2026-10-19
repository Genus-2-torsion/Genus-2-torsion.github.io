/* family_3_9_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E x E1) with J(Q)_tors exactly [3,9].
   Found 2026-10-19 (torsion_jac notes/isopair_2026_10_19.md).

   Source.  The Howe-Leprevost-Poonen (2,2)-gluing (HLP 2000, Prop. 4) of the Kubert X_1(9) curve
       E = E_9(t) : y^2 + (1 - c)xy - by = x^3 - bx^2,  c = t^2(t - 1),  b = c(t^2 - t + 1)     (P = (0,0) of order 9)
   with its 3-isogenous curve E1 = E/<3P>, along the isomorphism E[2] -> E1[2] induced by the 3-isogeny phi.
   phi(P) has order 3, and the kernel of the dual isogeny is mu_3, so E1[3] = <phi(P)> + mu_3.  The gluing has
   coefficients in Q(t); after removing square factors and scaling x it is the curve f_t below, f_t = F(x^2) with F a
   cubic over Q(t), and phi1 : (x,y) -> (x^2, y) maps C_t onto Ea : Y^2 = F(X), phi2 : (x,y) -> (1/x^2, y/x^3) onto
   Eb : Y^2 = X^3 F(1/X); below, Ea is Q(t)-isomorphic to E and Eb to E1.
   Argument.
   (0) E, P, E1 = E/<3P>, phi; mu_3 inside E1[3] (a second linear factor of the 3-division polynomial of E1 whose points
       are defined over Q(t)(sqrt -3) and swapped by its conjugation).
   (1) Splitting and gluing.  Ea = E and Eb = E1 over Q(t); phi1^* + phi2^* : Ea x Eb -> J is an isogeny over Q with
       kernel of order 4, so J[3] = Ea[3] + Eb[3] as Galois modules.  The gluing is along phi on E[2].
   (2) The marked group over Q(t).  Da = phi1^*(P') has order 9, Db = phi2^*(phi(P)') has order 3, and Db is not in
       <3 Da>, so G = <Da, Db> = [3,9] is marked over Q(t) and embeds in J_t(Q)_tors.
   (3) Exactness on a residue class.  For every 7-integral t = 2 (mod 7) the curve C_t has the reduction of C_2 at 7,
       with J(F_7) = [3,3,9] of order 81.  So J_t(Q)_tors is a 3-group embedding in [3,3,9].  By (1),
       J_t(Q)[3] = Ea(Q)[3] + Eb(Q)[3], and each summand is cyclic (Weil pairing; Q has no primitive cube root of unity),
       so J_t(Q)_tors has 3-rank <= 2.  A group of 3-rank <= 2 containing [3,9] and different from it contains [9,9] or
       [3,27], and neither embeds in [3,3,9].  Hence J_t(Q)_tors is exactly [3,9].  (#J(F_p) = #E(F_p)^2 is always
       divisible by 81, so a gcd of point counts cannot certify this.)
   (4) Infinitely many.  The G2-invariants of the members differ, so t -> M_2 is nonconstant with finite fibres.
   (5) Sanity: TorsionSubgroup at t = 2, 9, -5 is [3,9].
   Run (from the repository root):
       magma -b pipeline/magma/family_3_9_split_certificate.m > data/knowledge/sources/family_3_9_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
f := (-t^9 + 6*t^8 - 12*t^7 + 9*t^6 - 3*t^4 + t^3)*x^6 + (t^12 - 12*t^11 + 54*t^10 - 113*t^9 + 99*t^8 - 21*t^6 - 36*t^5 + 27*t^4 + 13*t^3 - 12*t^2 + 1)*x^4 + (-2*t^12 + 24*t^11 - 108*t^10 + 181*t^9 + 72*t^8 - 540*t^7 + 447*t^6 + 72*t^5 - 189*t^4 + 19*t^3 + 24*t^2 - 2)*x^2 + t^12 - 12*t^11 + 54*t^10 - 131*t^9 + 207*t^8 - 216*t^7 + 141*t^6 - 36*t^5 - 27*t^4 + 31*t^3 - 12*t^2 + 1;
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

// ---------------------------------------------------------------- (0) source: E = E_9(t), E1 = E/<3P>
cK := t^2*(t - 1); dK := t^2 - t + 1; bK := cK*dK;
E := EllipticCurve([Ft| 1 - cK, -bK, -bK, 0, 0]);
P9 := E![0, 0];
assert 9*P9 eq E!0 and 3*P9 ne E!0;
E1, phi := IsogenyFromKernel(E, x - (3*P9)[1]);
assert Degree(phi) eq 3;
dp := DefiningPolynomials(phi);
P9b := E1![ Evaluate(d, Eltseq(P9)) : d in dp ];
assert P9b ne E1!0 and 3*P9b eq E1!0;
// mu_3 = kernel of the dual isogeny, inside E1[3]
b1, b2, b3, b4, b6 := Explode(aInvariants(E1));
lins := [ h[1] : h in Factorization(DivisionPolynomial(E1, 3)) | Degree(h[1]) eq 1 ];
x0s := [ -Coefficient(h, 0)/Coefficient(h, 1) : h in lins ];
assert #x0s eq 2 and P9b[1] in x0s;
x0 := [ z : z in x0s | z ne P9b[1] ][1];
bb := b1*x0 + b3; dsc := bb^2 + 4*(x0^3 + b2*x0^2 + b4*x0 + b6);
okmu, w0 := IsSquare(-3*dsc);
assert okmu and dsc ne 0;
print "(0) E = E_9(t) : y^2 + (1 - c)xy - by = x^3 - bx^2, c = t^2(t - 1), b = c(t^2 - t + 1), P = (0,0) of order 9;";
print "    E1 = E/<3P>, phi : E -> E1 of degree 3, phi(P) has order 3; the 3-division polynomial of E1 has a second";
printf "    linear factor x - x0 with -3 (disc of the y-quadratic at x0) = (%o)^2, so E1[3] = <phi(P)> + mu_3\n", w0;

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
Da := Pull1(J, ma(P9), ca); Db := Pull2(J, mb(P9b), cb);
assert 9*Da eq J!0 and 3*Da ne J!0 and Db ne J!0 and 3*Db eq J!0;
assert Db notin [J!0, 3*Da, 6*Da];
print "(2) over Q(t): Da = phi1^*(P') has order 9, Db = phi2^*(phi(P)') has order 3 and is not in <3 Da>,";
print "    so G = <Da, Db> = [3,9] is marked over Q(t)";

// ---------------------------------------------------------------- (3) exactness on the residue class t = 2 mod 7
p := 7; r0 := 2;
members := [Q| 2, 9, -5, 16, 23, 9/8];
sanity := [Q| 2, 9, -5];
fbar := ReductionClass(f, p, r0);
for m in members do
    assert Valuation(m, p) ge 0 and GF(p)!m eq GF(p)!r0;
    assert PolynomialRing(GF(p))!SpecPoly(f, m) eq fbar;
end for;
n, chi := FrobData(fbar);
inv := Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
assert inv eq [3, 3, 9] and n eq 81;
assert not EmbedsIn([9, 9], inv) and not EmbedsIn([3, 27], inv) and EmbedsIn([3, 9], inv);
printf "(3) p = %o, t = %o mod %o: reduction %o (squarefree, degree %o); J(F_%o) = %o, chi = %o\n", p, r0, p, fbar, Degree(fbar), p, inv, chi;
print "    J_t(Q)_tors is a 3-group of 3-rank <= 2 (J(Q)[3] = Ea(Q)[3] + Eb(Q)[3], each cyclic) containing [3,9],";
print "    and neither [9,9] nor [3,27] embeds in [3,3,9]: every 7-integral t = 2 mod 7 has J_t(Q)_tors exactly [3,9]";

// ---------------------------------------------------------------- (4) nonconstant moduli map, (5) sanity
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } eq #members;
printf "(4) the G2-invariants of the members %o are pairwise distinct\n", members;

for m in sanity do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(5) TorsionSubgroup at t = %o: %o\n", m, invm;
    assert invm eq [3, 9];
end for;
printf "time %o s\n", Cputime(T0);
print "FAMILY_3_9_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E x E1) and J(Q)_tors exactly [3,9]";
quit;
