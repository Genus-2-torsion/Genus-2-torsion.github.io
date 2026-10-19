/* family_21_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E7 x E3) with J(Q)_tors exactly [21].
   Found 2026-10-19 (torsion_jac notes/hlpP2_2026_10_19.md; the family is the one of
   notes/z21_leprevost_and_hlp_2026_07_18.md, "Compact HLP gluing family").

   Source.  Howe-Leprevost-Poonen 2000, Lemma 11 and Prop. 4 (the Z/21 entry of their Table 1, restricted to
   the line on which E3 is the partner given by Lemma 11).  E7 = E_7(t) is the Kubert (Tate normal form) curve
       E_7(t) : y^2 + (1 - c)xy - by = x^3 - bx^2,   b = t^3 - t^2,  c = t^2 - t,   P7 = (0,0) of order 7,
   with short model y^2 = x^3 + Ax + B, A = -27 c4, B = -54 c6.  Lemma 11: with q = B^2/A^3 the curve
       E3 : y^2 = x^3 + (x - q)^2,   P3 = (0, q) of order 3,
   has 2-torsion Galois-isomorphic to that of E7 (beta = -alpha^2/A + (B/A^2) alpha - 1 for alpha^3 + A alpha + B = 0).
   Prop. 4 for this isomorphism gives, after removing a square,
       C_t : y^2 = f_t(x) = A^7 - 3A^5 B x^2 + (A^6 + 3A^3 B^2) x^4 - A B^3 x^6.
   f_t = F(x^2), so phi1 : (x,y) -> (x^2, y) maps C_t onto Ea : Y^2 = F(X) and phi2 : (x,y) -> (1/x^2, y/x^3)
   onto Eb : Y^2 = X^3 F(1/X); below, Ea is Q(t)-isomorphic to E3 and Eb to E7.

   Argument.
   (1) Splitting.  Wherever f_t is squarefree of degree 6, phi1^* + phi2^* : Ea x Eb -> J is an isogeny over Q
       (the involutions (x,y) -> (-x,y) and (x,y) -> (-x,-y) compose to the hyperelliptic involution), so J is
       isogenous over Q to E3 x E7.
   (2) The marked group over Q(t).  D7 = phi2^*(image of P7) has order 7 and D3 = phi1^*(image of P3) has order 3
       (Mumford arithmetic over Q(t)), so D21 = D7 + D3 has order 21 and Z/21 embeds in J_t(Q)_tors at every
       smooth member.
   (3) Exactness on a residue class.  For every 5-integral t = 2 (mod 5) the curve C_t has the reduction of C_2
       at 5, with #J(F_5) = 21.  Reduction is injective on J_t(Q)_tors (p odd, good reduction), so J_t(Q)_tors
       is exactly [21].
   (4) Infinitely many.  The G2-invariants of the members differ, so t -> M_2 is nonconstant with finite fibres,
       and the residue class gives infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup at t = 2, 7, -3 is [21].
   Run (from the repository root):
       magma -b pipeline/magma/family_21_split_certificate.m > data/knowledge/sources/family_21_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
SetQuitOnError(true);                   // fail closed: a failed assertion stops the run before the marker
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

function MonEC(G)                       // Y^2 = G(X), G a cubic, as a Weierstrass curve (X' = lc X, Y' = lc Y)
    cs := [ Coefficient(G, i) : i in [0..3] ];
    return EllipticCurve([0, cs[3], 0, cs[2]*cs[4], cs[1]*cs[4]^2]), cs[4];
end function;
// pullbacks of points (X, Y) of the Weierstrass models (scaled by c) to the Jacobian of y^2 = F(x^2)
function Pull1(JJ, Pt, c) PP := Parent((JJ!0)[1]); xx := PP.1; XX := Pt[1]/c; YY := Pt[2]/c; return elt< JJ | xx^2 - XX, PP!YY >; end function;
function Pull2(JJ, Pt, c) PP := Parent((JJ!0)[1]); xx := PP.1; XX := Pt[1]/c; YY := Pt[2]/c; return elt< JJ | xx^2 - 1/XX, (YY/XX)*xx >; end function;

Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);

// ---------------------------------------------------------------- (0) source: E_7(t), the Lemma 11 partner E3, the HLP sextic
bK := t^3 - t^2; cK := t^2 - t;
E7 := EllipticCurve([Ft| 1 - cK, -bK, -bK, 0, 0]);
P7 := E7![0, 0];
assert 7*P7 eq E7!0 and P7 ne E7!0;
c4, c6 := Explode(cInvariants(E7));
A := -27*c4; B := -54*c6;
E7s := EllipticCurve([Ft| 0, 0, 0, A, B]);
assert IsIsomorphic(E7, E7s);
q := B^2/A^3;
E3 := EllipticCurve([Ft| 0, 1, 0, -2*q, q^2]);             // y^2 = x^3 + (x - q)^2
P3 := E3![0, q];
assert 3*P3 eq E3!0 and P3 ne E3!0;
// Lemma 11: the root map alpha -> beta is a Galois-equivariant bijection E7[2] -> E3[2]
Kc<al> := ext< Ft | x^3 + A*x + B >;
be := -al^2/A + (B/A^2)*al - 1;
assert be^3 + (be - q)^2 eq 0;
f := A^7 - 3*A^5*B*x^2 + (A^6 + 3*A^3*B^2)*x^4 - A*B^3*x^6;
assert Degree(f) eq 6 and IsSquarefree(f);
assert Discriminant(f) eq 64*A^52*B^3*(4*A^3 + 27*B^2)^2;
printf "A = %o\nB = %o\n", A, B;
print "C_t : y^2 = f_t(x) = A^7 - 3A^5 B x^2 + (A^6 + 3A^3 B^2) x^4 - A B^3 x^6";

// ---------------------------------------------------------------- (1) the two elliptic quotients
F := P![ Coefficient(f, 2*i) : i in [0..3] ];
Fr := P![ Coefficient(F, 3 - i) : i in [0..3] ];
Ea, ca := MonEC(F); Eb, cb := MonEC(Fr);
oka, ma := IsIsomorphic(E3, Ea);
okb, mb := IsIsomorphic(E7, Eb);
assert oka and okb;
print "(1) phi1 = (x^2, y) : C_t -> Ea : Y^2 = F(X), phi2 = (1/x^2, y/x^3) : C_t -> Eb : Y^2 = X^3 F(1/X);";
print "    Ea is Q(t)-isomorphic to E3 and Eb to E7 = E_7(t): J is isogenous over Q to E3 x E7";

// ---------------------------------------------------------------- (2) the marked group over Q(t)
J := Jacobian(HyperellipticCurve(f));
D7 := Pull2(J, mb(P7), cb);
D3 := Pull1(J, ma(P3), ca);
assert 7*D7 eq J!0 and D7 ne J!0 and 3*D3 eq J!0 and D3 ne J!0;
D21 := D7 + D3;
assert 21*D21 eq J!0 and 7*D21 ne J!0 and 3*D21 ne J!0;
print "(2) over Q(t): D7 = phi2^*(P7) has order 7, D3 = phi1^*(P3) has order 3, D21 = D7 + D3 has order 21";

// ---------------------------------------------------------------- (3) exactness on the residue class t = 2 mod 5
p := 5; r0 := 2;
members := [Q| 2, 7, -3, 12, 17];
fbar := ReductionClass(f, p, r0);
for m in members do
    assert GF(p)!m eq GF(p)!r0;
    assert PolynomialRing(GF(p))!SpecPoly(f, m) eq fbar;
end for;
n, chi := FrobData(fbar);
inv := Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
assert n eq 21 and inv eq [21];
printf "(3) p = %o, t = %o mod %o: reduction %o (squarefree, degree %o); #J(F_%o) = %o, J(F_%o) = %o, chi = %o\n", p, r0, p, fbar, Degree(fbar), p, n, p, inv, chi;
print "    hence every 5-integral t = 2 mod 5 has J_t(Q)_tors exactly [21]";

// ---------------------------------------------------------------- (4) nonconstant moduli map, (5) sanity
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } eq #members;
printf "(4) the G2-invariants of the members %o are pairwise distinct\n", members;
for m in [Q| 2, 7, -3] do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(5) TorsionSubgroup at t = %o: %o\n", m, invm;
    assert invm eq [21];
end for;
printf "time %o s\n", Cputime(T0);
print "FAMILY_21_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E1 x E2) and J(Q)_tors exactly [21]";
quit;
