/* family_15_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E5 x E3) with J(Q)_tors exactly [15].
   Found 2026-10-19 (torsion_jac notes/hlpP2_2026_10_19.md).

   Source.  Howe-Leprevost-Poonen 2000, Lemma 11 and Prop. 4: the construction of their Z/21 entry (Table 1,
   section 3.6) with the Z/7 curve replaced by a Z/5 curve.  E5 = E_5(t) is the Kubert (Tate normal form) curve
       E_5(t) : y^2 + (1 - c)xy - by = x^3 - bx^2,   b = c = t,   P5 = (0,0) of order 5,
   with short model y^2 = x^3 + Ax + B, A = -27 c4, B = -54 c6.  Lemma 11: with q = B^2/A^3 the curve
       E3 : y^2 = x^3 + (x - q)^2,   P3 = (0, q) of order 3,
   has 2-torsion Galois-isomorphic to that of E5 (beta = -alpha^2/A + (B/A^2) alpha - 1 for alpha^3 + A alpha + B = 0).
   Prop. 4 for this isomorphism gives, after removing a square,
       C_t : y^2 = f_t(x) = A^7 - 3A^5 B x^2 + (A^6 + 3A^3 B^2) x^4 - A B^3 x^6.
   f_t = F(x^2), so phi1 : (x,y) -> (x^2, y) maps C_t onto Ea : Y^2 = F(X) and phi2 : (x,y) -> (1/x^2, y/x^3)
   onto Eb : Y^2 = X^3 F(1/X); below, Ea is Q(t)-isomorphic to E3 and Eb to E5.

   Argument.
   (1) Splitting.  Wherever f_t is squarefree of degree 6, phi1^* + phi2^* : Ea x Eb -> J is an isogeny over Q
       (the involutions (x,y) -> (-x,y) and (x,y) -> (-x,-y) compose to the hyperelliptic involution), so J is
       isogenous over Q to E3 x E5.
   (2) The marked group over Q(t).  D5 = phi2^*(image of P5) has order 5 and D3 = phi1^*(image of P3) has order 3
       (Mumford arithmetic over Q(t)), so D15 = D5 + D3 has order 15 and Z/15 embeds in J_t(Q)_tors at every
       smooth member.
   (3) Exactness on a residue class.  For every 7- and 11-integral t = 4 (mod 77) the curve C_t has the
       reductions of C_4 at 7 and 11, with #J(F_7) = 60 (J(F_7) = [2,30]) and #J(F_11) = 225 (J(F_11) = [15,15]).
       Reduction is injective on J_t(Q)_tors (p odd, good reduction), so #J_t(Q)_tors divides gcd(60, 225) = 15,
       and J_t(Q)_tors is exactly [15].
   (4) Infinitely many.  The G2-invariants of the members differ, so t -> M_2 is nonconstant with finite fibres,
       and the residue class gives infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup at t = 4, 81, -73 is [15].
   Run (from the repository root):
       magma -b pipeline/magma/family_15_split_certificate.m > data/knowledge/sources/family_15_split_certificate.log */
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

// ---------------------------------------------------------------- (0) source: E_5(t), the Lemma 11 partner E3, the HLP sextic
bK := t; cK := t;
E5 := EllipticCurve([Ft| 1 - cK, -bK, -bK, 0, 0]);
P5 := E5![0, 0];
assert 5*P5 eq E5!0 and P5 ne E5!0;
c4, c6 := Explode(cInvariants(E5));
A := -27*c4; B := -54*c6;
E5s := EllipticCurve([Ft| 0, 0, 0, A, B]);
assert IsIsomorphic(E5, E5s);
q := B^2/A^3;
E3 := EllipticCurve([Ft| 0, 1, 0, -2*q, q^2]);             // y^2 = x^3 + (x - q)^2
P3 := E3![0, q];
assert 3*P3 eq E3!0 and P3 ne E3!0;
// Lemma 11: the root map alpha -> beta is a Galois-equivariant bijection E5[2] -> E3[2]
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
okb, mb := IsIsomorphic(E5, Eb);
assert oka and okb;
print "(1) phi1 = (x^2, y) : C_t -> Ea : Y^2 = F(X), phi2 = (1/x^2, y/x^3) : C_t -> Eb : Y^2 = X^3 F(1/X);";
print "    Ea is Q(t)-isomorphic to E3 and Eb to E5 = E_5(t): J is isogenous over Q to E3 x E5";

// ---------------------------------------------------------------- (2) the marked group over Q(t)
J := Jacobian(HyperellipticCurve(f));
D5 := Pull2(J, mb(P5), cb);
D3 := Pull1(J, ma(P3), ca);
assert 5*D5 eq J!0 and D5 ne J!0 and 3*D3 eq J!0 and D3 ne J!0;
D15 := D5 + D3;
assert 15*D15 eq J!0 and 5*D15 ne J!0 and 3*D15 ne J!0;
print "(2) over Q(t): D5 = phi2^*(P5) has order 5, D3 = phi1^*(P3) has order 3, D15 = D5 + D3 has order 15";

// ---------------------------------------------------------------- (3) exactness on the residue class t = 4 mod 77
primes := [7, 11]; t0 := 4; modulus := 77;
members := [Q| 4, 81, -73, 158, 235];
nJ := [];
for p in primes do
    r0 := Integers()!(GF(p)!t0);
    fbar := ReductionClass(f, p, r0);
    for m in members do
        assert GF(p)!m eq GF(p)!r0;
        assert PolynomialRing(GF(p))!SpecPoly(f, m) eq fbar;
    end for;
    n, chi := FrobData(fbar);
    inv := Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
    assert &*inv eq n;
    Append(~nJ, n);
    printf "(3) p = %o, t = %o mod %o: reduction %o (squarefree, degree %o); #J(F_%o) = %o, J(F_%o) = %o, chi = %o\n", p, r0, p, fbar, Degree(fbar), p, n, p, inv, chi;
end for;
assert GCD(nJ) eq 15 and nJ eq [60, 225];
printf "    gcd of #J(F_p) = %o: every 7- and 11-integral t = %o mod %o has J_t(Q)_tors exactly [15]\n", GCD(nJ), t0, modulus;

// ---------------------------------------------------------------- (4) nonconstant moduli map, (5) sanity
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } eq #members;
printf "(4) the G2-invariants of the members %o are pairwise distinct\n", members;
for m in [Q| 4, 81, -73] do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(5) TorsionSubgroup at t = %o: %o\n", m, invm;
    assert invm eq [15];
end for;
printf "time %o s\n", Cputime(T0);
print "FAMILY_15_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E1 x E2) and J(Q)_tors exactly [15]";
quit;
