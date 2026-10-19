/* family_2_4_4_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E_s x E_s) with J(Q)_tors exactly [2,4,4].
   Found 2026-10-19 (torsion_jac notes/isopair_2026_10_19.md).

   Source.  A Howe-Leprevost-Poonen (2,2)-self-gluing (HLP 2000, Prop. 4; the construction of their Corollary 6) of a
   curve with Z/2 x Z/4 with itself, along an automorphism of its 2-torsion, on the locus where a Kummer coincidence
   makes the gain of the gluing contain a point of order 4.  The curve is
       C_s : y^2 = f_s(x) = -(s^4 + 1)(x^2 + 1)(x^2 - s^4)(s^4 x^2 - 1),
   f_s = F(x^2) with F(X) = -(s^4 + 1)(X + 1)(X - s^4)(s^4 X - 1) = X^3 F(1/X), so phi1 : (x,y) -> (x^2, y) and
   phi2 : (x,y) -> (1/x^2, y/x^3) both map C_s onto E_s : Y^2 = F(X), and J is isogenous over Q to E_s x E_s.
   E_s is Q(s)-isomorphic to the Kubert curve E_{2,4}(t) : y^2 + xy - by = x^3 - bx^2, b = t^2 - 1/16, at
   t = m/(2(m^2 + 1)), m = (s + 1)/(s - 1), so E_s(Q(s)) contains Z/2 x Z/4; P4 below has 2 P4 = (-1, 0).  The gluing
   identifies (R, 0) on the first copy with (1/R, 0) on the second: it fixes 2 P4 and swaps (s^4, 0) and (s^-4, 0).
   On E_{2,4}(t) the Kummer class of P4 is invariant under this swap iff 1 - 16t^2 is a square; at t = m/(2(m^2 + 1))
   it is ((m^2 - 1)/(m^2 + 1))^2.  Then the halves of P4 are defined over a quadratic extension and differ there by
   2 P4, and the class W = phi1^*(Q') + phi2^*(Q'), 2Q = P4, is rational of order 4.  The census curve 2.4.4.b is the
   member s = 2 (and s = 1/2, -2, -1/2).
   Argument.
   (0) E_s, its isomorphism with E_{2,4}(t(s)), and P4 of order 4 with 2 P4 = (-1, 0).
   (1) Splitting and gluing.  F(X) = X^3 F(1/X), so Ea = Eb = E_s; phi1^* + phi2^* : E_s x E_s -> J is an isogeny
       over Q with kernel of order 4.  The gluing map (R, 0) -> (1/R, 0) on E_s[2] is not the identity, so it is not
       induced by an automorphism of E_s (+-1 act trivially on E_s[2]).
   (2) The marked group over Q(s).  D1 = phi1^*(P4') has order 4 with 2 D1 = [x^2 + 1, 0].  The halving quartic of P4
       splits into two quadratics; for Q a half of P4 over the quadratic field, W = phi1^*(Q') + phi2^*(Q') has
       Mumford coefficients in Q(s) and order 4.  With c = [x^2 - s^4, 0], the 32 classes i D1 + j W + k c are
       distinct, so G = <D1, W, c> = Z/4 x Z/4 x Z/2 = [2,4,4] is marked over Q(s) and embeds in J_s(Q)_tors.
   (3) Exactness on a residue class: every s = 86 (mod 287 = 7 * 41), integral at 7 and 41, has the reductions of C_86
       at 7 and 41, where J(F_7) = [2,4,8] (a 2-group) and J(F_41) = [4,4,12,12].  J_s(Q)_tors embeds in both, and the
       largest group that embeds in both is [2,4,4], so J_s(Q)_tors is exactly [2,4,4].  (#J(F_p) = #E_s(F_p)^2 is
       always divisible by 64, so a gcd of point counts cannot certify this; the group structures are used instead.)
   (4) Infinitely many.  The G2-invariants of the members differ, so s -> M_2 is nonconstant with finite fibres.
   (5) Sanity: TorsionSubgroup at s = 86, -201, 373 is [2,4,4]; the member s = 2 (outside the class) is the census
       curve 2.4.4.b, and s = 4/3 (outside the class) has J(Q)_tors = [2,4,8] (E_s(Q) = Z/2 x Z/8 there).
   Run (from the repository root):
       magma -b pipeline/magma/family_2_4_4_split_certificate.m > data/knowledge/sources/family_2_4_4_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

Ft<s> := FunctionField(Q); P<x> := PolynomialRing(Ft);
f := -(s^4 + 1)*(x^2 + 1)*(x^2 - s^4)*(s^4*x^2 - 1);
assert &and[ Coefficient(f, i) eq 0 : i in [1, 3, 5] ] and Degree(f) eq 6 and IsSquarefree(f);
F := P![ Coefficient(f, 2*i) : i in [0..3] ];            // f = F(x^2)
printf "C_s : y^2 = f_s(x) = %o\n", f;
printf "    the irreducible factors of f_s over Q(s) have degrees %o\n", [ Degree(h[1]) : h in Factorization(f) ];

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

// ---------------------------------------------------------------- (0) source: E_s = E_{2,4}(t(s)) and P4
F := P![ Coefficient(f, 2*i) : i in [0..3] ];
Ea, ca := MonEC(F);                                   // E_s : Y^2 = F(X), as a Weierstrass model (X' = ca X)
mm := (s + 1)/(s - 1); tt := mm/(2*(mm^2 + 1)); bK := tt^2 - 1/16;
E24 := EllipticCurve([Ft| 1, -bK, -bK, 0, 0]);
assert 4*E24![0,0] eq E24!0 and 2*E24![0,0] ne E24!0 and #Roots(MonicCubic(E24)) eq 3;
ok24, m24 := IsIsomorphic(E24, Ea);
assert ok24 and IsSquare(1 - 16*tt^2);
P4 := Ea![ ca*(-s^4 - s^2 - 1)/s^2, ca*(-s^10 - s^8 - 2*s^6 - 2*s^4 - s^2 - 1)/s^2 ];
assert 4*P4 eq Ea!0 and 2*P4 eq Ea![-ca, 0];
print "(0) E_s : Y^2 = F(X) = -(s^4 + 1)(X + 1)(X - s^4)(s^4 X - 1) is Q(s)-isomorphic to E_{2,4}(t) : y^2 + xy - by = x^3 - bx^2,";
print "    b = t^2 - 1/16, at t = m/(2(m^2 + 1)), m = (s + 1)/(s - 1) (where 1 - 16t^2 is a square): E_s(Q(s)) contains Z/2 x Z/4;";
print "    P4 = (X, Y) = (-(s^4 + s^2 + 1)/s^2, -(s^10 + s^8 + 2s^6 + 2s^4 + s^2 + 1)/s^2) has order 4 and 2 P4 = (-1, 0)";

// ---------------------------------------------------------------- (1) the elliptic quotients and the gluing
Fr := P![ Coefficient(F, 3 - i) : i in [0..3] ];
assert Fr eq F;
print "(1) F(X) = X^3 F(1/X): phi1 = (x^2, y) and phi2 = (1/x^2, y/x^3) both map C_s onto E_s : Y^2 = F(X), so J is isogenous";
print "    over Q to E_s x E_s (phi1^* + phi2^* has kernel of order 4); the gluing identifies (R, 0) with (1/R, 0) on E_s[2]:";
print "    it fixes 2 P4 = (-1, 0) and swaps (s^4, 0) and (s^-4, 0), so it is not induced by an automorphism of E_s";

// ---------------------------------------------------------------- (2) the marked group over Q(s)
J := Jacobian(HyperellipticCurve(f));
D1 := Pull1(J, P4, ca);
assert 4*D1 eq J!0 and 2*D1 eq elt< J | x^2 + 1, 0 > and 2*D1 ne J!0;
a1, a2, a3, a4, a6 := Explode(aInvariants(Ea));
b2, b4, b6, b8 := Explode(bInvariants(Ea));
quart := x^4 - b4*x^2 - 2*b6*x - b8 - P4[1]*(4*x^3 + b2*x^2 + 2*b4*x + b6);      // x(2Q) = x(P4)
fq := Factorization(quart);
assert #fq eq 2 and &and[ Degree(h[1]) eq 2 : h in fq ];
K1<w> := ext< Ft | fq[1][1] >;
xQ := w; x4 := K1!P4[1]; y4m := K1!(-P4[2] - a1*P4[1] - a3);
cub := xQ^3 + a2*xQ^2 + a4*xQ + a6; L := a1*xQ + a3; S := 3*xQ^2 + 2*a2*xQ + a4;
yQ := (2*cub - y4m*L + S*(x4 - xQ))/(L + 2*y4m + a1*(x4 - xQ));     // the tangent at Q passes through -P4
EK := BaseChange(Ea, K1); Qpt := EK![xQ, yQ];
assert 2*Qpt eq EK!Eltseq(P4);
JK := BaseChange(J, K1);
Gk := Pull1(JK, Qpt, K1!ca) + Pull2(JK, Qpt, K1!ca);
assert &and[ IsCoercible(Ft, c) : c in Coefficients(Gk[1]) cat Coefficients(Gk[2]) ];
W := elt< J | P![ Ft!c : c in Coefficients(Gk[1]) ], P![ Ft!c : c in Coefficients(Gk[2]) ] >;
assert 4*W eq J!0 and 2*W ne J!0;
cc := elt< J | x^2 - s^4, 0 >;
assert cc ne J!0 and 2*cc eq J!0;
els := { i*D1 + j*W + k*cc : i in [0..3], j in [0..3], k in [0..1] };
assert #els eq 32;
print "(2) over Q(s): D1 = phi1^*(P4') has order 4, 2 D1 = [x^2 + 1, 0]; the halving quartic of P4 splits into two quadratics,";
print "    and for a half Q of P4 over the quadratic field the class W = phi1^*(Q') + phi2^*(Q') lies in J(Q(s)):";
printf "    W = (%o, %o)\n", W[1], W[2];
print "    W has order 4, and with c = [x^2 - s^4, 0] the 32 classes i D1 + j W + k c are distinct:";
print "    G = <D1, W, c> = Z/4 x Z/4 x Z/2 = [2,4,4] is marked over Q(s)";

// ---------------------------------------------------------------- (3) exactness on the residue class s = 86 mod 287
function StructureBound(invlist)        // the largest finite abelian group embedding in every group of the list
    ls := &meet[ Seqset(&cat[ PrimeDivisors(n) : n in inv | n gt 1 ]) : inv in invlist ];
    out := [Integers()| ];
    for l in ls do
        es := [ Reverse(Sort([ Valuation(n, l) : n in inv | n mod l eq 0 ])) : inv in invlist ];
        r := Min([ #e : e in es ]);
        out cat:= [ l^Min([ e[i] : e in es ]) : i in [1..r] ];
    end for;
    return Invariants(AbelianGroup(out));
end function;
primes := [7, 41]; residues := [2, 4];
assert CRT(residues, primes) eq 86;
members := [Q| 86, -201, 373, 86 + 287/2, 86 + 287/3];
sanity := [Q| 86, -201, 373];
invs := [];
for i in [1..2] do
    p := primes[i]; r0 := residues[i];
    fbar := ReductionClass(f, p, r0);
    for m in members do
        assert Valuation(m, p) ge 0 and GF(p)!m eq GF(p)!r0;
        assert PolynomialRing(GF(p))!SpecPoly(f, m) eq fbar;
    end for;
    n, chi := FrobData(fbar);
    inv := Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
    Append(~invs, inv);
    printf "(3) p = %o, s = %o mod %o: reduction %o (squarefree, degree %o); J(F_%o) = %o, chi = %o\n", p, r0, p, fbar, Degree(fbar), p, inv, chi;
end for;
assert invs eq [[2, 4, 8], [4, 4, 12, 12]];
bnd := StructureBound(invs);
assert bnd eq [2, 4, 4];
printf "    the largest group embedding in all J(F_p) has invariants %o = G\n", bnd;
print "    so every s = 86 mod 287, integral at 7 and 41, has J_s(Q)_tors exactly [2,4,4]";

// ---------------------------------------------------------------- (4) nonconstant moduli map, (5) sanity
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } eq #members;
printf "(4) the G2-invariants of the members %o are pairwise distinct\n", members;
C0 := HyperellipticCurve(RX![0,15,40,0,161,-240], RX![0,1]);
Cm := HyperellipticCurve(SpecPoly(f, 2));
assert IsIsomorphic(Cm, SimplifiedModel(C0));
print "    the member s = 2 is the census curve 2.4.4.b (isomorphic over Q)";
inv43 := TorsionInvariants(SpecPoly(f, 4/3));
assert inv43 eq [2, 4, 8];
print "    the member s = 4/3 (outside the class) has J(Q)_tors = [2,4,8]";
for m in sanity do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(5) TorsionSubgroup at s = %o: %o\n", m, invm;
    assert invm eq [2, 4, 4];
end for;
printf "time %o s\n", Cputime(T0);
print "FAMILY_2_4_4_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E_s x E_s) and J(Q)_tors exactly [2,4,4]";
quit;
