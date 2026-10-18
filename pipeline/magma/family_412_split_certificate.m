/* family_412_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E1 x E2) with J(Q)_tors exactly [4,12].
   Found 2026-10-18 (torsion_jac notes/f412_2026_10_18.md).

   Source.  The Howe-Leprevost-Poonen (2,2)-gluing (HLP 2000, Prop. 4) of the Kubert X_1(12) curve
       E2 = E_12(t) : y^2 + (1 - c)xy - by = x^3 - bx^2,  m = (3t - 3t^2 - 1)/(t - 1), f = m/(1 - t),
                      d = m + t, c = f(d - 1), b = cd,
   (torsion Z/12, P = (0,0) of order 12) with its 3-isogenous curve E1 = E2/<4P>, along the isomorphism
   E2[2] -> E1[2] induced by the 3-isogeny phi (the route of the [2,10] and [2,2,2,6] certificates).
   The HLP formula is evaluated over K = Q(t)(sqrt delta), delta the discriminant of the quadratic factor of
   the 2-division polynomial of E2; the result has coefficients in Q(t) and is, up to a square, with tau = t^2 - t,
       C_t : y^2 = f_t(x) = -(6tau^2 + 6tau + 1) ((2tau + 1)(4tau + 1)x^2 - (6tau + 1)) q(x) q(-x),
             q(x) = tau^2 x^2 - (6tau^2 + 6tau + 1)x + (3tau + 1)^2.
   f_t = F(x^2) with F a cubic over Q(t).  phi1 : (x,y) -> (x^2, y) maps C_t onto Ea : Y^2 = F(X),
   phi2 : (x,y) -> (1/x^2, y/x^3) onto Eb : Y^2 = X^3 F(1/X); below, Ea is Q(t)-isomorphic to E2 and Eb to E1.
   The census curve 4.12.a (Balakrishnan-Najman-Shnidman-Sutherland, table 2) is the member t = 1/3 (and 1 - t = 2/3).

   Argument.
   (1) Splitting.  Wherever f_t is squarefree of degree 6, phi1^* + phi2^* : Ea x Eb -> J is an isogeny over Q
       with kernel of order 4, so J is split over Q and J(Q)[l^inf] = Ea(Q)[l^inf] x Eb(Q)[l^inf] for odd l.
   (2) The marked group over Q(t).  D12 = phi1^*(Pa), Pa the image of P on Ea, has order 12.  The gain of the
       gluing: for Q in E2(Qbar) with 2Q = 3P, the class G4 = phi1^*(Q') + phi2^*(phi(Q)') (primes: images on Ea, Eb)
       does not depend on the choice of Q (the map Q -> phi1^*(Q') + phi2^*(phi(Q)') kills E2[2], which is the
       kernel of the gluing), so it is Galois-invariant.  It is computed over the quartic extension Q(t)(x(Q)),
       its Mumford coefficients lie in Q(t), and over Q(t) it has order 4 with 2 G4 = [q, 0] and 2 G4 != 6 D12.
       So G = <G4, D12> = [4,12] is marked over Q(t) and embeds in J_t(Q)_tors at every smooth member.
       (The pullbacks D12 and phi2^*(Pb), Pb the image of P on Eb, of order 4, generate only [2,12].)
   (3) mu_3 inside Eb[3].  The kernel of the dual isogeny E1 -> E2 is phi(E2[3]) = E2[3]/<4P>, which the Weil
       pairing identifies with Hom(<4P>, mu_3) = mu_3.  Checked directly on Eb: the 3-division polynomial of Eb has
       a linear factor x - x0 over Q(t), and -3 y0^2 is a square in Q(t), so the points (x0, +-y0) are defined
       over Q(t)(sqrt -3) and swapped by its conjugation.  If Eb(Q)[3] were nonzero, Eb[3] = Z/3 + mu_3 as a Galois
       module, and Frobenius at a good prime p = 1 (mod 3) would act trivially on Eb[3]; with the rational point 4P
       of Ea this gives J(F_p)[3] of rank >= 3.
   (4) Exactness on a residue class.  For every 7-integral t = 5 (mod 7) the curve C_t has the reduction of C_5
       at 7, with J(F_7) = [12,12].  Reduction is injective on J_t(Q)_tors.
       - The 2-part of J_t(Q)_tors lies in [4,4] and contains <G4, 3 D12> = [4,4], so it is [4,4].
       - The odd part is a 3-group of exponent 3 (#J(F_7) = 144).  J(Q)[3^inf] = Ea(Q)[3^inf] x Eb(Q)[3^inf] by (1);
         7 = 1 (mod 3) and J(F_7)[3] has rank 2, so Eb(Q)[3] = 0 by (3); Ea(Q)[3] is cyclic (Weil pairing).
         So J_t(Q)[3^inf] = <4 D12> = Z/3.
       Hence J_t(Q)_tors is exactly [4,12].  (#J(F_p) = #Ea(F_p)^2 is a square, so a gcd of point counts cannot
       certify this; the group structure is used instead.)  For those t the data x0, y0^2 of (3) also have
       7-integral, nonzero specializations (their denominators do not vanish mod 7 at t = 5).
   (5) Infinitely many.  The G2-invariants of the members differ, so t -> M_2 is nonconstant with finite
       fibres, and the residue class gives infinitely many isomorphism classes.
   (6) Sanity: TorsionSubgroup at t = 1/3 (census 4.12.a), 5, -2 is [4,12].
   Run (from the repository root):
       magma -b pipeline/magma/family_412_split_certificate.m > data/knowledge/sources/family_412_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
tau := t^2 - t; A := 6*tau^2 + 6*tau + 1;
q := tau^2*x^2 - A*x + (3*tau + 1)^2;
qm := tau^2*x^2 + A*x + (3*tau + 1)^2;
s2 := (2*tau + 1)*(4*tau + 1)*x^2 - (6*tau + 1);
f := -A*s2*q*qm;
assert &and[ Coefficient(f, i) eq 0 : i in [1, 3, 5] ] and Degree(f) eq 6 and IsSquarefree(f);
assert #Factorization(f) eq 3;                             // three irreducible quadratics over Q(t)
F := P![ Coefficient(f, 2*i) : i in [0..3] ];            // f = F(x^2)
printf "C_t : y^2 = f_t(x) = %o\n", f;
print "    = -(6tau^2+6tau+1)((2tau+1)(4tau+1)x^2 - (6tau+1)) q(x) q(-x), q(x) = tau^2 x^2 - (6tau^2+6tau+1)x + (3tau+1)^2, tau = t^2 - t";

// ---------------------------------------------------------------- (0) source: the HLP gluing of E2 = E_12(t) and E2/<4P>
mK := (3*t - 3*t^2 - 1)/(t - 1); fK := mK/(1 - t); dK := mK + t; cK := fK*(dK - 1); bK := cK*dK;
E2 := EllipticCurve([Ft| 1-cK, -bK, -bK, 0, 0]);
P12 := E2![0,0];
assert 12*P12 eq E2!0 and 6*P12 ne E2!0 and 4*P12 ne E2!0;
E1, phi := IsogenyFromKernel(E2, x - (4*P12)[1]);
assert Degree(phi) eq 3;
function MonicCubic(E)
    a1, a2, a3, a4, a6 := Explode(aInvariants(E));
    return x^3 + (a2 + a1^2/4)*x^2 + (a4 + a1*a3/2)*x + (a6 + a3^2/4);
end function;
g2 := MonicCubic(E2); g1 := MonicCubic(E1);
fa2 := Factorization(g2);
lin := [ g[1] : g in fa2 | Degree(g[1]) eq 1 ][1]; quad := [ g[1] : g in fa2 | Degree(g[1]) eq 2 ][1];
K<z> := ext< Ft | quad >; PK<X> := PolynomialRing(K);
be := [ K!(-Coefficient(lin, 0)), z, -Coefficient(quad, 1) - z ];      // roots: the x-coordinates of E2[2]
dp := DefiningPolynomials(phi);
a1, a2, a3, a4, a6 := Explode(aInvariants(E2));
function PhiX(xT) yT := -(a1*xT + a3)/2; vals := [ Evaluate(d, [xT, yT, 1]) : d in dp ]; return vals[1]/vals[3]; end function;
al := [ PhiX(b) : b in be ];                                           // x(phi(T_i)): the matching labels of E1[2]
assert &and[ Evaluate(PolynomialRing(K)!g1, a) eq 0 : a in al ] and #Seqset(al) eq 3;
aa1 := (al[3]-al[2])^2/(be[3]-be[2]) + (al[2]-al[1])^2/(be[2]-be[1]) + (al[1]-al[3])^2/(be[1]-be[3]);
bb1 := (be[3]-be[2])^2/(al[3]-al[2]) + (be[2]-be[1])^2/(al[2]-al[1]) + (be[1]-be[3])^2/(al[1]-al[3]);
aa2 := al[1]*(be[3]-be[2]) + al[2]*(be[1]-be[3]) + al[3]*(be[2]-be[1]);
bb2 := be[1]*(al[3]-al[2]) + be[2]*(al[1]-al[3]) + be[3]*(al[2]-al[1]);
Df := ((al[1]-al[2])*(al[1]-al[3])*(al[2]-al[3]))^2;
Dg := ((be[1]-be[2])*(be[1]-be[3])*(be[2]-be[3]))^2;
AA := Dg*aa1/aa2; BB := Df*bb1/bb2;
hHLP := -(AA*(al[2]-al[1])*(al[1]-al[3])*X^2 + BB*(be[2]-be[1])*(be[1]-be[3]))
        *(AA*(al[3]-al[2])*(al[2]-al[1])*X^2 + BB*(be[3]-be[2])*(be[2]-be[1]))
        *(AA*(al[1]-al[3])*(al[3]-al[2])*X^2 + BB*(be[1]-be[3])*(be[3]-be[2]));
hq := P![ Ft!c : c in Coefficients(hHLP) ];
ratio := LeadingCoefficient(hq)/LeadingCoefficient(f);
assert hq eq ratio*f;
okr, sqr := IsSquare(ratio);
assert okr;
print "(0) the HLP (Prop. 4) gluing of E2 = E_12(t) and E1 = E2/<4P> along phi on 2-torsion has coefficients in Q(t) and is (sqrt ratio)^2 f_t";

// ---------------------------------------------------------------- (1) the two elliptic quotients
function MonEC(G)
    cs := [ Coefficient(G, i) : i in [0..3] ];
    return EllipticCurve([0, cs[3], 0, cs[2]*cs[4], cs[1]*cs[4]^2]), cs[4];
end function;
Fr := P![ Coefficient(F, 3 - i) : i in [0..3] ];
Ea, ca := MonEC(F); Eb, cb := MonEC(Fr);
oka, ma := IsIsomorphic(E2, Ea);
okb, mb := IsIsomorphic(E1, Eb);
assert oka and okb;
print "(1) phi1 = (x^2, y) : C_t -> Ea : Y^2 = F(X), phi2 = (1/x^2, y/x^3) : C_t -> Eb : Y^2 = X^3 F(1/X);";
print "    Ea is Q(t)-isomorphic to E2 = E_12(t) and Eb to E1 = E2/<4P>: J is isogenous over Q to E1 x E2";

// ---------------------------------------------------------------- (2) the marked group over Q(t)
J := Jacobian(HyperellipticCurve(f));
// pullbacks of points (X, Y) of the Weierstrass models Ea, Eb (scaled by ca, cb) to the Jacobian
function Pull1(JJ, Pt, c) PP := Parent((JJ!0)[1]); xx := PP.1; XX := Pt[1]/c; YY := Pt[2]/c; return elt< JJ | xx^2 - XX, PP!YY >; end function;
function Pull2(JJ, Pt, c) PP := Parent((JJ!0)[1]); xx := PP.1; XX := Pt[1]/c; YY := Pt[2]/c; return elt< JJ | xx^2 - 1/XX, (YY/XX)*xx >; end function;
D12 := Pull1(J, ma(P12), ca);
assert 12*D12 eq J!0 and 6*D12 ne J!0 and 4*D12 ne J!0;
Pb := mb(E1![ Evaluate(d, Eltseq(P12)) : d in dp ]);                // the image of P on Eb, of order 4
D4b := Pull2(J, Pb, cb);
assert 4*D4b eq J!0 and 2*D4b ne J!0 and 2*D4b eq 6*D12;
assert D4b ne 3*D12 and D4b ne 9*D12;                    // so <D12, phi2^*(Pb)> = [2,12]: the pullbacks alone give only [2,12]
// the gain: Q with 2Q = T4 = 3P, over the quartic field Q(t)(x(Q))
T4 := 3*P12;
b2, b4, b6, b8 := Explode(bInvariants(E2));
quart := x^4 - b4*x^2 - 2*b6*x - b8 - T4[1]*(4*x^3 + b2*x^2 + 2*b4*x + b6);   // x(2Q) = x(T4)
fq := Factorization(quart);
assert #fq eq 1 and Degree(fq[1][1]) eq 4;
K1<w> := ext< Ft | fq[1][1] >;
xQ := w;
// y(Q) is linear in x(Q): the tangent at Q passes through -T4 (eliminate yQ^2 with the curve equation)
x4 := K1!T4[1]; y4m := K1!(-T4[2] - a1*T4[1] - a3);
cub := xQ^3 + a2*xQ^2 + a4*xQ + a6; L := a1*xQ + a3; S := 3*xQ^2 + 2*a2*xQ + a4;
yQ := (2*cub - y4m*L + S*(x4 - xQ))/(L + 2*y4m + a1*(x4 - xQ));
E2K := BaseChange(E2, K1); E1K := BaseChange(E1, K1);
EaK := BaseChange(Ea, K1); EbK := BaseChange(Eb, K1);
Qpt := E2K![xQ, yQ];
assert 2*Qpt eq E2K!Eltseq(T4);
vals := [ Evaluate(d, [Qpt[1], Qpt[2], K1!1]) : d in dp ];
phiQ := E1K![ vals[1]/vals[3], vals[2]/vals[3] ];
maK := Isomorphism(E2K, EaK, [ K1!c : c in IsomorphismData(ma) ]);
mbK := Isomorphism(E1K, EbK, [ K1!c : c in IsomorphismData(mb) ]);
JK := BaseChange(J, K1);
Gk := Pull1(JK, maK(Qpt), K1!ca) + Pull2(JK, mbK(phiQ), K1!cb);
assert &and[ IsCoercible(Ft, c) : c in Coefficients(Gk[1]) cat Coefficients(Gk[2]) ];
G4 := elt< J | P![ Ft!c : c in Coefficients(Gk[1]) ], P![ Ft!c : c in Coefficients(Gk[2]) ] >;
assert 4*G4 eq J!0 and 2*G4 ne J!0;
assert 2*G4 eq elt< J | q, 0 >;
assert 2*G4 ne 6*D12;                        // <G4> meets <D12> trivially
printf "(2) over Q(t): D12 = phi1^*(Pa) has order 12; phi2^*(Pb) has order 4, 2 phi2^*(Pb) = 6 D12, phi2^*(Pb) not in <D12>: together [2,12];\n";
printf "    the gain G4 = phi1^*(Q') + phi2^*(phi(Q)'), 2Q = 3P, computed over Q(t)(x(Q)) (quartic), lies in J(Q(t)):\n";
printf "    G4 = (%o, %o)\n", G4[1], G4[2];
print "    G4 has order 4, 2 G4 = [q, 0] != 6 D12, so G = <G4, D12> = [4,12] is marked over Q(t)";

// ---------------------------------------------------------------- (3) mu_3 inside Eb[3]
psi3b := DivisionPolynomial(Eb, 3);
lins := [ g[1] : g in Factorization(psi3b) | Degree(g[1]) eq 1 ];
assert #lins eq 1;
x0 := -Coefficient(lins[1], 0)/Coefficient(lins[1], 1);
_, ab2, _, ab4, ab6 := Explode(aInvariants(Eb));
y0sq := x0^3 + ab2*x0^2 + ab4*x0 + ab6;
okmu, w0 := IsSquare(-3*y0sq);
assert okmu and y0sq ne 0;
printf "(3) the 3-division polynomial of Eb has the single linear factor x - x0 over Q(t), and -3 y0^2 = (%o)^2:\n", w0;
print "    the kernel of the 3-isogeny Eb -> Ea is mu_3";

// ---------------------------------------------------------------- (4) exactness on the residue class t = 5 mod 7
p := 7; r0 := 5;
members := [Q| 1/3, 5, -2, 12, 19];
fbar := ReductionClass(f, p, r0);
for m in members do
    assert GF(p)!m eq GF(p)!r0;
    assert PolynomialRing(GF(p))!SpecPoly(f, m) eq fbar;
end for;
for c in [x0, y0sq, w0] do
    Nc, Dc := IntNumDen(c);
    assert GF(p)!Evaluate(Dc, r0) ne 0 and GF(p)!Evaluate(Nc, r0) ne 0;
end for;
n, chi := FrobData(fbar);
inv := Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
assert inv eq [12, 12] and n eq 144 and p mod 3 eq 1;
printf "(4) p = %o, t = %o mod %o: reduction %o (squarefree, degree %o); J(F_%o) = %o, chi = %o\n", p, r0, p, fbar, Degree(fbar), p, inv, chi;
print "    x0, y0^2 and w0 have 7-integral nonzero specializations on the class";
print "    the 2-part of J_t(Q)_tors is <G4, 3 D12> = [4,4]; the odd part is a 3-group of exponent 3";
print "    7 = 1 mod 3 and J(F_7)[3] has rank 2, so Eb(Q)[3] = 0 by (3), and J_t(Q)[3^inf] = <4 D12>";
print "    hence every 7-integral t = 5 mod 7 has J_t(Q)_tors exactly [4,12]";

// ---------------------------------------------------------------- (5) nonconstant moduli map, (6) sanity
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } eq #members;
printf "(5) the G2-invariants of the members %o are pairwise distinct\n", members;
C412a := HyperellipticCurve(RX![9, -27, 53, -62, 53, -27, 9], RX![0, 1, 1]);
assert G2Invariants(HyperellipticCurve(SpecPoly(f, 1/3))) eq G2Invariants(C412a);
print "    the member t = 1/3 is the census curve 4.12.a (equal G2-invariants)";
for m in [Q| 1/3, 5, -2] do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(6) TorsionSubgroup at t = %o: %o\n", m, invm;
    assert invm eq [4, 12];
end for;
printf "time %o s\n", Cputime(T0);
print "FAMILY_412_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E1 x E2) and J(Q)_tors exactly [4,12]";
quit;
