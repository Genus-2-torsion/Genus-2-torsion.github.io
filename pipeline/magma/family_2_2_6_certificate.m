/* family_2_2_6_certificate.m -- certificate that infinitely many pairwise non-isomorphic
   genus-2 curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,2,6].

   The family (torsion_jac/notes/simplefam_2026_10_18.md).  For a cubic h and a constant c, the curve
   y^2 = h^2 - c^2 = (h - c)(h + c) has div(y - h) = 3 oo+ - 3 oo-, so D3 = oo+ - oo- has order 3
   (the mechanism of the BNSS2026 [2,2,2,12] curves).  Take the two cubics
       F2 = (x - p)(x - q)(x - r),   F1 = x (x^2 - (p+q+r) x + (pq+pr+qr)),
   which differ by the constant F1 - F2 = pqr.  Then f = F1 F2 has the rational roots 0, p, q, r and the
   quadratic factor x^2 - (p+q+r) x + (pq+pr+qr), irreducible for generic p, q, r: factor type
   [1,1,1,1,2], rational 2-rank 3.  So J(Q) contains (Z/2)^3 x Z/3 = [2,2,6].
   One-parameter slice: p = 1, q = 2, r = t:
       f_t = x (x^2 - (t+3) x + (3t+2)) (x-1)(x-2)(x-t),   F1 - F2 = 2t.

   Argument.
   (1) Over Q(t): D3 = oo+ - oo- has order 3; T1 = [x(x-1), 0], T2 = [x(x-2), 0], T3 = [x(x-t), 0] have
       order 2 and their 7 nonzero sums are nonzero, so they generate (Z/2)^3; with D3 this marks
       [2,2,6] (order 24).
   (2) Residue certificate (family_cert_lib.m) at t = 3 mod 5 and t = 27 mod 29 (t = -2 mod 145):
       #J(F_5) = 24 = #G already, and the Frobenius polynomial at 29 is strict.  So every 5- and
       29-integral t = -2 mod 145 gives J(Q)_tors exactly [2,2,6] and a geometrically simple Jacobian.
   (3) The G2-invariants of two members differ, so infinitely many isomorphism classes.
   Control: the member t = -2 is Q-isomorphic to the census curve 2.2.6.a (checked below).
   Run (from the repository root):
       magma -b pipeline/magma/family_2_2_6_certificate.m > data/knowledge/sources/family_2_2_6_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// (0) the general identity over Q(p,q,r)
F3<p3,q3,r3> := FunctionField(Q, 3); P3<x3> := PolynomialRing(F3);
G2c := (x3-p3)*(x3-q3)*(x3-r3); G1c := x3*(x3^2 - (p3+q3+r3)*x3 + (p3*q3+p3*r3+q3*r3));
assert G1c - G2c eq p3*q3*r3;
print "over Q(p,q,r): x(x^2 - (p+q+r)x + (pq+pr+qr)) - (x-p)(x-q)(x-r) = pqr";

// (1) the slice p = 1, q = 2, r = t and the marked group over Q(t)
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
F2 := (x-1)*(x-2)*(x-t); F1 := x*(x^2 - (t+3)*x + (3*t+2));
assert F1 - F2 eq 2*t;
f := F1*F2;
assert Degree(f) eq 6 and IsSquarefree(f) and LeadingCoefficient(f) eq 1;
assert [ <Degree(g[1]), g[2]> : g in Factorization(f) ] eq [ <1,1>, <1,1>, <1,1>, <1,1>, <2,1> ];
printf "f_t = %o\n", f;
C := HyperellipticCurve(f); J := Jacobian(C);
// f = h^2 - t^2 with h = (F1 + F2)/2: div(y - h) = 3 oo+ - 3 oo-
h := (F1 + F2)/2; assert f eq h^2 - t^2;
D3 := (C![1, 1, 0]) - (C![1, -1, 0]);
assert D3 ne J!0 and 3*D3 eq J!0;
T := [ elt< J | x*(x-1), 0 >, elt< J | x*(x-2), 0 >, elt< J | x*(x-t), 0 > ];
assert &and[ Ti ne J!0 and 2*Ti eq J!0 : Ti in T ];
sums := { i*T[1] + j*T[2] + k*T[3] : i, j, k in [0,1] };
assert #sums eq 8;
grp := { S + e*D3 : S in sums, e in [0..2] };
assert #grp eq 24;
print "over Q(t): D3 = oo+ - oo- has order 3, [x(x-1),0], [x(x-2),0], [x(x-t),0] generate (Z/2)^3; the marked group has order 24 = [2,2,6]";

// (2), (3) residue certificate at t = 3 mod 5, t = 27 mod 29 (t = -2 mod 145)
g, st := ResidueCertificate(f, [5, 29], [3, 27], [-2, 143, -147, 288], 24, [2,2,6]);
for m in [143, -147] do
    printf "sanity: t = %o has TorsionSubgroup %o\n", m, TorsionInvariants(SpecPoly(f, m));
end for;

// control: t = -2 is the census curve 2.2.6.a
RX2<X2> := PolynomialRing(Q);
Ccensus := HyperellipticCurve(3*X2^5 - 4*X2^3 - X2^2 + X2, X2^2 + 1);   // 2.2.6.a: f = [0,1,-1,-4,0,3], h = [1,0,1]
assert IsIsomorphic(HyperellipticCurve(SpecPoly(f, -2)), Ccensus);
print "control: the member t = -2 is Q-isomorphic to census curve 2.2.6.a";
print "FAMILY_2_2_6_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,2,6]";
quit;
