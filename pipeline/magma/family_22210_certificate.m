/* family_22210_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,2,2,10].

   The family: Elkies' Clebsch-Klein family (Elkies 2024, Prop. 2.11; torsion_jac/notes/
   elkies22210_richelot.md)
       y^2 = x (x - r1^2)(x - r2^2)(x - r3^2)(x - r4^2)(x - r5^2),  sum r_i = 0,  sum r_i^3 = 0,
   i.e. (r1 : ... : r5) on the Clebsch diagonal cubic surface S.  Why the 5-torsion: e1 = e3 = 0 for the
   r_i, so prod (X - r_i) = X^5 + e2 X^3 + e4 X - e5 and prod (x - r_i^2) = x (x^2 + e2 x + e4)^2 - e5^2.
   Hence f = h^2 - e5^2 x with h = x (x^2 + e2 x + e4), div(y - h) = (0,0) + 2 oo+ - 3 oo-, and
   [oo+ - oo-] has order dividing 5.
   One-parameter slice: the tangent section of S at Elkies' point P0 = (1, -8, -7, 5, 9).  For v = v1 + u v2
   in the tangent plane (v1 = (7,0,0,-10,3), v2 = (0,56,0,-17,-39)), the line P0 + lambda v meets S again at
       r(u) = (sum v_i^3) P0 - 3 (sum P0_i v_i^2) v,
   a rational curve on S (the tangent section is a nodal cubic with node P0).

   Argument.
   (1) Over Q(u): r(u) lies on S, f_u is squarefree of degree 6 and monic, so oo+, oo- are rational, and
       D5 = [oo+ - oo-] has order 5.  The six rational roots give J[2] = (Z/2)^4: the classes
       T_i = [x (x - r_i^2), 0], i = 1..4, have all 15 nonzero sums nonzero.  So (Z/2)^4 x Z/5 =
       Z/2 x Z/2 x Z/2 x Z/10 is marked over Q(u).
   (2) Residue certificate (family_cert_lib.m) at u = 7 mod 11: #J(F_11) = 80 = #G, strict.  So every
       11-integral u = 7 mod 11 gives J(Q)_tors exactly [2,2,2,10] and a geometrically simple Jacobian.
   (3) The G2-invariants of two members differ, so infinitely many isomorphism classes.
   Run (from the census root):
       magma -b pipeline/magma/family_22210_certificate.m > data/knowledge/sources/family_22210_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// (1) the rational curve on the Clebsch surface and the family over Q(u)
P0 := [Q| 1, -8, -7, 5, 9];
assert &+P0 eq 0 and &+[ c^3 : c in P0 ] eq 0;
v1 := [Q| 7, 0, 0, -10, 3]; v2 := [Q| 0, 56, 0, -17, -39];
for v in [v1, v2] do       // tangent directions at P0
    assert &+v eq 0 and &+[ P0[i]^2*v[i] : i in [1..5] ] eq 0;
end for;
assert Rank(Matrix(Q, [P0, v1, v2])) eq 3;
Fu<u> := FunctionField(Q); P<x> := PolynomialRing(Fu);
v := [ v1[i] + u*v2[i] : i in [1..5] ];
S3 := &+[ v[i]^3 : i in [1..5] ];
S12 := &+[ P0[i]*v[i]^2 : i in [1..5] ];
r := [ S3*P0[i] - 3*S12*v[i] : i in [1..5] ];
assert &+r eq 0 and &+[ c^3 : c in r ] eq 0;
printf "r(u) = %o\n", r;
f := x*&*[ x - c^2 : c in r ];
assert Degree(f) eq 6 and IsMonic(f) and Discriminant(f) ne 0;
// the contact identity behind the 5-torsion
e2 := &+[ r[i]*r[j] : i, j in [1..5] | i lt j ];
e4 := &+[ r[i]*r[j]*r[k]*r[l] : i, j, k, l in [1..5] | i lt j and j lt k and k lt l ];
e5 := &*r;
assert f eq (x*(x^2 + e2*x + e4))^2 - e5^2*x;
C := HyperellipticCurve(f); J := Jacobian(C);
oo := PointsAtInfinity(C);
assert #oo eq 2;
D5 := J![oo[1], oo[2]];
assert D5 ne J!0 and 5*D5 eq J!0;
T := [ elt< J | x*(x - r[i]^2), 0 > : i in [1..4] ];
assert &and[ 2*T[i] eq J!0 : i in [1..4] ];
sums := [ &+[ T[i] : i in S ] : S in Subsets({1..4}) | #S gt 0 ];
assert #Seqset(sums) eq 15 and J!0 notin sums;
print "over Q(u): r(u) lies on the Clebsch surface, f = h^2 - e5^2 x, [oo+ - oo-] has order 5 and the root classes generate (Z/2)^4: [2,2,2,10] is marked";

// (2), (3) residue certificate at u = 7 mod 11
g, st := ResidueCertificate(f, [11], [7], [7, 18, -4], 80, [2,2,2,10]);
print "FAMILY_22210_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,2,2,10]";
quit;
