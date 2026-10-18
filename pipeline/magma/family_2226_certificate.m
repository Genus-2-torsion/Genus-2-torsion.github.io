/* family_2226_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,2,2,6].

   The family: the record's M(2,2,2,6) chart (census_infinity_record.md section 2.5;
   torsion_jac/notes/m2226_order6_doubling.md)
       y^2 = x (x+2s^2-sn)(x+2s^2+sm-2sn-mn)(x+2s^2+sm-sn-mn)(2x-mn)(2x+4s^2-4sn-mn),
   in the quintic model X = 1/x, Y = y/x^3 of that note:
       Y^2 = L1 L2 L3 L4 L5,  L1 = 1 + (2s^2-sn)X, L2 = 1 + (2s^2+sm-2sn-mn)X, L3 = 1 + (2s^2+sm-sn-mn)X,
       L4 = 2 - mnX, L5 = 2 + (4s^2-4sn-mn)X,
   where P = (0, 2) gives g = [P - oo] of order 6.  One-parameter slice: s = 7, m = 2, n free (the line
   through the record's fiber (s,m,n) = (7,2,1)).

   Argument.
   (1) Over Q(n): g = [X, 2] has order 6 (6g = 0, 2g != 0, 3g != 0).  The five linear factors give
       J[2] rational; the classes Ti = [X - root(Li), 0], i = 1..4, have all 15 nonzero sums nonzero, so
       they generate (Z/2)^4.  So (Z/2)^4 + <2g> = (Z/2)^4 x Z/3 = Z/2 x Z/2 x Z/2 x Z/6 is marked over Q(n).
   (2) Residue certificate (family_cert_lib.m) at n = 5 mod 13 and n = 2 mod 23 (n = 278 mod 299):
       #J(F_13) = 240, #J(F_23) = 432, gcd 48 = #G, strict at 23.  So every 13- and 23-integral
       n = 278 mod 299 gives J(Q)_tors exactly [2,2,2,6] and a geometrically simple Jacobian.
   (3) The G2-invariants of two members differ, so infinitely many isomorphism classes.
   Run (from the census root):
       magma -b pipeline/magma/family_2226_certificate.m > data/knowledge/sources/family_2226_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// (1) the family and the marked group over Q(n)
Fn<n> := FunctionField(Q); P<x> := PolynomialRing(Fn);
s := Fn!7; m := Fn!2;
L := [ 1 + (2*s^2 - s*n)*x, 1 + (2*s^2 + s*m - 2*s*n - m*n)*x, 1 + (2*s^2 + s*m - s*n - m*n)*x,
       2 - m*n*x, 2 + (4*s^2 - 4*s*n - m*n)*x ];
f := &*L;
assert Degree(f) eq 5 and Evaluate(f, 0) eq 4;
// it is the note's sextic with x = 1/X
Px<xx> := PolynomialRing(Fn);
sext := xx*(xx + 2*s^2 - s*n)*(xx + 2*s^2 + s*m - 2*s*n - m*n)*(xx + 2*s^2 + s*m - s*n - m*n)*(2*xx - m*n)*(2*xx + 4*s^2 - 4*s*n - m*n);
assert &and[ Evaluate(sext, 1/x0) eq Evaluate(f, x0)/x0^6 : x0 in [Fn!2, Fn!3, Fn!5, Fn!7, Fn!11, Fn!13, Fn!17] ];
J := Jacobian(HyperellipticCurve(f));
g6 := elt< J | x, 2 >;
assert 6*g6 eq J!0 and 2*g6 ne J!0 and 3*g6 ne J!0;
rts := [ -Coefficient(l, 0)/Coefficient(l, 1) : l in L ];
T := [ elt< J | x - rts[i], 0 > : i in [1..4] ];
assert &and[ 2*T[i] eq J!0 : i in [1..4] ];
sums := [ &+[ T[i] : i in S ] : S in Subsets({1..4}) | #S gt 0 ];
assert #Seqset(sums) eq 15 and J!0 notin sums;
print "over Q(n), s = 7, m = 2: g = [X, 2] has order 6 and the root classes generate (Z/2)^4: [2,2,2,6] is marked";

// (2), (3) residue certificate at n = 5 mod 13, n = 2 mod 23
g, st := ResidueCertificate(f, [13, 23], [5, 2], [278, -21, 577], 48, [2,2,2,6]);
print "FAMILY_2226_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,2,2,6]";
quit;
