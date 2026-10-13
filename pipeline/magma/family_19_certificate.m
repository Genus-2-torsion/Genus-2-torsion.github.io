/* family_19_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [19].

   The family F1 (derived 2026-10-13, torsion_jac/notes/z38_2026_10_13.md, section 3; census curve
   19.a is its member r = -1).  Quintic model C_r : y^2 = f_r with W = infinity and
       lambda = -r^4/(r-1)^3,  S = (lambda (x-1)^3 + x^4)/(x - r),
       A = (S - x(x-1)(x-r))/2,  B = (S + x(x-1)(x-r))/2,
       f_r = A^2 + lambda x (x-1)^4 = B^2 - x^5 (x-1).
   With P0 = (0, A(0)), P1 = (1, A(1)) and D_i = P_i - W, the two expressions for f_r give
   div(y - A) = P0 + 4 P1 - 5W and div(y - B) = 5 P0 + P1 - 6W, a relation lattice of
   determinant -19.

   Argument.
   (1) Over Q(r): f = B^2 - x^5(x-1) holds exactly, and D0 = [x, A(0)] satisfies 19 D0 = 0,
       D0 != 0 (Magma arithmetic over Q(r)).  On every smooth fiber D0 specializes to the class
       of P0 - W with P0 != W, which is nonzero of order dividing 19, so of order 19.  Hence
       J(Q) of C_r contains Z/19 whenever f_r is squarefree of degree 5.
   (2) The coefficients of f_r are in Z[1/2][r, 1/(r-1)].  For a rational r0 with
       r0 = 2 mod 3 (3-integral), r0 - 1 and 2 are 3-adic units, so f_{r0} = f_2 mod 3.
       f_2 mod 3 is squarefree of degree 5, so C_{r0} has good reduction at 3, the same as C_2.
   (3) #J(F_3) = 19 for C_2.  Torsion injects into J(F_3) at a prime of good reduction with
       e = 1 < 3 - 1, so #J(Q)_tors divides 19 and J(Q)_tors = Z/19 by (1).  The Frobenius
       polynomial at 3 is strict (irreducible, pi^n of degree 4 for n = 2..12), so J is
       geometrically simple.
   (4) The G2-invariants at r = 2 and r = 5 differ, so the moduli map r -> M_2 is nonconstant,
       with finite fibers: the r0 = 2 mod 3 give infinitely many pairwise non-isomorphic curves.
   Run:  magma -b pipeline/magma/family_19_certificate.m > data/knowledge/sources/family_19_certificate.log
*/
SetColumns(0);
Q := Rationals(); Z := Integers();
Fr<r> := FunctionField(Q); P<x> := PolynomialRing(Fr);
lam := -r^4/(r-1)^3;
S := (lam*(x-1)^3 + x^4) div (x - r);
assert S*(x - r) eq lam*(x-1)^3 + x^4;
A := (S - x*(x-1)*(x-r))/2; B := (S + x*(x-1)*(x-r))/2;
f := A^2 + lam*x*(x-1)^4;
assert Degree(A) le 2 and Degree(B) eq 3 and IsMonic(B) and Degree(f) eq 5;
assert f eq B^2 - x^5*(x-1);
print "over Q(r): f = A^2 + lambda x(x-1)^4 = B^2 - x^5(x-1)";

// (1) the 19-torsion section
J := Jacobian(HyperellipticCurve(f));
D0 := elt< J | x, Evaluate(A, 0) >;
D1 := elt< J | x - 1, Evaluate(A, 1) >;
assert D0 ne J!0 and 19*D0 eq J!0 and D1 eq -5*D0;
print "over Q(r): D0 = [x, A(0)] has order 19, and D1 = -5 D0";

// (2) the coefficients are 3-integral for r = 2 mod 3, and reduce to those of f_2
den := [ Factorization(Denominator(c)) : c in Coefficients(f) ];
assert &and[ &and[ t[1] eq r - 1 : t in d ] : d in den ];
R<X> := PolynomialRing(Q);
spec := func< g, r0 | R![ Evaluate(c, r0) : c in Coefficients(g) ] >;
f2 := spec(f, 2);
assert &and[ Denominator(c) mod 3 ne 0 : c in Coefficients(f2) ];
for c in Coefficients(f) do      // the numerators have coefficients in Z[1/2]
    assert &and[ Denominator(e) in {1, 2, 4} : e in Coefficients(Numerator(c)) ];
end for;
f2p := PolynomialRing(GF(3))!f2;
assert Degree(f2p) eq 5 and IsSquarefree(f2p);
for r0 in [-1, 5, 2/7, 8/7, 1/2, 11/4] do   // sanity: other members of the residue class reduce to f_2
    assert GF(3)!r0 eq 2;
    assert PolynomialRing(GF(3))!spec(f, r0) eq f2p;
end for;
printf "f_2 = %o; f_r = f_2 mod 3 is squarefree of degree 5 for every 3-integral r = 2 mod 3\n", f2;

// (3) the reduction at 3
Lp := LPolynomial(HyperellipticCurve(f2p) : Al := "Naive");
Tt := R.1;
chi := R!(Tt^4*Evaluate(R!Lp, 1/Tt));
assert Evaluate(Lp, 1) eq 19 and Evaluate(chi, 1) eq 19;
assert IsIrreducible(chi) and &and[ Degree(MinimalPolynomial(NumberField(chi).1^n)) eq 4 : n in [2..12] ];
printf "p = 3: #J(F_3) = %o, chi = %o (strict)\n", Evaluate(Lp, 1), chi;

// sanity: exact torsion of two members
for r0 in [2, 5] do
    C := IntegralModel(SimplifiedModel(HyperellipticCurve(spec(f, r0))));
    T := TorsionSubgroup(Jacobian(C));
    assert Invariants(T) eq [19];
    printf "r = %o: TorsionSubgroup = %o\n", r0, Invariants(T);
end for;

// positive control: r = -1 (also = 2 mod 3) is census curve 19.a
C19a := HyperellipticCurve(X^6 + X^5 + X^4 - 5*X^3 - 3*X^2 + 5*X - 2, X^2 + X + 1);
assert IsIsomorphic(HyperellipticCurve(spec(f, -1)), C19a);
print "r = -1: the fiber is census curve 19.a";

// (4) the moduli map is nonconstant
assert G2Invariants(HyperellipticCurve(spec(f, 2))) ne G2Invariants(HyperellipticCurve(spec(f, 5)));
print "G2-invariants at r = 2 and r = 5 differ";
print "FAMILY_19_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [19]";
quit;
