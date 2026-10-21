// geometric_26_39_checks.m
//
// Supplementary checks for the census proof "No rational points of order 26 or
// 39 on geometrically split abelian surfaces" (21 October 2026).  The proof does
// not depend on this script; it confirms two hand computations and tests two
// statements of the proof on a known example.
//
//   (A) Section 4 (the CM case).  For q = 2, 3 enumerate every integer pair
//       (a, b) such that P(T) = T^4 + a T^3 + b T^2 + q a T + q^2 has the form
//       (T^2 - c_1 T + q)(T^2 - c_2 T + q) with c_1, c_2 real in
//       [-2 sqrt(q), 2 sqrt(q)].  This is a necessary condition for P to be the
//       characteristic polynomial of Frobenius of an abelian surface over F_q.
//       Check that P(1) never equals 26 or 52 for q = 3, and that P(1) <= 25
//       (so P(1) is never a multiple of 39) for q = 2.
//
//   (B) Section 11 (remarks).  The census curve 19.b,
//       y^2 + (x^3 + x^2 + 1) y = x^2 + x, has Jacobian J_1(13) (LMFDB
//       169.a.169.1, Sato-Tate group E_6, so Gamma = C_6 in the notation of the
//       proof).  Check: J(Q)_tors = Z/19; at p = 2 the characteristic polynomial
//       of Frobenius has no root in F_13, so J[13]^ss has no one-dimensional
//       constituent; and #J(F_2) is prime to 3, so J[3]^ss has no trivial
//       constituent.  These are the conclusions of Proposition 9.1 and Lemma 10.1
//       for this surface.
//
// Every claimed fact is checked by an assert; the script ends with the marker
// ALL_CHECKS_PASSED.  Run from this directory:
//     magma -b geometric_26_39_checks.m < /dev/null > geometric_26_39_checks.log 2>&1

SetColumns(0);
SetMemoryLimit(16*10^9);
v1, v2, v3 := GetVersion();
printf "Magma V%o.%o-%o\n", v1, v2, v3;

// ---------------------------------------------------------------------------
// (A) Possible values of P(1) for abelian surfaces over F_2 and F_3.
// ---------------------------------------------------------------------------
// With h(X) = X^2 + a X + (b - 2q) we have P(T) = T^2 h(T + q/T), and the roots
// c_1, c_2 of h are real in [-2 sqrt q, 2 sqrt q] iff
//   a^2 - 4(b - 2q) >= 0,  a^2 <= 16 q,  and  h(2 sqrt q), h(-2 sqrt q) >= 0;
// the last condition is 2q + b >= 2|a| sqrt q, i.e. 2q + b >= 0 and
// (2q + b)^2 >= 4 a^2 q.  All conditions are exact integer inequalities.
// They force -2q <= b <= 6q, so the search range below is complete.
function PossibleCounts(q)
  vals := {};
  for a in [-8*q .. 8*q] do
    if a^2 gt 16*q then continue; end if;
    for b in [-2*q - 5 .. 6*q + 5] do
      if a^2 - 4*(b - 2*q) lt 0 then continue; end if;
      if 2*q + b lt 0 or (2*q + b)^2 lt 4*a^2*q then continue; end if;
      Include(~vals, 1 + a + b + q*a + q^2);
    end for;
  end for;
  return vals;
end function;

V3 := PossibleCounts(3);
V2 := PossibleCounts(2);
print "q = 3: possible values of P(1):", Sort(SetToSequence(V3));
print "q = 2: possible values of P(1):", Sort(SetToSequence(V2));
assert 26 notin V3 and 52 notin V3;
assert Max(V3) le 55;                           // Weil bound (1 + sqrt 3)^4 < 56
assert Max(V2) le 25;                           // so 39 never divides #A(F_2)
assert &and[v mod 39 ne 0 : v in V2];
print "CHECK no abelian surface over F_3 has 26 or 52 points; none over F_2 has a multiple of 39 points: ok";

// ---------------------------------------------------------------------------
// (B) The census curve 19.b = X_1(13).
// ---------------------------------------------------------------------------
R<x> := PolynomialRing(Rationals());
C := HyperellipticCurve(x^2 + x, x^3 + x^2 + 1);
print "CURVE", C;
J := Jacobian(SimplifiedModel(C));
TJ := TorsionSubgroup(J);
assert Invariants(TJ) eq [19];
print "J(Q)_tors = Z/19";
// The curve is isomorphic to the LMFDB model of 169.a.169.1.
C2 := HyperellipticCurve(x^5 + x^4, x^3 + x + 1);
assert IsIsomorphic(C, C2);
print "isomorphic to y^2 + (x^3 + x + 1) y = x^5 + x^4 (LMFDB 169.a.169.1)";
assert Discriminant(C) eq -169;                 // the model has good reduction at 2
L2 := LPolynomial(ChangeRing(C, GF(2)));        // L_2(T) = T^4 P_2(1/T)
P2 := R!Reverse(Coefficients(L2));
print "characteristic polynomial of Frobenius at 2:", P2;
assert Degree(P2) eq 4 and Coefficient(P2, 0) eq 4;
P2mod13 := PolynomialRing(GF(13))!P2;
assert #Roots(P2mod13) eq 0;
print "CHECK P_2 has no root mod 13, so J[13]^ss has no one-dimensional constituent: ok";
n2 := Evaluate(L2, 1);
print "#J(F_2) =", n2;
assert n2 eq 19 and n2 mod 3 ne 0;
print "CHECK #J(F_2) is prime to 3, so J[3]^ss has no trivial constituent: ok";

print "ALL_CHECKS_PASSED";
quit;
