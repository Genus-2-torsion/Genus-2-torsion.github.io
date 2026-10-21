// No rational point of order 63 on a geometrically split abelian surface over Q that is
// simple over Q (genus-2 torsion census, proof of 21 October 2026).
//
// This script checks every numerical fact used in the proof:
//   1. the residue classes modulo 7 of the discriminants of the trace-one rotations
//      (Section 7 of the proof);
//   2. the bound #A(F_2) <= 25 < 63 for abelian surfaces over F_2 (Section 4, the CM case);
//   3. the algebraic identities in Q(sqrt(-3)) used in Sections 9 and 10;
//   4. the field L = Q(zeta_9)^+ = Q[x]/(x^3-3x+1), and that Frob_5 generates Gal(L/Q);
//   5. the elliptic curve 162b1: its model, conductor, non-CM, a_5 = 0, and its torsion
//      Z/21 over L (Derickx-Etropolski-van Hoeij-Morrow-Zureick-Brown, Theorem A);
//   6. the Frobenius polynomial and the point count #A(F_5) = 21 of Section 10.
// Run as:  magma -b qsimple_63_arithmetic.m
// Nothing here depends on GRH or on any unproved hypothesis.

SetColumns(0);
SetMemoryLimit(16*10^9);
v1, v2, v3 := GetVersion();
printf "Magma version V%o.%o-%o\n", v1, v2, v3;

Qx<X> := PolynomialRing(Rationals());

// ---------------------------------------------------------------- 1. residues modulo 7
// A rotation of order n in {3,4,6} in PGL_2(Q) has a trace-one representative with
// characteristic polynomial X^2 - X + 1/c, where c = 1, 2, 3 respectively.
F7 := GF(7);
F7y<y> := PolynomialRing(F7);
discs := [Discriminant(X^2 - X + 1/c) : c in [1, 2, 3]];
assert discs eq [-3, -1, -1/3];
red := [F7!Numerator(d) / F7!Denominator(d) : d in discs];
assert red eq [F7!4, F7!6, F7!2];
assert IsSquare(red[1]) and IsSquare(red[3]);           // n = 3, 6: split at 7
assert not IsSquare(red[2]);                             // n = 4: unramified at 7
assert IsIrreducible(y^2 - y + F7!1/F7!2);               // X^2 - X + 1/2 mod 7
// The eigenvalue ratio of a rotation of order n is a primitive n-th root of unity; for
// n = 3, 6 these lie in F_7^* and reduction mod 7 is injective on mu_6 (7 = 1 mod 6).
assert #{x : x in F7 | x ne 0 and x^6 eq 1} eq 6;
// For V_4 in the unramified case the residue of a = R^2 is a non-square, so
// X^2 - a is irreducible mod 7 (the non-squares mod 7 are 3, 5, 6).
assert forall{a : a in [3, 5, 6] | not IsSquare(F7!a) and IsIrreducible(y^2 - F7!a)};
print "1. residues modulo 7: -3 = 4 and -1/3 = 2 are squares, -1 = 6 is not: OK";

// ---------------------------------------------------------------- 2. the CM case at 2
// (1 + sqrt 2)^4 = 17 + 12 sqrt 2 < 63, since 12^2 * 2 < (63 - 17)^2.
assert 12^2 * 2 lt (63 - 17)^2;
// Stronger: exhaustive list of P(1) for q-Weil polynomials T^4+aT^3+bT^2+qaT+q^2, q = 2.
// Write P(T) = (T^2 - b1 T + q)(T^2 - b2 T + q) with b1, b2 the roots of
// h(B) = B^2 + aB + (b - 2q); the condition is that both roots are real and lie in
// [-2 sqrt q, 2 sqrt q].  This is a necessary condition, so the list contains all values.
function WeilValues(q)
  vals := {};
  B := 4*Isqrt(4*q) + 8;
  for a in [-B..B] do
    for b in [-6*q-10..6*q+10] do
      if a^2 - 4*(b - 2*q) lt 0 then continue; end if;      // real roots
      if a^2 gt 16*q then continue; end if;                  // vertex in the interval
      // h(+-2 sqrt q) = 2q + b +- 2a sqrt q >= 0  <=>  2q + b >= 0 and (2q+b)^2 >= 4 a^2 q
      if 2*q + b lt 0 then continue; end if;
      if (2*q + b)^2 lt 4*a^2*q then continue; end if;
      Include(~vals, 1 + a + b + q*a + q^2);
    end for;
  end for;
  return vals;
end function;
V2 := WeilValues(2);
printf "possible #A(F_2) for an abelian surface A over F_2: %o\n", Sort(SetToSequence(V2));
assert Max(V2) eq 25;
assert 63 notin V2;
assert forall{n : n in V2 | n mod 63 ne 0};
print "2. no abelian surface over F_2 has 63 | #A(F_2): OK";

// ---------------------------------------------------------------- 3. identities in Q(sqrt -3)
// Gamma = C_6: R^2 - R + 1/3 = 0, u = 3R - 1, S = u^-2 R, T = 3S.
K6<R6> := NumberField(X^2 - X + 1/3);
u := 3*R6 - 1;
assert u^2 - u + 1 eq 0;                                 // u is a primitive 6th root of 1
assert u^3 eq -1 and u^6 eq 1;
assert R6 eq (1 + u)/3 and R6^2 eq u/3;
S := u^-2 * R6;
assert S^2 eq -1/3;
T := 3*S;
assert T eq u^-2 * (1 + u) and T^2 eq -3;
// T and zeta_3 - 1 = u^2 - 1 generate the same ideal of Z[u] (the prime above 3).
O6 := MaximalOrder(K6);
assert ideal<O6 | O6!T> eq ideal<O6 | O6!(u^2 - 1)>;
assert Norm(ideal<O6 | O6!T>) eq 3;
assert Norm(O6!(1 + u)) eq 3;
// det R^k = 3^-k > 0 for every element of Gamma = C_6.
assert Norm(R6) eq 1/3;
// The eigenvalue ratio of a trace-one rotation of order 3 (resp. 6) has order 3 (resp. 6).
K3a<R3> := NumberField(X^2 - X + 1);
assert R3^6 eq 1 and R3^3 eq -1;
assert (R3/(1 - R3))^3 eq 1 and R3/(1 - R3) ne 1;       // the roots are R3 and 1 - R3
r6 := R6/(1 - R6);
assert r6^6 eq 1 and r6^3 ne 1 and r6^2 ne 1;
// Gamma = C_3: R^2 - R + 1 = 0, so R = u is a primitive 6th root of unity, and
// u^(2k) * u^k = u^(3k) = (-1)^k is a scalar.
assert R3^2 - R3 + 1 eq 0 and R3^3 eq -1;
print "3. identities in Q(sqrt(-3)): u^2-u+1 = 0, S^2 = -1/3, T^2 = -3, (T) = (zeta_3 - 1): OK";

// ---------------------------------------------------------------- 4. the field L
K9<z9> := CyclotomicField(9);
assert MinimalPolynomial(z9 + 1/z9) eq X^3 - 3*X + 1;
L<al> := NumberField(X^3 - 3*X + 1);
assert IsNormal(L);
assert Discriminant(MaximalOrder(L)) eq 81;              // unramified outside 3
F5 := GF(5);
F5y<w5> := PolynomialRing(F5);
assert IsIrreducible(w5^3 - 3*w5 + 1);                   // 5 is inert in L
assert 5 mod 9 notin {1, 8};
print "4. L = Q(zeta_9)^+ = Q[x]/(x^3-3x+1) is cyclic cubic, unramified at 5, and 5 is inert: OK";

// ---------------------------------------------------------------- 5. the curve 162b1
E := EllipticCurve("162b1");
assert aInvariants(E) eq [1, -1, 1, -5, 5];              // y^2 + xy + y = x^3 - x^2 - 5x + 5
assert Conductor(E) eq 162 and Factorization(162) eq [<2, 1>, <3, 4>];
assert Denominator(jInvariant(E)) ne 1;                  // non-integral j: no CM
printf "j(162b1) = %o\n", jInvariant(E);
a5 := TraceOfFrobenius(E, 5);
assert a5 eq 0;
assert #ChangeRing(E, F5) eq 6;                           // 5 + 1 - a_5
EL := BaseChange(E, L);
P := EL![1 - 2*al, 2 - 2*al];
assert Order(P) eq 21;
// The point printed after Theorem A of that paper, (2a-3, 2a-2), is not on this model:
// no point of E(L) has x-coordinate 2a-3.  The proof uses only the torsion statements.
assert not IsPoint(EL, 2*al - 3);
TL, mTL := TorsionSubgroup(EL);
assert Invariants(TL) eq [21];
// Hence E(L) has exactly EulerPhi(21) = 12 points of order 21, i.e. 6 classes {+-P}, which
// give 6 distinct L-points of the fine moduli space Y_1(21), none of them rational.
assert #{t : t in TL | Order(t) eq 21} eq 12;
print "5. 162b1: conductor 162, non-CM, a_5 = 0, and E(L)_tors = Z/21 generated by (1-2a, 2-2a): OK";

// ---------------------------------------------------------------- 6. Frobenius at 5
// V_l(A) = V_l(E_1) (x) (Q(zeta_3) (x) Q_l)(xi^-1) with xi(Frob_5) a primitive cube root of
// unity; the characteristic polynomial of Frob_5 is N(X^2 - a_5 z X + 5 z^2).
K3<z> := CyclotomicField(3);
K3X<Y> := PolynomialRing(K3);
for zeta in [z, z^2] do
  f := (Y^2 - a5*zeta*Y + 5*zeta^2) * (Y^2 - a5*zeta^-1*Y + 5*zeta^-2);
  Pf := Qx![Rationals()!c : c in Coefficients(f)];
  assert Pf eq X^4 - 5*X^2 + 25;
  assert Evaluate(Pf, 1) eq 21;
  assert Norm(1 - a5*zeta + 5*zeta^2) eq 21;
end for;
assert 21 mod 63 ne 0;                                   // Z/63 does not embed in a group of order 21
print "6. #A(F_5) = N(1 + 5 zeta^2) = 21, which is not divisible by 63: OK";

print "ALL_CHECKS_PASSED";
quit;
