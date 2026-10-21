// No rational point of order 63 on a geometrically split abelian surface over Q that is
// simple over Q (genus-2 torsion census, proof of 21 October 2026).
//
// Consistency check, not part of the proof.  The census curve 21.b,
//     y^2 + (x^3 + x + 1) y = x^5 + 2x^4 + 2x^3 + x^2   (LMFDB 324.a.648.1),
// is geometrically split and simple over Q with torsion Z/21.  This script checks that
//   (a) for every prime 5 <= p < 300 its Frobenius characteristic polynomial equals the one
//       predicted in Section 10 for the twists of 162b1^2 by a cubic character of conductor 9,
//         (X^2 - a_p X + p)^2                         if p = +-1 mod 9,
//         N_{Q(zeta_3)/Q}(X^2 - a_p z X + p z^2)      otherwise (z a primitive cube root of 1),
//       with a_p = a_p(162b1);
//   (b) #J(F_5) = 21 and J(Q)_tors = Z/21.
// So the bound #A(F_5) = 21 of Section 10 is attained by a Jacobian in the census.
// Run as:  magma -b qsimple_63_census_21b.m
// Nothing here depends on GRH or on any unproved hypothesis.

SetColumns(0);
SetMemoryLimit(16*10^9);
v1, v2, v3 := GetVersion();
printf "Magma version V%o.%o-%o\n", v1, v2, v3;

Qx<x> := PolynomialRing(Rationals());
Zt<t> := PolynomialRing(Integers());
E := EllipticCurve("162b1");
assert aInvariants(E) eq [1, -1, 1, -5, 5];

C := HyperellipticCurve(x^5 + 2*x^4 + 2*x^3 + x^2, x^3 + x + 1);
D := Integers()!Discriminant(C);
printf "discriminant of the model of 21.b: %o\n", Factorization(D);
assert forall{p : p in PrimeDivisors(D) | p in {2, 3}};

K3<z> := CyclotomicField(3);
K3T<T> := PolynomialRing(K3);
function Predicted(p)
  ap := TraceOfFrobenius(E, p);
  if p mod 9 in {1, 8} then
    return (t^2 - ap*t + p)^2;
  end if;
  f := (T^2 - ap*z*T + p*z^2) * (T^2 - ap*z^2*T + p*z^4);
  return Zt![Integers()!(Rationals()!c) : c in Coefficients(f)];
end function;

n := 0;
for p in PrimesInInterval(5, 299) do
  Lp := LPolynomial(ChangeRing(C, GF(p)));
  assert Degree(Lp) eq 4;
  Pp := Zt!Reverse(Coefficients(Lp));                    // T^4 L(1/T)
  assert Pp eq Predicted(p);
  n +:= 1;
end for;
printf "(a) Frobenius polynomials agree with the prediction at all %o primes 5 <= p < 300\n", n;

L5 := LPolynomial(ChangeRing(C, GF(5)));
assert Evaluate(L5, 1) eq 21;
assert Zt!Reverse(Coefficients(L5)) eq t^4 - 5*t^2 + 25;
Cs := IntegralModel(SimplifiedModel(C));                 // y^2 = 4f + h^2, isomorphic to C over Q
J := Jacobian(Cs);
TJ := TorsionSubgroup(J);
assert Invariants(TJ) eq [21];
print "(b) #J(F_5) = 21 and J(Q)_tors = Z/21";

print "ALL_CHECKS_PASSED";
quit;
