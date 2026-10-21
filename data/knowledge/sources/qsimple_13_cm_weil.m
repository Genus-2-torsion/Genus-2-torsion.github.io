// qsimple_13_cm_weil.m
// Certificate for the CM case of the proof that a geometrically split genus-2
// Jacobian over Q with a rational point of order 13 is not simple over Q
// (proofs/qsimple-13-obstruction.html, Section 4).
//
// Checks:
//  (1) Among all Weil polynomials P(T) = T^4 + a T^3 + b T^2 + 2 a T + 4 of abelian
//      surfaces over F_2, the only one with 13 | P(1) is T^4 + 2T^3 + 2T^2 + 4T + 4.
//      The enumeration is exact: P = (T^2 - b1 T + 2)(T^2 - b2 T + 2) with real
//      |b_i| <= 2 sqrt 2, i.e. X^2 + a X + (b-4) has both roots in [-2 sqrt 2, 2 sqrt 2].
//  (2) That polynomial is irreducible, its root field is Q(zeta_12), and a = 2 is nonzero.
//  (3) Q(i) and Q(sqrt -3) contain only 4 and 6 roots of unity; 13 splits in both;
//      y^2 = x^3 - x and y^2 = x^3 + 1 have CM by Z[i] and Z[zeta_3] (conductors 32, 36)
//      and good ordinary reduction at 13.
SetColumns(0);
SetMemoryLimit(16*10^9);
v1, v2, v3 := GetVersion();
printf "Magma V%o.%o-%o\n", v1, v2, v3;

P<T> := PolynomialRing(Integers());

// (1) exact enumeration of Weil polynomials over F_2 (q = 2)
weil := [];
for a in [-6..6] do
  for b in [-20..30] do
    // X^2 + a X + (b-4) has both roots real and in [-2 sqrt 2, 2 sqrt 2]  <=>
    //   a^2 - 4(b-4) >= 0,  b + 4 >= 0,  (b+4)^2 >= 8 a^2,  a^2 <= 32
    if a^2 - 4*(b-4) ge 0 and b + 4 ge 0 and (b+4)^2 ge 8*a^2 and a^2 le 32 then
      Append(~weil, <a, b>);
    end if;
  end for;
end for;
// the search box contains all solutions: |a| <= 5 and -4 <= b <= a^2/4 + 4 <= 12
assert &and[ Abs(w[1]) le 5 and w[2] ge -4 and w[2] le 12 : w in weil ];
// independent floating-point cross-check of the root condition
CC := ComplexField(60);
for a in [-6..6] do
  for b in [-20..30] do
    f := T^4 + a*T^3 + b*T^2 + 2*a*T + 4;
    rts := Roots(f, CC);
    num := &and[ Abs(Abs(r[1]) - Sqrt(CC!2)) lt 10^-20 : r in rts ] and &+[ r[2] : r in rts ] eq 4;
    assert num eq (<a, b> in weil);
  end for;
end for;
values := { 5 + 3*w[1] + w[2] : w in weil };   // P(1) = 1 + a + b + 2a + 4
printf "number of Weil polynomials over F_2 (abelian surfaces): %o\n", #weil;
printf "values of P(1): %o\n", Sort(SetToSequence(values));
assert Max(values) le 33;
hits := [ w : w in weil | (5 + 3*w[1] + w[2]) mod 13 eq 0 ];
printf "Weil polynomials with 13 | P(1): %o\n", [ T^4 + w[1]*T^3 + w[2]*T^2 + 2*w[1]*T + 4 : w in hits ];
assert hits eq [ <2, 2> ];
P0 := T^4 + 2*T^3 + 2*T^2 + 4*T + 4;
assert Evaluate(P0, 1) eq 13;

// (2) the root field
assert IsIrreducible(P0);
K0 := NumberField(P0);
ok := IsIsomorphic(K0, CyclotomicField(12));
assert ok;
print "P0 = T^4 + 2T^3 + 2T^2 + 4T + 4 is irreducible, P0(1) = 13, root field = Q(zeta_12)";
assert Coefficient(P0, 3) eq 2;     // a = 2 <> 0, so P0(T) <> P0(-T)
assert P0 ne Evaluate(P0, -T);
// the imaginary quadratic subfields of Q(zeta_12) are Q(i) and Q(sqrt -3)
subs := [ s[1] : s in Subfields(K0, 2) ];
discs := { Discriminant(MaximalOrder(s)) : s in subs };
printf "quadratic subfields of Q(zeta_12): discriminants %o\n", discs;
assert discs eq { -4, -3, 12 };

// (3) the two CM fields and the curves E_M
for d in [-1, -3] do
  M := QuadraticField(d);
  OM := MaximalOrder(M);
  nmu := #TorsionUnitGroup(OM);
  fac := Factorization(13*OM);
  printf "Q(sqrt %o): %o roots of unity; 13 splits: %o\n", d, nmu, #fac eq 2 and &and[ f[2] eq 1 : f in fac ];
  assert nmu eq (d eq -1 select 4 else 6);
  assert #fac eq 2;
end for;
EI := EllipticCurve([-1, 0]);    // y^2 = x^3 - x
EJ := EllipticCurve([0, 1]);     // y^2 = x^3 + 1
for pair in [ <EI, -4, 32>, <EJ, -3, 36> ] do
  E := pair[1];
  okcm, D := HasComplexMultiplication(E);
  a13 := TraceOfFrobenius(E, 13);
  printf "%o: conductor %o, CM discriminant %o, a_13 = %o\n", aInvariants(E), Conductor(E), D, a13;
  assert okcm and D eq pair[2] and Conductor(E) eq pair[3];
  assert a13 mod 13 ne 0;          // good ordinary reduction at 13
end for;

print "ALL_CHECKS_PASSED";
quit;
