/* sutherland_2_22_check.m -- the [2,22] curve without real multiplication (A. V. Sutherland, personal
   communication, November 2026), cited in the note on the simple [2,22] row of the infinity column.

   C: y^2 + (x^2 + x)y = 9x^6 + 111x^5 + 279x^4 - 485x^3 - 282x^2 + 1836x + 1456   ([[f],[h]] as sent)
   Sutherland reports J(Q)_tors = [2,22] and End(J_Qbar) = Z.  Checked here:
     1. J(Q)_tors (Magma's TorsionSubgroup, on an integral simplified model);
     2. geometric simplicity: a strict prime (chi irreducible, no degree drop of pi^n for n <= 12,
        the criterion of verify_lib.m);
     3. no real multiplication defined over Q: if chi_p is irreducible, then J mod p is simple with
        End^0 = Q(pi), so End_Q(J) (x) Q, which commutes with Frobenius, embeds in Q(pi).  A real
        quadratic subfield of the quartic CM field Q(pi) is its maximal real subfield
        Q(pi + p/pi) = Q(sqrt(disc)), disc = c1^2 - 4(c2 - 2p) for chi = T^4 + c1 T^3 + c2 T^2 + p c1 T + p^2.
        Two good primes with irreducible chi and different squarefree parts of disc exclude RM over Q;
     4. the conductor (Magma's Conductor; its exponent at 2 uses Ogg's formula, without a guarantee
        when v_2(disc) >= 12) and the rigorous comparison with the census curve 2.22.a (19044) through
        the odd part.
   End(J_Qbar) = Z itself is not certified here; the Frobenius data are only consistent with it.
   Run:  magma -b pipeline/magma/sutherland_2_22_check.m < /dev/null > data/knowledge/sources/sutherland_2_22_check.log */
SetColumns(0);
SetMemoryLimit(8*10^9);
Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);
f := P![1456, 1836, -282, -485, 279, 111, 9];
h := P![0, 1, 1];
C := HyperellipticCurve(f, h);
printf "C: %o\n", C;

// 1. torsion
T := TorsionSubgroup(Jacobian(IntegralModel(SimplifiedModel(C))));
printf "torsion: %o\n", Invariants(T);

Cm := ReducedMinimalWeierstrassModel(C);
D := Z!Discriminant(Cm);
printf "reduced minimal model: %o\n", Cm;
printf "discriminant: %o = %o\n", D, Factorization(D);

function Chi(p)
  return P!Reverse(Coefficients(LPolynomial(ChangeRing(Cm, GF(p)))));
end function;

function IsStrict(chi)
  if Degree(chi) ne 4 or not IsIrreducible(chi) then return false; end if;
  K := NumberField(chi); pi := K.1;
  return forall{n : n in [2..12] | Degree(MinimalPolynomial(pi^n)) eq 4};
end function;

function SquarefreePart(d)
  return Sign(d) * &*[Z | t[1]^(t[2] mod 2) : t in Factorization(d)];
end function;

// 2. strict prime
good := [p : p in PrimesInInterval(3, 100) | D mod p ne 0];
strict := [p : p in good | IsStrict(Chi(p))];
printf "strict primes below 100: %o\n", strict;
printf "geometrically simple: %o\n", #strict gt 0;

// 3. real quadratic field Q(pi + p/pi) at the good primes with irreducible chi
fields := [];
for p in good do
  chi := Chi(p);
  if not IsIrreducible(chi) then printf "p = %o: chi reducible\n", p; continue; end if;
  c1 := Z!Coefficient(chi, 3); c2 := Z!Coefficient(chi, 2);
  s := SquarefreePart(c1^2 - 4*(c2 - 2*p));
  printf "p = %o: chi = %o, Q(pi + p/pi) = Q(sqrt(%o))\n", p, chi, s;
  Append(~fields, <p, s>);
end for;
quad := {t[2] : t in fields | t[2] ne 1};
printf "distinct real quadratic fields: %o\n", #quad;
printf "no RM over Q: %o\n", #quad ge 2;

// 4. conductor
N := Conductor(Cm);
odd := &*[Z | t[1]^t[2] : t in Factorization(N) | t[1] ne 2];
printf "conductor: %o = %o\n", N, Factorization(N);
printf "odd part %o > 19044 (conductor of 2.22.a): %o\n", odd, odd gt 19044;
print "SUTHERLAND_2_22_CHECK_DONE";
quit;
