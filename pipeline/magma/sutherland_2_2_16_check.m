/* sutherland_2_2_16_check.m -- End(J_Qbar) = Z for the geometrically simple [2,2,16] curve of
   A. V. Sutherland (personal communication, November 2026), census curve 2.2.16.b.

   C: y^2 + (x^2 + x)y = 360000x^6 - 11736600x^5 + 128899980x^4 - 496170642x^3 + 166739244x^2
                         + 519681459x - 284225007
   The torsion and the class are certified by verify_lib.m (data/logs/2.2.16.b.log).  Here:
   End(J_Qbar) = Z by the argument of sutherland_2_22_check.m, step 5: at a good prime p where chi_p is
   irreducible and Q(pi^n) = Q(pi) for every n (enough: n <= 30, since tau(pi)/pi would be a root of
   unity in a Galois closure of degree <= 8), End^0(J_{Fpbar}) is the quartic CM field K_p = Q(pi), and
   End^0(J_Qbar) embeds in it by reduction.  Two such primes whose fields are not isomorphic and have no
   common quadratic subfield leave End^0(J_Qbar) = Q.  Control: the RM census curve 2.22.a is not
   certified by the same test.
   Run:  magma -b pipeline/magma/sutherland_2_2_16_check.m < /dev/null > data/knowledge/sources/sutherland_2_2_16_check.log */
SetColumns(0);
SetMemoryLimit(8*10^9);
Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

function AbsolutelyStrict(chi)
  if Degree(chi) ne 4 or not IsIrreducible(chi) then return false; end if;
  K := NumberField(chi); pi := K.1;
  return forall{n : n in [2..30] | Degree(MinimalPolynomial(pi^n)) eq 4};
end function;
function QuadDiscs(K)   // the discriminants of the quadratic subfields of K
  return {Discriminant(Integers(F[1])) : F in Subfields(K, 2)};
end function;
function EndIsZ(Cm, D, label)
  ChiC := func<p | P!Reverse(Coefficients(LPolynomial(ChangeRing(Cm, GF(p)))))>;
  goodC := [p : p in PrimesInInterval(3, 200) | D mod p ne 0];
  absp := [p : p in goodC | AbsolutelyStrict(ChiC(p))];
  printf "%o: absolutely strict primes below 200: %o\n", label, absp;
  for i in [1..#absp] do
    for j in [i+1..#absp] do
      p := absp[i]; q := absp[j];
      Kp := NumberField(ChiC(p)); Kq := NumberField(ChiC(q));
      dp := QuadDiscs(Kp); dq := QuadDiscs(Kq);
      if IsEmpty(dp meet dq) and not IsIsomorphic(Kp, Kq) then
        printf "%o: p = %o, chi = %o, quadratic subfields of discriminant %o\n", label, p, ChiC(p), dp;
        printf "%o: q = %o, chi = %o, quadratic subfields of discriminant %o\n", label, q, ChiC(q), dq;
        printf "%o: K_p and K_q are not isomorphic and have no common quadratic subfield\n", label;
        return true;
      end if;
    end for;
  end for;
  printf "%o: quadratic subfields at the absolutely strict primes: %o\n", label,
    [<p, QuadDiscs(NumberField(ChiC(p)))> : p in absp];
  return false;
end function;

C := HyperellipticCurve(P![-284225007, 519681459, 166739244, -496170642, 128899980, -11736600, 360000], P![0, 1, 1]);
printf "C: %o\n", C;
Cm := ReducedMinimalWeierstrassModel(C);
D := Z!Discriminant(Cm);
printf "reduced minimal model: %o\n", Cm;
printf "discriminant: %o\n", Factorization(D);
printf "torsion: %o\n", Invariants(TorsionSubgroup(Jacobian(IntegralModel(SimplifiedModel(C)))));
printf "End(J_Qbar) = Z: %o\n", EndIsZ(Cm, D, "C");
C0 := ReducedMinimalWeierstrassModel(HyperellipticCurve(P![0,-6,12,-5,9,-3,1], P![0,1,1]));
printf "control 2.22.a (RM by Z[phi]) certified End = Z: %o (expected false)\n", EndIsZ(C0, Z!Discriminant(C0), "2.22.a");
print "SUTHERLAND_2_2_16_CHECK_DONE";
quit;
