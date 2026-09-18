/* family_22212_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,2,2,12].

   Ingredients (Balakrishnan-Najman-Shnidman-Sutherland, arXiv:2608.28543, Theorems 3.1 and 3.3):
   for P = [a:b:c:u:v] on the open surface S° (a^2+b^2+c^2 = u^2+v^2, a^4+b^4+c^4 = u^4+v^4, all
   coordinates nonzero with distinct squares) the curve C_P : y^2 = x(x-a^2)(x-b^2)(x-c^2)(x-u^2)(x-v^2)
   has [2,2,2,12] inside Jac(C_P)(Q); the elliptic curve E : y^2 = x^3 - 21x - 20 (conductor 288,
   rank 1, Q = (-3,4) of infinite order) maps to S by psi (below), psi(2Q) = P_1 = [120:143:266:218:241],
   and psi(E) is not contained in the degenerate locus.

   Argument (audit of 2026-09-18; verified here).  Let p in {29, 37, 71}.  The point 2Q has p-integral
   coordinates, the denominators of psi are p-adic units at 2Q, and the five coordinates of psi(2Q)
   are p-adic units.  If n = 2 mod ord(Q mod p) then nQ and 2Q have the same reduction mod p, hence
   congruent coordinates, so psi(nQ) = psi(2Q) as points of P^4(F_p) and C_{psi(nQ)} has the same
   reduction mod p as C_{P_1} (the curve depends on the projective point up to Q-isomorphism).  With
   M = lcm(16, 10, 14) = 560 and n = 2 + 560k this holds at all three primes.  C_{P_1} has good
   reduction at 29, 37, 71 with #J(F_29) = 864, #J(F_37) = 1344, gcd 96, and a strict Frobenius
   polynomial at 71.  So for every k with psi((2+560k)Q) in S° (all but finitely many k, and in fact
   every k, since the reductions are nondegenerate): #J(Q)_tors divides 96 and contains [2,2,2,12] of
   order 96, hence equals it; and J is geometrically simple (strict prime 71).  The G2-invariants at
   psi(2Q) and psi(4Q) differ, so the moduli map E -> M_2 is nonconstant and has finite fibers:
   infinitely many of the curves are pairwise non-isomorphic.
   Run:  magma -b pipeline/magma/family_22212_certificate.m > data/knowledge/sources/family_22212_certificate.log
*/
SetColumns(0);
Q := Rationals(); Z := Integers(); R<x> := PolynomialRing(Q);
E := EllipticCurve([0, 0, 0, -21, -20]);
assert Conductor(E) eq 288;
Qp := E![-3, 4];
assert Order(Qp) eq 0;
function Psi(P)
  xx := P[1]; yy := P[2];
  w := (3*xx + yy - 6)/(yy + 18);
  t := -(4*xx^2 - 7*xx + 12*yy + 16)/((xx - 8)*(3*xx + 2*yy + 12));
  return [3*t*(1-t), 3*t-1, 2*w*(3*t^2+1), 3*t^2-3*t+2, 6*t^2-3*t+1], [yy + 18, (xx - 8)*(3*xx + 2*yy + 12)];
end function;
function Normalise(v)   // projective point -> primitive integer coordinates
  d := LCM([Denominator(c) : c in v]); w := [Z!(c*d) : c in v]; g := GCD(w); return [c div g : c in w];
end function;
function CurveOf(P) return HyperellipticCurve(x*&*[x - c^2 : c in P]); end function;
function OnSurface(P) return P[1]^2+P[2]^2+P[3]^2 eq P[4]^2+P[5]^2 and P[1]^4+P[2]^4+P[3]^4 eq P[4]^4+P[5]^4; end function;
function IsStrict(chi)
  if Degree(chi) ne 4 or not IsIrreducible(chi) then return false; end if;
  K := NumberField(chi); pi := K.1;
  for n in [2..12] do if Degree(MinimalPolynomial(pi^n)) lt 4 then return false; end if; end for;
  return true;
end function;
P2 := 2*Qp; v, dens := Psi(P2);
P1 := Normalise(v); assert P1 eq [120, 143, 266, 218, 241] and OnSurface(P1);
C1 := CurveOf(P1); g1 := HyperellipticPolynomials(C1); D1 := Z!Discriminant(g1);
assert Invariants(TorsionSubgroup(Jacobian(C1))) eq [2, 2, 2, 12];
printf "psi(2Q) = P_1 = %o, J(Q)_tors = [2,2,2,12]\n", P1;
counts := [];
for p in [29, 37, 71] do
  assert D1 mod p ne 0;                                            // good reduction of C_{P_1}
  assert &and[Valuation(c, p) ge 0 : c in [P2[1], P2[2]]];         // 2Q p-integral
  assert &and[Valuation(c, p) eq 0 : c in dens];                   // denominators of psi are units at 2Q
  assert &and[Valuation(c, p) eq 0 : c in v];                      // coordinates of psi(2Q) are units
  o := Order(ChangeRing(E, GF(p))![-3, 4]);
  chi := R!Reverse(Coefficients(LPolynomial(ChangeRing(C1, GF(p)) : Al := "Naive")));
  Append(~counts, Z!Evaluate(chi, 1));
  printf "p = %o: ord(Q mod p) = %o, #J(F_p) = %o, chi = %o%o\n", p, o, Evaluate(chi, 1), chi, IsStrict(chi) select " (strict)" else "";
  if p eq 71 then assert IsStrict(chi); end if;
  assert o eq [16, 10, 14][Position([29, 37, 71], p)];
end for;
assert GCD(counts[1], counts[2]) eq 96 and LCM([16, 10, 14]) eq 560;
C2 := CurveOf(Normalise(Psi(4*Qp)));
assert G2Invariants(C1) ne G2Invariants(C2);
printf "gcd(#J(F_29), #J(F_37)) = 96 = #[2,2,2,12]; G2-invariants at psi(2Q) and psi(4Q) differ\n";
// sanity check of the congruence at k = 1
P562 := Normalise(Psi(562*Qp)); assert OnSurface(P562);
for p in [29, 37, 71] do
  i := [j : j in [1..5] | P1[j] mod p ne 0][1]; lam := (GF(p)!P562[i]) / GF(p)!P1[i];
  assert &and[GF(p)!P562[j] eq lam*GF(p)!P1[j] : j in [1..5]];
end for;
printf "psi(562 Q) = psi(2Q) in P^4(F_p) for p = 29, 37, 71 (sanity check of the congruence)\n";
printf "FAMILY_22212_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,2,2,12]\n";
quit;
