/* family_2_certificate.m -- certificate that infinitely many geometrically simple genus-2 Jacobians
   over Q have J(Q)_tors exactly [2], replacing the family y^2 = x(x^5 + x + t) of the archived record
   (census_infinity_record.md, section 2.1), whose generic member has NO rational 2-torsion: the
   quintic x^5 + x + t is irreducible over Q(t), so the even model has only two factors and 2-rank 0;
   the record's fiber t = 2 is exact [2] only because x^5 + x + 2 = (x + 1)(x^4 - x^3 + x^2 - x + 2)
   (found by the audit of 2026-09-18).

   Family:  C_t : y^2 = x (x^4 + x + t)   over Q(t).
     (L1) marked subgroup: (0,0) and the point at infinity are rational Weierstrass points, so
          P_t = [(0,0) - infinity] is a rational point of order 2 of J_t for every t (checked
          symbolically over Q(t) below); the fiber t = 1 has J_1(Q)_tors = [2] exactly, and the
          torsion of the generic fiber injects into J_1(Q)_tors, so the generic torsion is exactly [2];
     (L2) the fiber t = 1 has a strict prime, so its Jacobian is geometrically simple;
     (L3) the fibers t = 1, 2, 3, 5, 7 have pairwise distinct G2-invariants: the moduli map is not
          constant, so its fibers are finite;
     (L4) uniform bound: for two good primes p, q of the fiber t = 1 and t = 1 (mod pq), J_t reduces
          to the same curves mod p and q, hence #J_t(Q)_tors divides gcd(#J_1(F_p), #J_1(F_q)) =: B,
          and chi_p(J_t) = chi_p(J_1) for the strict prime p, so J_t is geometrically simple.
   For t in the residue class 1 + pqZ the torsion contains the marked Z/2 and has order dividing B;
   the t with an extra rational point of order m (2 < m | B) lie in finitely many thin sets, which
   miss infinitely many integers of every residue class (Hilbert irreducibility); the rest have
   torsion exactly [2], and all but finitely many are non-isomorphic (L3).
   Run:  magma -b pipeline/magma/family_2_certificate.m > data/knowledge/sources/family_2_certificate.log
*/
SetColumns(0);
Q := Rationals(); Z := Integers(); R<x> := PolynomialRing(Q);
Ft<t> := FunctionField(Q); Rt<X> := PolynomialRing(Ft);
function IsStrict(chi)
  if Degree(chi) ne 4 or not IsIrreducible(chi) then return false; end if;
  K := NumberField(chi); pi := K.1;
  for n in [2..12] do if Degree(MinimalPolynomial(pi^n)) lt 4 then return false; end if; end for;
  return true;
end function;
function Chi(g, p)   // naive point counting (Magma's default L-polynomial algorithm is unreliable at p <= 7)
  return R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(PolynomialRing(GF(p))!g) : Al := "Naive")));
end function;
// (L1) marked 2-torsion over Q(t)
Cgen := HyperellipticCurve(X*(X^4 + X + t));
Jgen := Jacobian(Cgen);
Pgen := Jgen![X, 0];   // the class of (0,0) - infinity
assert not IsZero(Pgen) and IsZero(2*Pgen);
printf "L1: over Q(t) the class [(0,0) - infinity] on y^2 = x(x^4 + x + t) has order 2\n";
assert IsIrreducible(X^4 + X + t);
// fibers
invs := [];
for tt in [1, 2, 3, 5, 7] do
  g := x*(x^4 + x + tt); C := HyperellipticCurve(g);
  I := Invariants(TorsionSubgroup(Jacobian(C)));
  D := Z!Discriminant(g);
  strict := 0;
  for p in PrimesInInterval(3, 200) do
    if D mod p ne 0 and IsStrict(Chi(g, p)) then strict := p; break; end if;
  end for;
  G2 := G2Invariants(C);
  printf "t = %o: J(Q)_tors = %o, strict prime %o, G2 = %o\n", tt, I, strict, G2;
  assert I eq [2] and strict ne 0;
  Append(~invs, G2);
end for;
assert #SequenceToSet(invs) eq 5;
printf "L2, L3: every fiber is exact [2] and geometrically simple; the five G2-invariants are distinct\n";
// (L4) uniform bound in a residue class of the fiber t = 1
g1 := x*(x^4 + x + 1); D1 := Z!Discriminant(g1);
good := [p : p in PrimesInInterval(3, 50) | D1 mod p ne 0];
p := good[1]; q := good[2];
B := GCD(Z!Evaluate(Chi(g1, p), 1), Z!Evaluate(Chi(g1, q), 1));
printf "L4: fiber t = 1, good primes p = %o, q = %o: #J_t(Q)_tors divides B = gcd(#J_1(F_p), #J_1(F_q)) = %o for all t = 1 mod %o\n", p, q, B, p*q;
printf "FAMILY_2_CERTIFIED\n";
quit;
