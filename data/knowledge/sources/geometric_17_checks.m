// geometric_17_checks.m
// Certificate for the computations used in the proof that no geometrically split abelian
// surface over Q has a rational point of order 17 (proofs/geometric-17-obstruction.html).
//
// Checks:
//  (A) Weil polynomials P(T) = T^4 + a T^3 + b T^2 + 2 a T + 4 of abelian surfaces over F_2:
//      exact enumeration; no value P(1) is divisible by 17 (Section 4).
//  (B) Residues modulo 17: -3 and -1/3 are non-squares, -1 = 4^2 is a square, the roots of
//      X^2 - X + 1/2 are 11 and 7 with ratio of order 4, and F_17^* has order 16
//      (Sections 5 and 6).
//  (C) The modular curve X_0(34) (Sections 6 and 7):
//      (C0) Phi_2 and Phi_17 as returned by ClassicalModularPolynomial are the classical
//           modular polynomials: Phi_N is monic of degree N+1 in the second variable, of
//           degree N+1 in the first, and Phi_N(j(q), j(q^N)) = O(q^m) with m > (N+1)^2,
//           where j(q) = E_4^3/Delta is computed from the definitions;
//      (C1) Magma's SmallModularCurve(34) is a nonsingular plane quartic X with
//           #X(F_p) = p + 1 - tr(T_p | S_2(34)) for all primes 3 <= p < 100, p <> 17
//           (a consistency check, not used in the proof);
//      (C2) identification of (X, j, w2) with (X_0(34), j, w_2):
//           (a) the function j of the database has degree 54 = [PSL_2(Z) : Gamma_0(34)];
//               w2 : (x:y:z) -> (z:-y:x) and w17 : (x:y:z) -> (x:-y:z) are involutions of X;
//               Phi_2(j, j o w2) = 0 in the function field of X;
//           (b) a cross-check with Magma's jInvariant on X at 5 points over number fields;
//           (c) Phi_17(j, j o w17) = 0 in the function field of X: the function vanishes at
//               more than 1944 = 18*54 + 18*54 geometric points of X at which j and j o w17
//               are regular, and its polar divisor has degree at most 1944;
//      (C3) the quotient: a nonsingular curve Yq of genus 1 and a morphism pi : X -> Yq of
//           degree 2 with pi o w2 = pi; Yq is isomorphic over Q to the elliptic curve
//           [1,-1,1,-6,-4] (17a2), whose torsion subgroup is Z/2 x Z/2 and whose 2-Selmer
//           group has order 4, so that its rank is 0; hence Yq has exactly 4 rational points;
//      (C4) the fibres of pi above them: two fibres consist of two rational points each,
//           and these four points are the poles of j (the cusps), of orders 34, 17, 2, 1;
//           each of the other two fibres is one pair of conjugate points defined over
//           Q(sqrt -15), at which j is not an algebraic integer.
//      Consequently a non-cuspidal point x of X_0(34) over a quadratic field K with
//      w_2(x) equal to the conjugate of x is defined over Q(sqrt -15), so K = Q(sqrt -15).
// Every value quoted in the text of the proof is asserted below, not only printed.
SetColumns(0);
SetMemoryLimit(16*10^9);
v1, v2, v3 := GetVersion();
printf "Magma V%o.%o-%o\n", v1, v2, v3;

// ---------------------------------------------------------------------------------------
// (A) Weil polynomials over F_2
P<T> := PolynomialRing(Integers());
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
printf "number of polynomials T^4 + a T^3 + b T^2 + 2 a T + 4 with all roots of absolute value sqrt 2: %o\n", #weil;
assert #weil eq 35;
printf "values of P(1): %o\n", Sort(SetToSequence(values));
assert Sort(SetToSequence(values)) eq [1..16] cat [19, 20, 25];
assert &and[ v mod 17 ne 0 : v in values ];
print "CHECK (A): no abelian surface over F_2 has a number of points divisible by 17";

// ---------------------------------------------------------------------------------------
// (B) residues modulo 17
F17 := GF(17);
assert not IsSquare(F17!-3);
assert not IsSquare(F17!-1/F17!3);
assert F17!4^2 eq F17!-1;
assert Order(F17!3) eq 16;        // F_17^* is cyclic of order 16
assert IsIrreducible(PolynomialRing(F17)![1, -1, 1]);        // X^2 - X + 1
assert IsIrreducible(PolynomialRing(F17)![1/F17!3, -1, 1]);  // X^2 - X + 1/3
rts17 := { t[1] : t in Roots(PolynomialRing(F17)![1/F17!2, -1, 1]) };   // X^2 - X + 1/2
assert rts17 eq { F17!11, F17!7 };
assert Order(F17!11/F17!7) eq 4;
assert F17!11/F17!7 eq F17!4;
print "CHECK (B): -3 and -1/3 are non-squares modulo 17, -1 = 4^2, X^2 - X + 1/2 has the roots 11 and 7 with ratio of order 4, and F_17^* has order 16";

// ---------------------------------------------------------------------------------------
// (C0) the modular polynomials
function CheckModularPolynomial(N)
  Phi := ClassicalModularPolynomial(N);
  assert Rank(Parent(Phi)) eq 2;
  assert BaseRing(Parent(Phi)) cmpeq Integers();
  assert Degree(Phi, 1) eq N + 1 and Degree(Phi, 2) eq N + 1;
  assert Coefficient(Phi, 2, N + 1) eq 1;       // monic in the second variable
  bound := (N + 1)^2;
  prec := 2*bound + 60;
  Lq<q> := LaurentSeriesRing(Rationals(), prec);
  // j = E_4^3/Delta, with Delta = q*prod(1-q^n)^24 and Euler's pentagonal number theorem
  eta := Lq!0;
  k := 0;
  while k*(3*k - 1) div 2 lt prec do
    eta +:= (-1)^k*q^(k*(3*k - 1) div 2);
    if k gt 0 and k*(3*k + 1) div 2 lt prec then
      eta +:= (-1)^k*q^(k*(3*k + 1) div 2);
    end if;
    k +:= 1;
  end while;
  eta +:= O(q^prec);
  Delta := q*eta^24;
  E4 := 1 + 240*&+[ DivisorSigma(3, n)*q^n : n in [1..prec-1] ] + O(q^prec);
  jq := E4^3/Delta;
  p0 := AbsolutePrecision(jq);
  assert p0 ge 2*bound;
  assert Valuation(jq) eq -1 and Coefficient(jq, -1) eq 1 and Coefficient(jq, 0) eq 744 and Coefficient(jq, 1) eq 196884;
  jm := jInvariant(q);                          // Magma's q-expansion, for comparison
  assert IsWeaklyZero(jq - jm) and AbsolutePrecision(jq - jm) ge 2*bound;
  jN := &+[ Coefficient(jq, k)*q^(N*k) : k in [-1..p0-1] | Coefficient(jq, k) ne 0 ] + O(q^(N*p0));
  pj := [ Lq!1 ]; pN := [ Lq!1 ];
  for k in [1..N+1] do Append(~pj, pj[k]*jq); Append(~pN, pN[k]*jN); end for;
  cs := Coefficients(Phi); ms := Monomials(Phi);
  g := &+[ cs[i]*pj[Degree(ms[i], 1) + 1]*pN[Degree(ms[i], 2) + 1] : i in [1..#cs] ];
  printf "N = %o: j(q) known modulo q^%o; Phi_N(j(q), j(q^N)) = O(q^%o); (N+1)^2 = %o\n", N, p0, AbsolutePrecision(g), bound;
  assert IsWeaklyZero(g) and AbsolutePrecision(g) gt bound;
  assert AbsolutePrecision(g) eq (N eq 2 select 72 else 402);
  return Phi;
end function;
Phi2 := CheckModularPolynomial(2);
Phi17 := CheckModularPolynomial(17);
print "CHECK (C0): Phi_2 and Phi_17 are the classical modular polynomials";

// ---------------------------------------------------------------------------------------
// (C) X_0(34)
N := 34;
X := SmallModularCurve(N);
P3<x,y,z> := CoordinateRing(AmbientSpace(X));
fX := DefiningPolynomial(X);
IX := Ideal(X);
print "model:", fX;
assert fX eq x^4 - y^4 + x^3*z + 3*x*y^2*z - 2*x^2*z^2 + x*z^3 + z^4;
assert IsNonsingular(X) and Degree(X) eq 4 and Genus(X) eq 3;

// (C1) model check
S := CuspidalSubspace(ModularSymbols(N, 2, +1));
assert Dimension(S) eq 3;
for p in PrimesInInterval(3, 100) do
  if p eq 17 then continue; end if;
  Xp := ChangeRing(X, GF(p));
  assert IsNonsingular(Xp);
  assert #Points(Xp) eq p + 1 - Trace(HeckeOperator(S, p));
end for;
print "CHECK (C1): nonsingular plane quartic, #X(F_p) = p+1-tr T_p for 3 <= p < 100, p <> 17";

// (C2a) the j-function and the Atkin-Lehner involutions
FX := FunctionField(X);
j := FX!jFunction(X, N);
assert Degree(j) eq 54;
w2 := AtkinLehnerInvolution(X, N, 2);
w17 := AtkinLehnerInvolution(X, N, 17);
w2eqs := DefiningPolynomials(w2);
w17eqs := DefiningPolynomials(w17);
assert w2eqs eq [z, -y, x] and w17eqs eq [x, -y, z];
assert Evaluate(fX, w2eqs) eq fX and Evaluate(fX, w17eqs) eq fX;     // both preserve X
assert [ Evaluate(e, w2eqs) : e in w2eqs ] eq [x, y, z];             // involutions
assert [ Evaluate(e, w17eqs) : e in w17eqs ] eq [x, y, z];
assert w2*w2 eq IdentityMap(X) and w17*w17 eq IdentityMap(X);
j2 := Pullback(w2, j);
assert Evaluate(Phi2, [j, j2]) eq 0;
assert j2 ne j;
print "CHECK (C2a): deg j = 54, w2 = (z:-y:x), w17 = (x:-y:z), and Phi_2(j, j o w2) = 0 in the function field of X";

// (C2b) cross-check with Magma's jInvariant at points over number fields
Qy<Y> := PolynomialRing(Rationals());
nchk := 0;
for xv in [2..40] do
  g := Evaluate(fX, [xv, Y, 1]);
  fac := [ t[1] : t in Factorization(g) | Degree(t[1]) ge 2 ];
  if #fac eq 0 then continue; end if;
  L<y0> := NumberField(fac[1]);
  XL := BaseChange(X, L);
  R := XL![xv, y0, 1];
  j0 := jInvariant(R, N);
  a2 := AtkinLehnerInvolution(XL, N, 2);
  a17 := AtkinLehnerInvolution(XL, N, 17);
  jj2 := jInvariant(a2(R), N);
  jj17 := jInvariant(a17(R), N);
  assert Evaluate(Phi2, [j0, jj2]) eq 0 and Evaluate(Phi17, [j0, jj2]) ne 0;
  assert Evaluate(Phi17, [j0, jj17]) eq 0 and Evaluate(Phi2, [j0, jj17]) ne 0;
  nchk +:= 1;
  if nchk eq 5 then break; end if;
end for;
assert nchk eq 5;
print "CHECK (C2b): Phi_2 and Phi_17 tests with Magma's jInvariant at 5 points over number fields";

// (C2c) Phi_17(j, j o w17) = 0 on X.
// j = jnum/jden with forms jnum, jden of degree 53 in x, y, z.  The function
// Phi_17(j, j o w17) has polar divisor at most 18*(j)_infty + 18*(j o w17)_infty, of degree
// 18*54 + 18*54 = 1944.  It vanishes at every point (x0 : y0 : 1) of X tested below; these
// are pairwise distinct geometric points (x0 rational, y0 a root of a squarefree quartic),
// counted with the degree of y0.
jr := ProjectiveFunction(j);
jnum := P3!Numerator(jr);
jden := P3!Denominator(jr);
assert IsHomogeneous(jnum) and IsHomogeneous(jden) and Degree(jnum) eq Degree(jden);
assert not (jden in IX);
assert FX!(jnum/jden) eq j;
// a second computation of deg j = 54, independent of Degree(j): on the nonsingular curve X
// the polar divisor of jnum/jden has degree 4*deg(jden) minus the length of the base scheme
// X cap {jnum = jden = 0}
Zb := Scheme(AmbientSpace(X), [fX, jnum, jden]);
assert Dimension(Zb) eq 0;
printf "deg jden = %o, length of the base scheme = %o\n", Degree(jden), Degree(Zb);
assert 4*Degree(jden) - Degree(Zb) eq 54;
npts := 0; nx := 0; xv := 1;
while npts le 2000 do
  xv +:= 1;
  h := Evaluate(fX, [xv, Y, 1]);
  assert Degree(h) eq 4;
  if not IsSquarefree(h) then continue; end if;
  for t in Factorization(h) do
    d := Degree(t[1]);
    if d eq 1 then
      L := Rationals(); y0 := -Coefficient(t[1], 0)/Coefficient(t[1], 1);
    else
      L := NumberField(t[1]); y0 := L.1;
    end if;
    R := [ L!xv, y0, L!1 ];
    assert Evaluate(fX, R) eq 0;
    Sw := [ Evaluate(e, R) : e in w17eqs ];
    assert Evaluate(fX, Sw) eq 0 and Sw[3] eq 1;
    dR := Evaluate(jden, R); dS := Evaluate(jden, Sw);
    if dR eq 0 or dS eq 0 then continue; end if;      // not used: j or j o w17 may have a pole
    assert Evaluate(Phi17, [ Evaluate(jnum, R)/dR, Evaluate(jnum, Sw)/dS ]) eq 0;
    npts +:= d;
  end for;
  nx +:= 1;
end while;
printf "Phi_17(j, j o w17) vanishes at %o geometric points of X (x0 = 2..%o, %o values used)\n", npts, xv, nx;
assert npts gt 1944;
assert npts eq 2004 and xv eq 502 and nx eq 501;
print "CHECK (C2c): Phi_17(j, j o w17) = 0 in the function field of X";

// (C3) the quotient
G := AutomorphismGroup(X, [w2]);
assert #G eq 2;
Yq, pi := CurveQuotient(G);
print "quotient curve:", Yq;
assert IsNonsingular(Yq) and Genus(Yq) eq 1;
pieqs := DefiningPolynomials(pi);
assert #pieqs eq 4 and not &and[ e in IX : e in pieqs ];
assert &and[ Evaluate(e, pieqs) in IX : e in DefiningPolynomials(Yq) ];    // pi maps X to Yq
pw := [ Evaluate(e, w2eqs) : e in pieqs ];
assert &and[ (pw[i]*pieqs[k] - pw[k]*pieqs[i]) in IX : i, k in [1..4] ];    // pi o w2 = pi
FY := FunctionField(Yq);
AY := AmbientSpace(Yq);
pb := [ Pullback(pi, FY!(AY.i/AY.4)) : i in [1..3] ];
assert &and[ Pullback(w2, t) eq t : t in pb ];                              // pi o w2 = pi, again
assert Degree(pi) eq 2;
print "CHECK (C3a): Yq is a nonsingular curve of genus 1 and pi : X -> Yq has degree 2 with pi o w2 = pi";
allY := PointSearch(Yq, 1000);
E, phi := EllipticCurve(Yq, allY[1]);
Em, psi := MinimalModel(E);
print "X_0(34)/w_2 is isomorphic to", aInvariants(Em), CremonaReference(Em);
assert aInvariants(Em) eq [1, -1, 1, -6, -4];
assert Conductor(Em) eq 17;
assert CremonaReference(Em) eq "17a2";
// phi : Yq -> E has degree 1 and both curves are nonsingular, so phi is an isomorphism;
// psi : E -> Em is the isomorphism onto the minimal model.
assert Degree(phi) eq 1;
assert IsIsomorphism(psi);
Tor, mT := TorsionSubgroup(Em);
assert Invariants(Tor) eq [2, 2];
// Rank 0.  The rational 2-torsion of Em is full, so Em(Q)/2Em(Q) has order 2^(rank+2), and
// it embeds into the 2-Selmer group, which has order 4.  (For a curve with full rational
// 2-torsion the 2-descent involves no class group of a number field.)
assert #TwoTorsionSubgroup(Em) eq 4;
assert #TwoSelmerGroup(Em) eq 4;
r, proved := Rank(Em);
assert r eq 0 and proved;
print "[1,-1,1,-6,-4]: conductor 17, torsion Z/2 x Z/2, 2-Selmer group of order 4, rank 0";
// phi*psi is an isomorphism from Yq onto the minimal model, whose group of rational points
// is its torsion subgroup (rank 0); so it suffices to find 4 rational points of Yq with
// distinct images.
comp := phi * psi;
imgs := [ comp(q) : q in allY ];
print "rational points of the quotient:", allY, " images:", imgs;
assert #allY eq 4 and #Set(imgs) eq #Tor;
print "CHECK (C3b): the quotient has exactly 4 rational points";

// (C4) the fibres
cuspTable := [ <[Rationals()| 0, 1, 1], 34>, <[Rationals()| -1, 1, 0], 17>, <[Rationals()| 0, -1, 1], 2>, <[Rationals()| 1, 1, 0], 1> ];
mps := {};
nCusp := 0; nQuad := 0; poleDegree := 0; poleOrders := [];
for Pt in allY do
  D := Pullback(pi, Place(Pt));
  assert Degree(D) eq 2;
  supp := Support(D);
  printf "fibre over the rational point %o of X_0(34)/w_2:\n", Pt;
  if #supp eq 2 then
    for pl in supp do
      assert Degree(pl) eq 1 and Valuation(D, pl) eq 1;
      assert Valuation(j, pl) lt 0;            // a pole of j: a cusp
      Pc := RepresentativePoint(pl);
      assert pi(Pc) eq Pt;
      assert w2(Pc) ne Pc and pi(w2(Pc)) eq Pt;
      printf "   rational point %o : cusp (pole of j of order %o), w2 maps it to %o\n", Pc, -Valuation(j, pl), w2(Pc);
      nCusp +:= 1;
      poleDegree +:= -Valuation(j, pl);
      Append(~poleOrders, -Valuation(j, pl));
      assert exists{ t : t in cuspTable | t[1] eq Eltseq(Pc) and t[2] eq -Valuation(j, pl) };
    end for;
  else
    assert #supp eq 1;
    pl := supp[1];
    assert Degree(pl) eq 2 and Valuation(D, pl) eq 1;
    K := ResidueClassField(pl);
    KK := AbsoluteField(NumberField(K));
    dK := Discriminant(MaximalOrder(KK));
    assert Valuation(j, pl) ge 0;
    jv := Evaluate(j, pl);
    mp := MinimalPolynomial(jv);
    printf "   one place of degree 2, residue field of discriminant %o, j has minimal polynomial %o\n", dK, mp;
    assert dK eq -15;
    assert Degree(mp) eq 2;
    Include(~mps, [ Rationals()!c : c in Coefficients(mp) ]);
    assert not &and[ IsIntegral(c) : c in Coefficients(mp) ];   // j is not an algebraic integer: no CM
    // the two conjugate values of j are 2-isogenous, as they must be (w2 exchanges the two points)
    rts := [ t[1] : t in Roots(mp, KK) ];
    assert #rts eq 2;
    assert Evaluate(Phi2, [rts[1], rts[2]]) eq 0;
    nQuad +:= 1;
  end if;
end for;
assert nCusp eq 4 and nQuad eq 2;
assert poleDegree eq 54;      // the four rational points are all the poles of j
assert Sort(poleOrders) eq [1, 2, 17, 34];
assert mps eq { [Rationals()| 25153757/8, -11779/4, 1], [Rationals()| 81029690385834182553917/2251799813685248, 7319387769191/17179869184, 1] };
print "CHECK (C4): the fibres are 2+2 rational cusps (all the poles of j, of orders 34, 17, 2, 1) and two pairs of conjugate points over Q(sqrt -15)";

// (C4') the same quadratic points, directly: over Q(s), s^2 = -15, the four points
// (x0 : y0 : 1) below lie on X, are not rational, are mapped by pi to the two rational points
// of the quotient that are not images of cusps, and w2 maps each of them to its conjugate.
K15<s> := QuadraticField(-15);
sig := hom< K15 -> K15 | -s >;
xf := FX!(x/z); yf := FX!(y/z);
nDirect := 0;
quadPts := {};
for Pt in allY do
  D := Pullback(pi, Place(Pt));
  supp := Support(D);
  if #supp eq 2 then continue; end if;
  pl := supp[1];
  assert Valuation(xf, pl) ge 0 and Valuation(yf, pl) ge 0;
  KK := AbsoluteField(NumberField(ResidueClassField(pl)));
  ok, iso := IsIsomorphic(KK, K15);
  assert ok;
  x0 := iso(KK!Evaluate(xf, pl)); y0 := iso(KK!Evaluate(yf, pl));
  for p in [ [x0, y0, K15!1], [sig(x0), sig(y0), K15!1] ] do
    assert Evaluate(fX, p) eq 0;
    Include(~quadPts, [p[1], p[2]]);
    assert not (p[1] in Rationals() and p[2] in Rationals());
    v := [ Evaluate(e, p) : e in pieqs ];
    c := Eltseq(Pt);
    assert v ne [0, 0, 0, 0] and &and[ v[i]*c[k] - v[k]*c[i] eq 0 : i, k in [1..4] ];   // pi(p) = Pt
    wp := [ p[3], -p[2], p[1] ];                  // w2(p)
    cp := [ sig(a) : a in p ];                    // the conjugate of p
    assert &and[ wp[i]*cp[k] - wp[k]*cp[i] eq 0 : i, k in [1..3] ];
    printf "   point (%o : %o : 1) over Q(s), s^2 = -15, above %o: w2 maps it to its conjugate\n", p[1], p[2], Pt;
    nDirect +:= 1;
  end for;
end for;
assert nDirect eq 4;
assert quadPts eq { [(7+s)/8, (1-s)/8], [(7-s)/8, (1+s)/8], [(7+s)/8, (-1+s)/8], [(7-s)/8, (-1-s)/8] };
print "CHECK (C4'): the four points over Q(sqrt -15), with w2(P) equal to the conjugate of P";

print "ALL_CHECKS_PASSED";
quit;
