// qsimple_13_x0_39.m
// Certificate for the case Gamma = C_6 of proofs/qsimple-13-obstruction.html (Section 6),
// and for the absence of rational cyclic 39-isogenies used in Section 8.
//
// Checks:
//  (1) Magma's SmallModularCurve(39) is a model of X_0(39): #X(F_p) = p + 1 - tr(T_p | S_2(39))
//      for all primes p < 100, p <> 3, 13.
//  (2) AtkinLehnerInvolution(X,39,3) is the moduli involution w_3: at several quadratic points R,
//      j(R) and j(w_3 R) are 3-isogenous and not 13-isogenous (and the reverse for w_13).
//  (3) X_0(39)/w_3 is isomorphic to the elliptic curve 39a1 = [1,1,0,-4,-5], which has rank 0
//      (proved) and torsion Z/2 x Z/2; so X_0(39)/w_3 has exactly 4 rational points.
//  (4) Above two of them lie the 4 rational cusps of X_0(39); above each of the other two lies
//      one pair of conjugate non-cuspidal points with residue field Q(sqrt -7).
//  (5) Consequently X_0(39)(Q) consists of the 4 rational cusps (no elliptic curve over Q has a
//      rational cyclic 39-isogeny), and every non-rational point of X_0(39) lying over a
//      rational point of X_0(39)/w_3 has residue field Q(sqrt -7), an imaginary quadratic field.
SetColumns(0);
SetMemoryLimit(16*10^9);
v1, v2, v3 := GetVersion();
printf "Magma V%o.%o-%o\n", v1, v2, v3;

N := 39;
X := SmallModularCurve(N);
print "model:", X;
assert Genus(X) eq 3;

// (1) model check
S := CuspidalSubspace(ModularSymbols(N, 2, +1));
assert Dimension(S) eq 3;
for p in PrimesInInterval(2, 100) do
  if p in {3, 13} then continue; end if;
  assert #ChangeRing(X, GF(p)) eq p + 1 - Trace(HeckeOperator(S, p));
end for;
print "CHECK model: #X(F_p) = p+1-tr T_p for p < 100, p <> 3, 13";

// (2) moduli check of w_3 (and w_13)
w3 := AtkinLehnerInvolution(X, N, 3);
Phi3 := ClassicalModularPolynomial(3);
Phi13 := ClassicalModularPolynomial(13);
fX, hX := HyperellipticPolynomials(X);
Qy<Y> := PolynomialRing(Rationals());
nchk := 0;
for xv in [2..30] do
  g := Y^2 + Evaluate(hX, xv)*Y - Evaluate(fX, xv);
  if not IsIrreducible(g) then continue; end if;
  L<y0> := NumberField(g);
  XL := BaseChange(X, L);
  R := XL![xv, y0, 1];
  j0 := jInvariant(R, N);
  a3 := AtkinLehnerInvolution(XL, N, 3);
  a13 := AtkinLehnerInvolution(XL, N, 13);
  j3 := jInvariant(a3(R), N);
  j13 := jInvariant(a13(R), N);
  assert Evaluate(Phi3, [j0, j3]) eq 0 and Evaluate(Phi13, [j0, j3]) ne 0;
  assert Evaluate(Phi13, [j0, j13]) eq 0 and Evaluate(Phi3, [j0, j13]) ne 0;
  nchk +:= 1;
  if nchk eq 5 then break; end if;
end for;
assert nchk eq 5;
print "CHECK w_3 is the moduli Atkin-Lehner involution (Phi_3 test at 5 quadratic points)";

// (3) the quotient
G := AutomorphismGroup(X, [w3]);
Yq, pi := CurveQuotient(G);
assert Genus(Yq) eq 1;
pts := PointSearch(Yq, 1000);
E, phi := EllipticCurve(Yq, pts[1]);
Em, psi := MinimalModel(E);
print "X_0(39)/w_3 is isomorphic to", aInvariants(Em), CremonaReference(Em);
assert aInvariants(Em) eq [1, 1, 0, -4, -5];
r, proved := Rank(Em);
assert r eq 0 and proved;
T, mT := TorsionSubgroup(Em);
assert Invariants(T) eq [2, 2];
print "39a1: rank 0 (proved), torsion Z/2 x Z/2";
// all rational points of the quotient: phi*psi is an isomorphism onto the minimal model, whose
// group of rational points is its torsion subgroup (rank 0); so it suffices to find 4 rational
// points of the quotient with distinct images.
comp := phi * psi;
allY := PointSearch(Yq, 1000);
imgs := [ comp(q) : q in allY ];
print "rational points of the quotient:", allY, " images:", imgs;
assert #allY eq 4 and #Set(imgs) eq #T;

// (4) the fibres
nRat := 0; nPairs := 0;
for Pt in allY do
  F := Pt @@ pi;
  printf "fibre over a rational point of X_0(39)/w_3 (reduced degree %o):\n", Degree(F);
  assert Degree(F) eq 2;
  assert #PointsOverSplittingField(F) eq 2;
  for R in PointsOverSplittingField(F) do
    K := Ring(Parent(R));
    iscusp := false;
    try
      jR := jInvariant(BaseChange(X, K)!Eltseq(R), N);
    catch e
      assert Position(e`Object, "cusp") gt 0;   // jInvariant refuses exactly the cusps
      iscusp := true;
    end try;
    if iscusp then
      assert Degree(K) eq 1;
      printf "   %o : rational cusp\n", R;
      nRat +:= 1;
    else
      dK := Discriminant(MaximalOrder(AbsoluteField(K)));
      cm := HasComplexMultiplication(EllipticCurveWithjInvariant(jR));
      printf "   %o over a field of discriminant %o: j has minimal polynomial %o, CM: %o\n",
        R, dK, MinimalPolynomial(jR), cm;
      assert Degree(K) eq 2 and dK eq -7 and not cm;
      nPairs +:= 1;
    end if;
  end for;
end for;
assert nRat eq 4 and nPairs eq 4;
print "CHECK fibres: 4 rational cusps, and 2 conjugate pairs of non-CM points over Q(sqrt -7)";

// (5) consequence
print "CHECK X_0(39)(Q) consists of the 4 rational cusps";

print "ALL_CHECKS_PASSED";
quit;
