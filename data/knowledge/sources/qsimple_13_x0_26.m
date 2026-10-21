// qsimple_13_x0_26.m
// Certificate for the case Gamma = C_4 of proofs/qsimple-13-obstruction.html (Section 7):
// the rational points of X_0(26)/w_2 and the points of X_0(26) above them.
//
// Checks:
//  (1) Magma's SmallModularCurve(26) is a model of X_0(26): #X(F_p) = p + 1 - tr(T_p | S_2(26))
//      for all primes 3 <= p < 100, p <> 13.
//  (2) AtkinLehnerInvolution(X,26,2) is the moduli involution w_2: at several quadratic points R,
//      j(R) and j(w_2 R) are 2-isogenous and not 13-isogenous (and the reverse for w_13).
//  (3) X_0(26)/w_2 is isomorphic to the elliptic curve 26a1 = [1,0,1,-5,-8], which has rank 0
//      (proved) and torsion Z/3; so X_0(26)/w_2 has exactly 3 rational points.
//  (4) Above these 3 points lie the 4 rational cusps of X_0(26) and one pair of conjugate points
//      over Q(sqrt 13) whose j-invariant is a root of the Hilbert class polynomial H_{-52}
//      (CM points). Every other point of X_0(26) above a rational point of X_0(26)/w_2 is
//      therefore a cusp or a CM point.
//  (5) Consequently X_0(26)(Q) consists of the 4 rational cusps (no non-CM elliptic curve over Q
//      has a rational cyclic 26-isogeny).
SetColumns(0);
SetMemoryLimit(16*10^9);
v1, v2, v3 := GetVersion();
printf "Magma V%o.%o-%o\n", v1, v2, v3;

N := 26;
X := SmallModularCurve(N);
print "model:", X;
assert Genus(X) eq 2;

// (1) model check
S := CuspidalSubspace(ModularSymbols(N, 2, +1));
assert Dimension(S) eq 2;
for p in PrimesInInterval(3, 100) do
  if p eq 13 then continue; end if;
  assert #ChangeRing(X, GF(p)) eq p + 1 - Trace(HeckeOperator(S, p));
end for;
print "CHECK model: #X(F_p) = p+1-tr T_p for 3 <= p < 100, p <> 13";

// (2) moduli check of w_2 (and w_13)
w2 := AtkinLehnerInvolution(X, N, 2);
w13 := AtkinLehnerInvolution(X, N, 13);
Phi2 := ClassicalModularPolynomial(2);
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
  a2 := AtkinLehnerInvolution(XL, N, 2);
  a13 := AtkinLehnerInvolution(XL, N, 13);
  j2 := jInvariant(a2(R), N);
  j13 := jInvariant(a13(R), N);
  assert Evaluate(Phi2, [j0, j2]) eq 0 and Evaluate(Phi13, [j0, j2]) ne 0;
  assert Evaluate(Phi13, [j0, j13]) eq 0 and Evaluate(Phi2, [j0, j13]) ne 0;
  nchk +:= 1;
  if nchk eq 5 then break; end if;
end for;
assert nchk eq 5;
print "CHECK w_2 is the moduli Atkin-Lehner involution (Phi_2 test at 5 quadratic points)";

// (3) the quotient
G := AutomorphismGroup(X, [w2]);
Yq, pi := CurveQuotient(G);
assert Genus(Yq) eq 1;
pts := PointSearch(Yq, 100);
E, phi := EllipticCurve(Yq, pts[1]);
Em, psi := MinimalModel(E);
print "X_0(26)/w_2 is isomorphic to", aInvariants(Em), CremonaReference(Em);
assert aInvariants(Em) eq [1, 0, 1, -5, -8];
r, proved := Rank(Em);
assert r eq 0 and proved;
T, mT := TorsionSubgroup(Em);
assert #T eq 3;
print "26a1: rank 0 (proved), torsion Z/3";
// all rational points of the quotient: preimages of the 3 torsion points under phi*psi
// all rational points of the quotient: phi*psi is an isomorphism onto the minimal model, whose
// group of rational points is its torsion subgroup (rank 0); so it suffices to find 3 rational
// points of the quotient with distinct images.
comp := phi * psi;
allY := PointSearch(Yq, 100);
imgs := [ comp(q) : q in allY ];
print "rational points of the quotient:", allY, " images:", imgs;
assert #allY eq 3 and #Set(imgs) eq #T;

// (4) the fibres
H52 := HilbertClassPolynomial(-52);
print "H_{-52} =", H52;
nRat := 0; nCM := 0;
for Pt in allY do
  F := Pt @@ pi;
  printf "fibre over a rational point of X_0(26)/w_2 (reduced degree %o):\n", Degree(F);
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
      mp := MinimalPolynomial(jR);
      printf "   %o over a field of degree %o: j has minimal polynomial %o\n", R, Degree(K), mp;
      assert Degree(K) eq 2 and Discriminant(MaximalOrder(AbsoluteField(K))) eq 13;
      assert Parent(mp)!H52 eq mp;     // CM by the order of discriminant -52
      nCM +:= 1;
    end if;
  end for;
end for;
assert nRat eq 4 and nCM eq 2;
print "CHECK fibres: 4 rational cusps and one conjugate pair over Q(sqrt 13) with j a root of H_{-52}";

// (5) consequence
print "CHECK X_0(26)(Q) consists of the 4 rational cusps";

print "ALL_CHECKS_PASSED";
quit;
