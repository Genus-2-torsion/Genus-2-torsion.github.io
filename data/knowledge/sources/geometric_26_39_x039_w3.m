// geometric_26_39_x039_w3.m
//
// Computation for the census proof "No rational points of order 26 or 39 on
// geometrically split abelian surfaces" (21 October 2026), Section 9.
//
// It checks the facts about the modular curve X_0(39) and its quotient by the
// Atkin-Lehner involution w_3 that the proof uses:
//   (1) the model of X_0(39) from Magma's small modular curve database passes
//       basic consistency checks (genus, degree and cusp widths of the j-map,
//       point counts against S_2(Gamma_0(39)));
//   (2) Magma's involution AtkinLehnerInvolution(X,39,3) is the moduli
//       involution w_3: Phi_3(j, j o w_3) = 0 identically, while
//       Phi_13(j, j o w_3), Phi_3(j, j o w_13), Phi_3(j, j o w_39) are nonzero;
//   (3) X_0(39)/w_3 is a genus-1 curve isomorphic over Q to the elliptic curve
//       39a1, which has rank 0 (unconditionally) and torsion Z/2 x Z/2, so the
//       quotient has exactly 4 rational points;
//   (4) the fibres of X_0(39) -> X_0(39)/w_3 over these 4 points: over two of
//       them the fibre consists of two rational cusps; over each of the other
//       two the fibre is a pair of conjugate points with residue field
//       Q(sqrt(-7)).  In particular X_0(39)(Q) consists of the 4 cusps, and a
//       non-rational point of X_0(39) over a quadratic field K whose image in
//       X_0(39)/w_3 is rational forces K = Q(sqrt(-7)).
//
// Every claimed fact is checked by an assert; the script ends with the marker
// ALL_CHECKS_PASSED.  Run from this directory:
//     magma -b geometric_26_39_x039_w3.m < /dev/null > geometric_26_39_x039_w3.log 2>&1

SetColumns(0);
SetMemoryLimit(16*10^9);
v1, v2, v3 := GetVersion();
printf "Magma V%o.%o-%o\n", v1, v2, v3;

// ---------------------------------------------------------------------------
// (1) The model of X_0(39).
// ---------------------------------------------------------------------------
X := SmallModularCurve(39);
print "MODEL", X;
assert Genus(X) eq 3;

j := jFunction(X, 39);
assert Degree(j) eq 56;                         // [PSL_2(Z) : Gamma_0(39)] = 56
cusps := [Cusp(X, 39, d) : d in [1, 3, 13, 39]];
assert #Set(cusps) eq 4;
poles := [-Valuation(j, Place(c)) : c in cusps];
print "pole orders of j at the four cusps:", poles;
assert Sort(poles) eq [1, 3, 13, 39];           // the cusp widths of Gamma_0(39)

M := CuspidalSubspace(ModularSymbols(39, 2, +1));
assert Dimension(M) eq 3;
NF := NewformDecomposition(M);
for p in [p : p in PrimesUpTo(100) | p notin {3, 13}] do
  np := #Points(ChangeRing(X, GF(p)));
  tr := &+[Integers()!Trace(HeckeOperator(f, p)) : f in NF];
  assert np eq p + 1 - tr;
end for;
print "CHECK #X(F_p) = p + 1 - Tr(T_p | S_2(Gamma_0(39))) for all primes p < 100, p != 3, 13: ok";

// ---------------------------------------------------------------------------
// (2) The Atkin-Lehner involutions and the moduli interpretation of w_3.
// ---------------------------------------------------------------------------
w3  := AtkinLehnerInvolution(X, 39, 3);
w13 := AtkinLehnerInvolution(X, 39, 13);
w39 := AtkinLehnerInvolution(X, 39, 39);
idX := IdentityMap(X);
assert w3*w3 eq idX and w3 ne idX;
assert w13*w13 eq idX and w13 ne idX;
assert w3*w13 eq w39;

j3  := Pullback(w3, j);
j13 := Pullback(w13, j);
j39 := Pullback(w39, j);
Phi3  := ClassicalModularPolynomial(3);
Phi13 := ClassicalModularPolynomial(13);
assert Evaluate(Phi3,  [j, j3])  eq 0;          // j o w_3 is 3-isogenous to j
assert Evaluate(Phi13, [j, j3])  ne 0;
assert Evaluate(Phi13, [j, j13]) eq 0;          // j o w_13 is 13-isogenous to j
assert Evaluate(Phi3,  [j, j13]) ne 0;
assert Evaluate(Phi3,  [j, j39]) ne 0;
assert Evaluate(Phi3,  [j, j])   ne 0;
print "CHECK Phi_3(j, j o w3) = 0 and Phi_13(j, j o w13) = 0 in the function field;";
print "      Phi_13(j, j o w3), Phi_3(j, j o w13), Phi_3(j, j o w39), Phi_3(j, j) are nonzero: ok";

Mw := Kernel(AtkinLehner(M, 3) - 1);
assert Dimension(Mw) eq 1;                      // so X_0(39)/w_3 has genus 1
print "dim S_2(Gamma_0(39))^{w_3 = +1} =", Dimension(Mw);

// ---------------------------------------------------------------------------
// (3) The quotient X_0(39)/w_3 and its rational points.
// ---------------------------------------------------------------------------
G := AutomorphismGroup(X, [w3]);
Q, pi := CurveQuotient(G);
print "QUOTIENT", Q;
assert Genus(Q) eq 1;
assert Degree(pi) eq 2;

pts := PointSearch(Q, 1000);
print "rational points found on X_0(39)/w_3:", pts;
assert #pts eq 4;

E, phi := EllipticCurve(Q, pts[1]);             // phi : Q -> E
ok, phiinv := IsInvertible(phi);
assert ok;                                      // phi is an isomorphism over Q
Em, psi := MinimalModel(E);
print "X_0(39)/w_3 is isomorphic over Q to", Em, "Cremona label", CremonaReference(Em);
assert aInvariants(Em) eq [1, 1, 0, -4, -5];    // 39a1
assert CremonaReference(Em) eq "39a1";

T, mT := TorsionSubgroup(Em);
assert Invariants(T) eq [2, 2];
// Full rational 2-torsion: the 2-division polynomial splits over Q, so the
// 2-Selmer group is computed from the etale algebra Q x Q x Q and no class
// group of a number field is involved.  #Sel^2 = 4 = #E(Q)[2] gives rank 0.
assert #Roots(DivisionPolynomial(Em, 2)) eq 3;
S2 := TwoSelmerGroup(Em);
assert #S2 eq 4;
r, proved := Rank(Em);
assert r eq 0 and proved;
print "39a1: torsion Z/2 x Z/2, #Sel^2 = 4, rank 0 (unconditional)";

imgs := [psi(phi(q)) : q in pts];
print "images of the 4 points in 39a1:", imgs;
assert #Set(imgs) eq 4;                         // all of 39a1(Q), so Q(Q) = pts
print "CHECK X_0(39)/w_3 has exactly 4 rational points: ok";

// ---------------------------------------------------------------------------
// (4) The fibres over the 4 rational points.
// ---------------------------------------------------------------------------
K<s> := QuadraticField(-7);
XK := BaseChange(X, K);
w3K := AtkinLehnerInvolution(XK, 39, 3);
sig := hom<K -> K | -s>;
cuspfibres := 0;
sqrtm7fibres := 0;
ratpts := {@ @};
for q in pts do
  F := q @@ pi;
  assert Dimension(F) eq 0;
  assert Degree(ReducedSubscheme(F)) eq 2;      // exactly two geometric points
  RF := RationalPoints(F);
  if #RF eq 2 then
    cuspfibres +:= 1;
    for R in RF do
      RX := X!Eltseq(R);
      assert RX in cusps;
      Include(~ratpts, RX);
    end for;
    printf "fibre over %o: two rational points, both cusps: %o\n", q, RF;
  else
    assert #RF eq 0;
    RK := RationalPoints(F, K);
    assert #RK eq 2;                            // a conjugate pair over Q(sqrt(-7))
    sqrtm7fibres +:= 1;
    printf "fibre over %o: one conjugate pair over Q(sqrt(-7)): %o\n", q, RK;
    for R in RK do
      RX := XK!Eltseq(R);
      assert RX notin [XK!Eltseq(c) : c in cusps];
      Rc := XK![sig(a) : a in Eltseq(RX)];
      assert Rc ne RX;                          // not a rational point
      assert w3K(RX) eq Rc;                     // w_3(R) = conjugate of R
      jR := jInvariant(RX, 39);
      jw := jInvariant(w3K(RX), 39);
      assert jw eq sig(jR);
      assert Evaluate(Phi3, [jR, jw]) eq 0;
      assert not HasComplexMultiplication(EllipticCurveWithjInvariant(jR));
      printf "   point %o: j = %o (not CM), w_3 maps it to its conjugate\n", RX, jR;
    end for;
  end if;
end for;
assert cuspfibres eq 2 and sqrtm7fibres eq 2;
assert #ratpts eq 4 and Set(ratpts) eq Set(cusps);
print "CHECK X_0(39)(Q) = the 4 cusps; the other two fibres are conjugate pairs over Q(sqrt(-7)): ok";
print "CONCLUSION: if K is a quadratic field and x in X_0(39)(K) is not rational with pi(x) rational, then K = Q(sqrt(-7)).";

print "ALL_CHECKS_PASSED";
quit;
