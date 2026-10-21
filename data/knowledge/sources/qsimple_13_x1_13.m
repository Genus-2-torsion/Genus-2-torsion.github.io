// qsimple_13_x1_13.m
// Certificate for the case "l = 2, image C_3, discriminant in d_K Q*^2" of
// proofs/qsimple-13-obstruction.html (Section 9.5(c1)).
//
// Checks:
//  (1) C : y^2 + (x^3+x^2+1) y = x^2 + x has genus 2, and the Tate normal form
//      E(b,c) : Y^2 + (1-c) X Y - b Y = X^3 - b X^2 with r = 1 - x y, s = 1 - x y/(y+1),
//      b = r s (r-1), c = s (r-1), has (0,0) of exact order 13 over the function field of C.
//      Hence C is a model of X_1(13) (a nonconstant map between curves of genus 2 has degree 1).
//  (2) j(E(b,c)) lies in Q(x), so it is invariant under the hyperelliptic involution
//      y -> -y - (x^3+x^2+1); its denominator is x^13 (x+1)^13 (x^3+4x^2+x-1).
//      A point of C with x = t in Q is defined over K = Q(sqrt D(t)),
//      D(t) = (t^3+t^2+1)^2 + 4(t^2+t), and D has no real root.
//  (3) j - 1728 = -D(t) g(t)^2 / den(j) with g in Q[t] of degree 18. Hence
//      (j(t) - 1728) D(t) is, modulo squares, -t(t+1)(t^3+4t^2+t-1), unless g(t) = 0 (then j = 1728).
//  (4) The genus-2 curve C1 : s^2 = -t(t+1)(t^3+4t^2+t-1) has Mordell-Weil rank 0,
//      unconditionally: the class group of the cubic field Q[t]/(t^3+4t^2+t-1) (discriminant 169)
//      is trivial with Proof := "Full", and #Sel_2(J) = #J(Q)[2] = 4. Chabauty0 gives
//      C1(Q) = {(-1,0), (0,0), infinity}, i.e. t in {-1, 0, infinity}: cusps (poles of j).
//  (5) Sanity check (not used in the proof): the point with t = -2 lies over K = Q(sqrt 17);
//      the curve E there has a K-point of order 13, j(E) in Q, E is K-isomorphic to its
//      conjugate, and its Q-model of minimal twist has a rational 13-isogeny.
// Magma's class group computations are rigorous unless SetClassGroupBounds("GRH") is called;
// this script never calls it.
SetColumns(0);
SetMemoryLimit(16*10^9);
v1, v2, v3 := GetVersion();
printf "Magma V%o.%o-%o\n", v1, v2, v3;

Q := Rationals();
P<t> := PolynomialRing(Q);

// (1) the model and the universal curve
Fx<x> := FunctionField(Q);
Py<Y> := PolynomialRing(Fx);
F<y> := FunctionField(Y^2 + (x^3+x^2+1)*Y - x^2 - x);
assert Genus(F) eq 2;
r := 1 - x*y; s := 1 - x*y/(y+1);
b := r*s*(r-1); c := s*(r-1);
E := EllipticCurve([1-c, -b, -b, 0, 0]);
P0 := E![0, 0];
assert IsZero(13*P0) and &and[ not IsZero(k*P0) : k in [1..12] ];
print "CHECK (0,0) has exact order 13 on the universal curve over the function field of C";

// (2) the j-invariant
ok, jx := IsCoercible(Fx, jInvariant(E));
assert ok;
Nj := Evaluate(Numerator(jx), t);
Dj := Evaluate(Denominator(jx), t);
c0 := LeadingCoefficient(Dj); Nj := Nj/c0; Dj := Dj/c0;
printf "j in Q(x): numerator degree %o, denominator %o\n", Degree(Nj), Factorization(Dj);
assert Degree(Nj) eq 42;
assert Dj eq t^13*(t+1)^13*(t^3+4*t^2+t-1);
Dk := (t^3+t^2+1)^2 + 4*(t^2+t);
assert Dk eq t^6 + 2*t^5 + t^4 + 2*t^3 + 6*t^2 + 4*t + 1;
assert #Roots(Dk, RealField(50)) eq 0;
print "CHECK D(t) = (t^3+t^2+1)^2 + 4(t^2+t) has no real root";

// (3) j - 1728
q := Nj - 1728*Dj;
quo, rem := Quotrem(q, Dk);
assert rem eq 0;
ok, g := IsSquare(-quo);
assert ok and Degree(g) eq 18;
print "CHECK j - 1728 = -D(t) g(t)^2 / den(j), deg g = 18";
f1 := -t*(t+1)*(t^3+4*t^2+t-1);
// (j - 1728) D = -D^2 g^2 / Dj, and Dj / (t(t+1)(t^3+4t^2+t-1)) = (t(t+1))^12 is a square
assert Dj eq (t*(t+1))^12 * (-f1);

// (4) the genus-2 curve C1
C1 := HyperellipticCurve(f1);
J := Jacobian(C1);
assert Genus(C1) eq 2;
Fc := NumberField(t^3+4*t^2+t-1);
Oc := MaximalOrder(Fc);
assert Discriminant(Oc) eq 169;
Cl := ClassGroup(Oc : Proof := "Full");
assert #Cl eq 1;
print "CHECK class group of Q[t]/(t^3+4t^2+t-1) (disc 169) is trivial (Proof := Full)";
Sel := TwoSelmerGroup(J);
J2 := TwoTorsionSubgroup(J);
printf "#Sel_2(J) = %o, #J(Q)[2] = %o\n", #Sel, #J2;
assert #Sel eq 4 and #J2 eq 4;
print "CHECK rank J(C1)(Q) = 0 (2-Selmer group equals the rational 2-torsion)";
Tor := TorsionSubgroup(J);
assert Invariants(Tor) eq [2, 2];
pts := Chabauty0(J);
print "C1(Q) =", pts;
assert Set(pts) eq { C1![1,0,0], C1![-1,0,1], C1![0,0,1] };
assert Evaluate(Dj, 0) eq 0 and Evaluate(Dj, -1) eq 0;
print "CHECK C1(Q) consists of t = -1, 0, infinity: cusps of X_1(13) (poles of j)";

// (5) sanity check at t = -2 (not used in the proof)
t0 := -2;
d0 := Evaluate(Dk, t0);
assert d0 eq 17;
K<w> := QuadraticField(17);
PK<Z> := PolynomialRing(K);
y0 := Roots(Z^2 + (t0^3+t0^2+1)*Z - t0^2 - t0)[1][1];
r0 := 1 - t0*y0; s0 := 1 - t0*y0/(y0+1);
bb := r0*s0*(r0-1); cc := s0*(r0-1);
EK := EllipticCurve([1-cc, -bb, -bb, 0, 0]);
assert Order(EK![0,0]) eq 13;
ok, j0 := IsCoercible(Q, jInvariant(EK));
assert ok and j0 eq Evaluate(jx, t0);
sig := hom< K -> K | -w >;
EKc := EllipticCurve([ sig(a) : a in aInvariants(EK) ]);
assert IsIsomorphic(EK, EKc);
E0 := MinimalTwist(EllipticCurveWithjInvariant(j0));
_, dcl := IsogenousCurves(E0);
printf "t = -2: K = Q(sqrt 17), j = %o, E0 = %o, isogeny class degree %o\n", j0, aInvariants(E0), dcl;
assert aInvariants(E0) eq [1, -1, 0, -139, 965];
assert dcl eq 13;   // E0 has a rational 13-isogeny
assert not HasComplexMultiplication(E0);

print "ALL_CHECKS_PASSED";
quit;
