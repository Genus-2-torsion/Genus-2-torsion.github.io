// qsimple_45_checks.m -- certificate for the census proof
// "No rational point of order 45 on a geometrically split abelian surface simple over Q"
// (14 November 2026).  Run:  magma -b qsimple_45_checks.m < /dev/null > qsimple_45_checks.log
// Every value quoted in the proof is ASSERTED below; the log ends with ALL_CHECKS_PASSED.
// A failed assertion stops the script before that line (magma -b still exits with status 0,
// so the log, not the exit status, is the verdict).  No randomized step affects the output.
//   Part A  Sections 1-3: torsion lists, residues modulo 3 and 5, the Weil bound.
//   Part B  Section 5: the flex model, the 3-division field and the dual 3-isogeny.
//   Part C  Section 6: the primitive cyclotomic equation (Proposition 5).
//   Part D  Section 7: X_1(15), the curve W and the four families of quadratic points.
//   Part E  Section 9: consistency checks (not used in the proof).
SetColumns(0);
SetMemoryLimit(16*10^9);
SetSeed(1);
// class groups are computed with the Minkowski bound (unconditionally) inside all descents
SetClassGroupBoundGenerators(map< PowerStructure(RngOrd) -> Integers() | O :-> Ceiling(MinkowskiBound(O)) >);
Q := Rationals();
S2S := func< s | "[" cat Join([Sprint(a) : a in s], ", ") cat "]" >;

// Explicit cube root of a nonzero rational function g, or false.  Factorization drops the
// leading unit, so the unit is recovered from the leading coefficients; h^3 = g is asserted.
function CubeRootRatFun(g)
  n := Numerator(g); dn := Denominator(g);
  okc, c := IsPower(Q!(LeadingCoefficient(n)/LeadingCoefficient(dn)), 3);
  if not okc then return false, _; end if;
  fn := Factorization(n); fd := Factorization(dn);
  if not &and[f[2] mod 3 eq 0 : f in fn cat fd] then return false, _; end if;
  P := Parent(n);
  r := (P!c)*&*([P!1] cat [f[1]^(f[2] div 3) : f in fn]);
  s := &*([P!1] cat [f[1]^(f[2] div 3) : f in fd]);
  h := r/s;
  assert h^3 eq g;
  return true, h;
end function;

//====================================================================================
print "PART A: torsion lists, residues and the Weil bound (Sections 1-3)";
// torsion groups of elliptic curves over quadratic fields (Kenku-Momose, Kamienny), as [n1, n2], n1 | n2
KMK := [[1, m] : m in [1..16] cat [18]] cat [[2, 2*m] : m in [1..6]] cat [[3, 3], [3, 6], [4, 4]];
assert #KMK eq 26;
hasorder := func< G, n | G[2] mod n eq 0 >;      // Z/n1 x Z/n2 with n1 | n2 has an element of order n iff n | n2
assert [G : G in KMK | hasorder(G, 45)] eq [];
assert [G : G in KMK | hasorder(G, 15)] eq [[1, 15]];
assert [G : G in KMK | hasorder(G, 9)] eq [[1, 9], [1, 18]];
Mazur := [[1, m] : m in [1..10] cat [12]] cat [[2, 2*m] : m in [1..4]];
assert [G : G in Mazur | hasorder(G, 15)] eq [];
print "A1 over quadratic fields no element of order 45, Z/15 is the only group with an element of order 15; over Q none of order 15: OK";
F5 := GF(5); F3 := GF(3); PF5<Y5> := PolynomialRing(F5); PF3<Y3> := PolynomialRing(F3);
assert not IsSquare(F5!(-3)) and not IsSquare(F5!(-1)/F5!3) and IsSquare(F5!(-1));
assert IsIrreducible(Y5^2 - Y5 + 1) and IsIrreducible(Y5^2 - Y5 + 1/F5!3);
rts := [rr[1] : rr in Roots(Y5^2 - Y5 + 1/F5!2)];
assert #rts eq 2 and Order(rts[1]/rts[2]) eq 4;
assert not IsSquare(F3!(-1)) and IsIrreducible(Y3^2 - Y3 + 1/F3!2);
RW := RealField(30);
assert (1 + Sqrt(RW!2))^4 lt 34 and Abs((1 + Sqrt(RW!2))^4 - (17 + 12*Sqrt(RW!2))) lt 10^-20;
print "A2 mod 5: -3, -1/3 non-squares, -1 a square, eigenvalue ratio of X^2-X+1/2 of order 4; mod 3: X^2-X+1/2 irreducible; (1+sqrt2)^4 = 17+12sqrt2 < 34: OK";

//====================================================================================
print "PART B: flex model, 3-division field and dual 3-isogeny (Section 5)";
FA<a1,a3> := FunctionField(Q, 2);
EA := EllipticCurve([FA | a1, 0, a3, 0, 0]);
assert Discriminant(EA) eq a3^3*(a1^3 - 27*a3);
RA<X> := PolynomialRing(FA);
hA := 3*X^3 + a1^2*X^2 + 3*a1*a3*X + 3*a3^2;
assert DivisionPolynomial(EA, 3) eq X*hA;
PA := EA![0, 0];
assert PA ne EA!0 and 3*PA eq EA!0 and -PA eq EA![0, -a3];
// the roots of h are -(e^2 + a1 e + a1^2)/9 with e^3 = a1^3 - 27 a3
RE<e> := PolynomialRing(FA);
AE := quo< RE | e^3 - (a1^3 - 27*a3) >;
eA := AE!e;
xe := -(eA^2 + a1*eA + a1^2)/9;
assert 3*xe^3 + a1^2*xe^2 + 3*a1*a3*xe + 3*a3^2 eq 0;
print "B1 y^2+a1xy+a3y=x^3: Disc = a3^3(a1^3-27a3), psi_3 = x(3x^3+a1^2x^2+3a1a3x+3a3^2), -P = (0,-a3), h(-(e^2+a1e+a1^2)/9) = 0 if e^3 = a1^3-27a3: OK";
// E' = E/<P> and the dual isogeny phihat : E' -> E
E2A, phiA := IsogenyFromKernel(EA, X);
assert aInvariants(E2A) eq [a1, 0, a3, -5*a1*a3, -a1^3*a3 - 7*a3^2];
phA := DualIsogeny(phiA);
assert Domain(phA) eq E2A and Codomain(phA) eq EA;
assert phiA*phA eq MultiplicationByMMap(EA, 3);      // phihat o phi = [3]
FE2 := FunctionField(E2A);
PE2 := Ambient(E2A);
xp := FE2!(PE2.1/PE2.3); yp := FE2!(PE2.2/PE2.3);
dpA := DefiningPolynomials(phA);
assert FE2!(dpA[1]/dpA[3]) eq (xp^3 - 9*a1*a3*xp - a1^3*a3 - 27*a3^2)/(3*xp + a1^2)^2;
assert FE2!(dpA[2]/dpA[3]) eq ((yp - 4*a3)/(3*xp + a1^2))^3;
print "B2 E/<P> = [a1,0,a3,-5a1a3,-a1^3a3-7a3^2]; phihat o phi = [3]; x o phihat = (x^3-9a1a3x-a1^3a3-27a3^2)/(3x+a1^2)^2, y o phihat = ((y-4a3)/(3x+a1^2))^3: OK";

//====================================================================================
print "PART C: the primitive cyclotomic equation (Section 6, Proposition 5)";
QT<T> := PolynomialRing(Q);
L<z> := NumberField(T^4 + T^3 + T^2 + T + 1);
OL := MaximalOrder(L);
tau := -(z^2 + z^3);
assert tau^2 - tau - 1 eq 0;
assert Discriminant(OL) eq 125;
RR30 := RealField(30);
MB := (Factorial(4)/4^4)*(4/Pi(RR30))^2*Sqrt(RR30!125);
assert MB lt 2 and MinkowskiBound(OL) lt 2;
assert ClassNumber(OL) eq 1;
UL, mUL := UnitGroup(OL);
assert Invariants(UL) eq [10, 0];
assert Order((OL!(-z))@@mUL) eq 10;
assert sub< UL | (OL!(-z))@@mUL, (OL!tau)@@mUL > eq UL;
print "C1 L = Q(zeta_5): Minkowski bound", Sprint(RealField(5)!MB), "< 2, class number 1, unit group = <-zeta> x <tau>: OK";
s2 := hom< L -> L | z^2 >; s3 := hom< L -> L | z^3 >;
assert s2(tau) eq -1/tau and s3(tau) eq -1/tau and tau*s2(tau)*s3(tau) eq 1/tau;
PAB<Av,Bv> := PolynomialRing(L, 2);
assert &*[Av + Bv*z^k : k in [1..4]] eq Av^4 - Av^3*Bv + Av^2*Bv^2 - Av*Bv^3 + Bv^4;
PZ<ai,bi> := PolynomialRing(Q, 2);
Phi10 := func< a, b | a^4 - a^3*b + a^2*b^2 - a*b^3 + b^4 >;
Phi5 := func< a, b | a^4 + a^3*b + a^2*b^2 + a*b^3 + b^4 >;
assert Phi5(ai, bi) eq Phi10(-ai, bi);
pi := 1 - z;
assert pi^4 eq 5*z^2/tau^2;
assert pi eq 5*tau*(z^4/(tau*pi))^3;
assert (1 - z)/(5*tau) eq ((2*z^3 - z^2 + z + 3)/5)^3;
print "C2 sigma_2(tau) = sigma_3(tau) = -1/tau; N(a+b*zeta) = Phi10(a,b); Phi5(a,b) = Phi10(-a,b); pi^4 = 5zeta^2/tau^2; pi = 5tau(zeta^4/(tau pi))^3: OK";
P2L := 2*OL;
assert IsPrime(P2L);
k16, red := ResidueClassField(P2L);
assert #k16 eq 16;
tb := red(OL!tau);
assert tb ne 1 and tb^3 eq 1;                      // tau mod 2 has order 3: not a cube in F_16^*
cubeclass := func< w | [j : j in [0..2] | (w/tb^j)^5 eq 1] >;   // the cubes of F_16^* are its 5-torsion
assert cubeclass(red(OL!1)) eq [0] and cubeclass(red(OL!z)) eq [0] and cubeclass(red(OL!(1 + z))) eq [1];
assert cubeclass(red(OL!5)) eq [0];
print "C3 2 is inert; modulo 2 and cubes: 1 ~ tau^0, zeta ~ tau^0, 1+zeta ~ tau^1, 5 ~ tau^0: OK";
F11 := Factorization(11*OL);
assert #F11 eq 4 and &and[Norm(f[1]) eq 11 : f in F11];
Ecyc := [ EllipticCurve([L | 0, 0, 4*z^3 + 3*z^2 + 7*z + 6, 0, -10*z^3 - 5*z^2 - 10*z]),
          EllipticCurve([L | 0, 0, z^3 + 7*z^2 + 13*z + 9, 0, -35*z^3 - 55*z^2 - 35*z]) ];
P2pl<tt,YY,WW> := ProjectiveSpace(L, 2);
P1 := ProjectiveSpace(Q, 1);
for j in [0, 1] do
  Ej := Ecyc[j + 1];
  C := Curve(P2pl, YY^3 - tau^j*(tt + z*WW)*(tt + z^2*WW)*(tt + z^3*WW));
  assert IsNonsingular(C) and Genus(C) eq 1;
  Em, mm := EllipticCurve(C, C![-z, 0, 1]);
  ok, iso := IsIsomorphic(Ej, Em);
  assert ok;
  assert Norm(Discriminant(Ej)) eq 3^12*5^12;      // good reduction at the four primes above 11
  lo, hi := RankBounds(Ej);
  assert lo eq 1 and hi eq 1;
  G, gm, b1, b2 := MordellWeilGroup(Ej);
  assert b1 and b2;                                // rank and generators proved: full Mordell-Weil group
  assert Invariants(G) eq [3, 0];
  fmap := iso*Inverse(mm);                         // E_j -> D_j
  pp := DefiningPolynomials(fmap);
  cov := map< Ej -> P1 | [pp[1], pp[3]] >;         // E_j -> D_j -> P^1, (t:Y:W) -> (t:W)
  pq := DefiningPolynomials(mm*Inverse(iso)*cov);
  assert pq[1]*C.3 - pq[2]*C.1 in Ideal(C);        // the composite D_j -> P^1 is (t:W)
  gmap := map< G -> Ej | a :-> gm(a) >;
  nn, vv, rr, cc := Chabauty(gmap, cov, 11);
  assert nn eq #vv and #cc[2] eq 0;                // every coset resolved, no unresolved neighbourhood
  assert rr in {1, 2};                             // saturation required at rr; b2 gives the full group
  covx := Extend(cov);
  params := {covx(gm(g)) : g in vv};
  expected := j eq 0 select {P1![0, 1], P1![1, 0]} else {P1![1, 1], P1![-1, 1]};
  assert params eq expected;
  print "C4 j =", j, ": D_j smooth of genus 1 with Jacobian E_j =", S2S(aInvariants(Ej)),
        "; Norm(Disc E_j) = 3^12*5^12; rank 1; E_j(L) = Z/3 + Z (proved); Chabauty at 11: upper bound",
        nn, "= number of points found, no unresolved cosets, saturation required at", rr,
        "; rational t-values", Sprint(params);
end for;
print "C5 Proposition 5: rational t-values are 0, oo on D_0 and 1, -1 on D_1: OK";

//====================================================================================
print "PART D: X_1(15) and its quadratic points (Section 7)";
// C1 the Tate normal form E_d
Fd<d> := FunctionField(Q);
Ed := EllipticCurve([Fd | 1 - d, -d, -d, 0, 0]);
assert Discriminant(Ed) eq d^5*(d^2 - 11*d - 1);
P5 := Ed![0, 0];
assert P5 ne Ed!0 and 5*P5 eq Ed!0;
Rd<XX> := PolynomialRing(Fd);
psi3 := 3*XX^4 + (d^2 - 6*d + 1)*XX^3 + 3*(d^2 - d)*XX^2 + 3*d^2*XX - d^3;
assert DivisionPolynomial(Ed, 3) eq psi3;
assert Evaluate(psi3, 0) eq -d^3 and Evaluate(psi3, d) eq d^5;
Fdu<dd,uu> := FunctionField(Q, 2);
zz := dd*uu/(uu - 1);
ps := func< w | 3*w^4 + (dd^2 - 6*dd + 1)*w^3 + 3*(dd^2 - dd)*w^2 + 3*dd^2*w - dd^3 >;
assert ps(dd*uu)/dd^3 eq (uu - 1)^2*(uu*zz^2 + 3*uu*(uu - 1)*zz + (uu - 1));
assert (uu - 1)*zz/uu eq dd and (uu - 1)*zz eq dd*uu;
print "D1 Disc(E_d) = d^5(d^2-11d-1); (0,0) has order 5; psi_3 as stated; psi_3(0) = -d^3, psi_3(d) = d^5; psi_3(du)/d^3 = (u-1)^2(uz^2+3u(u-1)z+u-1) for z = du/(u-1): OK";
// C2 the affine curve (8)
A3<u,zc,y> := AffineSpace(Q, 3);
f1 := u*zc^2 + 3*u*(u - 1)*zc + (u - 1);
dz := (u - 1)*zc;                                  // = u*d = x_3
f2 := u^2*y^2 + (u - dz)*u*dz*y - dz*u*y - (u^2*dz^3 - u*dz*dz^2);
FU := FieldOfFractions(CoordinateRing(A3));
dF := FU!dz/FU!u;
assert FU!f2 eq FU!u^2*(FU!y^2 + (1 - dF)*FU!dz*FU!y - dF*FU!y - FU!dz^3 + dF*FU!dz^2);
Xa := Curve(A3, [f1, f2]);
assert IsIrreducible(Xa) and Dimension(Xa) eq 1;
Sg := SingularSubscheme(Xa);
assert Dimension(Sg) eq 0;
spts := PointsOverSplittingField(Sg);
assert #spts eq 3;
sdv := [];
for p in spts do
  uv := p[1]; zv := p[2];
  assert uv ne 0;
  dv := (uv - 1)*zv/uv;
  assert dv eq 0 or dv^2 - 11*dv - 1 eq 0;        // every singular point is a cusp
  if dv eq 0 then assert [p[1], p[2], p[3]] eq [1, 0, 0]; end if;
  Append(~sdv, dv eq 0 select "(1,0,0), d = 0" else "d^2-11d-1 = 0");
end for;
print "D2 the curve (8) is irreducible; its 3 affine singular points are cusps:", S2S(sdv);
// C3 the elliptic curve W = 15a8 and the birational map W -> (8)
W := EllipticCurve([Q | 1, 1, 1, 0, 0]);
lo, hi := RankBounds(W);
assert lo eq 0 and hi eq 0;
TW := TorsionSubgroup(W);
assert Invariants(TW) eq [4];
OW := W!0; T1 := W![0, 0]; T2 := W![-1, 0]; T3 := W![0, -1];
assert 2*T1 eq T2 and 3*T1 eq T3 and 4*T1 eq OW;
FW := FunctionField(W);
PW := Ambient(W);
xw := FW!(PW.1/PW.3); yw := FW!(PW.2/PW.3);
// the moduli functions on W (the formulas displayed in Section 7)
function UofW(x, yy) return (x + 1)/(x^2 + x + 1); end function;
function ZofW(x, yy) return x*(x^2 - 1 - (x + 2)*yy)/((x + 1)*(x^2 + x + 1)); end function;
function DofW(x, yy) return x^3*(1 - x^2 + (x + 2)*yy)/((x + 1)^2*(x^2 + x + 1)); end function;
function Y3ofW(x, yy) return x^5*(x - 2*x^2 - 2*x^3 + (x^2 + 4*x + 1)*yy)/((x + 1)^2*(x^2 + x + 1)^3); end function;
uW := UofW(xw, yw); zW := ZofW(xw, yw); dW := DofW(xw, yw); y3W := Y3ofW(xw, yw);
x3W := dW*uW;
assert dW eq (uW - 1)*zW/uW and x3W eq (uW - 1)*zW;
assert Evaluate(f1, [uW, zW, y3W]) eq 0 and Evaluate(f2, [uW, zW, y3W]) eq 0;
xinv := ((uW - 1)^3*(uW + 1)^2 + uW^2*(uW - 1)^2*zW + uW^2*(-9*uW^3 + 3*uW^2 + 7*uW + 1)*y3W
          - uW^3*(3*uW + 2)*zW*y3W)/(uW*(uW - 1)^2*(uW^2 - uW - 1));
yinv := -((uW - 1)^4 + (uW - 1)^3*zW + 3*uW^2*(uW - 1)*y3W + uW^2*zW*y3W)/(uW*(uW - 1)^3);
assert xinv eq xw and yinv eq yw;                  // rho o pi = id: pi is birational
assert y3W^2 + (1 - dW)*x3W*y3W - dW*y3W eq x3W^3 - dW*x3W^2;
assert 3*x3W^4 + (dW^2 - 6*dW + 1)*x3W^3 + 3*(dW^2 - dW)*x3W^2 + 3*dW^2*x3W - dW^3 eq 0;
print "D3 W = [1,1,1,0,0]: rank 0, W(Q) = Z/4 = {O, T1 = (0,0), T2 = (-1,0) = 2T1, T3 = (0,-1) = 3T1}; (u,z,y3) maps W to (8) and (x,y) inverts it: OK";
// C4 the flex coefficient: translate (x3,y3) to (0,0), then shear y -> y + s x
b1c := 1 - dW; b2c := -dW; b3c := -dW; b4c := FW!0; b6c := FW!0; r := x3W; tt3 := y3W;
a3t := b3c + r*b1c + 2*tt3;
a3W := 2*y3W + (1 - dW)*x3W - dW;
assert a3t eq a3W and a3W ne 0;
sh := (b4c + 2*r*b2c - tt3*b1c + 3*r^2)/a3t;
a2n := b2c - sh*b1c + 3*r - sh^2;
a4n := b4c - sh*b3c + 2*r*b2c - (tt3 + r*sh)*b1c + 3*r^2 - 2*sh*tt3;
a6n := b6c + r*b4c + r^2*b2c + r^3 - tt3*b3c - tt3^2 - r*tt3*b1c;
assert a2n eq 0 and a4n eq 0 and a6n eq 0;
a1n := b1c + 2*sh;
assert a3W^3*(a1n^3 - 27*a3W) eq dW^5*(dW^2 - 11*dW - 1);   // Disc of the flex model = Disc(E_d)
deltaW := dW^2*(dW^2 - 11*dW - 1);
assert dW^5*(dW^2 - 11*dW - 1)/deltaW eq dW^3;     // Disc = d^3 * delta: same cube class
print "D4 after translation and shear E_d becomes y^2+a1xy+a3y=x^3 with a3 = 2y3+(1-d)x3-d != 0 and the same discriminant d^5(d^2-11d-1) = d^3*delta: OK";
// C5 the four rational points are cusps
plO := Place(OW); pl1 := Place(T1); pl2 := Place(T2); pl3 := Place(T3);
assert Valuation(dW, plO) eq -3 and Valuation(dW, pl2) eq -3;
assert Valuation(dW, pl1) eq 3 and Valuation(dW, pl3) eq 3;
print "D5 d has poles of order 3 at O and T2 and zeros of order 3 at T1 and T3: the rational points of W are cusps";
// C6 the four families
tW := [ yw/(xw + yw + 1), xw + 1, yw/xw, -(xw + yw + 1)/xw ];
DivOf := func< pl | Divisor(pl) >;
assert Divisor(tW[1]) eq 2*DivOf(pl1) - 2*DivOf(pl3);
assert Divisor(tW[2]) eq 2*DivOf(pl2) - 2*DivOf(plO);
assert Divisor(tW[3]) eq DivOf(pl1) + DivOf(pl2) - DivOf(pl3) - DivOf(plO);
assert Divisor(tW[4]) eq DivOf(pl2) + DivOf(pl3) - DivOf(pl1) - DivOf(plO);
assert Evaluate(tW[1], plO) eq 1 and Evaluate(tW[1], pl2) eq 1;
assert Evaluate(tW[2], pl1) eq 1 and Evaluate(tW[2], pl3) eq 1;
assert [Degree(tW[i]) : i in [1..4]] eq [2, 2, 2, 2];
ynegw := -yw - xw - 1;                             // the negation (x,y) -> (x,-y-x-1) of W
assert ynegw/xw eq tW[4] and ynegw/(xw + ynegw + 1) eq 1/tW[1] and xw + 1 eq tW[2];
print "D6 div t_1 = 2T1-2T3, div t_2 = 2T2-2O, div t_3 = T1+T2-T3-O, div t_4 = T2+T3-T1-O; t_1(O) = t_1(T2) = 1, t_2(T1) = t_2(T3) = 1; each t_j has degree 2; negation maps t_3 to t_4 and t_1 to 1/t_1: OK";
// C7 traces and norms in Q(W)/Q(t_j), in explicit quadratic extensions of Q(t)
Ft<t> := FunctionField(Q);
Pt := PolynomialRing(Q); tp := Pt.1;
RtZ<Zv> := PolynomialRing(Ft);
Bq := tp^2 - tp + 1; Q10 := tp^4 - tp^3 + tp^2 - tp + 1; Q5 := tp^4 + tp^3 + tp^2 + tp + 1;
F12 := (tp^2 - tp + 1)*(tp^4 - 11*tp^3 + 21*tp^2 - 11*tp + 1);
F34 := (tp^2 - tp + 1)^2*(tp + 1)*(tp^4 + tp^3 - 9*tp^2 + tp + 1);
Npaper := [
  -2*tp^10 - 17*tp^9 + 60*tp^8 - 30*tp^7 - 126*tp^6 + 228*tp^5 - 126*tp^4 - 30*tp^3 + 60*tp^2 - 17*tp - 2,
  2*tp^10 + 8*tp^9 - 24*tp^8 + 30*tp^7 - 45*tp^6 + 60*tp^5 - 45*tp^4 + 30*tp^3 - 24*tp^2 + 8*tp + 2,
  -tp^22 - 6*tp^21 - 7*tp^20 + 16*tp^19 + 21*tp^18 - 44*tp^17 - 37*tp^16 + 88*tp^15 + 38*tp^14 - 132*tp^13
    - 17*tp^12 + 152*tp^11 - 17*tp^10 - 132*tp^9 + 38*tp^8 + 88*tp^7 - 37*tp^6 - 44*tp^5 + 21*tp^4
    + 16*tp^3 - 7*tp^2 - 6*tp - 1,
  -tp^22 - 5*tp^21 - 2*tp^20 + 15*tp^19 - 8*tp^18 - 37*tp^17 + 55*tp^16 + 20*tp^15 - 119*tp^14 + 89*tp^13
    + 73*tp^12 - 170*tp^11 + 73*tp^10 + 89*tp^9 - 119*tp^8 + 20*tp^7 + 55*tp^6 - 37*tp^5 - 8*tp^4
    + 15*tp^3 - 2*tp^2 - 5*tp - 1 ];
Dpaper := [ Bq^3*Q10, Bq^3*Q10, tp^9*Q5, tp^9*Q5 ];
Dname := [ "B(t)^3 Phi10(t)", "B(t)^3 Phi10(t)", "t^9 Phi5(t)", "t^9 Phi5(t)" ];
Fpaper := [ F12, F12, F34, F34 ];
Qcyc := [ Q10, Q10, Q5, Q5 ];
// the generic point of W over Q(t_j): a quadratic relation and the other coordinate
qrel := [ Zv^2 - t/(1 - t)^2*Zv - t/(1 - t)^2,    // family 1: (1-t)^2 x^2 - t x - t = 0, y = t(x+1)/(1-t)
          Zv^2 + t*Zv - t*(t - 1)^2,              // family 2: y^2 + t y - t(t-1)^2 = 0,   x = t - 1
          Zv^2 + (1 - t - t^2)*Zv - t,            // family 3: x^2 + (1-t-t^2) x - t = 0,  y = t x
          Zv^2 + (1 - t - t^2)*Zv - t ];          // family 4: the same relation,          y = -(t+1) x - 1
gcof := [* *];
for j in [1..4] do
  q := qrel[j];
  assert not IsSquare(Discriminant(q));            // Q(W) is a quadratic extension of Q(t_j)
  Kj<Xj> := ext< Ft | q >;
  case j:
    when 1: xs := Xj; ys := t*(Xj + 1)/(1 - t); tv := ys/(xs + ys + 1);
    when 2: xs := Kj!(t - 1); ys := Xj; tv := xs + 1;
    when 3: xs := Xj; ys := t*Xj; tv := ys/xs;
    when 4: xs := Xj; ys := -(t + 1)*Xj - 1; tv := -(xs + ys + 1)/xs;
  end case;
  assert ys^2 + xs*ys + ys eq xs^3 + xs^2;           // the generic point lies on W
  assert tv eq Kj!t;                                // and has t_j-value t
  dK := DofW(xs, ys); x3K := dK*UofW(xs, ys); y3K := Y3ofW(xs, ys);
  a3K := 2*y3K + (1 - dK)*x3K - dK;
  delK := dK^2*(dK^2 - 11*dK - 1);
  Sj := Trace(a3K); Nmj := Norm(a3K);
  Rj := Sj^2/Nmj - 2;
  assert Rj eq Evaluate(Npaper[j], t)/Evaluate(Dpaper[j], t);   // R_j = N_j/Den_j exactly
  assert &and[Denominator(c) eq 1 : c in Coefficients(Npaper[j])];
  assert Degree(Npaper[j]) eq [10, 10, 22, 22][j] and Degree(Dpaper[j]) eq [10, 10, 13, 13][j];
  assert Coefficients(Npaper[j]) eq Reverse(Coefficients(Npaper[j]));   // palindromic
  assert Resultant(Qcyc[j], Npaper[j]) eq 25;
  okg, gj := CubeRootRatFun(Norm(delK)/Evaluate(Fpaper[j], t));
  assert okg;                                       // N(delta) = F_j * g_j^3
  Append(~gcof, gj);
  print "D7 family", j, ": Tr(a3)^2/N(a3) - 2 = N_j/(" cat Dname[j] cat ") with N_j as in the text (degree",
        Degree(Npaper[j]), ", integral, palindromic); Res(cyclotomic quartic, N_j) = 25; N(delta) = F_j * g^3 with g =", Sprint(gj);
end for;
assert Npaper[1] ne Npaper[2] and Npaper[3] ne Npaper[4];
assert Evaluate(F12, -1) eq 135 and Evaluate(F34, 1) eq -10;
assert not IsPower(Q!135, 3) and not IsPower(Q!(-10), 3);
print "D8 F12(-1) = 135 and F34(1) = -10, neither a rational cube";
// C9 the fibres over t = 0, 1, -1, oo
function CuspType(pl)
  if Valuation(dW, pl) lt 0 then return "pole of d"; end if;
  dv := Evaluate(dW, pl);
  if dv eq 0 then return "d = 0"; end if;
  if dv^2 - 11*dv - 1 eq 0 then return "d^2-11d-1 = 0"; end if;
  return "not a cusp";
end function;
excl := [ Q!(-1), Q!(-1), Q!1, Q!1 ];
exdisc := [ -15, -15, 5, 5 ];
exd := [ Q!(-1/8), Q!8, Q!(1/2), Q!(-2) ];
for j in [1..4] do
  f := tW[j];
  cuspvals := [];
  for t0 in [* Q!0, Q!1, Q!(-1), "oo" *] do
    if Type(t0) eq MonStgElt then
      fib := Poles(f); mult := [-Valuation(f, pl) : pl in fib];
    else
      fib := Zeros(f - t0); mult := [Valuation(f - t0, pl) : pl in fib];
    end if;
    assert &+[Degree(fib[i])*mult[i] : i in [1..#fib]] eq 2;
    if Type(t0) ne MonStgElt and t0 eq excl[j] then
      // the exceptional fibre: one non-cuspidal point of degree 2
      assert #fib eq 1 and Degree(fib[1]) eq 2 and CuspType(fib[1]) eq "not a cusp";
      pl := fib[1];
      K := ResidueClassField(pl);
      assert Discriminant(MaximalOrder(K)) eq exdisc[j];
      dv := Evaluate(dW, pl); x3v := Evaluate(x3W, pl); y3v := Evaluate(y3W, pl);
      assert dv eq exd[j];                           // d is rational at this point
      EK := EllipticCurve([K | 1 - dv, -dv, -dv, 0, 0]);
      assert Order(EK![0, 0]) eq 5 and Order(EK![x3v, y3v]) eq 3;
      assert Invariants(TorsionSubgroup(EK)) eq [15];
      sig := [g : g in Automorphisms(K) | g(K.1) ne K.1][1];
      a3v := 2*y3v + (1 - dv)*x3v - dv;
      assert a3v/sig(a3v) eq -1;                     // condition (4) holds here
      nD := Q!Norm(Discriminant(EK));
      assert not IsPower(nD, 3);                     // condition (3) fails
      assert IsPower(nD/Evaluate(Fpaper[j], t0), 3);
      gv := Evaluate(gcof[j], t0);
      assert gv ne 0 and gv eq [Q!5/768, Q!80/3, Q!(-5/8), Q!(-10)][j];   // finite and nonzero
      assert Evaluate(Npaper[j], t0)/Evaluate(Dpaper[j], t0) eq -2;       // R_j(t_0) = -2, so r = -2
      print "D9 family", j, ": over t =", Sprint(t0), "one non-cuspidal point, residue field of discriminant", exdisc[j],
            ", d =", Sprint(Q!dv), ", E_d(K)_tors = Z/15, a3/sigma(a3) = -1 (so (4) holds), N(Disc) =", Sprint(nD),
            "is not a cube (so (3) fails), N(Disc)/F_j(t) is a cube, g(t) =", Sprint(gv), ", R_j(t) = -2";
    else
      types := [CuspType(pl) : pl in fib];
      assert not "not a cusp" in types;
      Append(~cuspvals, Sprint(t0));
      print "D9 family", j, ": over t =", Sprint(t0), "places of degrees", S2S([Degree(pl) : pl in fib]),
            "with multiplicities", S2S(mult), ", all cusps:", S2S(types);
    end if;
  end for;
  assert cuspvals eq (j le 2 select ["0", "1", "oo"] else ["0", "-1", "oo"]);
end for;
print "D10 Proposition 6: no non-cuspidal quadratic point of X_1(15) satisfies both (3) and (4): OK";

//====================================================================================
print "PART E: consistency checks (Section 9; not used in the proof)";
Qx<xq> := PolynomialRing(Q);
// census curves 9.c, 18.c, 15.c: y^2 + h y = f; the census map of degree 3,
// (x,Y) -> (p/q, Y*hh/q^2) with Y = 2y + h, onto v^2 = c3 u^3 + c2 u^2 + c1 u + c0 over K
cases := [*
  <"9.c", "[9]", xq^2 + 1, xq^6 - 2*xq^5 + 2*xq^4 + 9*xq^3 - 12*xq^2 + 3*xq + 26, xq,
   ["-8*w", "90*w - 56", "36*w + 20", "4"], ["-1/2", "1", "1"], ["5*w - 1", "1/2*(-w - 8)", "-3*w - 1", "1"],
   ["1/4*(-11*w - 9)", "1/2*(-8*w - 3)", "1/2*(-3*w + 3)", "1"], [3], [true, true]>,
  <"18.c", "[18]", xq^2 - 6, -xq^5 + 5*xq^3 - 14*xq^2 - 16*xq + 10, xq^3 + xq,
   ["-20*w - 48", "-18*w - 41", "-2*w - 2", "1"], ["-2*w - 5", "-2*w + 2", "1"], ["2*w + 8", "2*w - 1", "-2", "1"],
   ["w + 1", "-6*w - 15", "-w - 3", "1"], [6], [true, true]>,
  <"15.c", "[15]", xq^2 - 3, xq^6 + 9*xq^5 + 38*xq^4 + 47*xq^3 + 72*xq^2 - 18*xq + 1, xq^2,
   ["-432*w + 864", "36*w + 45", "24*w - 12", "4"], ["0", "-w + 1", "1"], ["w - 2", "-7*w + 15", "1/2*(-4*w + 9)", "1"],
   ["-4*w + 7", "1/2*(5*w - 9)", "-w + 3", "1"], [15], [true, false]> *];
parents := [* *];
for cs in cases do
  K<w> := NumberField(cs[3]);
  KX<Xk> := PolynomialRing(K);
  ev := func< s | eval s >;
  c3, c2, c1, c0 := Explode([ev(s) : s in cs[6]]);
  pK := &+[ev(cs[7][i])*Xk^(i - 1) : i in [1..#cs[7]]];
  qK := &+[ev(cs[8][i])*Xk^(i - 1) : i in [1..#cs[8]]];
  hK := &+[ev(cs[9][i])*Xk^(i - 1) : i in [1..#cs[9]]];
  g := KX!(4*cs[4] + cs[5]^2);
  assert g*hK^2 eq qK*(c3*pK^3 + c2*pK^2*qK + c1*pK*qK^2 + c0*qK^3);   // the census map is well defined
  E := EllipticCurve([K | 0, c2, 0, c1*c3, c0*c3^2]);
  TE, mTE := TorsionSubgroup(E);
  assert Invariants(TE) eq cs[10];
  sig := [a : a in Automorphisms(K) | a(w) ne w][1];
  c3ok := IsPower(Q!Norm(Discriminant(E)), 3);
  P3 := [mTE(x) : x in TE | Order(x) eq 3];
  assert #P3 eq 2;
  c4ok := [IsPower(2*P[2]/sig(2*P[2]), 3) : P in P3];    // a3 = 2y(P) after moving P to (0,0)
  opp := [IsPower(2*P[2]*sig(2*P[2]), 3) : P in P3];
  assert [c3ok, c4ok[1]] eq cs[11] and c4ok[1] eq c4ok[2];
  assert opp eq [false, false];
  // type (N): E has no CM, and some split prime has a_p(E)^2 <> a_p(E^sigma)^2 = a_{sigma p}(E)^2
  assert not HasComplexMultiplication(E);
  OK := MaximalOrder(K);
  found := 0;
  for p in PrimesUpTo(200) do
    Fp := Factorization(p*OK);
    if #Fp eq 2 and &and[Valuation(Discriminant(E), f[1]) eq 0 and
                         &and[Valuation(c, f[1]) ge 0 : c in aInvariants(E)] : f in Fp] then
      if TraceOfFrobenius(E, Fp[1][1])^2 ne TraceOfFrobenius(E, Fp[2][1])^2 then found := p; break; end if;
    end if;
  end for;
  assert found gt 0;
  Append(~parents, <K, E>);
  print "E1", cs[1], cs[2], ": K of discriminant", Discriminant(MaximalOrder(K)), ", E(K)_tors =", S2S(cs[10]),
        ", E without CM and not geometrically isogenous to E^sigma (traces at the primes above", found,
        "), (3):", c3ok, ", (4):", c4ok[1], ", opposite sign a3*sigma(a3) a cube:", opp[1];
end for;
// family 1 over t = 2 gives the parent of 15.c: (3) holds and (4) fails
K := parents[3][1]; E15 := parents[3][2];
w := K.1;
assert w^2 eq 3;
isos := [];
for xx in [1 + w, 1 - w] do
  yy := -2*(xx + 1);
  assert yy^2 + xx*yy + yy eq xx^3 + xx^2 and yy/(xx + yy + 1) eq 2;
  dv := DofW(xx, yy); x3v := dv*UofW(xx, yy); y3v := Y3ofW(xx, yy);
  EK := EllipticCurve([K | 1 - dv, -dv, -dv, 0, 0]);
  assert Order(EK![x3v, y3v]) eq 3 and Invariants(TorsionSubgroup(EK)) eq [15];
  sig := hom< K -> K | -w >;
  a3v := 2*y3v + (1 - dv)*x3v - dv;
  assert IsPower(Q!Norm(Discriminant(EK)), 3) and not IsPower(a3v/sig(a3v), 3);
  Append(~isos, IsIsomorphic(EK, E15));
end for;
assert isos eq [false, true];
print "E2 family 1, t = 2: the point (1 - sqrt3, -2(2 - sqrt3)) of W gives E_d with E_d(K)_tors = Z/15, K-isomorphic to the elliptic curve of 15.c; (3) holds, (4) fails";
// the exceptional fibres come from elliptic curves over Q
for dv in [Q!(-1/8), Q!8, Q!(1/2), Q!(-2)] do
  E0 := EllipticCurve([Q | 1 - dv, -dv, -dv, 0, 0]);
  M := aInvariants(MinimalModel(E0));
  assert M eq (dv in {Q!(-1/8), Q!8} select [1, 1, 1, 22, -9] else [1, 1, 1, -3, 1]);
  assert Conductor(E0) eq 50;
  lab := CremonaReference(MinimalModel(E0));
  assert lab eq (dv in {Q!(-1/8), Q!8} select "50b2" else "50b1");
  print "E3 d =", Sprint(dv), ": E_d is defined over Q with minimal model", S2S(M), "(Cremona", lab cat ")";
end for;

print "ALL_CHECKS_PASSED";
quit;
