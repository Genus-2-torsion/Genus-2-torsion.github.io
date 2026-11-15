// Certificate for "Rational 3-torsion and the 2-primary torsion of geometrically split abelian surfaces"
// (census proof qsimple-96-obstruction).  Run: magma -b qsimple_96_checks.m < /dev/null
// Every claim is an assert; the run ends with ALL_CHECKS_PASSED.
// Parts: A (Lemma 1), B, C, D, F (Lemma 2), E, G (Lemma 3 and Remark 5), H (Lemma 7, Proposition 8 and Proposition 9),
// S (Remark 2).
SetColumns(0);
SetMemoryLimit(16*10^9);
Q := Rationals();
Z := Integers();
Qx<x> := PolynomialRing(Q);

// ---------------------------------------------------------------------------
// Part A.  #B(F_5) is never 96 for an abelian surface B over F_5.
// #B(F_5) = h(6) for h(y) = y^2 - s*y + t in Z[y] with both roots real in [-2 sqrt 5, 2 sqrt 5].
// Roots real and in the interval  <=>  s^2 >= 4t,  s^2 <= 80,  20 + t >= 0  and  20 s^2 <= (20 + t)^2.
vals := {};
for s in [-9..9] do for t in [-20..20] do
  if s^2 - 4*t ge 0 and s^2 le 80 and 20 + t ge 0 and 20*s^2 le (20 + t)^2 then
    Include(~vals, 36 - 6*s + t);
  end if;
end for; end for;
// |s| <= 4 sqrt 5 < 9 and |t| <= 20 for such h, so the loop is exhaustive.
assert not 96 in vals;
assert Max(vals) eq 100;
printf "A: possible values of h(6): %o\n", Sort(SetToSequence(vals));
// For the remark after Lemma 1: over F_7 the value 96 does occur (a product of curves with 8 and 12 points).
orders7 := { #EllipticCurve([GF(7) | a, b]) : a, b in GF(7) | 4*a^3 + 27*b^2 ne 0 };
assert 8 in orders7 and 12 in orders7;

// ---------------------------------------------------------------------------
// Part B.  Kubert's normal form for a point of order 12.
Qt<t> := FunctionField(Q);
function kubert_bc(tau)     // Kubert, Table 3, entry for Z/12
  m := (3*tau - 3*tau^2 - 1)/(tau - 1); f := m/(1 - tau); d := m + tau;
  c := f*(d - 1); b := c*d;
  return b, c;
end function;
function fam_bc(t)          // the same in factored form
  return t*(2*t - 1)*(2*t^2 - 2*t + 1)*(3*t^2 - 3*t + 1)/(t - 1)^4, -t*(2*t - 1)*(3*t^2 - 3*t + 1)/(t - 1)^3;
end function;
function cusp_pol(t)        // vanishes exactly at the finite values of t other than 1 at which E_t is singular
  return t*(2*t - 1)*(6*t^2 - 6*t + 1)*(3*t^2 - 3*t + 1)*(2*t^2 - 2*t + 1);
end function;
bK, cK := kubert_bc(t); bt, ct := fam_bc(t);
assert bK eq bt and cK eq ct;
Et := EllipticCurve([1 - ct, -bt, -bt, 0, 0]);
P := Et![0, 0];
assert 12*P eq Et!0 and 6*P ne Et!0 and 4*P ne Et!0;      // (0,0) has exact order 12
dd := bt/ct; ff := ct/(dd - 1);
assert (dd - ff)/(1 - ff) eq t;                            // the parameter is recovered from (b, c)
c4t, c6t := Explode(cInvariants(Et)); Dt := Discriminant(Et); jt := jInvariant(Et);
assert Dt eq t^12*(2*t - 1)^6*(6*t^2 - 6*t + 1)*(3*t^2 - 3*t + 1)^4*(2*t^2 - 2*t + 1)^3/(t - 1)^24;
assert Denominator(c4t*(t - 1)^12) eq 1 and Denominator(c6t*(t - 1)^18) eq 1;
assert Degree(Numerator(c4t)) eq 16 and Degree(Numerator(c6t)) eq 24;
// j has poles exactly at t = 0, 1/2, 1, oo and at the roots of the three quadratics
assert GCD(Numerator(c4t), Numerator(Dt)) eq 1 and Evaluate(Numerator(c4t), 1) ne 0;
assert Degree(Numerator(jt)) eq 48 and Degree(Denominator(jt)) eq 46;
print "B: Kubert family: (0,0) has exact order 12; t = (d-f)/(1-f); Delta = t^12 (2t-1)^6 (6t^2-6t+1) (3t^2-3t+1)^4 (2t^2-2t+1)^3 / (t-1)^24; j has poles of order 12, 6, 12, 2 at t = 0, 1/2, 1, oo and 1, 4, 3 at the roots of 6t^2-6t+1, 3t^2-3t+1, 2t^2-2t+1";

// Universality.  Tate normal form E(b,c): y^2 + (1-c)xy - by = x^3 - bx^2 with the point P = (0,0).
// Division polynomials psi_n in Z[a1,...,a6,x,y] (Silverman, Exercise III.3.7): for a point P != O of an elliptic
// curve, n P = O  <=>  psi_n(P) = 0.  Write psi_n = g_n for n odd and psi_n = psi_2 g_n for n even, with g_n a
// polynomial in x, and psi_2^2 = F = 4x^3 + b2 x^2 + 2 b4 x + b6.  Put f_n = g_n (n odd), f_n = F g_n (n even).
// Then n P = O  <=>  f_n(x(P)) = 0, and f_n has coefficients in Z[a1,...,a6].  The g_n are computed here from the
// recursion of the exercise, in Z[b,c][x] for the curve E(b,c), and compared with Magma's DivisionPolynomial.
Fbc<bb, cc> := FunctionField(Q, 2);
Pbc := Parent(Numerator(bb));
Ebc := EllipticCurve([Fbc | 1 - cc, -bb, -bb, 0, 0]);
Pbcx<xx> := PolynomialRing(Pbc);
pb := Pbc!Numerator(bb); pc := Pbc!Numerator(cc);
A1 := 1 - pc; A2 := -pb; A3 := -pb; A4 := Pbc!0; A6 := Pbc!0;
B2 := A1^2 + 4*A2; B4 := 2*A4 + A1*A3; B6 := A3^2 + 4*A6; B8 := A1^2*A6 + 4*A2*A6 - A1*A3*A4 + A2*A3^2 - A4^2;
assert [Fbc | B2, B4, B6, B8] eq bInvariants(Ebc);
FF2 := 4*xx^3 + B2*xx^2 + 2*B4*xx + B6;                   // psi_2^2
gdiv := AssociativeArray();
gdiv[0] := Pbcx!0; gdiv[1] := Pbcx!1; gdiv[2] := Pbcx!1;
gdiv[3] := 3*xx^4 + B2*xx^3 + 3*B4*xx^2 + 3*B6*xx + B8;
gdiv[4] := 2*xx^6 + B2*xx^5 + 5*B4*xx^4 + 10*B6*xx^3 + 10*B8*xx^2 + (B2*B8 - B4*B6)*xx + (B4*B8 - B6^2);
for n in [5..12] do
  m := n div 2;
  if IsOdd(n) then        // psi_{2m+1} = psi_{m+2} psi_m^3 - psi_{m-1} psi_{m+1}^3
    if IsEven(m) then gdiv[n] := FF2^2*gdiv[m+2]*gdiv[m]^3 - gdiv[m-1]*gdiv[m+1]^3;
    else gdiv[n] := gdiv[m+2]*gdiv[m]^3 - FF2^2*gdiv[m-1]*gdiv[m+1]^3; end if;
  else                    // psi_2 psi_{2m} = psi_m (psi_{m+2} psi_{m-1}^2 - psi_{m-2} psi_{m+1}^2)
    gdiv[n] := gdiv[m]*(gdiv[m+2]*gdiv[m-1]^2 - gdiv[m-2]*gdiv[m+1]^2);
  end if;
end for;
Fbcx := PolynomialRing(Fbc);
function PsiAtOrigin(n)
  fn := IsOdd(n) select gdiv[n] else FF2*gdiv[n];
  assert Degree(fn) eq (IsOdd(n) select (n^2 - 1) div 2 else (n^2 + 2) div 2);
  assert LeadingCoefficient(fn) eq (IsOdd(n) select n else 2*n);
  assert Fbcx![ Fbc!co : co in Coefficients(fn) ] eq Fbcx!DivisionPolynomial(Ebc, n);     // agrees with Magma
  return Evaluate(fn, 0);
end function;
function IrrFactors(g)      // the set of monic irreducible factors (the constant factor plays no role)
  return { Pbc | fa[1] : fa in Factorization(g) };
end function;
assert Denominator(Discriminant(Ebc)) eq 1 and pb in IrrFactors(Pbc!Numerator(Discriminant(Ebc)));   // b = 0: singular
assert IrrFactors(PsiAtOrigin(4)) eq { pb, pc };                      // 4P = 0  <=>  c = 0
assert IrrFactors(PsiAtOrigin(5)) eq { pb, pb - pc };                 // 5P = 0  <=>  b = c
assert IrrFactors(PsiAtOrigin(6)) eq { pb, pb - pc^2 - pc };          // 6P = 0  <=>  b = c^2 + c
Phi12 := pb^4 - 1/3*pb^3*pc^2 - 3*pb^3*pc + 10/3*pb^2*pc^2 + 1/3*pb*pc^4 - 5/3*pb*pc^3 + 1/3*pc^6 + 1/3*pc^4;
assert IrrFactors(PsiAtOrigin(12)) eq { pb, pc, pb - pc^2 - pc, Phi12 };
// so: P has exact order 12  <=>  Phi12(b,c) = 0 and b c (b - c^2 - c) != 0   (given a nonzero discriminant)
assert Evaluate(Phi12, [bt, ct]) eq 0;                                // the family lies on Phi12 = 0
// The inverse map T = (d - f)/(1 - f), d = b/c, f = c/(d - 1), as a rational function of (b, c).
dB := bb/cc; fB := cc/(dB - 1); TB := (dB - fB)/(1 - fB);
allowed := { pb, pc, pb - pc, pb - pc^2 - pc };
function OnlyAllowed(g)     // all irreducible factors of the polynomial g are among b, c, b - c, b - c^2 - c
  return IrrFactors(Pbc!g) subset allowed;
end function;
assert TB eq (bb^2 - bb*cc - cc^3)/(cc*(bb - cc^2 - cc));
assert OnlyAllowed(Denominator(TB));
assert TB - 1 eq (bb - cc)^2/(cc*(bb - cc^2 - cc));                    // so T is defined and T != 1 when c(b-c)(b-c^2-c) != 0
bT, cT := fam_bc(TB);
assert OnlyAllowed(Denominator(bT - bb)) and OnlyAllowed(Denominator(cT - cc));
assert IsDivisibleBy(Pbc!Numerator(bT - bb), Phi12) and IsDivisibleBy(Pbc!Numerator(cT - cc), Phi12);
// so: if Phi12(b,c) = 0 and c (b - c)(b - c^2 - c) != 0, then (b, c) = (b(T), c(T)) with T = T(b,c) != 1
print "B: universality: (0,0) has exact order 12 on E(b,c) iff Phi12(b,c) = 0 and b c (b - c^2 - c) != 0, and then (b,c) = (b(t), c(t)) for t = (d-f)/(1-f) != 1";

// ---------------------------------------------------------------------------
// Part C.  The 2-division polynomial of E_t and the condition for full 2-torsion.
e1 := (6*P)[1];
assert e1 eq t*(3*t^2 - 3*t + 1)/(t - 1);
b2, b4, b6, b8 := Explode(bInvariants(Et));
Rt<W> := PolynomialRing(Qt);
f2 := 4*W^3 + b2*W^2 + 2*b4*W + b6;
qq := f2 div (W - e1);
assert f2 eq qq*(W - e1) and Degree(qq) eq 2 and LeadingCoefficient(qq) eq 4;
qt := 48*t^4 - 96*t^3 + 80*t^2 - 32*t + 4;
assert qt eq 4*(6*t^2 - 6*t + 1)*(2*t^2 - 2*t + 1);
assert qt eq (12*(t - 1/2)^2 - 1)*(4*(t - 1/2)^2 + 1);
ht := (2*t - 1)^3*(2*t^2 - 2*t + 1)/(2*(t - 1)^6);
assert Discriminant(qq) eq qt*ht^2;                         // disc = q(t) h(t)^2
RS<S> := PolynomialRing(Qt);
FF<s> := ext< Qt | S^2 - qt >;                              // the function field of X: s^2 = q(t)
Bt := Coefficient(qq, 1);
e2 := (-Bt + ht*s)/8; e3 := (-Bt - ht*s)/8;
RX<V> := PolynomialRing(FF);
assert RX!f2 eq 4*(V - e1)*(V - e2)*(V - e3);
// every rational function of t used so far has a power of t - 1 as denominator, so the identities of Parts B and C
// specialize to every t != 1 in a field of characteristic 0
function PowerOfTminus1(f)
  den := Denominator(f);
  return den eq (Parent(den).1 - 1)^Degree(den);
end function;
assert forall{ f : f in [bt, ct, Dt, c4t, c6t, e1, Bt, Coefficient(qq, 0), ht] | PowerOfTminus1(f) };
assert Numerator(ht) eq (2*t - 1)^3*(2*t^2 - 2*t + 1)/2;      // h(t) != 0 when E_t is an elliptic curve
print "C: x(6P) = e1 = t(3t^2-3t+1)/(t-1); the other 2-division roots are (-B +- h s)/8 with s^2 = q(t) = 4(6t^2-6t+1)(2t^2-2t+1), h = (2t-1)^3 (2t^2-2t+1)/(2(t-1)^6)";

// ---------------------------------------------------------------------------
// Part D.  The curve X: s^2 = q(t).  Rank 0, X(Q) = Z/4, all four rational points are cusps.
q := 48*x^4 - 96*x^3 + 80*x^2 - 32*x + 4;
C := HyperellipticCurve(q);
EC, phi := EllipticCurve(C, C![0, 2]);
phiinv := Inverse(phi);
assert CremonaReference(EC) eq "24a4";
r, exact := Rank(EC);
assert r eq 0 and exact;
TC, mTC := TorsionSubgroup(EC);
assert Invariants(TC) eq [4];
XQ := [ C![0, 2], C![0, -2], C![1, 2], C![1, -2] ];
assert { phi(pt) : pt in XQ } eq { mTC(g) : g in TC };       // these are all rational points of X
assert not IsSquare(Q!48);                                   // the two points at infinity are not rational
assert phi(C![0, 2]) eq EC!0;
assert cusp_pol(Q!0) eq 0;                                   // t = 0: E_t is singular; t = 1: a pole of b and c
print "D: X = 24a4, rank 0, X(Q) = Z/4 = {(0,+-2), (1,+-2)}, all cusps (t = 0, 1)";

// ---------------------------------------------------------------------------
// Part E.  The four involutions iota_R: P -> R - P of X (R in X(Q)).  For each of them there is lambda in Q(X) with
//   c4(t o iota) = lambda^4 c4(t),  c6(t o iota) = lambda^6 c6(t),  lambda * iota^*(lambda) = 1.
ECF := BaseChange(EC, FF);
gen := ECF![ Evaluate(e, [FF!t, s, 1]) : e in DefiningEquations(phi) ];    // the generic point of X on the model EC
function IsSq(d)   // d in FF; returns true and a square root, or false
  cs := Eltseq(d); A := Qt!cs[1]; B := Qt!cs[2];
  if B eq 0 then
    fl, rr := IsSquare(A); if fl then return true, FF!rr; end if;
    fl, rr := IsSquare(A/qt); if fl then return true, rr*s; end if;
    return false, _;
  end if;
  fl, n := IsSquare(A^2 - qt*B^2);
  if not fl then return false, _; end if;
  for sg in [1, -1] do
    fl, xx := IsSquare((A + sg*n)/2);
    if fl and xx ne 0 then
      yy := B/(2*xx);
      if (xx + yy*s)^2 eq d then return true, xx + yy*s; end if;
    end if;
  end for;
  return false, _;
end function;
function EvalQt(f, a)   // f in Q(t), a in any field extension of Q
  return Evaluate(Numerator(f), a)/Evaluate(Denominator(f), a);
end function;
function EvalFF(f, a, b)   // f in FF evaluated at t = a, s = b
  cs := Eltseq(f);
  return EvalQt(Qt!cs[1], a) + EvalQt(Qt!cs[2], a)*b;
end function;
E1 := BaseChange(Et, FF);
c4a, c6a := Explode(cInvariants(E1));
lambdas := AssociativeArray(); iotas := AssociativeArray();
for g in TC do
  R := mTC(g);
  P2 := ECF!Eltseq(R) - gen;
  co := [ Evaluate(e, Eltseq(P2)) : e in DefiningEquations(phiinv) ];
  t2 := co[1]/co[3]; s2 := co[2]/co[3]^2;                   // (t, s) o iota_R
  assert s2^2 eq EvalQt(qt, t2);
  assert EvalFF(t2, t2, s2) eq FF!t and EvalFF(s2, t2, s2) eq s;         // iota_R is an involution
  b2_, c2_ := fam_bc(t2);
  E2 := EllipticCurve([1 - c2_, -b2_, -b2_, 0, 0]);
  assert jInvariant(E2) eq jInvariant(E1);                               // j o iota_R = j
  c4b, c6b := Explode(cInvariants(E2));
  fl, lam := IsSq((c6b/c4b)/(c6a/c4a));
  assert fl;
  assert c4b eq lam^4*c4a and c6b eq lam^6*c6a;
  assert lam*EvalFF(lam, t2, s2) eq 1;                                   // the cocycle condition
  lambdas[Eltseq(R)] := lam; iotas[Eltseq(R)] := <t2, s2>;
  printf "E: R = %o (order %o): t o iota_R = %o;  lambda_R = %o\n", R, Order(g), t2, lam;
end for;
// the formulas quoted in the proof of Lemma 3: (t o iota_R, lambda_R) for the four points R
w0 := s/(4*(6*t^2 - 6*t + 1));
l0 := FF!((3*t^2 - 3*t + 1)/(2*t^3)); l1 := (5*t^2 - 5*t + 1)*s/(4*t^3*(6*t^2 - 6*t + 1));
assert { <iotas[k][1], lambdas[k]> : k in Keys(iotas) } eq
       { <FF!t, FF!1>, <FF!(1 - t), FF!((t - 1)^3/t^3)>, <1/2 - w0, l0 - l1>, <1/2 + w0, l0 + l1> };
// R = (0, -2) gives the involution (t, s) -> (t, -s), and R = O gives t -> 1/2 - w0
assert iotas[Eltseq(phi(C![0, -2]))] eq <FF!t, -s> and iotas[Eltseq(EC!0)][1] eq 1/2 - w0;
// R = (1, 2) gives (t, s) -> (1 - t, s), and R = (1, -2), the point of order 2, gives t -> 1/2 + w0
assert iotas[Eltseq(phi(C![1, 2]))] eq <FF!(1 - t), s> and iotas[Eltseq(phi(C![1, -2]))][1] eq 1/2 + w0;
assert Order(phi(C![1, -2])) eq 2 and Order(phi(C![1, 2])) eq 4 and Order(phi(C![0, -2])) eq 4;

// ---------------------------------------------------------------------------
// Part F.  X over L = Q(zeta_12): rank 0, exactly 16 points, all cusps.
// rank X(L) = sum of the ranks over Q of the twists of X by 1, -1, -3, 3 (the quadratic subfields of L).
tw := AssociativeArray(); tw[-1] := "48a4"; tw[-3] := "72a1"; tw[3] := "144b1";
for d in [-1, -3, 3] do
  Ed := QuadraticTwist(EC, d);
  assert CremonaReference(Ed) eq tw[d];
  r, exact := Rank(Ed);
  assert r eq 0 and exact;
end for;
// 13 splits completely in L, EC has good reduction at 13, and the torsion of EC(L) injects into EC(F_13)
assert Z!Discriminant(EC) mod 13 ne 0;
assert #ChangeRing(EC, GF(13)) eq 16;
// 37 also splits completely in L, and 13 does not divide #EC(F_37): so EC(L) has no point of order 13, its torsion
// subgroup has order prime to 13 and injects into EC(F_13).  Hence #EC(L) <= 16.
assert Z!Discriminant(EC) mod 37 ne 0 and #ChangeRing(EC, GF(37)) mod 13 ne 0;
L<z> := CyclotomicField(12);
ECL := BaseChange(EC, L);
TL, mTL := TorsionSubgroup(ECL);
ptsL := { mTL(g) : g in TL };
assert #ptsL eq 16;                                          // so X(L) consists of exactly these 16 points
// the points among them with rational coordinates are the four points of Part D (a second proof of X(Q) = Z/4)
assert { ECL!Eltseq(mTC(g)) : g in TC } eq { pt : pt in ptsL | forall{ co : co in Eltseq(pt) | co in Q } };
CL := BaseChange(C, L);
phiinvL := map< ECL -> CL | [ Evaluate(e, [ECL.1, ECL.2, ECL.3]) : e in DefiningEquations(phiinv) ] >;
tvals := {* *};
for pt in ptsL do
  CP := phiinvL(pt);
  if CP[3] eq 0 then Include(~tvals, "oo"); continue; end if;
  tv := CP[1]/CP[3];
  assert tv eq 1 or cusp_pol(tv) eq 0;                       // a cusp
  Include(~tvals, Sprint(MinimalPolynomial(tv)));
end for;
assert Multiplicity(tvals, "oo") eq 2;
printf "F: X(Q(zeta_12)) has rank 0 and 16 points, all cusps; minimal polynomials of t (with multiplicity): %o\n", tvals;

// ---------------------------------------------------------------------------
// Part G.  One quadratic point in each family P + P^sigma = R (the line through -R of slope m; for R = O the
//          vertical line x = m): the curve E_t has torsion Z/2 x Z/12 over K and descends to Q as in Lemma 3.
a1, a2, a3, a4, a6 := Explode(aInvariants(EC));
function QuadPoint(R, m)
  if R eq EC!0 then
    pol := x^2 + (a1*m + a3)*x - (m^3 + a2*m^2 + a4*m + a6);          // in y
    assert IsIrreducible(pol);
    K<w> := NumberField(pol);
    return K, [K | m, w];
  end if;
  mR := -R; x0 := mR[1]; y0 := mR[2];
  ln := y0 + m*(x - x0);
  gg := ln^2 + a1*x*ln + a3*ln - (x^3 + a2*x^2 + a4*x + a6);
  pol := -gg div (x - x0);
  assert Degree(pol) eq 2 and IsIrreducible(pol);
  K<w> := NumberField(pol);
  return K, [K | w, y0 + m*(w - x0)];
end function;
for g in TC do
  R := mTC(g);
  K, xy := QuadPoint(R, 3);
  sig := Automorphisms(K)[2];
  assert sig(K.1) ne K.1;
  ECK := BaseChange(EC, K); CK := BaseChange(C, K);
  PK := ECK!xy; PKs := ECK![sig(xy[1]), sig(xy[2])];
  assert PK + PKs eq ECK!Eltseq(R);                                    // P + P^sigma = R
  phiinvK := map< ECK -> CK | [ Evaluate(e, [ECK.1, ECK.2, ECK.3]) : e in DefiningEquations(phiinv) ] >;
  CP := phiinvK(PK);
  assert CP[3] ne 0;
  t0 := CP[1]/CP[3]; s0 := CP[2]/CP[3]^2;
  assert s0^2 eq Evaluate(q, t0);
  assert t0 ne 1 and cusp_pol(t0) ne 0;                                // not a cusp
  io := iotas[Eltseq(R)];
  assert EvalFF(io[1], t0, s0) eq sig(t0) and EvalFF(io[2], t0, s0) eq sig(s0);   // P^sigma = iota_R(P)
  b0, c0 := fam_bc(t0);
  EK := EllipticCurve([1 - c0, -b0, -b0, 0, 0]);
  assert Invariants(TorsionSubgroup(EK)) eq [2, 12];
  u := EvalFF(lambdas[Eltseq(R)], t0, s0);
  assert u*sig(u) eq 1 and u ne -1;
  ww := 1 + u;                                                         // u = ww/sig(ww)
  assert u eq ww/sig(ww);
  c4K, c6K := Explode(cInvariants(EK));
  A4 := c4K*ww^4; A6 := c6K*ww^6;
  assert A4 in Q and A6 in Q;
  E0 := EllipticCurve([-27*Q!A4, -54*Q!A6]);
  assert IsIsomorphic(BaseChange(E0, K), EK);
  printf "G: R = %o, K = Q(sqrt(%o)): t = %o, torsion of E_t(K) = [2,12], E_t is the base change of the curve over Q with j = %o\n", R, SquarefreeFactorization(Z!Discriminant(MaximalOrder(K))), t0, jInvariant(E0);
end for;
// A test that does not use the involutions iota_R: K = Q(sqrt(13)), over which X has positive rank.  For the point
// P0 with t = 2 sqrt(13) - 7 and the points n P0 + T (1 <= n <= 3, T in X(Q)): t is not rational, the curve E_t has
// torsion Z/2 x Z/12 over K, j(E_t) is rational, and E_t is the base change of an elliptic curve over Q.
K13<w13> := QuadraticField(13);
sig13 := hom< K13 -> K13 | -w13 >;
EC13 := BaseChange(EC, K13); C13 := BaseChange(C, K13);
phi13 := map< C13 -> EC13 | [ Evaluate(e, [C13.1, C13.2, C13.3]) : e in DefiningEquations(phi) ] >;
phiinv13 := map< EC13 -> C13 | [ Evaluate(e, [EC13.1, EC13.2, EC13.3]) : e in DefiningEquations(phiinv) ] >;
t13 := 2*w13 - 7;
fl, s13 := IsSquare(Evaluate(q, t13));
assert fl;
P0 := phi13(C13![t13, s13]);
assert Order(P0) eq 0;                                                 // a point of infinite order
n13 := 0;
for n in [1..3] do for g in TC do
  PT := n*P0 + EC13!Eltseq(mTC(g));
  CP := phiinv13(PT);
  assert CP[3] ne 0;
  t0 := CP[1]/CP[3];
  assert not t0 in Q and t0 ne 1 and cusp_pol(t0) ne 0;
  b0, c0 := fam_bc(t0);
  EK := EllipticCurve([1 - c0, -b0, -b0, 0, 0]);
  assert Invariants(TorsionSubgroup(EK)) eq [2, 12];
  jK := jInvariant(EK);
  assert jK in Q and jK ne 0 and jK ne 1728;
  E1 := EllipticCurveFromjInvariant(Q!jK);
  c4K, c6K := Explode(cInvariants(EK)); d4, d6 := Explode(cInvariants(E1));
  delta := (c6K/d6)/(c4K/d4);                                          // E_t is the twist of E1 by delta
  fl, m := IsSquare(Q!Norm(delta));
  assert fl;                                                           // so delta is rational up to a square of K
  if delta/m eq -1 then m := -m; end if;
  rr := m/(Q!Norm(1 + delta/m));                                       // delta = rr * ((1 + delta/m))^2
  assert delta eq rr*(1 + delta/m)^2;
  assert IsIsomorphic(BaseChange(QuadraticTwist(E1, rr), K13), EK);
  n13 +:= 1;
end for; end for;
assert n13 eq 12;
print "G: K = Q(sqrt(13)): the 12 points n P0 + T (P0 with t = 2 sqrt(13) - 7, 1 <= n <= 3, T in X(Q)) give curves E_t with torsion [2,12], rational j-invariant, and E_t is a base change from Q";
// The contrast of Remark 5 (not used in the proof): the analogue of Lemma 3 for Z/2 x Z/10 is false.  Over
// K = Q(sqrt(17)) the curve E(b,c) with b = (5 - 3 sqrt(17))/72, c = (sqrt(17) - 7)/24 (the parameter
// f = (9 + sqrt(17))/16 of Kubert's family for Z/10: d = f^2/(f - (f-1)^2), c = f(d-1), b = cd) has torsion
// Z/2 x Z/10, its j-invariant is rational and it is isomorphic over K to its conjugate, but it is not a base
// change from Q.  Indeed, if E is isomorphic over K to (E_0)_K with E_0 over Q, then c4(E) = v^4 c4(E_0) and
// c6(E) = v^6 c6(E_0) for some v in K^*, so c6(E)/c4(E) is v^2 times a rational number and its norm is a square.
K17<w17> := QuadraticField(17);
sig17 := hom< K17 -> K17 | -w17 >;
f10 := (9 + w17)/16;
d10 := f10^2/(f10 - (f10 - 1)^2); c10 := f10*(d10 - 1); b10 := c10*d10;
assert b10 eq (5 - 3*w17)/72 and c10 eq (w17 - 7)/24;
E10 := EllipticCurve([1 - c10, -b10, -b10, 0, 0]);
assert Order(E10![0, 0]) eq 10;
assert Invariants(TorsionSubgroup(E10)) eq [2, 10];
j10 := jInvariant(E10);
assert j10 eq 141420761/9216;
assert IsIsomorphic(EllipticCurve([ sig17(a) : a in aInvariants(E10) ]), E10);
c4w, c6w := Explode(cInvariants(E10));
assert c4w ne 0 and c6w ne 0;
n10 := Q!Norm(c6w/c4w);
assert not IsSquare(n10) and IsSquare(17*n10);
print "G: contrast, K = Q(sqrt(17)): E(b,c) with b = (5 - 3 sqrt(17))/72, c = (sqrt(17) - 7)/24 has torsion [2,10], j = 141420761/9216 and is isomorphic to its conjugate, but the norm of c6/c4 is 17 times a square, so it is not a base change from Q";

// ---------------------------------------------------------------------------
// Part H.  Congruences for the Frobenius polynomials of Sections 6 and 7 (Lemma 7, Propositions 8 and 9).
// Rotations: the three polynomials of [P612, Lemma 2] (elements of order n = 3, 4, 6), with integers a, u, b.
PH<X, av, pv> := PolynomialRing(Z, 3);
rotpol := AssociativeArray();
rotpol[3] := X^4 + av*X^3 + (av^2 - pv)*X^2 + av*pv*X + pv^2;
rotpol[4] := X^4 - 2*av*X^3 + 2*av^2*X^2 - 2*pv*av*X + pv^2;
rotpol[6] := X^4 + 3*av*X^3 + (3*av^2 + pv)*X^2 + 3*av*pv*X + pv^2;
// their values at X = 1 are the norm forms displayed in [P612, Lemma 2]
assert Evaluate(rotpol[3], [1, av, pv]) eq (1 + av)^2 - (1 + av)*(av + pv) + (av + pv)^2;
assert Evaluate(rotpol[4], [1, av, pv]) eq (1 - av)^2 + (pv - av)^2;
assert Evaluate(rotpol[6], [1, av, pv]) eq (1 + 2*av + pv)^2 - (1 + 2*av + pv)*(av + pv) + (av + pv)^2;
for n in [3, 4, 6] do
  for a in [0..31] do for p in [1..31 by 2] do
    if Evaluate(rotpol[n], [1, a, p]) mod 32 eq 0 then assert p mod 8 eq 1; end if;
  end for; end for;
end for;
print "H: rotations of order 3, 4, 6: if 32 divides P(1) and p is odd, then p = 1 mod 8";
// Involutions (Lemma 7): with beta1 beta2 = delta, rho^2 = r^2 c and r^2 c delta = s p,
//   (X^2 - rho^2 beta1^2)(X^2 - rho^2 beta2^2) = (X^2 + s p)^2 - c (r (beta1 + beta2))^2 X^2.
PI<Xi, r, c, be1, be2> := PolynomialRing(Q, 5);
lhs := (Xi^2 - r^2*c*be1^2)*(Xi^2 - r^2*c*be2^2);
sp := r^2*c*be1*be2;
assert lhs eq (Xi^2 + sp)^2 - c*(r*(be1 + be2))^2*Xi^2;
// c squarefree (so c is not 0 mod 4) and not 1 mod 8:  32 | x^2 - c y^2  =>  4 | x.
for c0 in [0..31] do
  if c0 mod 4 eq 0 or c0 mod 8 eq 1 then continue; end if;
  for xx in [0..31] do for yy in [0..31] do
    if (xx^2 - c0*yy^2) mod 32 eq 0 then assert xx mod 4 eq 0; end if;
  end for; end for;
end for;
// c not 1 mod 3:  3 | x^2 - c y^2  =>  3 | x.
for c0 in [0, 2] do for xx in [0..2] do for yy in [0..2] do
  if (xx^2 - c0*yy^2) mod 3 eq 0 then assert xx mod 3 eq 0; end if;
end for; end for; end for;
print "H: involutions: if c is not a square in Q_2 then 32 | x^2 - c y^2 forces 4 | x; if c is not a square in Q_3 then 3 | x^2 - c y^2 forces 3 | x";

// ---------------------------------------------------------------------------
// Part S.  The census curve 2.2.2.6.c (Remark 2; not used in the proof): its Jacobian has rational torsion
// subgroup [2,2,2,6], and it is simple over Q because its Frobenius polynomial at 29 is irreducible.
fS := 6993*x^5 - 25655829*x^4 + 14380189111*x^3 + 16095094765893*x^2 + 2757085929534360*x + 1355999975081200;
CS := HyperellipticCurve(fS, x^2 + x);
assert Z!Discriminant(CS) mod 29 ne 0;
chi29 := Reverse(LPolynomial(ChangeRing(CS, GF(29))));
assert chi29 eq x^4 - 26*x^2 + 841 and IsIrreducible(chi29);
ordsS := [ Z!Evaluate(LPolynomial(ChangeRing(CS, GF(pp))), 1) : pp in [13, 19, 29, 41, 43] ];
assert forall{ pp : pp in [13, 19, 29, 41, 43] | Z!Discriminant(CS) mod pp ne 0 };
assert GCD(ordsS) eq 48;                                               // so #J(Q)_tors divides 48
JS := Jacobian(SimplifiedModel(CS));
TS := TorsionSubgroup(JS);
assert Invariants(TS) eq [2, 2, 2, 6];
print "S: census curve 2.2.2.6.c: rational torsion subgroup [2,2,2,6]; Frobenius polynomial at 29 is x^4 - 26x^2 + 841, irreducible";

print "ALL_CHECKS_PASSED";
quit;
