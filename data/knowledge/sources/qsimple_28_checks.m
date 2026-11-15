// qsimple_28_checks.m  (Magma V2.30-1)
// Certificate for the census proof qsimple_28_obstruction.tex: an abelian surface A over Q with a principal
// polarization defined over Q, simple over Q and geometrically split, with a rational point of order 7, has
// #A(Q)[2^oo] <= 2; hence [28] and [2,14] do not occur for Jacobians in this class.
// Every check is an assertion; a failed assertion stops the run before the final marker.
// Parts:  A  supersingular Weil polynomials (Lemma 1);  B  X1(14) = 14a4, moduli map, diamonds, cusps (Lemma 4);
//         C  the six families of quadratic points, norms of the 2-division invariants, exceptional values;
//         D  rational points of the condition curves (Propositions 13 and 14);
//         E  3-adic input of Proposition 11;
//         V  validation, not used in the proof (Remarks 1-3).
// No call to SetClassGroupBounds: all class group computations are unconditional.
// Final marker: ALL_28_CHECKS_PASSED.
SetColumns(0); SetMemoryLimit(16*10^9);
Q := Rationals(); Z := Integers();
QT<T> := PolynomialRing(Q);
Qm<m> := PolynomialRing(Q);

//////////////////////////////////////////////////////////////////////////////
// Part A.  Weil polynomials of supersingular abelian surfaces over F_p, p in {5,7,11,13,17,19,23}:
// P(T) = T^4 + a1 T^3 + a2 T^2 + p a1 T + p^2, real polynomial h(y) = y^2 + a1 y + a2 - 2p with roots in
// [-2 sqrt p, 2 sqrt p].  Supersingular <=> all roots of P are sqrt(p) times roots of unity.  We list them,
// compare with the classification of Lemma 1 (h in {y^2 - k p : 0 <= k <= 4}, plus y^2 +- 5y + 5 for p = 5),
// and check that 28 | P(1) forces p = +-1 mod 7 (p = 7 included).
printf "Part A\n";
function SuperSingularList(p)
  L := [];
  B := Isqrt(16*p);
  for a1 in [-B..B] do
    for a2 in [-6*p..6*p] do
      // real Weil polynomial h(y) = y^2 + a1 y + (a2 - 2p): real roots in [-2 sqrt p, 2 sqrt p]
      c := a2 - 2*p;
      disc := a1^2 - 4*c;
      if disc lt 0 then continue; end if;
      // roots (-a1 +- sqrt disc)/2 within [-2sqrt p, 2 sqrt p]: |root| <= 2 sqrt p  <=> h(+-2sqrt p) >= 0 and |a1| <= 4 sqrt p
      // roots in [-2 sqrt p, 2 sqrt p]  <=>  a1^2 <= 16p,  4p + c >= 0  and  (4p + c)^2 >= 4 a1^2 p
      if a1^2 gt 16*p or 4*p + c lt 0 or (4*p + c)^2 lt 4*a1^2*p then continue; end if;
      P := T^4 + a1*T^3 + a2*T^2 + p*a1*T + p^2;
      // supersingular: p | a1 and p | a2 (slopes 1/2 for abelian surfaces: Newton polygon of P)
      if a1 mod p eq 0 and a2 mod p eq 0 then
        // and P is a Weil polynomial: roots of absolute value sqrt p (guaranteed by the real-root condition)
        Append(~L, <a1, a2, Z!Evaluate(P, 1)>);
      end if;
    end for;
  end for;
  return L;
end function;
for p in [5, 7, 11, 13, 17, 19, 23] do
  L := SuperSingularList(p);
  hs := [ <x[1], x[2] - 2*p> : x in L ];  // (a1, c) with h = y^2 + a1 y + c
  expected := [ <0, -k*p> : k in [0..4] ] cat (p eq 5 select [<5, 5>, <-5, 5>] else []);
  // every listed h must be expected (the list may contain non-isogeny-class-realizable ones: Lemma 1 only needs "contained in")
  assert forall{h : h in hs | h in expected};
  vals := {x[3] : x in L};
  assert vals subset {(p + 1)^2 - k*p : k in [0..4]} join (p eq 5 select {11, 71} else {});
  if p eq 5 then assert {11, 71} subset vals; end if;
  if p eq 7 then assert vals eq {36, 43, 50, 57, 64}; end if;
  printf "A p = %o (p mod 7 = %o): supersingular real Weil polynomials (a1, a2 - 2p) = %o, values P(1) = %o\n", p, p mod 7, hs, vals;
  if p mod 7 notin {1, 6} then
    assert forall{v : v in vals | v mod 28 ne 0};
  end if;
end for;
// for Q(sqrt(-7)): 5 is inert
assert KroneckerSymbol(-7, 5) eq -1;
printf "A_DONE\n";

//////////////////////////////////////////////////////////////////////////////
// Part B.  X1(14) = 14a4 : Y^2 + XY + Y = X^3 - X.
printf "Part B\n";
EC := EllipticCurve([1, 0, 1, -1, 0]);
assert CremonaReference(CremonaDatabase(), EC) eq "14a4";
S2 := TwoSelmerGroup(EC);
TC, mTC := TorsionSubgroup(EC);
assert Invariants(TC) eq [6];
assert #S2 eq 2;                         // 2-Selmer group of order 2 = E(Q)[2]: rank 0
lo, hi := RankBounds(EC); assert hi eq 0;
ratpts := {mTC(g) : g in TC};
assert ratpts eq {EC!0, EC![0,0], EC![0,-1], EC![1,0], EC![1,-2], EC![-1,0]};
assert Order(EC![0,0]) eq 3 and Order(EC![-1,0]) eq 2 and Order(EC![1,0]) eq 6;
FE<X, Y> := FunctionField(EC);
function KubertData(t)
  b := t^3 - t^2; c := t^2 - t;
  a1 := 1 - c; a2 := -b; a3 := -b; a4 := 0*t; a6 := 0*t;
  b2 := a1^2 + 4*a2; b4 := 2*a4 + a1*a3; b6 := a3^2 + 4*a6; b8 := a1^2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3^2 - a4^2;
  c4 := b2^2 - 24*b4; c6 := -b2^3 + 36*b2*b4 - 216*b6;
  Dl := -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6;
  return b, c, b2, b4, b6, c4, c6, Dl;
end function;
tt := 1 - X*(X - 1)^2/Y^2;
sg := -(X - 1)^2/(2*Y);
b, c, b2, b4, b6, c4, c6, Dl := KubertData(tt);
ee := b/(1 - c - 2*tt*sg);
// e is a root of the 2-division polynomial 4x^3 + b2 x^2 + 2 b4 x + b6 of E_7(t)
assert 4*ee^3 + b2*ee^2 + 2*b4*ee + b6 eq 0;
jj := c4^3/Dl;
// the map t : X1(14) -> X1(7) (the t-line) has degree 3, j(E_7(t)) has degree 24 in t, and j has degree
// 3 * 24 = 72 on X1(14)
assert Degree(Denominator(Divisor(tt))) eq 3;
Qtl<tl> := FunctionField(Q);
_, _, _, _, _, c4l, _, Dll := KubertData(tl);
assert Max(Degree(Numerator(c4l^3/Dll)), Degree(Denominator(c4l^3/Dll))) eq 24;
assert Degree(Denominator(Divisor(jj))) eq 72;
// the plane model of the pairs (t, e): G(t, e) = 4e^3 + b2 e^2 + 2 b4 e + b6 = 0; its singular points are degenerate
P2<t0, e0> := PolynomialRing(Q, 2);
bb, cc, bb2, bb4, bb6, cc4, cc6, DDl := KubertData(t0);
G := 4*e0^3 + bb2*e0^2 + 2*bb4*e0 + bb6;
// G is irreducible over Q; as deg t = 3 = deg_e G and Q(X) = Q(t, e) (below), G(t, e)/4 is the minimal
// polynomial of e over Q(t) (also over Qbar(t)), so G = 0 is geometrically irreducible
faG := Factorization(G);
assert #faG eq 1 and faG[1][2] eq 1 and Degree(G, e0) eq 3;
A2 := AffineSpace(Q, 2);
Cte := Curve(A2, G);
assert IsAbsolutelyIrreducible(Cte);
// the singular subscheme: show t is degenerate there (b = 0 or Delta = 0)
Isg := ideal<P2 | G, Derivative(G, t0), Derivative(G, e0)>;
assert bb*DDl in Radical(Isg);
printf "B singular locus of the (t,e)-model lies over b*Delta_7(t) = 0\n";
// (t, e) <- (X, Y) is birational: t and e generate the function field (degree of Q(t,e) in FE is 1)
// check: X and Y are rational functions of t and sg, with sg = (1 - c - b/e)/(2t), via the explicit formulas
// X = Xb(u, sg), Y = Yb(u, sg) with u = t:
u := tt;
Xb := (u^3 - 4*u*sg^2 - u^2 + 2*u*sg + 4*sg^2 - u - 2*sg + 1)/(-u^3 + 3*u^2 - 3*u + 1);
Yb := (2*u^2*sg - 8*sg^3 + 4*sg^2 - 2*sg)/(-u^3 + 3*u^2 - 3*u + 1);
assert Xb eq X and Yb eq Y;
assert sg eq (1 - c - b/ee)/(2*tt);
printf "B (X, Y) are rational functions of (t, e): the map X1(14) -> {G = 0} is birational\n";
// cusps: the six rational points are poles of j
for C in ratpts do assert Valuation(jj, Place(C)) lt 0; end for;
// widths: O, (0,0), (0,-1) have width 7; (1,0), (1,-2), (-1,0) have width 14; w(C) = w(-C)
wid := func<C | -Valuation(jj, Place(C))>;
assert [wid(C) : C in [EC!0, EC![0,0], EC![0,-1]]] eq [7,7,7];
assert [wid(C) : C in [EC![1,0], EC![1,-2], EC![-1,0]]] eq [14,14,14];
for C in ratpts do assert wid(C) eq wid(-C); end for;
// non-cuspidal points with Y = 0 do not exist: Y = 0 forces X in {0, 1, -1}
assert Roots(QT!(T^3 - T)) eq [<-1,1>, <0,1>, <1,1>];
// diamond operators: translation by the 3-torsion point (0,0) preserves the K-isomorphism class of E_7(t)
P := EC(FE)![X, Y];
Pt := P + EC(FE)![0, 0];
ttau := 1 - Pt[1]*(Pt[1] - 1)^2/Pt[2]^2;
_, _, _, _, _, c4t, c6t, Dlt := KubertData(ttau);
assert c4t^3/Dlt eq jj;
assert IsSquare(c4t*c6/(c4*c6t));
// [-1] is not a diamond operator
Pm := -P;
tm := 1 - Pm[1]*(Pm[1] - 1)^2/Pm[2]^2;
_, _, _, _, _, c4m, c6m, Dlm := KubertData(tm);
assert c4m^3/Dlm ne jj;
// no point of X1(14) of degree <= 2 has j in {0, 1728}: the zeros of j and of j - 1728 are places of degree
// 6, 18 and 12, 24
for val in [0, 1728] do
  Zv := Numerator(Divisor(jj - val));
  dg := {Degree(pr[1]) : pr in Decomposition(Zv)};
  assert dg eq (val eq 0 select {6, 18} else {12, 24});
  printf "B zeros of j - %o: places of degrees %o\n", val, dg;
end for;
printf "B_DONE\n";

//////////////////////////////////////////////////////////////////////////////
// Part C.  Families of quadratic points.  For R in X1(14)(Q):
//   R = O:     Q = (m, y), y^2 + (m+1) y - (m^3 - m) = 0          (Q^sigma = -Q)
//   R != O:   line through -R with slope m; x(Q) root of q_{R,m}    (Q + Q^sigma = R)
// Dp = discriminant of the quadratic factor of the 2-division polynomial (up to 16), bp = f'(e).
printf "Part C\n";
F<mf> := FunctionField(Q);
PF<XX> := PolynomialRing(F);
function FamilyField(R)
  if R eq EC!0 then
    q := XX^2 + (mf + 1)*XX - (mf^3 - mf);
    K<w> := ext<F | q>;
    return q, K, K!mf, w;
  end if;
  P0 := -R; xr := P0[1]; yr := P0[2];
  cub := (yr + mf*(XX - xr))^2 + XX*(yr + mf*(XX - xr)) + (yr + mf*(XX - xr)) - XX^3 + XX;
  qq, rr := Quotrem(cub, XX - xr);
  assert rr eq 0;
  qq := qq/LeadingCoefficient(qq);
  K<w> := ext<F | qq>;
  return qq, K, w, yr + mf*(w - xr);
end function;
function TwoDivisionInvariants(x0, y0)
  t := 1 - x0*(x0 - 1)^2/y0^2;
  s := -(x0 - 1)^2/(2*y0);
  b := t^3 - t^2; c := t^2 - t;
  a1 := 1 - c; a2 := -b; a3 := -b;
  b2 := a1^2 + 4*a2; b4 := a1*a3; b6 := a3^2;
  e := b/(1 - c - 2*t*s);
  be := 4*e + b2; ga := be*e + 2*b4;
  assert ga*e + b6 eq 0;
  Dp := be^2 - 16*ga;
  bp := 4*e^2 + be*e + ga;
  return Dp, bp, t, s, e;
end function;
function SqClassQm(r)
  num := Qm!Numerator(r); den := Qm!Denominator(r);
  fa := Factorization(num*den);
  g := &*[Qm | f[1] : f in fa | IsOdd(f[2])];
  lc := LeadingCoefficient(num*den);
  cst := SquarefreeFactorization(Z!(Numerator(lc)*Denominator(lc)));
  return cst*g;
end function;
function FacWithUnit(g)
  // Factorization drops the leading unit; recover it and print unit * (factor)^e * ...
  fac := Factorization(g);
  unit := g div &*[Parent(g) | t[1]^t[2] : t in fac];
  assert unit*&*[Parent(g) | t[1]^t[2] : t in fac] eq g;
  return Sprint(unit) cat &cat[" * (" cat Sprint(t[1]) cat ")" cat (t[2] gt 1 select "^" cat Sprint(t[2]) else "") : t in fac];
end function;
function RationalRootsOf(r)
  num := Qm!Numerator(r); den := Qm!Denominator(r);
  return {x[1] : x in Roots(num) cat Roots(den)};
end function;
expD := AssociativeArray(); expB := AssociativeArray();
expD[EC!0]       := -(m^3 - 2*m^2 - m + 1);          expB[EC!0]       := (m + 1)*(m^3 - 9*m^2 - m + 1);
expD[EC![0,-1]]  := -(m + 1)*(m^3 + 2*m^2 - m - 1);  expB[EC![0,-1]]  := (m + 2)*(m^3 + 2*m^2 - 8*m - 8);
expD[EC![0,0]]   := -m*(m^3 + m^2 - 2*m - 1);        expB[EC![0,0]]   := (m - 1)*(m^3 + m^2 - 9*m - 1);
expD[EC![1,-2]]  := m*(m - 1)*(m^3 + 3*m^2 - 4*m + 1);   expB[EC![1,-2]] := -m*(m - 1)*(m^3 - 4*m^2 + 3*m + 1);
expD[EC![1,0]]   := (m + 1)*(m + 2)*(m^3 + 7*m^2 + 14*m + 7); expB[EC![1,0]] := -(m + 1)*(m + 2)*(m^3 - 7*m - 7);
expD[EC![-1,0]]  := m*(m + 1)*(m^3 + 5*m^2 + 6*m + 1);   expB[EC![-1,0]] := -m*(m + 1)*(m^3 - 2*m^2 - m + 1);
// d = disc q_{R,m} for the two families used in the text: A (R = O) and B (R = (-1,0))
expd := AssociativeArray();
expd[EC!0] := (m + 1)*(4*m^2 - 3*m + 1);  expd[EC![-1,0]] := m^4 + 2*m^3 + 7*m^2 + 6*m + 1;
excep := AssociativeArray();
for R in ratpts do
  q, K, xq, yq := FamilyField(R);
  assert yq^2 + xq*yq + yq eq xq^3 - xq;
  Dp, bp, t, s, e := TwoDivisionInvariants(xq, yq);
  ND := Norm(Dp); Nb := Norm(bp); d := Discriminant(q);
  if IsDefined(expd, R) then assert d eq Evaluate(expd[R], mf); end if;
  assert SqClassQm(ND) eq expD[R];
  assert SqClassQm(Nb) eq expB[R];
  // exceptional parameter values: rational roots of numerators/denominators of d, N(Dp), N(bp), N(t), N(1-t),
  // N(e), N(Delta), and of the coordinates of t, s, e, Dp, bp in the basis 1, w
  ex := RationalRootsOf(d) join RationalRootsOf(ND) join RationalRootsOf(Nb) join RationalRootsOf(Norm(t)) join
        RationalRootsOf(Norm(1 - t)) join RationalRootsOf(Norm(e));
  _, _, _, _, _, _, _, Dl7 := KubertData(t);
  ex join:= RationalRootsOf(Norm(Dl7));
  for el in [t, s, e, Dp, bp] do
    for co in Eltseq(el) do if co ne 0 then ex join:= RationalRootsOf(co); end if; end for;
  end for;
  excep[R] := ex;
  assert Denominator(d) eq 1;
  printf "C family R = %o: q = %o, d = %o = %o; square classes N(Dp) ~ %o, N(bp) ~ %o; exceptional m: %o\n", R, q, Qm!Numerator(d), FacWithUnit(Qm!Numerator(d)), expD[R], expB[R], ex;
end for;
// the exceptional values of the two families used in the text
assert excep[EC!0] eq {Q| -1, 0, 1/2, 1} and excep[EC![-1,0]] eq {Q| -1, -1/2, 0};
// Exceptional values: either q_{R,m0} is reducible over Q (no quadratic point), or m0 gives a degenerate curve,
// or we check the conditions directly on the specialized quadratic point.
function SpecializedPoint(R, m0)
  // returns ok, K, x0, y0 for the quadratic point of family R at m0 (ok = false if q reducible)
  if R eq EC!0 then
    q0 := QT!(T^2 + (m0 + 1)*T - (m0^3 - m0));
    if not IsIrreducible(q0) then return false, 0, 0, 0; end if;
    K0<w0> := NumberField(q0);
    return true, K0, K0!m0, w0;
  end if;
  P0 := -R; xr := P0[1]; yr := P0[2];
  cub := (yr + m0*(T - xr))^2 + T*(yr + m0*(T - xr)) + (yr + m0*(T - xr)) - T^3 + T;
  q0, r0 := Quotrem(cub, T - xr);
  assert r0 eq 0;
  q0 := q0/LeadingCoefficient(q0);
  if not IsIrreducible(q0) then return false, 0, 0, 0; end if;
  K0<w0> := NumberField(q0);
  return true, K0, w0, yr + m0*(w0 - xr);
end function;
function SqClassQ(x)
  x := Q!x; return SquarefreeFactorization(Z!(Numerator(x)*Denominator(x)));
end function;
genuine := {};
for R in ratpts do
  for m0 in excep[R] do
    ok, K0, x0, y0 := SpecializedPoint(R, m0);
    if not ok then printf "C exceptional R=%o m=%o: q reducible\n", R, m0; continue; end if;
    // nondegenerate quadratic point: test directly
    t0 := 1 - x0*(x0 - 1)^2/y0^2;
    b0 := t0^3 - t0^2; c0 := t0^2 - t0;
    E0 := EllipticCurve([1 - c0, -b0, -b0, 0, 0]);
    f0 := Polynomial([b0^2, 2*(-b0*(1 - c0)), (1 - c0)^2 - 4*b0, 4]);
    rts := Roots(f0);
    assert #rts eq 1;
    e0 := rts[1][1];
    qd := f0 div (Parent(f0).1 - e0);
    Dp0 := Discriminant(qd); bp0 := Evaluate(qd, e0);
    d0 := SqClassQ(Discriminant(DefiningPolynomial(K0)));
    nd := SqClassQ(Norm(Dp0)); nb := SqClassQ(Norm(bp0));
    printf "C exceptional R=%o m=%o: quadratic point over Q(sqrt %o), class N(Dp) = %o, N(bp) = %o\n", R, m0, d0, nd, nb;
    // (N1) must fail, or (family A) (N2) must fail
    assert nd ne 1 or (R in {EC!0, EC![0,0], EC![0,-1]} and nb notin {1, d0});
    Include(~genuine, <R[1], R[2], R[3], m0, d0, nd, nb>);
  end for;
end for;
// exactly three genuine quadratic points occur at exceptional values, and (N1) fails at all of them
assert genuine eq {<Q!0, Q!1, Q!0, Q!1/2, 3, -2, -39>, <Q!-1, Q!0, Q!1, Q!-1/2, -7, 14, 14>,
                   <Q!1, Q!0, Q!1, Q!-3, 17, 2, 26>};
printf "C genuine points at exceptional values (R, m, d, N(Dp), N(bp) classes): %o\n", genuine;
printf "C_DONE\n";

//////////////////////////////////////////////////////////////////////////////
// Part D.  Rational points of the condition curves.
printf "Part D\n";
function Genus2AllPoints(f)
  H := HyperellipticCurve(f);
  assert Genus(H) eq 2;
  J := Jacobian(H);
  lo, hi := RankBounds(J);
  if hi eq 0 then
    return Chabauty0(J), 0;
  end if;
  assert lo eq 1 and hi eq 1;
  pts := RationalPoints(H : Bound := 1000);
  gens := [];
  for p1 in pts do for p2 in pts do
    try
      Pj := J![p1, p2];
      if Order(Pj) eq 0 then Append(~gens, Pj); end if;
    catch err ; end try;
  end for; end for;
  assert #gens gt 0;
  return Chabauty(gens[1]), 1;
end function;
// D1: family B (R in the coset of (-1,0)): y^2 = class of N(Dp) has genus 2 and rank 0
expms := AssociativeArray();
expms[EC![1,-2]] := {Q| 0, 1};  expms[EC![1,0]] := {Q| -2, -1};  expms[EC![-1,0]] := {Q| -1, 0};
for R in [EC![1,-2], EC![1,0], EC![-1,0]] do
  pts, r := Genus2AllPoints(expD[R]);
  assert r eq 0;
  ms := {p[1]/p[3] : p in pts | p[3] ne 0};
  assert ms eq expms[R] and #pts eq 3 and #{p : p in pts | p[3] eq 0} eq 1;
  printf "D family B, R = %o: (N1) curve y^2 = %o, rank 0, finite m: %o\n", R, expD[R], ms;
  for m0 in ms do
    ok := SpecializedPoint(R, m0);
    assert not ok;          // q reducible: no quadratic point
  end for;
end for;
// D2: family A (R = O).  (N1) & (N2a):  gD(m) and (m+1) c(m) both in Q*^2, c(m) = m^3 - 9m^2 - m + 1.
// Since (m+1)^4 gD(tau m) = (m+1) c(m) with tau(m) = (1-m)/(1+m), this says: m and tau(m) both give points of
// C_D : Y^2 = H(x) = x^3 + 2x^2 - x - 1 (x = -m), an elliptic curve of rank 1 with trivial torsion.
// The 2-descent map C_D(Q)/2C_D(Q) -> L*/L*^2 (L = Q(alpha), H(alpha) = 0), x - alpha, has image {1, kappa0}.
// Hence (m + alpha)((alpha - 1) m + alpha + 1)(m + 1) lies in delta L*^2, delta in {1, kappa0}, and m is the
// image of an L-point of E_delta : delta v^2 = (alpha - 1)(m + alpha)(m + (alpha+1)/(alpha-1))(m + 1).
cb := m^3 - 9*m^2 - m + 1;
gD := expD[EC!0];
FQ<mq> := FunctionField(Q);
assert (mq + 1)^4*Evaluate(gD, (1 - mq)/(1 + mq)) eq (mq + 1)*Evaluate(cb, mq);
Hx := T^3 + 2*T^2 - T - 1;
assert Evaluate(Hx, -m) eq gD;
CD := EllipticCurve(Hx);
lo, hi := RankBounds(CD); assert lo eq 1 and hi eq 1;
assert #TorsionSubgroup(CD) eq 1;
P0 := CD![-1, -1];
assert Order(P0) eq 0;
L<al> := NumberField(Hx);
assert IsAbelian(L) and Discriminant(Integers(L)) eq 49;     // L = Q(zeta_7)^+
kappa0 := L!(P0[1]/P0[3]) - al;
assert not IsSquare(kappa0);               // P0 is not in 2 C_D(Q), so the image of C_D(Q) is {1, kappa0}
// the product of the descent images of (-m, .) and (-tau(m), .), times (m + 1)^2
FL<mL> := FunctionField(L);
assert (-mL - al)*(-(1 - mL)/(1 + mL) - al)*(mL + 1)^2 eq (mL + al)*((al - 1)*mL + al + 1)*(mL + 1);
assert -(al + 1)/(al - 1) eq 2*al^2 + 6*al + 3;
PL<XL> := PolynomialRing(L);
P1Q := ProjectiveSpace(Q, 1);
aL := al - 1;
rL := [ -al, -(al + 1)/(al - 1), L!(-1) ];
mrat := {Q|}; ninf := 0;
for delta in [L!1, kappa0] do
  Ed := EllipticCurve(&*[XL - delta*aL*ri : ri in rL]);    // V^2 = prod (M - delta a r_i), m = M/(delta a)
  Ecov := map<Ed -> P1Q | [Ed.1, delta*aL*Ed.3]>;
  lo, hi := RankBounds(Ed);
  Td, mTd := TorsionSubgroup(Ed);
  printf "D E_delta, delta = %o: rank bounds over L %o %o, torsion %o\n", delta, lo, hi, Invariants(Td);
  if hi eq 0 then
    pts := [mTd(g) : g in Td];
  else
    assert lo eq hi and hi eq 1;
    fi, Ad, mw := PseudoMordellWeilGroup(Ed);
    assert fi;                              // image of finite odd index in E_delta(L)
    V, Rch := Chabauty(mw, Ecov);
    gens := [mw(Ad.i) : i in [1..Ngens(Ad)]];
    for p in PrimeDivisors(Rch) do
      if p eq 2 then continue; end if;      // odd index
      S := Saturation(gens, p);
      // p-saturated: the saturation generates the same group
      assert Set(S) eq Set(gens);
    end for;
    pts := [mw(v) : v in V];
    printf "D E_delta, delta = %o: elliptic Chabauty: %o points, R = %o\n", delta, #pts, Rch;
  end if;
  mvals := {L|}; ni := 0;
  for P in pts do
    im := Ecov(P);
    if im[2] eq 0 then ni +:= 1; ninf +:= 1; continue; end if;
    mv := im[1]/im[2];
    Include(~mvals, L!mv);
    if mv in Q then Include(~mrat, Q!mv); end if;
  end for;
  printf "D E_delta, delta = %o: m-values of the points found: %o and %o point(s) with m = oo\n", delta, Sort([Sprint(x) : x in mvals]), ni;
  if delta eq 1 then
    // E_1: rank 1 < 3 = [L:Q]; elliptic Chabauty returns R = 10 and the two points with m = oo, -1
    assert lo eq 1 and hi eq 1 and Rch eq 10 and #pts eq 2 and ni eq 1 and mvals eq {L!(-1)};
  else
    // E_kappa0: rank 0, torsion (Z/2)^2, points with m in {oo, -1, -alpha, 2 alpha^2 + 6 alpha + 3}
    assert hi eq 0 and Invariants(Td) eq [2, 2] and ni eq 1 and mvals eq {L!(-1), -al, 2*al^2 + 6*al + 3};
  end if;
end for;
printf "D (N1)&(N2a): finite rational m-values from E_delta(L): %o (plus m = oo)\n", mrat;
assert mrat subset {Q!(-1)};
assert not SpecializedPoint(EC!0, Q!(-1));     // m = -1: q reducible (the cusp (-1,0)); m = oo is the cusp O
// (N1) & (N2b): genus-2 quotient y^2 = (4m^2 - 3m + 1)(m^3 - 9m^2 - m + 1)
fb := (4*m^2 - 3*m + 1)*cb;
pts, r := Genus2AllPoints(fb);
ms := {p[1]/p[3] : p in pts | p[3] ne 0};
assert r eq 1 and ms eq {Q| 0} and #{p : p in pts | p[3] eq 0} eq 1;
printf "D (N1)&(N2b): genus-2 curve of rank %o, finite m: %o\n", r, ms;
for m0 in ms do
  ok := SpecializedPoint(EC!0, m0);
  assert not ok;
end for;
// also the discriminant of family A
assert (m + 1)^2 + 4*(m^3 - m) eq (m + 1)*(4*m^2 - 3*m + 1);
printf "D_DONE\n";

//////////////////////////////////////////////////////////////////////////////
// Part E.  3-adic input on family A (R = O): J1 = j(E_Q), J2 = j(E_{-Q}) = j(E_{Q^sigma}).
printf "Part E\n";
q, K, xq, yq := FamilyField(EC!0);
function JofXY(x0, y0)
  t := 1 - x0*(x0 - 1)^2/y0^2;
  _, _, _, _, _, c4, _, Dl := KubertData(t);
  return c4^3/Dl;
end function;
y2 := -yq - mf - 1;
J1 := JofXY(xq, yq); J2 := JofXY(xq, y2);
assert J1*J2 eq Norm(J1);
NJ := F!(J1*J2); TJ := F!(J1 + J2);
rho := F!((TJ^2 - 2*NJ)/NJ);
P2Z<aa, bb> := PolynomialRing(Z, 2);
function Homog(r)
  num := Qm!Numerator(r); den := Qm!Denominator(r);
  D := Max(Degree(num), Degree(den));
  l := LCM([Denominator(x) : x in Coefficients(num) cat Coefficients(den)]);
  num := num*l; den := den*l;
  g := GCD([Z!x : x in Coefficients(num) cat Coefficients(den)]);
  num := num/g; den := den/g;
  A := &+[P2Z | (Z!Coefficient(num, i))*aa^i*bb^(D - i) : i in [0..Degree(num)]];
  B := &+[P2Z | (Z!Coefficient(den, i))*aa^i*bb^(D - i) : i in [0..Degree(den)]];
  return A, B;
end function;
A, B := Homog(NJ);
Cr, Dr := Homog(rho);
for pt in [[1,0], [0,1], [1,1], [2,1]] do
  assert Z!Evaluate(A, pt) mod 3 ne 0;      // numerator of N(J) is a 3-adic unit
  assert Z!Evaluate(B, pt) mod 3 eq 0;      // denominator of N(J) is divisible by 3
  assert Z!Evaluate(Dr, pt) mod 3 ne 0;     // denominator of rho is a 3-adic unit
end for;
printf "E for every m in P^1(Q): v_3(N(J)(m)) < 0 and v_3(rho(m)) >= 0\n";
printf "E_DONE\n";
//////////////////////////////////////////////////////////////////////////////
// Part V.  Validation, not used in the proof (Remarks 1-3).
printf "Part V\n";
// V1 (Remark 1).  The point of family A at m = 2: K = Q(sqrt 33), E = E_7(t(x)).
okV, KV, xV, yV := SpecializedPoint(EC!0, Q!2);
assert okV and SqClassQ(Discriminant(DefiningPolynomial(KV))) eq 33;
tV := 1 - xV*(xV - 1)^2/yV^2; bV := tV^3 - tV^2; cV := tV^2 - tV;
EV := EllipticCurve([1 - cV, -bV, -bV, 0, 0]);
assert Invariants(TorsionSubgroup(EV)) eq [14];
assert not IsIntegral(jInvariant(EV));               // j(E) is not integral, so E has no CM
fV := Polynomial([bV^2, -2*bV*(1 - cV), (1 - cV)^2 - 4*bV, 4]);
rtsV := Roots(fV); assert #rtsV eq 1;
qdV := fV div (Parent(fV).1 - rtsV[1][1]);
assert IsSquare(Q!Norm(Discriminant(qdV)));          // N(Delta_E) is a square: biquadratic 2-division field
EVm := MinimalModel(EV);
dec31 := Decomposition(Integers(KV), 31);
assert #dec31 eq 2 and forall{dd : dd in dec31 | dd[2] eq 1 and Valuation(Discriminant(EVm), dd[1]) eq 0};
ap31 := [TraceOfFrobenius(EVm, dd[1]) : dd in dec31];
assert Set(ap31) eq {4, -10};                        // |a_P| <> |a_P'|: E is not geometrically isogenous to E^sigma
printf "V1 m = 2: K = Q(sqrt 33), E(K)_tors = Z/14, j(E) not integral, N(Delta_E) is a square, a_P = %o at the two primes above 31\n", ap31;
// V2 (Remark 2).  The LMFDB curves 2.0.7.1-196.2-a3 and -a4 over K = Q(a), a^2 - a + 2 = 0.
K7<a7> := NumberField(T^2 - T + 2);
assert Discriminant(Integers(K7)) eq -7;
sg7 := hom<K7 -> K7 | 1 - a7>;
assert sg7(a7) ne a7;
Ea3 := EllipticCurve([K7 | 1, a7 - 1, a7 + 1, 6 - 5*a7, 2*a7 - 1]);
Ea4 := EllipticCurve([K7 | 1, -a7, a7, 2 + 4*a7, 2 - 3*a7]);
Ea3s := EllipticCurve([sg7(x) : x in aInvariants(Ea3)]);
Ta3, mTa3 := TorsionSubgroup(Ea3);
assert Invariants(Ta3) eq [14];
assert not IsIntegral(jInvariant(Ea3));             // j(E) is not integral, so E has no CM
T2a3 := [mTa3(g) : g in Ta3 | Order(g) eq 2];
assert #T2a3 eq 1;
PK7<xk7> := PolynomialRing(K7);
E1a3 := IsogenyFromKernel(Ea3, xk7 - T2a3[1][1]/T2a3[1][3]);
assert IsIsomorphic(E1a3, Ea3s) and IsIsomorphic(Ea4, Ea3s);   // E/<T> = E^sigma = a4: degree 2
okB, KB, xB, yB := SpecializedPoint(EC![-1,0], Q!(-1/2));
assert okB;
okI, iso7 := IsIsomorphic(KB, K7);
assert okI;
ptsB := [<iso7(xB), iso7(yB)>, <sg7(iso7(xB)), sg7(iso7(yB))>];
assert EC(K7)![ptsB[1][1], ptsB[1][2]] + EC(K7)![ptsB[2][1], ptsB[2][2]] eq EC(K7)![-1, 0];
isoa3 := [];
for xy in ptsB do
  tB := 1 - xy[1]*(xy[1] - 1)^2/xy[2]^2; bB := tB^3 - tB^2; cB := tB^2 - tB;
  Append(~isoa3, IsIsomorphic(EllipticCurve([1 - cB, -bB, -bB, 0, 0]), Ea3));
end for;
assert true in isoa3;                               // E = E_7(t(x)), x of family B at m = -1/2
aiE := aInvariants(Ea3);
b2E := aiE[1]^2 + 4*aiE[2]; b4E := 2*aiE[4] + aiE[1]*aiE[3]; b6E := aiE[3]^2 + 4*aiE[5];
fE := Polynomial([b6E, 2*b4E, b2E, 4]);
rtsE := Roots(fE); assert #rtsE eq 1;
qdE := fE div (Parent(fE).1 - rtsE[1][1]);
assert SqClassQ(Norm(Discriminant(qdE))) eq 14 and SqClassQ(Norm(Evaluate(qdE, rtsE[1][1]))) eq 14;
dec3 := Decomposition(Integers(K7), 3);
assert #dec3 eq 1 and Norm(dec3[1][1]) eq 9 and Valuation(jInvariant(Ea3), dec3[1][1]) eq 0;
printf "V2 2.0.7.1-196.2-a3: torsion Z/14, no CM, E/<T> = E^sigma = a4, family B at m = -1/2, N(Delta_E) ~ 14, N(beta_E) ~ 14, v_3(j) = 0\n";
// V3 (Remark 2).  The census curve 14.c (LMFDB 9604.b3): #J(F_p) = #(Res_{K/Q} E)(F_p) for E = a3.
Qx<xx> := PolynomialRing(Q);
C14 := HyperellipticCurve(3*xx^6 - 27*xx^5 + 47*xx^4 + 64*xx^3 - 7*xx^2 + 9*xx + 15, xx^2 + xx);
J14 := Jacobian(C14);
Em3 := MinimalModel(Ea3);
tested := [];
for p in PrimesInInterval(3, 120) do
  if Valuation(Z!Discriminant(C14), p) ne 0 then continue; end if;
  decp := Decomposition(Integers(K7), p);
  if exists{dd : dd in decp | Valuation(Discriminant(Em3), dd[1]) ne 0} then continue; end if;
  assert #BaseChange(J14, GF(p)) eq &*[#Reduction(Em3, dd[1]) : dd in decp];
  Append(~tested, p);
end for;
assert #tested eq 27 and tested eq [p : p in PrimesInInterval(5, 119) | p ne 7];
printf "V3 14.c: #J(F_p) = #(Res E)(F_p) at the %o primes p < 120 of good reduction for both: %o\n", #tested, tested;
// V4 (Remark 3).  Family A at m = a/b, |a| <= 6, 1 <= b <= 6: v_P(j(E)) = v_P(j(E^sigma)) is a negative multiple
// of 7 at every P | 3.
nV4 := 0; vals3 := {};
for aN in [-6..6] do for bD in [1..6] do
  if GCD(aN, bD) ne 1 then continue; end if;
  ok4, K4, x4, y4 := SpecializedPoint(EC!0, Q!aN/bD);
  if not ok4 then continue; end if;
  t4 := 1 - x4*(x4 - 1)^2/y4^2; b4v := t4^3 - t4^2; c4v := t4^2 - t4;
  j4 := jInvariant(EllipticCurve([1 - c4v, -b4v, -b4v, 0, 0]));
  sg4 := [s : s in Automorphisms(K4) | s(K4.1) ne K4.1][1];
  for dd in Decomposition(Integers(K4), 3) do
    v4 := Valuation(j4, dd[1]);
    assert v4 lt 0 and v4 mod 7 eq 0 and v4 eq Valuation(sg4(j4), dd[1]);
    Include(~vals3, v4);
  end for;
  nV4 +:= 1;
end for; end for;
assert nV4 eq 44;
printf "V4 family A: %o quadratic points with |a|, b <= 6; valuations of j at the primes above 3: %o\n", nV4, vals3;
printf "V_DONE\n";
printf "ALL_28_CHECKS_PASSED\n";
quit;
