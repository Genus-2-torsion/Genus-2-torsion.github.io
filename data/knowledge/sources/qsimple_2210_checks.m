// Certificate for "Rational 5-torsion and the 2-rank of geometrically split Jacobians"
// (census proof qsimple-2210-obstruction).  Run: magma -b qsimple_2210_checks.m < /dev/null
// Every claim is an assert; the run ends with ALL_CHECKS_PASSED.
SetColumns(0);
SetMemoryLimit(16*10^9);
Q := Rationals();
Z := Integers();

// ---------------------------------------------------------------------------
// Part A.  #B(F_3) is never 40 for an abelian surface B over F_3.
// #B(F_3) = h(4) for h(y) = y^2 - s*y + t in Z[y] with both roots real in [-2 sqrt 3, 2 sqrt 3].
// Roots real and in the interval  <=>  s^2 >= 4t,  48 + 4t >= 0  and  192 s^2 <= (48 + 4t)^2.
vals := {};
for s in [-7..7] do for t in [-12..12] do
  if s^2 - 4*t ge 0 and 48 + 4*t ge 0 and 192*s^2 le (48 + 4*t)^2 then
    Include(~vals, 16 - 4*s + t);
  end if;
end for; end for;
// |s| <= 4 sqrt 3 < 7 and |t| <= 12 for such h, so the loop is exhaustive.
assert not 40 in vals;
assert Max(vals) eq 49;
printf "A: possible values of h(4): %o\n", Sort(SetToSequence(vals));

// ---------------------------------------------------------------------------
// Part B.  The curve X: y^2 = x^3 - 2x^2 + 1 (x = 2t, y = s): rank 0, X(Q) = Z/6, all cusps.
X := EllipticCurve([0,-2,0,0,1]);
assert CremonaReference(X) eq "20a2";
r, exact := Rank(X);
assert r eq 0 and exact;
TX, mTX := TorsionSubgroup(X);
assert Invariants(TX) eq [6];
XQ := { mTX(g) : g in TX };
assert XQ eq { X!0, X![0,1], X![0,-1], X![1,0], X![2,1], X![2,-1] };
function kubert_ab(t)   // Kubert's Z/10 normal form y^2 + (1-c)xy - by = x^3 - bx^2
  c := -t*(t-1)*(2*t-1)/(t^2-3*t+1);
  b := t^3*(t-1)*(2*t-1)/(t^2-3*t+1)^2;
  return [1-c, -b, -b, 0, 0];
end function;
function disc_ai(ai)
  a1,a2,a3,a4,a6 := Explode(ai);
  b2 := a1^2 + 4*a2; b4 := 2*a4 + a1*a3; b6 := a3^2 + 4*a6;
  b8 := a1^2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3^2 - a4^2;
  return -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6;
end function;
for x0 in [0, 1, 2] do        // t = 0, 1/2, 1: the Kubert curve is singular (cusps)
  assert disc_ai(kubert_ab(Q!x0/2)) eq 0;
end for;
print "B: X = 20a2, rank 0, X(Q) = Z/6 = {O, (0,+-1), (1,0), (2,+-1)}, all cusps (t = oo, 0, 1/2, 1)";

// ---------------------------------------------------------------------------
// Part C.  The Kubert family: point of order 10, the 2-division roots, the full 2-torsion condition.
Qt<t> := FunctionField(Q);
Et := EllipticCurve(kubert_ab(t));
P := Et![0,0];
assert 10*P eq Et!0 and 5*P ne Et!0 and 2*P ne Et!0;
e1 := (5*P)[1];
assert e1 eq -t^3*(t-1)/(t^2-3*t+1);
b2,b4,b6,b8 := Explode(bInvariants(Et));
Rt<W> := PolynomialRing(Qt);
f2 := 4*W^3 + b2*W^2 + 2*b4*W + b6;
qq := f2 div (W - e1);
assert f2 eq qq*(W - e1);
Bt := (-8*t^4+12*t^3-4*t+1)/(t^2-3*t+1)^2;
assert qq eq 4*W^2 + Bt*W + Coefficient(qq,0);
assert Discriminant(qq) eq (2*t-1)^5*(4*t^2-2*t-1)/(t^2-3*t+1)^4;   // Lemma 7(b): (2t-1)(4t^2-2t-1) times a nonzero square
ait := kubert_ab(t);
assert (1-ait[1])^2/(-ait[2] - (1-ait[1])) eq t;               // the parameter is recovered as t = c^2/(b - c)
// over the function field of X the other two roots are e2,e3 = (-B +- rt)/8, rt = (2t-1)^2 s/(t^2-3t+1)^2
PtS<SS> := PolynomialRing(Qt);
FX<s> := ext< Qt | SS^2 - (2*t-1)*(4*t^2-2*t-1) >;
rt := (2*t-1)^2*s/(t^2-3*t+1)^2;
e2 := (-Bt + rt)/8; e3 := (-Bt - rt)/8;
RX<V> := PolynomialRing(FX);
assert RX!f2 eq 4*(V - e1)*(V - e2)*(V - e3);
print "C: Kubert family checked: (0,0) has order 10, t = c^2/(b-c), x(5P) = e1, disc = (2t-1)^5(4t^2-2t-1)/(t^2-3t+1)^4, e2,e3 = (-B +- rt)/8 over s^2 = (2t-1)(4t^2-2t-1)";

// ---------------------------------------------------------------------------
// Part D.  The five families of quadratic points with P + P^sigma = R, R in X(Q) - {O}:
//          the line through -R of rational slope m.  Square classes of d = disc(q) and N_ij.
Qm<m> := FunctionField(Q);
PZ<z> := PolynomialRing(Q);
PX<XX> := PolynomialRing(Qm);
function sqclass(f)    // square class of f in Q(m)^*: <squarefree integer, monic squarefree polynomial>
  f := Qm!f;
  n := Numerator(f)*Denominator(f);
  fa := Factorization(n);
  g := PZ!1;
  for ff in fa do if IsOdd(ff[2]) then g *:= PZ!Coefficients(ff[1]); end if; end for;
  u := Q!LeadingCoefficient(n);
  un := Numerator(u)*Denominator(u);
  return <Sign(un)*Squarefree(Abs(un)), g>, n;
end function;
// R : < common class g_R of the two equal norms, class of the third norm ("-1" or "-d"), d = disc q_{R,m} >
expected := AssociativeArray();
expected[[0,1]]  := < <1, z^3 - 2*z + 1>,         "-1", z^4 + 4*z^2 - 8*z + 4 >;
expected[[2,-1]] := < <-1, z^3 - 4*z^2 + 4*z - 1>, "-d", z^4 - 8*z^2 + 8*z >;
expected[[1,0]]  := < <1, z^4 - 6*z^2 + 5>,        "-1", z^4 - 2*z^2 + 5 >;
expected[[2,1]]  := < <1, z^3 + 4*z^2 + 4*z + 1>,  "-d", z^4 - 8*z^2 - 8*z >;
expected[[0,-1]] := < <-1, z^3 - 2*z - 1>,         "-1", z^4 + 4*z^2 + 8*z + 4 >;
famdata := [];
specials := AssociativeArray();
for R in XQ do
  if R eq X!0 then continue; end if;
  key := [Z!R[1], Z!R[2]];
  S := -R;
  F := XX^3 - 2*XX^2 + 1 - (Qm!S[2] + m*(XX - Qm!S[1]))^2;
  q, rem := Quotrem(F, XX - Qm!S[1]);
  assert rem eq 0 and Degree(q) eq 2 and LeadingCoefficient(q) eq 1;
  K<xx> := ext< Qm | q >;
  tt := xx/2;
  ss := Qm!S[2] + m*(xx - Qm!S[1]);
  assert ss^2 eq (2*tt-1)*(4*tt^2-2*tt-1);
  den := tt^2 - 3*tt + 1;
  E1 := -tt^3*(tt-1)/den;
  BB := (-8*tt^4+12*tt^3-4*tt+1)/den^2;
  RT := (2*tt-1)^2*ss/den^2;
  es := [ E1, (-BB + RT)/8, (-BB - RT)/8 ];
  // the e_i are the 2-division roots of the specialized Kubert curve
  EK := EllipticCurve(kubert_ab(tt));
  c2,c4,c6 := Explode(bInvariants(EK));
  PK<WW> := PolynomialRing(K);
  assert 4*WW^3 + c2*WW^2 + 2*c4*WW + c6 eq 4*&*[ WW - e : e in es ];
  N := AssociativeArray(); Ncl := AssociativeArray(); Npoly := AssociativeArray();
  for pr in [[1,2],[1,3],[2,3]] do
    N[pr] := Norm(es[pr[1]] - es[pr[2]]);
    Ncl[pr], Npoly[pr] := sqclass(N[pr]);
  end for;
  ex := expected[key];
  assert Discriminant(q) eq Evaluate(ex[3], m);   // the column d of the table
  prs := [[1,2],[1,3],[2,3]];
  pairs := [ [a,b] : a, b in prs | a lt b and Ncl[a] eq Ncl[b] ];
  assert #pairs eq 1;                      // exactly two of the three norms share a square class
  same := pairs[1];
  assert Ncl[same[1]] eq ex[1];
  third := [ pr : pr in prs | not pr in same ][1];
  if ex[2] eq "-1" then
    assert Ncl[third] eq <-1, PZ!1>;
  else
    assert IsSquare(-N[third]/Discriminant(q));   // N_third = -d * square
  end if;
  printf "D: R = %o: q = %o, d = %o; N%o%o ~ N%o%o ~ %o*(%o); N%o%o ~ %o\n", key, q, Discriminant(q),
     same[1][1], same[1][2], same[2][1], same[2][2], ex[1][1], ex[1][2], third[1], third[2], ex[2];
  // special slopes: rational roots of every numerator and denominator involved, and of Norm(t^2-3t+1)
  polys := [ Numerator(Qm!Discriminant(q)), Denominator(Qm!Discriminant(q)), Numerator(Qm!Norm(den)), Denominator(Qm!Norm(den)) ];
  for pr in [[1,2],[1,3],[2,3]] do polys cat:= [ Numerator(N[pr]), Denominator(N[pr]) ]; end for;
  sp := {};
  for f in polys do for rr in Roots(PZ!Coefficients(f)) do Include(~sp, rr[1]); end for; end for;
  specials[key] := sp;
  Append(~famdata, <key, R, ex[1]>);
end for;

// ---------------------------------------------------------------------------
// Part E.  The curves y^2 = g_R(m): rank 0, all rational points found; their slopes.
function slope_status(R, mm)   // direct computation at one rational slope, from scratch
  S := -R;
  q := (z^3 - 2*z^2 + 1 - (S[2] + mm*(z - S[1]))^2) div (z - S[1]);
  if not IsIrreducible(q) then return "rational"; end if;
  K<w> := NumberField(q);
  tt := w/2;
  if tt^2 - 3*tt + 1 eq 0 then return "cusp"; end if;
  ai := kubert_ab(tt);
  if disc_ai(ai) eq 0 then return "singular"; end if;
  EK := EllipticCurve(ai);
  c2,c4,c6 := Explode(bInvariants(EK));
  PK<WW> := PolynomialRing(K);
  rts := [ rr[1] : rr in Roots(4*WW^3 + c2*WW^2 + 2*c4*WW + c6) ];
  assert #rts eq 3;
  d := Discriminant(q);
  NN := [ Norm(rts[1]-rts[2]), Norm(rts[1]-rts[3]), Norm(rts[2]-rts[3]) ];
  sq := func< a | a ne 0 and IsSquare(Q!a) >;
  ok := (sq(NN[1]) and (sq(NN[2]) or sq(NN[3]) or sq(d*NN[2]*NN[3]))) or
        (sq(NN[2]) and (sq(NN[1]) or sq(NN[3]) or sq(d*NN[1]*NN[3]))) or
        (sq(NN[3]) and (sq(NN[1]) or sq(NN[2]) or sq(d*NN[1]*NN[2])));
  return ok select "PASSES" else "fails";
end function;
// R : < Cremona label of the Jacobian of C_R, finite slopes of the points of C_R(Q) >, as stated in Prop. 8
expectedE := AssociativeArray();
expectedE[[0,1]]  := < "40a3", {Q| 0, 1} >;
expectedE[[0,-1]] := < "40a3", {Q| -1, 0} >;
expectedE[[2,-1]] := < "20a2", {Q| 0, 1, 2} >;
expectedE[[2,1]]  := < "20a2", {Q| -2, -1, 0} >;
expectedE[[1,0]]  := < "40a3", {Q| -1, 1} >;
allslopes := AssociativeArray();
for fd in famdata do
  key := fd[1]; R := fd[2]; cl := fd[3];
  g := cl[1]*cl[2];
  C := HyperellipticCurve(g);
  assert Genus(C) eq 1;
  pts := Points(C : Bound := 200);
  E, phi := EllipticCurve(C, pts[1]);
  rk, exact := Rank(E);
  assert rk eq 0 and exact;
  nT := #TorsionSubgroup(E);
  assert #pts eq nT;              // C is isomorphic to E over Q, so C(Q) = pts
  sl := { pp[1]/pp[3] : pp in pts | pp[3] ne 0 };
  allslopes[key] := sl;
  lab := CremonaReference(MinimalModel(E));
  assert lab eq expectedE[key][1] and sl eq expectedE[key][2];
  printf "E: R = %o: y^2 = %o is %o (rank 0, #C(Q) = %o), finite slopes %o\n", key, g, lab, nT, sl;
  for mm in sl do assert slope_status(R, mm) eq "rational"; end for;
end for;

// ---------------------------------------------------------------------------
// Part F.  Special slopes: direct computation; each gives a reducible q, a cusp (t^2-3t+1 = 0) or a singular curve.
for fd in famdata do
  key := fd[1]; R := fd[2];
  st := [ <mm, slope_status(R, mm)> : mm in Sort(SetToSequence(specials[key])) ];
  printf "F: R = %o: special slopes %o\n", key, st;
  assert forall{ x : x in st | x[2] in {"rational", "cusp", "singular"} };
end for;

// ---------------------------------------------------------------------------
// Part G.  Sanity: at m = 7/3 each family gives a quadratic point P with P + P^sigma = R
//          and a Kubert curve with E(K)_tors = Z/2 x Z/10; the generic square classes agree.
for fd in famdata do
  key := fd[1]; R := fd[2]; cl := fd[3];
  S := -R; mm := 7/3;
  q := (z^3 - 2*z^2 + 1 - (S[2] + mm*(z - S[1]))^2) div (z - S[1]);
  K<w> := NumberField(q);
  XK := BaseChange(X, K);
  Pt := XK![w, S[2] + mm*(w - S[1])];
  sig := [ a : a in Automorphisms(K) | a(w) ne w ][1];
  Ps := XK![sig(Pt[1]), sig(Pt[2])];
  assert Pt + Ps eq XK!Eltseq(R);
  EK := EllipticCurve(kubert_ab(w/2));
  assert Invariants(TorsionSubgroup(EK)) eq [2,10];
  assert slope_status(R, mm) eq "fails";
end for;
print "G: sample slope 7/3: P + P^sigma = R, E(K)_tors = [2,10], conditions fail";
print "ALL_CHECKS_PASSED";
quit;
