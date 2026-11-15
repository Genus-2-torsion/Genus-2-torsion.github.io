// qsplit_conditions_checks.m
//
// Certificate for the census proof "Necessary conditions for genus-2 Jacobians split over Q with
// torsion [25], [32], [2,28], [2,2,20] or [2,2,2,12]" (15 November 2026).
//
// Magma V2.30-1.  Run with   magma -b qsplit_conditions_checks.m
// The script reads no files.  Part H uses Magma's CremonaDatabase (all elliptic curves over Q of
// conductor below 500000).  Every computational fact used in the text is asserted, not only printed.
// The run is certified only by the final line ALL_QSPLIT_CONDITIONS_CHECKS_PASSED: Magma -b exits
// with status 0 even after a failed assertion, which stops the script before that line.
// No result depends on the generalized Riemann hypothesis or on a bounded search used as a proof of
// nonexistence; the database computation of Part H is used only relative to the completeness of
// Cremona's tables below conductor 500000, as stated in the text.
//
// Parts (numbering of the text):
//   A  Proposition 3: the torsion of the two elliptic curves and the gluing degrees
//   B  Lemma 9: universal families for Z/2 x Z/6 and Z/2 x Z/4, and the valuation lemma
//   C  Theorem 10 ([2,2,2,12]): local data at 2, 3, 5, 7, 11, 13
//   D  Theorem 11 ([2,28]): local data at 3 and 5
//   E  Theorem 5(b), Theorems 12 and 13 ([32]): lattice identities, traces, the matrix lemma
//   F  Lemmas 14 and 16 ([25] at the prime 5): Tate normal form, Morton's formulas, unit roots
//   G  Lemma 15 and Theorem 17 ([25]): the trace identity and the conductor dichotomy
//   H  Lemma 18, Theorem 19 and Remarks 3 and 4 ([25]): nonsplit primes and Cremona's tables
//   P  Remark 5: products with the product polarization (gluing degree 1)
//   V  validation only, not used in the proofs

SetColumns(0);
SetMemoryLimit(16*10^9);
Q := Rationals(); Z := Integers();

// ================================================================================================
print "PART A: parent torsion (Proposition 3)";
Mazur := [[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[12],[2,2],[2,4],[2,6],[2,8]];
function AG(inv)
  return AbelianGroup([x : x in inv | x gt 1]);
end function;
function EmbedsIn(H, G)
  if #G mod #H ne 0 then return false; end if;
  for S in Subgroups(G) do
    if #S`subgroup eq #H and Invariants(S`subgroup) eq Invariants(H) then return true; end if;
  end for;
  return false;
end function;
function OneSide(G, T1, T2)
  // is there a subgroup S of G isomorphic to T1 with G/S embedding in T2?
  for S in Subgroups(G) do
    if #S`subgroup eq #T1 and Invariants(S`subgroup) eq Invariants(T1) then
      if EmbedsIn(quo<G | S`subgroup>, T2) then return true; end if;
    end if;
  end for;
  return false;
end function;
targets := [ <"[25]",[25]>, <"[32]",[32]>, <"[2,28]",[2,28]>, <"[2,2,20]",[2,2,20]>,
             <"[2,2,2,12]",[2,2,2,12]> ];
expected := AssociativeArray();
expected["[25]"]       := {{[5]}};
expected["[32]"]       := {{[4],[8]}, {[8]}};
expected["[2,28]"]     := {{[7],[2,4]}};
expected["[2,2,20]"]   := {{[10],[2,4]}};
expected["[2,2,2,12]"] := {{[2,6],[2,4]}};
for tg in targets do
  G := AG(tg[2]); found := {};
  for i in [1..#Mazur] do for j in [i..#Mazur] do
    T1 := AG(Mazur[i]); T2 := AG(Mazur[j]);
    if EmbedsIn(T1, G) and EmbedsIn(T2, G) and OneSide(G, T1, T2) and OneSide(G, T2, T1) then
      Include(~found, {Mazur[i], Mazur[j]});
    end if;
  end for; end for;
  print "PARENT_PAIRS", tg[1], found;
  assert found eq expected[tg[1]];
end for;
// gluing-degree consequences: if p | n then E(Q)[p] ~ F(Q)[p]
function pPart(inv, m) return [Gcd(x, m) : x in inv | Gcd(x, m) gt 1]; end function;
assert pPart([7],2) ne pPart([2,4],2) and pPart([7],7) ne pPart([2,4],7);
assert pPart([10],2) ne pPart([2,4],2) and pPart([10],5) ne pPart([2,4],5);
assert pPart([2,6],2) eq pPart([2,4],2) and pPart([2,6],3) ne pPart([2,4],3);
assert pPart([2,6],4) eq [2,2] and pPart([2,4],4) eq [2,4];   // E(Q)[4] vs F(Q)[4]
assert pPart([8],8) ne pPart([4],8);                          // C8 vs C4 at level 8
// [32]: in Z/32 the elements of order 4 are the 8u with u odd, and 8u+8v lies in 16Z/32
G32 := AbelianGroup([32]); g := G32.1;
assert {x : x in G32 | Order(x) eq 4} eq {(8*u)*g : u in [1,3,5,7]};
assert forall{<u,v> : u, v in [1,3,5,7] | Order((8*u+8*v)*g) le 2};
assert forall{m : m in [1,3,5,7,9,11,13,15] | Order((8*m)*g) eq 4};
// the product of two Mazur groups never has exponent divisible by 32 or a point of order 25
assert forall{<a,b> : a, b in Mazur | Exponent(DirectSum(AG(a),AG(b))) mod 32 ne 0 and
                                      Exponent(DirectSum(AG(a),AG(b))) mod 25 ne 0};
print "PART_A_DONE";

// ================================================================================================
print "PART B: universal families and the valuation lemma (Lemma 9)";
K<t> := FunctionField(Q);
A := (t+3)*(t-1)^3; B := (t-3)*(t+1)^3;
assert A - B eq 16*t;
E := EllipticCurve([K|0,-A-B,0,A*B,0]);                 // E_t : y^2 = x(x-A)(x-B)
P3 := E![t^4-2*t^2+1, 4*(t^2-1)^2];
assert P3 ne E!0 and 3*P3 eq E!0;                        // a Q(t)-rational point of order 3
assert A ne B and A ne 0 and B ne 0;                     // full rational 2-torsion
jE := jInvariant(E);
assert jE eq 256*(A^2-A*B+B^2)^3/(A^2*B^2*(A-B)^2);
nE := Numerator(jE); dE := Denominator(jE);
assert GCD(nE,dE) eq 1 and Degree(nE) eq 24 and Degree(dE) eq 18;
PZ<s0> := PolynomialRing(Q);
discE := 16*A^2*B^2*(A-B)^2;
assert Discriminant(E) eq discE;
assert {r[1] : r in Roots(Evaluate(Numerator(discE), s0))} eq {Q| 0, 1, -1, 3, -3};
assert &+[r[2] : r in Roots(Evaluate(Numerator(discE), s0))] eq Degree(Numerator(discE));
assert {r[1] : r in Roots(Evaluate(dE, s0))} eq {Q| 0, 1, -1, 3, -3};     // poles of j: these and infinity
// Kubert's universal family for Z/2 x Z/6 and the parameter change t = (s-3)/2
Ks<s> := FunctionField(Q);
cK := (10-2*s)/(s^2-9); bK := cK + cK^2;
EK := EllipticCurve([Ks|1-cK,-bK,-bK,0,0]);
PK := EK![0,0];
assert 6*PK eq EK!0 and 2*PK ne EK!0 and 3*PK ne EK!0;
b2,b4,b6 := Explode(bInvariants(EK));
PX<X> := PolynomialRing(Ks);
assert forall{f : f in Factorization(4*X^3+b2*X^2+2*b4*X+b6) | Degree(f[1]) eq 1};   // full 2-torsion
jK := jInvariant(EK);
tt := (s-3)/2;
assert jK eq Evaluate(nE, tt)/Evaluate(dE, tt);
print "E_t: Z/2 x Z/6 over Q(t), deg j = 24; Kubert's family has j(s) = j(E_{(s-3)/2})";
// the degree of X_1(2,6) -> X(1) and X_1(2,4) -> X(1): index of {+-1} H in SL_2(Z/N)
function IndexPM(N)
  total := 0; pm := 0;
  for a,b,c,d in [0..N-1] do
    if (a*d-b*c) mod N ne 1 then continue; end if;
    total +:= 1;
    if c eq 0 and ((a eq 1 and d eq 1) or (a eq N-1 and d eq N-1)) and b mod 2 eq 0 then pm +:= 1; end if;
  end for;
  return total div pm;
end function;
assert IndexPM(6) eq 24 and IndexPM(4) eq 12;
// F_t : y^2 = x(x+1)(x+t^2)
F := EllipticCurve([K|0,t^2+1,0,t^2,0]);
P4 := F![t, t*(t+1)];
assert 2*P4 ne F!0 and 4*P4 eq F!0;
jF := jInvariant(F);
assert jF eq 256*(t^4-t^2+1)^3/(t^4*(t^2-1)^2);
// Kubert's family for Z/2 x Z/4 and the parameter change t = (4v+1)/(4v-1)
Kv<v> := FunctionField(Q);
bV := v^2-1/16;
E24 := EllipticCurve([Kv|1,-bV,-bV,0,0]);
P24 := E24![0,0];
assert 4*P24 eq E24!0 and 2*P24 ne E24!0;
b2,b4,b6 := Explode(bInvariants(E24));
PY<Y> := PolynomialRing(Kv);
assert forall{f : f in Factorization(4*Y^3+b2*Y^2+2*b4*Y+b6) | Degree(f[1]) eq 1};
tv := (4*v+1)/(4*v-1);
assert jInvariant(E24) eq Evaluate(Numerator(jF), tv)/Evaluate(Denominator(jF), tv);
print "F_t: Z/2 x Z/4 over Q(t); Kubert's family has j(v) = j(F_{(4v+1)/(4v-1)})";
// the identities of the hand proof
P2<a,b> := PolynomialRing(Z,2);
Ah := (a+3*b)*(a-b)^3; Bh := (a-3*b)*(a+b)^3;
Ch := a^4-6*a^2*b^2-3*b^4; Dh := 8*a*b^3;
assert Ah - Bh eq 16*a*b^3 and Ah eq Ch + Dh and Bh eq Ch - Dh;
assert Ah^2 - Ah*Bh + Bh^2 eq Ch^2 + 3*Dh^2;
// machine proof of the valuation lemma on all of P^1(Q_p): a disc t0 + p^k Z_p is accepted when
// v(N) is constant on it and v(D) is larger there; otherwise it is split into p subdiscs
NE := 256*(Ah^2-Ah*Bh+Bh^2)^3;  DE := Ah^2*Bh^2*(Ah-Bh)^2;
NF := 256*(a^4-a^2*b^2+b^4)^3;  DF := a^4*b^4*(a^2-b^2)^2;
assert IsHomogeneous(NE) and IsHomogeneous(DE) and TotalDegree(NE) eq TotalDegree(DE);
assert IsHomogeneous(NF) and IsHomogeneous(DF) and TotalDegree(NF) eq TotalDegree(DF);
assert jE eq Evaluate(NE,[t,K!1])/Evaluate(DE,[t,K!1]);
assert jF eq Evaluate(NF,[t,K!1])/Evaluate(DF,[t,K!1]);
Zx<x> := PolynomialRing(Z);
function V(c, p)
  if c eq 0 then return 10^9; end if;
  return Valuation(c, p);
end function;
function Accept(fN, fD, t0, k, p)
  gN := Evaluate(fN, t0 + p^k*x); gD := Evaluate(fD, t0 + p^k*x);
  cN := Coefficients(gN); cD := Coefficients(gD);
  vN0 := V(cN[1], p);
  if vN0 ge 10^9 then return false, 0; end if;
  for i in [2..#cN] do if V(cN[i],p) le vN0 then return false, 0; end if; end for;
  return Min([V(c,p) : c in cD]) gt vN0, vN0;
end function;
function Search(fN, fD, t0, k0, p, maxdepth)
  stack := [<t0,k0>]; leaves := 0; deepest := k0; vals := {Z|};
  while #stack gt 0 do
    nd := stack[#stack]; Prune(~stack);
    ok, vN := Accept(fN, fD, nd[1], nd[2], p);
    if ok then
      leaves +:= 1; deepest := Max(deepest, nd[2]); Include(~vals, vN);
    else
      assert nd[2] lt maxdepth;
      for r in [0..p-1] do Append(~stack, <nd[1] + r*p^nd[2], nd[2]+1>); end for;
    end if;
  end while;
  return leaves, deepest, vals;
end function;
cases := [ <"E_t", NE, DE, 2, {Z|8,32}, {Z|8}>, <"E_t", NE, DE, 3, {Z|0,6}, {Z|0}>,
           <"F_t", NF, DF, 3, {Z|0}, {Z|0}> ];
for cs in cases do
  p := cs[4];
  l1, d1, v1 := Search(Evaluate(cs[2],[x,Zx!1]), Evaluate(cs[3],[x,Zx!1]), 0, 0, p, 60);   // t in Z_p
  l2, d2, v2 := Search(Evaluate(cs[2],[Zx!1,x]), Evaluate(cs[3],[Zx!1,x]), 0, 1, p, 60);   // 1/t in pZ_p
  print "VALUATION_LEMMA", cs[1], "p =", p, ":", l1, "discs (depth", d1, ") and", l2, "discs (depth", d2, ")";
  assert v1 eq cs[5] and v2 eq cs[6];
end for;
print "PART_B_DONE";

// ================================================================================================
print "PART C: [2,2,2,12] local data (Theorem 10)";
expE := AssociativeArray(); expF := AssociativeArray(); expDiff := AssociativeArray();
expE[5] := {Z|}; expE[7] := {Z|-4}; expE[11] := {Z|0}; expE[13] := {Z|2};
expF[5] := {Z|-2}; expF[7] := {Z|0}; expF[11] := {Z|-4,4}; expF[13] := {Z|-2,6};
expDiff[5] := {Z|4,8}; expDiff[7] := {Z|4,8,12}; expDiff[11] := {Z|4,8,12,16}; expDiff[13] := {Z|4,8,12,16};
for p in [5,7,11,13] do
  k := GF(p); sq := {k| y^2 : y in k};
  R<T> := PolynomialRing(k);
  eTr := {Z|}; fTr := {Z|}; eTr2 := {Z|}; fTr2 := {Z|};
  for a0 in k do for b0 in k do
    if 4*a0^3+27*b0^2 eq 0 then continue; end if;
    // method 1: Legendre sums, roots, and the halving criterion for a point of order 4
    N := 1;
    for y in k do w := y^3+a0*y+b0; if w eq 0 then N +:= 1; elif w in sq then N +:= 2; end if; end for;
    rts := [r[1] : r in Roots(T^3+a0*T+b0)];
    if #rts eq 3 then
      if N mod 3 eq 0 then Include(~eTr, p+1-N); end if;
      if exists{i : i in [1..3] | forall{j : j in [1..3] | j eq i or (rts[i]-rts[j]) in sq}} then
        Include(~fTr, p+1-N);
      end if;
    end if;
    // method 2: the group structure
    C := EllipticCurve([k|a0,b0]); Gp := AbelianGroup(C); ii := Invariants(Gp);
    twrank := #[m : m in ii | m mod 2 eq 0];
    if twrank eq 2 and #Gp mod 3 eq 0 then Include(~eTr2, p+1-#Gp); end if;
    if twrank eq 2 and Exponent(Gp) mod 4 eq 0 then Include(~fTr2, p+1-#Gp); end if;
  end for; end for;
  hasseE := {Z| y : y in [-2*p..2*p] | y^2 le 4*p and (p+1-y) mod 12 eq 0};
  hasseF := {Z| y : y in [-2*p..2*p] | y^2 le 4*p and (p+1-y) mod 8 eq 0};
  assert hasseE eq expE[p] and hasseF eq expF[p];     // used in the proof (torsion injection + Hasse)
  assert eTr eq hasseE and fTr eq hasseF and eTr2 eq hasseE and fTr2 eq hasseF;   // and sharp
  multF := {Z| p+1, -(p+1)};
  multE := ((2*(p+1)) mod 3 eq 0) select {Z| p+1, -(p+1)} else {Z| p+1};
  assert (p in {7,13}) eq (multE eq {Z|p+1});          // E split at 7 and 13 by its 3-torsion
  diffs := {Z| AbsoluteValue(y-z) : y in hasseE, z in hasseF join multF} join
           {Z| AbsoluteValue(y-z) : y in hasseF, z in hasseE join multE};
  assert 0 notin diffs and diffs eq expDiff[p];
  assert forall{dd : dd in diffs | Set(PrimeDivisors(dd)) subset {2,3}};
  // l = p: ordinary good reduction has unramified constituent a_p mod p, multiplicative +-1;
  // supersingular good reduction (a_p = 0 mod p) has irreducible E[p]
  unE := {k!y : y in hasseE | y mod p ne 0}; ssE := exists{y : y in hasseE | y mod p eq 0};
  unF := {k!y : y in hasseF | y mod p ne 0}; ssF := exists{y : y in hasseF | y mod p eq 0};
  assert unE meet {k|1,-1} eq {} and unF meet {k|1,-1} eq {} and unE meet unF eq {};
  assert ssE eq (p eq 11) and ssF eq (p eq 7);
  uE := AssociativeArray(); uF := AssociativeArray();
  uE[5] := {k|}; uE[7] := {k|3}; uE[11] := {k|}; uE[13] := {k|2};
  uF[5] := {k|3}; uF[7] := {k|}; uF[11] := {k|4,-4}; uF[13] := {k|6,11};
  assert unE eq uE[p] and unF eq uF[p];
  print "p =", p, " E good traces", hasseE, " F good traces", hasseF, " differences", diffs;
end for;
// local torsion bounds at 5, 7, 11, 13 (formal group torsion-free, E0/E1 of order p+1 when nonsplit,
// component group of order at most 2 when nonsplit, at most 4 when additive)
twoPart := func<m | 2^Valuation(m,2)>;
assert 12 gt 4 and 8 gt 4;                                 // additive reduction impossible
for p in [5,7,11,13] do
  nsF := twoPart(2*(p+1)) ge 8; nsE := (2*(p+1)) mod 3 eq 0;
  assert nsF eq (p in {7,11}) and nsE eq (p in {5,11});    // F split at 5, 13; E split at 7, 13
end for;
// equal splitting characters at 5 (E), 7 (F), 11 (both) for l != p: opposite characters differ by 2(p+1)
assert forall{p : p in [5,7,11] | Set(PrimeDivisors(2*(p+1))) subset {2,3}};
assert 7 in PrimeDivisors(2*(13+1));      // the trace argument would fail at 13; there torsion decides
// p = 2
k2 := GF(2); counts := {Z|};
for a1,a2,a3,a4,a6 in [0,1] do
  ok, EE := IsEllipticCurve([k2|a1,a2,a3,a4,a6]);
  if not ok then continue; end if;
  Include(~counts, 1 + #[ <y,z> : y,z in k2 | z^2+a1*y*z+a3*z eq y^3+a2*y^2+a4*y+a6 ]);
end for;
assert Max(counts) le 5 and {N : N in counts | N mod 4 eq 0} eq {4};   // good F at 2: #F(F_2) = 4
assert 2+1-4 eq -1 and {AbsoluteValue(-1-3), AbsoluteValue(-1+3)} eq {4,2};
assert 2*twoPart(3)*2 eq 4;                                // nonsplit F at 2: 2-primary torsion <= 4 < 8
assert PrimeDivisors(3-(-3)) eq [2,3];                     // E nonsplit vs F split at 2
// p = 3
assert (2*(3+1)) mod 3 ne 0;                               // nonsplit E at 3: no point of order 3
assert PrimeDivisors(4-(-4)) eq [2];                       // F nonsplit vs E split at 3
assert &*[2,3,5,7,11,13] eq 30030 and 30030^2 eq 901800900;
print "PART_C_DONE";

// ================================================================================================
print "PART D: [2,28] local data (Theorem 11)";
function Counts(p)
  k := GF(p); res := {Z|};
  for a1,a2,a3,a4,a6 in [0..p-1] do
    ok, EE := IsEllipticCurve([k|a1,a2,a3,a4,a6]);
    if not ok then continue; end if;
    Include(~res, 1 + #[ <y,z> : y,z in k | z^2+a1*y*z+a3*z eq y^3+a2*y^2+a4*y+a6 ]);
  end for;
  return res;
end function;
for p in [3,5] do
  S := Counts(p);
  t7 := {Z| p+1-N : N in S | N mod 7 eq 0};
  t8 := {Z| p+1-N : N in S | N mod 8 eq 0};
  h7 := {Z| y : y in [-2*p..2*p] | y^2 le 4*p and (p+1-y) mod 7 eq 0};
  h8 := {Z| y : y in [-2*p..2*p] | y^2 le 4*p and (p+1-y) mod 8 eq 0};
  assert t7 eq h7 and t8 eq h8;
  if p eq 3 then assert Max(S) eq 7 and h7 eq {Z|-3} and h8 eq {Z|}; end if;
  if p eq 5 then assert h7 eq {Z|-1} and h8 eq {Z|-2}; end if;
  print "p =", p, " traces with 7 | N:", h7, " traces with 8 | N:", h8;
end for;
// p = 3: E good is supersingular; comparisons with F multiplicative
assert (-3) mod 3 eq 0;
assert {(-3)-4, (-3)+4} eq {-7, 1};                        // only l = 7, and 7 does not divide n
assert 7 gt 4 and (2*(3+1)) mod 7 ne 0;                    // E bad at 3 is split
assert 4-(-4) eq 8 and GF(3)!1 ne GF(3)!(-1);              // F split at 3
// p = 5
assert twoPart(5+1)*2 eq 4;                                // F nonsplit at 5: 2-primary torsion <= 4 < 8
assert (-2) mod 5 ne 0 and (-1) mod 5 ne 0;                // good reductions are ordinary
assert (-1)-(-2) eq 1 and (-1)-6 eq -7 and 6-(-2) eq 8;    // trace differences for l != 5
assert GF(5)!(-1) eq GF(5)!4 and GF(5)!(-2) eq GF(5)!3 and #{GF(5)|4,3,1} eq 3;   // l = 5
assert (2*(5+1)) mod 7 ne 0;                               // E bad at 5 is split
assert 15^2 eq 225;
print "PART_D_DONE";

// ================================================================================================
print "PART E: [32] (Theorem 5(b), Theorems 12 and 13)";
// Theorem 5(b), step 1: primitive vectors of order 4 fixed by c = [[1,1],[0,-1]] are odd multiples of e1
R4 := Integers(4);
c1 := Matrix(R4,2,2,[1,1,0,-1]);
fixed := [ w : w in [Vector(R4,[y,z]) : y,z in [0..3]] | w*Transpose(c1) eq w and
           not (Z!w[1] mod 2 eq 0 and Z!w[2] mod 2 eq 0) ];
assert Set(fixed) eq {Vector(R4,[1,0]), Vector(R4,[3,0])};
Pq<a,b0,c0,d> := PolynomialRing(Q,4);
gm := Matrix(Pq,2,2,[a,b0,c0,d]);
Pc := Matrix(Pq,2,2,[1,-1/2,0,1]);                          // e+ = e1, e- = e2 - e1/2
gn := Pc^-1*gm*Pc;
assert gn eq Matrix(Pq,2,2,[a+c0/2, b0+(d-a)/2-c0/4, c0, d-c0/2]);
assert Pc^-1*Matrix(Pq,2,2,[1,1,0,-1])*Pc eq Matrix(Pq,2,2,[1,0,0,-1]);
Pz<y1,y2,y3,y4> := PolynomialRing(Q,4);
sub := [1+4*y1, y2, 4*y3, 1+2*y4];                          // a = 1 mod 4, c0 = 0 mod 4, d odd
for ent in Eltseq(gn) do
  assert forall{cf : cf in Coefficients(Evaluate(ent, sub)) | IsIntegral(cf)};
end for;
// step 2: the corner identities giving 2 R0 in R
Pg<g11,g12,g21,g22> := PolynomialRing(Q,4);
Gm := Matrix(Pg,2,2,[g11,g12,g21,g22]);
cd := Matrix(Pg,2,2,[1,0,0,-1]); I2 := IdentityMatrix(Pg,2);
ep := (I2+cd)/2; em := (I2-cd)/2;
assert (I2+cd)*Gm - g11*(I2+cd) eq 2*ep*Gm*em;
assert (I2-cd)*Gm - g22*(I2-cd) eq 2*em*Gm*ep;
// Theorem 13: S c = chi(c) c S and det S = -1 mod 4 force S = I mod 2 (b = 0: 8 solutions, b = 1: 4 solutions)
nsol := [0,0];
for bb in [0..1] do
  cb := Matrix(R4,2,2,[1,bb,0,-1]);
  sgn := (bb mod 2 eq 0) select R4!1 else R4!(-1);
  for s11,s12,s21,s22 in [0..3] do
    Sm := Matrix(R4,2,2,[s11,s12,s21,s22]);
    if Determinant(Sm) ne R4!(-1) or Sm*cb ne sgn*cb*Sm then continue; end if;
    nsol[1 + (bb mod 2)] +:= 1;
    assert ChangeRing(ChangeRing(Sm,Z),Integers(2)) eq IdentityMatrix(Integers(2),2);
  end for;
end for;
print "MATRIX_LEMMA: solutions for b = 0 (chi(c) = +1) and b = 1 (chi(c) = -1):", nsol;
assert nsol eq [8, 4];
cb := Matrix(R4,2,2,[1,1,0,-1]);
sols := [ Matrix(R4,2,2,[s11,s12,s21,s22]) : s11,s12,s21,s22 in [0..3] |
          Determinant(Matrix(R4,2,2,[s11,s12,s21,s22])) eq R4!(-1) and
          Matrix(R4,2,2,[s11,s12,s21,s22])*cb eq -cb*Matrix(R4,2,2,[s11,s12,s21,s22]) ];
assert forall{Sm : Sm in sols | Sm[2,1] eq 2*Sm[1,1] and Sm[2,2] eq -Sm[1,1] and IsOdd(Z!Sm[1,1]) and IsEven(Z!Sm[1,2])};
// the product model (P,Q) -> (P+Q, P-Q): kernel on E[2] x E[2] is the diagonal, pull-back 2 x product
M := Matrix(Z,2,2,[1,1,1,-1]);
assert Transpose(M)*M eq 2*IdentityMatrix(Z,2) and Determinant(M)^2 eq 4;
assert [ <y,z> : y,z in [0,1] | (y+z) mod 2 eq 0 and (y-z) mod 2 eq 0 ] eq [<0,0>,<1,1>];
// torsion over quadratic fields: the only group in the list containing Z/2 x Z/8 is Z/2 x Z/8
quad := [[m] : m in [1..16]] cat [[18]] cat [[2,2*m] : m in [1..6]] cat [[3,3],[3,6],[4,4]];
assert [ T : T in quad | #T eq 2 and T[1] mod 2 eq 0 and T[2] mod 8 eq 0 ] eq [[2,8]];
assert Exponent(AG([2,8])) eq 8;
// Theorem 12: on the coset chi = -1, 2(1 + det) = 0 mod 8 iff det = 3 mod 4; p = 3 mod 8 gives 8 mod 16
assert forall{dt : dt in [1,3,5,7] | ((2*(1+dt)) mod 8 eq 0) eq (dt mod 4 eq 3)};
assert forall{p : p in [3,11,19,43,59,67,83] | (1+p) mod 8 eq 4 and (2*(1+p)) mod 16 eq 8};
assert forall{r : r in [0..7] | ((2*(4+8*r)) mod 16) eq 8};
print "PART_E_DONE";

// ================================================================================================
print "PART F: [25] at the prime 5 (Lemmas 14 and 16)";
Pd<dd> := PolynomialRing(Z);
Ed := EllipticCurve([FieldOfFractions(Pd)|1-dd,-dd,-dd,0,0]);
c4,c6 := Explode(cInvariants(Ed)); Dl := Discriminant(Ed);
c4 := Pd!c4; c6 := Pd!c6; Dl := Pd!Dl;
assert Dl eq dd^5*(dd^2-11*dd-1);
assert c4 eq dd^4-12*dd^3+14*dd^2+12*dd+1;
assert c6 eq -dd^6+18*dd^5-75*dd^4-75*dd^2-18*dd-1;
assert Ed![0,0] ne Ed!0 and 5*(Ed![0,0]) eq Ed!0;
// (x,y) -> (x/d^2 - 1/d, y/d^3 - x/d^2) maps E_d isomorphically to E_{-1/d}, for every d != 0
Kd<dv> := FunctionField(Q);
Edv := EllipticCurve([Kd|1-dv,-dv,-dv,0,0]); ev := -1/dv;
assert IsIsomorphic(Edv, EllipticCurve([Kd|1-ev,-ev,-ev,0,0]));
Pxy<xx,yy> := PolynomialRing(Kd,2);
fd := yy^2+(1-dv)*xx*yy-dv*yy-xx^3+dv*xx^2;
XX := xx/dv^2-1/dv; YY := yy/dv^3-xx/dv^2;
fe := YY^2+(1-ev)*XX*YY-ev*YY-XX^3+ev*XX^2;
assert fe eq fd/dv^6;
// v_5(d) > 0: c4 is a unit and v(Delta) = 5 v(d); v_5(d) < 0 is reduced to this case by d -> -1/d
assert Evaluate(c4,0) eq 1 and Coefficient(Dl,5) eq -1 and Coefficient(Dl,6) eq -11 and Coefficient(Dl,7) eq 1;
assert Evaluate(c4,0) mod 5 eq 1;
// the reduction of E_d with 5 | d is y^2 + xy = x^3, whose node at (0,0) has tangents y(y+x) = 0
P5xy<xr,yr> := PolynomialRing(GF(5),2);
fr := yr^2 + (1-0)*xr*yr - 0*yr - (xr^3 - 0*xr^2);
assert fr eq yr^2 + xr*yr - xr^3;
assert Evaluate(fr,[0,0]) eq 0 and Evaluate(Derivative(fr,xr),[0,0]) eq 0 and Evaluate(Derivative(fr,yr),[0,0]) eq 0;
assert {f[1] : f in Factorization(yr^2+xr*yr)} eq {yr, xr+yr};
// d a unit at 5
assert dd^2-11*dd-1 - (dd-3)^2 eq -5*(dd+2);
Rk<kk> := PolynomialRing(Z);
assert Evaluate(dd^2-11*dd-1, 3+5*kk) eq 25*(kk^2-kk-1);
assert kk^2-kk-1 eq (kk-3)^2+5*(kk-2);
assert Evaluate(dd^2-11*dd-1, 18+25*kk) eq 125*(1+5*kk+5*kk^2);
assert Content(Evaluate(c4,3+5*kk)) mod 5 eq 0 and Content(Evaluate(c6,3+5*kk)) mod 5 eq 0;
assert Content(Evaluate(c6,18+25*kk)) mod 25 eq 0;            // type III: v(c6) >= 2, so j = 1728 mod 5
assert (Evaluate(c4,3) div 5) mod 5 ne 0 and Evaluate(c4,3) mod 5 eq 0;
assert 2*68-11 eq 125 and (2*3-11) mod 5 eq 0;                // the exceptional disk d = 68 mod 125
// the smooth fibres over F_5: d = 1, 4 give 5 points and d = 2 gives 10; d = 0, 3 are singular
for d0 in [0..4] do
  ok, Ec := IsEllipticCurve([GF(5)|1-d0,-d0,-d0,0,0]);
  assert ok eq (d0 in {1,2,4});
  if ok then assert #Ec eq (d0 eq 2 select 10 else 5); end if;
end for;
// j = 1728 is ordinary in characteristic 5
for a0 in [1..4] do assert TraceOfFrobenius(EllipticCurve([GF(5)|a0,0])) mod 5 ne 0; end for;
// unit roots: alpha = 1 + 5c mod 25 with alpha^2 - a_5 alpha + 5 = 0; c = 4 for 5 points, 3 for 10
cval := AssociativeArray();
for N in [5,10] do
  a5 := 6-N;
  lifts := [y : y in [0..24] | (y^2-a5*y+5) mod 25 eq 0 and y mod 5 eq 1];
  assert #lifts eq 1;
  cval[N] := ((lifts[1]-1) div 5) mod 5;
end for;
assert cval[5] eq 4 and cval[10] eq 3;
lam := func<cE, cF | (cF*Modinv(cE,5)) mod 5>;
assert lam(cval[5],cval[10]) eq 2 and lam(cval[10],cval[5]) eq 3;
// Morton's formula: in Q(zeta_5)(u), with d(u) = -(eps^5 u^5 - eps^-5)/(u^5+1), u^5 = -(d-alpha)/(d-beta)
// and x(u) = N(u)/D(u) is a root of the 5-division polynomial of E_{d(u)}
K5<z> := CyclotomicField(5);
sq5 := z-z^2-z^3+z^4;
assert sq5^2 eq 5;
eps := (-1+sq5)/2;
Ru<u> := PolynomialRing(K5);
Nu := -eps^4*&*[u-z^j*(1+z)^2 : j in [0..3]];
Du := (u+z^2)*(u+z^3)*(u+1)^2;
assert Degree(Nu) eq 4 and Degree(Du) eq 4 and Degree(GCD(Nu,Du)) eq 0;
Fu := FieldOfFractions(Ru);
bu := (eps^5*u^5-eps^(-5))/(u^5+1); du := -Fu!bu;
alpha := (11+5*sq5)/2; beta := (11-5*sq5)/2;
assert alpha*beta eq -1 and alpha+beta eq 11;                 // the roots of d^2 - 11d - 1
assert u^5 eq -(du-alpha)/(du-beta);
Eu := EllipticCurve([Fu|1-du,-du,-du,0,0]);
xu := Fu!(Nu/Du);
assert Evaluate(DivisionPolynomial(Eu,5),xu) eq 0;
assert xu ne 0 and xu ne du;
// the valuation computation: for d = 3 + 5k, q(d) = (2k-1-sqrt5)/(2k-1+sqrt5)
Pk<kq> := PolynomialRing(K5);
qd := ((3+5*kq)-alpha)/((3+5*kq)-beta);
assert qd eq (2*kq-1-sq5)/(2*kq-1+sq5);
assert qd - 1 eq -2*sq5/(2*kq-1+sq5);
Pw<ww> := PolynomialRing(K5);
qw := (5*ww-sq5)/(5*ww+sq5);                                   // 2k - 1 = 5w
assert -qw eq (1-sq5*ww)/(1+sq5*ww) and -qw-1 eq -2*sq5*ww/(1+sq5*ww);
print "PART_F_DONE";

// ================================================================================================
print "PART G: [25] trace identity and conductor dichotomy (Lemma 15, Theorem 17)";
Pt<a,b,c,d,w> := PolynomialRing(Z,5);
Mt := Matrix(Pt,2,2,[1+5*a,b,5*c,d]);
assert Trace(Mt) - 1 - Determinant(Mt) eq 5*(a*(1-d)+b*c);
assert (a+b*w)*(1-d) + b*(c+(d-1)*w) eq a*(1-d)+b*c;
hc := AssociativeArray();
hc[2] := [5]; hc[3] := [5]; hc[5] := [5,10]; hc[7] := [5,10]; hc[11] := [10,15]; hc[13] := [10,15,20];
for p in [2,3,5,7,11,13] do
  L := [N : N in [1..4*p] | N mod 5 eq 0 and (N-(p+1))^2 le 4*p];
  assert L eq hc[p] and Max(L) lt 25;                          // a good count is never divisible by 25
  assert forall{N : N in L | Valuation(N,5) eq 1};             // equal counts give exponent 5
end for;
assert forall{p : p in [2,3,7,11,13] | (2*(p+1)) mod 5 ne 0};  // bad reduction there is split
assert (2*(5+1)) mod 5 ne 0;                                   // no nonsplit reduction at 5
ratio := func<n, m | ((m div 5)*Modinv(n div 5, 5)) mod 5>;    // #F = ratio * #E mod 25
assert {ratio(10,15), ratio(15,10)} eq {4};
assert {ratio(5,10), ratio(10,5)} eq {2,3};
assert 2*3*11 eq 66 and 2*3*5*7 eq 210;
print "PART_G_DONE";

// ================================================================================================
print "PART H: [25] nonsplit obstruction and Cremona's tables (Lemma 18, Theorem 19)";
// an infinite family excluded only by Lemma 18
assert (Z!Evaluate(c4,2)) mod 19 ne 0;
assert not IsSquare(GF(19)!(-Evaluate(c6,2)));
assert (78^2-11*78-1) mod 361 eq 171 and 171 eq 9*19;
assert 78 mod 19 eq 2 and 2166 mod 361 eq 0 and 10830 mod 361 eq 0;
assert 2244 eq 78+2166 and 10830 eq 5*2166;
assert 2244 mod 6 eq 0 and 10830 mod 6 eq 0 and 2244 mod 5 eq 4 and 10830 mod 5 eq 0;
assert Valuation(19+1,5) eq 1;
assert Evaluate(c4,0) mod 2 eq 1 and Evaluate(c4,0) mod 3 eq 1;
for d0 in [78, 2244] do
  Ec := MinimalModel(EllipticCurve([Q|1-d0,-d0,-d0,0,0]));
  L2 := LocalInformation(Ec,2); L3 := LocalInformation(Ec,3); L19 := LocalInformation(Ec,19);
  assert L2[3] eq 1 and L2[6] and L3[3] eq 1 and L3[6];        // split multiplicative at 2 and 3
  assert L19[2] eq 1 and L19[3] eq 1 and not L19[6];           // nonsplit I1 at 19
  if d0 eq 2244 then assert Conductor(Ec) mod 5 ne 0 and TraceOfFrobenius(Ec,5) eq 1; end if;
  print "d =", d0, LocalInformation(Ec,2), LocalInformation(Ec,3), L19;
end for;
// the complete scan of Cremona's tables over conductors divisible by 6 and not by 25
D := CremonaDatabase();
assert LargestConductor(D) eq 499999;
ncurves := 0; nclasses := 0; C5 := []; C10 := [];
for N in [6..LargestConductor(D) by 6] do
  if N mod 25 eq 0 then continue; end if;
  nclasses +:= NumberOfIsogenyClasses(D,N);
  for Ec in EllipticCurves(D,N) do
    ncurves +:= 1;
    inv := Invariants(TorsionSubgroup(Ec));
    if inv eq [5] then Append(~C5, <N, CremonaReference(D,Ec)>);
    elif inv eq [10] then Append(~C10, <N, CremonaReference(D,Ec)>); end if;
  end for;
end for;
print "SCAN", ncurves, "curves in", nclasses, "classes;", #C5, "with torsion Z/5 and", #C10, "with torsion Z/10";
assert ncurves eq 1472509 and nclasses eq 980367 and #C5 eq 464 and #C10 eq 32;
function Witnesses(Ec)
  W := [];
  for p in PrimeDivisors(Conductor(Ec)) do
    if p mod 5 ne 4 or p mod 25 eq 24 then continue; end if;   // v_5(p+1) = 1
    Lp := LocalInformation(Ec,p);
    if Lp[3] eq 1 and not Lp[6] and Lp[2] mod 5 ne 0 then Append(~W, p); end if;
  end for;
  return W;
end function;
wit := [Witnesses(EllipticCurve(D,r[2])) : r in C5];
assert #[w : w in wit | #w gt 0] eq 242;
dich := [i : i in [1..#C5] | C5[i][1] mod 66 eq 0 or C5[i][1] mod 210 eq 0];
surv := [i : i in dich | #wit[i] eq 0];
assert #dich eq 155 and #surv eq 76;
mins := Min([C5[i][1] : i in surv]);
assert mins eq 8778 and [C5[i][2] : i in surv | C5[i][1] eq mins] eq ["8778x1"];
small := [i : i in dich | C5[i][1] lt 8778];
assert [C5[i][2] : i in small] eq ["1254k1","1914o1","1914p1"] and [wit[i] : i in small] eq [[19],[29],[29]];
for i in small do
  Ec := EllipticCurve(D,C5[i][2]); p := wit[i][1]; Lp := LocalInformation(Ec,p);
  assert Lp[2] eq 1 and Lp[3] eq 1 and not Lp[6] and Valuation(p+1,5) eq 1;
  cls := [Fc : Fc in EllipticCurves(D,C5[i][1]) | IsIsogenous(Ec,Fc)];
  assert #cls eq 2;
  other := [Fc : Fc in cls | not IsIsomorphic(Fc,Ec)][1];
  ok, phi := IsIsogenous(Ec, other);
  assert ok and Degree(phi) eq 5;
  assert #TorsionSubgroup(other) eq 1;
  Lo := LocalInformation(other,p); assert Lo[2] eq 5 and not Lo[6];
  assert #Roots(PolynomialRing(Q)!DivisionPolynomial(Ec,2)) eq 0 and #TwoTorsionSubgroup(Ec) eq 1;
  print "LOW_PARENT", C5[i][2], "witness", p, Lp, "; class partner", CremonaReference(D,other), "with trivial torsion and", Lo;
end for;
// the trace test on pairs: #E(F_p)/5 mod 5 at the primes 7 <= p <= 59
prs := [p : p in PrimesUpTo(59) | p gt 5];
function ClassLabel(str)
  while str[#str] in "0123456789" do str := str[1..#str-1]; end while;
  return str;
end function;
function Vals(lab)
  Ec := EllipticCurve(D,lab); N := Conductor(Ec);
  return [N mod p eq 0 select -1 else ((p+1-TraceOfFrobenius(Ec,p)) div 5) mod 5 : p in prs];
end function;
function PairTest(Lst)
  vals := [Vals(r[2]) : r in Lst];
  noniso := 0; same := []; maxp := 0;
  for i in [1..#Lst] do for j in [i+1..#Lst] do
    possible := {1,2,3,4}; first := 0;
    for m in [1..#prs] do
      a := vals[i][m]; b := vals[j][m];
      if a eq -1 or b eq -1 then continue; end if;
      possible := {k : k in possible | b eq (k*a) mod 5};
      if IsEmpty(possible) and first eq 0 then first := prs[m]; end if;
    end for;
    if ClassLabel(Lst[i][2]) eq ClassLabel(Lst[j][2]) then
      Append(~same, <Lst[i][2], Lst[j][2], possible>);
    else
      noniso +:= 1; assert IsEmpty(possible); maxp := Max(maxp, first);
    end if;
  end for; end for;
  return noniso, same, maxp;
end function;
// a curve paired with itself passes only with lambda = 1
selfok := true;
for r in C5 do
  vv := Vals(r[2]);
  possible := {k : k in {1,2,3,4} | forall{m : m in [1..#prs] | vv[m] eq -1 or vv[m] eq (k*vv[m]) mod 5}};
  selfok := selfok and possible eq {1};
end for;
assert selfok;
n1, same1, m1 := PairTest(C5);
print "PAIRS Z/5:", n1, "pairs in different classes, all excluded (largest prime needed", m1, "); same class:", same1;
assert n1 eq 107415 and m1 eq 59 and same1 eq [<"165066d1","165066d2",{1}>];
n2, same2, m2 := PairTest(C5 cat C10);
print "PAIRS Z/5 and Z/10:", n2, "pairs in different classes, all excluded; same class:", #same2;
assert n2 eq 122743 and #same2 eq 17 and forall{r : r in same2 | r[3] eq {1}};
// 66c1 satisfies (i)-(iii) and has no witness prime, but has torsion Z/10
E66 := EllipticCurve(D,"66c1");
assert Invariants(TorsionSubgroup(E66)) eq [10] and Conductor(E66) eq 66 and Witnesses(E66) eq [];
assert forall{p : p in [2,3,11] | LocalInformation(E66,p)[3] eq 1 and LocalInformation(E66,p)[6]};
assert TraceOfFrobenius(E66,5) eq 6-10;
assert 8778*500000 eq 4389000000;
print "PART_H_DONE";

// ================================================================================================
print "PART P: products with the product polarization (Remark 5)";
E7 := EllipticCurve(D,"26b1"); F24 := EllipticCurve(D,"15a1");
E26 := EllipticCurve(D,"30a2"); E10 := EllipticCurve(D,"66c1");
assert Invariants(TorsionSubgroup(E7)) eq [7] and Invariants(TorsionSubgroup(F24)) eq [2,4];
assert Invariants(TorsionSubgroup(E26)) eq [2,6] and Invariants(TorsionSubgroup(E10)) eq [10];
assert Invariants(DirectSum(TorsionSubgroup(E7), TorsionSubgroup(F24))) eq [2,28];
assert Invariants(DirectSum(TorsionSubgroup(E26), TorsionSubgroup(F24))) eq [2,2,2,12];
assert Invariants(DirectSum(TorsionSubgroup(E10), TorsionSubgroup(F24))) eq [2,2,20];
assert Conductor(E7) eq 26 and Conductor(F24) eq 15 and Conductor(E26) eq 30;
assert Conductor(E7)*Conductor(F24) eq 390 and Conductor(E26)*Conductor(F24) eq 450;
assert 450 mod 901800900 ne 0 and 390 mod 225 ne 0;
print "PRODUCTS: 26b1 x 15a1 has torsion [2,28]; 30a2 x 15a1 has torsion [2,2,2,12] and conductor 450; 66c1 x 15a1 has torsion [2,2,20]";
print "PART_P_DONE";

// ================================================================================================
print "PART V: validation only (not used in the proofs)";
// (V1) the trace relation of Theorem 12 on a pair with Z/8 and Z/4 related by the twist chi_{-1}
uq := Q!(13/18); bq := (2*uq-1)*(uq-1); cq := bq/uq;
EV := EllipticCurve([Q|1-cq,-bq,-bq,0,0]); FV := QuadraticTwist(EV,-1);
assert Invariants(TorsionSubgroup(EV)) eq [8] and Invariants(TorsionSubgroup(FV)) eq [4];
NEV := Conductor(EV); NFV := Conductor(FV); cnt := 0;
for p in PrimesUpTo(2000) do
  if p eq 2 or NEV mod p eq 0 or NFV mod p eq 0 then continue; end if;
  ae := TraceOfFrobenius(EV,p); af := TraceOfFrobenius(FV,p);
  assert (ae-(1+p)) mod 8 eq 0 and af eq KroneckerSymbol(-1,p)*ae;
  if p mod 8 eq 3 then assert (ae-af) mod 16 eq 8; cnt +:= 1; end if;
end for;
assert cnt eq 76;
// (V2) the valuation lemma on a box of parameters
for bb in [1..60] do for aa in [-60..60] do
  if Gcd(aa,bb) ne 1 then continue; end if;
  t0 := aa/bb;
  if t0 notin {Q|0,1,-1} then
    assert Valuation(256*(t0^4-t0^2+1)^3/(t0^4*(t0^2-1)^2),3) lt 0;
  end if;
  if t0 in {Q|0,1,-1,3,-3} then continue; end if;
  A0 := (t0+3)*(t0-1)^3; B0 := (t0-3)*(t0+1)^3;
  j0 := 256*(A0^2-A0*B0+B0^2)^3/(A0^2*B0^2*(A0-B0)^2);
  assert Valuation(j0,2) lt 0 and Valuation(j0,3) lt 0;
end for; end for;
// (V3) reduction types at 5 of sample Tate normal forms, and the fifth roots of q(d) over Q_5(zeta_5)
types := ["I0","I0","II","I5","III","II","III","III"];
for i in [1..8] do
  d0 := [1,2,3,5,18,48,68,193][i];
  Lc := LocalInformation(MinimalModel(EllipticCurve([Q|1-d0,-d0,-d0,0,0])),5);
  assert Sprint(Lc[5]) eq types[i];
end for;
P5 := pAdicField(5,80); B5<v5> := PolynomialRing(P5);
L5<pi5> := ext<P5|((v5+1)^5-1) div v5>;
zz := 1+pi5; ss := zz-zz^2-zz^3+zz^4;
assert Valuation(ss^2-5) ge 60;
aa5 := (11+5*ss)/2; bb5 := (11-5*ss)/2;
C5x<T5> := PolynomialRing(L5);
for d0 in [1,2,3,5,18,48,68,193] do
  qq := (d0-aa5)/(d0-bb5);
  degs := Sort([Degree(f[1]) : f in Factorization(T5^5-qq)]);
  additive := d0 mod 5 eq 3;
  exceptional := additive and Valuation(2*d0-11,5) ge 3;
  assert (degs eq [1,1,1,1,1]) eq ((not additive) or exceptional);
  assert ((not additive) or exceptional) or degs eq [5];
  printf "d = %o: additive %o, exceptional %o, fifth roots of q(d) of degrees %o\n", d0, additive, exceptional, degs;
end for;
print "PART_V_DONE";

print "ALL_QSPLIT_CONDITIONS_CHECKS_PASSED";
quit;
