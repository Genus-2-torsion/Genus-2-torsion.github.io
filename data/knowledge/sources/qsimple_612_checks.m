// Certificate for the census proof qsimple_612_obstruction.tex (15 November 2026):
// "Rational torsion containing (Z/3)^2 forces 2-primary torsion of order at most 4 on principally
// polarized abelian surfaces that are geometrically split and simple over Q".
// Magma V2.30-1.  Run:  magma -b qsimple_612_checks.m < /dev/null > qsimple_612_checks.log
// Every check is an assertion.  A failed assertion stops the run before the final line
// ALL_612_CHECKS_PASSED (Magma's exit status is 0 even then, so only that line certifies the run).
//   Part A: Theorem 1 (the image of G_Q in GL(A[3]) and its quotients).
//   Part B: Lemma 2 (characteristic polynomials of Frobenius elements of order 3, 4 and 6).
//   Part C: Lemma 2 (i) and (ii) (valuations, by residues).
//   Part D: Theorem 4 (supersingular Weil polynomials and their values at 1).
//   Part E: Theorem 8 (the inertia element).
//   Part F: the Corollary (the three groups).
//   Part V: Remark 3 (validation, not used in the proof).
SetColumns(0); SetMemoryLimit(16*10^9);

// =============================================================================================
// Part A. Theorem 1.
// =============================================================================================
// A1. In a symplectic basis (w1, w2, w1', w2') of A[3] with W = <w1, w2>, the Weil pairing has
// Gram matrix Jf = [[0, I], [-I, 0]].  A matrix g = [[I, X], [0, lam I]] (lam = +-1) is a symplectic
// similitude with multiplier lam, i.e. g^t Jf g = lam Jf, if and only if X is symmetric.  So the
// group U of Theorem 1(b) embeds in the symmetric matrices: #U <= 27.
F3 := GF(3);
I2 := IdentityMatrix(F3, 2); Z2 := ZeroMatrix(F3, 2, 2);
Jf := BlockMatrix(2, 2, [Z2, I2, -I2, Z2]);
Blk := func< X, lam | BlockMatrix(2, 2, [I2, X, Z2, lam*I2]) >;
nsym := 0;
for e in CartesianPower(F3, 4) do
  X := Matrix(F3, 2, 2, [e[1], e[2], e[3], e[4]]);
  for lam in [F3!1, F3!-1] do
    g := Blk(X, lam);
    assert (Transpose(g)*Jf*g eq lam*Jf) eq (Transpose(X) eq X);
  end for;
  if Transpose(X) eq X then nsym +:= 1; end if;
end for;
assert nsym eq 27;
print "A1: similitudes [[I,X],[0,lam I]] <=> X symmetric; 27 symmetric X";

// A2. The group M of all such similitudes (order 54) and its subgroups G that contain an element
// with lam = -1 (the image of G_Q, since the mod-3 cyclotomic character is onto {+-1}).  For each
// conjugacy class of such G, with U = G meet U0 (lam = 1):  #G = 2*#U with #U a power of 3 at most
// 27; every element of G outside U is an involution inverting U; U is the only subgroup of index 2;
// every quotient of G is trivial or generalized dihedral (U' x| C2, C2 acting by -1); and a quotient
// isomorphic to a group in the list of (I3) is C1, C2 or D6.
function IsGenDih3(H)
  if #H eq 1 then return true; end if;
  if #H mod 2 ne 0 then return false; end if;
  m := #H div 2;
  if m ne 1 and not (IsPrimePower(m) and PrimeDivisors(m) eq [3]) then return false; end if;
  Us := [ N`subgroup : N in NormalSubgroups(H) | N`order eq m ];
  for U in Us do
    if IsAbelian(U) and Exponent(U) in {1,3} then
      t := [ x : x in H | x notin U ][1];
      if Order(t) eq 2 and forall{ y : y in Generators(U) | t*y*t^-1 eq y^-1 } then return true; end if;
    end if;
  end for;
  return false;
end function;
GammaList := [ CyclicGroup(1), CyclicGroup(2), CyclicGroup(3), CyclicGroup(4), CyclicGroup(6),
               AbelianGroup(GrpPerm, [2,2]), DihedralGroup(3), DihedralGroup(4), DihedralGroup(6) ];
GammaNames := [ "C1", "C2", "C3", "C4", "C6", "V4", "D6", "D8", "D12" ];
GL4 := GL(4, F3);
Esym := [ Matrix(F3, 2, 2, [1,0,0,0]), Matrix(F3, 2, 2, [0,1,1,0]), Matrix(F3, 2, 2, [0,0,0,1]) ];
M := sub< GL4 | [ GL4!Blk(X, F3!1) : X in Esym ] cat [ GL4!Blk(Z2, F3!-1) ] >;
U0 := sub< M | [ GL4!Blk(X, F3!1) : X in Esym ] >;
assert #M eq 54 and #U0 eq 27 and IsAbelian(U0) and Exponent(U0) eq 3;
assert forall{ g : g in M | Transpose(Matrix(g))*Jf*Matrix(g) eq Matrix(g)[3,3]*Jf };
fM, Mp := CosetAction(M, sub< M | >);
assert #Mp eq 54;
U0p := fM(U0);
nG := 0; seen := {};
for S in Subgroups(Mp) do
  G := S`subgroup;
  if G subset U0p then continue; end if;
  nG +:= 1;
  U := G meet U0p;
  assert #G eq 2*#U and #U in {1, 3, 9, 27};
  assert forall{ t : t in G | t in U or (Order(t) eq 2 and forall{ u : u in Generators(U) | t*u*t^-1 eq u^-1 }) };
  idx2 := [ H`subgroup : H in Subgroups(G : OrderEqual := #G div 2) ];
  assert #idx2 eq 1 and idx2[1] eq U;
  for N in NormalSubgroups(G) do
    Q := quo< G | N`subgroup >;
    assert IsGenDih3(Q);
    for i in [1..#GammaList] do
      if IsIsomorphic(Q, GammaList[i]) then Include(~seen, GammaNames[i]); end if;
    end for;
  end for;
end for;
assert nG eq 28;
assert seen eq {"C1", "C2", "D6"};
printf "A2: %o classes of subgroups G of M with lam(G) = {+-1}; quotients in the list of (I3): %o\n", nG,
       Sort(SetToSequence(seen));

// A3. The abstract groups U x| C2 (U elementary abelian of order 3^k, k <= 3, C2 acting by -1):
// #G = 2*3^k, every quotient is trivial or of the same form, and the members of the list of (I3) of
// this form are exactly C1, C2 and D6.
for k in [0..3] do
  if k eq 0 then G := CyclicGroup(2); else
    U := AbelianGroup([3 : i in [1..k]]);
    Ap := AutomorphismGroup(U);
    inv := hom< U -> U | [ -U.i : i in [1..k] ] >;
    C2g := CyclicGroup(GrpAb, 2);
    G := SemidirectProduct(U, C2g, hom< C2g -> Ap | [ Ap!inv ] >);
  end if;
  G := PCGroup(G);
  assert #G eq 2*3^k;
  for N in NormalSubgroups(G) do
    assert IsGenDih3(quo< G | N`subgroup >);
  end for;
end for;
ok := [ GammaNames[i] : i in [1..#GammaList] | IsGenDih3(GammaList[i]) ];
assert ok eq ["C1", "C2", "D6"];
print "A3: quotients of U x| C2 (#U <= 27) are generalized dihedral; list members of this form:", ok;
print "PART_A_OK";

// =============================================================================================
// Part B. Lemma 2: N(X^2 - tX + p eps) for the three orders, with zeta = zeta_3 and i = sqrt(-1).
// =============================================================================================
Qz<z> := CyclotomicField(3);
R<a,b,p,u> := PolynomialRing(Qz, 4);
S<X> := PolynomialRing(R);
conjQ := hom< Qz -> Qz | z^2 >;
function ConjR(r)
  cs, ms := CoefficientsAndMonomials(r);
  return &+[ R | conjQ(cs[i])*ms[i] : i in [1..#cs] ];
end function;
function ConjS(f)
  return S![ ConjR(c) : c in Coefficients(f) ];
end function;
// order 3: eps = zeta, t = a zeta^2
f3 := X^2 - a*z^2*X + p*z;
P3 := f3*ConjS(f3);
assert P3 eq X^4 + a*X^3 + (a^2-p)*X^2 + a*p*X + p^2;
assert Evaluate(f3, 1) eq (1+a) + (a+p)*z;
assert Evaluate(P3, 1) eq (1+a)^2 - (1+a)*(a+p) + (a+p)^2;     // = N((1+a)+(a+p)zeta)
// order 6: eps = -zeta^2 (a primitive 6th root of unity), t = b(-2-zeta) = b zeta sqrt(-3)
assert (-z^2)^6 eq 1 and (-z^2)^2 ne 1 and (-z^2)^3 ne 1;
f6 := X^2 - b*(-2-z)*X + p*(-z^2);
P6 := f6*ConjS(f6);
assert P6 eq X^4 + 3*b*X^3 + (3*b^2+p)*X^2 + 3*b*p*X + p^2;
assert Evaluate(f6, 1) eq (1+2*b+p) + (b+p)*z;
assert Evaluate(P6, 1) eq (1+2*b+p)^2 - (1+2*b+p)*(b+p) + (b+p)^2;
// the other primitive 6th root of unity, eps = -zeta, t = b(1-zeta): the same shape with b -> -b
f6b := X^2 - b*(1-z)*X + p*(-z);
P6b := f6b*ConjS(f6b);
assert P6b eq X^4 - 3*b*X^3 + (3*b^2+p)*X^2 - 3*b*p*X + p^2;
// square roots used for gamma_g: zeta = (zeta^2)^2, -zeta^2 = (zeta sqrt(-3))^2/3, sqrt(-3) = 1 + 2 zeta
s3 := 1 + 2*z; assert s3^2 eq -3;
assert z eq (z^2)^2;
assert -z^2 eq (z*s3)^2/3;  assert z*s3 eq -2 - z;
assert -z eq (z^2*s3)^2/3;  assert z^2*s3 eq 1 - z;
// order 4: eps = i, t = u(1+i); eps = -i gives the conjugate factor and the same polynomial
Qi<ii> := CyclotomicField(4);
Ri<pp,uu> := PolynomialRing(Qi, 2);
Si<Y> := PolynomialRing(Ri);
conjI := hom< Qi -> Qi | -ii >;
function ConjRi(r)
  cs, ms := CoefficientsAndMonomials(r);
  return &+[ Ri | conjI(cs[i])*ms[i] : i in [1..#cs] ];
end function;
function ConjSi(f)
  return Si![ ConjRi(c) : c in Coefficients(f) ];
end function;
assert ii^2 eq -1;
f4 := Y^2 - uu*(1+ii)*Y + ii*pp;
P4 := f4*ConjSi(f4);
assert P4 eq Y^4 - 2*uu*Y^3 + 2*uu^2*Y^2 - 2*pp*uu*Y + pp^2;
assert ConjSi(f4) eq Y^2 - uu*(1-ii)*Y - ii*pp;
assert Evaluate(f4, 1) eq (1-uu) + ii*(pp-uu);
assert Evaluate(P4, 1) eq (1-uu)^2 + (pp-uu)^2;
// square roots: i = (1+i)^2/2, -i = (1-i)^2/2
assert ii eq (1+ii)^2/2 and -ii eq (1-ii)^2/2;
print "PART_B_OK";

// =============================================================================================
// Part C. Lemma 2 (i), (ii).  N3(a,p) = N((1+a)+(a+p)zeta), N6(b,p) = N((1+2b+p)+(b+p)zeta),
// N4(u,p) = (1-u)^2 + (p-u)^2; the residues mod 9 (resp. 16) depend only on the residues of the
// arguments mod 9 (resp. 16).
// =============================================================================================
N3 := func< A, P | (1+A)^2 - (1+A)*(A+P) + (A+P)^2 >;
N6 := func< B, P | (1+2*B+P)^2 - (1+2*B+P)*(B+P) + (B+P)^2 >;
N4 := func< U, P | (1-U)^2 + (P-U)^2 >;
// 3-adic, p = 2 mod 3: v_3(N3) <= 1, 3 does not divide N6 or N4
for A in [0..8] do for P in [0..8] do
  if P mod 3 eq 2 then
    assert N3(A,P) mod 9 ne 0;
    assert N6(A,P) mod 3 ne 0;
    assert N4(A,P) mod 3 ne 0;
  end if;
end for; end for;
// 2-adic, p = 3 mod 4: v_2(N3), v_2(N6) in {0, 2}, and 8 does not divide N4
for A in [0..15] do for P in [0..15] do
  if P mod 4 eq 3 then
    for N in [N3(A,P), N6(A,P)] do
      assert (N mod 2 eq 1) or (N mod 8 eq 4);
    end for;
    assert N4(A,P) mod 8 ne 0;
  end if;
end for; end for;
// the hypotheses on p are needed
assert exists{ <A,P> : A in [0..8], P in [0..8] | P mod 3 eq 1 and N3(A,P) mod 9 eq 0 };
assert exists{ <A,P> : A in [0..15], P in [0..15] | P mod 4 eq 1 and N3(A,P) mod 16 eq 0 };
print "PART_C_OK";

// =============================================================================================
// Part D. Theorem 4.
// D1. Supersingular Weil polynomials of abelian surfaces over F_p, p <= 23:
// P = X^4 + c1 X^3 + c2 X^2 + p c1 X + p^2 with all complex roots of absolute value sqrt(p) iff
// h(y) = y^2 + c1 y + (c2 - 2p) has two real roots in [-2 sqrt p, 2 sqrt p]; supersingular means
// p | c1 and p | c2.
// =============================================================================================
function IsWeil(pp, c1, c2)
  c := c2 - 2*pp; D := c1^2 - 4*c;
  if D lt 0 then return false; end if;
  if c1^2 gt 16*pp then return false; end if;
  s := 4*pp + c;
  if s lt 0 then return false; end if;
  return s^2 ge 4*pp*c1^2;
end function;
SSP1 := func< P | [ (P+1)^2, P^2+P+1, P^2+1, P^2-P+1, (P-1)^2 ] >;   // the values 1 + kp + p^2, |k| <= 2
SSVALS := AssociativeArray();
VALS := AssociativeArray();
for pp in [2,3,5,7,11,13,17,19,23] do
  ss := [];
  B := Isqrt(16*pp);
  for c1 in [-B..B] do for c2 in [-6*pp-1..6*pp+1] do
    if IsWeil(pp, c1, c2) and c1 mod pp eq 0 and c2 mod pp eq 0 then Append(~ss, <c1,c2>); end if;
  end for; end for;
  SSVALS[pp] := ss;
  VALS[pp] := { 1 + c[1] + c[2] + pp*c[1] + pp^2 : c in ss };
  printf "p=%o supersingular (c1,c2): %o ; P(1) values %o\n", pp, ss, VALS[pp];
  if pp ge 7 then
    assert Set(ss) eq { <0, k*pp> : k in [-2..2] };
    assert VALS[pp] eq Set(SSP1(pp));                            // D2: the list of Theorem 4
  end if;
end for;
assert Set(SSVALS[5]) eq { <0, k*5> : k in [-2..2] } join { <5,15>, <-5,15> };
assert VALS[5] eq Set(SSP1(5)) join {11, 71};
assert VALS[5] eq {11, 16, 21, 26, 31, 36, 71} and VALS[7] eq {36, 43, 50, 57, 64};
// (for p >= 17, p | c1 and |c1| <= 4 sqrt(p) < p force c1 = 0 and then c2 = kp, |k| <= 2)
print "D1, D2: supersingular Weil polynomials and their values at 1 for p <= 23";
// D3. 72 never divides a supersingular value at 1 when p = 5, 7 mod 12 (residues mod 72*12 suffice
// for the five values 1 + kp + p^2), nor for p = 5 (all seven values).
for r in [0..72*12-1] do
  if r mod 12 in {5,7} then
    for v in SSP1(r) do assert v mod 72 ne 0; end for;
  end if;
end for;
assert forall{ v : v in VALS[5] | v mod 72 ne 0 };
// D4. For p prime to 6: p^2 + 1 is prime to 3; p^2 +- p + 1 are odd and not divisible by 9.  So
// 9 | P(1) or 12 | P(1) forces P(1) = (p +- 1)^2.
for r in [0..36*9-1] do
  if GCD(r, 6) eq 1 then
    v := SSP1(r);
    assert v[3] mod 3 ne 0;
    assert v[2] mod 9 ne 0 and v[4] mod 9 ne 0;
    assert v[2] mod 2 eq 1 and v[4] mod 2 eq 1;
  end if;
end for;
// D5. With pi^2 = -p (resp. p), the group of F_p-points is killed by p + 1 (resp. p - 1): no point of
// order 12 if p = 5, 7 mod 12, and no point of order 9 if p is not +-1 mod 9.
for r in [0..35] do if r mod 12 in {5,7} then
  assert ((r mod 3 eq 2) and (r+1) mod 12 ne 0) or ((r mod 3 eq 1) and (r-1) mod 12 ne 0);
end if; end for;
for r in [0..8] do if GCD(r,3) eq 1 and r mod 9 notin {1,8} then
  assert (r+1) mod 9 ne 0 and (r-1) mod 9 ne 0;
end if; end for;
// D6. p = 5 (inert in Q(sqrt(-3))) and p = 7 (inert in Q(i)): 9 | P(1) or 12 | P(1) only for P(1) = 36,
// which is (p+1)^2 for p = 5 and (p-1)^2 for p = 7; and 5, 7 are not +-1 mod 9, and are 5, 7 mod 12.
assert { v : v in VALS[5] | v mod 9 eq 0 or v mod 12 eq 0 } eq {36} and (5+1)^2 eq 36;
assert { v : v in VALS[7] | v mod 9 eq 0 or v mod 12 eq 0 } eq {36} and (7-1)^2 eq 36;
assert KroneckerSymbol(-3, 5) eq -1 and KroneckerSymbol(-4, 7) eq -1;
assert 5 mod 9 notin {1,8} and 7 mod 9 notin {1,8} and 5 mod 12 eq 5 and 7 mod 12 eq 7;
print "PART_D_OK";

// =============================================================================================
// Part E. Theorem 8.
// =============================================================================================
// E1. F = Q(s), s^2 = a.  For t = y s (an element with gamma(t) = -t):
// N_{F/Q}(X^2 - tX + 1) = (X^2 - tX + 1)(X^2 + tX + 1) = (X^2+1)^2 - a y^2 X^2.
Ra<aa,yy> := PolynomialRing(Rationals(), 2);
Sa<ss> := PolynomialRing(Ra);
Qa<sb> := quo< Sa | ss^2 - aa >;
Ta<XX> := PolynomialRing(Qa);
lhs := (XX^2 - yy*sb*XX + 1)*(XX^2 + yy*sb*XX + 1);
rhs := (XX^2+1)^2 - aa*yy^2*XX^2;
assert lhs eq Ta!rhs;
print "E1: norm identity";
// E2. P_c = X^4 + (2-c)X^2 + 1 = (X^2+1)^2 - cX^2.  For c in [-200, 200], P_c is a product of
// cyclotomic polynomials iff c in {0,1,2,3,4}; the factorizations for c = 0..4; and c in a Z^2 with
// a > 1 squarefree and 0 < c <= 4 forces (a,c) in {(2,2), (3,3)}.
Z := Integers(); Px<x> := PolynomialRing(Z);
function IsCycloProduct(f)
  fac := Factorization(f);
  for t in fac do
    if not exists{ n : n in [1..60] | t[1] eq CyclotomicPolynomial(n) } then return false; end if;
  end for;
  return true;
end function;
good := [ c : c in [-200..200] | IsCycloProduct(x^4 + (2-c)*x^2 + 1) ];
assert good eq [0,1,2,3,4];
assert x^4+2*x^2+1 eq CyclotomicPolynomial(4)^2;
assert x^4+x^2+1 eq CyclotomicPolynomial(3)*CyclotomicPolynomial(6);
assert x^4+1 eq CyclotomicPolynomial(8);
assert x^4-x^2+1 eq CyclotomicPolynomial(12);
assert x^4-2*x^2+1 eq CyclotomicPolynomial(1)^2*CyclotomicPolynomial(2)^2;
pairs := [ <A,c> : A in [2..200], c in [1..4] | IsSquarefree(A) and c mod A eq 0 and IsSquare(c div A) ];
assert pairs eq [ <2,2>, <3,3> ];
print "E2: cyclotomic analysis";
// E3. Powers: X^2 = -1 mod Phi_4, X^4 = -1 mod Phi_8, X^6 = -1 mod Phi_12, so for semisimple rho(g) with
// characteristic polynomial Phi_4^2, Phi_8, Phi_12 the element h = g, g^2, g^3 has rho(h)^2 = -1.
assert (x^2+1) mod CyclotomicPolynomial(4) eq 0;
assert (x^4+1) mod CyclotomicPolynomial(8) eq 0;
assert (x^6+1) mod CyclotomicPolynomial(12) eq 0;
for pr in [<4,2>, <8,4>, <12,6>] do
  Mc := CompanionMatrix(CyclotomicPolynomial(pr[1]));
  assert Mc^pr[2] eq -IdentityMatrix(Z, Nrows(Mc));
end for;
print "E3: power identities";
// E4. T_2 is free of rank 2 over Z_2[i] (text), so (h acting on A[2^infinity]) is the standard
// model i = [[0,-1],[1,0]] + [[0,-1],[1,0]]; modulo 2^k (k <= 12) the kernel of i - 1 on (Z/2^k)^4 has
// order 4 and is killed by 2.
Jm := Matrix(Z, 2, 2, [0,-1,1,0]);
I4 := DiagonalJoin(Jm, Jm);
assert I4^2 eq -IdentityMatrix(Z, 4);
for k in [1..12] do
  G := AbelianGroup([2^k, 2^k, 2^k, 2^k]);
  Mk := I4 - IdentityMatrix(Z, 4);
  f := hom< G -> G | [ G!Eltseq(ColumnSubmatrix(Mk, j, 1)) : j in [1..4] ] >;
  Kf := Kernel(f);
  assert #Kf eq 4;
  assert &and[ 2*g eq Kf!0 : g in Kf ];
end for;
print "E4: fixed points of h on A[2^infinity]";
// E5. (h - 1)(h + 1) = h^2 - 1 = -2, so h - 1 is invertible on T_l for odd l.
assert (I4 - IdentityMatrix(Z, 4))*(I4 + IdentityMatrix(Z, 4)) eq -2*IdentityMatrix(Z, 4);
print "E5: (h-1)(h+1) = -2";
print "PART_E_OK";

// =============================================================================================
// Part F. The Corollary: [6,12], [3,24] and [2,6,6] have order 72, contain (Z/3)^2 (their 3-torsion
// subgroup has order 9) and have 2-primary part of order 8.
// =============================================================================================
for inv in [ [6,12], [3,24], [2,6,6] ] do
  Ab := AbelianGroup(inv);
  assert #Ab eq 72;
  assert #{ g : g in Ab | 3*g eq Ab!0 } eq 9;
  assert #{ g : g in Ab | 8*g eq Ab!0 } eq 8;
end for;
print "PART_F_OK";

// =============================================================================================
// Part V. Validation, not used in the proof (Remark 3).
// =============================================================================================
// V1. The census curve 2.2.2.2.d (LMFDB 1521.b1): y^2 + (x^2+x) y = 3x^5 + 16x^4 + 17x^3 + x^2 - 3x.
// Its discriminant is 3^6 13^6, so it has good reduction at every prime p not in {3, 13}; its
// rational torsion subgroup is (Z/2)^4.  At the primes 5 <= p < 200, p != 13: if p splits in Q(sqrt 13)
// the L-polynomial is a square; if p is inert, it is (1 - pT^2)^2 + 3 b_p^2 T^2 with b_p >= 0, a square
// iff b_p = 0.  b_7 = b_11 = b_19 = 2 and b_37 = 4; in particular the L-polynomials at 7 and 11 are not
// squares.
Pq<t> := PolynomialRing(Rationals());
CV := HyperellipticCurve(3*t^5 + 16*t^4 + 17*t^3 + t^2 - 3*t, t^2 + t);
assert Discriminant(CV) eq 3^6*13^6;
assert Invariants(TorsionSubgroup(Jacobian(SimplifiedModel(CV)))) eq [2,2,2,2];
PT<T> := PolynomialRing(Z);
bvals := AssociativeArray();
nsplit := 0; ninert := 0;
for q in PrimesInInterval(5, 199) do
  if q eq 13 then continue; end if;
  Lq := PT!LPolynomial(ChangeRing(CV, GF(q)));
  issq := IsSquare(Lq);
  if KroneckerSymbol(13, q) eq 1 then
    nsplit +:= 1;
    assert issq;
  else
    ninert +:= 1;
    d := Coefficient(Lq, 2) + 2*q;
    assert Coefficient(Lq, 1) eq 0 and Coefficient(Lq, 3) eq 0 and Coefficient(Lq, 4) eq q^2;
    assert d ge 0 and d mod 3 eq 0 and IsSquare(d div 3);
    bq := Isqrt(d div 3); bvals[q] := bq;
    assert Lq eq (1 - q*T^2)^2 + 3*bq^2*T^2;
    assert issq eq (bq eq 0);
  end if;
end for;
assert [ bvals[q] : q in [7, 11, 19, 37] ] eq [2, 2, 2, 4];
assert KroneckerSymbol(-3, 7) eq 1 and KroneckerSymbol(-39, 11) eq 1;   // 7 splits in Q(sqrt -3), 11 in Q(sqrt -39)
printf "V1: 2.2.2.2.d: %o split and %o inert primes checked; b_p at 7, 11, 19, 37: %o\n", nsplit, ninert,
       [ bvals[q] : q in [7, 11, 19, 37] ];
// V2. The LMFDB curves 3969.d.250047.1 and 25600.f.512000.1: rational torsion subgroups [2,6] and [2,2,2].
C1 := HyperellipticCurve(t - 4*t^3 + 5*t^4 - 3*t^5, 1 + t + t^2);
assert Invariants(TorsionSubgroup(Jacobian(SimplifiedModel(C1)))) eq [2,6];
C2c := HyperellipticCurve(2*t + 5*t^2 - t^3 - 5*t^4 + 2*t^5);
assert Invariants(TorsionSubgroup(Jacobian(C2c))) eq [2,2,2];
print "V2: torsion of 3969.d.250047.1 and 25600.f.512000.1";
print "PART_V_OK";
print "ALL_612_CHECKS_PASSED";
quit;
