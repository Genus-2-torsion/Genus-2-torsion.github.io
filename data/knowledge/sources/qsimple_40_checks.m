// Certificate for the census proof "Rational 5-torsion and the 2-primary torsion of geometrically split
// abelian surfaces" (proofs/qsimple-40-obstruction.html): part K is the input (I10) used in Lemma 4';
// parts S and J check the sharpness examples of Remark 4.  Run: magma -b qsimple_40_checks.m < /dev/null
// Every claim is an assert; the run ends with ALL_T5SQ_CHECKS_PASSED.
// The reused computational inputs (Lemma 1, Lemma 7 and Proposition 8 of the [2,2,10] proof) are certified
// by the published qsimple_2210_checks.m.
SetColumns(0);
SetMemoryLimit(16*10^9);
Q := Rationals();
PZ<T> := PolynomialRing(Integers());

// ---------------------------------------------------------------------------
// Part K.  Input (I10): Kenku-Momose, Theorem (0.1).  X_1(20) = X_1(1,20) has positive genus,
// and (1,20) is not among the exceptions (1,11),(1,13),(1,14),(1,15),(1,16),(1,18),(2,10),(2,12).
// Hence Y_1(20)(k) is empty for every quadratic field k.  (X_1(2,10), an exception, has genus 1.)
g20 := Genus(Gamma1(20));
assert g20 eq 3;
exceptions := { [1,11], [1,13], [1,14], [1,15], [1,16], [1,18], [2,10], [2,12] };
assert not [1,20] in exceptions;
assert [2,10] in exceptions;
printf "K: genus X_1(20) = %o > 0 and (1,20) is not a Kenku-Momose exception: no point of order 20 over quadratic fields\n", g20;

// ---------------------------------------------------------------------------
// Part S.  Sharpness.  E = LMFDB 2.2.13.1-36.1-a3 over K = Q(sqrt 13) = Q(w), w^2 - w - 3 = 0.
// Claims: E(K)_tors = Z/2 x Z/10; E has no CM (j not integral); E is K-isomorphic to E^sigma;
// 7 is inert in K, E has good reduction at 7 O_K, and T^4 - a T^2 + 49 (a = a_{7O_K}(E)), the
// Frobenius polynomial at 7 of A = Res_{K/Q} E, is irreducible over Q, so A is simple over Q.
// Then A is of type (S) with Gamma = C_2, principally polarized, and A(Q)_tors = E(K)_tors = Z/2 x Z/10.
K<w> := NumberField(T^2 - T - 3);
OK := MaximalOrder(K);
assert Discriminant(OK) eq 13;
E := EllipticCurve([K | 1 + w, 0, 1 + w, -115 + 48*w, 218 - 96*w]);
TE := TorsionSubgroup(E);
assert Invariants(TE) eq [2, 10];
j := jInvariant(E);
assert j eq 476379541/236196;                 // rational, not integral: E has no CM
assert not IsIntegral(Q!j);
sig := [ a : a in Automorphisms(K) | a(w) ne w ][1];
Es := EllipticCurve([ sig(a) : a in aInvariants(E) ]);
assert IsIsomorphic(E, Es);
fac := Factorization(7*OK);
assert #fac eq 1 and fac[1][2] eq 1;          // 7 is inert in K
P7 := fac[1][1];
assert Valuation(Discriminant(E), P7) eq 0;   // good reduction at 7 O_K (integral model)
assert forall{ a : a in aInvariants(E) | IsIntegral(a) };
a7 := TraceOfFrobenius(Reduction(E, P7));
assert a7 eq 10;
chi7 := T^4 - a7*T^2 + 49;
assert IsIrreducible(chi7);
assert Evaluate(chi7, 1) eq 40;               // #A(F_7) = #E(F_49) = 40
assert #Reduction(E, P7) eq 40;
printf "S: E = 2.2.13.1-36.1-a3: E(K)_tors = Z/2 x Z/10, j = %o (no CM), E ~= E^sigma over K, chi_7(Res E) = %o irreducible\n", j, chi7;
// An inert prime with 8 not dividing #E(F_{p^2}) (consistent with Proposition 8 of the [2,2,10] proof):
fac11 := Factorization(11*OK);
assert #fac11 eq 1 and fac11[1][2] eq 1;
assert Valuation(Discriminant(E), fac11[1][1]) eq 0;
assert #Reduction(E, fac11[1][1]) eq 100;     // 100 = 4 mod 8
print "S: #E(F_121) = 100, not divisible by 8";
// The second curve of Remark 3: LMFDB 2.2.17.1-36.1-a4 over Q(sqrt 17) = Q(v), v^2 - v - 4 = 0, with torsion
// Z/2 x Z/10; 11 is inert in Q(sqrt 17) and #E(F_121) = 100.
K2<v> := NumberField(T^2 - T - 4);
OK2 := MaximalOrder(K2);
assert Discriminant(OK2) eq 17;
E2 := EllipticCurve([K2 | 1, -v, v, -270 - 174*v, 2840 + 1819*v]);
assert Invariants(TorsionSubgroup(E2)) eq [2, 10];
assert forall{ a : a in aInvariants(E2) | IsIntegral(a) };
fac11b := Factorization(11*OK2);
assert #fac11b eq 1 and fac11b[1][2] eq 1;    // 11 is inert in Q(sqrt 17)
assert Valuation(Discriminant(E2), fac11b[1][1]) eq 0;
assert #Reduction(E2, fac11b[1][1]) eq 100;
print "S: E = 2.2.17.1-36.1-a4 over Q(sqrt 17): E(K)_tors = Z/2 x Z/10, 11 inert, #E(F_121) = 100, not divisible by 8";

// ---------------------------------------------------------------------------
// Part J.  Sharpness by a Jacobian.  C : y^2 = 63x^6 + 6x^5 + 184x^4 + 184x^2 - 6x + 63 over Q.
// Claims: J(Q)_tors = Z/2 x Z/10; J is simple over Q; over K = Q(sqrt 13), C is isomorphic to the
// Howe-Leprevost-Poonen curve y^2 = h(X) glued from E (Part S) and E^sigma along psi((al_i,0)) = (sigma(al_i),0),
// h is even, and its two elliptic quotients (X,y) -> (X^2,y), (X,y) -> (1/X^2, y/X^3) are K-isomorphic
// to E and to E^sigma, which are K-isomorphic and have no CM.  Hence J_K ~ E^2, all geometric
// endomorphisms are defined over K, and J is of type (S) with Gamma = C_2 (Gamma != 1 as J is Q-simple).
PQ<x> := PolynomialRing(Q);
fC := 63*x^6 + 6*x^5 + 184*x^4 + 184*x^2 - 6*x + 63;
C := HyperellipticCurve(fC);
TJ := TorsionSubgroup(Jacobian(C));
assert Invariants(TJ) eq [2, 10];
assert Factorization(Integers()!Discriminant(C)) eq [<2,26>, <3,10>, <13,6>];
L7 := LPolynomial(ChangeRing(C, GF(7)));
assert L7 eq Parent(L7)![1, 0, -10, 0, 49];
assert IsIrreducible(L7);
printf "J: C : y^2 = %o; J(Q)_tors = Z/2 x Z/10; L_7(T) = %o irreducible, so J is simple over Q\n", fC, L7;
// the HLP sextic from E over K
function HLPSextic(al, be)
    KK := Parent(al[1]); PKK<X> := PolynomialRing(KK);
    a1 := (al[3]-al[2])^2/(be[3]-be[2]) + (al[2]-al[1])^2/(be[2]-be[1]) + (al[1]-al[3])^2/(be[1]-be[3]);
    b1 := (be[3]-be[2])^2/(al[3]-al[2]) + (be[2]-be[1])^2/(al[2]-al[1]) + (be[1]-be[3])^2/(al[1]-al[3]);
    a2 := al[1]*(be[3]-be[2]) + al[2]*(be[1]-be[3]) + al[3]*(be[2]-be[1]);
    b2 := be[1]*(al[3]-al[2]) + be[2]*(al[1]-al[3]) + be[3]*(al[2]-al[1]);
    assert a1 ne 0 and a2 ne 0 and b1 ne 0 and b2 ne 0;
    Df := ((al[1]-al[2])*(al[1]-al[3])*(al[2]-al[3]))^2;
    Dg := ((be[1]-be[2])*(be[1]-be[3])*(be[2]-be[3]))^2;
    A := Dg*a1/a2; B := Df*b1/b2;
    return -(A*(al[2]-al[1])*(al[1]-al[3])*X^2 + B*(be[2]-be[1])*(be[1]-be[3]))
           *(A*(al[3]-al[2])*(al[2]-al[1])*X^2 + B*(be[3]-be[2])*(be[2]-be[1]))
           *(A*(al[1]-al[3])*(al[3]-al[2])*X^2 + B*(be[1]-be[3])*(be[3]-be[2]));
end function;
PK<X> := PolynomialRing(K);
ea1, ea2, ea3, ea4, ea6 := Explode(aInvariants(E));
cub := X^3 + (ea2 + ea1^2/4)*X^2 + (ea4 + ea1*ea3/2)*X + (ea6 + ea3^2/4);
al := [ r[1] : r in Roots(cub) ];
assert #al eq 3;
h := HLPSextic(al, [ sig(r) : r in al ]);
assert forall{ i : i in [1,3,5] | Coefficient(h, i) eq 0 };          // h is even
CK := HyperellipticCurve(h);
assert IsIsomorphic(BaseChange(C, K), CK);
ht := PK![ Coefficient(h, 2*i) : i in [0..3] ];                        // h(X) = ht(X^2)
htr := PK![ Coefficient(ht, 3-i) : i in [0..3] ];                      // U^3 ht(1/U)
function cubic_curve(g)     // the elliptic curve V^2 = g(U), g a cubic over K
  c3 := LeadingCoefficient(g);
  return EllipticCurve([0, Coefficient(g,2), 0, Coefficient(g,1)*c3, Coefficient(g,0)*c3^2]);
end function;
E1 := cubic_curve(ht); E2 := cubic_curve(htr);
assert IsIsomorphic(E1, E) or IsIsomorphic(E1, Es);
assert IsIsomorphic(E2, E) or IsIsomorphic(E2, Es);
assert IsIsomorphic(E1, E2);
print "J: over K, C ~= y^2 = h(X) (HLP gluing of E and E^sigma), h even; both quotients are K-isomorphic to E ~= E^sigma (no CM): type (S), Gamma = C_2";
print "ALL_T5SQ_CHECKS_PASSED";
quit;
