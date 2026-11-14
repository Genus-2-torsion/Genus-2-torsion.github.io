/* family_40_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [40].

   The family (derived 2026-11-14, record torsion_jac/notes/simple40_2026_11_14.md): a relation-lattice locus in the
   sense of Leprevost and Elkies, with a rational Weierstrass point W at infinity, two rational non-Weierstrass
   points P1, P2 and a rational Weierstrass point E, and the relations
       f = A^2 + lam x^4 (x - e),  A(e) = 0:   div(y - A)  = 4 P1 + E - 5W,   P1 = (0, A(0)), E = (e, 0),
       f = A'^2 + lam (x - 1)^5:               div(y - A') = 5 P2 - 5W,       P2 = (1, A'(1)).
   Subtracting, (A - A')(A + A') = lam R4 with R4 = (x-1)^5 - x^4 (x-e), so R4 must split as 2 + 2.  The curve of
   such splittings has genus 0; with the parameter k:
       al0 = (2k+5)/(k^3+5k^2+10k+10),  al1 = k al0,  kap = (10 al0^2 - 5 al0 - al1)/(al0^2 al1),  e = 5 + kap,
       g1 = x^2 + al1 x + al0,  g2 = x^2 + (10/kap - al1) x - 1/(kap al0),   R4 = kap g1 g2,
       lam = -g1(e)/(kap g2(e)),  A = (g1 + lam kap g2)/2,  A' = (lam kap g2 - g1)/2,  f_k = A^2 + lam x^4 (x - e).

   Argument.
   (1) Over Q(k): the two expressions for f agree, A(e) = 0, D1 = [x, A(0)] satisfies 4 D1 = [x - e, 0] != 0
       (so D1 has order 8), D2 = [x - 1, A'(1)] has order 5, and D1 + D2 has order 40 (Mumford arithmetic over Q(k)).
       So <D1, D2> = Z/40 is marked over Q(k).
   (2) Residue certificate (family_cert_lib.m) at k = 6 mod 7 and k = 16 mod 17, i.e. k = -1 mod 119:
       #J(F_7) = 40 and #J(F_17) = 320, gcd 40 = #G, and the Frobenius at 17 is strict.  So every 7- and
       17-integral k = -1 mod 119 gives J(Q)_tors exactly [40] and a geometrically simple Jacobian.
   (3) The G2-invariants of two members differ, so infinitely many isomorphism classes.
   The family does not contain census curve 40.a (Elkies 2002): 40.a has no rational point P with P - W of order 8
   for either of its rational Weierstrass points W.
   Run (from the census root):
       magma -b pipeline/magma/family_40_certificate.m > data/knowledge/sources/family_40_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// (1) the family and the marked group over Q(k)
Fk<k> := FunctionField(Q); P<x> := PolynomialRing(Fk);
al0 := (2*k + 5)/(k^3 + 5*k^2 + 10*k + 10);
al1 := k*al0;
kap := (10*al0^2 - 5*al0 - al1)/(al0^2*al1);
e := 5 + kap;
g1 := x^2 + al1*x + al0;
g2 := x^2 + (10/kap - al1)*x - 1/(kap*al0);
assert (x-1)^5 - x^4*(x - e) eq kap*g1*g2;
lam := -Evaluate(g1, e)/(kap*Evaluate(g2, e));
A  := (g1 + lam*kap*g2)/2;
A1 := (lam*kap*g2 - g1)/2;
f := A^2 + lam*x^4*(x - e);
assert f eq A1^2 + lam*(x-1)^5;
assert Evaluate(A, e) eq 0 and Degree(f) eq 5;
print "over Q(k): f = A^2 + lam x^4 (x - e) = A'^2 + lam (x - 1)^5 with A(e) = 0";
J := Jacobian(HyperellipticCurve(f));
D1 := elt< J | x, Evaluate(A, 0) >;
D2 := elt< J | x - 1, Evaluate(A1, 1) >;
DE := elt< J | x - e, 0 >;
assert DE ne J!0 and 2*DE eq J!0 and 4*D1 eq DE;
assert D2 ne J!0 and 5*D2 eq J!0;
D40 := D1 + D2;
assert 20*D40 ne J!0 and 8*D40 ne J!0 and 40*D40 eq J!0;
print "over Q(k): D1 = [x, A(0)] has order 8 (4 D1 = [x - e, 0]), D2 = [x - 1, A'(1)] has order 5, D1 + D2 has order 40";

// (2), (3) residue certificate at k = 6 mod 7, k = 16 mod 17
g, st := ResidueCertificate(f, [7, 17], [6, 16], [-1, 118, -120, 117/2], 40, [40]);
print "FAMILY_40_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [40]";
quit;
