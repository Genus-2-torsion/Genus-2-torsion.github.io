/* family_33_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [3,3].

   The family: Bruin-Flynn-Testa, Descent via (3,3)-isogeny on Jacobians of genus 2 curves,
   Acta Arith. 165 (2014), Theorem 6 (arXiv:1401.0580): with parameters r, s, t,
       H1 = x^2 + r x + t,  lambda1 = 4s,     H2 = x^2 + x + r,  lambda2 = 4st,
       G1 = (s - st - 1)x^3 + 3s(r-t)x^2 + 3sr(r-t)x - st^2 + sr^3 + t,
       G2 = (s - st + 1)x^3 + 3s(r-t)x^2 + 3sr(r-t)x - st^2 + sr^3 - t,
       F = G1^2 + lambda1 H1^3 = G2^2 + lambda2 H2^3.
   (Formulas as stated in the paper, transcribed 2026-10-18; the identity is checked below over Q(r,s,t).)
   One-parameter slice used here: r = 2, s = 3, t free.

   Argument.
   (1) Over Q(r,s,t) the two expressions for F agree.  On the slice, div(y - Gi) = 3 Di - 3(oo+ + oo-)
       with Di = {Hi = 0, y = Gi}, so Ti = [Hi, Gi mod Hi] has order 3 (checked: Ti != 0, 3 Ti = 0), and
       T2 != T1, T2 != -T1, so <T1, T2> = (Z/3)^2 is marked over Q(t).
   (2) Residue certificate (family_cert_lib.m) at t = 4 mod 5 and t = 6 mod 13 (t = 19 mod 65):
       #J(F_5) = 9, #J(F_13) = 198, gcd 9 = #G, strict at 13.  So every 5- and 13-integral t = 19 mod 65
       gives J(Q)_tors exactly [3,3] and a geometrically simple Jacobian.
   (3) The G2-invariants of two members differ, so infinitely many isomorphism classes.
   Run (from the census root):
       magma -b pipeline/magma/family_33_certificate.m > data/knowledge/sources/family_33_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// the BFT identity over Q(r,s,t)
F3<r3,s3,t3> := FunctionField(Q, 3); P3<x3> := PolynomialRing(F3);
BFT := func< r, s, tt, xx | [ xx^2 + r*xx + tt, 4*s, xx^2 + xx + r, 4*s*tt,
    (s - s*tt - 1)*xx^3 + 3*s*(r-tt)*xx^2 + 3*s*r*(r-tt)*xx - s*tt^2 + s*r^3 + tt,
    (s - s*tt + 1)*xx^3 + 3*s*(r-tt)*xx^2 + 3*s*r*(r-tt)*xx - s*tt^2 + s*r^3 - tt ] >;
H1, l1, H2, l2, G1, G2 := Explode(BFT(r3, s3, t3, x3));
assert G1^2 + l1*H1^3 eq G2^2 + l2*H2^3;
print "over Q(r,s,t): G1^2 + 4s H1^3 = G2^2 + 4st H2^3 (BFT Theorem 6)";

// (1) the slice r = 2, s = 3 and the marked group over Q(t)
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
H1, l1, H2, l2, G1, G2 := Explode(BFT(Ft!2, Ft!3, t, x));
f := G1^2 + l1*H1^3;
assert f eq G2^2 + l2*H2^3 and Degree(f) eq 6;
J := Jacobian(HyperellipticCurve(f));
T1 := elt< J | H1, G1 mod H1 >;
T2 := elt< J | H2, G2 mod H2 >;
assert T1 ne J!0 and 3*T1 eq J!0 and T2 ne J!0 and 3*T2 eq J!0;
assert T2 ne T1 and T2 ne -T1;
print "over Q(t), r = 2, s = 3: T1 = [H1, G1] and T2 = [H2, G2] have order 3 and T2 is not in <T1>: (Z/3)^2 is marked";

// (2), (3) residue certificate at t = 4 mod 5, t = 6 mod 13
g, st := ResidueCertificate(f, [5, 13], [4, 6], [19, 84, -46], 9, [3,3]);
print "FAMILY_33_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [3,3]";
quit;
