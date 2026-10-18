/* family_210_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,10].

   The family: the record's Elkies 5-torsion chart with three rational Weierstrass points
   (census_infinity_record.md section 2.7; torsion_jac/code/claude_census_infty_round6b.m), slice (t, 2, 3).
   Elkies (Contemp. Math. 796 (2024), Thm 2.3 / Cor 2.4): for a quadratic Q the curve
       y^2 = g,  g = 4Q^3 + (x^2 - 6x + 1)Q^2 + 2x(x-1)Q + x^2
   carries the 5-torsion class [D - K], D = {Q = 0, y = x} (on Q = 0, g = x^2).  A parameter tau gives the
   Weierstrass point lam(tau) with q(tau) = Q(lam)/lam: q = tau^2/(tau^2-1), m = tau/(tau^2-1),
   lam = (q-1)(-(2q-1) + 2m)/q.  Q is the quadratic interpolating Q(lam_i) = q_i lam_i for tau = t, 2, 3.

   Argument.
   (1) Over Q(t): g has the three rational roots lam_1, lam_2, lam_3 (checked).  D5 = [u, x mod u],
       u = Q/lc(Q), satisfies 5 D5 = 0, D5 != 0.  T12 = [(x-lam1)(x-lam2), 0] and T13 = [(x-lam1)(x-lam3), 0]
       are distinct nonzero 2-torsion classes, so <T12, T13> = (Z/2)^2 and <D5, T12, T13> = Z/2 x Z/10
       is marked over Q(t).
   (2) Residue certificate (family_cert_lib.m) at t = 1 mod 5 and t = 14 mod 17 (t = 31 mod 85):
       #J(F_5) = 40, #J(F_17) = 260, gcd 20 = #G, strict at 17.  So every 5- and 17-integral t = 31 mod 85
       gives J(Q)_tors exactly [2,10] and a geometrically simple Jacobian.
   (3) The G2-invariants of two members differ, so infinitely many isomorphism classes.
   Run (from the census root):
       magma -b pipeline/magma/family_210_certificate.m > data/knowledge/sources/family_210_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// (1) the family and the marked group over Q(t)
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
function WPdata(tau)
    q := tau^2/(tau^2-1); m := tau/(tau^2-1);
    assert m^2 eq q*(q-1);
    return (q-1)*(-(2*q-1) + 2*m)/q, q;
end function;
lams := []; qs := [];
for tau in [t, Ft!2, Ft!3] do lam, q := WPdata(tau); Append(~lams, lam); Append(~qs, q); end for;
assert #Seqset(lams) eq 3;
Qq := Interpolation(lams, [ qs[i]*lams[i] : i in [1..3] ]);
assert Degree(Qq) eq 2;
f := 4*Qq^3 + (x^2 - 6*x + 1)*Qq^2 + 2*x*(x - 1)*Qq + x^2;
assert Degree(f) eq 6 and &and[ Evaluate(f, l) eq 0 : l in lams ];
J := Jacobian(HyperellipticCurve(f));
u := Qq/LeadingCoefficient(Qq);
assert (f - x^2) mod u eq 0;      // the points {Q = 0, y = x} lie on the curve
D5 := elt< J | u, x mod u >;
assert D5 ne J!0 and 5*D5 eq J!0;
T12 := elt< J | (x - lams[1])*(x - lams[2]), 0 >;
T13 := elt< J | (x - lams[1])*(x - lams[3]), 0 >;
assert T12 ne J!0 and T13 ne J!0 and T12 ne T13 and 2*T12 eq J!0 and 2*T13 eq J!0;
print "over Q(t): g has the rational roots lam(t), lam(2), lam(3); D5 = [Q, x] has order 5; T12, T13 generate (Z/2)^2: [2,10] is marked";

// (2), (3) residue certificate at t = 1 mod 5, t = 14 mod 17
g, st := ResidueCertificate(f, [5, 17], [1, 14], [31, 116, -54], 20, [2,10]);
print "FAMILY_210_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,10]";
quit;
