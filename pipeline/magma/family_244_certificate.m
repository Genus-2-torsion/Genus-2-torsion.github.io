/* family_244_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,4,4].

   The family.  Nicholls' three-parameter family Lambda_334 (C. Nicholls, DPhil thesis, Oxford 2018,
   Prop. 5.9.6), in the form used by BNSS for their [8,8] lane (Genus2TorsionNotebook,
   code/claude_ov_88verify.m, function MemberSTV; transcribed here, not run):
       A = s^2 - t^4 + t^2,  u = (-s^2 A v^2 - 2 A v - 1)/(-s^2 t A v^2 + t),
       a = A/(1 - t^2),  b = A/(u^2 s^2 + 1 - t^2),  c = t^2,
       d2 = A (s^2 u^2 + t^4 - 2t^2 + 1)(s^4 u^2 - s^2 t^2 u^2 + s^2 u^2 - t^6 + 3t^4 - 3t^2 + 1),
       l1 = (-a+b+c-1)x^2 + (2a - 2bc)x + (abc - ab - ac + bc),  l2 = -x^2 + bc,  l3 = x^2 - a,
       C : y^2 = f = d2 l1 l2 l3 / lc(l1 l2 l3).
   bc is a square, so l2 splits.  l3 splits exactly when a = w^2, and this locus is rational:
       w = (m^2 t^2 - t^4 + t^2)/((m^2 + t^2 - 1) t),  s = m (w - t)
   (BNSS, notes/claude_ov_88_2026_07_25.md, logged 27 such members with torsion [2,4,4]; the rational
   parametrisation is from the review of 2026-10-14, ~/random_research/bnss_repos/review_A_notes.md, 1.1).
   This certificate uses the slice t = m = 2 (s = -24/7, w = 2/7, a = 4/49, c = 4), with parameter v.

   Argument.
   (1) Over Q(v): l2 = -(x - r)(x + r) and l3 = (x - w)(x + w) split, so f has factor type [1,1,1,1,2].
       Halving identities: for R = l3 and R = l2, with F = f/R, there are a linear l, a constant c and a
       monic quadratic u in Q(v)[x] with F = c u^2 + R l^2.  Then (y - Rl)(y + Rl) = c R u^2, so
       div(y - Rl) = (zeros of R) + 2 D_u - 3 D_inf and D = [u, Rl mod u] satisfies 2D = [R].  The
       identities are found by solving for (l, c, u) over Q(v) (a resultant in the two coefficients of l);
       the proof is the Mumford check over Q(v): 2 D3 = [l3], 2 D2 = [l2], 4 D3 = 4 D2 = 0.
       With T = [(x - w)(x - r)] (a third 2-torsion class), the 32 sums a D3 + b D2 + c T
       (a, b in Z/4, c in Z/2) are distinct, so Z/2 x Z/4 x Z/4 is marked over Q(v).
   (2) Residue certificate (family_cert_lib.m) at v = 2 mod 19 and v = 2 mod 31, i.e. every 19- and
       31-integral v = 2 mod 589: #J(F_19) = 288, #J(F_31) = 1088, gcd 32 = #[2,4,4]; the Frobenius at 31
       is strict, so these Jacobians are geometrically simple.  The member v = 2, i.e. (t,m,v) = (2,2,2),
       has TorsionSubgroup [2,4,4].  (The library's log lines call the family parameter t; here it is v.)
   (3) The G2-invariants of the members differ, so v -> M_2 is nonconstant and the class gives
       infinitely many isomorphism classes.
   Run (from the repository root):
       magma -b pipeline/magma/family_244_certificate.m > family_244_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// (0) the slice t = m = 2 of the l3-split locus of Lambda_334, over Q(v)
Fv<v> := FunctionField(Q); P<x> := PolynomialRing(Fv);
tt := Fv!2; mm := Fv!2;
w := (mm^2*tt^2 - tt^4 + tt^2)/((mm^2 + tt^2 - 1)*tt);
s := mm*(w - tt);
AA := s^2 - tt^4 + tt^2;
uu := (-s^2*AA*v^2 - 2*AA*v - 1)/(-s^2*tt*AA*v^2 + tt);
a := AA/(1 - tt^2); b := AA/(uu^2*s^2 + 1 - tt^2); c := tt^2;
d2 := AA*(s^2*uu^2 + tt^4 - 2*tt^2 + 1)*(s^4*uu^2 - s^2*tt^2*uu^2 + s^2*uu^2 - tt^6 + 3*tt^4 - 3*tt^2 + 1);
l1 := (-a + b + c - 1)*x^2 + (2*a - 2*b*c)*x + (a*b*c - a*b - a*c + b*c);
l2 := -x^2 + b*c; l3 := x^2 - a;
f1 := l1*l2*l3;
f := d2*f1/LeadingCoefficient(f1);
assert Degree(f) eq 6;
assert a eq w^2;
ok, r := IsSquare(b*c); assert ok;
assert l3 eq (x - w)*(x + w) and l2 eq -(x - r)*(x + r);
printf "slice t = m = 2: s = %o, w = %o, a = %o, c = %o; over Q(v): l3 = (x - w)(x + w), l2 = -(x - r)(x + r) with r = %o\n",
    s, w, a, c, r;
J := Jacobian(HyperellipticCurve(f));

// (1) halves of [R] for R | f quadratic: solve F = c u^2 + R l^2, l = L0 + L1 x, over Q(v)
R2<L0,L1> := PolynomialRing(Fv, 2);
PX<XX> := PolynomialRing(R2);
function Halves(R)
    F := f div R; assert F*R eq f;
    G := PX!F - (PX!R)*(L0 + L1*XX)^2;
    g := [ Coefficient(G, i) : i in [0..4] ];
    // G = g4 (X^2 + u1 X + u0)^2  <=>  e1 = e2 = 0
    e1 := 8*g[5]^2*g[2] - 4*g[5]*g[4]*g[3] + g[4]^3;
    e2 := 64*g[5]^3*g[1] - (4*g[5]*g[3] - g[4]^2)^2;
    res := UnivariatePolynomial(Resultant(e1, e2, L0));
    out := [];
    for h in Factorization(res) do
        if Degree(h[1]) ne 1 then continue; end if;
        lam1 := -Coefficient(h[1], 0)/Coefficient(h[1], 1);
        q1 := UnivariatePolynomial(Evaluate(e1, [L0, R2!lam1]));
        q2 := UnivariatePolynomial(Evaluate(e2, [L0, R2!lam1]));
        for rt in Roots(GCD(q1, q2)) do
            l := rt[1] + lam1*x;
            Gx := F - R*l^2;
            cc := LeadingCoefficient(Gx);
            if cc eq 0 or Degree(Gx) ne 4 then continue; end if;
            okk, uD := IsSquare(Gx/cc);
            if not okk then continue; end if;
            assert F eq cc*uD^2 + R*l^2;
            Append(~out, < elt< J | uD, (R*l) mod uD >, l, cc >);
        end for;
    end for;
    return out;
end function;

T3 := elt< J | x^2 - a, 0 >;
T2 := elt< J | x^2 - b*c, 0 >;
Tx := elt< J | (x - w)*(x - r), 0 >;
assert T3 ne J!0 and T2 ne J!0 and Tx ne J!0 and #{ T3, T2, Tx, T3 + T2, T3 + Tx, T2 + Tx, T3 + T2 + Tx } eq 7;
print "over Q(v): [l3], [l2], [(x - w)(x - r)] are independent 2-torsion classes";
H3 := Halves(l3); H2 := Halves(l2);
assert #H3 ge 1 and #H2 ge 1;
D3 := H3[1][1]; D2 := H2[1][1];
printf "halving identity for l3: l = %o, c = %o\n", H3[1][2], H3[1][3];
printf "D3 = %o\n", D3;
printf "halving identity for l2: l = %o, c = %o\n", H2[1][2], H2[1][3];
printf "D2 = %o\n", D2;
assert 2*D3 eq T3 and 2*D2 eq T2 and 4*D3 eq J!0 and 4*D2 eq J!0;
print "over Q(v): 2 D3 = [l3], 2 D2 = [l2], 4 D3 = 4 D2 = 0";
sums := { i*D3 + j*D2 + k*Tx : i, j in [0..3], k in [0..1] };
assert #sums eq 32;
print "over Q(v): the 32 sums i D3 + j D2 + k T are distinct: [2,4,4] is marked";

// (2) + (3) residue certificate, sanity at v = 2, distinct G2-invariants
g, st := ResidueCertificate(f, [19, 31], [2, 2], [Q| 2, 591, -587, 591/590], 32, [2,4,4]);
print "FAMILY_244_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,4,4]";
quit;
