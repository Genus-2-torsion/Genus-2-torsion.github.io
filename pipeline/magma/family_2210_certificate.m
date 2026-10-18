/* family_2210_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,2,10].

   The family: the record's Elkies 5-torsion chart with four rational Weierstrass points
   (census_infinity_record.md section 2.7; torsion_jac/code/claude_census_infty_round7d.m).  Elkies'
   curve y^2 = g, g = 4Q^3 + (x^2 - 6x + 1)Q^2 + 2x(x-1)Q + x^2, has the 5-torsion class [D - K],
   D = {Q = 0, y = x}.  Here Q = q x + kappa Cq(x), Cq(x) = q^2 x^2 + 2q(2q-1)(q-1) x + (q-1)^2, with
   q = q(s): both roots of Cq are rational Weierstrass points; lam1 = lam(t1), lam2 = lam(3) are two more
   (lam(tau) = (q-1)(-(2q-1) + 2m)/q, q = tau^2/(tau^2-1), m = tau/(tau^2-1)).  The two conditions
   Q(lam_i) = q_i lam_i fix kappa and cut out (for t2 = 3) a nodal cubic in (t1, s), parametrized by
       t1 = -(45 + 3 tau^2)/(9 - 5 tau^2),  s = tau t1.

   Argument.
   (1) Over Q(tau): the two values of kappa agree, g has the four rational roots lam1, lam2 and the roots
       rho1, rho2 of Cq (checked).  D5 = [u, x mod u], u = Q/lc(Q), has order 5.  The classes
       [(x-lam1)(x-w), 0], w = lam2, rho1, rho2, are 2-torsion and all 7 nonzero sums of them are nonzero,
       so they generate (Z/2)^3; with D5 this marks Z/2 x Z/2 x Z/10 over Q(tau).
   (2) Residue certificate (family_cert_lib.m) at tau = 2 mod 5 and tau = 5 mod 13 (tau = 57 mod 65):
       #J(F_5) = 40, #J(F_13) = 120, gcd 40 = #G, strict at 13.  So every 5- and 13-integral
       tau = 57 mod 65 gives J(Q)_tors exactly [2,2,10] and a geometrically simple Jacobian.
   (3) The G2-invariants of two members differ, so infinitely many isomorphism classes.
   Run (from the census root):
       magma -b pipeline/magma/family_2210_certificate.m > data/knowledge/sources/family_2210_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// (1) the family and the marked group over Q(tau)
Ft<tau> := FunctionField(Q); P<x> := PolynomialRing(Ft);
function WPq(v) return v^2/(v^2-1), v/(v^2-1); end function;
t1 := -(45 + 3*tau^2)/(9 - 5*tau^2);
s := tau*t1;
q1, m1 := WPq(t1); lam1 := (q1-1)*(-(2*q1-1) + 2*m1)/q1;
q2, m2 := WPq(Ft!3); lam2 := (q2-1)*(-(2*q2-1) + 2*m2)/q2;
q, m := WPq(s);
Cq := P![ (q-1)^2, 2*q*(2*q-1)*(q-1), q^2 ];
kap1 := (q1 - q)*lam1/Evaluate(Cq, lam1);
kap2 := (q2 - q)*lam2/Evaluate(Cq, lam2);
assert kap1 eq kap2;
Qq := q*x + kap1*Cq;
assert Evaluate(Qq, lam1) eq q1*lam1 and Evaluate(Qq, lam2) eq q2*lam2;
f := 4*Qq^3 + (x^2 - 6*x + 1)*Qq^2 + 2*x*(x - 1)*Qq + x^2;
rho := [ rt[1] : rt in Roots(Cq) ];
assert #rho eq 2;
W := [lam1, lam2] cat rho;
assert Degree(f) eq 6 and #Seqset(W) eq 4 and &and[ Evaluate(f, w) eq 0 : w in W ];
J := Jacobian(HyperellipticCurve(f));
u := Qq/LeadingCoefficient(Qq);
assert (f - x^2) mod u eq 0;
D5 := elt< J | u, x mod u >;
assert D5 ne J!0 and 5*D5 eq J!0;
T := [ elt< J | (x - W[1])*(x - W[j]), 0 > : j in [2..4] ];
assert &and[ 2*T[i] eq J!0 : i in [1..3] ];
sums := [ &+[ T[i] : i in S ] : S in Subsets({1..3}) | #S gt 0 ];
assert #Seqset(sums) eq 7 and J!0 notin sums;
print "over Q(tau): kappa is consistent, g has 4 rational roots, D5 = [Q, x] has order 5 and the root-pair classes generate (Z/2)^3: [2,2,10] is marked";

// (2), (3) residue certificate at tau = 2 mod 5, tau = 5 mod 13
g, st := ResidueCertificate(f, [5, 13], [2, 5], [57, 122, -8], 40, [2,2,10]);
print "FAMILY_2210_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,2,10]";
quit;
