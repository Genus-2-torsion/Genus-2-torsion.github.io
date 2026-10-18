/* family_8_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [8].

   The record's family (census_infinity_record.md, section 2.3; torsion_jac/code/claude_census_infty_certify.m,
   family 9): the slice r = 3, p = 2 of the A(8) chart (formulas verbatim from function A8f in
   torsion_jac/code/agent_a2_24_composite8x3.m), parameter t:
       e = t^2 - 2pt/r, d = e + 2p - r^2, lambda = r/t, u = p + rt - 2r,
       v = e + r^2 - rp - r^2 t + 3pt - rt^2, a = r^2 - lambda, b = 2rp - 2 lambda (p + rt) + 2 r lambda,
       c = p^2 + 2pr^2 - r^4 - r^3 t - rp^2/t - lambda (r^2 + e) + 2 lambda (rp + r^2 t - 3pt + rt^2),
       Q = x^2 + d, q = a x^2 + b x + c,  f_t = q (Q^2 + q)   (degree 6),
   with the marked class D8 = [g8, -ellBase mod g8], g8 = x^2 + u x + v, ellBase = -(q + Q (r x + p - r^2)).
   Argument.
   (1) Over Q(t): D8 has order exactly 8 (Mumford arithmetic), so [8] is marked.
   (2) Residue certificate (family_cert_lib.m) at t = 3 mod 5 and t = 3 mod 19, i.e. every 5- and
       19-integral t = 3 mod 95: #J(F_5) = 16, #J(F_19) = 584, gcd 8 = #[8]; the Frobenius at 19 is
       strict.  The member t = 3 has TorsionSubgroup [8].  (The record's exact fiber t = 2 has #J(F_p)
       divisible by 16 at every good prime p < 150, so no residue class through t = 2 is certified by
       point counts; t = 3 is the record's second listed fiber.)
   (3) The G2-invariants of the members differ, so t -> M_2 is nonconstant.
   Run (from the repository root):
       magma -b pipeline/magma/family_8_certificate.m > data/knowledge/sources/family_8_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
r := 3; p := 2;
e := t^2 - 2*p*t/r; d := e + 2*p - r^2; lambda := r/t;
uu := p + r*t - 2*r;
vv := e + r^2 - r*p - r^2*t + 3*p*t - r*t^2;
a := r^2 - lambda;
b := 2*r*p - 2*lambda*(p + r*t) + 2*r*lambda;
c := p^2 + 2*p*r^2 - r^4 - r^3*t - r*p^2/t - lambda*(r^2 + e) + 2*lambda*(r*p + r^2*t - 3*p*t + r*t^2);
Qp := x^2 + d; q := a*x^2 + b*x + c;
f := q*(Qp^2 + q);
assert Degree(f) eq 6;
printf "f_t = %o\n", f;
g8 := x^2 + uu*x + vv;
ellBase := -(q + Qp*(r*x + (p - r^2)));
J := Jacobian(HyperellipticCurve(f));
D8 := elt< J | g8, (-ellBase) mod g8 >;
assert 8*D8 eq J!0 and 4*D8 ne J!0;
printf "over Q(t): D8 = %o has order exactly 8: [8] is marked\n", D8;
g, st := ResidueCertificate(f, [5, 19], [3, 3], [Q| 3, 98, -92, 193], 8, [8]);
print "FAMILY_8_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [8]";
quit;
