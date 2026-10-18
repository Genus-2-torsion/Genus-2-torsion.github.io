/* family_2_certificate.m -- residue-class certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [2].

   It replaces the family y^2 = x(x^5 + x + t) of the archived record (census_infinity_record.md, section 2.1),
   whose generic member has NO rational 2-torsion: the quintic x^5 + x + t is irreducible over Q(t), so the even
   model has only two factors and 2-rank 0; the record's fiber t = 2 is exact [2] only because
   x^5 + x + 2 = (x + 1)(x^4 - x^3 + x^2 - x + 2) (found by the audit of 2026-09-18).

   Family: C_t : y^2 = x (x^4 + x + t) over Q(t).  (0,0) and the point at infinity are rational Weierstrass
   points, so [(0,0) - oo] is a rational point of order 2 over Q(t) (checked below), hence on every smooth fiber.
   Certificate (family_cert_lib.m): every 3- and 5-integral t = 2 (mod 15) has the reductions of t = 2 at 3 and
   at 5, with #J(F_3) = 22 and #J(F_5) = 20 (gcd 2) and a strict Frobenius at 3, so J_t(Q)_tors = [2] and J_t is
   geometrically simple; members with different G2-invariants give infinitely many isomorphism classes.
   (The version of 2026-09-18 used the class t = 1 (mod 15), where the gcd is 4, so exactness there needed a
   Hilbert-irreducibility step; corrected by the audit of 2026-10-18.)
   Run (from the repository root):
       magma -b pipeline/magma/family_2_certificate.m > data/knowledge/sources/family_2_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
f := x*(x^4 + x + t);
printf "f_t = %o\n", f;
J := Jacobian(HyperellipticCurve(f));
T := elt< J | x, 0 >;
assert T ne J!0 and 2*T eq J!0;
print "over Q(t): [(0,0) - oo] has order 2 on y^2 = x(x^4 + x + t)";
assert IsIrreducible(x^4 + x + t);
g, st := ResidueCertificate(f, [3, 5], [2, 2], [Q| 2, 17, -13, 32, 1/8], 2, [2]);
print "FAMILY_2_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2]";
quit;
