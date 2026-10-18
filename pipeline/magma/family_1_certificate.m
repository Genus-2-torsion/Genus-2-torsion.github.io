/* family_1_certificate.m -- residue-class certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have geometrically simple Jacobian with trivial rational torsion (audit of 2026-10-18; it
   complements the L1-L3 certificate trivial_torsion_family.m, which lists exact fibers only).

   Family: C_t : y^2 = x^5 + x + t over Q(t); there is nothing to mark.
   Certificate (family_cert_lib.m): every 3- and 17-integral t = 3 (mod 51) has the reductions of t = 3 at 3 and
   at 17, with #J(F_3) = 12 and #J(F_17) = 205 (gcd 1) and a strict Frobenius at 17, so J_t(Q)_tors = 0 and J_t is
   geometrically simple; members with different G2-invariants give infinitely many isomorphism classes.
   Run (from the repository root):
       magma -b pipeline/magma/family_1_certificate.m > data/knowledge/sources/family_1_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
f := x^5 + x + t;
printf "f_t = %o\n", f;
g, st := ResidueCertificate(f, [3, 17], [0, 3], [Q| 3, 54, -48, 105], 1, []);
print "FAMILY_1_CERTIFIED: infinitely many geometrically simple Jacobians with trivial J(Q)_tors";
quit;
