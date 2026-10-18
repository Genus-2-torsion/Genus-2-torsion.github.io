/* family_6_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [6].

   NOT the record's family.  The record (census_infinity_record.md, section 3) uses the plain M(12) chart
   T = a x^2 - x + r, y^2 + (x-r)(T+1) y = a x^2 T(T+1), whose class 2P - K (P = (0,0)) has order 6 and
   whose fibers are generically exact [6].  But that chart admits no residue-class certificate: on the
   slices r = 3 and a = 1, every good reduction found has 12 | #J(F_p), and on r = 3 the marked
   2-torsion class 3(2P - K) is divisible by 2 in J(F_p) in every good residue class
   (torsion_jac/notes/oldcert_B_2026_10_18.md).  So no good prime can cap the 2-part at Z/2.
   This certificate uses the contact-6 construction of the inventory's [6] row
   (torsion_jac/notes/contact6_m36.md) without the normalization h(0) = 1, which forces a root at x = 0:
       h = t + x + x^2 + x^3,  C_t : y^2 = f_t = h^2 - (x-1)^6 = (h - (x-1)^3)(h + (x-1)^3),
   of degree 5 with factorization type [2,3] over Q(t).  div(y - h) = 6P - 6 oo with P = (1, h(1)), so
   D6 = [x - 1, h(1)] has order 6 (checked over Q(t)).
   Certificate: t = 1 mod 3 gives #J(F_3) = 12; t = 5 mod 11 gives #J(F_11) = 102 with a strict
   Frobenius; so every 3- and 11-integral t = 16 mod 33 has torsion exactly [6] and a geometrically
   simple Jacobian (family_cert_lib.m).
   Run:  magma -b pipeline/magma/family_6_certificate.m > data/knowledge/sources/family_6_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
h := t + x + x^2 + x^3;
f := h^2 - (x-1)^6;
assert Degree(f) eq 5 and IsSquarefree(f);
printf "f_t = %o\n", f;
J := Jacobian(HyperellipticCurve(f));
D6 := elt< J | x - 1, Evaluate(h, 1) >;
assert 6*D6 eq J!0 and 2*D6 ne J!0 and 3*D6 ne J!0;
print "over Q(t): D6 = [x - 1, h(1)] has order 6";
g, st := ResidueCertificate(f, [3, 11], [1, 5], [16, -17, 49, 1/31], 6, [6]);
print "FAMILY_6_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [6]";
quit;
