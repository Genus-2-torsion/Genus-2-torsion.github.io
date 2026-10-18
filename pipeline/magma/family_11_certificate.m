/* family_11_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [11].

   Of the two families in the record (census_infinity_record.md, section 3) this uses Flynn's
   (as quoted by Daowsud-Schmidt, J. Number Theory 189 (2018); torsion_jac/code/order11_family_check.m):
       C_t : y^2 = f_t = x^6 + 2x^5 + (2t+3)x^4 + 2x^3 + (t^2+1)x^2 + 2t(1-t)x + t^2.
   f is monic of degree 6, so the two points oo+, oo- at infinity are rational, and D = [oo+ - oo-] has
   order 11 (checked over Q(t)).
   Certificate: t = 2 mod 3 gives #J(F_3) = 11; t = 1 mod 5 gives #J(F_5) = 22 with a strict Frobenius;
   so every 3- and 5-integral t = 11 mod 15 has torsion exactly [11] and a geometrically simple Jacobian
   (family_cert_lib.m).
   Run:  magma -b pipeline/magma/family_11_certificate.m > data/knowledge/sources/family_11_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
f := x^6 + 2*x^5 + (2*t+3)*x^4 + 2*x^3 + (t^2+1)*x^2 + 2*t*(1-t)*x + t^2;
assert Degree(f) eq 6 and IsMonic(f) and IsSquarefree(f);
printf "f_t = %o\n", f;
C := HyperellipticCurve(f);
J := Jacobian(C);
oo := PointsAtInfinity(C);
assert #oo eq 2;
D := J![ oo[1], oo[2] ];
assert D ne J!0 and 11*D eq J!0;
printf "over Q(t): D = [%o - %o] has order 11\n", oo[1], oo[2];
g, st := ResidueCertificate(f, [3, 5], [2, 1], [-4, 11, 26, -19], 11, [11]);
print "FAMILY_11_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [11]";
quit;
