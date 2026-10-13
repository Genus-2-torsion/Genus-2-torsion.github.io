/* family_228_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,2,8].

   The odd M_1(8,2,2) family of torsion_jac/notes/infinite_families_inventory.md (lines 429-451):
       f_{u,v} = ((1-u)x + 1)((1-v)x + 1)((u+v+2)x + 1)(-x^2 + (u^2 v - u^2 + u v^2 - u - v^2 - v - 2)x
                 - (u^2 + uv + v^2 + u + v + 1)),
   with Q = (-1, uv(u+v+1)).  On the slice v = 3, over Q(u): [Q - oo] has order 8, the three rational
   roots give three independent 2-torsion classes, and 4[Q - oo] is one of their combinations, so
   [2,2,8] is marked.  (Certificate found in the notes review of 2026-10-13.)
   Certificate: u = 4 mod 7 gives #J(F_7) = 64; u = 11 mod 53 gives #J(F_53) = 3296 with a strict
   Frobenius; gcd 32 = #[2,2,8].  So every u = 11 mod 371 (7- and 53-integral) has torsion exactly
   [2,2,8] and a geometrically simple Jacobian (family_cert_lib.m).
   Run:  magma -b pipeline/magma/family_228_certificate.m > data/knowledge/sources/family_228_certificate.log */
SetColumns(0);
load "pipeline/magma/family_cert_lib.m";
Fu<u> := FunctionField(Q); P<x> := PolynomialRing(Fu);
v := 3;
quad := -x^2 + (u^2*v - u^2 + u*v^2 - u - v^2 - v - 2)*x - (u^2 + u*v + v^2 + u + v + 1);
lins := [ (1-u)*x + 1, (1-v)*x + 1, (u+v+2)*x + 1 ];
f := &*lins * quad;
assert Degree(f) eq 5 and Evaluate(f, -1) eq (u*v*(u+v+1))^2;
J := Jacobian(HyperellipticCurve(f));
Qc := elt< J | x + 1, u*v*(u+v+1) >;
assert 8*Qc eq J!0 and 4*Qc ne J!0;
Ts := [ elt< J | l/LeadingCoefficient(l), 0 > : l in lins ];
combos := { a*Ts[1] + b*Ts[2] + c*Ts[3] : a, b, c in [0,1] };
assert #combos eq 8 and 4*Qc in combos;
print "over Q(u), v = 3: [Q - oo] has order 8, the rational roots give (Z/2)^3, and 4[Q - oo] lies in it: [2,2,8] is marked";
g, st := ResidueCertificate(f, [7, 53], [4, 11], [11, 382, -360, 753], 32, [2,2,8]);
print "FAMILY_228_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,2,8]";
quit;
