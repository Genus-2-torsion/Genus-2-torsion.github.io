/* family_3_9_certificate.m -- certificate that infinitely many pairwise non-isomorphic
   genus-2 curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [3,9].

   The family (torsion_jac/notes/simple39_2026_10_19.md; found 2026-10-19 on the pencil
   U9(d, 2/3, 1/2) of the [9]-chart of torsion_jac/notes/m318_2026_10_13.md, with t = 27(1 - d)/4).
       C_t : y^2 = f_t(x) = (t^2 x^3 + 3(t + 1) x + 1)^2 + 32 t^2 x^3,
       disc_x(f_t) = 2^21 3^6 t^16 (t - 2)^4 (4t + 1)^2,
   carries the class D = [x^2 + (3x + 1)/t, (2t - 12)x - 4] of order 9 over Q(t).  With
       B = t^2 x^3 + 6t^2 x^2 + (9t + 3)x + 2t + 1,   (4t + 1) f_t = B^2 + 4t^2 (t x^2 - x - 1)^3,
   so after the base change t = (s^2 - 1)/4 (4t + 1 = s^2) the divisor of y - B/s is 3 times a divisor
   over t x^2 - x - 1, and S = [x^2 - (x + 1)/t, 2s(2x + 1)] is a second rational class of order 3.

   Argument.
   (1) Over Q(s): 9D = 0, 3D != 0, 3S = 0, S != 0, S not in <3D>; the 27 sums iD + jS are distinct,
       so [3,9] is marked.
   (2) Residue certificate (family_cert_lib.m) at s = 5 mod 7 and s = 5 mod 17 (s = 5 mod 119):
       #J(F_7) = 81, #J(F_17) = 351, gcd 27 = #G, strict at 17.  So every 7- and 17-integral
       s = 5 mod 119 gives a geometrically simple Jacobian with J(Q)_tors exactly [3,9].
   (3) The G2-invariants of two members differ, so infinitely many isomorphism classes.
   Run (from the repository root):
       magma -b pipeline/magma/family_3_9_certificate.m > data/knowledge/sources/family_3_9_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// (0) the pencil over Q(t): the order-9 class and the contact identity
Ft<t> := FunctionField(Q); Pt<xt> := PolynomialRing(Ft);
ft := (t^2*xt^3 + 3*(t + 1)*xt + 1)^2 + 32*t^2*xt^3;
Bt := t^2*xt^3 + 6*t^2*xt^2 + (9*t + 3)*xt + 2*t + 1;
assert (4*t + 1)*ft eq Bt^2 + 4*t^2*(t*xt^2 - xt - 1)^3;
assert Discriminant(ft) eq 2^21*3^6*t^16*(t - 2)^4*(4*t + 1)^2;
Jt := Jacobian(HyperellipticCurve(ft));
Dt := elt< Jt | xt^2 + (3*xt + 1)/t, (2*t - 12)*xt - 4 >;
assert 9*Dt eq Jt!0 and 3*Dt ne Jt!0;
print "over Q(t): f_t = (t^2 x^3 + 3(t+1)x + 1)^2 + 32t^2 x^3, disc 2^21 3^6 t^16 (t-2)^4 (4t+1)^2; D has order 9; (4t+1) f_t = B^2 + 4t^2 (t x^2 - x - 1)^3";

// (1) the base change t = (s^2 - 1)/4 and the marked group over Q(s)
Fs<s> := FunctionField(Q); P<x> := PolynomialRing(Fs);
ts := (s^2 - 1)/4;
f := (ts^2*x^3 + 3*(ts + 1)*x + 1)^2 + 32*ts^2*x^3;
assert Degree(f) eq 6 and Discriminant(f) ne 0;
J := Jacobian(HyperellipticCurve(f));
D := elt< J | x^2 + (3*x + 1)/ts, (2*ts - 12)*x - 4 >;
S := elt< J | x^2 - (x + 1)/ts, 2*s*(2*x + 1) >;
assert 9*D eq J!0 and 3*D ne J!0;
assert 3*S eq J!0 and S ne J!0 and S ne 3*D and S ne -3*D;
grp := { i*D + j*S : i in [0..8], j in [0..2] };
assert #grp eq 27;
print "over Q(s), t = (s^2 - 1)/4: D (order 9) and S = [x^2 - (x+1)/t, 2s(2x+1)] (order 3, not in <3D>) generate a group of order 27 = [3,9]";

// (2), (3) residue certificate at s = 5 mod 7, s = 5 mod 17 (s = 5 mod 119)
g, st := ResidueCertificate(f, [7, 17], [5, 5], [Q| 5, 124, -114, 243], 27, [3,9]);
for m in [124, -114] do
    printf "sanity: s = %o has TorsionSubgroup %o\n", m, TorsionInvariants(SpecPoly(f, m));
end for;
print "FAMILY_3_9_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [3,9]";
quit;
