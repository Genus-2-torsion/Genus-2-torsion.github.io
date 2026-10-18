/* family_30_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [30].

   The family: the record's simultaneous contact-5/contact-6 family
   (torsion_jac/notes/contact5_contact6_order30_family.md), branch -1.  Odd quintics with
       f = h6^2 - (x-1)^6 = h5^2 - K x^5,  h6 = x^3 + A x^2 + B x + C,  h5 = e x^2 + d x + c,
   on the genus-0 branch
       t = (5R^2 - 20R + 19)/(R^2 - 5),  Y = -2(5R^2 - 22R + 25)/(R^2 - 5)  (Y^2 = 5t^2 - 6t + 5),
       u = t^3,  s = t^5 + t^4 + (5/2)t^3 + (1/2)t - t(t - 1/2)(t + 1)Y,
       q = (15u^5 + 90u^4 + 20u^3 s - 6u^2 s^2 + 231u^3 + 2u^2 s - 15u s^2 + 90u^2 - 20u s + 15u - 2s)
           /(u^6 + 6u^4 s - 2u^4 + 15u^3 s - u s^3 + u^2),
       A = (s+q)/2, e = (s-q)/2, B = (15 - s q)/2, C = (u^2+1)/(2u), c = (u^2-1)/(2u),
       d = (B C + 3)/c, K = -2A - 6.

   Argument.
   (1) Over Q(R): the two expressions for f agree, so div(y - h5) = 5 (0,c) - 5 oo and
       div(y - h6) = 6 (1, h6(1)) - 6 oo.  D5 = [x, c] has order 5 and D6 = [x - 1, h6(1)] has order 6
       (checked: 5 D5 = 0 != D5, 6 D6 = 0, 2 D6 != 0, 3 D6 != 0), so D5 + D6 has order 30.
   (2) Residue certificate (family_cert_lib.m) at R = 0 mod 11 and R = 6 mod 13 (R = 110 mod 143):
       #J(F_11) = 150, #J(F_13) = 120, gcd 30 = #G, strict at 13.  So every 11- and 13-integral
       R = 110 mod 143 gives J(Q)_tors exactly [30] and a geometrically simple Jacobian.
   (3) The G2-invariants of two members differ, so infinitely many isomorphism classes.
   Run (from the census root):
       magma -b pipeline/magma/family_30_certificate.m > data/knowledge/sources/family_30_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// (1) the family and the marked class over Q(R)
FR<R> := FunctionField(Q); P<x> := PolynomialRing(FR);
t := (5*R^2 - 20*R + 19)/(R^2 - 5);
Y := -2*(5*R^2 - 22*R + 25)/(R^2 - 5);
assert Y^2 eq 5*t^2 - 6*t + 5;
u := t^3;
s := t^5 + t^4 + (5/2)*t^3 + (1/2)*t - t*(t - 1/2)*(t + 1)*Y;
q := (15*u^5 + 90*u^4 + 20*u^3*s - 6*u^2*s^2 + 231*u^3 + 2*u^2*s - 15*u*s^2 + 90*u^2 - 20*u*s + 15*u - 2*s)
     /(u^6 + 6*u^4*s - 2*u^4 + 15*u^3*s - u*s^3 + u^2);
A := (s + q)/2; e := (s - q)/2; B := (15 - s*q)/2;
CC := (u^2 + 1)/(2*u); c := (u^2 - 1)/(2*u);
d := (B*CC + 3)/c; K := -2*A - 6;
h6 := x^3 + A*x^2 + B*x + CC; h5 := e*x^2 + d*x + c;
f := h6^2 - (x - 1)^6;
assert f eq h5^2 - K*x^5 and Degree(f) eq 5;
J := Jacobian(HyperellipticCurve(f));
D5 := elt< J | x, c >;
D6 := elt< J | x - 1, Evaluate(h6, 1) >;
assert D5 ne J!0 and 5*D5 eq J!0;
assert 6*D6 eq J!0 and 2*D6 ne J!0 and 3*D6 ne J!0;
print "over Q(R): f = h6^2 - (x-1)^6 = h5^2 - K x^5; D5 = [x, c] has order 5 and D6 = [x-1, h6(1)] has order 6, so D5 + D6 has order 30";

// (2), (3) residue certificate at R = 0 mod 11, R = 6 mod 13
g, st := ResidueCertificate(f, [11, 13], [0, 6], [110, 253, -33], 30, [30]);
print "FAMILY_30_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [30]";
quit;
