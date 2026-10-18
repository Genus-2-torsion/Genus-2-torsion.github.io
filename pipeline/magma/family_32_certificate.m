/* family_32_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [32].

   The family: the record's "reconstructed Elkies order-32 component" (torsion_jac/notes/
   elkies32_reconstruction.md).  Curves y^2 = h^2 - a^2 x^5 (x+1), h = a x^3 + b x^2 + c x + 1, with
   D = [(0,1) - oo]; div(y - h) = 5 (0,1) + (-1, h(-1)) - 6 oo, and on the component 16 D = [(r,0) - oo].
   The note gives the component as a genus-0 plane curve of degree 8 in (z, r) = (a/p, r), p = a - b + c - 1,
   with c = (-a^4 + 4a^3 p - 8a^2 p^2 + 8a p^3 + 8a^2 p - 16p^3)/(4a^2 p), b = a + c - 1 - p, and records
   only an implicit form.  A parametrization (found 2026-10-18 from Magma's Parametrization at the printed
   point, moved by a Moebius change of parameter) is
       z = 2(w^4 - w^3 + w^2 - w + 1)/(w^2 (w^2 - w + 1)),   r = 1/(w^4 - 1),
       p = -w^5 (w^6 - 2w^5 + 4w^4 - 4w^3 + 4w^2 - 2w + 1)/(w^4 - w^3 + w^2 - w + 1),
   with p the root of f(r) = 0 that passes through the printed point; w = 2 is Elkies' printed curve
   y^2 = (15x - 1)(1056x^4 + 156183x^3 + 26297x^2 + 649x - 121).

   Argument.
   (1) Over Q(w): (z, r) satisfies the component equation of the note, f(r) = 0, and D = [x, 1] satisfies
       16 D = [x - r, 0] != 0 and 32 D = 0, so D has order 32.
   (2) Residue certificate (family_cert_lib.m) at w = 2 mod 7 and w = 3 mod 13 (w = 16 mod 91):
       #J(F_7) = 64, #J(F_13) = 160, gcd 32 = #G, strict at 7.  So every 7- and 13-integral w = 16 mod 91
       gives J(Q)_tors exactly [32] and a geometrically simple Jacobian.
   (3) The G2-invariants of two members differ, so infinitely many isomorphism classes.
   Run (from the census root):
       magma -b pipeline/magma/family_32_certificate.m > data/knowledge/sources/family_32_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// (1) the family and the marked class over Q(w)
Fw<w> := FunctionField(Q); P<x> := PolynomialRing(Fw);
z := 2*(w^4 - w^3 + w^2 - w + 1)/(w^2*(w^2 - w + 1));
r := 1/(w^4 - 1);
Comp := func< Z, Rr | 3*Z^4*Rr^4 + 9*Z^4*Rr^3 - 16*Z^3*Rr^4 + 10*Z^4*Rr^2 - 56*Z^3*Rr^3 + 32*Z^2*Rr^4
    + 5*Z^4*Rr - 72*Z^3*Rr^2 + 144*Z^2*Rr^3 - 64*Z*Rr^4 + Z^4 - 40*Z^3*Rr + 208*Z^2*Rr^2 - 224*Z*Rr^3
    + 80*Rr^4 - 8*Z^3 + 120*Z^2*Rr - 288*Z*Rr^2 + 160*Rr^3 + 24*Z^2 - 160*Z*Rr + 160*Rr^2 - 32*Z + 80*Rr + 16 >;
assert Comp(z, r) eq 0;
p := -w^5*(w^6 - 2*w^5 + 4*w^4 - 4*w^3 + 4*w^2 - 2*w + 1)/(w^4 - w^3 + w^2 - w + 1);
a := z*p;
c := (-a^4 + 4*a^3*p - 8*a^2*p^2 + 8*a*p^3 + 8*a^2*p - 16*p^3)/(4*a^2*p);
b := a + c - 1 - p;
assert p eq a - b + c - 1;
h := a*x^3 + b*x^2 + c*x + 1;
f := h^2 - a^2*x^5*(x + 1);
assert Degree(f) eq 5 and Evaluate(f, r) eq 0;
printf "a = %o\nb = %o\nc = %o\n", a, b, c;
J := Jacobian(HyperellipticCurve(f));
D := elt< J | x, 1 >;
D16 := 16*D;
assert D16 eq elt< J | x - r, 0 > and D16 ne J!0 and 2*D16 eq J!0;
print "over Q(w): (z, r) lies on the reconstructed component, f(r) = 0, and D = [x, 1] has 16 D = [x - r, 0] != 0, 32 D = 0: order 32";

// positive control: w = 2 is Elkies' printed curve
X := RX.1;
Celk := HyperellipticCurve((15*X - 1)*(1056*X^4 + 156183*X^3 + 26297*X^2 + 649*X - 121));
assert G2Invariants(HyperellipticCurve(SpecPoly(f, 2))) eq G2Invariants(Celk);
print "w = 2: the fiber has the G2-invariants of Elkies' printed order-32 curve";

// (2), (3) residue certificate at w = 2 mod 7, w = 3 mod 13
g, st := ResidueCertificate(f, [7, 13], [2, 3], [16, 107, -75], 32, [32]);
print "FAMILY_32_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [32]";
quit;
