/* family_23_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [23].

   The family: the Kuru-Sadek construction (arXiv:2410.14455), genus-2 specialization, as transcribed in
   torsion_jac/code/order23_kuru_sadek_check.m (alpha read with beta^((g+1)/2), which reproduces their
   t = 2 example):
       beta = (t^2+1)^2/(4t^2), sbeta = (t^2+1)/(2t), s = (t^2-1)/(2t),
       alpha = beta - s^5/(beta sbeta), lambda = (alpha-1)^4/((alpha-beta)^2 alpha),
       A = (x^3 (x-alpha)^2 - (x-1)((x-1)^4 - lambda (x-beta)^2 x)) / (2 (x-alpha)(x-beta)),
       f_t = A^2 - lambda x^4 (x-1).

   Argument.
   (1) Over Q(t): A is a polynomial (exact division), f_t has degree 5, and D = [x, A(0)] satisfies
       23 D = 0, D != 0, so it has order 23.
   (2) Residue certificate (family_cert_lib.m) at t = 3 mod 7 and t = 4 mod 13 (t = 17 mod 91):
       #J(F_7) = 46, #J(F_13) = 115, gcd 23 = #G, strict at 7.  So every 7- and 13-integral t = 17 mod 91
       gives J(Q)_tors exactly [23] and a geometrically simple Jacobian.
   (3) The G2-invariants of two members differ, so infinitely many isomorphism classes.
   Run (from the census root):
       magma -b pipeline/magma/family_23_certificate.m > data/knowledge/sources/family_23_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// (1) the family and the marked class over Q(t)
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
beta := (t^2 + 1)^2/(4*t^2);
sbeta := (t^2 + 1)/(2*t);
s := (t^2 - 1)/(2*t);
assert sbeta^2 eq beta;
alpha := beta - s^5/(beta*sbeta);
lambda := (alpha - 1)^4/((alpha - beta)^2*alpha);
num := x^3*(x - alpha)^2 - (x - 1)*((x - 1)^4 - lambda*(x - beta)^2*x);
den := 2*(x - alpha)*(x - beta);
A, rem := Quotrem(num, den);
assert rem eq 0;
f := A^2 - lambda*x^4*(x - 1);
assert Degree(f) eq 5;
J := Jacobian(HyperellipticCurve(f));
D := elt< J | x, Evaluate(A, 0) >;
assert D ne J!0 and 23*D eq J!0;
print "over Q(t): A = (x^3(x-alpha)^2 - (x-1)((x-1)^4 - lambda(x-beta)^2 x))/(2(x-alpha)(x-beta)) is a polynomial and D = [x, A(0)] has order 23";

// (2), (3) residue certificate at t = 3 mod 7, t = 4 mod 13
g, st := ResidueCertificate(f, [7, 13], [3, 4], [17, 108, -74], 23, [23]);
print "FAMILY_23_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [23]";
quit;
