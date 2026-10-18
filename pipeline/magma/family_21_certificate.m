/* family_21_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [21].

   The family: Leprevost's 1991 one-parameter family, as transcribed in
   torsion_jac/code/z21_leprevost_family_verify.m (notes/z21_leprevost_and_hlp_2026_07_18.md):
       p2 = t^14+4t^13+19t^12+32t^11+113t^10+188t^9+379t^8+448t^7+379t^6+188t^5+113t^4+32t^3+19t^2+4t+1,
       p1 = t^10+4t^9+17t^8+24t^7+54t^6+56t^5+54t^4+24t^3+17t^2+4t+1,
       p0 = t^6+4t^5+15t^4+16t^3+15t^2+4t+1,
       A21 = p2 x^2 - 2(t^2+1)^2 p1 x + (t^2+1)^4 p0,
       k21 = 64 t^4 (t+1)^2 (t^2+1)^3 (t^4+2t^3+6t^2+2t+1)^3,
       f_t = A21^2 - k21 x^3 (x-1)^2.

   Argument.
   (1) Over Q(t): D = [x, A21(0)] satisfies 21 D = 0, 3 D != 0, 7 D != 0, so it has order 21.
   (2) Residue certificate (family_cert_lib.m) at t = 1 mod 5: #J(F_5) = 21 = #G, strict.  So every
       5-integral t = 1 mod 5 gives J(Q)_tors exactly [21] and a geometrically simple Jacobian.
   (3) The G2-invariants of two members differ, so infinitely many isomorphism classes.
   Run (from the census root):
       magma -b pipeline/magma/family_21_certificate.m > data/knowledge/sources/family_21_certificate.log */
SetColumns(0); SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// (1) the family and the marked class over Q(t)
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
p2 := t^14+4*t^13+19*t^12+32*t^11+113*t^10+188*t^9+379*t^8+448*t^7+379*t^6+188*t^5+113*t^4+32*t^3+19*t^2+4*t+1;
p1 := t^10+4*t^9+17*t^8+24*t^7+54*t^6+56*t^5+54*t^4+24*t^3+17*t^2+4*t+1;
p0 := t^6+4*t^5+15*t^4+16*t^3+15*t^2+4*t+1;
A21 := p2*x^2 - 2*(t^2+1)^2*p1*x + (t^2+1)^4*p0;
k21 := 64*t^4*(t+1)^2*(t^2+1)^3*(t^4+2*t^3+6*t^2+2*t+1)^3;
f := A21^2 - k21*x^3*(x-1)^2;
assert Degree(f) eq 5;
J := Jacobian(HyperellipticCurve(f));
D := elt< J | x, Evaluate(A21, 0) >;
assert 21*D eq J!0 and 3*D ne J!0 and 7*D ne J!0;
print "over Q(t): D = [x, A21(0)] has order 21";

// (2), (3) residue certificate at t = 1 mod 5
g, st := ResidueCertificate(f, [5], [1], [1, 6, -4, 6/11], 21, [21]);
print "FAMILY_21_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [21]";
quit;
