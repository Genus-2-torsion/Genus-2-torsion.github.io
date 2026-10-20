/* leprevost_3_9_split.m -- every member of Leprevost's [3,9] family C_t is geometrically split.

   Leprevost 1995 (J. Theor. Nombres Bordeaux 7, pp. 302-304) gives the family
       C_t : y^2 = f_t(x) = x^6 + 6(t-1)x^5 + 3(t-1)(3t-5)x^4 - (18t^2-34t+18)x^3 + 3(t-1)(5t-3)x^2 - 6t(t-1)x + t^2
   whose Jacobian contains [3,9] over Q(t).  He shows that J(C_3) is simple over Q, and writes that he could not
   decide whether J(C_t) is absolutely simple (p. 303-304).  Platonov-Zhgun-Petrunin 2013 (Theorem 5) show that J(C_3)
   is isogenous over Qbar to a product of elliptic curves.  Here this is proved for the whole family.

   (0) Transcription: f_3 is C_3 as reprinted by PZP2013 (p. 321); members have torsion [3,9].
   (1) Over K(t), K = Q(sqrt -3): f_t = q q' with conjugate cubics q, q'.  With p = x^2 - x and
       E(u) := Res_y(q'(y), u q(y) - p(y)) (a cubic in u), the identity q^3 E(p/q) = c q' hh^2 holds with an explicit
       cubic hh and c in K(t)^*.  So (x, y) -> (u, v) = (p/q, y hh/q^2) maps C_t onto the curve c v^2 = E(u), of degree
       3 (gcd(p, q) = 1), and E has three distinct roots: the target has genus 1, and J(C_t) is isogenous over K(t) to
       a product of two elliptic curves.
   (2) Specialisation: f_t is squarefree exactly for t not in {0, 1} (and t^2 - t + 1 != 0).  The identity is an
       identity in K[t][x]; the conditions gcd(p, q) = 1, deg E = 3, disc E != 0, c != 0 fail only at t in {0, 1} or
       at irrational t.  So for every rational t != 0, 1 the curve C_t covers an elliptic curve over K with degree 3.
   (3) Checks of (1) at members: for t = 3, 5, -1/2 and the primes p = 1 (mod 3) below 200 of good reduction, the
       L-polynomial of the target elliptic curve mod p divides that of C_t mod p.
   Run (from the repository root):
       magma -b pipeline/magma/leprevost_3_9_split.m > data/knowledge/sources/leprevost_3_9_split.log */
SetColumns(0); SetMemoryLimit(16*10^9);
Q := Rationals(); Z := Integers();

// (0) transcription
Qx<z> := PolynomialRing(Q);
fq := func< t0 | z^6+6*(t0-1)*z^5+3*(t0-1)*(3*t0-5)*z^4-(18*t0^2-34*t0+18)*z^3+3*(t0-1)*(5*t0-3)*z^2-6*t0*(t0-1)*z+t0^2 >;
assert fq(3) eq z^6 + 12*z^5 + 24*z^4 - 78*z^3 + 72*z^2 - 36*z + 9;
print "(0) f_3 = x^6 + 12x^5 + 24x^4 - 78x^3 + 72x^2 - 36x + 9 is C_3 as reprinted by PZP2013, p. 321";
for t0 in [Q | 3, 5, -1/2] do
    f0 := fq(t0); d := LCM([Z | Denominator(c) : c in Coefficients(f0)]);
    T := Invariants(TorsionSubgroup(Jacobian(HyperellipticCurve(d^2*f0))));
    assert T eq [3, 9];
    printf "(0) t = %o: J(Q)_tors = %o\n", t0, T;
end for;

// (1) the cover over K(t)
K<w> := QuadraticField(-3);
FK<t> := FunctionField(K); PK<x> := PolynomialRing(FK);
f := x^6+6*(t-1)*x^5+3*(t-1)*(3*t-5)*x^4-(18*t^2-34*t+18)*x^3+3*(t-1)*(5*t-3)*x^2-6*t*(t-1)*x+t^2;
q  := x^3 + (3*t - w - 3)*x^2 + ((-w - 3)*t + (3*w + 3)/2)*x + (w + 1)/2*t;
qc := x^3 + (3*t + w - 3)*x^2 + ((w - 3)*t + (-3*w + 3)/2)*x + (-w + 1)/2*t;
assert f eq q*qc;
p  := x^2 - x;
hh := x^3 + (-w - 3)/2*x^2 + w*t*x + (-w + 1)/2*t;
PU<u> := PolynomialRing(FK : Global := false);
P2<yy> := PolynomialRing(PU);
E := Resultant(P2![PU | Coefficient(qc, k) : k in [0..3]], u*P2![PU | Coefficient(q, k) : k in [0..3]] - P2![PU | Coefficient(p, k) : k in [0..2]]);
assert Degree(E) eq 3;
lhs := &+[PK | Coefficient(E, k)*p^k*q^(3 - k) : k in [0..3]];     // q^3 E(p/q)
c := LeadingCoefficient(lhs);
assert lhs eq c*qc*hh^2;
assert c ne 0;
assert GCD(p, q) eq 1;
dE := Discriminant(E);
assert dE ne 0;
printf "(1) over Q(sqrt -3)(t): f_t = q q', and q^3 E(p/q) = c q' hh^2 with p = x^2 - x, c = %o\n", c;
printf "(1) q = %o\n(1) hh = %o\n(1) E(u) = %o\n", q, hh, E;
printf "(1) gcd(p, q) = 1, deg E = 3, disc E = %o != 0: (x, y) -> (p/q, y hh/q^2) maps C_t with degree 3 onto the genus-1 curve c v^2 = E(u)\n", Factorization(Numerator(dE));
jE := jInvariant(EllipticCurve(HyperellipticCurve(c^(-1)*E)));
printf "(1) j-invariant of the target: %o\n", jE;

// (2) the exceptional parameters
PZt<tt> := PolynomialRing(K);
function Polys(a)   // numerator and denominator of an element of K(t) as polynomials over K
    return [PZt | Evaluate(Numerator(a), tt), Evaluate(Denominator(a), tt)];
end function;
bad := Polys(c) cat Polys(LeadingCoefficient(E)) cat Polys(dE) cat Polys(Resultant(p, q));
for k in [0..3] do bad cat:= [Polys(Coefficient(E, k))[2]]; end for;
for k in [0..3] do bad cat:= [Polys(Coefficient(q, k))[2], Polys(Coefficient(hh, k))[2]]; end for;
rts := {Q | };
for b in bad do
    if b eq 0 then continue; end if;
    for r in Roots(b) do ok, r0 := IsCoercible(Q, r[1]); if ok then Include(~rts, r0); end if; end for;
end for;
assert rts subset {Q | 0, 1};
printf "(2) disc f_t = %o; rational t where the cover may degenerate: %o (where f_t is not squarefree)\n", Factorization(Numerator(Discriminant(f))), rts;
print "(2) so for every rational t != 0, 1 the Jacobian of C_t is isogenous over Q(sqrt -3) to a product of two elliptic curves";

// (3) numerical check at members
for t0 in [Q | 3, 5, -1/2] do
    f0 := fq(t0);
    nprimes := 0;
    for pr in [r : r in PrimesInInterval(7, 200) | r mod 3 eq 1] do
        if Valuation(Numerator(Discriminant(f0)), pr) gt 0 or Valuation(Denominator(t0), pr) gt 0 or Valuation(Numerator(t0), pr) gt 0 then continue; end if;
        Fp := GF(pr); s := Sqrt(Fp!(-3));
        ev := func< a | Evaluate(Numerator(a), t0)/Evaluate(Denominator(a), t0) >;
        // reduce the target c v^2 = E(u) mod pr with w -> s
        Ecoeffs := [ ev(Coefficient(E, k)) : k in [0..3] ]; c0 := ev(c);
        Ered := [ Fp!(Eltseq(a)[1]) + Fp!(Eltseq(a)[2])*s : a in Ecoeffs ];
        cred := Fp!(Eltseq(c0)[1]) + Fp!(Eltseq(c0)[2])*s;
        Pp<X> := PolynomialRing(Fp);
        Ep := Pp!Ered;
        if cred eq 0 or Discriminant(Ep) eq 0 or Degree(Ep) ne 3 then continue; end if;
        QT := PolynomialRing(Q);
        LE := QT!LPolynomial(HyperellipticCurve(cred^(-1)*Ep));
        LC := QT!LPolynomial(HyperellipticCurve(PolynomialRing(Fp)!f0));
        assert LC mod LE eq 0;
        nprimes +:= 1;
    end for;
    assert nprimes ge 5;
    printf "(3) t = %o: at %o primes p = 1 (mod 3) below 200 the L-polynomial of the target divides that of C_t mod p\n", t0, nprimes;
end for;
print "LEPREVOST_3_9_SPLIT_DONE: every member C_t (t rational, t != 0, 1) of Leprevost's [3,9] family is geometrically split";
quit;
