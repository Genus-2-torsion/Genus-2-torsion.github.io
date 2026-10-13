/* family_2214_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,2,14].

   The contact-7 chart h = 1 - (7/2)x + a x^2 + b x^3, f = (h^2 + (x-1)^7)/x^2: x^2 f - h^2 = (x-1)^7 gives
   div(xy - h) = 7P - 7W with P = (1, h(1)), so [P - oo] has order 7.  With
       V(z) = 2z^5 + 4z^4 + 6z^3 + 8z^2 + 10z + 5,  A(z) = V(z)/(2(z+1)^2),
       b = (A(u) - A(s1))/(s1^2 - u^2),  a = A(s1) - b(1 - s1^2),
   f has the rational roots 1 - s1^2 and 1 - u^2 (torsion_jac/notes/contact7_two_root_surface_2026_07_11.md).
   A third root 1 - s3^2 is rational exactly on the collinearity locus G(s1, u, s3) = 0 of the points
   (1 - z^2, A(z)), z = s1, u, s3.  Its fiber over s1 = -1/2 is a plane sextic C of genus 1, birational to
   E : y^2 = x^3 - 12x (conductor 2^5 3^2, rank 1, P = (-2, 4)).  (Found in the notes review of 2026-10-13.)

   Argument.
   (1) Over Q(u), s1 = -1/2: [P - oo] has order 7 and the roots 3/4, 1 - u^2 give (Z/2)^2, so [2,14] is
       marked.  Over the function field of C, the third root 1 - s3^2 is a root of f_u: [2,2,14] is marked
       over Q(C).
   (2) A tuple (N1 : N2 : N3) of integral polynomials representing the birational map E -> C, with N3(P)
       a unit at 17 and 29.  For R = (1 + n k)P with n = lcm(ord(P mod 17), ord(P mod 29)), R reduces to P
       at 17 and 29.  So (u, s3) = (N1/N3, N2/N3)(R) is 17- and 29-integral and congruent to its value at P.
   (3) Residue certificate on f over Q(u) at u = u(P) mod 17 and mod 29 (family_cert_lib.m):
       #J(F_17) = 392 and #J(F_29) = 1120, gcd 56 = #[2,2,14], strict at 17.  The three roots stay distinct
       mod 17, so for every k >= 1 the curve C_{u(R)} has J(Q)_tors exactly [2,2,14] and a geometrically
       simple Jacobian.
   (4) u is a nonconstant function on E, so the u(R) are infinitely many distinct values.  The moduli map
       u -> M_2 is nonconstant (G2 at u(P) and u((1+n)P) differ), so these curves are infinitely many
       isomorphism classes.
   Run:  magma -b pipeline/magma/family_2214_certificate.m > data/knowledge/sources/family_2214_certificate.log */
SetColumns(0);
load "pipeline/magma/family_cert_lib.m";
V := func< z | 2*z^5 + 4*z^4 + 6*z^3 + 8*z^2 + 10*z + 5 >;
Az := func< z | V(z)/(2*(z+1)^2) >;
s1 := -1/2;

// (1) the family over Q(u)
Fu<u> := FunctionField(Q); P<x> := PolynomialRing(Fu);
b := (Az(u) - Az(Fu!s1))/(s1^2 - u^2);
a := Az(Fu!s1) - b*(1 - s1^2);
h := 1 - (7/2)*x + a*x^2 + b*x^3;
num := h^2 + (x-1)^7;
f := num div x^2; assert f*x^2 eq num and Degree(f) eq 5;
assert Evaluate(f, 1 - s1^2) eq 0 and Evaluate(f, 1 - u^2) eq 0;
J := Jacobian(HyperellipticCurve(f));
D7 := elt< J | x - 1, Evaluate(h, 1) >;
assert D7 ne J!0 and 7*D7 eq J!0;
print "over Q(u): [P - oo] has order 7 and f has the roots 3/4, 1 - u^2: [2,14] is marked";

// the collinearity locus and its fiber C over s1 = -1/2
P3<z1,z2,z3> := PolynomialRing(Q, 3);
M := Matrix(P3, 3, 3, [ 2*(z1+1)^2, 2*(z1+1)^2*(1-z1^2), V(z1),
                        2*(z2+1)^2, 2*(z2+1)^2*(1-z2^2), V(z2),
                        2*(z3+1)^2, 2*(z3+1)^2*(1-z3^2), V(z3) ]);
facs := Factorization(Determinant(M));
G := [ g[1] : g in facs | TotalDegree(g[1]) gt 2 ];
assert #G eq 1; G := G[1];
R3<XX,YY,WW> := PolynomialRing(Q, 3);
gh := Homogenization(Evaluate(G, [R3!s1, XX, YY]), WW);
PP := ProjectiveSpace(R3);
C := Curve(PP, gh);
assert Genus(C) eq 1;
// over the function field of C, 1 - s3^2 is a root of f_u (u = X/W, s3 = Y/W)
KC := FunctionField(C);
uC := KC!(XX/WW); s3C := KC!(YY/WW);
fC := &+[ Evaluate(Coefficient(f, i), uC)*(1 - s3C^2)^i : i in [0..5] ];
assert fC eq 0;
print "over Q(C): 1 - s3^2 is a third root of f_u, so [2,2,14] is marked over Q(C)";

// (2) E and the map E -> C
pt0 := C![-10, -10/7, 1];
Eg, phi := EllipticCurve(C, Places(pt0)[1]);
E, iso := MinimalModel(Eg);
assert aInvariants(E) eq [0, 0, 0, -12, 0];
r1, r2 := RankBounds(E); assert r1 eq 1 and r2 eq 1;
Pt := E![-2, 4]; assert Order(Pt) eq 0;
ok, phinv := IsInvertible(phi); assert ok;   // Eg -> C
psi := Inverse(iso)*phinv;                   // E -> Eg -> C
alts := [ DefiningPolynomials(psi) ];
try alts := AllDefiningPolynomials(psi); catch e; end try;
primes := [17, 29];
good := [];
PZ3 := PolynomialRing(Z, 3);
for tup in alts do
    // clear denominators: integral polynomials
    k := LCM([ Denominator(c) : c in &cat[ Coefficients(g) : g in tup ] ]);
    T := [ PZ3!(k*g) : g in tup ];
    vals := [ Evaluate(g, [Z!(-2), Z!4, Z!1]) : g in T ];
    if &and[ vals[3] mod p ne 0 : p in primes ] then Append(~good, T); end if;
end for;
assert #good ge 1;
T := good[1];
vals := [ Q!Evaluate(g, [-2, 4, 1]) : g in T ];
u0 := vals[1]/vals[3]; s30 := vals[2]/vals[3];
assert psi(Pt) eq C![u0, s30, 1];
printf "E = %o (rank 1), P = %o; psi(P) = (u, s3) = (%o, %o)\n", aInvariants(E), Pt, u0, s30;
n := LCM([ Order(ChangeRing(E, GF(p))!Pt) : p in primes ]);
printf "ord(P mod 17), ord(P mod 29) = %o; n = %o\n", [ Order(ChangeRing(E, GF(p))!Pt) : p in primes ], n;
printf "for R = (1 + n k)P, R = P mod 17 and mod 29, so u(R) = u(P) and s3(R) = s3(P) mod 17 and mod 29\n";
// sanity at k = 1
R1 := (1 + n)*Pt;
vR := [ Q!Evaluate(g, [R1[1], R1[2], R1[3]]) : g in T ];
u1 := vR[1]/vR[3]; s31 := vR[2]/vR[3];
assert psi(R1) eq C![u1, s31, 1];
for p in primes do assert GF(p)!u1 eq GF(p)!u0 and GF(p)!s31 eq GF(p)!s30; end for;
printf "sanity: u((1+n)P) has %o digits and = u(P) mod 17, 29\n", #Sprint(Numerator(u1));

// (3) residue certificate and distinct roots mod 17
roots0 := [ 1 - s1^2, 1 - u0^2, 1 - s30^2 ];
assert #{ GF(17)!r : r in roots0 } eq 3;
g, st := ResidueCertificate(f, primes, [ Z!(GF(p)!u0) : p in primes ], [u0, u1], 56, [2,2,14]);
print "FAMILY_2214_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,2,14]";
quit;
