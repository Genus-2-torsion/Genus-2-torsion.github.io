/* family_2244_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,2,4,4].

   The record's family (census_infinity_record.md, section 2.5): the rational curve in A(2,2,4,4) of
   BNSS2026 (Proposition "rational curve in A(2,2,4,4)"; torsion_jac/code/claude_census_infty_certify.m,
   family 7):
       a = (s-2)(s^2-s+1/2)(s^2-3/2)(s^4+s^2+9/4),   b = (s^2-2s+1/2)(s^2-2s+9/2)(s^2-3/2)(s^2+1),
       c = -(s-2)(2s^2-2s+1)(s^4+s^2+9/4),            d = (s-2)(2s^2-2s+1)(2s^2-3)(s^2+1),
       C_s : y^2 = f_s = x(x+A)(x+B)(x+C)(x+D),  (A,B,C,D) = (a^2,b^2,c^2,d^2)   (degree 5).
   Argument.
   (1) Lemma 2244(b) of BNSS2026: with s_i the elementary symmetric polynomials in a, b, c, d,
       D0 = [x^2 - s2 x + s4, (s1 s2 - s3) x - s1 s4] satisfies 2 D0 = [x, 0].
       Lemma 2244(c): (C-A)(D-A), (C-B)(D-B), (C-A)(C-B), (D-A)(D-B) are squares in Q(s) (checked), and
       then [R] = [(x+A)(x+B), 0] is divisible by 2.  An explicit half: with F = x(x+C)(x+D), find
       U = k x^2 + U1 x + U0 and l = k x + l0 in Q(s)[x] with F = R l^2 - U^2.  Then f - (R l)^2 = -k^2 R u^2
       (u = U/k), so div(y - R l) = (-A,0) + (-B,0) + 2(P1 + P2) - 6 oo, and D1 = [u, R l mod u] satisfies
       2 D1 = [R].  Evaluated at the five roots of f the identity becomes five linear equations in
       (U0, U1, k, l0), one system for each choice of signs of the square roots; they are solved over
       Q(s) and a solution satisfying the identity is kept.  Mumford checks over Q(s): 2 D0 = [x, 0],
       2 D1 = [R], and the 64 sums i D0 + j D1 + k [x+A, 0] + m [x+C, 0] (i, j in Z/4; k, m in Z/2) are
       distinct, so [2,2,4,4] is marked.
   (2) Residue certificate (family_cert_lib.m) at s = 3 mod 37 and s = 3 mod 43, i.e. every 37- and
       43-integral s = 3 mod 1591: #J(F_37) = 1216, #J(F_43) = 1856, gcd 64 = #[2,2,4,4]; the Frobenius
       at 37 is strict.  The record's fiber s = 3 lies in the class and has TorsionSubgroup [2,2,4,4].
   (3) The G2-invariants of the members differ, so s -> M_2 is nonconstant.
   (The library's log lines call the family parameter t; here it is s.)
   Run (from the repository root):
       magma -b pipeline/magma/family_2244_certificate.m > data/knowledge/sources/family_2244_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
Fs<s> := FunctionField(Q); P<x> := PolynomialRing(Fs);
a := (s-2)*(s^2-s+1/2)*(s^2-3/2)*(s^4+s^2+9/4);
b := (s^2-2*s+1/2)*(s^2-2*s+9/2)*(s^2-3/2)*(s^2+1);
c := -(s-2)*(2*s^2-2*s+1)*(s^4+s^2+9/4);
d := (s-2)*(2*s^2-2*s+1)*(2*s^2-3)*(s^2+1);
A := a^2; B := b^2; C := c^2; D := d^2;
f := x*(x+A)*(x+B)*(x+C)*(x+D);
assert Degree(f) eq 5;
J := Jacobian(HyperellipticCurve(f));
// (1a) the half of [x, 0] (Lemma 2244(b))
aa := [a, b, c, d];
s1 := &+aa; s2 := &+[ aa[i]*aa[j] : i, j in [1..4] | i lt j ];
s3 := &+[ aa[i]*aa[j]*aa[k] : i, j, k in [1..4] | i lt j and j lt k ]; s4 := &*aa;
D0 := elt< J | x^2 - s2*x + s4, (s1*s2 - s3)*x - s1*s4 >;
T0 := elt< J | x, 0 >;
assert T0 ne J!0 and 2*D0 eq T0;
print "over Q(s): D0 = [x^2 - s2 x + s4, (s1 s2 - s3) x - s1 s4] has 2 D0 = [x, 0] != 0";
// (1b) the half of [(x+A)(x+B), 0] (Lemma 2244(c))
okA, yA := IsSquare((C-A)*(D-A)); okB, yB := IsSquare((C-B)*(D-B));
okC, wC := IsSquare((C-A)*(C-B)); okD, wD := IsSquare((D-A)*(D-B));
assert okA and okB and okC and okD;
print "over Q(s): (C-A)(D-A), (C-B)(D-B), (C-A)(C-B), (D-A)(D-B) are squares";
R := (x+A)*(x+B); F := x*(x+C)*(x+D);
assert f eq R*F;
sols := [];
for e2, e3, e4, e5 in [1, -1] do
    // U(-A) = a yA, U(-B) = e2 b yB, U(0) = e3 a b l(0), U(-C) = e4 wC l(-C), U(-D) = e5 wD l(-D)
    rows := [ [Fs| 1, -A, A^2, 0], [Fs| 1, -B, B^2, 0], [Fs| 1, 0, 0, -e3*a*b],
              [Fs| 1, -C, C^2 + e4*wC*C, -e4*wC], [Fs| 1, -D, D^2 + e5*wD*D, -e5*wD] ];
    rhs := [Fs| a*yA, e2*b*yB, 0, 0, 0];
    ok, sol := IsConsistent(Transpose(Matrix(Fs, 5, 4, &cat rows)), Vector(Fs, rhs));
    if ok then
        U0, U1, k, l0 := Explode(Eltseq(sol));
        if k ne 0 then
            U := k*x^2 + U1*x + U0; l := k*x + l0;
            if F eq R*l^2 - U^2 then Append(~sols, <U/k, l>); end if;
        end if;
    end if;
end for;
assert #sols ge 1;
u1 := sols[1][1]; l := sols[1][2];
assert F eq R*l^2 - (LeadingCoefficient(l)*u1)^2;
D1 := elt< J | u1, (R*l) mod u1 >;
TR := elt< J | R, 0 >;
assert TR ne J!0 and 2*D1 eq TR;
printf "halving identity F = R l^2 - (k u)^2 with l = %o\n", l;
printf "D1 = %o\n", D1;
print "over Q(s): 2 D1 = [(x+A)(x+B), 0] != 0";
TA := elt< J | x + A, 0 >; TC := elt< J | x + C, 0 >;
assert #{ i*D0 + j*D1 + k*TA + m*TC : i, j in [0..3], k, m in [0..1] } eq 64;
print "over Q(s): the 64 sums i D0 + j D1 + k [x+A, 0] + m [x+C, 0] are distinct: [2,2,4,4] is marked";
// (2) + (3)
g, st := ResidueCertificate(f, [37, 43], [3, 3], [Q| 3, 1594, -1588], 64, [2,2,4,4]);
print "FAMILY_2244_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,2,4,4]";
quit;
