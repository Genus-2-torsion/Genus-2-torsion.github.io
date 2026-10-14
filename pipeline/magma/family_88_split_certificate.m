/* family_88_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have geometrically split Jacobian with J(Q)_tors exactly [8,8].

   Source of the construction.  BNSS's conjugation self-gluing of the Kubert X_1(8) curve
       E_t : y^2 + (1 - c)xy - by = x^3 - bx^2,  b = (2t - 1)(t - 1),  c = b/t,  P = (0,0) of order 8,
   glued to itself along the anti-isometry of E_t[2] that fixes the rational 2-torsion point and swaps
   the conjugate pair.  The halving condition is t(1 - t) = square, i.e. t = a^2/(a^2 + b^2)
   (Genus2TorsionNotebook, notes/claude_split_2266_1212_session_2026_08_13.md, section 2.1, and
   product/code/lane_1212.m, lane_1212_sym.m; there exact [8,8] is asserted for 59 values of t; that
   code was read, not run).  Here t = m^2/(m^2 + 1).

   The model (derived for this certificate).  With N = m^8 - 4m^6 - 2m^4 - 4m^2 + 1,
       F(u) = -N (u + 1)(16 m^8 (u + 1)^2 - N^2 u),     E'_m : Y^2 = F(u),
       C_m : y^2 = f(x) = F(x^2) = -N (x^2 + 1)(4m^4 x^2 - N x + 4m^4)(4m^4 x^2 + N x + 4m^4).
   F is palindromic (u^3 F(1/u) = F(u)), so
       phi1 : (x, y) -> (x^2, y)   and   phi2 : (x, y) -> (1/x^2, y/x^3)
   are two maps C_m -> E'_m of degree 2.  The induced gluing of E'[2] is (u,0) -> (1/u,0): it fixes
   u = -1 and swaps the two roots of 16m^8(u + 1)^2 - N^2 u (their product is 1), which is BNSS's
   conjugation self-gluing.  E'_m is isomorphic to E_t over Q(m) by an explicit affine map, checked below.

   Argument.
   (1) Splitting.  phi1 is a nonconstant morphism C_m -> E'_m.  Wherever f_m is squarefree of degree 6,
       F_m is a squarefree cubic with F_m(0) != 0, so E'_m is an elliptic curve and J(C_m) is isogenous
       over Q to E'_m x E'_m (phi1^* + phi2^* : E' x E' -> J is a (2,2)-isogeny).  So every member has
       a geometrically split Jacobian (in fact split over Q).
   (2) Containment over Q(m).  Let P' = (u_P, Y_P) be the image of P on E'.  Then
       A = phi1^*(P' - O) = [x^2 - u_P, Y_P]  and  B = phi2^*(P' - O) = [x^2 - 1/u_P, (Y_P/u_P) x].
       By Mumford arithmetic over Q(m): 8A = 8B = 0, 4A = 4B != 0.  A half D of A + B is found by
       solving the halving identity f - V^2 = c u_D^2 u_T (T = A + B, V = -v_T + u_T (L0 + L1 x)) over
       Q(m); the proof is the check 2D = A + B.  Then 8D = 0, 4D != 0 and 4D != 4A, so
       <A, D> = Z/8 x Z/8 (the 64 sums iA + jD are distinct).
   (3) Exactness on a residue class.  The gcd of #J(F_p) cannot certify [8,8] here.  E_t/<4P> has full
       rational 2-torsion and a point of order 8 over Q(m) (checked below), so 16 | #E'(F_p) and
       256 | #J(F_p) = #E'(F_p)^2 at every good p.  The certificate uses the group structure instead.
       Reduction J_m(Q)_tors -> J(F_p) is an injective homomorphism (p odd, good reduction), so
       J_m(Q)_tors embeds in every J(F_p).  For p-integral m = 2 mod 13 the curve has the reduction of
       m = 2 at 13, with J(F_13) = Z/16 x Z/16: so J_m(Q)_tors is a 2-group of 2-rank <= 2.  For m = 2
       mod 337, J(F_337) = Z/8 x Z/8 x Z/40 x Z/40, whose 2-primary part has exponent 8.  A group that
       embeds in both has order <= 2^(2*3) = 64, and it contains <A, D>.  So every 13- and 337-integral
       m = 2 mod 4381 gives J(Q)_tors = Z/8 x Z/8 exactly.  As a check on (1), the Frobenius polynomial
       of J at 13 and at 337 is the square of that of E'.
   (4) The G2-invariants of the members differ, so m -> M_2 is nonconstant with finite fibers, and the
       residue class gives infinitely many isomorphism classes.  The member m = 2 (t = 4/5) has
       TorsionSubgroup [8,8].
   Run (from the repository root):
       magma -b pipeline/magma/family_88_split_certificate.m > family_88_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

Fm<m> := FunctionField(Q); P<x> := PolynomialRing(Fm);
N := m^8 - 4*m^6 - 2*m^4 - 4*m^2 + 1;
F := -N*(x + 1)*(16*m^8*(x + 1)^2 - N^2*x);           // E' : Y^2 = F(u), u written as x
f := Evaluate(F, x^2);
assert f eq -N*(x^2 + 1)*(4*m^4*x^2 - N*x + 4*m^4)*(4*m^4*x^2 + N*x + 4*m^4);
assert Degree(F) eq 3 and Reverse(Coefficients(F)) eq Coefficients(F);
quad2 := 16*m^8*(x + 1)^2 - N^2*x;
assert #Factorization(quad2) eq 1 and Coefficient(quad2, 0) eq Coefficient(quad2, 2);
print "C_m: y^2 = f(x) = F(x^2); F is palindromic, so phi1(x,y) = (x^2, y) and phi2(x,y) = (1/x^2, y/x^3) map C_m to E': Y^2 = F(u)";
print "E'[2]: u = -1 and the two roots of the irreducible palindromic quadratic 16m^8(u+1)^2 - N^2 u (swapped by u -> 1/u)";

// the Kubert curve E_t on the conic t = m^2/(m^2+1), and the isomorphism E_t -> E'
tK := m^2/(m^2 + 1);
assert tK*(1 - tK) eq (m/(m^2 + 1))^2;
bK := (2*tK - 1)*(tK - 1); cK := bK/tK;
E := EllipticCurve([Fm| 1 - cK, -bK, -bK, 0, 0]);
P8 := E![0, 0];
assert 8*P8 eq E!0 and 4*P8 ne E!0;
a1, a2, a3, a4, a6 := Explode(aInvariants(E));
gE := x^3 + a2*x^2 + a4*x + a6 + (a1*x + a3)^2/4;      // Y^2 = gE(X) with Y = y + (a1 X + a3)/2
c1 := (-1/4*m^12 + 1/2*m^10 + 9/4*m^8 + 3*m^6 + 9/4*m^4 + 1/2*m^2 - 1/4)/m^4;
c0 := (-1/4*m^8 + m^6 + 1/2*m^4 - 1/4)/m^2;
kk := (1/2*m^22 - 5/2*m^20 - 9/2*m^18 + 21/2*m^16 + 45*m^14 + 79*m^12 + 79*m^10 + 45*m^8 + 21/2*m^6
       - 9/2*m^4 - 5/2*m^2 + 1/2)/m^2;
assert kk ne 0 and c1 ne 0 and kk^2*gE eq Evaluate(F, c1*x + c0);
print "E_t (Kubert X_1(8), t = m^2/(m^2+1), t(1-t) a square) -> E': (X, Y) -> (c1 X + c0, kk Y) is an isomorphism over Q(m)";
uP := c0; YP := kk*a3/2;                                // image of P = (0,0)
assert Evaluate(F, uP) eq YP^2;
printf "P' = (%o, %o)\n", uP, YP;

// (2) the marked (Z/8)^2 over Q(m)
J := Jacobian(HyperellipticCurve(f));
A := elt< J | x^2 - uP, YP >;
B := elt< J | x^2 - 1/uP, (YP/uP)*x >;
assert 8*A eq J!0 and 4*A ne J!0 and 8*B eq J!0 and 4*A eq 4*B;
print "over Q(m): A = phi1^*P' and B = phi2^*P' have order 8, and 4A = 4B";
R2<L0,L1> := PolynomialRing(Fm, 2);
PX<XX> := PolynomialRing(R2);
// rational halves D of a class T = [u_T, v_T] with deg u_T = 2: f - V^2 = c u_D^2 u_T, V = -v_T + u_T (L0 + L1 x)
function Halves(T)
    uT := T[1]; vT := T[2];
    assert Degree(uT) eq 2 and T[3] eq 2;
    V := -(PX!vT) + (PX!uT)*(L0 + L1*XX);
    G, rr := Quotrem(PX!f - V^2, PX!uT); assert rr eq 0;
    g := [ Coefficient(G, i) : i in [0..4] ];
    e1 := 8*g[5]^2*g[2] - 4*g[5]*g[4]*g[3] + g[4]^3;
    e2 := 64*g[5]^3*g[1] - (4*g[5]*g[3] - g[4]^2)^2;
    res := UnivariatePolynomial(Resultant(e1, e2, L0));
    out := [];
    for h in Factorization(res) do
        if Degree(h[1]) ne 1 then continue; end if;
        lam1 := -Coefficient(h[1], 0)/Coefficient(h[1], 1);
        q1 := UnivariatePolynomial(Evaluate(e1, [L0, R2!lam1]));
        q2 := UnivariatePolynomial(Evaluate(e2, [L0, R2!lam1]));
        for rt in Roots(GCD(q1, q2)) do
            Vx := -vT + uT*(rt[1] + lam1*x);
            Gx := (f - Vx^2) div uT;
            cc := LeadingCoefficient(Gx);
            if cc eq 0 or Degree(Gx) ne 4 then continue; end if;
            okk, uD := IsSquare(Gx/cc);
            if not okk then continue; end if;
            assert f - Vx^2 eq cc*uD^2*uT;
            Append(~out, elt< J | uD, Vx mod uD >);
        end for;
    end for;
    return out;
end function;
H := Halves(A + B);
assert #H ge 1;
D := H[1];
printf "D = %o\n", D;
assert 2*D eq A + B;
assert 8*D eq J!0 and 4*D ne J!0 and 4*D ne 4*A;
assert #{ i*A + j*D : i, j in [0..7] } eq 64;
print "over Q(m): 2D = A + B, and the 64 sums iA + jD are distinct: [8,8] is marked";

// why #J(F_p) cannot certify [8,8]: E_t/<4P> has Z/2 x Z/8 over Q(m)
E2, iso2 := IsogenyFromKernel(E, x - (4*P8)[1]);
b1, b2, b3, b4, b6 := Explode(aInvariants(E2));
assert &and[ Degree(h[1]) eq 1 : h in Factorization(x^3 + (b1^2/4 + b2)*x^2 + (b1*b3/2 + b4)*x + (b3^2/4 + b6)) ];
assert #DivisionPoints(iso2(P8), 2) ge 1 and 4*iso2(P8) eq E2!0 and 2*iso2(P8) ne E2!0;
print "E_t/<4P> has full rational 2-torsion and a point of order 8 over Q(m): 16 | #E'(F_p), 256 | #J(F_p) at every good p";

// (3) exactness from the group structures of J(F_13) and J(F_337) on the class m = 2 mod 4381
primes := [13, 337]; residues := [2, 2];
members := [Q| 2, 4383, -4379, 4383/4382];
odd := []; rk := []; ex := [];
for i in [1..#primes] do
    p := primes[i]; r0 := residues[i];
    fbar := ReductionClass(f, p, r0);
    for mm in members do
        assert GF(p)!mm eq GF(p)!r0;
        assert PolynomialRing(GF(p))!SpecPoly(f, mm) eq fbar;
    end for;
    n, chi := FrobData(fbar);
    Ap := AbelianGroup(Jacobian(HyperellipticCurve(fbar)));
    inv := Invariants(Ap);
    assert #Ap eq n;
    // E' mod p (monic model V^2 = U^3 + c2 U^2 + c1 c3 U + c0 c3^2) and its Frobenius polynomial
    Fp := GF(p); Pp := PolynomialRing(Fp);
    Fb := Pp![ Fp!Evaluate(Numerator(cf), r0)/Fp!Evaluate(Denominator(cf), r0) : cf in Coefficients(F) ];
    k3 := Coefficient(Fb, 3);
    Ebar := EllipticCurve(Pp![ Coefficient(Fb, 0)*k3^2, Coefficient(Fb, 1)*k3, Coefficient(Fb, 2), 1 ]);
    aE := p + 1 - #Ebar;
    assert chi eq (X^2 - aE*X + p)^2;
    v2 := Valuation(n, 2);
    Append(~odd, n div 2^v2);
    Append(~rk, #[ q : q in inv | IsEven(q) ]);
    Append(~ex, Max([ Valuation(q, 2) : q in inv ]));
    printf "p = %o, m = %o mod %o: reduction %o (squarefree, degree %o), J(F_%o) = %o (order %o), chi = %o = (%o)^2, #E'(F_%o) = %o\n",
        p, r0, p, fbar, Degree(fbar), p, inv, n, chi, X^2 - aE*X + p, p, #Ebar;
end for;
godd := GCD(odd); r2 := Min(rk); e2 := Min(ex);
printf "odd parts %o (gcd %o); 2-ranks %o (min %o); 2-primary exponents 2^%o (min 2^%o)\n", odd, godd, rk, r2, ex, e2;
assert godd eq 1 and r2 eq 2 and e2 eq 3;
print "J_m(Q)_tors embeds in both: a 2-group of 2-rank <= 2 and exponent <= 8, so of order <= 64 = #<A, D>";
inv0 := TorsionInvariants(SpecPoly(f, members[1]));
printf "sanity: m = %o (t = %o) has TorsionSubgroup %o\n", members[1], members[1]^2/(members[1]^2 + 1), inv0;
assert inv0 eq [8,8];

// (4) the moduli map is nonconstant
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, mm))) : mm in members } ge 2;
printf "G2-invariants of the members %o are not all equal\n", members;
print "FAMILY_88_SPLIT_CERTIFIED: infinitely many geometrically split Jacobians (split over Q, isogenous to E x E) with J(Q)_tors exactly [8,8]";
quit;
