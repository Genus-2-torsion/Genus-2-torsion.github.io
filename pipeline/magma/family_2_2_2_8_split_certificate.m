/* family_2_2_2_8_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E x E2_m) with J(Q)_tors exactly [2,2,2,8].
   Found 2026-10-19 (torsion_jac notes/g33ell_2026_10_19.md).

   Source.  Bielliptic curves with six rational Weierstrass points and one elliptic quotient fixed.  Let
   E = 210e2 (torsion Z/2 x Z/8) in the form Y^2 = g(u) = (u - e1)(u - e2)(u - e3), (e1, e2, e3) = (-36, 31/4, 28),
   Y = y + (a1 u + a3)/2.  A curve y^2 = c (x^2 - r1^2)(x^2 - r2^2)(x^2 - r3^2) has the degree-2 map
   (x, y) -> (x^2, y) onto Y^2 = c prod (u - r_i^2); this quotient is isomorphic to E when an affine map
   u -> mu u + k sends r_i^2 to e_i and c = mu (square).  The first condition is the conic
   r3^2 - r1^2 = lambda (r2^2 - r1^2), lambda = (e3 - e1)/(e2 - e1) = 256/175, with the rational point (1 : 1 : 1);
   the lines through it give
       r1 = m^2 - lambda,  r2 = m^2 - 2m + lambda,  r3 = -m^2 + 2 lambda m - lambda,
       C_m : y^2 = f(x) = c (x^2 - r1^2)(x^2 - r2^2)(x^2 - r3^2),  c = (e2 - e1)(r2^2 - r1^2),
   with phi1 : (x, y) -> (mu x^2 + k, mu y/s) onto E (s = r2^2 - r1^2, mu = (e2 - e1)/s, k = e1 - mu r1^2;
   certified by the identity s^2 g(mu x^2 + k) = mu^2 f(x) over Q(m)) and phi2 : (x, y) -> (1/x^2, y/x^3) onto
   E2_m : Y^2 = w^3 F(1/w), F(u) = c prod (u - r_i^2).  C_m is the (2,2)-gluing (Howe-Leprevost-Poonen 2000) of
   E and E2_m along (e_i, 0) -> (1/r_i^2, 0) (phi1^*(e_i, 0) = phi2^*(1/r_i^2, 0) = [x^2 - r_i^2, 0]), with all of
   J[2] rational.  The 2-Kummer images of E2_m on its 2-torsion equal those of E (checked over Q(m) in the note),
   so a 2-torsion point of E2_m is divisible by 2 over Q exactly when the corresponding point of E is;
   here E2_m(Q) contains Z/2 x Z/4 at every member, so #J(F_p) is always divisible by 128 = 2 #G and the group
   structures are used.

   Argument.
   (1) Splitting.  The identity holds over Q(m).  At every member of the class below f_m is squarefree of
       degree 6 (its reduction is), so s and mu are finite and nonzero, phi1 and phi2 are morphisms over Q of
       degree 2 onto E and E2_m, and J is isogenous over Q to E x E2_m (phi1^* + phi2^* is an isogeny with
       kernel of order 4).
   (2) The marked group over Q(m), by Mumford arithmetic: the 15 classes [(x - w)(x - w'), 0] for the pairs of
       the six rational Weierstrass points +-r_i, and 0, are 16 distinct classes killed by 2;
       A = phi1^*(P8) = [x^2 - (u(P8) - k)/mu, Y(P8) s/mu] for P8 = (4, -62) has order 8,
       and together they generate 64 classes: [2,2,2,8] is marked.
       The marked sections specialize injectively on smooth fibres, so [2,2,2,8] embeds in J_m(Q)_tors.
   (3) Exactness on a residue class.  For every 11- and 13-integral m = 8/5 (mod 143), C_m has the reductions of
       C_{8/5} at 11 and 13, with J(F_11) = [2,2,2,16] and
       J(F_13) = [2,2,8,8] (the Frobenius polynomials are chi_E chi_E2, as the splitting predicts).
       J_m(Q)_tors embeds in all of them, and the largest group that does is [2,2,2,8]: it is exactly [2,2,2,8].
   (4) The G2-invariants of the members differ, so m -> M_2 is nonconstant with finite fibres, and the class
       gives infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup [2,2,2,8] at m = 8/5 and m = 723/5; the member m = 8/5 is Q-isomorphic to the
       census curve 2.2.2.8.b, with E2_{8/5} = 21a1 (torsion Z/2 x Z/4).
       Small members (m = n/d, |n| <= 12, d <= 6) all have J(Q)_tors = [2,2,2,8] (see the note).
   Run (from the repository root):
       magma -b pipeline/magma/family_2_2_2_8_split_certificate.m > data/knowledge/sources/family_2_2_2_8_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

// The largest finite abelian group that embeds in every group with the given invariants: for each prime l,
// the sorted l-primary exponents of a subgroup are bounded componentwise by those of the group.
function StructureBound(invlist)
    ls := &meet[ Seqset(&cat[ PrimeDivisors(n) : n in inv | n gt 1 ]) : inv in invlist ];
    out := [Integers()| ];
    for l in ls do
        es := [ Reverse(Sort([ Valuation(n, l) : n in inv | n mod l eq 0 ])) : inv in invlist ];
        r := Min([ #e : e in es ]);
        out cat:= [ l^Min([ e[i] : e in es ]) : i in [1..r] ];
    end for;
    return Invariants(AbelianGroup(out));
end function;
function OrderUpTo(P, n)
    for k in [1..n] do if k*P eq Parent(P)!0 then return k; end if; end for; return -1;
end function;
// Y^2 = F(u), F a cubic with F(0) != 0, as a Weierstrass curve (monic model V^2 = U^3 + b2 U^2 + b1 k U + b0 k^2, U = k u, V = k Y)
function CubicCurve(F)
    k := Coefficient(F, 3);
    return EllipticCurve([0, Coefficient(F, 2), 0, Coefficient(F, 1)*k, Coefficient(F, 0)*k^2]);
end function;

// ---------------------------------------------------------------- (0) the fixed curve E
E := EllipticCurve([1,0,0,-1070,7812]);
assert Conductor(E) eq 210 and Invariants(TorsionSubgroup(E)) eq [2,8];
Pgen := E![4,-62]; ngen := 8;
assert Order(Pgen) eq ngen;
a1, a2, a3, a4, a6 := Explode(aInvariants(E));
u := X;                                                  // RX = Q[X] from family_cert_lib.m
g := u^3 + a2*u^2 + a4*u + a6 + (a1*u + a3)^2/4;          // E : Y^2 = g(u), Y = y + (a1 u + a3)/2
eps := [Q| -36, 31/4, 28]; e1, e2, e3 := Explode(eps);
assert g eq (u - e1)*(u - e2)*(u - e3);
printf "E = 210e2 = %o, torsion %o; Y^2 = g(u) = (u - e1)(u - e2)(u - e3) with (e1, e2, e3) = %o; Pgen = %o of order %o\n",
    aInvariants(E), Invariants(TorsionSubgroup(E)), eps, Pgen, ngen;

// ---------------------------------------------------------------- (1) the family over Q(m)
Fm<m> := FunctionField(Q); P<x> := PolynomialRing(Fm);
lam := (e3 - e1)/(e2 - e1);
r1 := m^2 - lam; r2 := m^2 - 2*m + lam; r3 := -m^2 + 2*lam*m - lam;      // lines through (1 : 1 : 1) on the conic below
assert r3^2 - r1^2 eq lam*(r2^2 - r1^2);
sQ := r2^2 - r1^2; mu := (e2 - e1)/sQ; kk := e1 - mu*r1^2; c := (e2 - e1)*sQ;
assert c eq mu*sQ^2;
assert [ (e - kk)/mu : e in eps ] eq [ r1^2, r2^2, r3^2 ];
F := c*(x - r1^2)*(x - r2^2)*(x - r3^2);
f := Evaluate(F, x^2);
assert f eq c*(x - r1)*(x + r1)*(x - r2)*(x + r2)*(x - r3)*(x + r3);
printf "(1) lambda = (e3 - e1)/(e2 - e1) = %o; (r1 : r2 : r3) = (%o : %o : %o) is on the conic r3^2 - r1^2 = lambda (r2^2 - r1^2)\n", lam, r1, r2, r3;
printf "(1) C_m : y^2 = f(x) = F(x^2), F(u) = c (u - r1^2)(u - r2^2)(u - r3^2), c = (e2 - e1)(r2^2 - r1^2) = %o\n", c;
gP := P![ Fm!cf : cf in Coefficients(g) ];
assert sQ^2*Evaluate(gP, mu*x^2 + kk) eq mu^2*f;
printf "(1) phi1 : (x, y) -> (mu x^2 + k, mu y/s) maps C_m onto E (Y^2 = g(u)): s^2 g(mu x^2 + k) = mu^2 f(x) over Q(m), with s = r2^2 - r1^2, mu = (e2 - e1)/s = %o, k = e1 - mu r1^2 = %o\n", mu, kk;
F2 := P!Reverse(Coefficients(F));                          // F2(w) = w^3 F(1/w)
assert Degree(F2) eq 3 and Coefficients(F2) eq Reverse(Coefficients(F));
print "(1) phi2 : (x, y) -> (1/x^2, y/x^3) maps C_m onto E2_m : Y^2 = F2(w) = w^3 F(1/w); so J(C_m) is isogenous over Q to E x E2_m";
print "    the six Weierstrass points (+-r_i, 0) of C_m are rational over Q(m)";

// ---------------------------------------------------------------- (2) the marked group over Q(m)
J := Jacobian(HyperellipticCurve(f));
lin := [x - r1, x + r1, x - r2, x + r2, x - r3, x + r3];
T2 := { J!0 } join { elt< J | lin[i]*lin[j], 0 > : i, j in [1..6] | i lt j };
assert #T2 eq 16 and &and[ 2*t eq J!0 : t in T2 ];
uP := Pgen[1]; YP := Pgen[2] + (a1*Pgen[1] + a3)/2;
A := elt< J | x^2 - (uP - kk)/mu, YP*sQ/mu >;              // phi1^*(Pgen) - phi1^*(O)
assert OrderUpTo(A, 64) eq ngen;
Gm := { i*A + t : i in [0..ngen-1], t in T2 };
printf "(2) over Q(m): the 15 pair classes [(x - w)(x - w'), 0] and 0 are 16 classes killed by 2 (J[2] is rational); A = phi1^*(Pgen) = %o has order %o; together they generate %o classes\n", A, OrderUpTo(A, 64), #Gm;
assert #Gm eq 64;
print "    so J(Q(m)) contains G = [2,2,2,8]";

// ---------------------------------------------------------------- (3) exactness on the residue class m = 8/5 mod [11, 13]
primes := [11, 13]; m0 := 8/5; M := &*primes;
members := [Q| m0, m0 + M, m0 - M, m0 + 2*M];
invs := [];
for p in primes do
    assert Integers()!Conductor(E) mod p ne 0;
    r0 := Integers()!(GF(p)!m0);
    fbar := ReductionClass(f, p, r0);
    for mm in members do
        assert (GF(p)!mm) eq (GF(p)!r0);
        assert PolynomialRing(GF(p))!SpecPoly(f, mm) eq fbar;
    end for;
    n, chi := FrobData(fbar);
    Ap := AbelianGroup(Jacobian(HyperellipticCurve(fbar)));
    assert #Ap eq n;
    Fp := GF(p); Pp := PolynomialRing(Fp);
    F2bar := Pp![ Fp!Evaluate(Numerator(cf), r0)/Fp!Evaluate(Denominator(cf), r0) : cf in Coefficients(F2) ];
    E2bar := CubicCurve(F2bar);
    aE := p + 1 - #ChangeRing(E, Fp); aE2 := p + 1 - #E2bar;
    assert chi eq (X^2 - aE*X + p)*(X^2 - aE2*X + p);
    Append(~invs, Invariants(Ap));
    printf "(3) p = %o, m = %o mod %o: reduction %o (squarefree, degree %o); J(F_p) = %o (order %o), chi = (%o)(%o): #E(F_p) = %o, #E2(F_p) = %o\n",
        p, r0, p, fbar, Degree(fbar), Invariants(Ap), n, X^2 - aE*X + p, X^2 - aE2*X + p, p + 1 - aE, p + 1 - aE2;
end for;
bnd := StructureBound(invs);
printf "(3) the largest group embedding in all J(F_p) has invariants %o\n", bnd;
assert bnd eq [2,2,2,8];
printf "(3) every p-integral m = %o mod %o (p in %o) has J(Q)_tors = G = [2,2,2,8]\n", m0, M, primes;

// ---------------------------------------------------------------- (4) infinitely many isomorphism classes, (5) sanity
G2s := {@ G2Invariants(HyperellipticCurve(SpecPoly(f, mm))) : mm in members @};
assert #G2s ge 2;
printf "(4) the G2-invariants of the members %o are not all equal: m -> M_2 is nonconstant, so the class gives infinitely many isomorphism classes\n", members;
for mm in members[1..2] do
    fm := SpecPoly(f, mm);
    invm := TorsionInvariants(fm);
    E2m := CubicCurve(PolynomialRing(Q)![ Evaluate(cf, mm) : cf in Coefficients(F2) ]);
    printf "(5) m = %o: TorsionSubgroup %o; E2_m = %o (minimal model), conductor %o, torsion %o\n", mm, invm, aInvariants(MinimalModel(E2m)), Conductor(E2m), Invariants(TorsionSubgroup(E2m));
    assert invm eq [2,2,2,8];
end for;
Cc := HyperellipticCurve(PolynomialRing(Q)![-45,42,30,-23,-9,3,1], PolynomialRing(Q)![0,1,1]);
assert IsIsomorphic(HyperellipticCurve(SpecPoly(f, 8/5)), Cc);
printf "(5) the member m = %o is Q-isomorphic to the census curve 2.2.2.8.b\n", 8/5;
printf "time %o s\n", Cputime(T0);
print "FAMILY_2_2_2_8_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E x E2_m) and J(Q)_tors exactly [2,2,2,8]";
quit;
