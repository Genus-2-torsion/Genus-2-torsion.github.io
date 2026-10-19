/* family_2_2_24_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E x F_u) with J(Q)_tors exactly [2,2,24].
   Found 2026-10-19 (torsion_jac notes/hlpEC_2026_10_19.md).

   Source.  HLP 2000, section 3.7 (the [2,2,24] entry of Table 1): the (2,2)-gluing (HLP Prop. 4) of
       E = E_{2,6}(241/81)   (Kubert X_1(2,6) curve, HLP Table 3, torsion Z/2 x Z/6, conductor 2.3.5.11.41.61)
   with the Kubert X_1(2,8) curves F_u = E_{2,8}(u) (torsion Z/2 x Z/8, P8 = (0,0)), along the isomorphism
   psi(T1) = T2', psi(T2) = T3', psi(T3) = T1' of the 2-torsion points of HLP Table 5 (T1' = 4 P8).  HLP's gain
   condition 2(t - 3) = -y^2, (t + 3)(t - 5) = (8u^2 - 1)(8u^2 + 8u + 1) z^2 is taken at y = 2/9 (t = 241/81), where
   it becomes the genus-1 curve
       D : s^2 = dpol(u) = -41 (8u^2 - 1)(8u^2 + 8u + 1),
   with Jacobian EM : y^2 = x^3 - 1681x (conductor 53792; HLP: rank 2).  On D a point of E x F_u that is not
   rational maps to a new rational 2-torsion point of the Jacobian (HLP Prop. 12).
   Model.  With f(z) = prod (z - x(T_i)) (E : y^2 = f(z), z = x, y' = y + (a1 x + a3)/2) and phi in Q(u)[z] the quadratic
   with phi(x(T_i)) = x(T'_psi(i)), phi = c0 + c1 z + c2 z^2, k = c1 + c2 (x(T1) + x(T2) + x(T3)), the HLP gluing is,
   up to a square (closed form of Prop. 4, see the note),
       C_u : y^2 = h(x) = c2 Fh(k x^2 - 1, c2 x^2),   Fh(X, Z) = Z^3 f(X/Z).
   phi2 : (x, y) -> (1/x^2, y/x^3) maps C_u onto Y^2 = c2^4 f((k - U)/c2), which is E, and phi1 : (x, y) -> (x^2, y)
   onto Ea : Y^2 = Ha(U) (h = Ha(x^2)), checked to be F_u over Q(u).

   Argument.
   (1) Splitting.  phi1^* + phi2^* : Ea x Eb -> J is an isogeny over Q wherever h is squarefree of degree 6.
   (2) The marked group.  Over Q(u): phi1^*, phi2^* are injective (prime degree, C not of genus 1);
       D6 = phi2^*(P6) (P6 = (0,0) of E, order 6) and D8 = phi1^*(P8) have orders 6 and 8, D6 + D8 has order 24, and
       the 2-torsion classes [x^2 - U_i, 0] (U_i = 1/(k - c2 x(T_i)), the roots of Ha) are not all in <D6 + D8>:
       so <D6 + D8, [x^2 - U_i, 0]> = [2,24] (the image of E(Q)_tors x F_u(Q)_tors).  The gain: U_2 is a square in
       Q(u) and U_3 is dpol(u) times a square, so on D both sqrt(U_2) and sqrt(U_3) are rational and
       T = [(x - sqrt U_2)(x - sqrt U_3), 0] is a rational 2-torsion class; it corresponds to a Galois-stable pair of
       Weierstrass points different from the pairs {+-sqrt U_i}, so T is not in the [2,24] (whose 2-torsion is
       {0, [x^2 - U_i, 0]}).  Hence [2,2,24] embeds in J(Q) at every point of D where h is squarefree of degree 6.
       Mumford arithmetic at the members below confirms it.
   (3) Exactness on a residue class.  #J(F_p) = #E(F_p) #F_u(F_p) is divisible by 192, so the group structure is
       used.  For every 13- and 47-integral u = 1/3 (mod 13*47), C_u has the reductions of C_{1/3} at 13 and 47, with
       J(F_13) = [2,2,48] (2-primary part [2,2,16], 3-primary part Z/3) and J(F_47) = [2,2,24,24] (2-primary
       exponent 8).  The largest group embedding in both is [2,2,24], so J(Q)_tors is exactly [2,2,24].
   (4) Infinitely many members.  The point B over u = -5/6 is sent to O on EM, and u = NX/NZ in the coordinates of EM;
       the member u = 1/3 lies over R0, and P = (-9, -120) has infinite order.  With m = lcm of the orders of P mod 13
       and 47, R_k = R0 + m k P reduces to R0 mod 13 and 47, NZ(R0) is a unit there, so u(R_k) is integral and = 1/3
       at both primes: each R_k gives a member of the class.  The R_k are distinct and u has degree 2 on D, so there
       are infinitely many u; G2-invariants at two points of D differ, so u -> M_2 is nonconstant on D with finite
       fibres: infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup at u = 1/3 and u = 5/27 is [2,2,24].  (u = 1/3, -5/6, -1/7, -5/14 give the same curve.)
       The member u = 1/3 (HLP's point) is the census curve 2.2.24.a (equal G2-invariants, checked below).
   Run (from the repository root):
       magma -b pipeline/magma/family_2_2_24_split_certificate.m > data/knowledge/sources/family_2_2_24_split_certificate.log */
SetColumns(0);
SetQuitOnError(true);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

function Kubert(b, c) return EllipticCurve([1 - c, -b, -b, 0, 0]); end function;
function MonEC(G)
    cs := [ Coefficient(G, i) : i in [0..3] ];
    return EllipticCurve([0, cs[3], 0, cs[2]*cs[4], cs[1]*cs[4]^2]), cs[4];
end function;
function OrderUpTo(P, n)
    for k in [1..n] do if k*P eq Parent(P)!0 then return k; end if; end for; return -1;
end function;
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

// ---------------------------------------------------------------- (0) E = E_{2,6}(241/81)
tE := Q!241/81;
E := Kubert((-2*tE^3 + 14*tE^2 - 22*tE + 10)/((tE + 3)^2*(tE - 3)^2), (-2*tE + 10)/((tE + 3)*(tE - 3)));
assert Invariants(TorsionSubgroup(E)) eq [2, 6] and Conductor(E) eq 2*3*5*11*41*61;
P6 := E![0, 0]; assert Order(P6) eq 6;
aE := aInvariants(E);
al := [ (-2*tE + 10)/((tE + 3)*(tE - 3)), (-tE^3 + 7*tE^2 - 11*tE + 5)/(4*(tE + 3)*(tE - 3)^2), (-2*tE^2 + 4*tE - 2)/((tE + 3)^2*(tE - 3)) ];
fE := &*[ X - a : a in al ];
b2, b4, b6 := Explode(bInvariants(E)[1..3]);
assert fE eq X^3 + b2/4*X^2 + b4/2*X + b6/4;
zP := 0; wP := aE[3]/2; assert Evaluate(fE, zP) eq wP^2;
printf "(0) E = E_{2,6}(241/81) = %o, conductor %o, torsion [2,6]; x(T1), x(T2), x(T3) = %o\n", aE, Factorization(Conductor(E)), al;

// ---------------------------------------------------------------- (1) the family over Q(u)
Fu<u> := FunctionField(Q); PZ<z> := PolynomialRing(Fu);
Fuu := Kubert((16*u^3 + 16*u^2 + 6*u + 1)/(8*u^2 - 1)^2, (16*u^3 + 16*u^2 + 6*u + 1)/(2*u*(4*u + 1)*(8*u^2 - 1)));
P8 := Fuu![0, 0]; assert OrderUpTo(P8, 20) eq 8;
be := [ (16*u^3 + 12*u^2 + 2*u)/(8*u^2 - 1)^2, (32*u^3 + 24*u^2 + 8*u + 1)/(16*u^2*(8*u^2 - 1)), (-32*u^4 - 32*u^3 - 12*u^2 - 2*u)/((4*u + 1)^2*(8*u^2 - 1)) ];
B2, B4, B6 := Explode(bInvariants(Fuu)[1..3]);
gF := z^3 + B2/4*z^2 + B4/2*z + B6/4;
assert gF eq &*[ z - bb : bb in be ];
assert (4*P8)[1] eq be[1];                                          // T1' = 4 P8
sig := [2, 3, 1];                                                   // psi(T_i) = T'_sig(i)
phi := &+[ be[sig[i]]*&*[ (z - al[j])/(al[i] - al[j]) : j in [1..3] | j ne i ] : i in [1..3] ];
assert &and[ Evaluate(phi, al[i]) eq be[sig[i]] : i in [1..3] ];
c2 := Coefficient(phi, 2); c1 := Coefficient(phi, 1); k := c1 + c2*&+al;
P<x> := PolynomialRing(Fu);
fEu := P![ Fu!c : c in Coefficients(fE) ];
Fh := func< XX, ZZ | &+[ Coefficient(fEu, i)*XX^i*ZZ^(3 - i) : i in [0..3] ] >;
f := c2*Fh(k*x^2 - 1, c2*x^2);
assert Degree(f) eq 6 and IsSquarefree(f) and &and[ Coefficient(f, i) eq 0 : i in [1, 3, 5] ];
Ha := P![ Coefficient(f, 2*i) : i in [0..3] ];
Hb := P![ Coefficient(f, 6 - 2*i) : i in [0..3] ];
assert Hb eq c2^4*Evaluate(fEu, (k - x)/c2);
Ea, la := MonEC(Ha);
c4a, c6a := Explode(cInvariants(Ea)); c4f, c6f := Explode(cInvariants(Fuu));
lam := (c6a/c6f)/(c4a/c4f);
assert c4a eq lam^2*c4f and c6a eq lam^3*c6f and IsSquare(lam);
printf "(1) C_u : y^2 = f_u(x) = c2 Fh(k x^2 - 1, c2 x^2), coefficients of degrees %o in u;\n", [ Degree(Numerator(c)) : c in Coefficients(f) ];
print "    phi2 = (1/x^2, y/x^3) maps C_u onto Y^2 = c2^4 f((k - U)/c2) = E, phi1 = (x^2, y) onto Ea, isomorphic to F_u over Q(u)";

// ---------------------------------------------------------------- (2) the marked group over Q(u), and the gain on D
J := Jacobian(HyperellipticCurve(f)); PJ := Parent((J!0)[1]); xj := PJ.1;
Ub := k - c2*zP; Yb := c2^2*wP;
D6 := elt< J | xj^2 - 1/Ub, (Yb/Ub)*xj >;
ok, iso := IsIsomorphic(Fuu, Ea); assert ok;
Q8 := iso(P8);
D8 := elt< J | xj^2 - Q8[1]/la, PJ!(Q8[2]/la) >;
assert OrderUpTo(D6, 30) eq 6 and OrderUpTo(D8, 30) eq 8 and OrderUpTo(D6 + D8, 60) eq 24;
Us := [ 1/(k - c2*a) : a in al ];
assert &and[ Evaluate(Ha, U0) eq 0 : U0 in Us ];
T2s := [ elt< J | xj^2 - U0, 0 > : U0 in Us ];
assert &and[ 2*T eq J!0 and T ne J!0 : T in T2s ];
assert exists{ T : T in T2s | T ne 12*(D6 + D8) };
print "(2) over Q(u): D6 = phi2^*(P6) has order 6, D8 = phi1^*(P8) order 8, D6 + D8 order 24, and some [x^2 - U_i, 0] is not";
print "    12 (D6 + D8): <D6 + D8, [x^2 - U_i, 0]> = [2,24] is marked over Q(u)";
PU<V> := PolynomialRing(Q);
dpol := -41*(8*u^2 - 1)*(8*u^2 + 8*u + 1);
ok2, r2 := IsSquare(Us[2]); assert ok2;
ok3, r3 := IsSquare(Us[3]*dpol); assert ok3;
printf "(2) the gain: U_2 = (%o)^2 and U_3 dpol(u) = (%o)^2 in Q(u), dpol = %o\n", r2, r3, dpol;
print "    so on D : s^2 = dpol(u), sqrt U_2 = r2 and sqrt U_3 = r3/s are rational, and T = [(x - r2)(x - r3/s), 0] is a rational";
print "    2-torsion class outside the [2,24]: [2,2,24] embeds in J(Q) at every good point of D";

// ---------------------------------------------------------------- the base D and its Jacobian EM
dq := PU!Evaluate(Numerator(dpol), V);
H := HyperellipticCurve(dq);
B := H![-5, 164, 6];                                               // u = -5/6
EH, mH := EllipticCurve(H, B);
EM, mM := MinimalModel(EH);
assert aInvariants(EM) eq [0, 0, 0, -1681, 0];
iso2 := mH*mM; inv := Inverse(iso2);
dp := DefiningPolynomials(inv);
Pinf := EM![-9, -120];
assert Order(Pinf) eq 0;
R0e := iso2(H![1, 41, 3]);                                          // over u = 1/3
im0 := inv(R0e); assert im0[1]/im0[3] eq 1/3;
printf "(4) D : s^2 = %o; EM = %o (rank bounds %o), B over u = -5/6 -> O, P = %o, the member u = 1/3 lies over R0 = %o\n",
    dq, aInvariants(EM), [RankBounds(EM)], Pinf, R0e;

// ---------------------------------------------------------------- (3) exactness on the residue class u = 1/3 mod 13, 47
primes := [13, 47]; ures := Q!1/3;
invs := [];
for p in primes do
    r0 := Integers()!(GF(p)!ures);
    fbar := ReductionClass(f, p, r0);
    n, chi := FrobData(fbar);
    A := AbelianGroup(Jacobian(HyperellipticCurve(fbar)));
    assert #A eq n;
    Append(~invs, Invariants(A));
    printf "(3) p = %o, u = %o mod %o: reduction %o (squarefree, degree %o); #J(F_p) = %o, J(F_p) = %o\n", p, r0, p, fbar, Degree(fbar), n, Invariants(A);
end for;
bnd := StructureBound(invs);
printf "(3) the largest group embedding in all J(F_p) has invariants %o\n", bnd;
assert bnd eq [2, 2, 24];

// ---------------------------------------------------------------- (4) members R_k = R0 + m k P
ords := [ Order(ChangeRing(EM, GF(p))!Eltseq(Pinf)) : p in primes ];
mstep := LCM(ords);
cden := LCM([ Denominator(c) : c in Coefficients(dp[1]) cat Coefficients(dp[3]) ]);
NX := cden*dp[1]; NZ := cden*dp[3];
vz := Evaluate(NZ, [R0e[1], R0e[2], 1]);
for p in primes do
    assert Valuation(Integers()!Discriminant(EM), p) eq 0;
    assert Valuation(R0e[1], p) ge 0 and Valuation(R0e[2], p) ge 0 and Valuation(vz, p) eq 0;
    assert &and[ Valuation(c, p) ge 0 : c in Coefficients(NX) cat Coefficients(NZ) ];
end for;
printf "(4) P has orders %o mod %o; R_k = R0 + %o k P reduces to R0, NZ(R0) is a unit at these primes, so u(R_k) = 1/3 mod p\n", ords, primes, mstep;
R1 := R0e + mstep*Pinf;
u1 := Evaluate(NX, [R1[1], R1[2], 1])/Evaluate(NZ, [R1[1], R1[2], 1]);
assert IsSquare(Evaluate(dq, u1)) and &and[ Valuation(u1 - ures, p) ge 1 : p in primes ];
printf "(4) sanity: u(R_1) has %o digits, lies on D and is = 1/3 mod %o\n", #Sprint(u1), primes;

// ---------------------------------------------------------------- (2), (4), (5) at members
function Member(u0)
    fu := SpecPoly(f, u0);
    C := HyperellipticCurve(fu); Ju := Jacobian(C); PP := Parent((Ju!0)[1]); xu := PP.1;
    sp := func< c | Evaluate(c, u0) >;
    d6 := elt< Ju | xu^2 - 1/sp(Ub), PP!(sp(Yb)/sp(Ub))*xu >;
    Fu0 := Kubert(sp((16*u^3 + 16*u^2 + 6*u + 1)/(8*u^2 - 1)^2), sp((16*u^3 + 16*u^2 + 6*u + 1)/(2*u*(4*u + 1)*(8*u^2 - 1))));
    Ha0 := RX![ Coefficient(fu, 2*i) : i in [0..3] ];
    Ea0, la0 := MonEC(Ha0);
    ok, iso0 := IsIsomorphic(Fu0, Ea0); assert ok;
    q8 := iso0(Fu0![0, 0]);
    d8 := elt< Ju | xu^2 - q8[1]/la0, PP!(q8[2]/la0) >;
    oks, s0 := IsSquare(Evaluate(dq, u0)); assert oks;
    a2 := sp(r2); a3 := sp(r3)/s0;
    assert Evaluate(fu, a2) eq 0 and Evaluate(fu, a3) eq 0;
    tt := elt< Ju | (xu - a2)*(xu - a3), 0 >;
    t2 := [ elt< Ju | xu^2 - sp(U0), 0 > : U0 in Us ];
    G := { i*(d6 + d8) + e*tt + (j eq 0 select Ju!0 else t2[j]) : i in [0..23], e in [0, 1], j in [0..3] };
    return fu, #G, OrderUpTo(d6 + d8, 60);
end function;
G2s := [];
for u0 in [Q| 1/3, 5/27] do
    fu, nG, o24 := Member(u0);
    Append(~G2s, G2Invariants(HyperellipticCurve(fu)));
    invm := TorsionInvariants(fu);
    printf "(2),(5) u = %o: D6 + D8 has order %o, with T and the [x^2 - U_i, 0] it generates %o classes; TorsionSubgroup %o\n", u0, o24, nG, invm;
    assert o24 eq 24 and nG eq 96 and invm eq [2, 2, 24];
end for;
assert G2s[1] ne G2s[2];
C2224a := HyperellipticCurve(RX![0, -158457693276225, -25400003349776, 19670508870560, -10404361990400, 4237596000000]);
assert G2s[1] eq G2Invariants(C2224a);
print "(5) the member u = 1/3 has the G2-invariants of the census curve 2.2.24.a";
print "(4) the G2-invariants at u = 1/3 and u = 5/27 (points of D) differ: u -> M_2 is nonconstant on D, so the R_k give infinitely many isomorphism classes";
printf "time %o s\n", Cputime(T0);
print "FAMILY_2_2_24_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E x F_u) and J(Q)_tors exactly [2,2,24]";
quit;
