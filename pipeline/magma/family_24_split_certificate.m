/* family_24_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E x F_m) with J(Q)_tors exactly [24].
   Found 2026-10-19 (torsion_jac notes/hlpEC_2026_10_19.md).

   Source.  The Howe-Leprevost-Poonen (2,2)-gluing (HLP 2000, Prop. 4 and section 3.5, the case N = 8, N' = 6,
   which HLP leave as "subsumed") of the fixed Kubert X_1(8) curve
       E = E_8(3/2) = 21a3 : y^2 + (1 - c)xy - by = x^3 - bx^2, b = 1, c = 2/3  (torsion Z/8, P_E = (0,0)),
   with the Kubert X_1(6) curves F_s = E_6(s) (b = s^2 + s, c = s, (0,0) of order 6).  Both have exactly one
   rational 2-torsion point, so E[2] = F_s[2] as Galois modules exactly when the quadratic fields of their other
   2-torsion points agree: (s + 1)(9s + 1) = 7 w^2 (Delta_8(3/2) = 7 mod squares).  This conic has the point
   (s, w) = (-1, 0), and w = m(s + 1) gives s = (1 - 7m^2)/(7m^2 - 9).
   Model.  Write E as y^2 = f(z) = z^3 + pE z + qE (z = x + sh).  The isomorphism psi : E[2] -> F_s[2] is given by
   the Tschirnhausen transformation phi(z) = c1 z + c2 (z^2 + 2pE/3), with c1, c2 in Q(m) below: the characteristic
   polynomial of phi(a), a a root of f, is the depressed 2-division cubic z^3 + a1(s) z + a0(s) of F_s (checked over
   Q(m)).  (The other choice of psi is m -> -m.)  The HLP gluing along psi is, up to a square,
       C_m : y^2 = h(x) = c2 Fh(c1 x^2 - 1, c2 x^2),   Fh(X, Z) = X^3 + pE X Z^2 + qE Z^3
   (closed form of HLP Prop. 4, see the note).  phi2 : (x, y) -> (1/x^2, y/x^3) maps C_m onto Eb : Y^2 = c2^4 f((c1 - U)/c2),
   which is E, and phi1 : (x, y) -> (x^2, y) maps C_m onto Ea : Y^2 = Ha(U) (h = Ha(x^2)), checked to be F_s over Q(m).

   Argument.
   (1) Splitting.  phi1^* + phi2^* : Ea x Eb -> J is an isogeny over Q wherever h is squarefree of degree 6, so J is
       isogenous over Q to E x F_s.
   (2) The marked group over Q(m).  phi1, phi2 have prime degree 2 and C is not of genus 1, so they do not factor
       through isogenies and phi1^*, phi2^* are injective.  D8 = phi2^*(P_E) has order 8 and D6 = phi1^*(Q6), Q6 the
       image of (0,0) of F_s, has order 6.  Mumford arithmetic over Q(m): 4 D8 = 3 D6 (the rational 2-torsion points
       are identified by the gluing) and D8 + D6 has order 24.  So Z/24 is marked over Q(m).
   (3) Exactness on a residue class.  #J(F_p) = #E(F_p) #F_s(F_p) is divisible by 48, so the group structure is used.
       For every 5- and 17-integral m = 56 (mod 85), C_m has the reductions of C_56 at 5 and 17, with
       J(F_5) = Z/48 (2-primary part cyclic, 3-primary part Z/3) and J(F_17) = Z/24 x Z/24 (2-primary exponent 8).
       J_m(Q)_tors embeds in both, and the largest group that does is Z/24: J_m(Q)_tors is exactly Z/24.
   (4) Infinitely many.  The members m = 56, 141, -29 have pairwise distinct G2-invariants, so m -> M_2 is
       nonconstant with finite fibres, and the residue class gives infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup at m = 56, 141, -29 is [24].  The member m = -1 (outside the class) is the census curve
       48.a (Howe 2015; Platonov-Petrunin), with TorsionSubgroup [48]: there the gain of HLP Prop. 12 halves the point of
       order 24.  The census curve 24.b is not a member (its Jacobian is isogenous to curves of conductor 15 and 20, and
       every member has the factor 21a3).
   Run (from the repository root):
       magma -b pipeline/magma/family_24_split_certificate.m > data/knowledge/sources/family_24_split_certificate.log */
SetColumns(0);
SetQuitOnError(true);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

function Kubert(b, c) return EllipticCurve([1 - c, -b, -b, 0, 0]); end function;
function DepCubic(E)               // Y^2 = z^3 + A z + B, z = x + sh, Y = y + (a1 x + a3)/2
    b2, b4, b6 := Explode(bInvariants(E)[1..3]);
    A2 := b2/4; A1 := b4/2; A0 := b6/4;
    return A1 - A2^2/3, A0 - A1*A2/3 + 2*A2^3/27, A2/3;
end function;
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

// ---------------------------------------------------------------- (0) the fixed curve E = E_8(3/2) = 21a3
E := Kubert(1, 2/3);
assert Conductor(E) eq 21 and Invariants(TorsionSubgroup(E)) eq [8];
PE := E![0, 0]; assert Order(PE) eq 8;
pE, qE, shE := DepCubic(E);
assert pE eq -1873/3888 and qE eq 80569/629856;
aE := aInvariants(E); zP := shE; wP := aE[3]/2;
fE := X^3 + pE*X + qE;
assert Evaluate(fE, zP) eq wP^2;
printf "(0) E = E_8(3/2) = %o = 21a3, torsion Z/8; f(z) = %o = %o\n", aE, fE, Factorization(fE);

// ---------------------------------------------------------------- (1) the family over Q(m)
Fm<m> := FunctionField(Q);
sm := (1 - 7*m^2)/(7*m^2 - 9);
assert (sm + 1)*(9*sm + 1) eq 7*(m*(sm + 1))^2;
c1 := (-184/27*m^4 + 368/63*m^2 + 3260/1323*m + 184/441)/(m^4 - 18/7*m^2 + 81/49);
c2 := (-16*m^4 + 96/7*m^2 + 368/49*m + 48/49)/(m^4 - 18/7*m^2 + 81/49);
Fs := Kubert(sm^2 + sm, sm);
a1s, a0s, shs := DepCubic(Fs);
R2<u1, u2> := PolynomialRing(Q, 2);
Ma := Matrix(Q, 3, 3, [0, 1, 0, 0, 0, 1, -qE, -pE, 0]);
MM := u1*ChangeRing(Ma, R2) + u2*(ChangeRing(Ma*Ma, R2) + (2*pE/3)*IdentityMatrix(R2, 3));
cp := CharacteristicPolynomial(MM);
assert Coefficient(cp, 2) eq 0;
assert Evaluate(Coefficient(cp, 1), [c1, c2]) eq a1s and Evaluate(Coefficient(cp, 0), [c1, c2]) eq a0s;
print "(1) s = (1 - 7m^2)/(7m^2 - 9) satisfies (s + 1)(9s + 1) = 7 (m(s + 1))^2; over Q(m) the characteristic polynomial of";
print "    phi(a) = c1 a + c2 (a^2 + 2pE/3) is the depressed 2-division cubic of F_s = E_6(s): F_s[2] = E[2] via (a_i,0) -> (phi(a_i),0)";
P<x> := PolynomialRing(Fm);
Fh := func< XX, ZZ | XX^3 + pE*XX*ZZ^2 + qE*ZZ^3 >;
f := c2*Fh(c1*x^2 - 1, c2*x^2);
assert Degree(f) eq 6 and IsSquarefree(f) and &and[ Coefficient(f, i) eq 0 : i in [1, 3, 5] ];
Ha := P![ Coefficient(f, 2*i) : i in [0..3] ];
Hb := P![ Coefficient(f, 6 - 2*i) : i in [0..3] ];
assert Hb eq c2^4*(((c1 - x)/c2)^3 + pE*(c1 - x)/c2 + qE);
Ea, la := MonEC(Ha);
c4a, c6a := Explode(cInvariants(Ea)); c4f, c6f := Explode(cInvariants(Fs));
lam := (c6a/c6f)/(c4a/c4f);
assert c4a eq lam^2*c4f and c6a eq lam^3*c6f and IsSquare(lam);
printf "(1) C_m : y^2 = f_m(x) = c2 Fh(c1 x^2 - 1, c2 x^2), coefficients of degrees %o in m;\n", [ Degree(Numerator(c)) : c in Coefficients(f) ];
print "    phi2 = (1/x^2, y/x^3) maps C_m onto Y^2 = c2^4 f((c1 - U)/c2) = E, phi1 = (x^2, y) onto Ea, isomorphic to F_s over Q(m)";

// ---------------------------------------------------------------- (2) the marked group over Q(m)
J := Jacobian(HyperellipticCurve(f));
PJ := Parent((J!0)[1]); xj := PJ.1;
Ub := c1 - c2*zP; Yb := c2^2*wP;                                   // P_E on Eb
D8 := elt< J | xj^2 - 1/Ub, (Yb/Ub)*xj >;
ok, iso := IsIsomorphic(Fs, Ea); assert ok;
Q6 := iso(Fs![0, 0]);
D6 := elt< J | xj^2 - Q6[1]/la, PJ!(Q6[2]/la) >;
assert OrderUpTo(D8, 30) eq 8 and OrderUpTo(D6, 30) eq 6;
assert 4*D8 eq 3*D6 and OrderUpTo(D8 + D6, 50) eq 24;
print "(2) over Q(m): D8 = phi2^*(P_E) has order 8, D6 = phi1^*((0,0) of F_s) has order 6, 4 D8 = 3 D6, and D8 + D6 has order 24:";
print "    Z/24 is marked over Q(m)";

// ---------------------------------------------------------------- (3) exactness on the residue class m = 56 mod 85
primes := [5, 17]; mres := 56;
members := [Q| 56, 141, -29];
invs := [];
for p in primes do
    r0 := Integers()!(GF(p)!mres);
    fbar := ReductionClass(f, p, r0);
    for mv in members do
        assert GF(p)!mv eq GF(p)!r0;
        assert PolynomialRing(GF(p))!SpecPoly(f, mv) eq fbar;
    end for;
    n, chi := FrobData(fbar);
    A := AbelianGroup(Jacobian(HyperellipticCurve(fbar)));
    assert #A eq n;
    Append(~invs, Invariants(A));
    printf "(3) p = %o, m = %o mod %o: reduction %o (squarefree, degree %o); #J(F_p) = %o, J(F_p) = %o\n",
        p, r0, p, fbar, Degree(fbar), n, Invariants(A);
end for;
bnd := StructureBound(invs);
printf "(3) the largest group embedding in all J(F_p) has invariants %o\n", bnd;
assert bnd eq [24];
printf "(3) every 5- and 17-integral m = %o mod %o has J(Q)_tors exactly Z/24\n", mres, &*primes;

// ---------------------------------------------------------------- (4), (5)
G2s := [ G2Invariants(HyperellipticCurve(SpecPoly(f, mv))) : mv in members ];
assert #Seqset(G2s) eq #members;
printf "(4) the G2-invariants of the members %o are pairwise distinct\n", members;
for mv in members do
    invm := TorsionInvariants(SpecPoly(f, mv));
    printf "(5) TorsionSubgroup at m = %o (s = %o): %o\n", mv, Evaluate(sm, mv), invm;
    assert invm eq [24];
end for;
C48a := HyperellipticCurve(RX![3, 0, 9, 0, -4], RX![0, 1, 0, 1]);
assert G2Invariants(HyperellipticCurve(SpecPoly(f, -1))) eq G2Invariants(C48a) and TorsionInvariants(SpecPoly(f, -1)) eq [48];
print "(5) remark: the member m = -1 (not in the class) is the census curve 48.a, with TorsionSubgroup [48] (the gain of HLP Prop. 12)";
printf "time %o s\n", Cputime(T0);
print "FAMILY_24_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E x F_m) and J(Q)_tors exactly [24]";
quit;
