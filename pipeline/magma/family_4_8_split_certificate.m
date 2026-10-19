/* family_4_8_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E1 x E2) with J(Q)_tors exactly [4,8].
   Found 2026-10-19 (torsion_jac notes/pow2_2026_10_19.md).

   Source.  The Howe-Leprevost-Poonen (2,2)-gluing (HLP 2000, Prop. 4) of
       E1 = E_{2,4}(s) : y^2 + xy - b y = x^3 - b x^2,  b = s^2 - 1/16,  P1 = (0,0) of order 4,  s = u/4,
       E2 = 210e2 : y^2 + xy = x^3 - 1070 x + 7812  (torsion Z/2 x Z/8),  P2 = (244, 3658) of order 8,
   along the isomorphism alpha : E1[2] -> E2[2] that sends the special 2-torsion point 2 P1 of E1 to a
   NON-special 2-torsion point of E2 (not 4 P2).  In the monic models Y^2 = (X - e1)(X - e2)(X - e3)
   (Y = y + (x - b)/2, resp. Y = y + x/2) the roots are
       E1: (s^2 - 1/16, -s/2 - 1/8, s/2 - 1/8)   (special root x(2 P1) first),
       E2: (-36, 31/4, 28)                       (x(4 P2) = 28 last),
   and alpha matches them in this order: the "sp -> ns" gluing of [2,4] x [2,8] with gain g = 1 (torsion_jac notes
   gain2_loci_2026_10_15.md, section 2: [2,4] x [2,8], g = 1, sp -> ns gives [4,8]).  No condition on s is needed:
   the gluing exists over Q(s) because both curves have full rational 2-torsion.
   Up to a square, the HLP curve is (derived here, checked below)
       C_u : y^2 = f_u(x) = -(81u^2 + 862u + 81)(81x^2 - u)(1024x^2 + (u - 1)^2)(700x^2 + (u + 1)^2),
   f_u = F(x^2) with F a cubic over Q(u), so
       phi1 : (x, y) -> (x^2, y)       maps C_u onto Ea : Y^2 = F(X),
       phi2 : (x, y) -> (1/x^2, y/x^3) maps C_u onto Eb : Y^2 = X^3 F(1/X),
   and Ea, Eb are Q(u)-isomorphic to E1(u/4), E2 in some order (checked below with IsIsomorphic over Q(u)).

   Argument.
   (1) Splitting.  Wherever f_u is squarefree of degree 6, phi1^* + phi2^* : Ea x Eb -> J is an isogeny over Q
       with kernel of order 4, so J is isogenous over Q to E1 x E2.
   (2) The marked group over Q(u).  A = pullback of P1 (order 4) and B = pullback of P2 (order 8); by Mumford
       arithmetic over Q(u) the 32 sums iA + jB (0 <= i < 4, 0 <= j < 8) are distinct, so <A, B> = Z/4 x Z/8.
       (With alpha sp -> sp instead, 2A = 4B and one gets [2,2,8]; this is the difference.)  The marked torsion
       sections specialize injectively on smooth fibres, so [4,8] embeds in J_u(Q)_tors at every smooth member.
   (3) Exactness on a residue class.  #J(F_p) = #E1(F_p) #E2(F_p) is divisible by 8 * 16 = 128 at every good p,
       so a gcd of point counts cannot certify [4,8]; the group structure is used (as in family_18_split and
       family_88_split).  Reduction J_u(Q)_tors -> J(F_p) is injective (p odd, good reduction).  For
       u = 4 mod 11, u = 8 mod 13, u = 10 mod 41 (u p-integral) the curve has a fixed reduction at each prime with
           J(F_11) = [2,4,16] (2-part partition (4,2,1)),  J(F_13) = [8,16] ((4,3)),  J(F_41) = [8,8,24] ((3,3,3)).
       The 2-part H of J_u(Q)_tors embeds in all three, so its partition is bounded termwise by (4,2,1), (4,3) and
       (3,3,3), i.e. by their termwise minimum (3,2), the partition of [4,8]; the odd parts 1, 1, 3 have gcd 1.
       As H contains <A, B> = [4,8], J_u(Q)_tors = [4,8] exactly for every such u.
   (4) Infinitely many.  The G2-invariants of the members differ, so u -> M_2 is nonconstant with finite fibres
       and the residue class gives infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup at two members is [4,8].
   Run (from the repository root):
       magma -b pipeline/magma/family_4_8_split_certificate.m > family_4_8_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

function Sqf(q)
    q := Rationals()!q; n := Numerator(q)*Denominator(q);
    return Sign(n)*&*([1] cat [p[1] : p in Factorization(Abs(n)) | p[2] mod 2 eq 1]);
end function;
function MonEC(G)       // Y^2 = c3 X^3 + c2 X^2 + c1 X + c0  ->  monic model (x, y) = (c3 X, c3 Y)
    c := [ Coefficient(G, i) : i in [0..3] ];
    return EllipticCurve([0, c[3], 0, c[2]*c[4], c[1]*c[4]^2]), c[4];
end function;
function MonicCubic(E, x)
    a1, a2, a3, a4, a6 := Explode(aInvariants(E));
    return x^3 + (a2 + a1^2/4)*x^2 + (a4 + a1*a3/2)*x + (a6 + a3^2/4);
end function;

// ---------------------------------------------------------------- (0) the two factors and the HLP gluing
Fu<u> := FunctionField(Q); Px<x> := PolynomialRing(Fu);
s := u/4; b1 := s^2 - 1/16;
E1 := EllipticCurve([Fu| 1, -b1, -b1, 0, 0]);
P1 := E1![0, 0];
assert 4*P1 eq E1!0 and 2*P1 ne E1!0;
al := [s^2 - 1/16, -s/2 - 1/8, s/2 - 1/8];
assert &and[ Evaluate(MonicCubic(E1, x), a) eq 0 : a in al ] and (2*P1)[1] eq al[1];
E2 := EllipticCurve([Q| 1, 0, 0, -1070, 7812]);
assert Conductor(E2) eq 210 and Invariants(TorsionSubgroup(E2)) eq [2, 8];
P2 := E2![244, 3658];
assert Order(P2) eq 8;
be := [Fu| -36, 31/4, 28];
assert &and[ Evaluate(MonicCubic(ChangeRing(E2, Fu), x), c) eq 0 : c in be ] and (4*P2)[1] eq 28;
print "E1 = E_{2,4}(u/4), P1 = (0,0) of order 4, monic roots (s^2-1/16, -s/2-1/8, s/2-1/8), special first";
print "E2 = 210e2 = [1,0,0,-1070,7812], torsion [2,8], P2 = (244,3658) of order 8, monic roots (-36, 31/4, 28), x(4P2) = 28";
print "alpha: roots matched in this order (special root of E1 -> -36, non-special on E2)";
a1 := (al[3]-al[2])^2/(be[3]-be[2]) + (al[2]-al[1])^2/(be[2]-be[1]) + (al[1]-al[3])^2/(be[1]-be[3]);
bb1 := (be[3]-be[2])^2/(al[3]-al[2]) + (be[2]-be[1])^2/(al[2]-al[1]) + (be[1]-be[3])^2/(al[1]-al[3]);
a2 := al[1]*(be[3]-be[2]) + al[2]*(be[1]-be[3]) + al[3]*(be[2]-be[1]);
bb2 := be[1]*(al[3]-al[2]) + be[2]*(al[1]-al[3]) + be[3]*(al[2]-al[1]);
Df := ((al[1]-al[2])*(al[1]-al[3])*(al[2]-al[3]))^2;
Dg := ((be[1]-be[2])*(be[1]-be[3])*(be[2]-be[3]))^2;
AA := Dg*a1/a2; BB := Df*bb1/bb2;
hHLP := -(AA*(al[2]-al[1])*(al[1]-al[3])*x^2 + BB*(be[2]-be[1])*(be[1]-be[3]))
        *(AA*(al[3]-al[2])*(al[2]-al[1])*x^2 + BB*(be[3]-be[2])*(be[2]-be[1]))
        *(AA*(al[1]-al[3])*(al[3]-al[2])*x^2 + BB*(be[1]-be[3])*(be[3]-be[2]));
f := -(81*u^2 + 862*u + 81)*(81*x^2 - u)*(1024*x^2 + (u - 1)^2)*(700*x^2 + (u + 1)^2);
ratio := Fu!(hHLP/f);
oksq, rt := IsSquare(ratio);
assert oksq;
print "C_u : y^2 = f_u(x) = -(81u^2 + 862u + 81)(81x^2 - u)(1024x^2 + (u - 1)^2)(700x^2 + (u + 1)^2)";
printf "the HLP (Prop. 4) curve is (%o)^2 f_u: C_u is Q(u)-isomorphic to the gluing\n", rt;

// ---------------------------------------------------------------- (1) the two elliptic quotients
F := Px![ Coefficient(f, 2*i) : i in [0..3] ];
Fr := Px![ Coefficient(F, 3 - i) : i in [0..3] ];
assert Evaluate(F, x^2) eq f;
Ea, ca := MonEC(F); Eb, cb := MonEC(Fr);
ok1, m1 := IsIsomorphic(E1, Eb);
E2u := ChangeRing(E2, Fu);
ok2, m2 := IsIsomorphic(E2u, Ea);
assert ok1 and ok2;
print "(1) phi1 = (x^2, y): C_u -> Ea : Y^2 = F(X) ~ E2 and phi2 = (1/x^2, y/x^3): C_u -> Eb : Y^2 = X^3 F(1/X) ~ E1(u/4) over Q(u):";
print "    J is isogenous over Q to E1 x E2";

// ---------------------------------------------------------------- (2) the marked [4,8] over Q(u)
J := Jacobian(HyperellipticCurve(f));
function Pull1(JJ, Pt, c) PP := Parent((JJ!0)[1]); xx := PP.1; XX := Pt[1]/c; YY := Pt[2]/c; return elt< JJ | xx^2 - XX, PP!YY >; end function;
function Pull2(JJ, Pt, c) PP := Parent((JJ!0)[1]); xx := PP.1; XX := Pt[1]/c; YY := Pt[2]/c; return elt< JJ | xx^2 - 1/XX, (YY/XX)*xx >; end function;
A := Pull2(J, m1(P1), cb);
B := Pull1(J, m2(E2u!Eltseq(P2)), ca);
assert 4*A eq J!0 and 2*A ne J!0 and 8*B eq J!0 and 4*B ne J!0;
GG := { i*A + j*B : i in [0..3], j in [0..7] };
assert #GG eq 32;
assert 2*A ne 4*B;
print "(2) over Q(u): A = phi2^*(P1) has order 4, B = phi1^*(P2) has order 8, and the 32 sums iA + jB are distinct:";
print "    <A, B> = Z/4 x Z/8 is marked over Q(u) (2A != 4B: the gluing sends 2 P1 to a non-special point)";
// the gluing kernel: pullbacks of matched 2-torsion points cancel
T1s := [ E1![a, -(a - b1)/2] : a in al ];
T2s := [ E2u![c, -c/2] : c in be ];
for i in [1..3] do assert Pull2(J, m1(T1s[i]), cb) + Pull1(J, m2(T2s[i]), ca) eq J!0; end for;
print "    pullbacks of matched 2-torsion points cancel: the gluing kernel is graph(alpha)";

// ---------------------------------------------------------------- (3) exactness on a residue class
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
primes := [11, 13, 41]; residues := [4, 8, 10];
ures := CRT(residues, primes); M := &*primes;
members := [Q| ures, ures - M, ures + M, (ures + M)/(1 + M)];
invs := [];
for i in [1..#primes] do
    p := primes[i]; r0 := residues[i];
    fbar := ReductionClass(f, p, r0);
    for m in members do
        assert Valuation(m, p) ge 0 and GF(p)!m eq GF(p)!r0;
        assert PolynomialRing(GF(p))!SpecPoly(f, m) eq fbar;
    end for;
    n, chi := FrobData(fbar);
    Ap := AbelianGroup(Jacobian(HyperellipticCurve(fbar)));
    assert #Ap eq n;
    Append(~invs, Invariants(Ap));
    printf "(3) p = %o, u = %o mod %o: reduction %o (squarefree, degree %o), J(F_%o) = %o (order %o), chi = %o\n",
        p, r0, p, fbar, Degree(fbar), p, Invariants(Ap), n, chi;
end for;
assert invs eq [[2,4,16], [8,16], [8,8,24]];
bnd := StructureBound(invs);
printf "(3) the largest group embedding in all J(F_p) has invariants %o\n", bnd;
assert bnd eq [4, 8];
printf "(3) every u = %o mod %o, p-integral for p in %o, has J(Q)_tors = <A, B> = [4,8]\n", ures, M, primes;

// ---------------------------------------------------------------- (4) infinitely many, (5) sanity
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } eq #members;
printf "(4) the G2-invariants of the members %o are pairwise distinct\n", members;
for m in members[1..2] do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(5) TorsionSubgroup at u = %o: %o\n", m, invm;
    assert invm eq [4, 8];
end for;
printf "time %o s\n", Cputime(T0);
print "FAMILY_4_8_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E1 x E2) and J(Q)_tors exactly [4,8]";
quit;
