/* family_2_16_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E1 x E2) with J(Q)_tors exactly [2,16].
   Found 2026-10-19 (torsion_jac notes/pow2_2026_10_19.md).

   Source.  The Howe-Leprevost-Poonen (2,2)-gluing (HLP 2000, Prop. 4) of E1 = 210e2 (torsion Z/2 x Z/8) with
   curves E2 with full rational 2-torsion, with gain g = 2 coming from the point of order 8 (torsion_jac notes
   gain2_loci_2026_10_15.md, section 2: [2,2] x [2,8] gives only census-known groups; here [2,16]).
   E1 : y^2 + xy = x^3 - 1070x + 7812, monic model Y^2 = (X + 36)(X - 31/4)(X - 28) (Y = y + x/2),
   P8 = (244, 3658) of order 8, 4 P8 = (28, -14).  The Kummer class of P8 is
       delta(P8) = (244 + 36, 244 - 31/4, 244 - 28) = (280, 945/4, 216) = (70, 105, 6) mod squares.
   A curve E2 : Y^2 = (X - e1)(X - e2)(X - e3) glued along alpha : (-36, 31/4, 28) -> (e1, e2, e3) has the gain
   (P8, (e3, 0)) iff the Kummer class of (e3, 0), (e3 - e1, e3 - e2, (e3 - e1)(e3 - e2)), equals delta(P8),
   i.e. e3 - e1 = 70 a^2 and e3 - e2 = 105 b^2.  With e3 = 0, s = a/b (and scaling by b^2, a quadratic twist
   by a square) this is
       E2(s) : Y^2 = X (X + 70 s^2)(X + 105),
   and the gluing is (derived here; HLP's formula, checked below)
       C_s : y^2 = f_s(x) = (128 - 27s^2)((8s^2 - 12)x^2 - 5)(35s^2 x^2 - 32)(140x^2 - 27).
   phi1 : (x, y) -> (x^2, y) maps C_s onto Ea : Y^2 = F(X) (f_s = F(x^2)), Q(s)-isomorphic to E2(s), and
   phi2 : (x, y) -> (1/x^2, y/x^3) onto Eb : Y^2 = X^3 F(1/X), Q(s)-isomorphic to E1 (checked below).

   Argument.
   (1) Splitting.  Wherever f_s is squarefree of degree 6, phi1^* + phi2^* : Ea x Eb -> J is an isogeny over Q
       with kernel of order 4, so J is isogenous over Q to E1 x E2(s).
   (2) The marked [2,16] over Q(s).  Over F = Q(sqrt 70, sqrt 105) the four halves of P8 on E1 are defined
       (DivisionPoints), and over F(s) the halves of T = (0,0) on E2(s) are (X, Y) = (r r', r r'(r + r')),
       r = +-sqrt(70) s, r' = +-sqrt(105) (r^2 = 0 - e1, r'^2 = 0 - e2).  For a half R1 of P8 and a half R2 of T
       the class D16 = pull(R1) + pull(R2) is computed by Mumford arithmetic over F(s); for suitable halves its
       coefficients lie in Q(s), and over Q(s): 16 D16 = 0, 8 D16 != 0.  With Z = the pullback of the 2-torsion
       point (-36, 18) of E1, the 32 sums i D16 + j Z are distinct, so <D16, Z> = Z/16 x Z/2 is marked over Q(s)
       and embeds in J_s(Q)_tors at every smooth member.
   (3) Exactness on a residue class.  #J(F_p) = #E1(F_p) #E2(F_p) is divisible by 16 * 4 = 64 at every good p,
       so the group structure is used.  For s = 1 mod 11 and s = 2 mod 13 (s 11- and 13-integral) the curve
       has fixed reductions with J(F_11) = [8,16] (2-part partition (4,3)) and J(F_13) = [2,2,80] ((4,1,1), odd
       part 5).  The termwise minimum of the 2-part partitions is (4,1), the partition of [2,16], and the odd
       parts 1 and 5 are coprime; J_s(Q)_tors embeds in both and contains <D16, Z>, so it is exactly [2,16].
   (4) Infinitely many.  The G2-invariants of the members differ, so s -> M_2 is nonconstant with finite fibres,
       and the residue class gives infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup at two members is [2,16].
   Run (from the repository root):
       magma -b pipeline/magma/family_2_16_split_certificate.m > family_2_16_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

function MonEC(G)       // Y^2 = c3 X^3 + c2 X^2 + c1 X + c0  ->  monic model (x, y) = (c3 X, c3 Y)
    c := [ Coefficient(G, i) : i in [0..3] ];
    return EllipticCurve([0, c[3], 0, c[2]*c[4], c[1]*c[4]^2]), c[4];
end function;
function Ord2(g)      // order of a class killed by 64, else -1
    z := Parent(g)!0; h := g; k := 1;
    while k le 64 do if h eq z then return k; end if; h := 2*h; k := 2*k; end while;
    return -1;
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

// ---------------------------------------------------------------- (0) E1, E2(s), the HLP gluing
E1 := EllipticCurve([Q| 1, 0, 0, -1070, 7812]);
assert Conductor(E1) eq 210 and Invariants(TorsionSubgroup(E1)) eq [2, 8];
P8 := E1![244, 3658];
assert Order(P8) eq 8 and 4*P8 eq E1![28, -14];
PQ := RX; z := X;
assert z^3 + z^2/4 - 1070*z + 7812 eq (z + 36)*(z - 31/4)*(z - 28);
dl := [ 244 - r : r in [-36, 31/4, 28] ];
assert dl eq [280, 945/4, 216];
print "E1 = 210e2, monic model Y^2 = (X + 36)(X - 31/4)(X - 28), P8 = (244, 3658), delta(P8) = (280, 945/4, 216) = (70, 105, 6) mod squares";
Fs<s> := FunctionField(Q); Px<x> := PolynomialRing(Fs);
al := [Fs| -36, 31/4, 28];
be := [Fs| -70*s^2, -105, 0];
assert IsSquare((be[3] - be[1])/(70*s^2)) and IsSquare(dl[1]/70) and IsSquare(dl[2]/105) and IsSquare((be[3] - be[2])/105);
E2 := EllipticCurve([Fs| 0, 70*s^2 + 105, 0, 7350*s^2, 0]);
print "E2(s) : Y^2 = X(X + 70s^2)(X + 105); alpha : (-36, 31/4, 28) -> (-70s^2, -105, 0)";
aa1 := (al[3]-al[2])^2/(be[3]-be[2]) + (al[2]-al[1])^2/(be[2]-be[1]) + (al[1]-al[3])^2/(be[1]-be[3]);
bb1 := (be[3]-be[2])^2/(al[3]-al[2]) + (be[2]-be[1])^2/(al[2]-al[1]) + (be[1]-be[3])^2/(al[1]-al[3]);
aa2 := al[1]*(be[3]-be[2]) + al[2]*(be[1]-be[3]) + al[3]*(be[2]-be[1]);
bb2 := be[1]*(al[3]-al[2]) + be[2]*(al[1]-al[3]) + be[3]*(al[2]-al[1]);
Df := ((al[1]-al[2])*(al[1]-al[3])*(al[2]-al[3]))^2;
Dg := ((be[1]-be[2])*(be[1]-be[3])*(be[2]-be[3]))^2;
AA := Dg*aa1/aa2; BB := Df*bb1/bb2;
hHLP := -(AA*(al[2]-al[1])*(al[1]-al[3])*x^2 + BB*(be[2]-be[1])*(be[1]-be[3]))
        *(AA*(al[3]-al[2])*(al[2]-al[1])*x^2 + BB*(be[3]-be[2])*(be[2]-be[1]))
        *(AA*(al[1]-al[3])*(al[3]-al[2])*x^2 + BB*(be[1]-be[3])*(be[3]-be[2]));
f := (128 - 27*s^2)*((8*s^2 - 12)*x^2 - 5)*(35*s^2*x^2 - 32)*(140*x^2 - 27);
oksq, rt := IsSquare(Fs!(hHLP/f));
assert oksq;
print "C_s : y^2 = f_s(x) = (128 - 27s^2)((8s^2 - 12)x^2 - 5)(35s^2 x^2 - 32)(140x^2 - 27)";
printf "the HLP (Prop. 4) curve is (%o)^2 f_s\n", rt;

// ---------------------------------------------------------------- (1) the two elliptic quotients
F3 := Px![ Coefficient(f, 2*k) : k in [0..3] ];
Fr := Px![ Coefficient(F3, 3 - k) : k in [0..3] ];
assert Evaluate(F3, x^2) eq f;
Ea, ca := MonEC(F3); Eb, cb := MonEC(Fr);
E1s := ChangeRing(E1, Fs);
okA, mA := IsIsomorphic(E2, Ea); okB, mB := IsIsomorphic(E1s, Eb);
assert okA and okB;
print "(1) phi1 = (x^2, y): C_s -> Ea ~ E2(s) and phi2 = (1/x^2, y/x^3): C_s -> Eb ~ E1 over Q(s): J is isogenous over Q to E1 x E2(s)";

// ---------------------------------------------------------------- (2) the marked [2,16] over Q(s)
J := Jacobian(HyperellipticCurve(f));
function Pull1(JJ, Pt, c) PP := Parent((JJ!0)[1]); xx := PP.1; XX := Pt[1]/c; YY := Pt[2]/c; return elt< JJ | xx^2 - XX, PP!YY >; end function;
function Pull2(JJ, Pt, c) PP := Parent((JJ!0)[1]); xx := PP.1; XX := Pt[1]/c; YY := Pt[2]/c; return elt< JJ | xx^2 - 1/XX, (YY/XX)*xx >; end function;
B8 := Pull2(J, mB(E1s!Eltseq(P8)), cb);
Zt := Pull2(J, mB(E1s![-36, 18]), cb);
assert Ord2(B8) eq 8 and Ord2(Zt) eq 2;
dA := IsomorphismData(mA); dB := IsomorphismData(mB);
im := mB(E1s!Eltseq(P8)); r, ss, tt, uu := Explode(dB);
conv := (uu^2*P8[1] + r eq im[1] and uu^3*P8[2] + ss*uu^2*P8[1] + tt eq im[2]) select 2 else 1;
if conv eq 1 then assert (P8[1]-r)/uu^2 eq im[1] and (P8[2] - ss*(P8[1]-r) - tt)/uu^3 eq im[2]; end if;
function ApplyIso(dat, pt, sub)
    r, ss, tt, uu := Explode([sub(c) : c in dat]);
    if conv eq 1 then return [(pt[1]-r)/uu^2, (pt[2] - ss*(pt[1]-r) - tt)/uu^3]; end if;
    return [uu^2*pt[1] + r, uu^3*pt[2] + ss*uu^2*pt[1] + tt];
end function;

F := AbsoluteField(NumberField([X^2 - 70, X^2 - 105]));
r70 := Sqrt(F!70); r105 := Sqrt(F!105);
FF<w> := FunctionField(F); PF<xF> := PolynomialRing(FF);
sub := func< c | Evaluate(Numerator(c), w)/Evaluate(Denominator(c), w) >;
JF := Jacobian(HyperellipticCurve(PF![ sub(c) : c in Coefficients(f) ]));
caF := sub(ca); cbF := sub(cb);
function PullF1(pt) im := ApplyIso(dA, pt, sub); return elt< JF | xF^2 - im[1]/caF, PF!(im[2]/caF) >; end function;
function PullF2(pt) im := ApplyIso(dB, pt, sub); XX := im[1]/cbF; YY := im[2]/cbF; return elt< JF | xF^2 - 1/XX, (YY/XX)*xF >; end function;
halvesP := DivisionPoints(ChangeRing(E1, F)!Eltseq(P8), 2);
assert #halvesP eq 4;
E2F := EllipticCurve([FF| 0, 70*w^2 + 105, 0, 7350*w^2, 0]);
halvesT := [];
for sb in [1, -1] do for sc in [1, -1] do
    rr := sb*r70*w; rc := sc*r105;
    assert rr^2 eq 0 - (-70*w^2) and rc^2 eq 0 - (-105);
    R2 := E2F![rr*rc, rr*rc*(rr + rc)];
    assert 2*R2 eq E2F![0, 0];
    Append(~halvesT, R2);
end for; end for;
PQs := PolynomialRing(Q);
function InQ(c)
    nu := Numerator(c); de := Denominator(c); lc := LeadingCoefficient(de);
    return &and[ IsCoercible(Q, a) : a in Coefficients(nu/lc) cat Coefficients(de/lc) ];
end function;
function ToQs(c)
    nu := Numerator(c); de := Denominator(c); lc := LeadingCoefficient(de);
    return Evaluate(PQs![Q!a : a in Coefficients(nu/lc)], s)/Evaluate(PQs![Q!a : a in Coefficients(de/lc)], s);
end function;
D16F := 0;
for R1 in halvesP do
    for R2 in halvesT do
        D := PullF2(R1) + PullF1(R2);
        if &and[ InQ(c) : c in Coefficients(D[1]) cat Coefficients(D[2]) ] then D16F := D; break R1; end if;
    end for;
end for;
assert D16F cmpne 0;
D16 := elt< J | Px![ ToQs(c) : c in Coefficients(D16F[1]) ], Px![ ToQs(c) : c in Coefficients(D16F[2]) ], D16F[3] >;
BT := Pull1(J, mA(E2![0, 0]), ca);
assert Ord2(D16) eq 16 and 2*D16 eq B8 + BT;
GG := { i*D16 + j*Zt : i in [0..15], j in [0..1] };
assert #GG eq 32;
printf "(2) over F(s), F = Q(sqrt 70, sqrt 105): D16 = pull(R1) + pull(R2), 2R1 = P8, 2R2 = (0,0) on E2(s), has coefficients in Q(s):\n     D16 = %o\n", D16;
print "    over Q(s): D16 has order 16, 2 D16 = pull(P8) + pull((0,0)), and with Z = pull((-36,18)) the 32 sums i D16 + j Z are distinct:";
print "    <D16, Z> = Z/16 x Z/2 is marked over Q(s)";

// ---------------------------------------------------------------- (3) exactness on a residue class
primes := [11, 13]; residues := [1, 2];
sres := CRT(residues, primes); M := &*primes;
members := [Q| sres, sres - M, sres + M, (sres + M)/(1 + M)];
invs := [];
for k in [1..#primes] do
    p := primes[k]; r0 := residues[k];
    fbar := ReductionClass(f, p, r0);
    for m in members do
        assert Valuation(m, p) ge 0 and GF(p)!m eq GF(p)!r0;
        assert PolynomialRing(GF(p))!SpecPoly(f, m) eq fbar;
    end for;
    n, chi := FrobData(fbar);
    Ap := AbelianGroup(Jacobian(HyperellipticCurve(fbar)));
    assert #Ap eq n;
    Append(~invs, Invariants(Ap));
    printf "(3) p = %o, s = %o mod %o: reduction %o (squarefree, degree %o), J(F_%o) = %o (order %o), chi = %o\n",
        p, r0, p, fbar, Degree(fbar), p, Invariants(Ap), n, chi;
end for;
assert invs eq [[8, 16], [2, 2, 80]];
bnd := StructureBound(invs);
printf "(3) the largest group embedding in all J(F_p) has invariants %o\n", bnd;
assert bnd eq [2, 16];
printf "(3) every s = %o mod %o, 11- and 13-integral, has J(Q)_tors = <D16, Z> = [2,16]\n", sres, M;

// ---------------------------------------------------------------- (4) infinitely many, (5) sanity
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } eq #members;
printf "(4) the G2-invariants of the members %o are pairwise distinct\n", members;
for m in members[1..2] do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(5) TorsionSubgroup at s = %o: %o\n", m, invm;
    assert invm eq [2, 16];
end for;
printf "time %o s\n", Cputime(T0);
print "FAMILY_2_16_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E1 x E2) and J(Q)_tors exactly [2,16]";
quit;
