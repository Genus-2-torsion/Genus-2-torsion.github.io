/* family_2_24_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E1 x E2) with J(Q)_tors exactly [2,24].
   Found 2026-10-19 (torsion_jac notes/hlpP2_2026_10_19.md).

   Source.  The Z/2 x Z/24 entry of Howe-Leprevost-Poonen 2000, Table 1 (section 3.4): the (2,2)-gluing
   (Prop. 4) of a curve with torsion Z/2 x Z/8 and a curve with torsion Z/2 x Z/6 (all 2-torsion rational, so
   every bijection of the 2-torsion is Galois), restricted to the line on which the first curve is fixed:
   E1 = 210e2 = [1,0,0,-1070,7812] (torsion Z/2 x Z/8, P8 = (244,3658), 4 P8 = (28,-14), T = (-36,18)) and
       E2 = E_{2,6}(t) : y^2 + (1 - c)xy - by = x^3 - bx^2,  c = (10 - 2t)/(t^2 - 9),  b = c + c^2,
   with P6 = (0,0) of order 6.  The gluing psi : E1[2] -> E2[2] sends 4 P8 to 3 P6 and the points with
   x = 31/4, -36 to those with x = -2(t - 1)^2/((t + 3)^2 (t - 3)), -(t - 1)^2 (t - 5)/(4 (t - 3)^2 (t + 3))
   (x-coordinates on y^2 = cubic, after completing the square).  The HLP sextic h has coefficients in Q(t), and
   x -> lambda x, y -> mu y (lambda, mu in Q(t)) takes it to the model
       C_t : y^2 = f_t(x)   (embedded below; f_t = F(x^2), coefficients of degree <= 13 in t).
   phi1 : (x,y) -> (x^2, y) maps C_t onto Ea : Y^2 = F(X), phi2 : (x,y) -> (1/x^2, y/x^3) onto
   Eb : Y^2 = X^3 F(1/X); below, Ea is Q(t)-isomorphic to E2 and Eb to E1.

   Argument.
   (1) Splitting.  Wherever f_t is squarefree of degree 6, phi1^* + phi2^* : Ea x Eb -> J is an isogeny over Q
       (the involutions (x,y) -> (-x,y) and (x,y) -> (-x,-y) compose to the hyperelliptic involution), so J is
       isogenous over Q to E1 x E2.
   (2) The marked group over Q(t).  D8 = phi2^*(P8) has order 8, D6 = phi1^*(P6) has order 6, 4 D8 = 3 D6 (the
       glued 2-torsion points), D24 = D8 + D6 has order 24, and DT = phi2^*(T) has order 2 and is not 12 D24
       (Mumford arithmetic over Q(t)).  So <D24, DT> = [2,24] embeds in J_t(Q)_tors at every smooth member.
   (3) Exactness on a residue class.  #J(F_p) = #E1(F_p) #E2(F_p) is always divisible by 192, and J(F_p) usually
       has 2-exponent >= 16 (over F_p the Kummer maps of E1 and E2 are both bijective onto the same group, so a
       half of some (P8, Q) is always F_p-rational in J), so the group structure at three primes is used.  For every 11-, 13- and 41-integral
       t = 24 (mod 5863), C_t has the reductions of C_24 at 11, 13 and 41, with J(F_11) = [4,48] (2-rank 2),
       J(F_13) = [2,2,48] (second 2-exponent 1) and J(F_41) = [2,2,24,24] (2-exponent 8; 3-rank 2, but
       3-exponent 3, and 3-rank 1 at 11).  J_t(Q)_tors embeds in all three (p odd, good reduction): its 2-part
       has rank <= 2 and exponents <= (8, 2), its odd part is <= Z/3, so it embeds in [2,24]: it is exactly
       [2,24].  (StructureBound: for every prime l the sorted l-primary exponents of a subgroup are bounded
       componentwise by those of the group.)
   (4) Infinitely many.  The G2-invariants of the members differ, so t -> M_2 is nonconstant with finite
       fibres, and the residue class gives infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup at t = 24, 5887 (in the class) and at t = 2, 7 (outside it) is [2,24].  (The member
       t = 21, outside the class, is the census curve 2.48.a, with torsion [2,48]: a 2-power gain of the gluing, see
       family_2_48_split_certificate.m.)
   Run (from the repository root):
       magma -b pipeline/magma/family_2_24_split_certificate.m > data/knowledge/sources/family_2_24_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
SetQuitOnError(true);                   // fail closed: a failed assertion stops the run before the marker
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

function MonEC(G)                       // Y^2 = G(X), G a cubic, as a Weierstrass curve (X' = lc X, Y' = lc Y)
    cs := [ Coefficient(G, i) : i in [0..3] ];
    return EllipticCurve([0, cs[3], 0, cs[2]*cs[4], cs[1]*cs[4]^2]), cs[4];
end function;
function MonicCubic(E, x)               // y^2 = MonicCubic(x) after completing the square
    a1, a2, a3, a4, a6 := Explode(aInvariants(E));
    return x^3 + (a2 + a1^2/4)*x^2 + (a4 + a1*a3/2)*x + (a6 + a3^2/4);
end function;
function HLPSextic(al, be)              // HLP 2000, Prop. 4, with psi(al_i) = be_i
    L := Parent(al[1]); PL<X> := PolynomialRing(L);
    aa1 := (al[3]-al[2])^2/(be[3]-be[2]) + (al[2]-al[1])^2/(be[2]-be[1]) + (al[1]-al[3])^2/(be[1]-be[3]);
    bb1 := (be[3]-be[2])^2/(al[3]-al[2]) + (be[2]-be[1])^2/(al[2]-al[1]) + (be[1]-be[3])^2/(al[1]-al[3]);
    aa2 := al[1]*(be[3]-be[2]) + al[2]*(be[1]-be[3]) + al[3]*(be[2]-be[1]);
    bb2 := be[1]*(al[3]-al[2]) + be[2]*(al[1]-al[3]) + be[3]*(al[2]-al[1]);
    Df := ((al[1]-al[2])*(al[1]-al[3])*(al[2]-al[3]))^2;
    Dg := ((be[1]-be[2])*(be[1]-be[3])*(be[2]-be[3]))^2;
    AA := Dg*aa1/aa2; BB := Df*bb1/bb2;
    return -(AA*(al[2]-al[1])*(al[1]-al[3])*X^2 + BB*(be[2]-be[1])*(be[1]-be[3]))
           *(AA*(al[3]-al[2])*(al[2]-al[1])*X^2 + BB*(be[3]-be[2])*(be[2]-be[1]))
           *(AA*(al[1]-al[3])*(al[3]-al[2])*X^2 + BB*(be[1]-be[3])*(be[3]-be[2]));
end function;
// pullbacks of points (X, Y) of the Weierstrass models (scaled by c) to the Jacobian of y^2 = F(x^2)
function Pull1(JJ, Pt, c) PP := Parent((JJ!0)[1]); xx := PP.1; XX := Pt[1]/c; YY := Pt[2]/c; return elt< JJ | xx^2 - XX, PP!YY >; end function;
function Pull2(JJ, Pt, c) PP := Parent((JJ!0)[1]); xx := PP.1; XX := Pt[1]/c; YY := Pt[2]/c; return elt< JJ | xx^2 - 1/XX, (YY/XX)*xx >; end function;
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
// f2 = h(lam x)/mu^2 for some lam, mu in the base field (h, f2 even sextics)
function SameUpToScaling(h, f2)
    okl, lam := IsSquare(Coefficient(f2, 2)*Coefficient(h, 0)/(Coefficient(f2, 0)*Coefficient(h, 2)));
    okm, mu := IsSquare(Coefficient(h, 0)/Coefficient(f2, 0));
    return okl and okm and &and[ Coefficient(f2, i)*mu^2 eq Coefficient(h, i)*lam^i : i in [0..6] ];
end function;

Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
f := (-10368*t^13 + 404352*t^12 - 6531840*t^11 + 58958080*t^10 - 353292672*t^9 + 1620756096*t^8 - 5595259392*t^7 + 8999155200*t^6 + 27131518080*t^5 - 192335486848*t^4 + 470196878592*t^3 - 580362151680*t^2 + 358681046400*t - 88035984000)*x^6 + (6561*t^12 - 236196*t^11 + 3424842*t^10 - 22566276*t^9 + 21810951*t^8 + 635091192*t^7 - 5459931668*t^6 + 35114636952*t^5 - 159717812289*t^4 + 350315165484*t^3 + 10188866922*t^2 - 1350417694356*t + 1633141477881)*x^4 + (-2827791*t^8 + 67866984*t^7 - 576869364*t^6 + 2364919416*t^5 - 7054491690*t^4 + 16618029336*t^3 + 145223799916*t^2 - 997272028536*t + 1517004081729)*x^2 + 293932800*t^4 - 3527193600*t^3 + 8817984000*t^2 - 89515238400*t + 246500755200;
assert &and[ Coefficient(f, i) eq 0 : i in [1, 3, 5] ] and Degree(f) eq 6 and IsSquarefree(f);
printf "C_t : y^2 = f_t(x) = %o\n", f;

// ---------------------------------------------------------------- (0) source: the HLP gluing of 210e2 and E_{2,6}(t)
E1 := EllipticCurve([Ft| 1, 0, 0, -1070, 7812]);
assert Conductor(EllipticCurve([Q| 1, 0, 0, -1070, 7812])) eq 210;
P8 := E1![244, 3658]; T1 := E1![28, -14]; T := E1![-36, 18];
assert 8*P8 eq E1!0 and 4*P8 ne E1!0 and 4*P8 eq T1 and 2*T eq E1!0 and T ne E1!0 and T ne T1;
cK := (10 - 2*t)/(t^2 - 9); bK := cK + cK^2;           // E_{2,6}(t): HLP Table 3, b = (-2t^3 + 14t^2 - 22t + 10)/((t+3)^2 (t-3)^2)
E2 := EllipticCurve([Ft| 1 - cK, -bK, -bK, 0, 0]);
P6 := E2![0, 0];
assert 6*P6 eq E2!0 and 2*P6 ne E2!0 and 3*P6 ne E2!0;
g1 := MonicCubic(E1, x); g2 := MonicCubic(E2, x);
al := [Ft| 28, 31/4, -36 ];               // x(4 P8), and the two other 2-torsion points of E1
be := [Ft| (10 - 2*t)/(t^2 - 9), -2*(t - 1)^2/((t + 3)^2*(t - 3)), -(t - 1)^2*(t - 5)/(4*(t - 3)^2*(t + 3)) ];
assert be[1] eq (3*P6)[1];                // psi(4 P8) = 3 P6: the special points are glued
assert &and[ Evaluate(g1, a) eq 0 : a in al ] and &and[ Evaluate(g2, b) eq 0 : b in be ] and #Seqset(be) eq 3;
h := P!HLPSextic(al, be);
assert SameUpToScaling(h, f);
print "(0) the HLP (Prop. 4) gluing of E1 = 210e2 and E2 = E_{2,6}(t) along psi (4 P8 -> 3 P6) has coefficients in Q(t),";
print "    f_t = h(lambda x)/mu^2 with lambda, mu in Q(t)";

// ---------------------------------------------------------------- (1) the two elliptic quotients
F := P![ Coefficient(f, 2*i) : i in [0..3] ];
Fr := P![ Coefficient(F, 3 - i) : i in [0..3] ];
Ea, ca := MonEC(F); Eb, cb := MonEC(Fr);
oka, ma := IsIsomorphic(E2, Ea);
okb, mb := IsIsomorphic(E1, Eb);
assert oka and okb;
print "(1) phi1 = (x^2, y) : C_t -> Ea : Y^2 = F(X), phi2 = (1/x^2, y/x^3) : C_t -> Eb : Y^2 = X^3 F(1/X);";
print "    Ea is Q(t)-isomorphic to E2 and Eb to E1 = 210e2: J is isogenous over Q to E1 x E2";

// ---------------------------------------------------------------- (2) the marked group over Q(t)
J := Jacobian(HyperellipticCurve(f));
D8 := Pull2(J, mb(P8), cb);
D6 := Pull1(J, ma(P6), ca);
DT := Pull2(J, mb(T), cb);
assert 8*D8 eq J!0 and 4*D8 ne J!0;
assert 6*D6 eq J!0 and 3*D6 ne J!0 and 2*D6 ne J!0 and 4*D8 eq 3*D6;
D24 := D8 + D6;
assert 24*D24 eq J!0 and 12*D24 ne J!0 and 8*D24 ne J!0;
assert 2*DT eq J!0 and DT ne J!0 and DT ne 12*D24;
print "(2) over Q(t): D8 = phi2^*(P8) has order 8, D6 = phi1^*(P6) has order 6, 4 D8 = 3 D6, D24 = D8 + D6 has order 24, DT = phi2^*(T) has order 2 and is not 12 D24: [2,24] is marked";

// ---------------------------------------------------------------- (3) exactness on the residue class t = 24 mod 5863
primes := [11, 13, 41]; t0 := 24; modulus := 5863;
members := [Q| 24, 5887, -5839, 11750];
invs := [];
for p in primes do
    r0 := Integers()!(GF(p)!t0);
    fbar := ReductionClass(f, p, r0);
    for m in members do
        assert Valuation(m, p) ge 0 and GF(p)!m eq GF(p)!r0;
        assert PolynomialRing(GF(p))!SpecPoly(f, m) eq fbar;
    end for;
    n, chi := FrobData(fbar);
    inv := Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
    assert &*inv eq n;
    Append(~invs, inv);
    printf "(3) p = %o, t = %o mod %o: reduction %o (squarefree, degree %o); J(F_%o) = %o, chi = %o\n", p, r0, p, fbar, Degree(fbar), p, inv, chi;
end for;
assert invs eq [[4, 48], [2, 2, 48], [2, 2, 24, 24]];
bnd := StructureBound(invs);
printf "    the largest group embedding in all J(F_p) has invariants %o\n", bnd;
assert bnd eq [2, 24];
printf "    hence every t = %o mod %o, integral at %o, has J_t(Q)_tors exactly [2,24]\n", t0, modulus, primes;

// ---------------------------------------------------------------- (4) nonconstant moduli map, (5) sanity
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } eq #members;
printf "(4) the G2-invariants of the members %o are pairwise distinct\n", members;
for m in [Q| 24, 5887] do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(5) TorsionSubgroup at t = %o (in the class): %o\n", m, invm;
    assert invm eq [2, 24];
end for;
for m in [Q| 2, 7] do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(5) TorsionSubgroup at t = %o (outside the class): %o\n", m, invm;
    assert invm eq [2, 24];
end for;
printf "time %o s\n", Cputime(T0);
print "FAMILY_2_24_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E1 x E2) and J(Q)_tors exactly [2,24]";
quit;
