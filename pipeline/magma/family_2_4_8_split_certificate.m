/* family_2_4_8_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E1 x E2) with J(Q)_tors exactly [2,4,8].
   Found 2026-10-19 (torsion_jac notes/hlpP2_2026_10_19.md).

   Source.  The Z/2 x Z/4 x Z/8 entry of Howe-Leprevost-Poonen 2000, Table 1 (section 3.4): the (2,2)-gluing
   (Prop. 4) of two curves with torsion Z/2 x Z/8 along a bijection of their (rational) 2-torsion that matches
   the points 4 P8 (the other choices give [8,8]), restricted to the line on which the first curve is fixed:
   E1 = 210e2 = [1,0,0,-1070,7812] (torsion Z/2 x Z/8, P8 = (244,3658), 4 P8 = (28,-14), T = (-36,18)) and
       E2 = E_{2,8}(t) : y^2 + (1 - c)xy - by = x^3 - bx^2,
            b = (16t^3 + 16t^2 + 6t + 1)/(8t^2 - 1)^2,  c = (16t^3 + 16t^2 + 6t + 1)/(2t(4t + 1)(8t^2 - 1)),
   with Q8 = (0,0) of order 8.  The gluing psi : E1[2] -> E2[2] sends 4 P8 to 4 Q8, and the points with
   x = 31/4, -36 to the points with the x-coordinates listed in be[2], be[3] below (on y^2 = cubic, after
   completing the square).  The HLP sextic h has coefficients in Q(t), and x -> lambda x, y -> mu y
   (lambda, mu in Q(t)) takes it to the model
       C_t : y^2 = f_t(x)   (embedded below; f_t = F(x^2), coefficients of degree <= 28 in t).
   phi1 : (x,y) -> (x^2, y) maps C_t onto Ea : Y^2 = F(X), phi2 : (x,y) -> (1/x^2, y/x^3) onto
   Eb : Y^2 = X^3 F(1/X); below, Ea is Q(t)-isomorphic to E2 and Eb to E1.

   Argument.
   (1) Splitting.  Wherever f_t is squarefree of degree 6, phi1^* + phi2^* : Ea x Eb -> J is an isogeny over Q
       (the involutions (x,y) -> (-x,y) and (x,y) -> (-x,-y) compose to the hyperelliptic involution), so J is
       isogenous over Q to E1 x E2.
   (2) The marked group over Q(t).  D8 = phi2^*(P8) and D8' = phi1^*(Q8) have order 8 and 4 D8 = 4 D8' (the
       glued 2-torsion points); D8' - D8 has order 4 and 2(D8' - D8) is not 4 D8, so <D8, D8'> = [4,8]; DT =
       phi2^*(T) has order 2 and lies outside <D8, D8'> (it is none of its three elements of order 2)
       (Mumford arithmetic over Q(t)).  So [2,4,8] embeds in J_t(Q)_tors at every smooth member.
   (3) Exactness on a residue class.  #J(F_p) = #E1(F_p) #E2(F_p) is always divisible by 256, so the group
       structure is used.  For every 23- and 41-integral t = 3 (mod 943), C_t has the reductions of C_3 at 23
       and 41, with J(F_23) = [2,2,4,32] (2-exponents 5,2,1,1; no odd part) and J(F_41) = [8,8,24]
       (2-exponents 3,3,3).  J_t(Q)_tors embeds in both (p odd, good reduction): it is a 2-group of 2-rank <= 3
       with exponents <= (3,2,1), so it embeds in [2,4,8]: it is exactly [2,4,8].  (StructureBound: for every
       prime l the sorted l-primary exponents of a subgroup are bounded componentwise by those of the group.)
   (4) Infinitely many.  The G2-invariants of the members differ, so t -> M_2 is nonconstant with finite
       fibres, and the residue class gives infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup at t = 3, 946 (in the class) and at t = 2, 5 (outside it) is [2,4,8].
   Run (from the repository root):
       magma -b pipeline/magma/family_2_4_8_split_certificate.m > data/knowledge/sources/family_2_4_8_split_certificate.log */
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
f := (-21743271936*t^28 - 152202903552*t^27 - 497377345536*t^26 - 1006985281536*t^25 - 1407743688704*t^24 - 1424064774144*t^23 - 1048389877760*t^22 - 531381092352*t^21 - 136342732800*t^20 + 44849299456*t^19 + 73012264960*t^18 + 44738224128*t^17 + 17073660160*t^16 + 3764438016*t^15 - 46024960*t^14 - 431994112*t^13 - 207600800*t^12 - 62016064*t^11 - 13420240*t^10 - 2181072*t^9 - 265919*t^8 - 23712*t^7 - 1464*t^6 - 56*t^5 - t^4)*x^6 + (6879707136*t^24 + 41278242816*t^23 + 113515167744*t^22 + 189191946240*t^21 + 213339930624*t^20 + 172767707136*t^19 + 105001795584*t^18 + 50333736960*t^17 + 19963194880*t^16 + 5740656640*t^15 - 690966016*t^14 - 2791641600*t^13 - 2384416192*t^12 - 1235939200*t^11 - 412298144*t^10 - 72060480*t^9 + 7284430*t^8 + 9541920*t^7 + 3470936*t^6 + 793768*t^5 + 126289*t^4 + 14080*t^3 + 1056*t^2 + 48*t + 1)*x^4 + (-723914496*t^16 - 2895657984*t^15 - 5067401472*t^14 - 5067401472*t^13 - 3055317408*t^12 - 931424832*t^11 + 116545392*t^10 + 269216784*t^9 + 143529985*t^8 + 51628448*t^7 + 17550328*t^6 + 6045464*t^5 + 1711777*t^4 + 336896*t^3 + 42112*t^2 + 3008*t + 94)*x^2 + 18370800*t^8 + 36741600*t^7 + 27556200*t^6 + 9185400*t^5 - 2480625*t^4 - 3628800*t^3 - 1360800*t^2 - 226800*t - 14175;
assert &and[ Coefficient(f, i) eq 0 : i in [1, 3, 5] ] and Degree(f) eq 6 and IsSquarefree(f);
printf "C_t : y^2 = f_t(x) = %o\n", f;

// ---------------------------------------------------------------- (0) source: the HLP gluing of 210e2 and E_{2,8}(t)
E1 := EllipticCurve([Ft| 1, 0, 0, -1070, 7812]);
assert Conductor(EllipticCurve([Q| 1, 0, 0, -1070, 7812])) eq 210;
P8 := E1![244, 3658]; T1 := E1![28, -14]; T := E1![-36, 18];
assert 8*P8 eq E1!0 and 4*P8 ne E1!0 and 4*P8 eq T1 and 2*T eq E1!0 and T ne E1!0 and T ne T1;
bK := (16*t^3 + 16*t^2 + 6*t + 1)/(8*t^2 - 1)^2; cK := (16*t^3 + 16*t^2 + 6*t + 1)/(2*t*(4*t + 1)*(8*t^2 - 1));   // E_{2,8}(t): HLP Table 3
E2 := EllipticCurve([Ft| 1 - cK, -bK, -bK, 0, 0]);
Q8 := E2![0, 0];
assert 8*Q8 eq E2!0 and 4*Q8 ne E2!0;
g1 := MonicCubic(E1, x); g2 := MonicCubic(E2, x);
al := [Ft| 28, 31/4, -36 ];               // x(4 P8), and the two other 2-torsion points of E1
be := [Ft| (1/4*t^3 + 3/16*t^2 + 1/32*t)/(t^4 - 1/4*t^2 + 1/64), (1/4*t^3 + 3/16*t^2 + 1/16*t + 1/128)/(t^4 - 1/8*t^2),
              (-1/4*t^4 - 1/4*t^3 - 3/32*t^2 - 1/64*t)/(t^4 + 1/2*t^3 - 1/16*t^2 - 1/16*t - 1/128) ];
assert be[1] eq (4*Q8)[1];                // psi(4 P8) = 4 Q8: the special points are glued
assert &and[ Evaluate(g1, a) eq 0 : a in al ] and &and[ Evaluate(g2, b) eq 0 : b in be ] and #Seqset(be) eq 3;
h := P!HLPSextic(al, be);
assert SameUpToScaling(h, f);
print "(0) the HLP (Prop. 4) gluing of E1 = 210e2 and E2 = E_{2,8}(t) along psi (4 P8 -> 4 Q8) has coefficients in Q(t),";
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
D8p := Pull1(J, ma(Q8), ca);
DT := Pull2(J, mb(T), cb);
assert 8*D8 eq J!0 and 4*D8 ne J!0 and 8*D8p eq J!0 and 4*D8p ne J!0 and 4*D8 eq 4*D8p;
D4 := D8p - D8;
assert 4*D4 eq J!0 and 2*D4 ne J!0 and 2*D4 ne 4*D8;      // <D8, D8'> = <D8> + <D4> = [4,8]
assert 2*DT eq J!0 and DT ne J!0 and DT notin { 4*D8, 2*D4, 4*D8 + 2*D4 };
print "(2) over Q(t): D8 = phi2^*(P8), D8' = phi1^*(Q8) have order 8, 4 D8 = 4 D8', D8' - D8 has order 4 with 2(D8' - D8) != 4 D8, DT = phi2^*(T) of order 2 outside <D8, D8'>: [2,4,8] is marked";

// ---------------------------------------------------------------- (3) exactness on the residue class t = 3 mod 943
primes := [23, 41]; t0 := 3; modulus := 943;
members := [Q| 3, 946, -940, 1889];
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
assert invs eq [[2, 2, 4, 32], [8, 8, 24]];
bnd := StructureBound(invs);
printf "    the largest group embedding in all J(F_p) has invariants %o\n", bnd;
assert bnd eq [2, 4, 8];
printf "    hence every t = %o mod %o, integral at %o, has J_t(Q)_tors exactly [2,4,8]\n", t0, modulus, primes;

// ---------------------------------------------------------------- (4) nonconstant moduli map, (5) sanity
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } eq #members;
printf "(4) the G2-invariants of the members %o are pairwise distinct\n", members;
for m in [Q| 3, 946] do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(5) TorsionSubgroup at t = %o (in the class): %o\n", m, invm;
    assert invm eq [2, 4, 8];
end for;
for m in [Q| 2, 5] do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(5) TorsionSubgroup at t = %o (outside the class): %o\n", m, invm;
    assert invm eq [2, 4, 8];
end for;
printf "time %o s\n", Cputime(T0);
print "FAMILY_2_4_8_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E1 x E2) and J(Q)_tors exactly [2,4,8]";
quit;
