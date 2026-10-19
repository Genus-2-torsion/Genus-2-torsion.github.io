/* family_2_48_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E1 x E2) with J(Q)_tors exactly [2,48].
   Found 2026-10-19 (torsion_jac notes/hlpP2_2026_10_19.md), as a by-product of the [2,24] certificate.

   Source.  The [2,24] family of family_2_24_split_certificate.m: the Howe-Leprevost-Poonen (2,2)-gluing
   (HLP 2000, Prop. 4) of E1 = 210e2 = [1,0,0,-1070,7812] (torsion Z/2 x Z/8, P8 = (244,3658), T = (-36,18)) and
   E2 = E_{2,6}(t) (torsion Z/2 x Z/6) along psi : alpha_i -> beta_i, alpha = (28, 31/4, -36), beta as in be[]
   below (psi(4 P8) = 3 P6, psi(T) = T3' := (beta_3, 0)), model C_t : y^2 = f_t(x).
   The gain.  By HLP 2000, Prop. 12, a point (P0, Q0) of E1 x E2 with 2P0 = a := P8 + T (order 8) and
   2Q0 = T3' maps to a rational point D of J = (E1 x E2)/graph(psi) iff the Kummer images agree under psi:
       (x(a) - alpha_k)(beta_3 - beta_k) is a square for k = 1, 2
   (the third component follows, the products of the three components being squares on both sides).  Over Q(t),
   (x(a) - alpha_1)(beta_3 - beta_1) = s1(t) h1(t)^2 and (x(a) - alpha_2)(beta_3 - beta_2) = s2(t) h2(t)^2 with
       s1 = 6t^2 - 12t - 90,   s2 = 15t^2 - 150t + 135,
   so the condition is the genus-1 curve X : Y1^2 = s1(t), Y2^2 = s2(t) (an intersection of two quadrics in
   P^3), at the t where h1, h2 are finite and nonzero.  The map
       phi : E' -> X,  E' = 4800c1 : y^2 = x^3 - x^2 + 67x + 237  (rank 1, E'(Q) = Z/2 x Z, P = (13, -56)),
   given by four quartics (embedded below; checked to map E' into X), sends T' = (-3, 0) to a point over t = 21.
   2D = image of (a, T3') = image of (P8 + T + T, 0) = D8 (T3' = psi(T) is glued to T), of order 8; so D has order
   16, 8 D = 4 D8 = image of (4 P8, 0) is not DT = image of (T, 0), and <D, DT, D6> = [2,48].

   Argument.
   (1) Splitting and (2) the marked [2,24] over Q(t): as in family_2_24_split_certificate.m (repeated below).
   (3) The gain locus.  The two Kummer products are s_k h_k^2 over Q(t) with the s_k above (checked), and the
       third is automatic; phi maps E' into X (the quadrics pull back to multiples of the Weierstrass cubic).
       At every member below the two Kummer products are checked to be nonzero rational squares.
   (4) Exactness.  For every 11- and 13-integral t = 21 (mod 143), C_t has the reductions of C_21 at 11 and 13,
       with J(F_11) = [4,48] (2-rank 2) and J(F_13) = [2,2,48] (second 2-exponent 1); J_t(Q)_tors embeds in both,
       hence in [2,48].
   (5) Infinitely many members in the class.  E' has good reduction at 11 and 13, and P has orders 8 and 4 there,
       so R_k = T' + 8kP reduces to T' mod 11 and 13.  The quartics are integral and the fourth is a unit at T'
       mod 11 and 13, so t(R_k) is 11- and 13-integral and = 21 (mod 143).  On that class h1, h2 are finite and
       nonzero and s1, s2 do not vanish (their values at 21 are units mod 11 and 13), so every R_k gives a member
       with the gain and J(Q)_tors exactly [2,48].  The R_k are distinct and t has finite degree on E', so
       infinitely many t occur; the G2-invariants at t(T'), t(T' + 8P), t(T' - 8P) differ: infinitely many
       isomorphism classes.
   (6) Sanity: TorsionSubgroup at t = 21 (in the class; Q-isomorphic to the census curve 2.48.a) and at
       t(P) = -85701/20879 (on the gain locus, outside the class) is [2,48].
   Run (from the repository root):
       magma -b pipeline/magma/family_2_48_split_certificate.m > data/knowledge/sources/family_2_48_split_certificate.log */
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

// ---------------------------------------------------------------- (3) the gain locus
aa := P8 + T;
assert 8*aa eq E1!0 and 4*aa ne E1!0;
v1 := [ aa[1] - al[k] : k in [1..3] ];
v2 := [ be[3] - be[1], be[3] - be[2], (be[3] - be[1])*(be[3] - be[2]) ];
assert IsSquare(&*v1) and IsSquare(&*v2);                // so the third components agree once the first two do
s1 := 6*t^2 - 12*t - 90; s2 := 15*t^2 - 150*t + 135;
ok1, h1 := IsSquare(v1[1]*v2[1]/s1); ok2, h2 := IsSquare(v1[2]*v2[2]/s2);
assert ok1 and ok2;
print "(3) Kummer products: (x(a) - alpha_k)(beta_3 - beta_k) = s_k h_k^2 over Q(t), s1 = 6t^2 - 12t - 90, s2 = 15t^2 - 150t + 135;";
print "    the third component is automatic (the products of the three components are squares)";
P3<TT, Y1, Y2, ZZ> := ProjectiveSpace(Q, 3);
q1X := Y1^2 - (6*TT^2 - 12*TT*ZZ - 90*ZZ^2);
q2X := Y2^2 - (15*TT^2 - 150*TT*ZZ + 135*ZZ^2);
X := Curve(P3, [q1X, q2X]);
assert Genus(X) eq 1;
Ep := EllipticCurve([Q| 0, -1, 0, 67, 237]);
assert Conductor(Ep) eq 4800;
PE<xe, ye, ze> := CoordinateRing(Ambient(Ep));
Wp := DefiningPolynomial(Ep);
dpz := [ 825440*xe^3*ye - 88761600*xe^3*ze + 6466000*xe^2*ye*ze - 255840*xe*ye^2*ze - 814247*ye^3*ze + 148401600*xe^2*ze^2 + 175754720*xe*ye*ze^2 + 127819380*ye^2*ze^2 - 9168427200*xe*ze^3 + 246404040*ye*ze^3 - 24527959200*ze^4,
         -3611520*xe^3*ye + 282816000*xe^3*ze + 21520320*xe^2*ye*ze - 639600*xe*ye^2*ze + 3637104*ye^3*ze + 910368000*xe^2*ze^2 - 226107840*xe*ye*ze^2 - 226782000*ye^2*ze^2 + 20653056000*xe*ze^3 - 738253440*ye*ze^3 + 46065888000*ze^4,
         28142400*xe^2*ye*ze - 1023360*xe*ye^2*ze + 31980*ye^3*ze - 1874880000*xe^2*ze^2 - 236853600*xe*ye*ze^2 - 94686480*ye^2*ze^2 - 596880000*xe*ze^3 - 2016662400*ye*ze^3 - 4086720000*ze^4,
         533*ye^3*ze + 11520900*ye^2*ze^2 - 1278360000*xe*ze^3 - 84693000*ye*ze^3 - 3515580000*ze^4 ];
for q in [q1X, q2X] do
    assert IsDivisibleBy(Evaluate(q, dpz), Wp);
end for;
print "    X : Y1^2 = s1(t), Y2^2 = s2(t) (genus 1); phi : E' = [0,-1,0,67,237] (4800c1) -> X by four integral quartics:";
print "    both quadrics pull back to multiples of the Weierstrass cubic";
Pg := Ep![13, -56]; Tp := Ep![-3, 0];
assert Order(Pg) eq 0 and 2*Tp eq Ep!0 and Tp ne Ep!0;
function TVal(R) v := Eltseq(R); return Evaluate(dpz[1], v)/Evaluate(dpz[4], v); end function;
assert Evaluate(dpz[4], Eltseq(Tp)) eq 319500000 and TVal(Tp) eq 21;
printf "    P = %o has infinite order, T' = %o has order 2; t(T') = %o\n", Pg, Tp, TVal(Tp);
function KumOK(t0)       // both Kummer products are nonzero rational squares at t = t0
    return &and[ (ev ne 0 and IsSquare(ev)) where ev := Evaluate(c, t0) : c in [v1[1]*v2[1], v1[2]*v2[2]] ];
end function;

// ---------------------------------------------------------------- (4) exactness at 11 and 13, (5) the members R_k = T' + 8kP
primes := [11, 13]; t0 := 21;
invs := [];
for p in primes do
    fbar := ReductionClass(f, p, t0);
    inv := Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
    n, chi := FrobData(fbar);
    Append(~invs, inv);
    printf "(4) p = %o, t = %o mod %o: reduction %o (squarefree, degree %o); J(F_%o) = %o, chi = %o\n", p, t0, p, fbar, Degree(fbar), p, inv, chi;
    assert Conductor(Ep) mod p ne 0;
    assert Integers()!Evaluate(dpz[4], Eltseq(Tp)) mod p ne 0;
    for c in [h1, h2, s1, s2] do Nc, Dc := IntNumDen(c); assert GF(p)!Evaluate(Nc, t0) ne 0 and GF(p)!Evaluate(Dc, t0) ne 0; end for;
end for;
assert invs eq [[4, 48], [2, 2, 48]];
bnd := StructureBound(invs);
printf "    the largest group embedding in all J(F_p) has invariants %o\n", bnd;
assert bnd eq [2, 48];
assert [ Order(ChangeRing(Ep, GF(p))!Pg) : p in primes ] eq [8, 4];
assert &and[ IsCoercible(Integers(), cc) : cc in &cat[ Coefficients(d) : d in dpz ] ];
print "(5) E' has good reduction at 11, 13 and P has orders 8, 4 there: R_k = T' + 8kP reduces to T' mod 11 and 13, the";
print "    quartics are integral and the fourth is a unit at T', so t(R_k) is 11- and 13-integral and = 21 mod 143;";
print "    h1, h2, s1, s2 are units at t = 21 mod 11 and 13, so the Kummer products are nonzero squares at every t(R_k)";
memb := [];
for k in [0, 1, -1] do
    tv := TVal(Tp + 8*k*Pg);
    for p in primes do
        assert Valuation(tv, p) ge 0 and GF(p)!tv eq GF(p)!t0;
        assert PolynomialRing(GF(p))!SpecPoly(f, tv) eq ReductionClass(f, p, t0);
    end for;
    assert KumOK(tv);
    Append(~memb, tv);
    printf "    k = %o: t(R_k) has %o digits, = 21 mod 143, both Kummer products are nonzero squares\n", k, #Sprint(tv);
end for;
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, tv))) : tv in memb } eq #memb;
print "    the G2-invariants at t(T'), t(T' + 8P), t(T' - 8P) are pairwise distinct: infinitely many isomorphism classes";

// ---------------------------------------------------------------- (6) sanity
C248a := HyperellipticCurve(PolynomialRing(Q)![1, -13, 19, 111, 19, -13, 1], PolynomialRing(Q)![0, 1, 1]);
assert IsIsomorphic(HyperellipticCurve(SpecPoly(f, 21)), C248a);
print "(6) the member t = 21 (in the class) is Q-isomorphic to the census curve 2.48.a";
for R in [Tp, Pg] do
    tv := TVal(R);
    assert KumOK(tv);
    invm := TorsionInvariants(SpecPoly(f, tv));
    printf "(6) t = %o: both Kummer products are squares; TorsionSubgroup %o\n", tv, invm;
    assert invm eq [2, 48];
end for;
printf "time %o s\n", Cputime(T0);
print "FAMILY_2_48_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E1 x E2) and J(Q)_tors exactly [2,48]";
quit;
