/* seed_credit_check.m -- which published curve is each seed curve whose credit names a paper?
   (credit errata of 2026-10-20 in data/knowledge/sources/PROVENANCE.md; the check of 2026-10-14 is
   attribution_check.m)

   Compares census curves (models as in data/curves/<id>.json) with the curves and one-parameter families
   printed in the cited sources, by Magma's IsIsomorphic over Q.  Every transcribed curve is checked against
   the torsion its source claims.  Equations of the sources (read from the page scans):
     PP2012a = Platonov–Petrunin, Dokl. Math. 85 (2012), 286–288: Thm 2 (p. 287) f14, f18; Thm 3 (p. 288)
               f11, f13, f17, f19, f23, f29; Thm 4 (p. 288) f28.
     PP2012b = Platonov–Petrunin, Dokl. Math. 86 (2012), 642–643: Thm 2 (p. 643) f36; Thm 3 f48,1, f48,2.
     PZP2013 = Platonov–Zhgun–Petrunin, Dokl. Math. 87 (2013), 318–321: p. 320 f28,2, f36 (and f48,1, f48,2).
     PP2015  = Platonov–Petrunin, Dokl. Math. 91 (2015), 220–221: p. 221 f28,1–f28,4, f36, f30,1–f30,4, f39.
     Elkies 2002 (web page g2_tors.html, "N=32: an example"): y^2 = (15x-1)(1056x^4+156183x^3+26297x^2+649x-121).
     BFT2014 = Bruin–Flynn–Testa, Acta Arith. 165 (2014), arXiv:1401.0580v1: Thm 6 (p. 7), C_rst : y^2 = G1^2 + 4s H1^3,
               H1 = x^2 + rx + t, G1 = (s-st-1)x^3 + 3s(r-t)x^2 + 3sr(r-t)x - st^2 + sr^3 + t; Examples 19-21
               (pp. 15-16) are C_{-3,-3,-3}, C_{-2,1,2}, C_{2,-1,-2}.  (Example 19 is printed with (12x^3 - 105)^2;
               Thm 6 gives (13x^3 - 105)^2, and only the latter has (Z/3)^2.)
     Leprevost 1995 = J. Theor. Nombres Bordeaux 7 (1995), 283-306: the [3,9] family C_t (p. 303, Thm 1.2.2)
               f_t = x^6 + 6(t-1)x^5 + 3(t-1)(3t-5)x^4 - (18t^2-34t+18)x^3 + 3(t-1)(5t-3)x^2 - 6t(t-1)x + t^2,
               and its member C_3 (p. 303).  P. 290: the order-21 curves of Leprevost 1991 (orders 15-21; his [2])
               all have a rational Weierstrass point.
     Flynn 1991 = Invent. Math. 106 (1991), 433-442: Result 1(b) (p. 435), g = 2, r = 1:
               Y^2 = (Y - X^3)(X^2 + w1 X + w0), D = (0,0) - oo of exact order 6, oo a rational Weierstrass point;
               Result 3 (p. 438), g = 2: Y^2 + t(X-1)Y = psi(X)^2 - t(X^4 + X^(r+1)), psi = X^(r+1) + ... + X^3,
               a divisor of exact order 11 - r (r = 1: order 10; r = 2: order 9).
     Leprevost 1997 = Manuscripta Math. 92 (1997), 47-63: C_{1,t} (p. 48), g = 2:
               y^2 = [x((x-1)^2 - x^2) + t((x-1)^2 + x^2)]^2 + 4t x^5, a rational point of order 20.
   A family search finds the rational t at which the family has the Igusa-Clebsch point of the census curve
   and tests isomorphism over Q there; each search is first run on a known member moved by x -> (2x+1)/(x-3).
   Run:  magma -b pipeline/magma/seed_credit_check.m < /dev/null > data/knowledge/sources/seed_credit_check.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
Q := Rationals();
P<x> := PolynomialRing(Q);
Tors := func<C | Invariants(TorsionSubgroup(Jacobian(SimplifiedModel(MinimalWeierstrassModel(C)))))>;

// census curves (data/curves/<id>.json)
CEN := [*
    <"5.a", [-6,11,-19,14,-9,1], [1]>,
    <"6.a", [3,21,46,26,-12,1], [0,1,1]>,
    <"7.a", [-2,3,0,-3,0,1], [0,0,0,1]>,
    <"9.a", [0,0,0,0,-1], [1,1,0,1]>,
    <"10.a", [7,16,0,-8,-2,1], [0,1,0,1]>,
    <"11.a", [0,0,1], [1,1,0,1]>,
    <"14.a", [0,1,1], [1,0,0,1]>,
    <"17.a", [1,3,1,1], [0,1,1,1]>,
    <"18.a", [0,2,4,0,-2], [1,0,0,1]>,
    <"19.a", [-2,5,-3,-5,1,1,1], [1,1,1]>,
    <"19.b", [0,1,1], [1,0,1,1]>,
    <"20.a", [0,0,0,0,0,-1], [1,1]>,
    <"21.a", [0,1,2,0,-1], [1,1,0,1]>,
    <"23.a", [1,-2,2,1,-2,-1,1], [0,0,1]>,
    <"28.a", [2,3,1,1,0,-1], [1,0,0,1]>,
    <"28.b", [4,10,-4,-1,5,-3,1], [0,1,1]>,
    <"29.a", [0,0,-1,2,0,-2,1], [1,1]>,
    <"30.a", [0,-1,16,-12,-10,4,2], [0,0,1,1]>,
    <"32.a", [-47,-49,49,45,1,-1,1], [0,1,1]>,
    <"36.a", [9,6,-14,13,30,3], [0,0,1]>,
    <"36.b", [9,-9,-1,-5,5,-1,1], [0,1,1]>,
    <"48.a", [3,0,9,0,-4], [0,1,0,1]>,
    <"2.28.a", [0,-3,4,34,-42,-3,9], [1,0,1]>,
    <"3.3.a", [3,3,7,2,4,0,1], [1]>,
    <"3.9.a", [326,-75,-230,-113,40,75,16], [1,0,0,1]>,
    <"3.9.c", [1,-3,-1,7,6,4,2], [0,1,1,1]> *];
CC := [ <c[1], HyperellipticCurve(P!c[2], P!c[3])> : c in CEN ];
Cen := func< id | [c[2] : c in CC | c[1] eq id][1] >;

// (1) printed curves: torsion (against the source's claim) and the census curves isomorphic to them over Q
function Crst(r, s, t)
    H1 := x^2 + r*x + t;
    G1 := (s - s*t - 1)*x^3 + 3*s*(r - t)*x^2 + 3*s*r*(r - t)*x - s*t^2 + s*r^3 + t;
    return G1^2 + 4*s*H1^3;
end function;
SRC := [*
 <"PP2012a f14 (claim: order 14)", (x^2+1)*(x^4+5*x^2+4*x+4)>,
 <"PP2012a f18 (claim: order 18)", (x^2-x+1)*(x^4-x^3+9*x^2+8*x-8)>,
 <"PP2012a f11 (claim: order 11)", x^6+2*x^4+2*x^3+5*x^2+2*x+1>,
 <"PP2012a f13 (claim: order 13)", (x+2)*(x^5-2*x^4-4*x^2+12*x-4)>,
 <"PP2012a f17 (claim: order 17)", (x+2)*(x^5+3*x^3+5*x+2)>,
 <"PP2012a f19 (claim: order 19)", x^6+2*x^5+x^4+2*x^3+6*x^2+4*x+1>,
 <"PP2012a f23 (claim: order 23)", 4*x^6+4*x^5-7*x^4-4*x^3+8*x^2+8*x+4>,
 <"PP2012a f29 (claim: order 29)", (x-1)*(x^5-x^4-4*x^3-12*x^2-12*x-4)>,
 <"PP2012a f28 (claim: order 28)", x^6+2*x^5-5*x^4-14*x^3-3*x^2+24*x+28>,
 <"PP2012b f36 (claim: order 36)", (x+2)*(2*x^2+x+2)*(2*x^3+9*x^2+5)>,
 <"PP2012b f48,1 (claim: order 48)", (x-2)*(x+2)*(x^4-10*x^2-3)>,
 <"PP2012b f48,2 (claim: order 48)", (x^2+x+2)*(x^4+3*x^3+21*x^2-27*x+18)>,
 <"PZP2013 f28,2 (claim: order 28)", (x+1)*(4*x^2-3*x+1)*(4*x^3+9*x^2-4*x+4)>,
 <"PZP2013 f36 (claim: order 36)", (x+2)*(3*x^2-6*x+8)*(3*x^3+6*x^2+3*x+4)>,
 <"PP2015 f28,1 (claim: order 28)", (9*x^4+84*x^3+229*x^2+142*x+32)*(x-2)*x>,
 <"PP2015 f28,2 (claim: order 28)", (4*x^3+21*x^2-27)*(4*x+5)*(x+2)*x>,
 <"PP2015 f28,3 (claim: order 28)", (4*x^4+8*x^3+24*x^2+4*x-7)*(x+2)*(x-2)>,
 <"PP2015 f28,4 (claim: order 28)", (2*x^3+13*x^2+8*x+2)*(x^2+4*x-6)*(2*x-3)>,
 <"PP2015 f36 (claim: order 36)", (3*x^4+28*x^3+88*x^2+56*x-52)*(3*x-2)*(x+2)>,
 <"PP2015 f30,1 (claim: order 30)", (16*x^4+56*x^3+57*x^2+11*x-8)*(x+3)*x>,
 <"PP2015 f30,2 (claim: order 30)", (16*x^4+72*x^3-111*x^2-117*x+144)*(x+3)*x>,
 <"PP2015 f30,3 (claim: order 30)", 2*(4*x^4+34*x^3+66*x^2-44*x-71)*(2*x-1)*(x+1)>,
 <"PP2015 f30,4 (claim: order 30)", (4*x^4+24*x^3+25*x^2-69*x+73)*(x^2-x+1)>,
 <"PP2015 f39 (claim: order 39)", (9*x^5+36*x^4-12*x^3+48*x^2+4*x+16)*(x+4)>,
 <"Elkies 2002 N=32 (claim: order 32)", (15*x-1)*(1056*x^4+156183*x^3+26297*x^2+649*x-121)>,
 <"BFT2014 Ex. 19, C_{-3,-3,-3} (claim: (Z/3)^2)", Crst(-3,-3,-3)>,
 <"BFT2014 Ex. 19 as printed, (12x^3-105)^2 - 12(x^2-3x-3)^3", (12*x^3-105)^2 - 12*(x^2-3*x-3)^3>,
 <"BFT2014 Ex. 20, C_{-2,1,2} (claim: (Z/3)^2)", Crst(-2,1,2)>,
 <"BFT2014 Ex. 21, C_{2,-1,-2} as printed (claim: (Z/3)^2)", 12*x^6+72*x^5+312*x^4+688*x^3+768*x^2+192*x+68>,
 <"Leprevost 1995 C_3 (claim: [3,9])", x^6+12*x^5+24*x^4-78*x^3+72*x^2-36*x+9> *];
for s in SRC do
    C := HyperellipticCurve(s[2]);
    g := G2Invariants(C);
    iso := [ c[1] : c in CC | G2Invariants(c[2]) eq g and IsIsomorphic(C, c[2]) ];
    printf "(1) %o: torsion %o; census curves isomorphic over Q (among those above): %o\n", s[1], Tors(C), iso;
end for;
printf "(1) BFT2014 Ex. 21: the printed curve is C_{2,-1,-2} of Thm 6: %o\n", IsIsomorphic(HyperellipticCurve(Crst(2,-1,-2)), HyperellipticCurve(12*x^6+72*x^5+312*x^4+688*x^3+768*x^2+192*x+68));
printf "(1) 3.3.a is the fiber (r,s,t) = (5/2,16/9,25/16) of BFT2014 Thm 6: %o (Thm 6 parametrises the curves with a rational (Z/3)^2 of trivial Weil pairing; that fiber is not printed in the paper)\n", IsIsomorphic(HyperellipticCurve(Crst(5/2,16/9,25/16)), Cen("3.3.a"));

// (2) one-parameter families
Ft<t> := FunctionField(Q);
Pt<X> := PolynomialRing(Ft);
function WeightedEq(a, b)   // the polynomials in t whose common rational roots give a = b in P(1,2,3,5)
    w := [1, 2, 3, 5];
    return [ Numerator(a[i]^w[j]*b[j]^w[i] - a[j]^w[i]*b[i]^w[j]) : i in [1..3], j in [2..4] | i lt j ];
end function;
function FiberSearch(F, C)   // F in Q(t)[X], the family y^2 = F; C a curve over Q
    ICt := IgusaClebschInvariants(F);
    ICc := [ Ft!a : a in IgusaClebschInvariants(C) ];
    eqs := [ e : e in WeightedEq(ICt, ICc) | e ne 0 ];
    g := GCD(eqs);
    rts := [ r[1] : r in Roots(g) ];
    out := [];
    for r in rts do
        cf := [ Evaluate(c, r) : c in Coefficients(F) ];
        Fr := P!cf;
        if Discriminant(Fr) eq 0 or Degree(Fr) lt 5 then continue; end if;
        if IsIsomorphic(HyperellipticCurve(Fr), C) then Append(~out, r); end if;
    end for;
    return rts, out;
end function;
function Member(F, r) return HyperellipticCurve(P![ Evaluate(c, r) : c in Coefficients(F) ]); end function;
function Moved(C) return Transformation(SimplifiedModel(C), [2, 1, 1, -3]); end function;   // x -> (2x+1)/(x-3)
fam := [*
 <"Flynn 1991 Result 3 (g = 2, r = 1; order 10)", 4*((X^3+X^2)^2 - t*(X^4+X^2)) + t^2*(X-1)^2, "10.a", -5/3>,
 <"Flynn 1991 Result 3 (g = 2, r = 2; order 9)", 4*(X^6 - t*(X^4+X^3)) + t^2*(X-1)^2, "9.a", 7/2>,
 <"Leprevost 1997 C_{1,t} (g = 2; order 20)", (X*((X-1)^2-X^2) + t*((X-1)^2+X^2))^2 + 4*t*X^5, "20.a", 3/5>,
 <"Leprevost 1995 [3,9] family C_t", X^6+6*(t-1)*X^5+3*(t-1)*(3*t-5)*X^4-(18*t^2-34*t+18)*X^3+3*(t-1)*(5*t-3)*X^2-6*t*(t-1)*X+t^2, "3.9.a", 11/4> *];
for fm in fam do
    F := fm[2]; t0 := fm[4];
    _, ctl := FiberSearch(F, Moved(Member(F, t0)));
    assert t0 in ctl;
    printf "(2) %o: torsion of the member t = %o is %o; control: the member t = %o, moved, is found at %o\n", fm[1], t0, Tors(Member(F, t0)), t0, ctl;
    rts, iso := FiberSearch(F, Cen(fm[3]));
    printf "(2) %o vs %o: rational t with the same Igusa-Clebsch point %o; isomorphic over Q at %o\n", fm[1], fm[3], rts, iso;
end for;

// (3) Flynn 1991 Result 1(b) (g = 2, r = 1): every member has a rational Weierstrass point W and a rational
//     point P with [P - W] of order 6.  Test: for each rational Weierstrass point W and each D of order 6,
//     is D + [W - oo] of the form [P - oo] (Mumford degree 1)?  (Model y^2 = F of odd degree.)
function Pairs(C)
    C1 := SimplifiedModel(C);
    F := HyperellipticPolynomials(C1);
    if Degree(F) eq 6 then
        rr := Roots(F);
        if #rr eq 0 then return -1; end if;   // no rational Weierstrass point
        a := rr[1][1];
        F := P!(x^6*Evaluate(F, a + 1/x));   // x -> a + 1/x, y -> y/x^3: the root a goes to infinity
        C1 := HyperellipticCurve(F);
    end if;
    assert Degree(F) eq 5;
    d := LCM([Denominator(c) : c in Coefficients(F)]);
    F := d^2*F; C1 := HyperellipticCurve(F);   // y -> d y: an integral model, as TorsionSubgroup needs
    J := Jacobian(C1);
    T, m := TorsionSubgroup(J);
    six := [ m(g) : g in T | Order(g) eq 6 ];
    Ws := [ J!0 ] cat [ elt< J | Parent(F).1 - r[1], 0 > : r in Roots(F) ];
    return #[ 0 : D in six, W in Ws | Degree((D + W)[1]) eq 1 ];
end function;
Fl := func< w1, w0 | HyperellipticCurve(-x^3*(x^2 + w1*x + w0), -(x^2 + w1*x + w0)) >;   // Y^2 = (Y - X^3)(X^2 + w1 X + w0)
ctl := Pairs(Moved(Fl(1, 3)));
assert ctl gt 0;
printf "(3) control: the member (w1, w0) = (1, 3) of Flynn 1991 Result 1(b), moved, has torsion %o and %o pairs (W, D)\n", Tors(Fl(1, 3)), ctl;
printf "(3) 6.a: torsion %o; pairs (W, D) with D of order 6 and D + [W - oo] = [P - oo]: %o, so 6.a is not a member of Result 1(b)\n", Tors(Cen("6.a")), Pairs(Cen("6.a"));

// (4) 21.a has no rational Weierstrass point (Leprevost 1995, p. 290: the order-21 curves of Leprevost 1991 have one)
f21, h21 := HyperellipticPolynomials(Cen("21.a"));
F21 := 4*f21 + h21^2;
printf "(4) 21.a: 4f + h^2 = %o has degree %o and rational roots %o: no rational Weierstrass point\n", F21, Degree(F21), [r[1] : r in Roots(F21)];
assert Degree(F21) eq 6 and #Roots(F21) eq 0;
print "SEED_CREDIT_CHECK_DONE";
quit;
