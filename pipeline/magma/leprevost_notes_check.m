/* leprevost_notes_check.m -- are the seed curves credited to Leprevost's C. R. notes curves of those notes?
   (second review of 2026-10-20 in data/knowledge/sources/PROVENANCE.md)

   The notes (C. R. Acad. Sci. Paris 313 (1991), 451-454 and 771-774; 316 (1993), 819-821) could not be read.  By
   their zbMATH Open summaries the 1991 notes give one-parameter families with a point of order 13, and 15, 17, 19,
   21, over Q(t), and the 1993 note gives one curve C_l for l = 21, 22, 23, 25, 26, 27, 29, two for l = 24, and the
   [3,9] family: the contents of the full versions, Leprevost, Manuscripta Math. 75 (1992), 303-326, and
   J. Theor. Nombres Bordeaux 7 (1995), 283-306.  Equations (read from the page scans):
     Leprevost 1995, p. 286 (Thm 1.2.1): f22 = (2x^2-2x+1)(2x^4-2x^3+x^2-4x+4), f24,1 = (2x^2-2x-1)(2x^4-10x^3+7x^2+4x-4),
               f24,2 = (x^2-x+1)(x^4-3x^3+8x^2-3x+1), f26 = (6x^2-6x+1)(6x^4-30x^3+49x^2-20x+4).
     Leprevost 1992, p. 316 (l = 13): y^2 = [(t+1)x^2 - 2x + 1]^2 - 4t x^3 (x-1)^2 ("celle obtenue en [3]", the 1991 note);
               p. 317 (l = 15): y^2 = [x^3 - (x-1)^2 (x-1-t)]^2 - 4t x^3 (x-1)^2, D = -4096 t^4 (t+1)^3 (16t^2 - 216t - 243);
               p. 319 (l = 17): y^2 = A^2 - 4t^2 (t-1)^3 x^3 (x-1)^2, A = (t^3+2t^2-3t+1)x^2 - t(2t^2+t-1)x + t^3;
               p. 322 (l = 19): y^2 = [a2 x^2 + a1 x + a0]^2 - lam x^3 (x-1)^2 with
                 a2 = -(t^10+6t^9+29t^8+168t^7+482t^6-348t^5+482t^4+168t^3+29t^2+6t+1),
                 a1 = 2(t+1)^2 (t^8+4t^7+20t^6+92t^5+22t^4+92t^3+20t^2+4t+1), a0 = -(t+1)^6 (t^4+14t^2+1),
                 lam = 1024 t^4 (t^2-1)^2 (t^4+14t^2+1)^2.
   (1) Each printed curve has a point of the claimed order; it is not isomorphic over Qbar to the census curve.
   (2) Each family has a point of the claimed order at two members; a family search (as in seed_credit_check.m: the
       t with the Igusa-Clebsch point of the census curve, then isomorphism over Q), first run on a known member moved
       by x -> (2x+1)/(x-3), finds no fiber for 15.a and 19.a (not even over Qbar: the common factor is 1, resp. supported
       on the singular members t = 0, +-1), and the fibers t = +-1 (13.a) and t = -1 (17.a).
   Run:  magma -b pipeline/magma/leprevost_notes_check.m < /dev/null > data/knowledge/sources/leprevost_notes_check.log */
SetColumns(0); SetMemoryLimit(16*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q);
TorsOrd := func< C | #TorsionSubgroup(Jacobian(SimplifiedModel(MinimalWeierstrassModel(C)))) >;
Cen := AssociativeArray();
Cen["13.a"] := HyperellipticCurve(P![0,0,-1,-1], P![1,1,1,1]);
Cen["15.a"] := HyperellipticCurve(P![0,-1,-1], P![1,1,1,1]);
Cen["17.a"] := HyperellipticCurve(P![1,3,1,1], P![0,1,1,1]);
Cen["19.a"] := HyperellipticCurve(P![-2,5,-3,-5,1,1,1], P![1,1,1]);
Cen["22.a"] := HyperellipticCurve(P![1,-1,-2,1], P![0,1,0,1]);
Cen["24.a"] := HyperellipticCurve(P![0,2,4,3,2], P![1,0,0,1]);
Cen["26.a"] := HyperellipticCurve(P![0,0,-3,-3,-3,-12], P![1,1]);

// (1) Leprevost 1995, p. 286
pr := [* <"C22", (2*x^2-2*x+1)*(2*x^4-2*x^3+x^2-4*x+4), 22, "22.a">,
         <"C24,1", (2*x^2-2*x-1)*(2*x^4-10*x^3+7*x^2+4*x-4), 24, "24.a">,
         <"C24,2", (x^2-x+1)*(x^4-3*x^3+8*x^2-3*x+1), 24, "24.a">,
         <"C26", (6*x^2-6*x+1)*(6*x^4-30*x^3+49*x^2-20*x+4), 26, "26.a"> *];
for c in pr do
    C := HyperellipticCurve(c[2]);
    n := TorsOrd(C);
    assert n mod c[3] eq 0;
    g1 := G2Invariants(C); g2 := G2Invariants(Cen[c[4]]);
    assert g1 ne g2;
    printf "(1) Leprevost 1995 %o: #J(Q)_tors = %o (a point of order %o); G2-invariants %o, of %o: %o, so not isomorphic over Qbar\n", c[1], n, c[3], g1, c[4], g2;
end for;

// (2) Leprevost 1992 families
Ft<t> := FunctionField(Q); Pt<X> := PolynomialRing(Ft);
function WeightedEq(a, b)
    w := [1, 2, 3, 5];
    return [ Numerator(a[i]^w[j]*b[j]^w[i] - a[j]^w[i]*b[i]^w[j]) : i in [1..3], j in [2..4] | i lt j ];
end function;
function FiberSearch(F, C)
    ICt := IgusaClebschInvariants(F);
    ICc := [ Ft!a : a in IgusaClebschInvariants(C) ];
    eqs := [ e : e in WeightedEq(ICt, ICc) | e ne 0 ];
    g := GCD(eqs);
    rts := [ r[1] : r in Roots(g) ];
    out := [];
    for r in rts do
        Fr := P![ Evaluate(c, r) : c in Coefficients(F) ];
        if Discriminant(Fr) eq 0 or Degree(Fr) lt 5 then continue; end if;
        if IsIsomorphic(HyperellipticCurve(Fr), C) then Append(~out, r); end if;
    end for;
    return rts, out, g;
end function;
function Member(F, r) return HyperellipticCurve(P![ Evaluate(c, r) : c in Coefficients(F) ]); end function;
function Moved(C) return Transformation(SimplifiedModel(C), [2, 1, 1, -3]); end function;   // x -> (2x+1)/(x-3)
a2 := -(t^10+6*t^9+29*t^8+168*t^7+482*t^6-348*t^5+482*t^4+168*t^3+29*t^2+6*t+1);
a1 := 2*(t+1)^2*(t^8+4*t^7+20*t^6+92*t^5+22*t^4+92*t^3+20*t^2+4*t+1);
a0 := -(t+1)^6*(t^4+14*t^2+1);
F15 := (X^3 - (X-1)^2*(X-1-t))^2 - 4*t*X^3*(X-1)^2;
assert LeadingCoefficient(F15)^2*Discriminant(F15) eq -4096*t^4*(t+1)^3*(16*t^2-216*t-243);
print "(2) l = 15: the discriminant is the printed D = -4096 t^4 (t+1)^3 (16t^2 - 216t - 243)";
fam := [* <"l = 13 (p. 316)", ((t+1)*X^2 - 2*X + 1)^2 - 4*t*X^3*(X-1)^2, 13, "13.a", [Q| 2, 3], 5/2, [Q| -1, 1]>,
          <"l = 15 (p. 317)", F15, 15, "15.a", [Q| 2, 3], 2/9, [Q| ]>,
          <"l = 17 (p. 319)", ((t^3+2*t^2-3*t+1)*X^2 - t*(2*t^2+t-1)*X + t^3)^2 - 4*t^2*(t-1)^3*X^3*(X-1)^2, 17, "17.a", [Q| 2, 3], 4/7, [Q| -1]>,
          <"l = 19 (p. 322)", (a2*X^2 + a1*X + a0)^2 - 1024*t^4*(t^2-1)^2*(t^4+14*t^2+1)^2*X^3*(X-1)^2, 19, "19.a", [Q| 2, 3], 5/3, [Q| ]> *];
for fm in fam do
    F := fm[2];
    for t0 in fm[5] do
        n := TorsOrd(Member(F, t0));
        assert n mod fm[3] eq 0;
        printf "(2) Leprevost 1992 %o: the member t = %o has #J(Q)_tors = %o (a point of order %o)\n", fm[1], t0, n, fm[3];
    end for;
    _, ctl := FiberSearch(F, Moved(Member(F, fm[6])));
    assert fm[6] in ctl;
    rts, iso, dg := FiberSearch(F, Cen[fm[4]]);
    assert Set(iso) eq Set(fm[7]);
    printf "(2) Leprevost 1992 %o: control (the member t = %o, moved) found at %o; %o: the Igusa-Clebsch equations have the common factor %o (rational roots %o), isomorphic over Q at %o\n", fm[1], fm[6], ctl, fm[4], Degree(dg) eq 0 select "1" else Sprint([<Sprint(Evaluate(a[1], t)), a[2]> : a in Factorization(dg)]), rts, iso;
end for;
print "LEPREVOST_NOTES_CHECK_DONE";
quit;
