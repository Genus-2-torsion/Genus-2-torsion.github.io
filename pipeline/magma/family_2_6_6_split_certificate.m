/* family_2_6_6_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves over Q
   have Jacobian split over Q (isogenous over Q to E1 x E2) with J(Q)_tors exactly [2,6,6].
   Found 2026-10-19 (torsion_jac notes/hlpES_2026_10_19.md).

   Source.  Howe-Leprevost-Poonen 2000, Section 3.7: glue (HLP Prop. 4) E1 = E_{2,6}(t) and E2 = E_{2,6}(u) (Tate
   normal form, HLP Table 3, P = (0,0) of order 6) along psi = identity on the 2-torsion points T1, T2, T3 labelled as
   in HLP Table 5.  The image of E1(Q) x E2(Q) is [6,6]; a half of (T1, T1') gives a new rational point of order 2
   (HLP Prop. 12) when 2(t - 3) = 2(u - 3) y^2 and (t + 3)(t - 5) = (u + 3)(u - 5) z^2 (HLP system (6)).  With
   t = 3 + (u - 3) y^2 this is the genus-1 curve ((u-3)y^2 + 6)((u-3)y^2 - 2) = (u+3)(u-5) z^2 over Q(u), with origin
   (y, z) = (1, 1) and HLP's point sigma = (-1, 1) of infinite order.  Here we take 2 sigma (computed below on a
   Weierstrass model of the generic fibre); it gives y(u) of degree 4 and t(u) of degree 9.  (The points y = +-1
   give t = u, the same curve, and are not used.)  C_u : y^2 = f_u(x) is the Q(u)-isomorphic model
   f(x) = hHLP(lam x)/mu^2 of the HLP Prop. 4 sextic; f = F(x^2) with F a cubic over Q(u).

   Argument.
   (1) Splitting.  phi1 : (x,y) -> (x^2, y) maps C_u onto Ea : Y^2 = F(X), phi2 : (x,y) -> (1/x^2, y/x^3) onto
       Eb : Y^2 = X^3 F(1/X); Ea is Q(u)-isomorphic to E2 and Eb to E1, so J is split over Q at every smooth member.
   (2) The marked group over Q(u).  f has four roots in Q(u) (and an irreducible quadratic factor), so the eight
       classes [(x - r_i)(x - r_j), 0], [quadratic, 0] and 0 form (Z/2)^3 in J(Q(u)); D3 = 2 phi2^*(P1) and
       D3' = 2 phi1^*(P2) have order 3 and D3' is not in <D3>.  So (Z/2)^3 x (Z/3)^2 = [2,6,6] is marked over Q(u).
   (3) Exactness on a residue class.  #J(F_p) is always divisible by 144, so the group structure is used.  Every
       11- and 17-integral u = 4 (mod 187) has the reductions of u = 4 at 11 and 17, with J(F_11) = [2,2,6,6] and
       J(F_17) = [2,6,24]; the largest group embedding in both is [2,6,6], so J_u(Q)_tors is exactly [2,6,6].
   (4) Infinitely many: the G2-invariants of the members differ.   (5) Sanity: TorsionSubgroup at u = 4 and 191.
   Run (from the repository root):
       magma -b pipeline/magma/family_2_6_6_split_certificate.m > data/knowledge/sources/family_2_6_6_split_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
SetQuitOnError(true);                          // fail closed: any failed assertion ends the run before the marker
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

// ---------------------------------------------------------------- helpers (shared by the hlpES certificates)
// HLP 2000, Table 3: Tate normal form y^2 + (1-c)xy - by = x^3 - bx^2, (0,0) of maximal order.
function KubertBC(N, t)
    case N:
        when 8: return 2*t^2 - 3*t + 1, (2*t^2 - 3*t + 1)/t;
        when 10: return (2*t^5 - 3*t^4 + t^3)/(t^2 - 3*t + 1)^2, (-2*t^3 + 3*t^2 - t)/(t^2 - 3*t + 1);
        when 12: return (12*t^6 - 30*t^5 + 34*t^4 - 21*t^3 + 7*t^2 - t)/(t - 1)^4, (-6*t^4 + 9*t^3 - 5*t^2 + t)/(t - 1)^3;
        when 24: return t^2 - 1/16, 0*t;                                                         // (2,4)
        when 26: return (-2*t^3 + 14*t^2 - 22*t + 10)/((t + 3)^2*(t - 3)^2), (-2*t + 10)/((t + 3)*(t - 3));   // (2,6)
        when 28: return (16*t^3 + 16*t^2 + 6*t + 1)/(8*t^2 - 1)^2, (16*t^3 + 16*t^2 + 6*t + 1)/(2*t*(4*t + 1)*(8*t^2 - 1)); // (2,8)
    end case;
    error "unknown level";
end function;
function TateCurve(b, c) return EllipticCurve([1 - c, -b, -b, 0, 0]); end function;
// y^2 = MonicCubic(E)(x) after completing the square in y (x unchanged)
function MonicCubic(E, P)
    a1, a2, a3, a4, a6 := Explode(aInvariants(E)); xx := P.1;
    return xx^3 + (a2 + a1^2/4)*xx^2 + (a4 + a1*a3/2)*xx + (a6 + a3^2/4);
end function;
// HLP 2000, Prop. 4: the sextic glued from labelled roots al (curve E) and be (curve F), psi((al_i,0)) = (be_i,0).
function HLPSextic(al, be)
    K := Parent(al[1]); PK<X> := PolynomialRing(K);
    a1 := (al[3]-al[2])^2/(be[3]-be[2]) + (al[2]-al[1])^2/(be[2]-be[1]) + (al[1]-al[3])^2/(be[1]-be[3]);
    b1 := (be[3]-be[2])^2/(al[3]-al[2]) + (be[2]-be[1])^2/(al[2]-al[1]) + (be[1]-be[3])^2/(al[1]-al[3]);
    a2 := al[1]*(be[3]-be[2]) + al[2]*(be[1]-be[3]) + al[3]*(be[2]-be[1]);
    b2 := be[1]*(al[3]-al[2]) + be[2]*(al[1]-al[3]) + be[3]*(al[2]-al[1]);
    assert a1 ne 0 and a2 ne 0 and b1 ne 0 and b2 ne 0;          // psi does not come from an isomorphism
    Df := ((al[1]-al[2])*(al[1]-al[3])*(al[2]-al[3]))^2;
    Dg := ((be[1]-be[2])*(be[1]-be[3])*(be[2]-be[3]))^2;
    A := Dg*a1/a2; B := Df*b1/b2;
    return -(A*(al[2]-al[1])*(al[1]-al[3])*X^2 + B*(be[2]-be[1])*(be[1]-be[3]))
           *(A*(al[3]-al[2])*(al[2]-al[1])*X^2 + B*(be[3]-be[2])*(be[2]-be[1]))
           *(A*(al[1]-al[3])*(al[3]-al[2])*X^2 + B*(be[1]-be[3])*(be[3]-be[2]));
end function;
// A smaller model: g(x) = f(lam x)/mu^2 with lam, mu in Q(t)* (so y^2 = g is Q(t)-isomorphic to y^2 = f), chosen
// prime by prime to minimise the valuations of the coefficients of the even sextic f.  Returns g, lam, mu.
function MinimizeEvenSextic(f)
    P := Parent(f); Ft := BaseRing(P); R := PolynomialRing(Q); RZ := PolynomialRing(Z);
    D := LCM([ Denominator(c) : c in Coefficients(f) ]);
    ff := f*D^2;
    a := [ R!Numerator(Coefficient(ff, 2*i)) : i in [0..3] ];
    assert &and[ Denominator(Coefficient(ff, 2*i)) eq 1 : i in [0..3] ];
    den := LCM([Z| Denominator(e) : e in Coefficients(c), c in a ]);
    a := [ RZ!(den^2*c) : c in a ];
    lam := Ft!1; mu := Ft!(1/(D*den));
    irr := &join[ {@ p[1] : p in Factorization(c) @} : c in a | c ne 0 ];
    ints := &join[ { p[1] : p in TrialDivision(Abs(Content(c)), 10000) | p[1] lt 10000 } : c in a | c ne 0 ];
    pis := [ RZ!p : p in irr | Degree(p) gt 0 ] cat [ RZ!p : p in Sort(Setseq(ints)) ];
    for pi in pis do
        v := [ a[i] eq 0 select 10^6 else Valuation(a[i], pi) : i in [1..4] ];
        best := <&+[ x : x in v | x lt 10^6 ], 0, 0>;
        for k in [-30..30] do
            w := [ v[i] + 2*(i-1)*k : i in [1..4] ];
            m := Min([ w[i] : i in [1..4] | a[i] ne 0 ]) div 2;
            tot := &+[ w[i] - 2*m : i in [1..4] | a[i] ne 0 ];
            if tot lt best[1] then best := <tot, k, m>; end if;
        end for;
        k := best[2]; m := best[3];
        if k ne 0 or m ne 0 then
            for i in [1..4] do
                if a[i] ne 0 then
                    e := 2*(i-1)*k - 2*m;
                    if e ge 0 then a[i] := a[i]*pi^e; else a[i] := ExactQuotient(a[i], pi^(-e)); end if;
                end if;
            end for;
            lam := lam*(Ft!R!pi)^k; mu := mu*(Ft!R!pi)^m;
        end if;
    end for;
    g := P![ Ft!R!a[1], 0, Ft!R!a[2], 0, Ft!R!a[3], 0, Ft!R!a[4] ];
    assert Evaluate(f, lam*P.1) eq mu^2*g;
    return g, lam, mu;
end function;
// the two elliptic quotients of y^2 = f(x) = G(x^2): phi1 = (x^2, y) onto Ea : Y^2 = G(X), phi2 = (1/x^2, y/x^3)
// onto Eb : Y^2 = X^3 G(1/X); Weierstrass models with X' = cX, Y' = cY (c the leading coefficient)
function MonEC(G)
    cs := [ Coefficient(G, i) : i in [0..3] ];
    return EllipticCurve([0, cs[3], 0, cs[2]*cs[4], cs[1]*cs[4]^2]), cs[4];
end function;
function Quotients(f)
    P := Parent(f);
    assert &and[ Coefficient(f, i) eq 0 : i in [1, 3, 5] ];
    G := P![ Coefficient(f, 2*i) : i in [0..3] ];
    Gr := P![ Coefficient(G, 3 - i) : i in [0..3] ];
    Ea, ca := MonEC(G); Eb, cb := MonEC(Gr);
    return Ea, ca, Eb, cb;
end function;
// pullbacks to J of points of Ea (through phi1) and Eb (through phi2)
function Pull1(JJ, Pt, c) PP := Parent((JJ!0)[1]); xx := PP.1; XX := Pt[1]/c; YY := Pt[2]/c; return elt< JJ | xx^2 - XX, PP!YY >; end function;
function Pull2(JJ, Pt, c) PP := Parent((JJ!0)[1]); xx := PP.1; XX := Pt[1]/c; YY := Pt[2]/c; return elt< JJ | xx^2 - 1/XX, (YY/XX)*xx >; end function;
// order of a torsion class (searching n <= nmax; 0 if not found)
function Ord(D, nmax)
    Z0 := Parent(D)!0; E := D;
    for n in [1..nmax] do
        if E eq Z0 then return n; end if;
        E := E + D;
    end for;
    return 0;
end function;
// the largest abelian group embedding in every group with the given invariant factors (as in family_18_split)
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
// specialize a class of J(Q(t)) at the rational t0 (curve y^2 = SpecPoly(f, t0))
function SpecClass(D, t0, Jt0)
    PQ := Parent((Jt0!0)[1]);
    ev := func< c | Evaluate(c, t0) >;
    return elt< Jt0 | PQ![ ev(c) : c in Coefficients(D[1]) ], PQ![ ev(c) : c in Coefficients(D[2]) ], D[3] >;
end function;
// the subgroup generated by gens (brute-force closure) and its invariant factors
function GeneratedGroup(gens)
    S := {@ Parent(gens[1])!0 @};
    for g in gens do
        h := g; newS := S;
        while not h in S do newS := newS join {@ s + h : s in S @}; h := h + g; end while;
        S := newS;
    end for;
    n := #S; Z0 := Parent(gens[1])!0; parts := [Integers()| ];
    for l in PrimeDivisors(n) do
        c := [1]; k := 0;
        repeat k +:= 1; Append(~c, #[ s : s in S | l^k*s eq Z0 ]); until c[k+1] eq l^Valuation(n, l);
        m := [ Valuation(c[i+1] div c[i], l) : i in [1..k] ];      // m[i] = number of cyclic factors of order >= l^i
        for j in [1..m[1]] do Append(~parts, l^#[ i : i in [1..k] | m[i] ge j ]); end for;
    end for;
    return n, Invariants(AbelianGroup(parts));
end function;
function SeqStr(S) return "[ " cat &cat[ Sprint(S[i]) cat (i lt #S select ", " else "") : i in [1..#S] ] cat " ]"; end function;
// the rational root r and the monic quadratic x^2 + p x + q of a cubic with exactly one rational root
function Split2(g)
    fa := Factorization(g);
    lin := [ h[1] : h in fa | Degree(h[1]) eq 1 ]; qu := [ h[1] : h in fa | Degree(h[1]) eq 2 ];
    assert #lin eq 1 and #qu eq 1;
    return -Coefficient(lin[1], 0)/Coefficient(lin[1], 1), Coefficient(qu[1], 1)/Coefficient(qu[1], 2), Coefficient(qu[1], 0)/Coefficient(qu[1], 2);
end function;
// HLP Prop. 4 for E1, E2 over Ft = BaseRing(P), each with one rational point of order 2 and the same quadratic field of
// 2-torsion: evaluated over K = Ft(w), w^2 = delta1; psi matches the rational roots and sends (-p1 + w)/2 to (-p2 + m w)/2,
// m^2 = delta2/delta1.  Returns the sextic, which has coefficients in Ft.
function GlueOnePoint(E1, E2, P)
    Ft := BaseRing(P);
    r1, p1, q1 := Split2(MonicCubic(E1, P)); r2, p2, q2 := Split2(MonicCubic(E2, P));
    del1 := p1^2 - 4*q1; del2 := p2^2 - 4*q2;
    okm, mm := IsSquare(del2/del1);
    assert okm;                                  // the 2-torsion of E1 and E2 are isomorphic Galois modules
    K<w> := ext< Ft | PolynomialRing(Ft).1^2 - del1 >;
    al := [ K!r1, (-p1 + w)/2, (-p1 - w)/2 ];
    be := [ K!r2, (-p2 + mm*w)/2, (-p2 - mm*w)/2 ];
    hK := HLPSextic(al, be);
    return P![ Ft!c : c in Coefficients(hK) ];
end function;
// HLP Prop. 4 for labelled rational roots (full rational 2-torsion)
function GlueFull(al, be, P)
    Ft := BaseRing(P);
    return P![ Ft!c : c in Coefficients(HLPSextic(al, be)) ];
end function;

// ---------------------------------------------------------------- (0) source: the section 2 sigma of HLP's system (6)
Fu<u> := FunctionField(Q); P<x> := PolynomialRing(Fu);
Py<yy> := PolynomialRing(Fu);
kk := (u + 3)*(u - 5);
Cy := HyperellipticCurve(kk*((u - 3)*yy^2 + 6)*((u - 3)*yy^2 - 2));   // Z^2 = kk^2 z^2 ... with Z = kk z
EW, phiW := EllipticCurve(Cy, Cy![1, kk, 1]);              // origin (y, z) = (1, 1)
sig := phiW(Cy![-1, kk, 1]);                               // sigma = (-1, 1)
R2 := Inverse(phiW)(2*sig);
assert R2[3] ne 0;
yv := R2[1]/R2[3];
assert yv eq (-u^4 + 12*u^3 - 126*u^2 + 348*u + 279)/(3*u^4 - 20*u^3 - 54*u^2 + 540*u - 981);
t := 3 + (u - 3)*yv^2;
oksq := IsSquare((t + 3)*(t - 5)*(u + 3)*(u - 5)); assert oksq;
printf "(0) on HLP's system (6): the section 2 sigma has y = %o, t = 3 + (u-3) y^2 =\n    %o;\n    (t-3)/(u-3) = y^2 and (t+3)(t-5)(u+3)(u-5) is a square\n", yv, t;
b1, c1 := KubertBC(26, t); E1 := TateCurve(b1, c1);
b2, c2 := KubertBC(26, u); E2 := TateCurve(b2, c2);
P1 := E1![0,0]; P2 := E2![0,0];
assert Ord(P1, 6) eq 6 and Ord(P2, 6) eq 6;
assert jInvariant(E1) ne jInvariant(E2);
function Roots26(z)       // HLP Table 5: x(T1), x(T2), x(T3) on E_{2,6}(z) (model with the square completed)
    return [ (-2*z+10)/((z+3)*(z-3)), (-z^3+7*z^2-11*z+5)/(4*(z+3)*(z-3)^2), (-2*z^2+4*z-2)/((z+3)^2*(z-3)) ];
end function;
al := Roots26(t); be := Roots26(u);
assert &and[ Evaluate(MonicCubic(E1, P), r) eq 0 : r in al ] and &and[ Evaluate(MonicCubic(E2, P), r) eq 0 : r in be ];
h0 := GlueFull(al, be, P);
f, lam, mu := MinimizeEvenSextic(h0);
assert Degree(f) eq 6 and IsSquarefree(f);
printf "(0) E1 = E_{2,6}(t), E2 = E_{2,6}(u) (different j), psi = identity on the Table 5 labels; HLP Prop. 4 gluing hHLP over Q(u);\n    C_u : y^2 = f_u(x) = hHLP(lam x)/mu^2 with\n    lam = %o\n    mu = %o\n", lam, mu;
printf "    f_u = %o\n", f;

// ---------------------------------------------------------------- (1) splitting
Ea, ca, Eb, cb := Quotients(f);
oka, ma := IsIsomorphic(E2, Ea); okb, mb := IsIsomorphic(E1, Eb);
assert oka and okb;
print "(1) Ea : Y^2 = F(X) is Q(u)-isomorphic to E2 = E_{2,6}(u), Eb : Y^2 = X^3 F(1/X) to E1 = E_{2,6}(t): J is split over Q";

// ---------------------------------------------------------------- (2) the marked group over Q(u)
J := Jacobian(HyperellipticCurve(f));
fa := Factorization(f);
lins := [ g[1] : g in fa | Degree(g[1]) eq 1 ]; quads := [ g[1] : g in fa | Degree(g[1]) eq 2 ];
assert #lins eq 4 and #quads eq 1;
lc := LeadingCoefficient(f);
two := [ J!0 ] cat [ elt< J | lins[i]*lins[j], 0 > : i, j in [1..4] | i lt j ] cat [ elt< J | quads[1], 0 > ];
assert #Seqset(two) eq 8 and &and[ 2*c eq J!0 : c in two ];
D3 := 2*Pull2(J, mb(P1), cb); D3p := 2*Pull1(J, ma(P2), ca);
assert Ord(D3, 3) eq 3 and Ord(D3p, 3) eq 3 and not (D3p in [J!0, D3, 2*D3]);
print "(2) over Q(u): f has four roots in Q(u); the eight classes from them are distinct and killed by 2 ((Z/2)^3);";
print "    D3 = 2 phi2^*(P1), D3' = 2 phi1^*(P2) have order 3 and D3' is not in <D3>: [2,6,6] = (Z/2)^3 x (Z/3)^2 is marked over Q(u)";

// ---------------------------------------------------------------- (3) exactness on the class u = 4 mod 187
primes := [11, 17]; res := [4, 4];
members := [Q| 4, 191, -183, 378, 565];
invs := [];
for i in [1..#primes] do
    p := primes[i]; r0 := res[i];
    fbar := ReductionClass(f, p, r0);
    for m in members do
        assert Valuation(m, p) ge 0 and GF(p)!m eq GF(p)!r0;
        assert PolynomialRing(GF(p))!SpecPoly(f, m) eq fbar;
    end for;
    n, chi := FrobData(fbar);
    inv := Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
    assert &*inv eq n;
    Append(~invs, inv);
    printf "(3) p = %o, u = %o mod %o: reduction %o (squarefree, degree %o), J(F_%o) = %o, chi = %o\n", p, r0, p, fbar, Degree(fbar), p, inv, chi;
end for;
bnd := StructureBound(invs);
printf "(3) the largest group embedding in all these J(F_p) is %o\n", bnd;
assert bnd eq [2, 6, 6];
print "(3) every 11- and 17-integral u = 4 mod 187 has J_u(Q)_tors = [2,6,6] exactly";

// ---------------------------------------------------------------- (4) infinitely many, (5) sanity
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } eq #members;
printf "(4) the G2-invariants of the members %o are pairwise distinct\n", members;
for m in [Q| 4, 191] do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(5) TorsionSubgroup at u = %o: %o\n", m, invm;
    assert invm eq [2, 6, 6];
end for;
printf "time %o s\n", Cputime(T0);
print "FAMILY_2_6_6_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E1 x E2) and J(Q)_tors exactly [2,6,6]";
quit;
