/* split_11_remark.m -- why the construction of the split [11] curve 11.b gives no family
   (remark for the open grade of split [11]; not a certificate).

   11.b (BNSS2026) is the (2,2)-gluing of an elliptic curve E over K = Q(sqrt 11) with a point of order 11
   to its Galois conjugate E^sigma.  Model (BNSS2026): X1(11) : r^2 - r A(s) + s = 0, A = s^3 - 3s^2 + 4s, a
   double cover of the s-line; E_t = E(b, c) in Tate normal form, b = r s (r - 1), c = s (r - 1), (0,0) of
   order 11.  For s in Q with r quadratic, K_s = Q(r) and sigma : r -> A(s) - r.  Every elliptic curve over a
   quadratic field with a point of order 11 is such an E_t (up to isomorphism): the quadratic points of the
   rank-0 elliptic curve X1(11) lie on the s-line up to a diamond operator, which does not change E.

   The gluing needs a G_K-isomorphism psi : E[2] -> E^sigma[2] with psi^sigma psi = id.  Writing psi as the
   Moebius map x -> mu + nu/(x - mu^sigma) between the (shifted) 2-division cubics, mu = m1 + m2 r with
   m1, m2, nu in Q, the conditions are (i) N(c0) = nu^3 and (ii) nu b = c0 a^sigma, where
   g(u) = f(u + mu) = u^3 + a u^2 + b u + c0.  Over Q(s) these have a 6-point solution scheme (one point for
   each bijection E[2] -> E^sigma[2]); eliminating gives an irreducible G6(nu) of degree 6, i.e. a curve Y over
   Q with a degree-6 map to the s-line whose rational points are the inputs of the construction.

   (1) Y over Q(s): G6 has degree 6 and is irreducible; H(S, V) is G6 with denominators cleared.
   (2) Modulo p = 101, 103, 107 the reduction of H is irreducible of degree 6 over F_p(s), with exact constant
       field F_p, and Y mod p has genus 8.  The geometric genus does not increase under reduction, so Y has
       genus >= 8 and Y(Q) is finite (Faltings): the construction gives finitely many curves.
   (3) s = 4/5 (BNSS's point; K = Q(sqrt 11)) has a rational nu; s = 2 has none (irreducible fibre).
   (4) Search: s = a/b with |a| <= HB, 1 <= b <= HB having a rational root nu of H(s, nu).
   (5) Cross-check of (4) on a smaller box without H: a G_K-isomorphism E[2] -> E^sigma[2] exists only if the
       conjugate 2-division cubic f^sigma has a root in K[x]/(f).
   Run (from the repository root):
       magma -b pipeline/magma/split_11_remark.m > data/knowledge/sources/split_11_remark.log */
SetColumns(0); SetMemoryLimit(16*10^9);
if not assigned HB then HB := 400; elif Type(HB) eq MonStgElt then HB := StringToInteger(HB); end if;
if not assigned HB2 then HB2 := 30; elif Type(HB2) eq MonStgElt then HB2 := StringToInteger(HB2); end if;
Q := Rationals();

// (1) the curve Y over Q(s)
FF<s> := FunctionField(Q);
A := s^3 - 3*s^2 + 4*s;
R<m1, m2, nu> := PolynomialRing(FF, 3);
// elements of Q(s)(r) as pairs [x0, x1] = x0 + x1 r, r^2 = A r - s; sigma: r -> A - r
function Mul(x, y) return [x[1]*y[1] - s*x[2]*y[2], x[1]*y[2] + x[2]*y[1] + A*x[2]*y[2]]; end function;
function Add(x, y) return [x[1]+y[1], x[2]+y[2]]; end function;
function Sc(c, x) return [c*x[1], c*x[2]]; end function;
function Sig(x) return [x[1] + A*x[2], -x[2]]; end function;
rr := [R|0, 1]; one := [R|1, 0];
bb := Mul(Mul(rr, Sc(s, one)), Add(rr, Sc(-1, one)));   // b = r s (r - 1)
cc := Mul(Sc(s, one), Add(rr, Sc(-1, one)));            // c = s (r - 1)
a1 := Add(one, Sc(-1, cc)); a3 := Sc(-1, bb); a2 := Sc(-1, bb);
// monic 2-division cubic of E(b, c): u^3 + (a2 + a1^2/4) u^2 + (a1 a3/2) u + a3^2/4
fA := Add(a2, Sc(1/4, Mul(a1, a1))); fB := Sc(1/2, Mul(a1, a3)); fC := Sc(1/4, Mul(a3, a3));
mu := [m1, m2]; mu2 := Mul(mu, mu); mu3 := Mul(mu2, mu);
ga := Add(fA, Sc(3, mu));                                // g(u) = f(u + mu) = u^3 + ga u^2 + gb u + gc
gb := Add(Add(fB, Sc(2, Mul(fA, mu))), Sc(3, mu2));
gc := Add(Add(Add(fC, Mul(fB, mu)), Mul(fA, mu2)), mu3);
e1 := Add(Mul(gc, Sig(gc)), [-nu^3, R!0]);               // (i)
e2 := Add(Sc(nu, gb), Sc(-1, Mul(gc, Sig(ga))));         // (ii), two components
assert e1[2] eq 0;
I := Saturation(ideal<R | e1[1], e2[1], e2[2]>, nu);
assert Dimension(I) eq 0 and QuotientDimension(I) eq 6;
G := Basis(EliminationIdeal(I, {nu}))[1];
PV<VV> := PolynomialRing(FF);
G6 := PV![MonomialCoefficient(G, nu^i) : i in [0..Degree(G, nu)]];
G6 := G6/LeadingCoefficient(G6);
assert Degree(G6) eq 6 and IsIrreducible(G6);
Pq<S, V> := PolynomialRing(Q, 2);
den := LCM([Denominator(c) : c in Coefficients(G6)]);
H := &+[Pq | Evaluate(Numerator(c*den), S)*V^(i-1) : i -> c in Coefficients(G6)];
H := H div GCD([Pq | Evaluate(Numerator(c*den), S) : c in Coefficients(G6)]);
printf "(1) Y over Q(s): solution scheme of degree 6, elimination polynomial G6(nu) of degree 6, irreducible over Q(s); H(S, V) has degree %o in S and %o in V\n", Degree(H, S), Degree(H, V);

// (2) genus of the reductions
for p in [101, 103, 107] do
    Fp := GF(p); Fs<t> := FunctionField(Fp); Pp<w> := PolynomialRing(Fs);
    Hp := &+[Pp | (Fp!(MonomialCoefficient(H, S^i*V^j)*LCM([Denominator(c) : c in Coefficients(H)])))*t^i*w^j : i in [0..Degree(H, S)], j in [0..6]];
    assert Degree(Hp) eq 6 and IsIrreducible(Hp);
    Fy := FunctionField(Hp/LeadingCoefficient(Hp));
    assert Degree(ExactConstantField(Fy)) eq 1;
    gp := Genus(Fy);
    assert gp eq 8;
    printf "(2) p = %o: H mod p irreducible of degree 6 over F_p(s), exact constant field F_%o, genus %o\n", p, p, gp;
end for;
print "(2) so Y has geometric genus >= 8 and Y(Q) is finite";

// (3) the fibres at s = 4/5 and s = 2
PQ<v> := PolynomialRing(Q);
cfs := [ [MonomialCoefficient(H, S^i*V^j) : i in [0..Degree(H, S)]] : j in [0..6] ];
Pt<tt> := PolynomialRing(Q);
cfp := [ Pt!c : c in cfs ];
FibreH := func< s0 | PQ![Evaluate(c, s0) : c in cfp] >;
h45 := FibreH(4/5);
assert [r[1] : r in Roots(h45)] eq [-32768/5166015625];
printf "(3) s = 4/5: fibre factors with degrees %o, rational nu = %o (K = Q(sqrt 11): the construction of 11.b)\n", [Degree(f[1]) : f in Factorization(h45)], -32768/5166015625;
h2 := FibreH(2);
assert Degree(h2) eq 6 and IsIrreducible(h2);
print "(3) s = 2: fibre irreducible of degree 6 (E is not 2-congruent to its conjugate over K)";

// (4) search on Y
hits := [];
for b in [1..HB] do for a in [-HB..HB] do
    if GCD(a, b) ne 1 then continue; end if;
    s0 := a/b;
    h := FibreH(s0);
    if h eq 0 then printf "(4) ZERO fibre at s = %o\n", s0; continue; end if;
    rts := Roots(h);
    if #rts gt 0 then
        g := Pt![s0, -(s0^3 - 3*s0^2 + 4*s0), 1];
        printf "(4) s = %o: rational nu %o, r quadratic over Q: %o\n", s0, [x[1] : x in rts], IsIrreducible(g);
        if IsIrreducible(g) then Append(~hits, s0); end if;
    end if;
end for; end for;
assert hits eq [4/5];
printf "(4) s = a/b with |a| <= %o, 1 <= b <= %o: the only point of Y with r quadratic is s = 4/5 (the others are cusps, r rational)\n", HB, HB;

// (5) cross-check without H: does f^sigma have a root in K[x]/(f)?
P<z> := PolynomialRing(Q);
function PsiExists(s0)
    A0 := s0^3 - 3*s0^2 + 4*s0;
    g := z^2 - A0*z + s0;
    if Discriminant(g) eq 0 or not IsIrreducible(g) then return false; end if;
    K<r> := NumberField(g);
    b0 := r*s0*(r - 1); c0 := s0*(r - 1);
    ai := [1 - c0, -b0, -b0, K!0, K!0];
    if Discriminant(EllipticCurve(ai)) eq 0 then return false; end if;
    E := EllipticCurve(ai);
    assert Order(E![0, 0]) eq 11;
    sig := hom< K -> K | A0 - r >;
    PK<X> := PolynomialRing(K);
    b2, b4, b6 := Explode(bInvariants(E));
    f := 4*X^3 + b2*X^2 + 2*b4*X + b6;
    fs := 4*X^3 + sig(b2)*X^2 + 2*sig(b4)*X + sig(b6);
    assert IsIrreducible(f);
    return #Roots(ChangeRing(fs, NumberField(f))) gt 0;
end function;
hits2 := [];
for a in [-HB2..HB2] do for b in [1..HB2] do
    if GCD(a, b) ne 1 then continue; end if;
    if PsiExists(a/b) then Append(~hits2, a/b); end if;
end for; end for;
assert hits2 eq [4/5];
printf "(5) s = a/b with |a| <= %o, 1 <= b <= %o: E[2] and E^sigma[2] are G_K-isomorphic only at s = 4/5\n", HB2, HB2;
printf "SPLIT_11_REMARK_DONE: the construction of 11.b gives finitely many curves; to height %o in s only 11.b\n", HB;
quit;
