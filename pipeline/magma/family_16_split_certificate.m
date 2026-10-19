/* family_16_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E1 x E2) with J(Q)_tors exactly [16].
   Found 2026-10-19 (torsion_jac notes/pow2_2026_10_19.md).

   Source.  The Howe-Leprevost-Poonen (2,2)-gluing (HLP 2000, Prop. 4) of E1 = 15a4 (torsion Z/8) with the
   curves E2 that have one rational 2-torsion point and the same 2-torsion field, with gain g = 2 (torsion_jac
   notes gain2_loci_2026_10_15.md, section 2: [2] x [8] with g = 2 gives [16]; the census curve 16.b is the
   gluing of 15a4 and 180a2).  In the monic model Y^2 = (X - 3/4)(X^2 + 2X + 37) of E1 (Y = y + (x + 1)/2) the
   2-torsion field is K = Q(i), the roots are e0 = 3/4, th1 = -1 + 6i and its conjugate, and P8 = (2, 6) has
   order 8.  Its Kummer class is delta(P8) = (x - e0, x - th1) = (N(beta), beta) with beta = 2 - th1 = 3 - 6i.
   A curve E2 : y^2 = x(x - th2)(x - th2'), th2 in K, glued along alpha : (e0, th1, th1') -> (0, th2, th2'),
   has the gain iff the Kummer class (th2 th2', -th2) of (0,0) on E2 equals delta(P8), i.e. iff
   -th2 = beta gamma^2 with gamma in K*.  With gamma = 1 + s i (gamma is defined up to Q*, and s = infinity is
   gamma = i) this gives, over Q(s),
       E2(s) : y^2 = x (x^2 + (6 + 24s - 6s^2) x + 45 (s^2 + 1)^2),     th2 = -(3 - 6i)(1 + s i)^2,
   and the gluing is (derived here; HLP's formula evaluated over K(s), checked below)
       C_s : y^2 = f_s(x) = (s^2 - 11s - 1)((s^2 + s - 1)x^2 + 1)(144(s^2 + 1)^2 x^4 + 24(11s^2 + 4s - 11)x^2 + 125).
   phi1 : (x, y) -> (x^2, y) maps C_s onto Ea : Y^2 = F(X) (f_s = F(x^2)), Q(s)-isomorphic to E2(s), and
   phi2 : (x, y) -> (1/x^2, y/x^3) onto Eb : Y^2 = X^3 F(1/X), Q(s)-isomorphic to E1.
   The member s = 0 (gamma = 1, E2 = 180a2) is the census curve 16.b.

   Argument.
   (1) Splitting.  Wherever f_s is squarefree of degree 6, phi1^* + phi2^* : Ea x Eb -> J is an isogeny over Q
       with kernel of order 4, so J is isogenous over Q to E1 x E2(s).
   (2) The marked Z/16 over Q(s).  Let F be the splitting field of z^4 - 6z^2 + 45 (degree 8; it contains i,
       sqrt(beta) and sqrt(beta')).  Over F the four halves of P8 on E1 are defined (DivisionPoints), and over
       F(s) the four halves of T = (0,0) on E2(s) are (x, Y) = (r r', r r'(r + r')) with r^2 = -th2 = beta gamma^2,
       r = +-sqrt(beta) gamma, r' = +-sqrt(beta') gamma'.  For a half R1 of P8 and a half R2 of T the class
       D16 = phi2^*(R1) + phi1^*(R2) is computed by Mumford arithmetic over F(s); for suitable halves its
       coefficients lie in Q(s) (Galois descent made explicit), and over Q(s): 16 D16 = 0 and 8 D16 != 0.
       So Z/16 is marked over Q(s), and it embeds in J_s(Q)_tors at every smooth member (the marked sections
       specialize injectively on smooth fibres).
   (3) Exactness on a residue class.  For s = 1 mod 7 and s = 1 mod 23 (s 7- and 23-integral) the curve has the
       reductions of s = 1, with #J(F_7) = 80 and #J(F_23) = 432; J_s(Q)_tors injects into both (p odd, good
       reduction), so its order divides gcd(80, 432) = 16 = #<D16>, and J_s(Q)_tors is exactly Z/16.
   (4) Infinitely many.  The G2-invariants of the members differ, so s -> M_2 is nonconstant with finite fibres,
       and the residue class gives infinitely many isomorphism classes.
   (5) Sanity: TorsionSubgroup at s = 1 and s = 162 is [16]; the member s = 0 is Q-isomorphic to the census
       curve 16.b (it lies outside the class); the member s = -2, also outside it, has torsion [48] (E2(-2) has
       a rational point of order 3), so the residue class is needed.
   Run (from the repository root):
       magma -b pipeline/magma/family_16_split_certificate.m > family_16_split_certificate.log */
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

// ---------------------------------------------------------------- (0) E1, E2(s) and the HLP gluing over K(s)
E1 := EllipticCurve([Q| 1, 1, 1, 35, -28]);
assert Conductor(E1) eq 15 and Invariants(TorsionSubgroup(E1)) eq [8];
P8 := E1![2, 6];
assert Order(P8) eq 8;
PQ := RX; z := X;
a1, a2, a3, a4, a6 := Explode(aInvariants(E1));
g1 := z^3 + (a2 + a1^2/4)*z^2 + (a4 + a1*a3/2)*z + (a6 + a3^2/4);
assert g1 eq (z - 3/4)*(z^2 + 2*z + 37);
print "E1 = 15a4 = [1,1,1,35,-28], torsion Z/8, P8 = (2,6); monic model Y^2 = (X - 3/4)(X^2 + 2X + 37), 2-torsion field Q(i)";
K<i> := QuadraticField(-1);
Ks<sK> := FunctionField(K); PK<XK> := PolynomialRing(Ks);
Fs<s> := FunctionField(Q); Px<x> := PolynomialRing(Fs);
PKs := PolynomialRing(K); PQs := PolynomialRing(Q);
function ToQ(c)         // an element of K(s) with rational coefficients, as an element of Q(s)
    nu := Numerator(c); de := Denominator(c); lc := LeadingCoefficient(de);
    nu := PKs!(nu/lc); de := PKs!(de/lc);
    assert &and[ IsCoercible(Q, a) : a in Coefficients(nu) cat Coefficients(de) ];
    return Evaluate(PQs![Q!a : a in Coefficients(nu)], s)/Evaluate(PQs![Q!a : a in Coefficients(de)], s);
end function;
th1 := -1 + 6*i; th1c := -1 - 6*i; beta := 2 - th1; betac := 2 - th1c;
assert beta eq 3 - 6*i;
gamK := 1 + sK*i; gamKc := 1 - sK*i;
th2 := -beta*gamK^2; th2c := -betac*gamKc^2;
E2 := EllipticCurve([Fs| 0, 6 + 24*s - 6*s^2, 0, 45*(s^2 + 1)^2, 0]);
assert ToQ(th2 + th2c) eq -(6 + 24*s - 6*s^2) and ToQ(th2*th2c) eq 45*(s^2 + 1)^2;
print "E2(s) : y^2 = x(x^2 + (6 + 24s - 6s^2)x + 45(s^2+1)^2) = x(x - th2)(x - th2'), th2 = -(3 - 6i)(1 + s i)^2";
al := [Ks| 3/4, th1, th1c]; be := [Ks| 0, th2, th2c];
aa1 := (al[3]-al[2])^2/(be[3]-be[2]) + (al[2]-al[1])^2/(be[2]-be[1]) + (al[1]-al[3])^2/(be[1]-be[3]);
bb1 := (be[3]-be[2])^2/(al[3]-al[2]) + (be[2]-be[1])^2/(al[2]-al[1]) + (be[1]-be[3])^2/(al[1]-al[3]);
aa2 := al[1]*(be[3]-be[2]) + al[2]*(be[1]-be[3]) + al[3]*(be[2]-be[1]);
bb2 := be[1]*(al[3]-al[2]) + be[2]*(al[1]-al[3]) + be[3]*(al[2]-al[1]);
Df := ((al[1]-al[2])*(al[1]-al[3])*(al[2]-al[3]))^2;
Dg := ((be[1]-be[2])*(be[1]-be[3])*(be[2]-be[3]))^2;
AA := Dg*aa1/aa2; BB := Df*bb1/bb2;
hK := -(AA*(al[2]-al[1])*(al[1]-al[3])*XK^2 + BB*(be[2]-be[1])*(be[1]-be[3]))
      *(AA*(al[3]-al[2])*(al[2]-al[1])*XK^2 + BB*(be[3]-be[2])*(be[2]-be[1]))
      *(AA*(al[1]-al[3])*(al[3]-al[2])*XK^2 + BB*(be[1]-be[3])*(be[3]-be[2]));
hHLP := Px![ ToQ(c) : c in Coefficients(hK) ];
f := (s^2 - 11*s - 1)*((s^2 + s - 1)*x^2 + 1)*(144*(s^2 + 1)^2*x^4 + 24*(11*s^2 + 4*s - 11)*x^2 + 125);
oksq, rt := IsSquare(Fs!(hHLP/f));
assert oksq;
print "C_s : y^2 = f_s(x) = (s^2 - 11s - 1)((s^2 + s - 1)x^2 + 1)(144(s^2+1)^2 x^4 + 24(11s^2 + 4s - 11)x^2 + 125)";
printf "the HLP (Prop. 4) curve, evaluated over Q(i)(s), has coefficients in Q(s) and is (%o)^2 f_s\n", rt;

// ---------------------------------------------------------------- (1) the two elliptic quotients
F3 := Px![ Coefficient(f, 2*k) : k in [0..3] ];
Fr := Px![ Coefficient(F3, 3 - k) : k in [0..3] ];
assert Evaluate(F3, x^2) eq f;
Ea, ca := MonEC(F3); Eb, cb := MonEC(Fr);
E1s := ChangeRing(E1, Fs);
ok1, m1 := IsIsomorphic(E1s, Eb); ok2, m2 := IsIsomorphic(E2, Ea);
assert ok1 and ok2;
print "(1) phi1 = (x^2, y): C_s -> Ea : Y^2 = F(X) ~ E2(s) and phi2 = (1/x^2, y/x^3): C_s -> Eb ~ E1 over Q(s): J is isogenous over Q to E1 x E2(s)";

// ---------------------------------------------------------------- (2) the marked Z/16 over Q(s)
J := Jacobian(HyperellipticCurve(f));
B8 := elt< J | x^2 - 1/(m1(E1s!Eltseq(P8))[1]/cb), ((m1(E1s!Eltseq(P8))[2]/cb)/(m1(E1s!Eltseq(P8))[1]/cb))*x >;
assert Ord2(B8) eq 8;
d1 := IsomorphismData(m1); d2 := IsomorphismData(m2);
im := m1(E1s!Eltseq(P8)); r, ss, tt, uu := Explode(d1);
conv := (uu^2*P8[1] + r eq im[1] and uu^3*P8[2] + ss*uu^2*P8[1] + tt eq im[2]) select 2 else 1;
if conv eq 1 then assert (P8[1]-r)/uu^2 eq im[1] and (P8[2] - ss*(P8[1]-r) - tt)/uu^3 eq im[2]; end if;
function ApplyIso(dat, pt, sub)
    r, ss, tt, uu := Explode([sub(c) : c in dat]);
    if conv eq 1 then return [(pt[1]-r)/uu^2, (pt[2] - ss*(pt[1]-r) - tt)/uu^3]; end if;
    return [uu^2*pt[1] + r, uu^3*pt[2] + ss*uu^2*pt[1] + tt];
end function;
F := SplittingField(z^4 - 6*z^2 + 45);
assert Degree(F) eq 8;
okI, iF := IsSquare(F!(-1)); okb, rb := IsSquare(F!(3 - 6*iF)); okbc, rbc := IsSquare(F!(3 + 6*iF));
assert okI and okb and okbc;
FF<w> := FunctionField(F); PF<xF> := PolynomialRing(FF);
sub := func< c | Evaluate(Numerator(c), w)/Evaluate(Denominator(c), w) >;
JF := Jacobian(HyperellipticCurve(PF![ sub(c) : c in Coefficients(f) ]));
caF := sub(ca); cbF := sub(cb);
function PullF1(pt) im := ApplyIso(d2, pt, sub); return elt< JF | xF^2 - im[1]/caF, PF!(im[2]/caF) >; end function;
function PullF2(pt) im := ApplyIso(d1, pt, sub); XX := im[1]/cbF; YY := im[2]/cbF; return elt< JF | xF^2 - 1/XX, (YY/XX)*xF >; end function;
halvesP := DivisionPoints(ChangeRing(E1, F)!Eltseq(P8), 2);
assert #halvesP eq 4;
gam := 1 + w*iF; gamc := 1 - w*iF;
th2F := -(3 - 6*iF)*gam^2; th2Fc := -(3 + 6*iF)*gamc^2;
E2F := EllipticCurve([FF| 0, sub(6 + 24*s - 6*s^2), 0, sub(45*(s^2 + 1)^2), 0]);
assert th2F + th2Fc eq -sub(6 + 24*s - 6*s^2) and th2F*th2Fc eq sub(45*(s^2 + 1)^2);
halvesT := [];
for sb in [1, -1] do for sc in [1, -1] do
    rr := sb*rb*gam; rc := sc*rbc*gamc;
    assert rr^2 eq -th2F and rc^2 eq -th2Fc;
    R2 := E2F![rr*rc, rr*rc*(rr + rc)];
    assert 2*R2 eq E2F![0, 0];
    Append(~halvesT, R2);
end for; end for;
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
assert Ord2(D16) eq 16;
BT := elt< J | x^2 - m2(E2![0, 0])[1]/ca, 0 >;         // phi1^*((0,0)), of order 2
assert Ord2(BT) eq 2 and 2*D16 eq B8 + BT;
printf "(2) over F(s), F = splitting field of z^4 - 6z^2 + 45 (degree 8): D16 = phi2^*(R1) + phi1^*(R2), 2R1 = P8, 2R2 = (0,0), has coefficients in Q(s):\n     D16 = %o\n", D16;
print "    over Q(s): 16 D16 = 0 and 8 D16 != 0, so Z/16 is marked; 2 D16 = phi2^*(P8) + phi1^*((0,0))";

// ---------------------------------------------------------------- (3) exactness on the class s = 1 mod 161
primes := [7, 23]; residues := [1, 1];
sres := CRT(residues, primes); M := &*primes;
members := [Q| 1, 1 + M, 1 - M, (1 + M)/(1 + 2*M)];
nJ := [];
for k in [1..#primes] do
    p := primes[k]; r0 := residues[k];
    fbar := ReductionClass(f, p, r0);
    for m in members do
        assert Valuation(m, p) ge 0 and GF(p)!m eq GF(p)!r0;
        assert PolynomialRing(GF(p))!SpecPoly(f, m) eq fbar;
    end for;
    n, chi := FrobData(fbar);
    Append(~nJ, n);
    printf "(3) p = %o, s = %o mod %o: reduction %o (squarefree, degree %o), #J(F_%o) = %o, chi = %o\n",
        p, r0, p, fbar, Degree(fbar), p, n, chi;
end for;
assert nJ eq [80, 432] and GCD(nJ) eq 16;
printf "(3) gcd of #J(F_p) = %o = #<D16>: every s = %o mod %o, 7- and 23-integral, has J(Q)_tors exactly Z/16\n", GCD(nJ), sres, M;

// ---------------------------------------------------------------- (4) infinitely many, (5) sanity
assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } eq #members;
printf "(4) the G2-invariants of the members %o are pairwise distinct\n", members;
for m in members[1..2] do
    invm := TorsionInvariants(SpecPoly(f, m));
    printf "(5) TorsionSubgroup at s = %o: %o\n", m, invm;
    assert invm eq [16];
end for;
C16b := HyperellipticCurve(RX![-125, 0, 97, 0, -26, 0, 2], RX![0, 1, 0, 1]);     // census 16.b
assert G2Invariants(HyperellipticCurve(SpecPoly(f, 0))) eq G2Invariants(C16b);
assert IsIsomorphic(HyperellipticCurve(SpecPoly(f, 0)), C16b);
E2_0 := EllipticCurve([Q| 0, 6, 0, 45, 0]);
assert CremonaReference(E2_0) eq "180a2";
print "(5) the member s = 0 (E2(0) = 180a2) is Q-isomorphic to the census curve 16.b (outside the class above)";
inv48 := TorsionInvariants(SpecPoly(f, -2));
assert inv48 eq [48] and Invariants(TorsionSubgroup(EllipticCurve([Q| 0, 6 - 48 - 24, 0, 45*25, 0]))) eq [6];
print "(5) outside the class, s = -2 has TorsionSubgroup [48] (E2(-2) has a rational point of order 3)";
printf "time %o s\n", Cputime(T0);
print "FAMILY_16_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E1 x E2) and J(Q)_tors exactly [16]";
quit;
