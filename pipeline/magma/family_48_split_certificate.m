/* family_48_split_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have Jacobian split over Q (isogenous over Q to E_s x F) with J(Q)_tors exactly [48].
   Found 2026-10-19 (torsion_jac notes/hlpEC_48_2026_10_19.md).

   Source.  Howe 2015 (Bull. London Math. Soc. 47, arXiv:1407.2654v3), Theorems 3.1, 3.2 and Remark 3.3, which use
   HLP 2000, Prop. 4 (gluing) and Prop. 12 (gain).
     F   : y^2 = g(x) = x^3 - 31x^2 + 256x = 42a1 (= E_8(1/3)), fixed, Q = (32, 96) of order 8, 4Q = (0,0);
     E_s : y^2 = f(x) = x(x^2 + a x + b),  a = -8(s^4 + 42s^2 - 147)/(s^2 + 63)^2,  b = 16(s^2 + 7)^3/(s^2 + 63)^3,
           P = (4(s^2 + 7)/(s^2 + 63), 224(s^2 + 7)/(s^2 + 63)^2) of order 6, 3P = (0,0);
     psi : E_s[2] -> F[2], (alpha_i, 0) -> (beta_i, 0), with alpha_2 = 4(s + r)^3(s - 3r)/(s^2 + 63)^2, beta_2 = (31 - 3r)/2,
           r^2 = -7 (alpha_1 = beta_1 = 0, alpha_3, beta_3 the conjugates).
     C^s_24 : d y^2 = x^6 + c4 x^4 + c2 x^2 + c0, with
           c4 = -31(s^4 + 42s^2 - (32200/93)s - 147),
           c2 = 2^8(s^8 + 84s^6 - (3472/3)s^5 + 1470s^4 - 48608s^3 + 53508s^2 + 170128s + 21609),
           c0 = (2^20 7/3) s (s^2 + 7)^3 (s^2 + 63),   d = s^4 + 42s^2 + (1736/3)s - 147
     (coefficients from the PDF, p. 5; the pdftotext rendering garbles them).  Below the model is used as
     Y^2 = f_s(x) = d (x^6 + c4 x^4 + c2 x^2 + c0) (Y = d y).
     The base is Howe's D : w^2 = z^3 + 14z^2 + 196z = 9408q1 (rank 2), with
           s(z, w) = -21((z^2 + 196)(z^2 + 56z + 196) + 32(z + 14)(z - 14)w)/(z^4 - 896z^3 - 24696z^2 - 175616z + 38416).

   Argument.
   (1) Gluing and splitting.  Over Q(s) the certificate rebuilds the HLP gluing of E_s and F along psi (the Tschirnhausen
       transformation phi in Q(s)[x] with phi(alpha_i) = beta_i, computed over Q(s)(r), and the closed form
       y^2 = c2 Fh(k x^2 - 1, c2 x^2) of HLP Prop. 4, Fh(X, Z) = Z^3 f(X/Z), k = c1 + c2(alpha_1 + alpha_2 + alpha_3)) and
       checks that f_s(x) is a square times this sextic at x/(8(s^2 + 63)) (Howe's substitution).  Independently,
       (x, y) -> (x^2, y) maps C onto Y^2 = d(U^3 + c4 U^2 + c2 U + c0) and (x, y) -> (1/x^2, y/x^3) onto
       Y^2 = d(c0 U^3 + c2 U^2 + c4 U + 1); these are checked to be isomorphic over Q(s) to F and to E_s, so J is
       isogenous over Q to E_s x F wherever f_s is squarefree of degree 6.
   (2) The marked group.
       - Over Q(s): the two quotient maps have prime degree 2 and C is not of genus 1, so the pullbacks are injective.
         D6 (pullback of P) has order 6, D8 (pullback of Q) has order 8, and T8 = (pullback of 3P) + D8, the image of
         (3P, Q), has order 8.  The sum D24 = (pullback of 2P) + D8, the image of (2P, Q), has order 24: this is Howe's
         Theorem 3.1.
       - The gain (Howe, Theorem 3.2).  Let L = Q[T]/(T^2 + 7), rho = T, A = 4(s + rho)^3(s - 3rho)/(s^2 + 63)^2 and
         B = (31 - 3rho)/2.  Since f and g each have one rational root, H^1(Q, E_s[2]) and H^1(Q, F[2]) are identified
         with the group of square classes of L^x through the components at alpha_2 and beta_2.  The Kummer image of 3P = (0,0) is -A, that of Q is
         32 - B, and psi matches the two components.
         Over the function field L(D) the certificate checks A(B - 32) = l^2 n^2, with
             l = (1 - rho)(s + rho)^2/(s^2 + 63),
             n = 6((7 + 11rho)w - 4z^2 + 784)/(8z^2 + 7(1 - 3rho)z + 1568).
         So at every point of D with a smooth member, iota(3P) = psi~(iota'(Q)), and by HLP Prop. 12 the image of a
         point (P0, Q0) with 2P0 = 3P, 2Q0 = Q is a rational point G of J with 2G = T8.
       - G has order 16, and G + (pullback of 2P) has order 48: Z/48 embeds in J(Q).
       - At the members below, G is found explicitly with TorsionSubgroup on an integral model (2G = T8 is checked by
         Mumford arithmetic), and TorsionSubgroup gives [48].
   (3) Exactness on a residue class.  #J(F_p) = #E_s(F_p) #F(F_p).  At p = 5, 6 | #E_s(F_5) <= 10 and 8 | #F(F_5) <= 10
       force #J(F_5) = 48.  For every point R of D with R = P2 = (7, -49) mod 5, s(R) is 5-integral and = s(P2) = 3
       (mod 5), C_{s(R)} has the reduction of C_3 at 5, and J(F_5) = Z/48 (cyclic).  So J(Q)_tors is exactly Z/48.
   (4) Infinitely many members.  P2 has order 8 in D(F_5), so R_k = (1 + 8k) P2 reduces to P2 mod 5.  The denominator of
       s does not vanish at z = 7 mod 5, so s(R_k) = 3 mod 5 and C_{s(R_k)} has the reduction of (3).  The R_k are
       distinct, and s has degree 16 on D (Remark 3.3), so they give infinitely many s.  The G2-invariants at s = 3 (P2)
       and s = 21/31 (P3) differ, so s -> M_2 is nonconstant, and the R_k give infinitely many isomorphism classes.
   (5) Sanity.  The members at P1, P2, P3 (s = -21, 3, 21/31) are Howe's Corollary 3.5 curves (equal G2-invariants), with
       TorsionSubgroup [48].  The census curve 48.a is not a member: every member has the isogeny factor 42a1, but J(48.a)
       is isogenous to curves of conductor 21 and 84.  Its torsion is [48] (TorsionSubgroup).
   Run (from the repository root):
       magma -b pipeline/magma/family_48_split_certificate.m > data/knowledge/sources/family_48_split_certificate.log */
SetColumns(0);
SetQuitOnError(true);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";
T0 := Cputime();

function MonEC(G)                                   // Y^2 = G(U), G a cubic, as a Weierstrass curve; returns it and lc(G)
    cs := [ Coefficient(G, i) : i in [0..3] ];
    return EllipticCurve([0, cs[3], 0, cs[2]*cs[4], cs[1]*cs[4]^2]), cs[4];
end function;
function OrderUpTo(P, n)
    for k in [1..n] do if k*P eq Parent(P)!0 then return k; end if; end for; return -1;
end function;
function IsoToSquare(E1, E2)                        // c4(E1) = lam^2 c4(E2), c6(E1) = lam^3 c6(E2), lam a square
    c4a, c6a := Explode(cInvariants(E1)); c4f, c6f := Explode(cInvariants(E2));
    lam := (c6a/c6f)/(c4a/c4f);
    return c4a eq lam^2*c4f and c6a eq lam^3*c6f and IsSquare(lam);
end function;

// ---------------------------------------------------------------- (0) the fixed curve F and the family E_s
F := EllipticCurve([0, -31, 0, 256, 0]);
assert Conductor(F) eq 42 and Invariants(TorsionSubgroup(F)) eq [8];
QF := F![32, 96]; assert Order(QF) eq 8 and 4*QF eq F![0, 0];
Fs<s> := FunctionField(Q); P<x> := PolynomialRing(Fs);
a := -8*(s^4 + 42*s^2 - 147)/(s^2 + 63)^2; b := 16*(s^2 + 7)^3/(s^2 + 63)^3;
Es := EllipticCurve([0, a, 0, b, 0]);
PE := Es![4*(s^2 + 7)/(s^2 + 63), 224*(s^2 + 7)/(s^2 + 63)^2];
assert OrderUpTo(PE, 12) eq 6 and 3*PE eq Es![0, 0];
printf "(0) F = %o = 42a1, Q = (32, 96) of order 8, 4Q = (0,0); E_s : y^2 = x(x^2 + a x + b) with P of order 6, 3P = (0,0)\n", aInvariants(F);

// ---------------------------------------------------------------- (1) Howe's C^s_24, and the HLP gluing
c4 := -31*(s^4 + 42*s^2 - (32200/93)*s - 147);
c2 := 2^8*(s^8 + 84*s^6 - (3472/3)*s^5 + 1470*s^4 - 48608*s^3 + 53508*s^2 + 170128*s + 21609);
c0 := (2^20*7/3)*s*(s^2 + 7)^3*(s^2 + 63);
d := s^4 + 42*s^2 + (1736/3)*s - 147;
fs := d*(x^6 + c4*x^4 + c2*x^2 + c0);
assert Degree(fs) eq 6 and IsSquarefree(fs);
f := x*(x^2 + a*x + b); g := x^3 - 31*x^2 + 256*x;
Pr<RR> := PolynomialRing(Fs); Kr<r> := ext< Fs | RR^2 + 7 >;
al := [ Kr!0, 4*(s + r)^3*(s - 3*r)/(s^2 + 63)^2, 4*(s - r)^3*(s + 3*r)/(s^2 + 63)^2 ];
be := [ Kr!0, (31 - 3*r)/2, (31 + 3*r)/2 ];
PKr<X> := PolynomialRing(Kr);
assert &and[ Evaluate(PKr![ Kr!c : c in Coefficients(f) ], al[i]) eq 0 and Evaluate(PKr![ Kr!c : c in Coefficients(g) ], be[i]) eq 0 : i in [1..3] ];
phi := &+[ be[i]*&*[ (X - al[j])/(al[i] - al[j]) : j in [1..3] | j ne i ] : i in [1..3] ];
cs := [ Coefficient(phi, i) : i in [0..2] ];
assert &and[ IsCoercible(Fs, c) : c in cs ];
cc0, cc1, cc2 := Explode([ Fs!c : c in cs ]);
k := cc1 + cc2*(-Coefficient(f, 2));
hHLP := cc2*&+[ Coefficient(f, i)*(k*x^2 - 1)^i*(cc2*x^2)^(3 - i) : i in [0..3] ];
hl := Evaluate(hHLP, x/(8*(s^2 + 63)));
ratio := LeadingCoefficient(fs)/LeadingCoefficient(hl);
assert fs eq ratio*hl and IsSquare(ratio);
print "(1) psi : (alpha_i,0) -> (beta_i,0) is given by phi in Q(s)[x]; the closed-form HLP gluing y^2 = c2 Fh(k x^2 - 1, c2 x^2), at";
print "    x/(8(s^2 + 63)) and times a square, is Howe's f_s = d(x^6 + c4 x^4 + c2 x^2 + c0): C^s_24 is the HLP gluing of E_s and F";
Ha := P![ Coefficient(fs, 2*i) : i in [0..3] ];                    // fs = Ha(x^2)
Hb := P![ Coefficient(fs, 6 - 2*i) : i in [0..3] ];                // x^6 fs(1/x) = Hb(x^2)
Ea, la := MonEC(Ha); Eb, lb := MonEC(Hb);
assert IsoToSquare(Ea, BaseChange(F, Fs)) and IsoToSquare(Eb, Es);
print "(1) (x,y) -> (x^2, y) maps C onto Y^2 = Ha(U), isomorphic over Q(s) to F; (x,y) -> (1/x^2, y/x^3) onto Y^2 = Hb(U), isomorphic to E_s";

// ---------------------------------------------------------------- (2) the marked group over Q(s)
J := Jacobian(HyperellipticCurve(fs)); PJ := Parent((J!0)[1]); xj := PJ.1;
FK := BaseChange(F, Fs); QK := FK!Eltseq(QF);
oka, isoa := IsIsomorphic(FK, Ea);  assert oka;
okb, isob := IsIsomorphic(Es, Eb); assert okb;
Pull1 := func< Pt | elt< J | xj^2 - Pt[1]/la, PJ!(Pt[2]/la) > >;                 // pullback along (x^2, y)
Pull2 := func< Pt | elt< J | xj^2 - lb/Pt[1], PJ!((Pt[2]/lb)*(lb/Pt[1]))*xj > >;  // pullback along (1/x^2, y/x^3)
D6 := Pull2(isob(PE)); D3 := Pull2(isob(2*PE)); T2 := Pull2(isob(3*PE));
D8 := Pull1(isoa(QK));
assert OrderUpTo(D6, 30) eq 6 and OrderUpTo(D8, 30) eq 8;
T8 := T2 + D8; D24 := D3 + D8;
assert OrderUpTo(T8, 30) eq 8 and OrderUpTo(D24, 50) eq 24;
print "(2) over Q(s): D6 = pullback of P has order 6, D8 = pullback of Q order 8, T8 = image of (3P, Q) order 8,";
print "    D24 = image of (2P, Q) order 24 (Howe, Theorem 3.1)";

// ---------------------------------------------------------------- (2) the gain over L(D): A(B - 32) = l^2 n^2
L<rho> := QuadraticField(-7);
FzL<zL> := FunctionField(L); PwL<WL> := PolynomialRing(FzL);
FD<wL> := FunctionField(WL^2 - (zL^3 + 14*zL^2 + 196*zL));
sD := -21*((zL^2 + 196)*(zL^2 + 56*zL + 196) + 32*(zL + 14)*(zL - 14)*wL)/(zL^4 - 896*zL^3 - 24696*zL^2 - 175616*zL + 38416);
AD := 4*(sD + rho)^3*(sD - 3*rho)/(sD^2 + 63)^2; BD := (31 - 3*rho)/2;
lD := (1 - rho)*(sD + rho)^2/(sD^2 + 63);
nD := 6*((7 + 11*rho)*wL - 4*zL^2 + 784)/(8*zL^2 + 7*(1 - 3*rho)*zL + 1568);
assert AD*(BD - 32) eq lD^2*nD^2;
// the components at alpha_2, beta_2: A = alpha_2 and B = beta_2 under rho = r
PL<XL> := PolynomialRing(L);
assert Evaluate(XL^3 - 31*XL^2 + 256*XL, BD) eq 0;          // B = beta_2 (rho = r, as in (1))
assert Evaluate(PolynomialRing(FD)![0, Evaluate(b, sD), Evaluate(a, sD), 1], AD) eq 0;   // A = alpha_2 is a root of f on D
print "(2) over L(D), L = Q(sqrt -7): A(B - 32) = l^2 n^2 (Howe, proof of Theorem 3.2), so iota(3P) = psi~(iota'(Q)) on D and by";
print "    HLP Prop. 12 there is a rational G with 2G = T8: G has order 16 and G + D3 has order 48, so Z/48 embeds in J(Q)";

// ---------------------------------------------------------------- the base D and the residue class
D := EllipticCurve([0, 14, 0, 196, 0]);
assert Conductor(D) eq 9408;
P1 := D![0, 0]; P2 := D![7, -49]; P3 := D![16, -104];
assert Order(P1) eq 2 and Order(P2) eq 0 and Order(P3) eq 0;
NsD := func< R | -21*((R[1]^2 + 196)*(R[1]^2 + 56*R[1] + 196) + 32*(R[1] + 14)*(R[1] - 14)*R[2]) >;
DsD := func< R | R[1]^4 - 896*R[1]^3 - 24696*R[1]^2 - 175616*R[1] + 38416 >;
sAt := func< R | NsD(R)/DsD(R) >;
assert sAt(P2) eq 3 and sAt(P3) eq 21/31;
printf "(4) D = %o = 9408q1 (rank bounds %o), P1 = (0,0) of order 2, P2 = (7,-49), P3 = (16,-104); s(P2) = 3, s(P3) = 21/31\n", aInvariants(D), [RankBounds(D)];

// ---------------------------------------------------------------- (3) exactness: p = 5, s = 3 mod 5
p := 5; r0 := 3;
fbar := ReductionClass(fs, p, r0);
n5, chi5 := FrobData(fbar);
A5 := AbelianGroup(Jacobian(HyperellipticCurve(fbar)));
assert #A5 eq n5 and Invariants(A5) eq [48];
printf "(3) p = 5, s = 3 mod 5: reduction %o (squarefree, degree %o), J(F_5) = %o: every member in the class has J(Q)_tors = Z/48\n", fbar, Degree(fbar), Invariants(A5);
// residues s mod 5 with good reduction all give Z/48 (#E_s(F_5) = 6, #F(F_5) = 8 are forced)
for rr in [1..4] do assert Invariants(AbelianGroup(Jacobian(HyperellipticCurve(ReductionClass(fs, p, rr))))) eq [48]; end for;
print "    (s = 1, 2, 3, 4 mod 5 all give good reduction and J(F_5) = Z/48)";

// ---------------------------------------------------------------- (4) members R_k = (1 + 8k) P2
o5 := Order(ChangeRing(D, GF(5))!Eltseq(P2)); assert o5 eq 8;
assert Valuation(Integers()!Discriminant(D), 5) eq 0 and GF(5)!DsD(P2) ne 0;
printf "(4) P2 has order %o in D(F_5); R_k = (1 + %o k) P2 reduces to P2 mod 5, s(R_k) is 5-integral and = 3 mod 5\n", o5, o5;
R1 := (1 + o5)*P2; s1 := sAt(R1);
assert Valuation(s1 - 3, 5) ge 1;
assert PolynomialRing(GF(5))!SpecPoly(fs, s1) eq fbar;
printf "(4) sanity: s(R_1) = s(9 P2) has %o digits, = 3 mod 5, with the same reduction at 5\n", #Sprint(s1);
assert G2Invariants(HyperellipticCurve(SpecPoly(fs, Q!3))) ne G2Invariants(HyperellipticCurve(SpecPoly(fs, Q!21/31)));
print "(4) the G2-invariants at s = 3 and s = 21/31 differ: s -> M_2 is nonconstant, so the R_k give infinitely many isomorphism classes";

// ---------------------------------------------------------------- (2), (5) members: explicit G, TorsionSubgroup, Howe's Corollary 3.5
howe := [ HyperellipticCurve(RX![0, -16, 8, 14, -5, -3, 1], RX![0, 1, 1]),
          HyperellipticCurve(RX![2, -6, 10, -11, 5, -1, 1], RX![0, 1, 1]),
          HyperellipticCurve(RX![80875050306064, -60387768796, 60403483004, -31429634, 15717859, -3651, 1217], RX![0, 1, 1]) ];
function MemberCheck(s0)
    sp := func< c | Evaluate(c, s0) >;
    f0 := SpecPoly(fs, s0);
    mu := LCM([ Denominator(c) : c in Coefficients(f0) ]);
    f1 := mu^2*f0;                                                     // integral model y1 = mu y
    J1 := Jacobian(HyperellipticCurve(f1)); P1x := Parent((J1!0)[1]); X1 := P1x.1;
    tr := func< Dv | elt< J1 | P1x![ sp(c) : c in Coefficients(Dv[1]) ], mu*P1x![ sp(c) : c in Coefficients(Dv[2]) ] > >;
    t8 := tr(T8); d3 := tr(D3);
    A, mA := TorsionSubgroup(J1);
    halves := [ mA(e) : e in A | 2*mA(e) eq t8 ];
    assert #halves ge 1;
    G := halves[1];
    return Invariants(A), OrderUpTo(t8, 20), OrderUpTo(G, 40), OrderUpTo(G + d3, 100), f1, #halves;
end function;
for i in [1..3] do
    R := [P1, P2, P3][i]; s0 := i eq 1 select Q!-21 else sAt(R);
    inv, o8, o16, o48, f1, nh := MemberCheck(s0);
    eqhowe := G2Invariants(HyperellipticCurve(f1)) eq G2Invariants(howe[i]);
    printf "(5) P%o (s = %o): TorsionSubgroup %o; T8 has order %o, a rational G with 2G = T8 has order %o (%o such G), G + D3 has order %o; Howe Cor. 3.5 curve: %o\n",
        i, s0, inv, o8, o16, nh, o48, eqhowe;
    assert inv eq [48] and o8 eq 8 and o16 eq 16 and o48 eq 48 and eqhowe;
end for;
// the census curve 48.a
C48a := HyperellipticCurve(RX![3, 0, 9, 0, -4], RX![0, 1, 0, 1]);
f48a := HyperellipticPolynomials(SimplifiedModel(C48a));
T48a := TorsionInvariants(f48a);
printf "(5) census 48.a: TorsionSubgroup %o; not a member (J(48.a) is isogenous to curves of conductor 21 and 84, every member has the factor 42a1)\n", T48a;
assert T48a eq [48];
printf "time %o s\n", Cputime(T0);
print "FAMILY_48_SPLIT_CERTIFIED: infinitely many genus-2 curves with Jacobian split over Q (isogenous to E_s x F) and J(Q)_tors exactly [48]";
quit;
