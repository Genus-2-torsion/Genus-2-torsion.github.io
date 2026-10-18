/* family_44_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [4,4].

   The record's family (census_infinity_record.md, section 2.6; torsion_jac/code/claude_census_infty_round4c.m):
   "Richelot partner #5" of the [2,2,4,4] family C_s : y^2 = f_s = x(x+A)(x+B)(x+C)(x+D) of
   family_2244_certificate.m.  Magma's RichelotIsogenousSurfaces lists the 15 (2,2)-isogenous surfaces in an
   order that depends on the factorization order of f over Q, so an index alone does not define a family.
   At the record's fibers s = 3, 5, 7, 9, 11, partner #5 is the quotient by the kernel
       K = { [(-C,0) - oo], [(0,0) - (-B,0)], [(-A,0) - (-D,0)] }
   (checked by Q-isomorphism of the codomains, torsion_jac notes/oldcert_A_2026_10_18.md); this certificate uses that kernel
   over Q(s).  K meets 2[2,2,4,4] = <[x, 0], [(x+A)(x+B), 0]> trivially, so the quotient carries [4,4].
   The codomain (Richelot; Cassels-Flynn, Prolegomena, ch. 9): G1 = x + C, G2 = x(x+B), G3 = (x+A)(x+D),
       H1 = G2'G3 - G2G3',  H2 = G3'G1 - G3G1',  H3 = G1'G2 - G1G2',  Del = det(coefficients of G1, G2, G3),
       C'_s : y^2 = f'_s = Del H1 H2 H3   (degree 6, coefficients in Q[s]).
   Argument.
   (1) Marking.  The classes D0, D1 of family_2244_certificate.m (2 D0 = [x, 0], 2 D1 = [(x+A)(x+B), 0])
       are pushed along the Richelot correspondence G1(x)H1(z) + G2(x)H2(z) = 0, y y' = Del G1(x)H1(z)(x - z),
       giving degree-4 divisors on C'_s and classes E1, E2 in J'(Q(s)).  (How E1, E2 are found plays no
       role in the proof.)  Mumford checks over Q(s): 4 E1 = 4 E2 = 0, and the 16 sums i E1 + j E2 are
       distinct, so [4,4] is marked.
   (2) Point counts cannot certify this family: J' is isogenous to J, so #J'(F_p) = #J(F_p) is divisible
       by #[2,2,4,4] = 64 at every good p.  The group structure does it.  Every 23-, 59- and 73-integral
       s = 5 mod 99061 (= 23*59*73) has the reductions of s = 5 at 23, 59 and 73, squarefree of degree 6;
       J'(Q)_tors embeds in each J'(F_p) (p odd, good reduction), so its 2-part has 2-rank <= 2 (p = 23),
       exponent dividing 4 (p = 59), and contains [4,4]: it is [4,4].  The odd parts of #J'(F_p) have gcd 1.
       The Frobenius at 73 is strict, so these Jacobians are geometrically simple.
   (3) The record's fiber s = 5 lies in the class and has TorsionSubgroup [4,4]; the G2-invariants of the
       members differ, so s -> M_2 is nonconstant.
   (The log lines call the family parameter t; here it is s.)
   Run (from the repository root):
       magma -b pipeline/magma/family_44_certificate.m > data/knowledge/sources/family_44_certificate.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
load "pipeline/magma/family_cert_lib.m";

// ---- local helpers: Richelot codomain and pushforward, group-structure residue certificate ----
function RichData(G1, G2, G3)
    K := BaseRing(Parent(G1));
    H1 := Derivative(G2)*G3 - G2*Derivative(G3);
    H2 := Derivative(G3)*G1 - G3*Derivative(G1);
    H3 := Derivative(G1)*G2 - G1*Derivative(G2);
    M := Matrix(K, 3, 3, [ Coefficient(G, i) : i in [0..2], G in [G1, G2, G3] ]);
    return H1, H2, H3, Determinant(M);
end function;
// image of the effective divisor {u(x) = 0, y = v(x)} under the correspondence: (U'(z), V'(z)), deg U' = 2 deg u
function PushEff(u, v, G1, G2, H1, H2, Del)
    Px := Parent(u); K := BaseRing(Px);
    Pz<z> := PolynomialRing(K);
    PXZ := PolynomialRing(Pz);
    lift := func< g | PXZ![ Pz!c : c in Coefficients(g) ] >;
    uz := lift(u);
    Phi := lift(G1)*Evaluate(H1, z) + lift(G2)*Evaluate(H2, z);
    U := Resultant(uz, Phi); U := U/LeadingCoefficient(U);
    assert Degree(U) eq 2*Degree(u);
    if Degree(u) eq 1 then
        Xz := Pz!(-Coefficient(u, 0)/Coefficient(u, 1));
    else
        rem := Phi mod uz;
        g, ia, _ := XGCD(Coefficient(rem, 1) mod U, U); assert g eq 1;
        Xz := (-Coefficient(rem, 0)*ia) mod U;
    end if;
    vX := Evaluate(lift(v), Xz) mod U;
    g, iv, _ := XGCD(vX, U); assert g eq 1;
    V := (Del*Evaluate(lift(G1), Xz)*Evaluate(H1, z)*(Xz - z)*iv) mod U;
    return Px!Eltseq(U), Px!Eltseq(V);
end function;
// exponents of 2 in the invariant factors (decreasing): the type of the 2-primary part
function TwoType(inv)
    return Reverse(Sort([ Valuation(n, 2) : n in inv | n mod 2 eq 0 ]));
end function;
// does the abelian 2-group of type a embed in the one of type b?
function EmbedsType(a, b)
    return #a le #b and &and[ a[i] le b[i] : i in [1..#a] ];
end function;
// all types of subgroups of the abelian 2-group of type b
function SubTypes(b)
    out := [ [Integers()| ] ];
    for i in [1..#b] do
        new := [];
        for a in out do
            if #a eq i - 1 then
                if i eq 1 then mx := b[1]; else mx := Min(b[i], a[i-1]); end if;
                for e in [1..mx] do Append(~new, a cat [e]); end for;
            end if;
        end for;
        out cat:= new;
    end for;
    return out;
end function;
function TypeToInv(a) return Reverse([ 2^e : e in a ]); end function;
// residue certificate by group structure, for a 2-group G (invariants Ginv)
function GroupStructureCertificate(f, primes, residues, members, Ginv)
    types := []; odds := []; anystrict := false;
    for i in [1..#primes] do
        p := primes[i]; r0 := residues[i];
        fbar := ReductionClass(f, p, r0);
        for m in members do
            assert (GF(p)!m) eq (GF(p)!r0);
            assert PolynomialRing(GF(p))!SpecPoly(f, m) eq fbar;
        end for;
        n, chi := FrobData(fbar);
        st := StrictChi(chi); anystrict := anystrict or st;
        inv := Invariants(AbelianGroup(Jacobian(HyperellipticCurve(fbar))));
        assert &*inv eq n;
        Append(~types, TwoType(inv)); Append(~odds, n div 2^Valuation(n, 2));
        printf "p = %o, t = %o mod %o: reduction %o (squarefree, degree %o), J(F_%o) = %o (order %o; 2-part %o), chi = %o%o\n",
            p, r0, p, fbar, Degree(fbar), p, inv, n, TypeToInv(TwoType(inv)), chi, st select " (strict)" else "";
    end for;
    G2 := TwoType(Ginv);
    cands := [ h : h in SubTypes(types[1]) | EmbedsType(G2, h) and &and[ EmbedsType(h, b) : b in types ] ];
    printf "abelian 2-groups containing %o that embed in every J(F_p)[2^oo]: %o\n", Ginv, [ TypeToInv(h) : h in cands ];
    assert cands eq [ G2 ];
    printf "odd parts of #J(F_p): %o, gcd %o\n", odds, GCD(odds);
    assert GCD(odds) eq 1 and anystrict and &*Ginv eq 2^Valuation(&*Ginv, 2);
    m0 := members[1];
    inv0 := TorsionInvariants(SpecPoly(f, m0));
    printf "sanity: t = %o has TorsionSubgroup %o\n", m0, inv0;
    assert inv0 eq Ginv;
    assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } ge 2;
    printf "G2-invariants of the members %o are not all equal\n", members;
    return true;
end function;

// ---- the source family and its marked classes (family_2244_certificate.m) ----
Fs<s> := FunctionField(Q); P<x> := PolynomialRing(Fs);
a := (s-2)*(s^2-s+1/2)*(s^2-3/2)*(s^4+s^2+9/4);
b := (s^2-2*s+1/2)*(s^2-2*s+9/2)*(s^2-3/2)*(s^2+1);
c := -(s-2)*(2*s^2-2*s+1)*(s^4+s^2+9/4);
d := (s-2)*(2*s^2-2*s+1)*(2*s^2-3)*(s^2+1);
A := a^2; B := b^2; C := c^2; D := d^2;
f := x*(x+A)*(x+B)*(x+C)*(x+D);
J := Jacobian(HyperellipticCurve(f));
aa := [a, b, c, d];
s1 := &+aa; s2 := &+[ aa[i]*aa[j] : i, j in [1..4] | i lt j ];
s3 := &+[ aa[i]*aa[j]*aa[k] : i, j, k in [1..4] | i lt j and j lt k ]; s4 := &*aa;
D0 := elt< J | x^2 - s2*x + s4, (s1*s2 - s3)*x - s1*s4 >;
assert 2*D0 eq elt< J | x, 0 >;
okA, yA := IsSquare((C-A)*(D-A)); okB, yB := IsSquare((C-B)*(D-B));
okC, wC := IsSquare((C-A)*(C-B)); okD, wD := IsSquare((D-A)*(D-B));
assert okA and okB and okC and okD;
R := (x+A)*(x+B); F := x*(x+C)*(x+D);
sols := [];
for e2, e3, e4, e5 in [1, -1] do
    rows := [ [Fs| 1, -A, A^2, 0], [Fs| 1, -B, B^2, 0], [Fs| 1, 0, 0, -e3*a*b],
              [Fs| 1, -C, C^2 + e4*wC*C, -e4*wC], [Fs| 1, -D, D^2 + e5*wD*D, -e5*wD] ];
    rhs := [Fs| a*yA, e2*b*yB, 0, 0, 0];
    ok, sol := IsConsistent(Transpose(Matrix(Fs, 5, 4, &cat rows)), Vector(Fs, rhs));
    if ok then
        U0, U1, k, l0 := Explode(Eltseq(sol));
        if k ne 0 then
            U := k*x^2 + U1*x + U0; l := k*x + l0;
            if F eq R*l^2 - U^2 then Append(~sols, <U/k, l>); end if;
        end if;
    end if;
end for;
assert #sols ge 1;
D1 := elt< J | sols[1][1], (R*sols[1][2]) mod sols[1][1] >;
assert 2*D1 eq elt< J | R, 0 >;
print "source (family_2244_certificate.m): 2 D0 = [x, 0] and 2 D1 = [(x+A)(x+B), 0] over Q(s)";

// ---- the codomain for the kernel {oo - (-C), 0 - (-B), (-A) - (-D)} ----
G1 := x + C; G2 := x*(x+B); G3 := (x+A)*(x+D);
assert G1*G2*G3 eq f;
H1, H2, H3, Del := RichData(G1, G2, G3);
assert Del ne 0;
fp := Del*H1*H2*H3;
assert Degree(fp) eq 6 and &and[ Denominator(cc) eq 1 : cc in Coefficients(fp) ];
printf "codomain f'_s = Del H1 H2 H3: degree 6, coefficients in Q[s] of degrees %o\n", [ Degree(Numerator(cc)) : cc in Coefficients(fp) ];
printf "H1 = %o\nH2 = %o\nH3 = %o\nDel = %o\n", H1, H2, H3, Del;
Jp := Jacobian(HyperellipticCurve(fp));
U, V := PushEff(D0[1], D0[2], G1, G2, H1, H2, Del); E1 := elt< Jp | U, V, 4 >;
U, V := PushEff(D1[1], D1[2], G1, G2, H1, H2, Del); E2 := elt< Jp | U, V, 4 >;
printf "E1 = %o\nE2 = %o\n", E1, E2;
assert 4*E1 eq Jp!0 and 4*E2 eq Jp!0;
assert 2*E1 eq elt< Jp | H3/LeadingCoefficient(H3), 0 > and 2*E2 eq elt< Jp | H1/LeadingCoefficient(H1), 0 >;
assert #{ i*E1 + j*E2 : i, j in [0..3] } eq 16;
print "over Q(s): 4 E1 = 4 E2 = 0, 2 E1 = [H3, 0], 2 E2 = [H1, 0], and the 16 sums i E1 + j E2 are distinct: [4,4] is marked on J'";
// the record's partner #5 at its fiber s = 5 (informational: Magma's ordering of partners may change)
f5 := SpecPoly(f, 5);
Rs := RichelotIsogenousSurfaces(Jacobian(HyperellipticCurve(f5)));
printf "fiber s = 5: C'_5 is Q-isomorphic to Magma's RichelotIsogenousSurfaces(J_5)[5] (the record's partner #5): %o\n",
    Type(Rs[5]) eq JacHyp and IsIsomorphic(HyperellipticCurve(SpecPoly(fp, 5)), Curve(Rs[5]));

// ---- residue certificate by group structure; sanity; moduli ----
ok := GroupStructureCertificate(fp, [23, 59, 73], [5, 5, 5], [Q| 5, 99066, -99056], [4,4]);
print "FAMILY_44_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [4,4]";
quit;
