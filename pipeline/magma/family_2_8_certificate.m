/* family_2_8_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2 curves
   over Q have geometrically simple Jacobian with J(Q)_tors exactly [2,8].

   The record's family (census_infinity_record.md, section 2.6; torsion_jac/code/claude_census_infty_round4c.m):
   "Richelot partner #4" of the [2,2,2,8] family C_t : y^2 = f_t = x(x+A)(x+B)(x+C)(x+D),
   (A,B,C,D) = (a^2,b^2,c^2,d^2), (a,b,c,d) = (-4t^2(t+1)/(t^2+t+1)^2, -t/(t+1), 1, t), of
   family_2228_certificate.m.  Magma's RichelotIsogenousSurfaces lists the 15 (2,2)-isogenous surfaces in an
   order that depends on the factorization order of f over Q, and the index does NOT define one family
   here: at the record's fibers t = 2, 3 partner #4 is the quotient by
       K = { [(-B,0) - oo], [(0,0) - (-A,0)], [(-C,0) - (-D,0)] },
   but at t = 5, 7 it is the quotient by { [(-A,0) - oo], [(0,0) - (-B,0)], [(-C,0) - (-D,0)] } (checked by
   Q-isomorphism of the codomains, torsion_jac notes/oldcert_A_2026_10_18.md).  This certificate uses the kernel K of the
   record's exact and strict fiber t = 2, over Q(t).  K does not contain [x, 0] = 4[P - oo], so the quotient
   carries [2,8].
   The codomain (Richelot; Cassels-Flynn, Prolegomena, ch. 9): G1 = x + B, G2 = x(x+A), G3 = (x+C)(x+D),
       H1 = G2'G3 - G2G3',  H2 = G3'G1 - G3G1',  H3 = G1'G2 - G1G2',  Del = det(coefficients of G1, G2, G3),
       C'_t : y^2 = f'_t = Del H1 H2 H3   (degree 6).
   Argument.
   (1) Marking.  The class [P - oo] of order 8 of family_2228_certificate.m is pushed along the Richelot
       correspondence G1(x)H1(z) + G2(x)H2(z) = 0, y y' = Del G1(x)H1(z)(x - z), giving a degree-2 divisor
       Gamma(P) on C'_t and the class E = [Gamma(P) - oo' - oo''] in J'(Q(t)).  (How E is found plays no role
       in the proof.)  Mumford checks over Q(t): E has order 8, 4E = [H3, 0], and with T = [H1, 0] the
       16 sums i E + j T are distinct, so [2,8] is marked.
   (2) Point counts cannot certify this family: J' is isogenous to J, so #J'(F_p) = #J(F_p) is divisible
       by #[2,2,2,8] = 64 at every good p.  The group structure does it.  Every 17-, 31- and 157-integral
       t = 2 mod 82739 (= 17*31*157) has the reductions of t = 2 at 17, 31 and 157, squarefree of degree 6;
       J'(Q)_tors embeds in each J'(F_p) (p odd, good reduction), and the only abelian 2-group containing
       [2,8] that embeds in J'(F_17)[2^oo], J'(F_31)[2^oo] and J'(F_157)[2^oo] is [2,8] (2-rank <= 2 at 17;
       the second invariant is 2 at 31; exponent 8 at 157).  The odd parts of #J'(F_p) have gcd 1.  The
       Frobenius at 31 is strict, so these Jacobians are geometrically simple.
   (3) The record's fiber t = 2 lies in the class and has TorsionSubgroup [2,8]; the G2-invariants of the
       members differ, so t -> M_2 is nonconstant.
   (Named family_2_8 because family_28 would be the name for the cyclic group [28].)
   Run (from the repository root):
       magma -b pipeline/magma/family_2_8_certificate.m > data/knowledge/sources/family_2_8_certificate.log */
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

// ---- the source family and its marked class (family_2228_certificate.m) ----
Ft<t> := FunctionField(Q); P<x> := PolynomialRing(Ft);
aa := [ -4*t^2*(t+1)/(t^2+t+1)^2, -t/(t+1), Ft!1, t ];
A := aa[1]^2; B := aa[2]^2; C := aa[3]^2; D := aa[4]^2;
f := x*(x+A)*(x+B)*(x+C)*(x+D);
J := Jacobian(HyperellipticCurve(f));
s1 := &+aa; s2 := &+[ aa[i]*aa[j] : i, j in [1..4] | i lt j ];
s3 := &+[ aa[i]*aa[j]*aa[k] : i, j, k in [1..4] | i lt j and j lt k ]; s4 := &*aa;
assert s2^2 eq 4*s4;
x0 := s2/2; y0 := (s1*s2 - s3)*x0 - s1*s4;
P8 := elt< J | x - x0, y0 >;
assert 8*P8 eq J!0 and 4*P8 eq elt< J | x, 0 > and 4*P8 ne J!0;
print "source (family_2228_certificate.m): [P - oo] has order 8 and 4[P - oo] = [x, 0] over Q(t)";

// ---- the codomain for the kernel {oo - (-B), 0 - (-A), (-C) - (-D)} ----
G1 := x + B; G2 := x*(x+A); G3 := (x+C)*(x+D);
assert G1*G2*G3 eq f;
H1, H2, H3, Del := RichData(G1, G2, G3);
assert Del ne 0;
fp := Del*H1*H2*H3;
assert Degree(fp) eq 6;
printf "codomain f'_t = Del H1 H2 H3 = %o\n", fp;
printf "H1 = %o\nH2 = %o\nH3 = %o\nDel = %o\n", H1, H2, H3, Del;
Jp := Jacobian(HyperellipticCurve(fp));
U, V := PushEff(P8[1], P8[2], G1, G2, H1, H2, Del); E := elt< Jp | U, V, 2 >;
printf "E = %o\n", E;
T := elt< Jp | H1/LeadingCoefficient(H1), 0 >;
assert 8*E eq Jp!0 and 4*E ne Jp!0 and 4*E eq elt< Jp | H3/LeadingCoefficient(H3), 0 >;
assert T ne Jp!0 and 2*T eq Jp!0 and T ne 4*E;
assert #{ i*E + j*T : i in [0..7], j in [0..1] } eq 16;
print "over Q(t): E has order 8, 4E = [H3, 0], T = [H1, 0] != 4E, and the 16 sums i E + j T are distinct: [2,8] is marked on J'";
// the record's partner #4 at its fiber t = 2 (informational: Magma's ordering of partners may change)
f2 := SpecPoly(f, 2);
Rs := RichelotIsogenousSurfaces(Jacobian(HyperellipticCurve(f2)));
printf "fiber t = 2: C'_2 is Q-isomorphic to Magma's RichelotIsogenousSurfaces(J_2)[4] (the record's partner #4): %o\n",
    Type(Rs[4]) eq JacHyp and IsIsomorphic(HyperellipticCurve(SpecPoly(fp, 2)), Curve(Rs[4]));

// ---- residue certificate by group structure; sanity; moduli ----
ok := GroupStructureCertificate(fp, [17, 31, 157], [2, 2, 2], [Q| 2, 82741, -82737], [2,8]);
print "FAMILY_2_8_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [2,8]";
quit;
