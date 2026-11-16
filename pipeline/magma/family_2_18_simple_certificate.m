/* 16 November 2026. Standalone certificate for infinitely many geometrically
   simple genus-2 Jacobians with exact rational torsion [2,18].
   Source construction: relation rows (3,2;0,0),(0,3;1,1).
   The helper library below is copied from the existing census family_cert_lib.m.
   No changes to the live census are made by this script. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);
/* family_cert_lib.m -- shared helpers for the residue-class family certificates
   (family_<group>_certificate.m).  Load from the repository root:
       load "pipeline/magma/family_cert_lib.m";

   The argument the certificates share.  Let C_t : y^2 = f_t with f in Q(t)[x] be a one-parameter family,
   and G a subgroup of J(Q(t)) marked over Q(t) (checked by Mumford arithmetic over Q(t)).
   - Reduction.  Let p be an odd prime and r0 in Z/p.  Write each coefficient of f as N(t)/D(t) with
     N, D in Z[t].  If D(r0) is not 0 mod p, then for every p-integral t = r0 mod p the coefficient is
     p-integral and congruent to N(r0)/D(r0).  So all those f_t have the same reduction fbar mod p.  If fbar
     is squarefree of degree 5 or 6, every such C_t has good reduction at p, and its reduction is the curve
     y^2 = fbar.
   - Torsion.  J_t(Q)_tors injects into J(F_p) (p odd and e = 1 < p - 1), so #J_t(Q)_tors divides #J(F_p).
     The marked torsion sections specialize injectively on smooth fibers (J[N] is finite etale over the
     smooth locus), so G embeds in J_t(Q)_tors.  If #J(F_p) = #G, or more generally the gcd over the chosen
     primes is #G, then J_t(Q)_tors is G.
   - Simplicity.  If the Frobenius polynomial chi of fbar is strict (irreducible, with pi^n of degree 4 for
     n = 2..12), then J_t is geometrically simple, since a decomposition of J_t over Qbar would reduce to one
     of J(F_p) over Fbar_p.
   - Infinitely many.  If two members of the residue class have different G2-invariants, the moduli map
     t -> M_2 is nonconstant, so it has finite fibers, and the infinitely many t in the class give
     infinitely many isomorphism classes.
   With several primes, use the residue class given by the Chinese remainder theorem. */

Q := Rationals(); Z := Integers();
RX<X> := PolynomialRing(Q);

// Strict Frobenius polynomial: irreducible, and pi^n has degree 4 for n = 2..12.
function StrictChi(chi)
    if not IsIrreducible(chi) then return false; end if;
    K := NumberField(chi);
    return &and[ Degree(MinimalPolynomial(K.1^n)) eq 4 : n in [2..12] ];
end function;

// Specialize f in Q(t)[x] at a rational t0.
function SpecPoly(f, t0)
    return RX![ Evaluate(c, t0) : c in Coefficients(f) ];
end function;

// c in Q(t): integer polynomials N, D in Z[t] with c = N/D.
function IntNumDen(c)
    n := Numerator(c); d := Denominator(c);
    k := LCM([ Denominator(e) : e in Coefficients(n) cat Coefficients(d) ]);
    PZ := PolynomialRing(Z);
    return PZ!(k*n), PZ!(k*d);
end function;

// The common reduction mod p of the f_t with t p-integral and t = r0 mod p.  Asserts that it exists
// (denominators are units) and is squarefree of degree 5 or 6.
function ReductionClass(f, p, r0)
    Fp := GF(p); Pp := PolynomialRing(Fp);
    cs := [];
    for c in Coefficients(f) do
        N, D := IntNumDen(c);
        dval := Fp!Evaluate(D, r0);
        assert dval ne 0;
        Append(~cs, (Fp!Evaluate(N, r0))/dval);
    end for;
    fbar := Pp!cs;
    assert Degree(fbar) ge 5 and IsSquarefree(fbar);
    return fbar;
end function;

// #J(F_p) and the Frobenius polynomial of y^2 = fbar (naive point counting: Magma's default
// LPolynomial is unreliable over very small fields).
function FrobData(fbar)
    Lp := LPolynomial(HyperellipticCurve(fbar) : Al := "Naive");
    chi := RX!(X^4*Evaluate(RX!Lp, 1/X));
    return Z!Evaluate(Lp, 1), chi;
end function;

// Integral simplified model, for TorsionSubgroup.
function TorsionInvariants(fQ)
    C := IntegralModel(SimplifiedModel(HyperellipticCurve(fQ)));
    return Invariants(TorsionSubgroup(Jacobian(C)));
end function;

/* Run the residue-class certificate for f in Q(t)[x] at the classes t = r_i mod p_i, and print it.
   members: rationals in the combined class, for the sanity checks; at least two.  The first is also
   checked with TorsionSubgroup.  Returns the gcd of the #J(F_p) and whether some prime is strict. */
function ResidueCertificate(f, primes, residues, members, Gorder, Ginv)
    nJ := []; anystrict := false;
    for i in [1..#primes] do
        p := primes[i]; r0 := residues[i];
        fbar := ReductionClass(f, p, r0);
        for m in members do
            assert (GF(p)!m) eq (GF(p)!r0);
            assert PolynomialRing(GF(p))!SpecPoly(f, m) eq fbar;
        end for;
        n, chi := FrobData(fbar);
        st := StrictChi(chi);
        anystrict := anystrict or st;
        Append(~nJ, n);
        printf "p = %o, t = %o mod %o: reduction %o (squarefree, degree %o), #J(F_%o) = %o, chi = %o%o\n",
            p, r0, p, fbar, Degree(fbar), p, n, chi, st select " (strict)" else "";
    end for;
    g := GCD(nJ);
    printf "gcd of #J(F_p) = %o; #G = %o\n", g, Gorder;
    assert g eq Gorder and anystrict;
    m0 := members[1];
    inv := TorsionInvariants(SpecPoly(f, m0));
    printf "sanity: t = %o has TorsionSubgroup %o\n", m0, inv;
    assert inv eq Ginv;
    assert #{ G2Invariants(HyperellipticCurve(SpecPoly(f, m))) : m in members } ge 2;
    printf "G2-invariants of the members %o are not all equal\n", members;
    return g, anystrict;
end function;

F<k>:=FunctionField(Q);P<x>:=PolynomialRing(F);
s:=(k^2+3)/4;r:=(k^2+3)/(k-1)^2;kap:=(s-1)^2/s^3;
X1:=x^3*(x-1)^2;X2:=(x-1)^3*(x-r)*(x-s);
R0:=x^3-(x-1)*(x-r)*(x-s);
D:=(x-1)^2;S:=kap*R0;A:=(D+S)/2;B:=(S-D)/2;
f:=A^2-kap*X1;
assert f eq B^2-kap*X2;
assert Evaluate(B,r) eq 0 and Evaluate(B,s) eq 0;
assert r ne s and Degree(f) eq 5 and Discriminant(f) ne 0;
J:=Jacobian(HyperellipticCurve(f));zero:=J!0;
P0:=elt<J|x,Evaluate(A,0)>;P1:=elt<J|x-1,Evaluate(A,1)>;
T1:=elt<J|x-r,0>;T2:=elt<J|x-s,0>;
assert 3*P0+2*P1 eq zero and 3*P1+T1+T2 eq zero;
assert 9*P0 eq zero and 3*P0 ne zero;
assert 2*T1 eq zero and 2*T2 eq zero;
assert T1 ne zero and T2 ne zero and T1+T2 ne zero;
D18:=P0+T1;
assert 18*D18 eq zero and 9*D18 ne zero and 6*D18 ne zero;
assert T2 ne 9*D18;
print "GENERIC_MARKING: C9 plus C2^2, hence C2 x C18";
print "RELATIONS: 3 P0 + 2 P1 = 0, 3 P1 + T1 + T2 = 0";
g,strict:=ResidueCertificate(f,[5,17],[0,0],[Q|0,85,-85,170],36,[2,18]);
// Direct absolute-simplicity check, with no bound on extension degrees.
Rz<z>:=PolynomialRing(Q);Fz:=FieldOfFractions(Rz);Rw<w>:=PolynomialRing(Fz);
chiw:=w^4+w^3+16*w^2+17*w+289;
scaled:=w^4+z*w^3+16*z^2*w^2+17*z^3*w+289*z^4;
h1:=289*z^4+527*z^3+800*z^2+527*z+289;
h2:=17*z^4+18*z^3+3*z^2+18*z+17;
assert IsIrreducible(RX![289,17,16,1,1]);
assert Resultant(chiw,scaled) eq 17^4*(z-1)^4*h1*h2^2;
assert IsIrreducible(h1) and IsIrreducible(h2);
// h1 and h2 are primitive irreducible integer polynomials, not monic.
// They cannot be cyclotomic. Thus the only root-of-unity root ratio is 1.
assert GCD([Z!c:c in Coefficients(h1)]) eq 1;
assert GCD([Z!c:c in Coefficients(h2)]) eq 1;
assert LeadingCoefficient(h1) ne 1 and LeadingCoefficient(h2) ne 1;
print "ABSOLUTE_SIMPLICITY: no distinct Frobenius root ratio is a root of unity";
C0:=ReducedMinimalWeierstrassModel(HyperellipticCurve(SpecPoly(f,0)));
print "SMALL_MEMBER",C0;
print "FAMILY_2_18_SIMPLE_CERTIFIED";
quit;
