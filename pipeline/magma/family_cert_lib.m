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
