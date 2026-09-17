/* trivial_torsion_family.m -- the L1-L3 certificate (see data/knowledge/sources/census_infinity_record.md,
   section 1) that infinitely many geometrically simple genus-2 Jacobians over Q have trivial rational
   torsion: the family y^2 = x^5 + x + t.
     L1  exactness: fibers with J(Q)_tors = {0} (Magma, TorsionSubgroup);
     L2  geometric simplicity: a strict prime per fiber (irreducible Frobenius polynomial, no degree
         drop of pi^n for n <= 12), as in pipeline/magma/verify_lib.m;
     L3  nonconstant modulus: the fibers have pairwise distinct G2-invariants, and the base is P^1.
   Run from the repository root:  magma -b pipeline/magma/trivial_torsion_family.m > data/knowledge/sources/trivial_torsion_family.log
*/
SetColumns(0);
SetMemoryLimit(4*10^9);
R<x> := PolynomialRing(Rationals());
Tp<T> := PolynomialRing(Rationals());
Z := Integers();

function IsStrict(chi)
    if Degree(chi) ne 4 or not IsIrreducible(chi) then return false; end if;
    K := NumberField(chi); pi := K.1;
    for n in [2..12] do
        if Degree(MinimalPolynomial(pi^n)) lt 4 then return false; end if;
    end for;
    return true;
end function;

function StrictPrime(g, pmax)
    D := Z!Discriminant(g);
    for p in PrimesInInterval(3, pmax) do
        if D mod p eq 0 or (Z!LeadingCoefficient(g)) mod p eq 0 then continue; end if;
        gp := PolynomialRing(GF(p))!g;
        if Degree(gp) lt 5 or not IsSquarefree(gp) then continue; end if;
        // Al := "Naive": Magma's default L-polynomial algorithm is wrong for some curves over F_3, F_5, F_7 (see verify_lib.m)
        chi := Tp!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp) : Al := "Naive")));
        if IsStrict(chi) then return p, chi; end if;
    end for;
    return 0, Tp!0;
end function;

invs := [];
t := 0;
while #invs lt 5 do
    t +:= 1;
    g := x^5 + x + t;
    if not IsIrreducible(g) then printf "t = %o: x^5 + x + t reducible, rational 2-torsion, skipped\n", t; continue; end if;
    C := HyperellipticCurve(g);
    I := Invariants(TorsionSubgroup(Jacobian(C)));
    p, chi := StrictPrime(g, 200);
    G2 := G2Invariants(C);
    printf "t = %o: y^2 = %o, J(Q)_tors = %o, strict prime %o, chi = %o, G2 = %o\n", t, g, I, p, chi, G2;
    if I ne [] or p eq 0 then printf "  (not a certificate fiber)\n"; continue; end if;
    Append(~invs, G2);
end while;
assert #SequenceToSet(invs) eq #invs;
printf "TRIVIAL_TORSION_FAMILY_CERTIFIED: %o fibers of y^2 = x^5 + x + t, all with trivial torsion, geometrically simple, pairwise distinct G2-invariants\n", #invs;
quit;
