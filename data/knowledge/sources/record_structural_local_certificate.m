// Exact controls for the September 13 structural audit.
// The uniform valuation and lattice proofs are in the consolidated note.
// No bounded search is used as a nonexistence certificate.
SetColumns(0); SetSeed(1); SetMemoryLimit(2*10^9);
SetClassGroupBoundGenerators(map< PowerStructure(RngOrd) -> Integers() |
    O :-> Ceiling(MinkowskiBound(O)) >);
Q:=Rationals(); R<x>:=PolynomialRing(Q);
polys:=[x*(x+3)*(x-1), x*(x-3)*(x+1),
    (x^2-1)*(x^2-9), -(x^2-1)*(x^2-9),
    x^4-6*x^2+1, -(x^4-6*x^2+1),
    -x*(x-1)*(3*x^2-3*x+1), x*(x^2-1)];
for i in [1..#polys] do
    C:=HyperellipticCurve(polys[i]);
    seed:=Points(C : Bound:=10)[1]; // choice of origin only
    E,mp:=EllipticCurve(C,seed);
    S:=TwoSelmerGroup(E); T2:=TwoTorsionSubgroup(E);
    assert #S eq #T2;
    printf "ELLIPTIC_DIRECT_DESCENT index=%o a_invariants=%o selmer_order=%o rational_two_torsion_order=%o rank_upper_bound=0\n",
        i,aInvariants(E),#S,#T2;
end for;
print "EIGHT_ELLIPTIC_DESCENTS_PASSED";

K<t>:=FunctionField(Q);
A:=(t+3)*(t-1)^3; B:=(t-3)*(t+1)^3;
D:=t^8-12*t^6+30*t^4+228*t^2+9;
assert A-B eq 16*t;
assert A^2-A*B+B^2 eq D;
jj:=D^3/(t^2*(t^2-9)^2*(t^2-1)^6);
assert jj eq 256*(A^2-A*B+B^2)^3/(A^2*B^2*(A-B)^2);
E:=EllipticCurve([K|0,-A-B,0,A*B,0]);
assert jInvariant(E) eq jj;
printf "FULL26_J_FORMULA %o\n",jj;
function ExpectedV3(u)
    h:=Valuation(u,3);
    if h lt 0 then return 6*h; end if;
    if h eq 0 then return -6*Valuation(u^2-1,3); end if;
    if h eq 1 then return -2*Valuation((u/3)^2-1,3); end if;
    return 2-2*h;
end function;
function ExpectedV2(u)
    h:=Valuation(u,2);
    if h lt 0 then return 6*h; end if;
    if h gt 0 then return -2*h; end if;
    if Valuation(u-1,2) lt 2 then u:=-u; end if;
    r:=Valuation(u-1,2); s:=Valuation(u+3,2);
    assert r ge 2 and s ge 2;
    if r eq 2 then assert s ge 3; else assert s eq 2; end if;
    return 16-6*r-2*s;
end function;
controls:={Q|a/b : a in [-25..25], b in [1..25] |
    a ne 0 and Abs(a) ne b and Abs(a) ne 3*b};
controls join:={Q|3^h*(1+3*u) : h in [-5..5], u in [-5..5] |
    3^h*(1+3*u) notin {Q|0,1,-1,3,-3}};
for u in controls do
    jv:=Evaluate(Numerator(jj),u)/Evaluate(Denominator(jj),u);
    assert Valuation(jv,3) eq ExpectedV3(u) and ExpectedV3(u) lt 0;
    assert Valuation(jv,2) eq ExpectedV2(u) and ExpectedV2(u) lt 0;
end for;
printf "FULL26_LOCAL_VALUATION_CONTROLS_PASSED count=%o\n",#controls;

// Formal projector for opposite quadratic signs of potentially
// multiplicative elliptic representations on inertia at 3.
Alg<z>:=quo<R|(x-1)^2*(x+1)^2>;
ep:=(z+z^-1+2)/4; em:=1-ep;
assert ep^2 eq ep and em^2 eq em and ep*em eq 0;
assert (z-1)^2*ep eq 0 and (z+1)^2*em eq 0;
assert (z+z^-1-2) eq (z-1)^2*z^-1;
g:=1+2*t;
assert (g+g^-1-2)/4 eq t^2/g;
assert (g+g^-1+2)/4 eq 1+t^2/g;
print "LOCAL_TWO_ADIC_PROJECTOR_IDENTITIES_PASSED";
print "STRUCTURAL_LOCAL_CERTIFICATE_DONE";
print "CERTIFIED structural local checks";
quit;
