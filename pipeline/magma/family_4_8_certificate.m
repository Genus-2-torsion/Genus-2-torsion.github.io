/* family_4_8_certificate.m -- certificate that infinitely many pairwise non-isomorphic genus-2
   curves over Q have geometrically simple Jacobian with J(Q)_tors exactly [4,8], over a RATIONAL base.
   (By-product of code/simple248_family_2_4_8_certificate.m; notes/simple248_2026_10_19.md.)

   Construction.  The M_1(8,2^w) odd model of torsion_jac/notes/m18_m14_halving.md (chart n = 1),
       C : y^2 = x A(x) B(x),  t = (2R^2 + (1 - w^2) R - 2w^2)/(4(w^2 - 1)),
   restricted to the rational curve
       w = (13R - 12)/(12R - 13),   R = 3(tau^2 - 4)/(2(tau^2 - 9)),
   on which both 2-descent conditions for halving T0 = [x, 0] hold identically in Q(tau).
   Argument.
   (1) Over Q(tau): P8 = [x + R, y0] (y0^2 = f(-R)) has order 8, the tangent half H (2H = T0) has order 4,
       and the 32 sums i P8 + j H are distinct: [4,8] is marked over Q(tau).
   (2) Residue certificate (family_cert_lib.m) at tau = 4 mod 13 and tau = 4 mod 41 (tau = 4 mod 533):
       #J(F_13) = 128, #J(F_41) = 1888, gcd 32 = #[4,8], with a strict prime.
   (3) The G2-invariants of the members differ, so infinitely many isomorphism classes.
   Run (from the census repository root):
       magma -b <this file> < /dev/null > <log> */
SetColumns(0); SetMemoryLimit(16*10^9);
SetQuitOnError(true);
load "pipeline/magma/family_cert_lib.m";

Ft<tau> := FunctionField(Q); Pt<xt> := PolynomialRing(Ft);
Rt := 3*(tau^2 - 4)/(2*(tau^2 - 9));
wt := (13*Rt - 12)/(12*Rt - 13);
tt := (2*Rt^2 + (1 - wt^2)*Rt - 2*wt^2)/(4*(wt^2 - 1));
At := xt^2 + (Rt^3 + 4*Rt^2*tt + Rt - 8*Rt*tt + 4*tt)*xt + Rt^4;
Bt := (Rt + 2 + 4*tt)*xt^2 + (Rt^2 + 4*Rt + 1 + 8*tt)*xt + (2*Rt^2 + Rt + 4*tt);
fT := xt*At*Bt;
assert Degree(fT) eq 5 and IsSquarefree(fT);

// (1) the marked group over Q(tau)
Jt := Jacobian(HyperellipticCurve(fT));
ok, y0 := IsSquare(Evaluate(fT, -Rt)); assert ok;
P8 := elt< Jt | xt + Rt, Pt!y0 >;
assert 8*P8 eq Jt!0 and 4*P8 ne Jt!0;
T0 := elt< Jt | xt, 0 >;
hT := At*Bt; c := [ Coefficient(hT, i) : i in [0..4] ];
V := -Rt^2*wt;
FU := 4*(c[4] - 2*c[5]*xt)*(c[2] - 2*c[5]*xt*V) - (c[3] - c[5]*(xt^2 + 2*V))^2;   // polynomial in U
found := false;
for rt in Roots(FU) do
    U0 := rt[1];
    okM, M := IsSquare(c[4] - 2*c[5]*U0); okN, N := IsSquare(c[2] - 2*c[5]*U0*V);
    if not (okM and okN) then continue; end if;
    for e in [1, -1] do
        if hT - xt*(M*xt + e*N)^2 eq c[5]*(xt^2 + U0*xt + V)^2 then
            uH := xt^2 + U0*xt + V; vH := (xt*(M*xt + e*N)) mod uH; found := true; break rt;
        end if;
    end for;
end for;
assert found;
H := elt< Jt | uH, vH >;
assert T0 ne Jt!0 and 2*T0 eq Jt!0 and 2*H eq T0;
assert #{ i*P8 + j*H : i in [0..7], j in [0..3] } eq 32;
print "(1) over Q(tau): P8 (order 8) and the half H of [x, 0] (order 4) generate a group of order 32 = [4,8]";

// (2), (3) residue certificate at tau = 4 mod 13 and mod 41 (tau = 4 mod 533)
g, st := ResidueCertificate(fT, [13, 41], [4, 4], [Q| 4, 537, -529, 1070], 32, [4,8]);
for m in [Q| 537, -529] do
    printf "sanity: tau = %o has TorsionSubgroup %o\n", m, TorsionInvariants(SpecPoly(fT, m));
end for;
print "FAMILY_4_8_CERTIFIED: infinitely many geometrically simple Jacobians with J(Q)_tors exactly [4,8], over the rational tau-line";
quit;
