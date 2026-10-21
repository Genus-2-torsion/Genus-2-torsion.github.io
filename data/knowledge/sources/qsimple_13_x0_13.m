// qsimple_13_x0_13.m
// Certificate for the case "l = 2, image C_3, discriminant a rational square" of
// proofs/qsimple-13-obstruction.html (Section 9.5(c2)), with a cross-check for l = 3.
//
// Checks:
//  (1) j13(u) = (u^2+5u+13)(u^4+7u^3+20u^2+19u+1)^3/u satisfies Phi_13(j13(u), j13(13/u)) = 0
//      identically, and has degree 14 = [SL_2(Z) : Gamma_0(13)]. Hence u -> (j13(u), j13(13/u))
//      identifies P^1 with X_0(13): every elliptic curve over Q with a rational 13-isogeny has
//      j = j13(u) for some u in Q^*.
//  (2) j13(u) - 1728 = (u^2+6u+13) h(u)^2 / u with h in Q[u] of degree 6. So for u in Q^* with
//      h(u) <> 0, j13(u) - 1728 is a square iff u^3 + 6u^2 + 13u is a square.
//  (3) s^2 = u^3 + 6u^2 + 13u is the elliptic curve 52a1 (via u = X - 2, s = Y), of rank 0
//      (proved) with torsion Z/2: its rational points are infinity and (0,0), i.e. u in {infinity, 0}.
//  (4) Cross-check (not needed for the proof): the curve w^3 = u^2(u^2+5u+13), on which j13(u)
//      is a cube, has genus 2 and a Jacobian of rank 0; Chabauty0 finds only the two points at
//      infinity of its hyperelliptic model, which lie over cusps of X_0(13).
SetColumns(0);
SetMemoryLimit(16*10^9);
v1, v2, v3 := GetVersion();
printf "Magma V%o.%o-%o\n", v1, v2, v3;

Q := Rationals();
Fu<u> := FunctionField(Q);
j13 := (u^2+5*u+13)*(u^4+7*u^3+20*u^2+19*u+1)^3/u;
assert Degree(Numerator(j13)) eq 14 and Denominator(j13) eq u;

// (1) the modular polynomial identity
Phi13 := ClassicalModularPolynomial(13);
val := Evaluate(Phi13, [j13, Evaluate(j13, 13/u)]);
assert val eq 0;
print "CHECK Phi_13(j13(u), j13(13/u)) = 0 identically; deg j13 = 14";

// (2) j - 1728
P<v> := PolynomialRing(Q);
N13 := Evaluate(Numerator(j13), v);
q := N13 - 1728*v;
quo, rem := Quotrem(q, v^2+6*v+13);
assert rem eq 0;
ok, h := IsSquare(quo);
assert ok and Degree(h) eq 6;
printf "CHECK j13 - 1728 = (u^2+6u+13) h(u)^2/u with h = %o\n", h;
// h(u) = 0 would give j = 1728 (CM); h has no rational root anyway
assert #Roots(h) eq 0;

// (3) the elliptic curve s^2 = u^3 + 6u^2 + 13u
E := EllipticCurve([0, 6, 0, 13, 0]);
Em := MinimalModel(E);
print "s^2 = u^3+6u^2+13u is", aInvariants(Em), CremonaReference(Em);
assert aInvariants(Em) eq [0, 0, 0, 1, -10];
ok, iso := IsIsomorphic(E, EllipticCurve([0, 0, 0, 1, -10]));
assert ok;
// the substitution u = X - 2, s = Y
PX<X> := PolynomialRing(Q);
assert Evaluate(v^3 + 6*v^2 + 13*v, X - 2) eq X^3 + X - 10;
r, proved := Rank(Em);
assert r eq 0 and proved;
T := TorsionSubgroup(Em);
assert #T eq 2;
pts := Points(E : Bound := 1000);
assert Set(pts) eq { E!0, E![0, 0] };
print "CHECK 52a1 has rank 0 (proved) and torsion Z/2; its points are infinity and (u,s) = (0,0)";

// (4) cross-check for l = 3: j13(u) is a cube iff u^2 (u^2+5u+13) is a cube
A2<uu, ww> := AffineSpace(Q, 2);
C3 := Curve(A2, ww^3 - uu^2*(uu^2 + 5*uu + 13));
assert Genus(C3) eq 2;
ok, H3, m3 := IsHyperelliptic(C3);
assert ok;
J3 := Jacobian(H3);
lo, hi := RankBounds(J3);
printf "L3: hyperelliptic model %o, rank bounds %o %o, torsion %o\n", H3, lo, hi, Invariants(TorsionSubgroup(J3));
assert lo eq 0 and hi eq 0;
pts3 := Chabauty0(J3);
print "L3 Chabauty0:", pts3;
assert #pts3 eq 2 and &and[ p[3] eq 0 : p in pts3 ];
for p in pts3 do
  try
    print "   preimage on w^3 = u^2(u^2+5u+13):", p @@ m3;
  catch e
    print "   (preimage not computed)";
  end try;
end for;

print "ALL_CHECKS_PASSED";
quit;
