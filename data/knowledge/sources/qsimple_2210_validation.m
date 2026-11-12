// Supplementary check (not used in the proof) for "Rational 5-torsion and the 2-rank of geometrically split Jacobians",
// Remark 2: for the eight LMFDB genus-2 curves with rational torsion (Z/2)^3 that are simple over Q and split over a
// quadratic field K, check #J(F_p) = #E(F_p^2) at inert p (product formula at split p) for the LMFDB elliptic factor E
// over K, find in the 2-isogeny class of E over K the curves with full K-rational 2-torsion (Lemma 4 predicts at least
// one), and check the conclusion of Lemma 6 for each of them (it must hold, since 8 | #J(F_p)).  Finally check the
// example after Lemma 6.  Every check is an assert; the run ends with VAL_DONE.
// Run: magma -b qsimple_2210_validation.m < /dev/null
SetColumns(0);
SetMemoryLimit(16*10^9);
Q := Rationals();
PQ<x> := PolynomialRing(Q);
data := [
 <"2457.b.199017.1", [0,-9,-13,1,4,-1], [0,1,1], [1,-1,1], "0,1;1,0;0,1;-5,0;4,-1">,
 <"2560.a.5120.1", [0,1,1,-2,-1,1], [], [1,0,1], "0,0;1,-1;0,0;1,0;1,-1">,
 <"2169.a.175689.1", [0,-1,-14,22,-9,1], [0,1,1], [1,-1,1], "1,0;0,1;0,1;5,-5;-2,-1">,
 <"2560.a.819200.1", [0,2,7,2,-7,2], [], [1,0,1], "1,1;-1,-1;0,0;0,-2;-1,2">,
 <"4096.b.65536.1", [0,-1,0,0,0,1], [], [-2,0,1], "0,0;-1,1;0,0;-4,-4;-6,-4">,
 <"4352.b.278528.1", [0,2,1,-4,-1,2], [], [1,0,1], "1,1;1,-1;1,1;-1,-1;-1,0">,
 <"25600.f.512000.1", [0,2,5,-1,-5,2], [], [1,0,1], "0,0;-1,1;0,0;3,6;5,-1">,
 <"8192.b.131072.1", [0,-4,6,0,-3,1], [], [1,0,1], "0,0;1,1;0,0;0,1;0,0"> ];
function twotors(E)
  K := BaseRing(E);
  a1,a2,a3,a4,a6 := Explode(aInvariants(E));
  PK<z> := PolynomialRing(K);
  b2,b4,b6 := Explode(bInvariants(E));
  return [ r[1] : r in Roots(4*z^3 + b2*z^2 + 2*b4*z + b6) ];
end function;
for D in data do
  lab := D[1];
  C0 := HyperellipticCurve(PQ!D[2], PQ!D[3]);
  C := SimplifiedModel(C0);
  J := Jacobian(C);
  K<a> := NumberField(PQ!D[4]);
  OK := Integers(K);
  ai := [ K![ StringToInteger(c) : c in Split(w, ",") ] : w in Split(D[5], ";") ];
  E0 := EllipticCurve(ai);
  printf "==== %o: J(Q)_tors = %o, K = Q(sqrt(%o))\n", lab, Invariants(TorsionSubgroup(J)), Squarefree(Integers()!Discriminant(PQ!D[4]));
  assert Invariants(TorsionSubgroup(J)) eq [2,2,2];
  bad := &*BadPrimes(C0) * Integers()!Norm(Discriminant(E0)) * 2;
  // (1) point counts: inert p -> #J(F_p) = #E0(F_{p^2}); split p -> product
  okI := true; okS := true; ni := 0; ns := 0;
  for p in PrimesInInterval(3, 400) do
    if bad mod p eq 0 then continue; end if;
    NJ := #BaseChange(J, GF(p));
    dec := Decomposition(OK, p);
    if #dec eq 1 and dec[1][2] eq 1 then
      P := dec[1][1]; k, mk := ResidueClassField(P);
      Ek := EllipticCurve([ mk(c) : c in aInvariants(E0) ]);
      ni +:= 1; if NJ ne #Ek then okI := false; end if;
    elif #dec eq 2 then
      prod := 1;
      for dd in dec do k, mk := ResidueClassField(dd[1]); prod *:= #EllipticCurve([ mk(c) : c in aInvariants(E0) ]); end for;
      ns +:= 1; if NJ ne prod then okS := false; end if;
    end if;
  end for;
  printf "  #J(F_p) = #E0(F_p^2) at %o inert p: %o; = product at %o split p: %o\n", ni, okI, ns, okS;
  assert okI and okS and ni gt 0 and ns gt 0;
  // (2) 2-isogeny BFS over K
  curves := [E0]; frontier := [E0];
  for depth in [1..6] do
    newf := [];
    for E in frontier do
      for r in twotors(E) do
        PK<z> := PolynomialRing(K);
        E2, phi := IsogenyFromKernel(E, z - r);
        E2 := WeierstrassModel(E2);
        if not exists{ F : F in curves | jInvariant(F) eq jInvariant(E2) and IsIsomorphic(F, E2) } then
          Append(~curves, E2); Append(~newf, E2);
        end if;
      end for;
    end for;
    frontier := newf;
    if #frontier eq 0 then break; end if;
  end for;
  full := [ E : E in curves | #twotors(E) eq 3 ];
  printf "  2-isogeny class over K: %o curves; with full K-rational 2-torsion: %o\n", #curves, #full;
  assert #full ge 1;
  // (3) the square-class conditions for the full-2-torsion curves
  d := Q!Squarefree(Integers()!Discriminant(PQ!D[4]));
  for E in full do
    rts := twotors(E);
    N := [ Q!Norm(rts[1]-rts[2]), Q!Norm(rts[1]-rts[3]), Q!Norm(rts[2]-rts[3]) ];
    sq := func< t | t ne 0 and IsSquare(t) >;
    // the conclusion of Lemma 6: some N_ij is a square, and N_ik or N_jk or d*N_ik*N_jk is a square
    two := (sq(N[1]) and sq(N[2])) or (sq(N[1]) and sq(N[3])) or (sq(N[2]) and sq(N[3]));
    dprod := (sq(N[1]) and sq(d*N[2]*N[3])) or (sq(N[2]) and sq(d*N[1]*N[3])) or (sq(N[3]) and sq(d*N[1]*N[2]));
    printf "    squarefree N12,N13,N23 = %o, d = %o: two norms squares %o, N_ij and d*N_ik*N_jk squares %o\n", [ Squarefree(Numerator(t)*Denominator(t))*Sign(t) : t in N ], d, two, dprod;
    assert two or dprod;
  end for;
end for;
// The example after Lemma 6: y^2 = x(x-1)(x-sqrt 2) over Q(sqrt 2) satisfies the conclusion of Lemma 6 only through
// d*N13*N23, and 8 divides #E(F_p^2) at every inert p < 400.
K<a> := NumberField(x^2 - 2);
OK := Integers(K);
E := EllipticCurve([0, -(1+a), 0, a, 0]);       // x(x-1)(x-a) = x^3 - (1+a)x^2 + a x
rts := [K!0, K!1, a];
assert Set(twotors(E)) eq Set(rts);
N := [ Q!Norm(rts[1]-rts[2]), Q!Norm(rts[1]-rts[3]), Q!Norm(rts[2]-rts[3]) ];
assert N eq [1, -2, -1];
assert IsSquare(2*N[2]*N[3]) and not IsSquare(N[2]) and not IsSquare(N[3]);
cnt := 0;
for p in PrimesInInterval(3, 400) do
  dec := Decomposition(OK, p);
  if #dec eq 1 and dec[1][2] eq 1 then
    k, mk := ResidueClassField(dec[1][1]);
    assert #EllipticCurve([ mk(c) : c in aInvariants(E) ]) mod 8 eq 0;
    cnt +:= 1;
  end if;
end for;
assert cnt gt 0;
printf "==== Lemma 6 example y^2 = x(x-1)(x-sqrt 2): N12,N13,N23 = %o, d = 2; 8 | #E(F_p^2) at all %o inert p < 400\n", N, cnt;
print "VAL_DONE";
quit;
