// Validation for the census proof qsimple-96-obstruction (not used in the proof).
// Run: magma -b qsimple_96_validation.m < /dev/null        (optional: B:=600 for another bound on p)
// Thirteen LMFDB curves whose Jacobians are of type (S), one for each possible Gamma and each sign of the class of
// an involution.  For every prime p < B of good reduction for the given model, the L-polynomial of the curve at p
// is asserted to have one of the shapes that [P612, Lemma 2] and Lemma 7 of the proof allow for that curve:
//   trivial Frobenius:        P = (X^2 - a X + p)^2
//   order 3, 4, 6:            the three polynomials of [P612, Lemma 2]
//   involution of class c:    P = (X^2 + s p)^2 - c t^2 X^2,  s = sign(c),  t an integer.
// The classes c of the involutions of each curve are an input of this script (they were read off from the data);
// the assertion is that no prime gives any other polynomial.  The run ends with VAL_DONE.
SetColumns(0);
SetMemoryLimit(16*10^9);
Q := Rationals();
Z := Integers();
Qx<x> := PolynomialRing(Q);
if not assigned B then B := 500; elif Type(B) eq MonStgElt then B := StringToInteger(B); end if;

// < LMFDB label, f, h, Sato-Tate group (LMFDB), Gamma, orders of the rotations, classes c of the involutions >
curves := [*
  < "576.a.576.1",       [0,-1,0,-1],       [1,1,1,1], "E_2",    "C2",  [],     [-1] >,
  < "25600.f.512000.1",  [0,2,5,-1,-5,2],   [0],       "J(E_1)", "C2",  [],     [2] >,
  < "3969.d.250047.1",   [0,1,0,-4,5,-3],   [1,1,1],   "J(E_1)", "C2",  [],     [3] >,
  < "324.a.648.1",       [0,0,1,2,2,1],     [1,1,0,1], "E_3",    "C3",  [3],    [] >,
  < "1296.a.20736.1",    [0,-1,4,-3,-1,1],  [0],       "E_3",    "C3",  [3],    [] >,
  < "256.a.512.1",       [0,-1,1,1,-3,2],   [1],       "E_4",    "C4",  [4],    [-1] >,
  < "169.a.169.1",       [0,0,0,0,1,1],     [1,1,0,1], "E_6",    "C6",  [3,6],  [-3] >,
  < "8192.b.131072.1",   [0,-4,6,0,-3,1],   [0],       "J(E_2)", "V4",  [],     [1,-1] >,
  < "102400.b.102400.1", [0,-1,0,-1,0,1],   [0],       "J(E_2)", "V4",  [],     [2,-1] >,
  < "262144.a.262144.1", [0,-1,0,-2,0,1],   [0],       "J(E_2)", "V4",  [],     [2,-1] >,
  < "2187.a.6561.1",     [-1],              [1,0,0,1], "J(E_3)", "D6",  [3],    [3] >,
  < "12800.c.128000.1",  [0,-1,3,0,-3,1],   [0],       "J(E_4)", "D8",  [4],    [-1,1,2] >,
  < "6075.a.18225.1",    [1,0,0,1],         [1,0,0,1], "J(E_6)", "D12", [3,6],  [-3,1,3] >
*];

// P(X) = X^4 + a1 X^3 + a2 X^2 + p a1 X + p^2.
function IsTrivialShape(a1, a2, p)
  return 4*a2 eq a1^2 + 8*p;
end function;
function IsRotationShape(n, a1, a2, p)        // [P612, Lemma 2]
  if n eq 3 then return a2 eq a1^2 - p; end if;                                   // a = a1
  if n eq 4 then return a1 mod 2 eq 0 and 2*a2 eq a1^2; end if;                   // u = -a1/2
  return a1 mod 3 eq 0 and 3*a2 eq a1^2 + 3*p;                                    // n = 6, b = a1/3
end function;
function IsInvolutionShape(c, a1, a2, p)      // Lemma 7: a1 = 0 and a2 = 2 s p - c t^2
  if a1 ne 0 then return false, _; end if;
  s := Sign(c);
  v := 2*s*p - a2;
  if v mod c ne 0 then return false, _; end if;
  fl, t := IsSquare(v div c);
  if not fl then return false, _; end if;
  return true, t;
end function;

for cv in curves do
  lab := cv[1]; rots := cv[6]; invs := cv[7];
  C := HyperellipticCurve(Qx!cv[2], Qx!cv[3]);
  disc := Z!Discriminant(C);
  count := AssociativeArray();
  keys := ["trivial"] cat [ "order " cat IntegerToString(n) : n in rots ] cat [ "involution c=" cat IntegerToString(c) : c in invs ];
  for k in keys do count[k] := 0; end for;
  nprimes := 0; nonzero_t := 0;
  for p in PrimesInInterval(3, B) do
    if disc mod p eq 0 then continue; end if;
    Lp := LPolynomial(ChangeRing(C, GF(p)));
    a1 := Z!Coefficient(Lp, 1); a2 := Z!Coefficient(Lp, 2);
    assert Z!Coefficient(Lp, 3) eq p*a1 and Z!Coefficient(Lp, 4) eq p^2;
    N := Z!Evaluate(Lp, 1);
    nprimes +:= 1;
    fits := [];
    if IsTrivialShape(a1, a2, p) then Append(~fits, "trivial"); end if;
    for n in rots do
      if IsRotationShape(n, a1, a2, p) then
        Append(~fits, "order " cat IntegerToString(n));
        if N mod 32 eq 0 then assert p mod 8 eq 1; end if;          // the congruence of Proposition 8
      end if;
    end for;
    for c in invs do
      fl, t := IsInvolutionShape(c, a1, a2, p);
      if fl then
        Append(~fits, "involution c=" cat IntegerToString(c));
        assert N eq (1 + Sign(c)*p)^2 - c*t^2;
        if t ne 0 then nonzero_t +:= 1; end if;
        if c mod 8 ne 1 and N mod 32 eq 0 then assert (1 + Sign(c)*p) mod 4 eq 0; end if;
        if c mod 3 ne 1 and N mod 3 eq 0 then assert (1 + Sign(c)*p) mod 3 eq 0; end if;
      end if;
    end for;
    assert #fits ge 1;                                               // no other polynomial occurs
    for k in fits do count[k] +:= 1; end for;
  end for;
  assert forall{ k : k in keys | count[k] ge 1 };                    // every allowed shape does occur
  printf "VAL %o  ST=%o  Gamma=%o  good primes below %o: %o;  fits:", lab, cv[4], cv[5], B, nprimes;
  for k in keys do printf "  %o: %o", k, count[k]; end for;
  printf "\n";
end for;
print "VAL_DONE";
quit;
