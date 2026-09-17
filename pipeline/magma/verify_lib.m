/* verify_lib.m -- Magma verification of a submitted genus-2 curve over Q for the census
   of rational torsion subgroups of genus-2 Jacobians (https://genus-2-torsion.github.io/).

   Loaded by a job file (written by pipeline/verify.py) that first defines:

     fcoeffs, hcoeffs : sequences of strings (rational numbers), coefficients in ASCENDING degree
                        of f and h in the model  y^2 + h(x) y = f(x)   (h may be [ "0" ])
     LogFile, OutFile : strings (paths); LogFile gets progress lines, OutFile one JSON object
     MemGB            : integer (optional), memory cap
     ClaimedGroup     : sequence of integers (optional), the torsion group the submitter claims
     StrictPrimeMax   : integer (optional, default 500): primes tried for the certificates
     RichelotDepth    : integer (optional, default 2)
     QuadDiscs        : sequence of integers (optional): extra squarefree d for Richelot over Q(sqrt d)
     Cover            : optional split certificate supplied by the submitter (see 5(e) below):
                        a tuple < fieldpoly (string in X, "" for Q), cubic (4 strings e3,e2,e1,e0, or
                                  [] when a Cremona label is given), p (strings), q (strings),
                                  hh (strings), note (string), cremona label (string, "" if none) >
                        coefficients ascending, elements of the field as polynomials in "w"
     ExistingModels   : list of < id, fcoeffs, hcoeffs > (optional): curves already in the census
                        with the same G2-invariants, tested for Q-isomorphism
     Generators       : optional list of < a (string in x), b (string in x), d (integer) >: generators of
                        J(Q)_tors in Mumford representation on the even model y^2 = 4f + h^2 (Magma's
                        SimplifiedModel), as printed by pipeline/magma/prepare_submission.m.  When they
                        are given, are torsion points, are independent, and the order of the group they
                        generate equals the gcd of #J(F_p) over good primes, J(Q)_tors is certified
                        without running TorsionSubgroup (the "generators" method); otherwise
                        TorsionSubgroup is run as usual.

   Every string that reaches `eval` was whitelisted by verify.py (digits, w/X, + - * / ^ ( ) spaces).

   What is verified (everything the certificate relies on is recomputed here):
     1. the model defines a smooth curve of genus 2 over Q; the even model y^2 = g, g = 4f + h^2
        (made integral by scaling) is the working model;
     2. J(Q)_tors, exactly: from submitted generators when their group's order equals the gcd of
        #J(F_p) over good primes (see Generators above), else by Stoll's algorithm as implemented
        in Magma (TorsionSubgroup); generators in Mumford representation on y^2 = 4f + h^2;
     3. geometric simplicity (Frobenius polynomials by brute-force point counting, see ChiAt):
        a good prime p with chi_p irreducible and no degree drop of pi^n for
        n <= 12 ("strict" prime; a geometrically split or CM surface shows a drop at every good
        prime, because a ratio of Frobenius eigenvalues that is a root of unity has order <= 12
        in the Galois closure of a quartic CM field) -- the criterion of the paper's
        verify_simple_certificates.m;
     4. Q-simplicity: a good prime with chi_p irreducible over Q (a Q-split surface has
        chi_p = chi_E * chi_E' at every good prime);
     5. splitness certificates: (a) an extra involution over Q (AutomorphismGroup) => split
        over Q; (b) Richelot (2,2)-isogenies over Q, depth <= RichelotDepth: a product of
        elliptic curves over Q => split over Q, an elliptic curve over a quadratic field
        (Weil restriction) => geometrically split; (c) more than one involution in the
        geometric automorphism group => geometrically split; (d) Richelot over quadratic
        fields Q(sqrt d); (e) a submitted map (x,y) -> (p/q, y*hh/q^2) from y^2 = g to a
        genus-1 curve v^2 = e3 u^3 + e2 u^2 + e1 u + e0 over a number field, checked by the
        polynomial identity (e3 p^3 + e2 p^2 q + e1 p q^2 + e0 q^3) q = g hh^2 with a
        nonzero Wronskian and nonzero discriminant of the cubic;
     6. the "split signature" (evidence only): at every good prime below 200 the strictness
        test fails;
     7. invariants: Igusa-Clebsch and G2-invariants, the discriminant of the working model.
   The reduced minimal Weierstrass model and the conductor are computed by conductor_lib.m in a
   separate job (they need the discriminant factored, which may not terminate quickly).
*/

if not assigned MemGB then MemGB := 16; end if;
if Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB * 10^9);
if not assigned StrictPrimeMax then StrictPrimeMax := 500; end if;
if not assigned RichelotDepth then RichelotDepth := 2; end if;
if not assigned QuadDiscs then QuadDiscs := []; end if;
if not assigned ClaimedGroup then ClaimedGroup := []; end if;
if not assigned ExistingModels then ExistingModels := [* *]; end if;
if not assigned Generators then Generators := [* *]; end if;
if not assigned MaxGeneratedOrder then MaxGeneratedOrder := 20000; end if;

T0 := Cputime();
procedure Log(s)
  Write(LogFile, Sprintf("[%o s] %o", RealField(6)!Cputime(T0), s));
end procedure;

// ---------------------------------------------------------------- tiny JSON writer
function JStr(s)
  t := "";
  for i in [1..#s] do
    c := s[i];
    if c eq "\"" then t cat:= "\\\"";
    elif c eq "\\" then t cat:= "\\\\";
    elif c eq "\n" then t cat:= "\\n";
    else t cat:= c; end if;
  end for;
  return "\"" cat t cat "\"";
end function;

function JVal(v)
  case Type(v):
    when MonStgElt: return JStr(v);
    when RngIntElt: return IntegerToString(v);
    when BoolElt:   return v select "true" else "false";
    when FldRatElt: return JStr(Sprint(v));
    when SeqEnum:   return "[" cat (#v eq 0 select "" else &cat[JVal(v[i]) cat (i lt #v select "," else "") : i in [1..#v]]) cat "]";
    when List:
      if #v gt 0 and &and[Type(v[i]) eq Tup and #v[i] eq 2 and Type(v[i][1]) eq MonStgElt : i in [1..#v]] then
        return "{" cat &cat[JStr(v[i][1]) cat ":" cat JVal(v[i][2]) cat (i lt #v select "," else "") : i in [1..#v]] cat "}";
      end if;
      return "[" cat (#v eq 0 select "" else &cat[JVal(v[i]) cat (i lt #v select "," else "") : i in [1..#v]]) cat "]";
  end case;
  return JStr(Sprint(v));
end function;

function JObj(pairs)
  return "{" cat (#pairs eq 0 select "" else &cat[JStr(pairs[i][1]) cat ":" cat JVal(pairs[i][2]) cat (i lt #pairs select "," else "") : i in [1..#pairs]]) cat "}";
end function;

procedure Emit(pairs)
  Write(OutFile, JObj(pairs) : Overwrite := true);
end procedure;

procedure Fail(reason, pairs)
  Log("FAIL: " cat reason);
  Emit([* <"ok", false>, <"error", reason> *] cat pairs);
end procedure;

Null := [* <"null", true> *];   // placeholder: verify.py maps {"null": true} to null

// ---------------------------------------------------------------- helpers
Q := Rationals(); Z := Integers();
R<x> := PolynomialRing(Q);
Tp := R;   // Frobenius polynomials live in the same ring (Magma caches PolynomialRing(Q)); printed with variable T

function CoeffStrings(f)   // polynomial over Q -> ascending coefficient strings ("0" for the zero polynomial)
  if f eq 0 then return ["0"]; end if;
  return [Sprint(c) : c in Coefficients(f)];
end function;

function PolyString(f, var)
  P<xx> := PolynomialRing(BaseRing(Parent(f)));
  s := Sprint(P!Eltseq(f));
  return SubstituteString(s, "xx", var);
end function;

function FactoredString(N)
  if N eq 0 then return "0"; end if;
  parts := [];
  if N lt 0 then Append(~parts, "-1"); end if;
  for t in Factorization(Abs(N)) do
    Append(~parts, t[2] eq 1 select Sprint(t[1]) else Sprintf("%o^%o", t[1], t[2]));
  end for;
  if #parts eq 0 then return "1"; end if;
  return &cat[parts[i] cat (i lt #parts select " * " else "") : i in [1..#parts]];
end function;

// Frobenius characteristic polynomial of Jac(y^2 = g) at p, or 0 if p is not admissible.
// The L-polynomial is computed by brute-force point counting (Al := "Naive"): Magma V2.29-4's default
// algorithm returns a wrong L-polynomial for a small proportion of genus-2 curves over F_3, F_5, F_7
// (found while seeding this census: for y^2 = -11x^6 - 16x^5 + 8x^4 + 30x^3 + 8x^2 - 16x - 11 over F_5
// the default gives 25T^4 + 4T^2 + 1, while #C(F_5) = 10, #C(F_25) = 36 give 25T^4 + 20T^3 + 13T^2 + 4T + 1).
// The result is cross-checked against Magma's point enumeration, the functional equation and the
// Weil bounds; on any inconsistency the prime is treated as inadmissible and the run is aborted.
function ChiAt(g, D, p)
  if p eq 2 or D mod p eq 0 then return Tp!0, ""; end if;
  if (Z!LeadingCoefficient(g)) mod p eq 0 then return Tp!0, ""; end if;
  gp := PolynomialRing(GF(p))!g;
  if Degree(gp) lt 5 or not IsSquarefree(gp) then return Tp!0, ""; end if;
  Cp := HyperellipticCurve(gp);
  L := LPolynomial(Cp : Al := "Naive");
  c := Coefficients(L);
  if #c ne 5 then return Tp!0, Sprintf("L-polynomial at p = %o has degree %o", p, #c - 1); end if;
  a1 := -c[2]; a2 := c[3];
  if c[4] ne -p*a1 or c[5] ne p^2 then return Tp!0, Sprintf("L-polynomial at p = %o violates the functional equation", p); end if;
  if Abs(a1) gt 4*Sqrt(p) or Abs(a2) gt 6*p then return Tp!0, Sprintf("L-polynomial at p = %o violates the Weil bounds", p); end if;
  if p + 1 - a1 ne #Points(Cp) then return Tp!0, Sprintf("L-polynomial at p = %o disagrees with the point count", p); end if;
  return Tp!Reverse(c), "";
end function;

// root-power strictness of chi: irreducible and no degree drop of pi^n for n <= 12
function IsStrict(chi)
  if Degree(chi) ne 4 or not IsIrreducible(chi) then return false; end if;
  K := NumberField(chi); pi := K.1;
  for n in [2..12] do
    if Degree(MinimalPolynomial(pi^n)) lt 4 then return false; end if;
  end for;
  return true;
end function;

// E an elliptic curve over a quadratic field K: is E isomorphic over K to the base change of an
// elliptic curve E0 over Q?  Then Res_{K/Q} E = E0 x E0^(D) over Q (D the discriminant of K).
// Method: j(E) must be rational; E is then a quadratic twist by some d in K of the curve E1/Q with
// that j-invariant (j <> 0, 1728), and E0 exists iff d lies in Q^* K^{*2}, i.e. d/sigma(d) = e^2
// with N(e) = 1 (then e = f/sigma(f) by Hilbert 90 and d/f^2 is rational).  The result is checked
// by IsIsomorphic at the end, so only the "yes" answer is used.
function DescendsToQ(E)
  K := BaseRing(E);
  if Type(K) eq FldRat or Degree(K) ne 2 or Type(BaseRing(K)) ne FldRat then return false, 0; end if;
  j := jInvariant(E);
  ok, jQ := IsCoercible(Q, j);
  if not ok or jQ eq 0 or jQ eq 1728 then return false, 0; end if;
  E1 := EllipticCurveFromjInvariant(jQ);
  E1K := BaseChange(E1, K);
  ok, d := IsQuadraticTwist(E, E1K);
  if not ok then return false, 0; end if;
  sigma := Automorphisms(K)[2];
  if sigma(K.1) eq K.1 then sigma := Automorphisms(K)[1]; end if;
  ok, e := IsSquare(d / sigma(d));
  if not ok or Norm(e) ne 1 then return false, 0; end if;
  f := e eq -1 select K.1 - sigma(K.1) else 1 + e;
  ok, q := IsCoercible(Q, d / f^2);
  if not ok then return false, 0; end if;
  E0 := QuadraticTwist(E1, q);
  if not IsIsomorphic(BaseChange(E0, K), E) then return false, 0; end if;
  return true, E0;
end function;

function SplitNode(L)   // a product of elliptic curves (SetCart) or a Weil restriction (CrvEll)?
  for s in L do
    if Type(s) eq SetCart or Type(s) eq CrvEll then return true, s; end if;
  end for;
  return false, 0;
end function;

function FieldName(K)   // Q, Q(sqrt(d)), or the defining polynomial
  if Type(K) eq FldRat then return "Q"; end if;
  if Degree(K) eq 2 and Type(BaseRing(K)) eq FldRat then
    d := Discriminant(MaximalOrder(K));
    d0, _ := SquarefreeFactorization(Abs(d));
    return Sprintf("Q(sqrt(%o))", Sign(d)*d0);
  end if;
  return "the field defined by " cat Sprint(DefiningPolynomial(K));
end function;

function DescribeNode(s)
  if Type(s) eq SetCart then
    return "product of two elliptic curves over " cat FieldName(BaseRing(s[1])) cat
      (Type(BaseRing(s[1])) eq FldRat select Sprintf(" (conductors %o and %o)", Conductor(s[1]), Conductor(s[2])) else "");
  elif Type(s) eq CrvEll then
    return "Weil restriction of an elliptic curve E over " cat FieldName(BaseRing(s)) cat " (the surface is E x E^sigma)";
  end if;
  return Sprint(Type(s));
end function;

// ---------------------------------------------------------------- 1. the curve
Log("start");
ok := true; msg := "";
try
  fpol := R![Q | eval s : s in fcoeffs];
  hpol := R![Q | eval s : s in hcoeffs];
catch e
  ok := false; msg := "could not parse the coefficients: " cat Sprint(e`Object);
end try;
if not ok then Fail(msg, [* *]); quit; end if;
g0 := 4*fpol + hpol^2;
if Degree(g0) notin {5, 6} then
  Fail(Sprintf("4f + h^2 has degree %o; a genus-2 model needs degree 5 or 6", Degree(g0)), [* *]); quit;
end if;
if not IsSquarefree(g0) then Fail("the curve is singular (4f + h^2 is not squarefree)", [* *]); quit; end if;
ok := true;
try
  C := HyperellipticCurve(fpol, hpol);
catch e
  ok := false; msg := Sprint(e`Object);
end try;
if not ok then Fail("HyperellipticCurve failed: " cat msg, [* *]); quit; end if;
if Genus(C) ne 2 then Fail(Sprintf("the curve has genus %o, not 2", Genus(C)), [* *]); quit; end if;
// integral even model y^2 = g used for reductions and certificates
dd := LCM([Z | Denominator(c) : c in Coefficients(g0)]);
g := R!(dd^2 * g0);
Cg := HyperellipticCurve(g);
D := Z!Discriminant(g);
Log(Sprintf("curve: y^2 + (%o) y = %o; working model y^2 = %o", hpol, fpol, g));

// ---------------------------------------------------------------- 2. torsion
Js := Jacobian(SimplifiedModel(C));
// upper bound: #J(Q)_tors divides gcd of #J(F_p) over good primes p >= 3 (reduction is injective on torsion)
gcdJ := 0; usedp := [];
chis := AssociativeArray();
procedure ChiInto(~chis, g, D, p)
  chi, msg := ChiAt(g, D, p);
  if msg ne "" then Fail("Magma point-count inconsistency: " cat msg, [* *]); quit; end if;
  chis[p] := chi;
end procedure;
for p in PrimesInInterval(3, 200) do
  ChiInto(~chis, g, D, p); chi := chis[p];
  if chi eq 0 then continue; end if;
  gcdJ := GCD(gcdJ, Z!Evaluate(chi, 1));
  Append(~usedp, p);
  if #usedp ge 25 then break; end if;
end for;
Log(Sprintf("#J(Q)_tors divides gcd #J(F_p) = %o (p = %o)", gcdJ, usedp));

tors_method := "";
gens := [* *];
inv := [];
ntors := 0;
gen_note := "";
if #Generators gt 0 then
  // fast path: submitted generators
  ok := true; msg := "";
  pts := [];
  try
    for gen in Generators do
      a := R!(eval gen[1]); b := R!(eval gen[2]); d := gen[3];
      if Type(d) eq MonStgElt then d := StringToInteger(d); end if;
      Append(~pts, elt<Js | a, b, d>);
    end for;
  catch e
    ok := false; msg := "a submitted generator is not a point of the Jacobian of y^2 = 4f + h^2: " cat Sprint(e`Object);
  end try;
  if ok then
    orders := [Order(P) : P in pts];
    if exists{o : o in orders | o eq 0} then
      ok := false; msg := "a submitted generator has infinite order";
    elif exists{o : o in orders | o eq 1} then
      ok := false; msg := "a submitted generator is the identity";
    end if;
  end if;
  if ok then
    N := &*orders;
    if N gt MaxGeneratedOrder then
      ok := false; msg := Sprintf("the submitted generators would generate a group of order %o; too large to check", N);
    end if;
  end if;
  if ok then
    // independence: the map Z/o_1 x ... x Z/o_r -> J is injective iff all N sums are distinct
    S := {};
    tuples := CartesianProduct([[0..o-1] : o in orders]);
    for t in tuples do
      Include(~S, &+[Js | t[i]*pts[i] : i in [1..#pts]]);
    end for;
    if #S ne N then
      ok := false; msg := Sprintf("the submitted generators are not independent (they generate a group of order %o, not %o)", #S, N);
    end if;
  end if;
  if not ok then
    Fail("generators rejected: " cat msg, [* *]); quit;
  end if;
  invg := Invariants(AbelianGroup(orders));
  if gcdJ eq N then
    inv := invg; ntors := N; tors_method := "generators";
    for i in [1..#pts] do
      Append(~gens, [* <"order", orders[i]>, <"mumford", [Sprint(pts[i][1]), Sprint(pts[i][2]), Sprint(pts[i][3])]> *]);
    end for;
    Log(Sprintf("J(Q)_tors = %o: the %o submitted generators are independent torsion points generating a group of order %o = gcd #J(F_p)", inv, #pts, N));
  else
    gen_note := Sprintf("the submitted generators generate a subgroup of order %o, but the point-count bound is %o; TorsionSubgroup was run", N, gcdJ);
    Log(gen_note);
  end if;
end if;
if tors_method eq "" then
  t1 := Cputime();
  Tg, mT := TorsionSubgroup(Js);
  inv := Invariants(Tg); ntors := #Tg; tors_method := "TorsionSubgroup";
  Log(Sprintf("J(Q)_tors = %o (TorsionSubgroup, %o s)", inv, RealField(6)!Cputime(t1)));
  for i in [1..Ngens(Tg)] do
    P := mT(Tg.i);
    Append(~gens, [* <"order", Order(Tg.i)>, <"mumford", [Sprint(P[1]), Sprint(P[2]), Sprint(P[3])]> *]);
  end for;
  if gcdJ ne 0 and gcdJ mod ntors ne 0 then
    Fail(Sprintf("internal inconsistency: #J(Q)_tors = %o does not divide gcd #J(F_p) = %o", ntors, gcdJ), [* *]); quit;
  end if;
end if;
claim_ok := #ClaimedGroup eq 0 or [Z | c : c in ClaimedGroup] eq inv;

// ---------------------------------------------------------------- 3. geometric simplicity (strict prime)
strictp := 0; strictchi := Tp!0; qsimplep := 0; qsimplechi := Tp!0;
nadm := 0; nred := 0;
stricts := AssociativeArray();
for p in PrimesInInterval(3, StrictPrimeMax) do
  if not IsDefined(chis, p) then ChiInto(~chis, g, D, p); end if;
  chi := chis[p];
  if chi eq 0 then continue; end if;
  nadm +:= 1;
  if IsIrreducible(chi) then
    if qsimplep eq 0 then qsimplep := p; qsimplechi := chi; end if;
    stricts[p] := IsStrict(chi);
    if stricts[p] then strictp := p; strictchi := chi; break; end if;
  else
    nred +:= 1; stricts[p] := false;
  end if;
end for;
if strictp ne 0 then
  Log(Sprintf("geometrically simple: strict prime %o, chi = %o", strictp, PolyString(strictchi, "T")));
else
  Log(Sprintf("no strict prime below %o (%o admissible primes, %o with reducible chi)", StrictPrimeMax, nadm, nred));
end if;
if qsimplep ne 0 then Log(Sprintf("Q-simple: chi irreducible at p = %o", qsimplep)); end if;

// split signature: every good prime < 200 fails the strictness test
sig_all_fail := true; sig_n := 0; sig_red := 0;
for p in PrimesInInterval(3, 200) do
  if not IsDefined(chis, p) then ChiInto(~chis, g, D, p); end if;
  chi := chis[p];
  if chi eq 0 then continue; end if;
  sig_n +:= 1;
  if not IsIrreducible(chi) then sig_red +:= 1; end if;
  if not IsDefined(stricts, p) then stricts[p] := IsStrict(chi); end if;
  if stricts[p] then sig_all_fail := false; end if;
end for;

// ---------------------------------------------------------------- 4. splitness certificates
cert := Null; geom_split := false; q_split := false;
function SetCert(kind, over, detail)
  return [* <"kind", kind>, <"over", over>, <"detail", detail> *];
end function;

if strictp eq 0 then
  // (a) an extra involution over Q
  try
    A, mA := AutomorphismGroup(C);
    invols := [a : a in A | Order(a) eq 2];
    if #invols gt 1 then
      // the hyperelliptic involution acts trivially on the x-line; any other involution has a
      // genus-1 quotient (Riemann-Hurwitz: an involution of a genus-2 curve with quotient of
      // genus 0 is the hyperelliptic one)
      Amb := Ambient(C);
      extra := [];
      for a in invols do
        eqs := DefiningPolynomials(mA(a));
        if eqs[1]*Amb.3 - Amb.1*eqs[3] ne 0 then Append(~extra, a); end if;
      end for;
      if #extra gt 0 then
        eqs := DefiningPolynomials(mA(extra[1]));
        es := [SubstituteString(SubstituteString(SubstituteString(Sprint(e), "$.1", "x"), "$.2", "y"), "$.3", "z") : e in eqs];
        geom_split := true; q_split := true;
        cert := SetCert("involution", "Q", Sprintf("an involution over Q other than the hyperelliptic one, (x : y : z) -> (%o : %o : %o) on the model y^2 + h y = f; its quotient is a genus-1 curve over Q", es[1], es[2], es[3]));
        Log("split over Q: " cat cert[3][2]);
      end if;
    end if;
  catch e
    Log("AutomorphismGroup failed: " cat Sprint(e`Object));
  end try;
  // (b) Richelot isogenies over Q: collect the degenerate codomains at depth <= RichelotDepth and
  //     prefer a product over Q, then a Weil restriction of a curve that descends to Q, then any
  //     Weil restriction (geometrically split only)
  if not q_split then
    try
      J := Jacobian(SimplifiedModel(C));
      L := RichelotIsogenousSurfaces(J);
      nodes := [* <1, u> : u in L | Type(u) eq SetCart or Type(u) eq CrvEll *];
      if RichelotDepth ge 2 then
        for u in L do
          if Type(u) ne JacHyp then continue; end if;
          for v in RichelotIsogenousSurfaces(u) do
            if Type(v) eq SetCart or Type(v) eq CrvEll then Append(~nodes, <2, v>); end if;
          end for;
        end for;
      end if;
      best := 0; bestrank := 0;   // 3: product over Q, 2: descending Weil restriction, 1: Weil restriction
      bestE0 := 0;
      for nd in nodes do
        depth := nd[1]; s := nd[2];
        if Type(s) eq SetCart then r := 3; E0 := 0;
        else
          desc, E0 := DescendsToQ(s);
          r := desc select 2 else 1;
        end if;
        if r gt bestrank then bestrank := r; best := nd; bestE0 := E0; end if;
        if r eq 3 then break; end if;
      end for;
      if bestrank gt 0 then
        depth := best[1]; s := best[2];
        geom_split := true;
        how := depth eq 1 select "(2,2)-isogenous over Q to a " else "a chain of two (2,2)-isogenies over Q reaches a ";
        if bestrank eq 3 then
          q_split := true;
          cert := SetCert("richelot", "Q", how cat DescribeNode(s));
        elif bestrank eq 2 then
          q_split := true;
          cert := SetCert("richelot", "Q", how cat DescribeNode(s) cat Sprintf("; E is the base change of the elliptic curve %o over Q (conductor %o), so the surface is isogenous over Q to E0 x E0^(D), D the discriminant of the quadratic field", aInvariants(bestE0), Conductor(bestE0)));
        else
          cert := SetCert("richelot", "Q", how cat DescribeNode(s));
        end if;
        Log((q_split select "split over Q: " else "geometrically split: ") cat cert[3][2]);
      end if;
    catch e
      Log("RichelotIsogenousSurfaces over Q failed: " cat Sprint(e`Object));
    end try;
  end if;
  // (c) geometric automorphism group with more than one involution
  if not geom_split then
    try
      G := GeometricAutomorphismGroup(C);
      ninv := #[a : a in G | Order(a) eq 2];
      if ninv gt 1 then
        geom_split := true;
        cert := SetCert("geometric_involution", "Qbar", Sprintf("the geometric automorphism group has order %o (%o involutions), so the curve is bielliptic over Qbar", #G, ninv));
        Log("geometrically split: " cat cert[3][2]);
      end if;
    catch e
      Log("GeometricAutomorphismGroup failed: " cat Sprint(e`Object));
    end try;
  end if;
  // (d) Richelot over quadratic fields
  if not geom_split then
    ds := [];
    for t in Factorization(g) do
      if Degree(t[1]) eq 2 then
        d := SquarefreeFactorization(Z!Numerator(Discriminant(t[1]))*Z!Denominator(Discriminant(t[1])));
        if d ne 1 and d notin ds then Append(~ds, d); end if;
      end if;
    end for;
    d0 := SquarefreeFactorization(D);
    if d0 ne 1 and d0 notin ds then Append(~ds, d0); end if;
    for d in QuadDiscs cat [-1, 2, -2, 3, -3, 5, -5, 6, -6, 7, -7, 10, -10, 11, -11, 13, -13] do
      if d notin ds then Append(~ds, d); end if;
    end for;
    for d in ds do
      if geom_split then break; end if;
      try
        K := QuadraticField(d);
        JK := Jacobian(SimplifiedModel(HyperellipticCurve(PolynomialRing(K)!g)));
        LK := RichelotIsogenousSurfaces(JK);
        found, s := SplitNode(LK);
        if not found and RichelotDepth ge 2 then
          for u in LK do
            if Type(u) ne JacHyp then continue; end if;
            found, s := SplitNode(RichelotIsogenousSurfaces(u));
            if found then break; end if;
          end for;
        end if;
        if found then
          geom_split := true;
          cert := SetCert("richelot", Sprintf("Q(sqrt(%o))", d), Sprintf("a chain of at most %o (2,2)-isogenies over Q(sqrt(%o)) reaches a %o", RichelotDepth, d, DescribeNode(s)));
          Log("geometrically split: " cat cert[3][2]);
        end if;
      catch e
        Log(Sprintf("Richelot over Q(sqrt(%o)) failed: %o", d, e`Object));
      end try;
    end for;
  end if;
end if;
// (e) a submitted cover (checked even when a certificate was already found, as it may be over Q)
cover_checked := false; cover_ok := false; cover_msg := "";
if assigned Cover then
  cover_checked := true;
  try
    fp := Cover[1];
    if fp eq "" then
      KK := Q; PK<w> := PolynomialRing(Q); over := "Q";
    else
      PX<X> := PolynomialRing(Q);
      KK<w> := NumberField(PX!(eval fp)); over := Sprint(DefiningPolynomial(KK));
    end if;
    PKx<xx> := PolynomialRing(KK);
    if #Cover[2] eq 0 and Cover[7] ne "" then
      // the target is a curve of the Cremona database: use its Weierstrass model v^2 = u^3 + a4 u + a6
      Ew := WeierstrassModel(EllipticCurve(CremonaDatabase(), Cover[7]));
      ai := aInvariants(Ew);
      cub := [KK | 1, 0, ai[4], ai[5]];
    else
      cub := [KK | eval s : s in Cover[2]];         // e3, e2, e1, e0
    end if;
    pv := PKx![KK | eval s : s in Cover[3]];
    qv := PKx![KK | eval s : s in Cover[4]];
    hv := PKx![KK | eval s : s in Cover[5]];
    e3, e2, e1, e0 := Explode(cub);
    NN := e3*pv^3 + e2*pv^2*qv + e1*pv*qv^2 + e0*qv^3;
    disc_cubic := Discriminant(PKx![e0, e1, e2, e3]);
    wr := pv*Derivative(qv) - Derivative(pv)*qv;
    if e3 eq 0 or disc_cubic eq 0 then
      cover_msg := "the target cubic is not the equation of a smooth genus-1 curve";
    elif wr eq 0 then
      cover_msg := "the map is constant (zero Wronskian)";
    elif NN*qv ne (PKx!g)*hv^2 then
      cover_msg := "the identity (e3 p^3 + e2 p^2 q + e1 p q^2 + e0 q^3) q = g hh^2 fails for y^2 = g, g = " cat Sprint(g);
    else
      cover_ok := true;
      geom_split := true;
      if fp eq "" then q_split := true; end if;
      cert := SetCert("cover", over, Sprintf("a nonconstant map of degree %o from y^2 = g to the genus-1 curve v^2 = (%o)u^3 + (%o)u^2 + (%o)u + (%o) over %o%o, (x,y) -> (p/q, y*hh/q^2) with p = %o, q = %o, hh = %o%o", Max(Degree(pv), Degree(qv)), e3, e2, e1, e0, over, Cover[7] eq "" select "" else " (the elliptic curve " cat Cover[7] cat ")", PolyString(pv, "x"), PolyString(qv, "x"), PolyString(hv, "x"), Cover[6] eq "" select "" else "; " cat Cover[6]));
      Log((q_split select "split over Q: " else "geometrically split: ") cat cert[3][2]);
    end if;
  catch e
    cover_msg := "could not check the cover: " cat Sprint(e`Object);
  end try;
  if not cover_ok then Log("submitted cover rejected: " cat cover_msg); end if;
end if;

// the class
if strictp ne 0 then cls := "simple";
elif q_split then cls := "qsplit";
elif geom_split and qsimplep ne 0 then cls := "gsplit";
elif geom_split then cls := "split_undecided_over_Q";
else cls := "undecided"; end if;
if strictp ne 0 and (geom_split or q_split) then
  Fail("internal inconsistency: both a simplicity and a splitness certificate were found", [* *]); quit;
end if;
if q_split and qsimplep ne 0 then
  Fail("internal inconsistency: split over Q but chi_p irreducible at some good prime", [* *]); quit;
end if;
Log("class: " cat cls);

// ---------------------------------------------------------------- 5. invariants
// (the reduced minimal model and the conductor need the discriminant factored; they are computed
//  by conductor_lib.m in a separate job with its own time limit)
IC := IgusaClebschInvariants(C);
G2 := G2Invariants(C);
disc_str := Sprint(D);
nW := #Roots(g) + (IsOdd(Degree(g)) select 1 else 0);

// Q-isomorphism with existing census curves of the same G2-invariants
dups := [];
for em in ExistingModels do
  try
    Ce := HyperellipticCurve(R![Q | eval s : s in em[2]], R![Q | eval s : s in em[3]]);
    if G2Invariants(Ce) eq G2 and IsIsomorphic(C, Ce) then Append(~dups, em[1]); end if;
  catch e
    Log("isomorphism test failed for " cat em[1] cat ": " cat Sprint(e`Object));
  end try;
end for;

v1, v2, v3 := GetVersion();
ver := Sprintf("%o.%o-%o", v1, v2, v3);

// ---------------------------------------------------------------- 6. output
Emit([*
  <"ok", true>,
  <"magma_version", ver>,
  <"curve", [*
      <"f", CoeffStrings(fpol)>, <"h", CoeffStrings(hpol)>,
      <"f_str", PolyString(fpol, "x")>, <"h_str", PolyString(hpol, "x")>,
      <"g", CoeffStrings(g)>, <"g_str", PolyString(g, "x")>, <"g_scale", Sprint(dd)>,
      <"discriminant", disc_str>, <"discriminant_digits", #Sprint(Abs(D))>,
      <"igusa_clebsch", [Sprint(c) : c in IC]>, <"g2_invariants", [Sprint(c) : c in G2]>,
      <"rational_weierstrass_points", nW>
  *]>,
  <"torsion", [*
      <"invariants", inv>, <"order", ntors>, <"generators", gens>,
      <"method", tors_method>, <"generators_note", gen_note>, <"generators_submitted", #Generators>,
      <"point_count_gcd", gcdJ>, <"point_count_primes", usedp>,
      <"claimed", ClaimedGroup>, <"claim_matches", claim_ok>
  *]>,
  <"simplicity", [*
      <"geometrically_simple", strictp ne 0>, <"strict_prime", strictp>,
      <"chi", strictp ne 0 select PolyString(strictchi, "T") else "">,
      <"primes_tried_below", StrictPrimeMax>, <"admissible_primes", nadm>, <"reducible_chi", nred>
  *]>,
  <"q_simple", [* <"certified", qsimplep ne 0>, <"prime", qsimplep>,
                  <"chi", qsimplep ne 0 select PolyString(qsimplechi, "T") else ""> *]>,
  <"split", [*
      <"geometrically_split", geom_split>, <"split_over_Q", q_split>,
      <"certificate", cert>,
      <"cover_submitted", cover_checked>, <"cover_accepted", cover_ok>, <"cover_error", cover_msg>,
      <"signature", [* <"all_good_primes_fail_strictness", sig_all_fail>, <"primes_checked", sig_n>,
                       <"reducible_chi", sig_red> *]>
  *]>,
  <"class", cls>,
  <"duplicates_over_Q", dups>,
  <"cputime", RealField(6)!Cputime(T0)>
*]);
Log("VERIFY_DONE");
