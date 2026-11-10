/* conductor_by_prime.m -- as conductor_lib.m (same job variables, same JSON output), but the conductor is
   computed prime by prime: N = prod p^Conductor(Cm, p) over the primes dividing the minimal discriminant,
   and the exponent at 2 on an isomorphic model x -> (aX + b)/(cX + d) given by the job variable
   Model2 := [a, b, c, d] (default [0, 1, 1, 0], i.e. x -> 1/X).  Magma's regular model at 2, which its
   conductor exponent at 2 needs, can take hours on one model and seconds on an isomorphic one; the
   exponent itself does not depend on the model.  The model is checked to be isomorphic to Cm first.
   Used by hand through `fill_conductors.py --reuse` for curves where conductor_lib.m times out at 2
   (2.2.16.b: no result in 4 h on the reduced minimal model, 2.4 s on x -> 1/X).  The same caveat as
   in conductor_lib.m: the exponent at 2 uses Ogg's formula, without a correctness guarantee when
   v_2(disc) >= 12; odd exponents are rigorous. */
if not assigned Model2 then Model2 := [0, 1, 1, 0]; end if;
if not assigned MemGB then MemGB := 8; end if;
if Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB * 10^9);
T0 := Cputime();
procedure Log(s) Write(LogFile, Sprintf("[%o s] %o", RealField(6)!Cputime(T0), s)); end procedure;
function JStr(s)
  t := "";
  for i in [1..#s] do
    c := s[i];
    if c eq "\"" then t cat:= "\\\""; elif c eq "\\" then t cat:= "\\\\"; elif c eq "\n" then t cat:= "\\n"; else t cat:= c; end if;
  end for;
  return "\"" cat t cat "\"";
end function;
function JList(v) return "[" cat (#v eq 0 select "" else &cat[JStr(v[i]) cat (i lt #v select "," else "") : i in [1..#v]]) cat "]"; end function;
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
Q := Rationals(); Z := Integers();
R<x> := PolynomialRing(Q);
fpol := R![Q | eval s : s in fcoeffs];
hpol := R![Q | eval s : s in hcoeffs];
C := HyperellipticCurve(fpol, hpol);
Cm := ReducedMinimalWeierstrassModel(C);
fm, hm := HyperellipticPolynomials(Cm);
Dm := Z!Discriminant(Cm);
Log(Sprintf("minimal model y^2 + (%o) y = %o, discriminant %o", hm, fm, Dm));
dfact := FactoredString(Dm);
Log("discriminant factored: " cat dfact);
// the minimal model first (ok: false), so that it is kept if the conductor fails or times out
mm := "\"minimal_model\":{\"f\":" cat JList([Sprint(c) : c in Coefficients(fm)]) cat ",\"h\":" cat JList([Sprint(c) : c in Coefficients(hm)]) cat "}," cat
  "\"minimal_discriminant\":" cat JStr(Sprint(Dm)) cat ",\"minimal_discriminant_factored\":" cat JStr(dfact);
Write(OutFile, "{\"ok\":false," cat mm cat "}" : Overwrite := true);
Log("MINIMAL_MODEL_DONE");
N := 1;
for t in Factorization(Abs(Dm)) do
  p := t[1];
  if p eq 2 then
    C2 := Transformation(Cm, Model2);
    error if not IsIsomorphic(C2, Cm), "Model2 does not give a model isomorphic to the curve";
    e := Conductor(C2, 2);
    Log(Sprintf("p = 2: exponent %o, on the model x -> (%o X + %o)/(%o X + %o)", e, Model2[1], Model2[2], Model2[3], Model2[4]));
  else
    e := Conductor(Cm, p);
    Log(Sprintf("p = %o: exponent %o", p, e));
  end if;
  N *:= p^e;
end for;
note := (Valuation(Dm, 2) ge 12 select "the exponent of 2 in the conductor was computed with Ogg's formula, for which Magma gives no correctness guarantee when v_2(disc) >= 12; odd exponents are rigorous" else "")
  cat Sprintf("%ocomputed prime by prime (conductor_by_prime.m), the exponent at 2 on the isomorphic model x -> (%o X + %o)/(%o X + %o)", Valuation(Dm, 2) ge 12 select "; " else "", Model2[1], Model2[2], Model2[3], Model2[4]);
Log(Sprintf("conductor %o = %o", N, FactoredString(N)));
Write(OutFile, "{\"ok\":true," cat mm cat
  ",\"conductor\":" cat JStr(Sprint(N)) cat ",\"conductor_factored\":" cat JStr(FactoredString(N)) cat ",\"conductor_note\":" cat JStr(note) cat
  ",\"cputime\":" cat Sprint(RealField(6)!Cputime(T0)) cat "}" : Overwrite := true);
Log("CONDUCTOR_DONE");
