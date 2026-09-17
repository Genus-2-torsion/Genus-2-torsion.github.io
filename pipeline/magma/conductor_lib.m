/* conductor_lib.m -- reduced minimal Weierstrass model, factored discriminant and conductor of a
   genus-2 curve y^2 + h(x) y = f(x) over Q.  Run by pipeline/verify.py as a separate Magma job
   under a time limit, because the conductor needs the discriminant factored.  The job defines
   fcoeffs, hcoeffs (ascending coefficient strings), LogFile, OutFile, MemGB (optional).
   Magma's conductor exponent at 2 uses Ogg's formula when v_2(disc) >= 12 and prints a warning
   that it gives no correctness guarantee; the certificate records this as conductor_note.  Odd
   exponents are rigorous (Liu / Ogg-Saito). */
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
N := Conductor(Cm);
note := Valuation(Dm, 2) ge 12 select "the exponent of 2 in the conductor was computed with Ogg's formula, for which Magma gives no correctness guarantee when v_2(disc) >= 12; odd exponents are rigorous" else "";
Log(Sprintf("conductor %o = %o", N, FactoredString(N)));
Write(OutFile, "{\"ok\":true,\"minimal_model\":{\"f\":" cat JList([Sprint(c) : c in Coefficients(fm)]) cat ",\"h\":" cat JList([Sprint(c) : c in Coefficients(hm)]) cat "}," cat
  "\"minimal_discriminant\":" cat JStr(Sprint(Dm)) cat ",\"minimal_discriminant_factored\":" cat JStr(dfact) cat
  ",\"conductor\":" cat JStr(Sprint(N)) cat ",\"conductor_factored\":" cat JStr(FactoredString(N)) cat ",\"conductor_note\":" cat JStr(note) cat
  ",\"cputime\":" cat Sprint(RealField(6)!Cputime(T0)) cat "}" : Overwrite := true);
Log("CONDUCTOR_DONE");
