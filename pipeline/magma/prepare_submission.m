/* prepare_submission.m -- run this on YOUR machine before submitting a curve to the census
   https://genus-2-torsion.github.io/ : it computes what the verifier on Mordell would otherwise
   have to compute, and prints the submission as JSON.

   Usage (Magma V2.28-9 or later):
     magma -b f:="x^5 - 8*x^4 + 16*x^3 - x" h:="x" pipeline/magma/prepare_submission.m
   or edit the two lines below and run  magma -b pipeline/magma/prepare_submission.m
   f and h are polynomials in x with rational coefficients: the curve is y^2 + h(x) y = f(x)
   (deg f <= 6, deg h <= 3; h may be "0").

   What it prints:
     * J(Q)_tors (Magma's TorsionSubgroup) and its generators in Mumford representation
       <a(x), b(x), d> on the even model y^2 = 4f + h^2 -- paste them into the submission: if
       the order of the group they generate equals the gcd of #J(F_p) over good primes (also
       printed), the verifier certifies the torsion from the generators without any search;
     * the gcd bound, and whether it is sharp;
     * a strict prime if one is found below 500 (then the Jacobian is geometrically simple), else
       whether the "split signature" holds (every good prime below 200 fails the strictness test);
     * the conductor, when the discriminant has at most 60 digits (Magma needs to factor it);
     * the submission JSON (paste it into the "bulk submission" format, or its fields into the
       GitHub issue form).
   Frobenius polynomials are computed with LPolynomial(C : Al := "Naive"): Magma's default
   algorithm returns wrong L-polynomials for some genus 2 curves over F_3, F_5, F_7.
*/
if not assigned f then f := "x^5 - 8*x^4 + 16*x^3 - x"; end if;
if not assigned h then h := "x"; end if;
SetColumns(0);
Q := Rationals(); Z := Integers();
R<x> := PolynomialRing(Q);
fpol := R!(eval f); hpol := R!(eval h);
C := HyperellipticCurve(fpol, hpol);
error if Genus(C) ne 2, "the curve does not have genus 2";
g0 := 4*fpol + hpol^2;
dd := LCM([Z | Denominator(c) : c in Coefficients(g0)]);
g := R!(dd^2*g0);
D := Z!Discriminant(g);
printf "curve: y^2 + (%o) y = %o\n", hpol, fpol;
printf "even model (Magma's SimplifiedModel): y^2 = %o\n", g0;

function LPoly(gp)  // reliable L-polynomial (see the header)
  return LPolynomial(HyperellipticCurve(gp) : Al := "Naive");
end function;
function IsStrict(chi)
  if Degree(chi) ne 4 or not IsIrreducible(chi) then return false; end if;
  K := NumberField(chi); pi := K.1;
  for n in [2..12] do
    if Degree(MinimalPolynomial(pi^n)) lt 4 then return false; end if;
  end for;
  return true;
end function;

// torsion
J := Jacobian(SimplifiedModel(C));
t0 := Cputime();
T, mT := TorsionSubgroup(J);
inv := Invariants(T);
printf "J(Q)_tors = %o (order %o; %o s)\n", inv, #T, RealField(4)!Cputime(t0);
gens := [];
for i in [1..Ngens(T)] do
  P := mT(T.i);
  Append(~gens, <Order(T.i), Sprint(P[1]), Sprint(P[2]), P[3]>);
  printf "  generator of order %o: <%o, %o, %o>\n", Order(T.i), P[1], P[2], P[3];
end for;

// gcd bound
gcdJ := 0; used := [];
for p in PrimesInInterval(3, 200) do
  if D mod p eq 0 or (Z!LeadingCoefficient(g)) mod p eq 0 then continue; end if;
  gp := PolynomialRing(GF(p))!g;
  if Degree(gp) lt 5 or not IsSquarefree(gp) then continue; end if;
  gcdJ := GCD(gcdJ, Z!Evaluate(LPoly(gp), 1)); Append(~used, p);
  if #used ge 25 then break; end if;
end for;
printf "gcd of #J(F_p) over p = %o: %o -- %o\n", used, gcdJ,
  gcdJ eq #T select "sharp: the verifier will certify the torsion from the generators alone"
                else "not sharp: the verifier will run TorsionSubgroup itself";

// class
strictp := 0; strictchi := 0; allfail := true; n200 := 0;
for p in PrimesInInterval(3, 500) do
  if D mod p eq 0 or (Z!LeadingCoefficient(g)) mod p eq 0 then continue; end if;
  gp := PolynomialRing(GF(p))!g;
  if Degree(gp) lt 5 or not IsSquarefree(gp) then continue; end if;
  chi := R!Reverse(Coefficients(LPoly(gp)));
  st := IsStrict(chi);
  if p lt 200 then n200 +:= 1; if st then allfail := false; end if; end if;
  if st then strictp := p; strictchi := chi; break; end if;
end for;
cls := "";
if strictp ne 0 then
  cls := "simple";
  printf "geometrically simple: strict prime %o, chi = %o\n", strictp, strictchi;
elif allfail and n200 gt 0 then
  cls := "split";
  printf "no strict prime below 500 and every good prime below 200 fails the strictness test: the Jacobian is probably geometrically split\n";
  printf "  (if it splits only through an odd-degree isogeny or over a field of degree > 2, supply a map to an elliptic curve -- see the submit page)\n";
else
  printf "no strict prime below 500 (inconclusive)\n";
end if;

// conductor
cond := "";
if #Sprint(Abs(D)) le 60 then
  try
    Cm := ReducedMinimalWeierstrassModel(C);
    N := Conductor(Cm); cond := Sprint(N);
    printf "conductor %o = %o\n", N, Factorization(N);
  catch e
    printf "conductor not computed: %o\n", e`Object;
  end try;
else
  printf "conductor not computed (discriminant has %o digits)\n", #Sprint(Abs(D));
end if;

// JSON
function JS(s) return "\"" cat s cat "\""; end function;
fs := [Sprint(c) : c in Coefficients(fpol)]; hs := hpol eq 0 select ["0"] else [Sprint(c) : c in Coefficients(hpol)];
lst := func<L | "[" cat (#L eq 0 select "" else &cat[JS(L[i]) cat (i lt #L select ", " else "") : i in [1..#L]]) cat "]">;
gj := "[" cat (#gens eq 0 select "" else &cat[ "[" cat JS(gens[i][2]) cat ", " cat JS(gens[i][3]) cat ", " cat Sprint(gens[i][4]) cat "]" cat (i lt #gens select ", " else "") : i in [1..#gens]]) cat "]";
printf "\n---- submission JSON (fill in the last five fields) ----\n";
printf "{\n  \"schema\": \"genus-2-torsion/submission/1\",\n  \"coeffs\": [%o, %o],\n  \"group\": %o,\n  \"class\": \"%o\",\n  \"generators\": %o,\n", lst(fs), lst(hs), Sprint(inv), cls, gj;
printf "  \"conductor\": \"%o\",\n  \"submitter\": \"\", \"affiliation\": \"\", \"github\": \"\",\n  \"discoverer\": \"\", \"year\": null, \"reference\": \"\", \"notes\": \"\"\n}\n", cond;
printf "---- generators field for the issue form ----\n%o\n", gj;
quit;
