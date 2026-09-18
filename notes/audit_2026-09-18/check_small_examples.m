SetColumns(0);R<x>:=PolynomialRing(Rationals());
for which in [1,2] do
 f:=which eq 1 select x*(x^4+x+1) else (x^2-1)*(x^2-4)*(x^2-49);gcd:=0;
 for p in PrimesInInterval(3,43) do
  if (Integers()!Discriminant(f)) mod p eq 0 then continue;end if;
  c:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(PolynomialRing(GF(p))!f):Al:="Naive")));
  jp:=Integers()!Evaluate(c,1);gcd:=GCD(gcd,jp);
  printf "SMALL_EXAMPLE which=%o p=%o chi=%o #J=%o gcd=%o\n",which,p,c,jp,gcd;
 end for;
end for;
quit;
