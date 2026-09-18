SetColumns(0);
SetMemoryLimit(8*10^9);
R<x> := PolynomialRing(Rationals());
function Strict(c)
  if not IsIrreducible(c) then return false; end if;
  K<a>:=NumberField(c);
  return &and[Degree(MinimalPolynomial(a^j)) eq 4 : j in [2..12]];
end function;
procedure Check(tag, f, h, want, p)
  if Discriminant(4*f+h^2) eq 0 then
    printf "FAMILY_SINGULAR %o\n",tag;
    return;
  end if;
  gg:=4*f+h^2; dd:=LCM([Denominator(c):c in Coefficients(gg)]);
  C:=HyperellipticCurve(dd^2*gg);
  got:=Invariants(TorsionSubgroup(Jacobian(SimplifiedModel(C))));
  str:="not requested";
  if p ne 0 then
    g:=4*f+h^2;
    if (Numerator(Discriminant(g))*Denominator(Discriminant(g))) mod p eq 0 then
      str:="BAD REDUCTION";
    else
      chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(PolynomialRing(GF(p))!g):Al:="Naive")));
      str:=Sprint(Strict(chi));
    end if;
  end if;
  printf "FAMILY %o expected=%o got=%o strict(%o)=%o\n",tag,want,got,p,str;
end procedure;
for t in [1,2,3,4,5] do
  Check(Sprintf("[2] t=%o",t),x*(x^5+x+t),R!0,[2],t eq 2 select 17 else 0);
end for;
Check("[2,2] t=3",x*(x-1)*(x-2)*(x^3+x+3),R!0,[2,2],11);
Check("[2,2,2] t=3",x*(x-1)*(x-2)*(x^2+x+3),R!0,[2,2,2],29);
Check("[2,2,2,2] t=5",x*(x-1)*(x-2)*(x-3)*(x-5),R!0,[2,2,2,2],11);
Check("[3] t=2",(x^3+x+1)^2-2*(x^2+1)^3,R!0,[3],13);
aa:=[Rationals()|1,2,3,-3,5,1/2,7,-5]; pp:=[17,17,13,11,23,11,19,11];
for i in [1..#aa] do
  a:=aa[i]; Check(Sprintf("[10] a=%o",a),(1+a*x+x^2)^2-(a+2)^2*x^5,R!0,[10],pp[i]);
end for;
aa:=[Rationals()|1,2,3,-1,-2,1/2,5,-3]; pp:=[19,17,19,19,17,11,11,19];
for i in [1..#aa] do
  s:=aa[i]; f:=x*((x^2+3*x+1+2*s^2)^2-8*s^2*(x+1)^2);
  Check(Sprintf("[4] s=%o",s),f,R!0,[4],pp[i]);
end for;
aa:=[Rationals()|3,5,7,-2,1/2,-4,11,13]; pp:=[17,11,11,13,29,11,13,37];
for i in [1..#aa] do
  t:=aa[i]; Check(Sprintf("[2,4] t=%o",t),x*(x^2+(2*t+1)*x+t^2)*(x^2+3*x+1),R!0,[2,4],pp[i]);
end for;
aa:=[Rationals()|1,3,-2,5,1/2,-4,7]; pp:=[13,19,17,13,29,17,17];
for i in [1..#aa] do
  t:=aa[i]; Check(Sprintf("[2,2,4] t=%o",t),x*(x+4)*(x+9)*(x^2+(2*t+1)*x+t^2),R!0,[2,2,4],pp[i]);
end for;
Check("[2,2,2,4] t=5",x*(x+1)*(x+4)*(x+9)*(x+25),R!0,[2,2,2,4],59);
t:=Rationals()!2; a:=-4*t^2*(t+1)/(t^2+t+1)^2; b:=-t/(t+1);c:=1;d:=t;
Check("[2,2,2,8] t=2",x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2),R!0,[2,2,2,8],31);
aa:=[<25,-26,-15>,<7,2,1>,<3,-2,5>,<11,4,-3>,<5,7,2>]; pp:=[23,47,41,23,67];
for i in [1..#aa] do
  s,m,n:=Explode(aa[i]);
  f:=x*(x+2*s^2-s*n)*(x+2*s^2+s*m-2*s*n-m*n)*(x+2*s^2+s*m-s*n-m*n)*(2*x-m*n)*(2*x+4*s^2-4*s*n-m*n);
  Check(Sprintf("[2,2,2,6] %o",aa[i]),f,R!0,[2,2,2,6],pp[i]);
end for;
for t in [2,3,4] do
  f:=x^6+2*x^5+(2*t+3)*x^4+2*x^3+(t^2+1)*x^2+2*t*(1-t)*x+t^2;
  Check(Sprintf("[11] Flynn t=%o",t),f,R!0,[11],0);
end for;
print "FAMILY_AUDIT_COMPLETE";
quit;
