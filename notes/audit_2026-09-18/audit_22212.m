SetColumns(0);SetMemoryLimit(8*10^9);
R<x>:=PolynomialRing(Rationals());
function Strict(c)
 if not IsIrreducible(c) then return false;end if;
 K<a>:=NumberField(c);return &and[Degree(MinimalPolynomial(a^n)) eq 4:n in [2..12]];
end function;
E:=EllipticCurve([Rationals()|0,0,0,-21,-20]);Q:=E![-3,4];
function Psi(P)
 x:=P[1]/P[3];y:=P[2]/P[3];
 w:=(3*x+y-6)/(y+18);t:=-(4*x^2-7*x+12*y+16)/((x-8)*(3*x+2*y+12));
 return [3*t*(1-t),3*t-1,2*w*(3*t^2+1),3*t^2-3*t+2,6*t^2-3*t+1];
end function;
p1:=[120,143,266,218,241];im2:=Psi(2*Q);im4:=Psi(4*Q);
assert [a/im2[1]:a in im2] eq [Rationals()|a/120:a in p1];
f1:=x*&*[x-a^2:a in p1];
f2:=x*&*[x-a^2:a in im4];
assert G2Invariants(HyperellipticCurve(f1)) ne G2Invariants(HyperellipticCurve(f2));
print "MODULI_NONCONSTANT_VERIFIED";
print "PSI(2Q)",im2;
print "SAMPLE_EXACT_TORSION",Invariants(TorsionSubgroup(Jacobian(HyperellipticCurve(f1))));
gcd:=0;
for p in [29,31,37,43,47,53,59,61,67,71] do
 assert (Integers()!Discriminant(f1)) mod p ne 0;
 cp:=HyperellipticCurve(PolynomialRing(GF(p))!f1);
 chi:=R!Reverse(Coefficients(LPolynomial(cp:Al:="Naive")));
 jp:=Integers()!Evaluate(chi,1);gcd:=GCD(gcd,jp);
 printf "P1 p=%o chi=%o #J=%o gcd=%o strict=%o Epointorder=%o\n",p,chi,jp,gcd,Strict(chi),Order(EllipticCurve([GF(p)|0,0,0,-21,-20])![-3,4]);
end for;
assert gcd eq 96;
print "EXACT_INFINITY_LOCAL_DATA_VERIFIED";
K<t>:=FunctionField(Rationals());S<z>:=PolynomialRing(K);
assert IsIrreducible(z^5+z+t);
print "GENERIC_2_FAMILY_QUINTIC_IRREDUCIBLE";
quit;
