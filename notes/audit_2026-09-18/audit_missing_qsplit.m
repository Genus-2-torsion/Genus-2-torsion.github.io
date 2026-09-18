SetColumns(0);R<x>:=PolynomialRing(Rationals());
for bc in [<2,4>,<2,5>,<2,7>,<3,5>] do
 b,c:=Explode(bc);
 f:=(x^2-1)*(x^2-b^2)*(x^2-c^2);
 C:=HyperellipticCurve(f);T:=TorsionSubgroup(Jacobian(C));
 E1:=EllipticCurve(HyperellipticCurve((x-1)*(x-b^2)*(x-c^2)));
 D:=HyperellipticCurve(x*(x-1)*(x-b^2)*(x-c^2));E2:=EllipticCurve(D,D![0,0,1]);
 printf "QSPLIT_EXAMPLE b,c=%o,%o f=%o torsion=%o N(E1)=%o N(E2)=%o N(J)=%o\n",b,c,f,Invariants(T),Conductor(E1),Conductor(E2),Conductor(E1)*Conductor(E2);
end for;
quit;
