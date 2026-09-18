SetColumns(0);R<x>:=PolynomialRing(Rationals());C:=HyperellipticCurve(R![26,41,-44,-30,25,-4],R![1,1]);
A:=AutomorphismGroup(C);print "AUT_Q",GroupName(A),#A;
L:=RichelotIsogenousSurfaces(Jacobian(SimplifiedModel(C)));
for obj in L do
 if Type(obj) eq CrvEll then
  K:=BaseRing(obj); print "ELLIPTIC_FACTOR",DefiningPolynomial(K),aInvariants(obj),"j",jInvariant(obj);
 elif Type(obj) eq SetCart then print "Q_PRODUCT",Conductor(obj[1]),Conductor(obj[2]);end if;
end for;
quit;
