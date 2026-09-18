SetColumns(0);R<x>:=PolynomialRing(Rationals());
function Psi(P)
 xx:=P[1]/P[3];yy:=P[2]/P[3];
 w:=(3*xx+yy-6)/(yy+18);t:=-(4*xx^2-7*xx+12*yy+16)/((xx-8)*(3*xx+2*yy+12));
 return [3*t*(1-t),3*t-1,2*w*(3*t^2+1),3*t^2-3*t+2,6*t^2-3*t+1];
end function;
for p in [29,37,71] do
 Ep:=EllipticCurve([GF(p)|0,0,0,-21,-20]);Qp:=Ep![-3,4];
 assert 560*Qp eq Ep!0;
 im:=Psi(2*Qp);
 assert [a/im[1]:a in im] eq [GF(p)|a/120:a in [120,143,266,218,241]];
 printf "PROGRESSION_MOD_%o_VERIFIED order(Q)=%o im2=%o\n",p,Order(Qp),im;
end for;
for t in [1,2,3] do
 f:=x*(x^4+x+t);C:=HyperellipticCurve(f);T:=TorsionSubgroup(Jacobian(C));
 printf "REPLACEMENT_2_FAMILY t=%o torsion=%o G2=%o\n",t,Invariants(T),G2Invariants(C);
end for;
quit;
