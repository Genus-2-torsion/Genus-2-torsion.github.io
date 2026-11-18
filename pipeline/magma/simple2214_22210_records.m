SetQuitOnError(true); SetColumns(0);
Q:=Rationals(); Z:=Integers(); R<x>:=PolynomialRing(Q);
h:=x^2+x; f14:=R![-35,165,-183,-95,123,69,9]; C14:=HyperellipticCurve(f14,h);
f10:=R![-7791,-7511,3544,3144,-402,-246,36]; C10:=HyperellipticCurve(f10,h);
function Chi(C,p)
 assert Valuation(Discriminant(C),p) eq 0; Cp:=ChangeRing(C,GF(p));
 return R!Reverse(Coefficients(LPolynomial(Cp:Al:="Naive")));
end function;
procedure CheckStrict(C,p,expected)
 chi:=Chi(C,p);assert chi eq expected and IsIrreducible(chi);
 K<t>:=FunctionField(Q); S<T>:=PolynomialRing(K);
 ratio:=R!Numerator(Resultant(Evaluate(chi,T),Evaluate(chi,t*T)));
 ratio:=ExactQuotient(ratio,(x-1)^4);assert Degree(ratio) eq 12;
 assert &and[GCD(ratio,R!CyclotomicPolynomial(i)) eq 1:i in [1..42]];
 print "STRICT",p,chi;
end procedure;
// Three rational finite branch points and infinity give three independent
// 2-torsion classes. The contact at P=(1,81/40) gives 7-torsion.
k:=1-7/2*x+7/40*x^2+87/20*x^3;
f:=(x-16/25)*(x-7/16)*(x+8)*(x^2+5*x-15/4);
assert x^2*f eq k^2+(x-1)^7; C:=HyperellipticCurve(f); assert IsIsomorphic(C,C14);
J:=Jacobian(C); D:=elt<J|x-1,81/40>; assert D ne J!0 and 7*D eq J!0;
T:=[elt<J|x-r,0>:r in [16/25,7/16,-8]];
assert #{i*T[1]+j*T[2]+l*T[3]:i,j,l in [0,1]} eq 8;
assert &and[2*A eq J!0:A in T]; assert Evaluate(Chi(C14,7),1) eq 56;
CheckStrict(C14,17,x^4+4*x^3+30*x^2+68*x+289);
// Elkies' Clebsch-Klein surface. The infinity difference gives 5-torsion.
r:=[2,4,10,-7,-9];assert &+r eq 0 and &+[a^3:a in r] eq 0;
f:=x*&*[x-a^2:a in r];C:=HyperellipticCurve(f);assert IsIsomorphic(C,C10);
J:=Jacobian(C);inf:=PointsAtInfinity(C);D:=J![inf[1],inf[2]];
assert D ne J!0 and 5*D eq J!0;
T:=[elt<J|x*(x-r[i]^2),0>:i in [1..4]];assert &and[2*A eq J!0:A in T];
assert #{i*T[1]+j*T[2]+k0*T[3]+l*T[4]:i,j,k0,l in [0,1]} eq 16;
counts:=[Z!Evaluate(Chi(C10,p),1):p in [23,29,53,61]];
assert counts eq [640,960,2400,4240] and GCD(counts) eq 80;
CheckStrict(C10,61,x^4+8*x^3+22*x^2+488*x+3721);
// Reduction at 2 in the charts where x or 1/x is finite. Both Jacobians
// have semistable toric dimension 2, as detailed in the accompanying note.
P<u,v>:=PolynomialRing(GF(2),2);
red:=func<f| &+[GF(2)!Coefficient(f,i)*u^i:i in [0..6]]>;
assert v^2+(u^2+u)*v+red(f14) eq (v+1+u^2+u^3)*(v+1+u+u^3);
assert red(f10) eq u+1;
E:=v^2+(u^2+u)*v+red(f10);
assert Evaluate(E,[u,v+1]) eq v^2+(u^2+u)*v+u^2;
Einf:=v^2+(u+u^2)*v+u^5+u^6;
assert Einf eq v^2+(u^2+u)*v+&+[GF(2)!Coefficient(f10,i)*u^(6-i):i in [0..6]];
// The two charts have an ordinary double point at u=0. The only other
// possible affine singular point is (1,0) where the first derivative is 1.
assert Evaluate(Derivative(E,1),[GF(2)!1,GF(2)!0]) eq 1;
assert Evaluate(E,[GF(2)!0,GF(2)!1]) eq 0 and Evaluate(Einf,[GF(2)!0,GF(2)!0]) eq 0;
assert Evaluate(E,[u,v+1])-(v^2+u*v+u^2) eq u^2*v;
assert Einf-(v^2+u*v) eq u^2*v+u^5+u^6;
// If the latter special fiber were geometrically reducible, the monic
// quadratic in v would have a root g in the integral closure of k[u] in
// k(u), hence a polynomial. Comparing leading terms gives deg g <= 2.
// Its coefficients satisfy the exact inconsistent system below.
A<a,b,c>:=PolynomialRing(GF(2),3);
I:=ideal<A|a^2+1,a+1,b^2+b+a,b+c,c^2+c>;
assert 1 in I;
assert Conductor(C14) eq 152100 and Conductor(C10) eq 22406283900;
print "CONDUCTOR_FACTORS",Factorization(152100),Factorization(22406283900);
print "SIMPLE_2214_AND_22210_CERTIFIED";quit;
