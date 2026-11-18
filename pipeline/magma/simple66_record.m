SetQuitOnError(true); SetColumns(0);
Q:=Rationals(); Z:=Integers(); R<x>:=PolynomialRing(Q);
f:=-8*x^5+256*x^4-479*x^3+130*x^2+1000*x;
h:=x^2+x; C:=HyperellipticCurve(f,h); J:=Jacobian(HyperellipticCurve(4*f+h^2));
u:=Q!3; r:=4/(u^3+3*u+2);
q:=3*x^2+6*x+4*r-1;
b:=8*(r-2)*x^3+3*r*(3*r-8)*x^2+6*r^2*(r-3)*x+r^2*(r^2-6*r+1);
assert IsIsomorphic(C,HyperellipticCurve(-q*b));
// Exact rational divisors on this model and an independent size upper bound.
A:=elt<J|x^2-x,60*x>; B:=elt<J|x^2-11*x-50,270*x+450>;
assert 6*A eq J!0 and 6*B eq J!0;
assert #{i*A+j*B:i,j in [0..5]} eq 36;
print "MUMFORD_GENERATORS", A,B;
function Chi(p)
 Cbar:=ChangeRing(C,GF(p)); assert IsNonsingular(Cbar);
 return R!Reverse(Coefficients(LPolynomial(Cbar:Al:="Naive")));
end function;
assert Evaluate(Chi(7),1) eq 36;
chi:=Chi(61); assert chi eq x^4-4*x^3+54*x^2-244*x+3721 and IsIrreducible(chi);
print "JACOBIAN_ORDER_7",Evaluate(Chi(7),1),"FROBENIUS_61",chi;
// The 12 distinct-index ratios are roots of a polynomial of degree 12.
// A root of unity therefore has phi(n)<=12; these n are all <=42.
K<t>:=FunctionField(Q); S<T>:=PolynomialRing(K);
f1:=Evaluate(chi,T); f2:=Evaluate(chi,t*T);
ratio:=R!Numerator(Resultant(f1,f2)); ratio:=ExactQuotient(ratio,(x-1)^4);
assert Degree(ratio) eq 12;
assert &and[GCD(ratio,R!CyclotomicPolynomial(i)) eq 1:i in [1..42]];
print "FROBENIUS_61_NO_ROOT_OF_UNITY_RATIOS";
// The special fiber in P(1,1,3) consists of graphs of two cubic sections.
H<X,Z0,Y>:=PolynomialRing(GF(2),3);
fb:=&+[GF(2)!Coefficient(f,i)*X^i*Z0^(6-i):i in [0..6]];
hb:=&+[GF(2)!Coefficient(h,i)*X^i*Z0^(3-i):i in [0..3]];
assert Y^2+hb*Y-fb eq (Y+X*Z0^2)*(Y+X^2*Z0);
assert X*Z0^2-X^2*Z0 eq X*Z0*(X+Z0);
// Three transverse rational nodes, no other intersections, so the semistable
// conductor exponent is the first Betti number 3-2+1=2.
N:=Conductor(C); assert N eq 11337300; assert (N div 4) eq 2834325;
print "ODD_CONDUCTOR_FACTORIZATION",Factorization(Z!(N div 4));
print "SIMPLE_66_CONDUCTOR_11337300_CERTIFIED";
quit;
