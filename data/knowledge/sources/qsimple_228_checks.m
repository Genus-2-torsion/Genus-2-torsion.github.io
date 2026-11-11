// Exact arithmetic for the proposed [2,28] Q-simple geometrically split exclusion.
// Mathematical representation and polarization arguments are in proof.md.
// No analytic-rank or BSD assumption is used.
SetColumns(0);
SetSeed(1);
Q:=Rationals(); Z:=Integers();
X:=SmallModularCurve(21); F<x,y>:=FunctionField(X);
assert aInvariants(X) eq [1,0,0,-4,-1];
assert Genus(X) eq 1;
j:=jFunction(X,21);
assert Degree(j) eq 32;
cusps:=[Cusp(X,21,d):d in [1,3,7,21]];
assert #Set(cusps) eq 4;
poles:=Sort([-Valuation(j,Place(c)):c in cusps]);
assert poles eq [1,3,7,21];
M:=CuspidalSubspace(ModularSymbols(21,2,+1));
assert Dimension(M) eq 1;
for p in [p:p in PrimesUpTo(43)|p notin {3,7}] do
 assert #Points(ChangeRing(X,GF(p))) eq p+1-Z!Trace(HeckeOperator(M,p));
end for;
w:=AtkinLehnerInvolution(X,21,3);
assert w*w eq IdentityMap(X);
Phi3:=ClassicalModularPolynomial(3);
assert Evaluate(Phi3,[j,Pullback(w,j)]) eq 0;
assert Evaluate(Phi3,[j,Pullback(AtkinLehnerInvolution(X,21,7),j)]) ne 0;
assert Evaluate(Phi3,[j,Pullback(AtkinLehnerInvolution(X,21,21),j)]) ne 0;
print "MODULAR_MODEL_JMAP_AND_W3_VERIFIED";

rb:=TwoPowerIsogenyDescentRankBound(X,X![-2,1]:Cutoff:=1,MaxSteps:=6);
assert rb eq 0;
TX,mTX:=TorsionSubgroup(X); assert Invariants(TX) eq [2,4];
jvalues:={Q|};
for a in TX do
 P:=mTX(a); jj:=Evaluate(j,P);
 if Type(jj) ne Infty then Include(~jvalues,Q!jj); end if;
end for;
assert jvalues eq {Q|-189613868625/128,3375/2,-1159088625/2097152,-140625/8};
for jj in jvalues do
 E:=EllipticCurveFromjInvariant(jj); D2:=DivisionPolynomial(E,2);
 assert IsIrreducible(D2) and not IsSquare(Discriminant(D2));
end for;
print "RATIONAL_X021_CLASSES_ALL_S3_AT_2",jvalues;

A,m:=AlgorithmicFunctionField(F);
t:=(y-x-3)/(x+1); H:=t*(x+2);
assert IsSquare(m((j-1728)/H));
assert Pullback(w,t) eq t;
q:=t/(x+2);
assert Pullback(w,H)/H eq q^2;
assert q*Pullback(w,q) eq 1;
assert H*(1+q)^2 eq t*(t^2+5*t+7);
assert H*(1-q)^2 eq t*(t^2+t+7);
assert x^2-(t^2+3*t+3)*x-t^2-6*t-10 eq 0;
R<z>:=PolynomialRing(A); L:=ext<A|z^2-m(j-1728)>;
assert Genus(L) eq 3;
print "DISCRIMINANT_COVER_AND_QUOTIENT_IDENTITIES_DONE";

for b in [1,5] do
 E:=EllipticCurve([Q|0,b,0,7,0]);
 rb:=TwoPowerIsogenyDescentRankBound(E,E![0,0]:Cutoff:=1,MaxSteps:=6);
 assert rb eq 0;
 T,tm:=TorsionSubgroup(E);
 pts:=[tm(a):a in T];
 ts:={Q|P[1]/P[3]:P in pts|P[3] ne 0};
 assert #{P:P in pts|P[3] eq 0} eq 1;
 if b eq 1 then
  assert Invariants(T) eq [6] and ts eq {Q|0,1,7};
 else
  assert Invariants(T) eq [2] and ts eq {Q|0};
 end if;
 print "QUOTIENT_ALGEBRAIC_RANK_ZERO",b,"TORSION",Invariants(T),"FINITE_T",ts;
end for;

// The t=0 and t=infinity fibers of X0(21) are exactly its four cusps.
assert Set([Evaluate(t,P):P in cusps]) eq {0,Infinity()};

K<a>:=QuadraticField(13);
for tt in [1,7] do
 xx:=tt eq 1 select (7+3*a)/2 else (73+21*a)/2;
 yy:=tt*(xx+1)+xx+3;
 P:=X(K)![xx,yy,1]; jj:=Evaluate(j,P);
 expected:=tt eq 1 select 6565897584689283072+1821052339310297088*a
                         else 75497472+20938752*a;
 assert jj eq expected;
 E:=EllipticCurveFromjInvariant(jj);
 k:=GF(29); roots:=Roots(Polynomial([k|-13,0,1])); assert #roots eq 2;
 for rr in roots do
  mp:=hom<K->k|rr[1]>;
  Ep:=EllipticCurve([mp(c):c in aInvariants(E)]);
  assert Discriminant(Ep) ne 0;
  ap:=Z!(30-#Ep);
  assert Abs(ap) eq 6;
  assert #Roots(DivisionPolynomial(Ep,2)) eq 3;
  assert (30-ap)^2 mod 7 ne 0 and (30+ap)^2 mod 7 ne 0;
  print "FIBER",tt,"SQRT13_MOD29",rr[1],"TRACE",ap,"FULL_TWO",true;
 end for;
end for;
assert 576 mod 7 ne 0 and 1296 mod 7 ne 0;
print "VERIFIED_228_MODULAR_CERTIFICATE_DONE";
quit;
