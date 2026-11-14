// family_6_6_certificate.m -- exact infinitude of geometrically simple [6,6] (14 November 2026).
// Bruin-Flynn-Testa's (Z/3)^2 family (BFT 2014, Theorem 6) on the locus s = -4/(r-2)^2, t = r^2/4, where the
// quintic is -q_r*b_r; on the rational curve r = 4/(u^3+3u+2) the cubic b_r has a rational root, which adds a
// second rational 2-torsion class, so [6,6] is marked over Q(u). Every u = 1/2 + 259k (k integral) has
// #J(F_7) = 36 with the marked classes generating J(F_7), and absolutely simple reduction at 37 (#J(F_37) = 1044),
// so J(Q)_tors is exactly [6,6] and J is geometrically simple; the G2-invariants vary. u = 1/2 is census 6.6.a.
// Record: torsion_jac/notes/simple66_2026_11_14.md (with the proof and an independent SymPy check).
// Run (from the repository root):
//     magma -b pipeline/magma/family_6_6_certificate.m < /dev/null > data/knowledge/sources/family_6_6_certificate.log
SetColumns(0);SetMemoryLimit(4*10^9);
Q:=Rationals();K<r>:=FunctionField(Q);Rx<x>:=PolynomialRing(K);
q:=3*x^2+6*x+4*r-1;
b:=8*(r-2)*x^3+3*r*(3*r-8)*x^2+6*r^2*(r-3)*x+r^2*(r^2-6*r+1);
f:=-q*b;
s:=-4/(r-2)^2;t:=r^2/4;scale:=(r-2)^2/r;
h1:=x^2+r*x+t;h2:=x^2+x+r;
g1:=(s-s*t-1)*x^3+3*s*(r-t)*x^2+3*s*r*(r-t)*x-s*t^2+s*r^3+t;
g2:=g1+2*(x^3-t);
g1:=scale*g1;g2:=scale*g2;
lam1:=4*s*scale^2;lam2:=lam1*t;
assert f eq g1^2+lam1*h1^3 and f eq g2^2+lam2*h2^3;
assert IsSquarefree(f) and Degree(f) eq 5;
print "BFT identities over Q(r) verified; f_r=-q_r b_r is smooth generically";
print "g1",g1,"g2",g2,"lambda1",lam1,"lambda2",lam2;

U<u>:=FunctionField(Q);P<X>:=PolynomialRing(U);
ru:=4/(u^3+3*u+2);root:=(u^2-2*u-1)/(u*(u^3+3*u+2));
sub:=hom<K->U|ru>;
mapP:=func<pol|P![sub(cc):cc in Coefficients(pol)]>;
F:=mapP(f);Q1:=mapP(q)/3;H1:=mapP(h1);H2:=mapP(h2);G1:=mapP(g1);G2:=mapP(g2);
assert Evaluate(F,root) eq 0 and IsSquarefree(F);
Q2:=ExactQuotient(F,(X-root)*Q1);Q2:=Q2/LeadingCoefficient(Q2);
assert Degree(Q2) eq 2;
J:=Jacobian(HyperellipticCurve(F));
A:=elt<J|H1,G1 mod H1>;B:=elt<J|H2,G2 mod H2>;
S:=elt<J|X-root,0>;T:=elt<J|Q1,0>;
assert 3*A eq J!0 and A ne J!0 and 3*B eq J!0 and B ne J!0;
assert B ne A and B ne -A;
assert 2*S eq J!0 and S ne J!0 and 2*T eq J!0 and T ne J!0 and S ne T;
assert #{i*A+j*B+k*S+l*T:i,j in [0..2],k,l in [0..1]} eq 36;
D:=elt<J|X+ru/2,Evaluate(G1,-ru/2)>;
assert 2*D eq A and 3*D eq T;
W:=B+S;assert #{i*D+j*W:i,j in [0..5]} eq 36;
print "Two order-6 generators: D=[(-r/2,g1(-r/2))-infinity] and W=B+S; 3D=T";
print "Over Q(u), two independent order-3 and two independent order-2 classes generate [6,6]";
print "r(u)",ru,"rational branch point",root;

R<z>:=PolynomialRing(Q);
spec:=func<pol,u0|R![Evaluate(cc,u0):cc in Coefficients(pol)]>;
f0:=spec(F,1/2);C0:=HyperellipticCurve(f0);
seed:=HyperellipticCurve(1872*z^5-3000*z^4+6969*z^3-1691*z^2+4875*z);
assert IsIsomorphic(C0,seed);
print "u=1/2 is Q-isomorphic to census curve 6.6.a";
for pol in [F,H1,H2,G1,G2,Q1,P!(X-root)] do
 for cc in Coefficients(pol) do
  for p in [7,37] do
   num:=Numerator(cc);den:=Denominator(cc);
   multiplier:=LCM([Denominator(c):c in Coefficients(num)] cat [Denominator(c):c in Coefficients(den)]);
   assert GF(p)!Evaluate(multiplier*den,1/2) ne 0;
  end for;
 end for;
end for;
print "All displayed coefficients and torsion representatives are regular at u=1/2 modulo 7 and 37";
counts:=[];
for p in [7,37] do
 Rp<xx>:=PolynomialRing(GF(p));fp:=Rp!f0;
 assert Degree(fp) eq 5 and IsSquarefree(fp);
 L:=LPolynomial(HyperellipticCurve(fp));chi:=R!(z^4*Evaluate(R!L,1/z));
 n:=Evaluate(L,1);Append(~counts,Integers()!n);
 if p eq 7 then
  J7:=Jacobian(HyperellipticCurve(fp));
  to7:=func<pol|Rp!spec(pol,1/2)>;
  A7:=elt<J7|to7(H1),to7(G1) mod to7(H1)>;B7:=elt<J7|to7(H2),to7(G2) mod to7(H2)>;
  S7:=elt<J7|xx-GF(p)!Evaluate(root,1/2),0>;T7:=elt<J7|to7(Q1),0>;
  assert #{i*A7+j*B7+k*S7+l*T7:i,j in [0..2],k,l in [0..1]} eq 36;
  print "The four marked classes generate all 36 points modulo 7";
 end if;
 printf "p=%o #J=%o chi=%o\n",p,n,chi;
 if p eq 37 then
  assert chi eq z^4-10*z^3+54*z^2-370*z+1369;
  assert IsIrreducible(chi);
  Nk<pi>:=NumberField(chi);
  assert &and[Degree(MinimalPolynomial(pi^j)) eq 4:j in [2..12]];
  print "Geometric simplicity certificate: pi^j has degree 4 for 1<=j<=12";
  V<zeta>:=FunctionField(Q);VV<Y>:=PolynomialRing(V);
  pp:=Y^4-10*Y^3+54*Y^2-370*Y+1369;
  rr:=Resultant(pp,Evaluate(pp,zeta*Y));
  c1:=1369*zeta^4+296*zeta^3-1746*zeta^2+296*zeta+1369;
  c2:=37*zeta^4+20*zeta^3+66*zeta^2+20*zeta+37;
  assert rr eq 37^4*(zeta-1)^4*c1*c2^2;
  ZZ<Zz>:=PolynomialRing(Q);
  for cc in [ZZ!c1,ZZ!c2] do
   for n in [1,2,3,4,5,6,8,10,12] do assert GCD(cc,CyclotomicPolynomial(n)) eq 1;end for;
  end for;
  print "Absolute simplicity independently certified: no quotient of distinct Frobenius roots is a root of unity";
 end if;
end for;
assert counts eq [36,1044] and GCD(counts) eq 36;

function IntegralCurve(pol)
 den:=LCM([Denominator(c):c in Coefficients(pol)]);
 return HyperellipticCurve(den^2*pol);
end function;
C1:=HyperellipticCurve(spec(F,1/2+259));
assert G2Invariants(C0) ne G2Invariants(C1);
print "The geometric moduli map is nonconstant: u=1/2 and u=519/2 have distinct G2 invariants";
for u0 in [Q!1/2,Q!519/2] do
 ff:=spec(F,u0);JJ:=Jacobian(HyperellipticCurve(ff));
 AA:=elt<JJ|spec(H1,u0),spec(G1,u0) mod spec(H1,u0)>;
 BB:=elt<JJ|spec(H2,u0),spec(G2,u0) mod spec(H2,u0)>;
 SS:=elt<JJ|z-Evaluate(root,u0),0>;TT:=elt<JJ|spec(Q1,u0),0>;
 assert #{i*AA+j*BB+k*SS+l*TT:i,j in [0..2],k,l in [0..1]} eq 36;
 for p in [7,37] do assert PolynomialRing(GF(p))!ff eq PolynomialRing(GF(p))!f0;end for;
 assert Invariants(TorsionSubgroup(Jacobian(IntegralCurve(ff)))) eq [6,6];
 printf "Independent TorsionSubgroup check u=%o: [6,6]\n",u0;
end for;
print "FAMILY_66_SIMPLE_CERTIFIED: every u=1/2+259*k, k integral, gives exact [6,6] and geometric simplicity; infinitely many geometric isomorphism classes";
quit;
