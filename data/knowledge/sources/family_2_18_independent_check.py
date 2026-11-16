"""Independent Sage check of the simple [2,18] family, 16 November 2026.
Uses polynomial identities, direct finite-field point counts, and PARI.
The generic torsion proof is by the two displayed principal divisors.
"""
from sage.all import QQ, GF, PolynomialRing, NumberField, pari
R=PolynomialRing(QQ,'k');k=R.gen();K=R.fraction_field();P=PolynomialRing(K,'x');x=P.gen()
s=(k*k+3)/4;r=(k*k+3)/(k-1)**2;lam=(s-1)**2/s**3
s,r,lam=map(K,[s,r,lam])
X=x**3*(x-1)**2;Y=(x-1)**3*(x-r)*(x-s)
A=((x-1)**2+lam*(x**3-(x-1)*(x-r)*(x-s)))/2
B=(lam*(x**3-(x-1)*(x-r)*(x-s))-(x-1)**2)/2
f=A*A-lam*X
assert f==B*B-lam*Y
assert B(r)==B(s)==0 and A(1)==B(1)!=0 and A(0)!=0
assert f.degree()==5 and f.discriminant()!=0 and r!=s
assert lam*r**3==(r-1)**2 and lam*s**3==(s-1)**2
print('SYMBOLIC_DIVISOR_IDENTITIES_VERIFIED')
f0=PolynomialRing(QQ,'x')([c(0) for c in f])
for p in [5,17]:
 counts=[]
 for power in [1,2]:
  F=GF(p**power,'a');fp=PolynomialRing(F,'x')(f0)
  assert fp.gcd(fp.derivative()).degree()==0 and fp.degree()==5
  n=1
  for a in F:
   b=fp(a)
   n+=1 if b==0 else 2 if b.is_square() else 0
  counts.append(n)
 n1,n2=counts; a1=n1-p-1;a2=(a1*a1-p*p-1+n2)//2
 Rt=PolynomialRing(QQ,'T');T=Rt.gen();chi=T**4+a1*T**3+a2*T**2+p*a1*T+p*p
 assert chi(1)==(n1*n1+n2)//2-p
 print('DIRECT_COUNTS',p,counts,'J_ORDER',chi(1),'FROBENIUS',chi)
 if p==5:assert chi(1)==36
 else:
  assert chi==T**4+T**3+16*T**2+17*T+289 and chi(1)==324 and chi.is_irreducible()
  L=NumberField(chi,'pi');pi=L.gen()
  assert all((pi**n).minpoly().degree()==4 for n in range(2,13))
  print('STRICT_FROBENIUS_INDEPENDENTLY_VERIFIED')
  # A direct proof for every extension degree, without a cutoff on n.
  Rz=PolynomialRing(QQ,'z');z=Rz.gen();Rw=PolynomialRing(Rz.fraction_field(),'w');w=Rw.gen()
  chiw=w**4+w**3+16*w**2+17*w+289
  scaled=w**4+z*w**3+16*z**2*w**2+17*z**3*w+289*z**4
  h1=289*z**4+527*z**3+800*z**2+527*z+289
  h2=17*z**4+18*z**3+3*z**2+18*z+17
  assert chiw.resultant(scaled)==17**4*(z-1)**4*h1*h2**2
  assert h1.is_irreducible() and h2.is_irreducible()
  assert not h1.monic().is_cyclotomic() and not h2.monic().is_cyclotomic()
  print('NO_DISTINCT_FROBENIUS_ROOT_RATIO_IS_A_ROOT_OF_UNITY')
print('SAGE_INDEPENDENT_CHECK_PASSED')
