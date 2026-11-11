"""Exact finite controls for the uniform order-256 trace proof.

These checks supplement the written identities; they do not bound q.
"""
from fractions import Fraction as Q
from itertools import product


def det(a):
    return a[0]*a[3]-a[1]*a[2]


def mul(a,b):
    return (a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
            a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3])


def inv(a):
    d=det(a)
    return (Q(a[3],d), Q(-a[1],d), Q(-a[2],d), Q(a[0],d))


F2=list(product(range(2), repeat=4))
ann=[]
for k in F2:
    if all((mul(m,k)[0]+mul(m,k)[3]) % 2 == 0 for m in F2):
        ann.append(k)
assert ann==[(0,0,0,0)]
print('TRACE_PAIRING_NONDEGENERATE')

gl8=sum(det(a)%2==1 for a in product(range(8),repeat=4))
h8=sum(det(a)%2==1 and a[0]==1 and a[2]==0 and a[1]%2==0 and a[3]%2==1
       for a in product(range(8),repeat=4))
assert (gl8,h8,gl8//h8)==(1536,16,96)
print(f'ELLIPTIC_STABILIZER_INDEX gl2_mod8={gl8} stabilizer={h8} index={gl8//h8}')

det_checks=0
for q in [2,4,8,16,32,64,128,256,512,1024]:
    for a,v in product(F2,repeat=2):
        g=(1+2*a[0],2*a[1],2*a[2],1+2*a[3])
        h=tuple(g[i]+q*v[i] for i in range(4))
        assert (det(h)-det(g))%q==0
        assert ((det(h)-det(g))//q-(v[0]+v[3]))%2==0
        det_checks+=1
print(f'DETERMINANT_TRACE_CONTROLS checks={det_checks}')

D=(8,0,0,2)
I=(1,0,0,1)
p_list=[p for p in product(range(4),repeat=4) if det(p)%2]
matrix_checks=0
fractional_checks=0
for a,b,c,d in product(range(4),repeat=4):
    delta=a*d-b*c
    if not delta%2:
        continue
    L=(Q(a),Q(c,4),Q(4*b),Q(d))
    expected=(Q(8*d,delta),Q(-2*c,delta),Q(-8*b,delta),Q(2*a,delta))
    assert mul(D,inv(L))==expected
    for p in p_list:
        N=mul(inv(p),L)
        dni=mul(D,inv(N))
        assert dni==mul(expected,p)
        assert all(x.denominator%2 and x.numerator%2==0 for x in dni)
        matrix_checks+=1
        fractional_checks+=int(c%4!=0)
print(f'RATIONAL_MATRIX_CONTROLS checks={matrix_checks} fractional_upper_right={fractional_checks}')
print('ORDER256_TRACE_CERTIFICATE_DONE')
